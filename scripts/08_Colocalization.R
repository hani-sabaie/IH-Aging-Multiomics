# ===== Standard coloc.abf analysis: UKB and FinnGen =====

rm(list = ls(all.names = TRUE))
gc()

library(data.table)
library(dplyr)
library(coloc)
library(purrr)

# -------------------------------------------------------------------------
# Locate repository root
# -------------------------------------------------------------------------

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (length(file_arg) == 1) {
  script_path <- sub("^--file=", "", file_arg)
  script_dir <- dirname(normalizePath(script_path))
  repo_root <- normalizePath(file.path(script_dir, ".."))
} else {
  repo_root <- normalizePath(".")
}

# -------------------------------------------------------------------------
# Input paths
#
# Defaults assume the repository-local data/ structure documented in
# DATA_SOURCES.md. Environment variables allow use of locally stored
# third-party data without copying them into the repository.
# -------------------------------------------------------------------------

eqtl_file <- Sys.getenv(
  "COLOC_EQTL_FILE",
  unset = file.path(repo_root, "data", "GTEx_v8", "myLiteCisEqtl.txt")
)

ukb_file <- Sys.getenv(
  "COLOC_UKB_GWAS",
  unset = file.path(repo_root, "data", "GWAS", "UKB", "UKB_gwas_for_smr.txt")
)

finngen_file <- Sys.getenv(
  "COLOC_FINNGEN_GWAS",
  unset = file.path(repo_root, "data", "GWAS", "FinnGen", "Finn_gwas_for_smr.txt")
)

outdir <- file.path(
  repo_root,
  "processed_results",
  "08_colocalization"
)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

for (f in c(eqtl_file, ukb_file, finngen_file)) {
  if (!file.exists(f)) {
    stop("Required input file not found: ", f)
  }
}

# -------------------------------------------------------------------------
# GTEx v8 subcutaneous-adipose eQTL data
# SMAD3 region +/- 1 Mb, matching the original analysis
# -------------------------------------------------------------------------

smad3_start <- 67357940
smad3_end   <- 67487507
cis_window  <- 1e6

eqtl <- fread(
  eqtl_file,
  select = c(
    "SNP", "BP", "A1", "A2",
    "Probe", "Gene", "b", "SE", "p"
  ),
  data.table = FALSE
)

eqtl_region <- eqtl %>%
  filter(
    BP > (smad3_start - cis_window),
    BP < (smad3_end + cis_window)
  ) %>%
  transmute(
    SNP,
    A1_eqtl = A1,
    A2_eqtl = A2,
    beta_eqtl = b,
    se_eqtl = SE,
    p_eqtl = p,
    gene = Gene
  )

# -------------------------------------------------------------------------
# Function for one GWAS cohort
# -------------------------------------------------------------------------

run_coloc <- function(
    gwas_file,
    study,
    n_total,
    n_cases,
    n_eqtl = 663) {

  message("Reading ", study, " GWAS...")

  gwas <- fread(
    gwas_file,
    select = c("SNP", "A1", "A2", "freq", "b", "se", "p"),
    data.table = FALSE
  ) %>%
    mutate(
      MAF_gwas = ifelse(freq <= 0.5, freq, 1 - freq)
    ) %>%
    transmute(
      SNP,
      A1_gwas = A1,
      A2_gwas = A2,
      beta_gwas = b,
      se_gwas = se,
      p_gwas = p,
      MAF_gwas
    )

  df_coloc <- inner_join(eqtl_region, gwas, by = "SNP")

  # Keep SNPs whose alleles can be aligned
  same_alleles <-
    (
      df_coloc$A1_eqtl == df_coloc$A1_gwas &
      df_coloc$A2_eqtl == df_coloc$A2_gwas
    ) |
    (
      df_coloc$A1_eqtl == df_coloc$A2_gwas &
      df_coloc$A2_eqtl == df_coloc$A1_gwas
    )

  df_coloc <- df_coloc[same_alleles, , drop = FALSE]

  # Flip GWAS beta when alleles are reversed
  flip <-
    df_coloc$A1_eqtl == df_coloc$A2_gwas &
    df_coloc$A2_eqtl == df_coloc$A1_gwas

  df_coloc$beta_gwas[flip] <- -df_coloc$beta_gwas[flip]

  genes <- unique(df_coloc$gene)

  results <- lapply(genes, function(g) {

    tmp <- df_coloc %>%
      filter(gene == g) %>%
      arrange(SNP) %>%
      distinct(SNP, .keep_all = TRUE) %>%
      filter(
        !is.na(beta_gwas),
        !is.na(se_gwas),
        !is.na(beta_eqtl),
        !is.na(se_eqtl),
        !is.na(MAF_gwas),
        p_eqtl > 0,
        p_gwas > 0,
        MAF_gwas > 0,
        MAF_gwas < 1
      )

    if (nrow(tmp) < 1) {
      return(NULL)
    }

    D1 <- list(
      snp = tmp$SNP,
      beta = tmp$beta_gwas,
      varbeta = tmp$se_gwas^2,
      type = "cc",
      s = n_cases / n_total,
      N = n_total
    )

    D2 <- list(
      snp = tmp$SNP,
      beta = tmp$beta_eqtl,
      varbeta = tmp$se_eqtl^2,
      type = "quant",
      N = n_eqtl
    )

    coloc.abf(
      dataset1 = D1,
      dataset2 = D2,
      MAF = tmp$MAF_gwas
    )
  })

  names(results) <- genes

  coloc_sum <- map_dfr(
    results,
    function(res) {

      if (is.null(res)) {
        return(NULL)
      }

      snp_tab <- res$results
      snp_tab <- snp_tab[order(-snp_tab$SNP.PP.H4), ]

      data.frame(
        nsnps = unname(res$summary["nsnps"]),
        pph0 = unname(res$summary["PP.H0.abf"]),
        pph1 = unname(res$summary["PP.H1.abf"]),
        pph2 = unname(res$summary["PP.H2.abf"]),
        pph3 = unname(res$summary["PP.H3.abf"]),
        pph4 = unname(res$summary["PP.H4.abf"]),
        top_snp = snp_tab$snp[1],
        top_snp_pph4 = snp_tab$SNP.PP.H4[1]
      )
    },
    .id = "gene"
  )

  coloc_sum %>%
    mutate(
      study = study,
      tissue = "Adipose Subcutaneous",
      .before = gene
    )
}

# -------------------------------------------------------------------------
# Run both discovery and replication cohorts
# -------------------------------------------------------------------------

ukb_coloc <- run_coloc(
  gwas_file = ukb_file,
  study = "UKB",
  n_total = 371810,
  n_cases = 28707
)

finngen_coloc <- run_coloc(
  gwas_file = finngen_file,
  study = "FinnGen",
  n_total = 207653,
  n_cases = 17096
)

coloc_all <- bind_rows(
  ukb_coloc,
  finngen_coloc
)

# PP.H4 > 0.75: colocalization criterion used in the manuscript
coloc_sig <- coloc_all %>%
  filter(pph4 > 0.75)

# -------------------------------------------------------------------------
# Save author-generated processed results
# -------------------------------------------------------------------------

all_file <- file.path(
  outdir,
  "COLOC_all_studies.csv"
)

sig_file <- file.path(
  outdir,
  "COLOC_all_studies_sig.csv"
)

fwrite(coloc_all, all_file)
fwrite(coloc_sig, sig_file)

cat("Total COLOC rows:", nrow(coloc_all), "\n")
cat("PP.H4 > 0.75 rows:", nrow(coloc_sig), "\n")
cat("All results:", all_file, "\n")
cat("Significant results:", sig_file, "\n")