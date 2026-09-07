# ===== Reproducible coloc.susie analysis: UKB and FinnGen =====

rm(list = ls(all.names = TRUE))
gc()

library(data.table)
library(dplyr)
library(susieR)
library(coloc)
library(Matrix)
library(tibble)
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
# Inputs
# -------------------------------------------------------------------------

eqtl_file <- Sys.getenv(
  "COLOC_EQTL_FILE",
  unset = file.path(
    repo_root,
    "data",
    "GTEx_v8",
    "myLiteCisEqtl.txt"
  )
)

ukb_gwas <- Sys.getenv(
  "SUSIE_UKB_GWAS",
  unset = file.path(
    repo_root,
    "data",
    "GWAS",
    "UKB",
    "UKB_SMAD3_250kb.txt"
  )
)

fin_gwas <- Sys.getenv(
  "SUSIE_FINNGEN_GWAS",
  unset = file.path(
    repo_root,
    "data",
    "GWAS",
    "FinnGen",
    "Finn_SMAD3_250kb.txt"
  )
)

ukb_bim <- Sys.getenv(
  "SUSIE_UKB_BIM",
  unset = file.path(
    repo_root,
    "data",
    "1000Genomes",
    "UKB_SMAD3_1000G_EUR.bim"
  )
)

ukb_ld <- Sys.getenv(
  "SUSIE_UKB_LD",
  unset = file.path(
    repo_root,
    "data",
    "1000Genomes",
    "UKB_SMAD3_1000G_EUR_LD.ld"
  )
)

fin_bim <- Sys.getenv(
  "SUSIE_FINNGEN_BIM",
  unset = file.path(
    repo_root,
    "data",
    "1000Genomes",
    "Finn_SMAD3_1000G_EUR.bim"
  )
)

fin_ld <- Sys.getenv(
  "SUSIE_FINNGEN_LD",
  unset = file.path(
    repo_root,
    "data",
    "1000Genomes",
    "Finn_SMAD3_1000G_EUR_LD.ld"
  )
)

required_files <- c(
  eqtl_file,
  ukb_gwas,
  fin_gwas,
  ukb_bim,
  ukb_ld,
  fin_bim,
  fin_ld
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  stop(
    "Required coloc.susie input file(s) not found:\n",
    paste(missing_files, collapse = "\n")
  )
}

# -------------------------------------------------------------------------
# Outputs
# -------------------------------------------------------------------------

processed_dir <- file.path(
  repo_root,
  "processed_results",
  "08_colocalization"
)

source_dir <- file.path(
  repo_root,
  "source_data",
  "figure_source_data",
  "Figure_12"
)

figdir <- file.path(
  repo_root,
  "outputs",
  "coloc_susie"
)

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# Cohort definitions
# -------------------------------------------------------------------------

cohorts <- list(
  UKB = list(
    gwas = ukb_gwas,
    bim = ukb_bim,
    ld = ukb_ld,
    ncase = 28707,
    ncontrol = 343103,
    smr_top_snp = "rs7181556"
  ),
  FinnGen = list(
    gwas = fin_gwas,
    bim = fin_bim,
    ld = fin_ld,
    ncase = 17096,
    ncontrol = 190557,
    smr_top_snp = "rs11315136"
  )
)

N_eqtl <- 663
sdY_eqtl <- 1

cat("coloc.susie input preflight passed.\n")
cat("UKB total N:", cohorts$UKB$ncase + cohorts$UKB$ncontrol, "\n")
cat("FinnGen total N:", cohorts$FinnGen$ncase + cohorts$FinnGen$ncontrol, "\n")
cat("GTEx eQTL N:", N_eqtl, "\n")
cat("UKB SMR top SNP:", cohorts$UKB$smr_top_snp, "\n")
cat("FinnGen SMR top SNP:", cohorts$FinnGen$smr_top_snp, "\n")

# =========================================================================
# Prepare SMAD3 eQTL data once
# =========================================================================

eqtl <- fread(eqtl_file)

eqtl <- eqtl[
  BP > (67357940 - 1e6) &
  BP < (67487507 + 1e6)
]

setnames(
  eqtl,
  old = c("Chr", "se"),
  new = c("CHR", "SE"),
  skip_absent = TRUE
)

eqtl <- eqtl[
  ,
  .(SNP, CHR, BP, A1, A2, Gene, b, SE, p)
]

eqtl[, `:=`(
  b = as.numeric(b),
  SE = as.numeric(SE),
  p = as.numeric(p),
  BP = as.numeric(BP)
)]

eqtl <- eqtl[
  complete.cases(eqtl) &
  Gene == "SMAD3"
]

eqtl[, `:=`(
  varbeta = SE^2,
  z = b / SE
)]

# =========================================================================
# coloc.susie workflow for one GWAS cohort
# =========================================================================

run_coloc_susie <- function(
    study,
    gwas_file,
    bim_file,
    ld_file,
    ncase,
    ncontrol,
    smr_top_snp) {

  cat("\n===== Running coloc.susie:", study, "=====\n")

  ntotal <- ncase + ncontrol
  case_fraction <- ncase / ntotal

  # -----------------------------------------------------------------------
  # GWAS
  # -----------------------------------------------------------------------

  gwas <- fread(gwas_file)
  gwas$rsID <- gwas$SNP

  gwas <- gwas %>%
    arrange(rsID, p) %>%
    group_by(rsID) %>%
    slice(1) %>%
    ungroup()

  gwas <- gwas[
    complete.cases(gwas) &
    gwas$b != 0,
  ]

  gwas <- gwas %>%
    mutate(
      s = case_fraction,
      varbeta = se^2
    ) %>%
    dplyr::select(
      rsID, A1, A2, freq,
      b, se, p, pos,
      s, varbeta
    )

  # -----------------------------------------------------------------------
  # BIM + LD
  # -----------------------------------------------------------------------

  bim <- fread(
    bim_file,
    col.names = c(
      "CHR", "rsID", "CM", "BP_ld",
      "A1_ld", "A2_ld"
    )
  )

  dup_bim <- duplicated(bim$rsID)
  bim2 <- bim[!dup_bim, ]

  ld_matrix <- as.matrix(
    read.table(
      ld_file,
      header = FALSE,
      check.names = FALSE
    )
  )

  ld_matrix <- ld_matrix[
    !dup_bim,
    !dup_bim,
    drop = FALSE
  ]

  rownames(ld_matrix) <- bim2$rsID
  colnames(ld_matrix) <- bim2$rsID

  # Remove invalid LD rows/columns
  diag_vals <- diag(ld_matrix)

  bad_snps <- names(diag_vals)[
    diag_vals <= 0 |
    is.na(diag_vals)
  ]

  good_snps <- setdiff(
    rownames(ld_matrix),
    bad_snps
  )

  ld_matrix_clean <- ld_matrix[
    good_snps,
    good_snps,
    drop = FALSE
  ]

  bim_clean <- bim2[
    bim2$rsID %in% good_snps,
  ]

  # -----------------------------------------------------------------------
  # Align GWAS alleles with LD reference
  # -----------------------------------------------------------------------

  merged <- gwas %>%
    inner_join(
      bim_clean %>%
        dplyr::select(
          rsID,
          A1_ld,
          A2_ld
        ),
      by = "rsID"
    )

  same <- (
    merged$A1 == merged$A1_ld &
    merged$A2 == merged$A2_ld
  )

  flipped <- (
    merged$A1 == merged$A2_ld &
    merged$A2 == merged$A1_ld
  )

  keep <- same | flipped
  merged <- merged[keep, ]

  merged$beta_aligned <- merged$b
  merged$freq_aligned <- merged$freq

  flip_idx <- which(flipped[keep])

  merged$beta_aligned[flip_idx] <-
    -merged$beta_aligned[flip_idx]

  merged$freq_aligned[flip_idx] <-
    1 - merged$freq_aligned[flip_idx]

  # MAF filter
  merged$MAF <- pmin(
    merged$freq_aligned,
    1 - merged$freq_aligned
  )

  merged <- merged[
    merged$MAF >= 0.01,
  ]

  # -----------------------------------------------------------------------
  # Final GWAS + LD SNP set
  # -----------------------------------------------------------------------

  final_snps <- intersect(
    rownames(ld_matrix_clean),
    merged$rsID
  )

  ld_matrix_final <- ld_matrix_clean[
    final_snps,
    final_snps,
    drop = FALSE
  ]

  ld_pd <- as.matrix(
    nearPD(
      ld_matrix_final,
      corr = TRUE
    )$mat
  )

  rownames(ld_pd) <- final_snps
  colnames(ld_pd) <- final_snps

  merged_final <- merged[
    merged$rsID %in% final_snps,
  ]

  merged_final <- merged_final[
    match(
      final_snps,
      merged_final$rsID
    ),
  ]

  stopifnot(
    all(rownames(ld_pd) == merged_final$rsID)
  )

  # -----------------------------------------------------------------------
  # GWAS + LD + SMAD3 eQTL common SNPs
  # -----------------------------------------------------------------------

  common_snps <- intersect(
    merged_final$rsID,
    eqtl$SNP
  )

  R_common <- ld_pd[
    common_snps,
    common_snps,
    drop = FALSE
  ]

  gwas_sub <- data.frame(
    SNP = merged_final$rsID,
    p = merged_final$p,
    beta = merged_final$beta_aligned,
    se = merged_final$se,
    varbeta = merged_final$se^2,
    pos = merged_final$pos,
    s = case_fraction,
    stringsAsFactors = FALSE
  )

  gwas_sub <- gwas_sub[
    gwas_sub$SNP %in% common_snps,
  ]

  gwas_sub <- gwas_sub[
    match(
      common_snps,
      gwas_sub$SNP
    ),
  ]

  eqtl_sub <- eqtl[
    eqtl$SNP %in% common_snps,
  ]

  eqtl_sub <- eqtl_sub[
    match(
      common_snps,
      eqtl_sub$SNP
    ),
  ]

  stopifnot(
    all(rownames(R_common) == common_snps),
    all(gwas_sub$SNP == common_snps),
    all(eqtl_sub$SNP == common_snps)
  )

  cat(
    study,
    "common GWAS+LD+eQTL SNPs:",
    length(common_snps),
    "\n"
  )

  # -----------------------------------------------------------------------
  # GWAS SuSiE
  # -----------------------------------------------------------------------

  dataset1 <- list(
    beta = setNames(
      gwas_sub$beta,
      gwas_sub$SNP
    ),
    varbeta = setNames(
      gwas_sub$varbeta,
      gwas_sub$SNP
    ),
    s = case_fraction,
    type = "cc",
    snp = gwas_sub$SNP,
    LD = R_common,
    position = gwas_sub$pos,
    N = ntotal
  )

  coloc:::check_dataset(
    dataset1,
    1
  )

  S3 <- runsusie(
    dataset1,
    nref = 503,
    coverage = 0.90
  )

  # -----------------------------------------------------------------------
  # eQTL SuSiE
  # -----------------------------------------------------------------------

  dataset2 <- list(
    beta = setNames(
      eqtl_sub$b,
      eqtl_sub$SNP
    ),
    varbeta = setNames(
      eqtl_sub$varbeta,
      eqtl_sub$SNP
    ),
    N = N_eqtl,
    sdY = sdY_eqtl,
    type = "quant",
    snp = eqtl_sub$SNP,
    LD = R_common,
    position = eqtl_sub$BP
  )

  coloc:::check_dataset(
    dataset2,
    2
  )

  S4 <- runsusie(
    dataset2,
    nref = 503,
    coverage = 0.90
  )

  # -----------------------------------------------------------------------
  # coloc.susie
  # -----------------------------------------------------------------------

  coloc_res <- coloc.susie(
    S3,
    S4
  )

  coloc_summary <- as.data.table(
    coloc_res$summary
  )

  # PP.H4 is a Bayesian posterior probability, not a P value.
  # PP.H4 > 0.75 is used as the colocalization-support criterion;
  # frequentist multiple-testing correction is therefore not applicable.
  coloc_supported <- coloc_summary[
    PP.H4.abf > 0.75
  ]

  # Add study for combined processed output
  coloc_summary[, study := study]
  coloc_supported[, study := study]

  setcolorder(
    coloc_summary,
    c(
      "study",
      setdiff(
        names(coloc_summary),
        "study"
      )
    )
  )

  setcolorder(
    coloc_supported,
    c(
      "study",
      setdiff(
        names(coloc_supported),
        "study"
      )
    )
  )

  # -----------------------------------------------------------------------
  # Figure 12 SNP-level source data
  # -----------------------------------------------------------------------

  locus_source <- data.table(
    study = study,
    SNP = common_snps,
    position = gwas_sub$pos,
    GWAS_p = gwas_sub$p,
    eQTL_p = eqtl_sub$p,
    GWAS_beta = gwas_sub$beta,
    GWAS_se = gwas_sub$se,
    eQTL_beta = eqtl_sub$b,
    eQTL_se = eqtl_sub$SE
  )

  # Local LD correlation with SMR top SNP
  if (smr_top_snp %in% common_snps) {

    locus_source[, LD_r_with_SMR_top :=
      as.numeric(
        R_common[
          smr_top_snp,
          common_snps
        ]
      )
    ]

  } else {

    locus_source[, LD_r_with_SMR_top := NA_real_]
  }

  locus_source[, SMR_top_SNP :=
    SNP == smr_top_snp
  ]

  # Mark SNPs from PP.H4-supported coloc.susie components.
  hit_snps <- unique(c(
    coloc_supported$hit1,
    coloc_supported$hit2
  ))

  hit_snps <- hit_snps[
    !is.na(hit_snps)
  ]

  locus_source[, coloc_susie_hit :=
    SNP %in% hit_snps
  ]

  fwrite(
    locus_source,
    file.path(
      source_dir,
      paste0(
        "Figure12_",
        study,
        "_locus_source_data.csv"
      )
    )
  )

  fwrite(
    coloc_summary,
    file.path(
      source_dir,
      paste0(
        "Figure12_",
        study,
        "_coloc_susie_summary_source_data.csv"
      )
    )
  )

  # -----------------------------------------------------------------------
  # Processed result
  # -----------------------------------------------------------------------

  fwrite(
    coloc_supported,
    file.path(
      processed_dir,
      paste0(
        study,
        "_coloc_susie_reproduced.csv"
      )
    )
  )

  cat(
    study,
    "PP.H4-supported coloc.susie rows:",
    nrow(coloc_supported),
    "\n"
  )

  if (nrow(coloc_supported) > 0) {
    print(coloc_supported)
  }

  invisible(
    list(
      S3 = S3,
      S4 = S4,
      coloc = coloc_res,
      supported = coloc_supported,
      locus_source = locus_source
    )
  )
}

# =========================================================================
# Run both cohorts
# =========================================================================

ukb_coloc_susie <- run_coloc_susie(
  study = "UKB",
  gwas_file = cohorts$UKB$gwas,
  bim_file = cohorts$UKB$bim,
  ld_file = cohorts$UKB$ld,
  ncase = cohorts$UKB$ncase,
  ncontrol = cohorts$UKB$ncontrol,
  smr_top_snp = cohorts$UKB$smr_top_snp
)

finngen_coloc_susie <- run_coloc_susie(
  study = "FinnGen",
  gwas_file = cohorts$FinnGen$gwas,
  bim_file = cohorts$FinnGen$bim,
  ld_file = cohorts$FinnGen$ld,
  ncase = cohorts$FinnGen$ncase,
  ncontrol = cohorts$FinnGen$ncontrol,
  smr_top_snp = cohorts$FinnGen$smr_top_snp
)

cat(
  "\nAll coloc.susie analyses completed.\n"
)