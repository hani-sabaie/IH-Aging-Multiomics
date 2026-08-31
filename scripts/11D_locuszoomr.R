# ===== Reproducible locuszoom / eQTL visualization for Figure 13A-D =====

rm(list = ls(all.names = TRUE))
gc()

library(data.table)
library(locuszoomr)
library(EnsDb.Hsapiens.v75)
library(ggplot2)

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

source_dir <- file.path(
  repo_root,
  "source_data",
  "figure_source_data",
  "Figure_13"
)

figdir <- file.path(
  repo_root,
  "outputs",
  "Figure_13"
)

dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# LDlink API token
# -------------------------------------------------------------------------

ldlink_token <- Sys.getenv("LDLINK_TOKEN")

if (!nzchar(ldlink_token)) {
  stop(
    "LDLINK_TOKEN is not set. ",
    "Set it as an environment variable before running this script."
  )
}

# -------------------------------------------------------------------------
# Figure input data
# -------------------------------------------------------------------------

ukb_file <- file.path(
  source_dir,
  "Figure13AC_UKB_GWAS_locus_source_data.tsv"
)

finngen_file <- file.path(
  source_dir,
  "Figure13BD_FinnGen_GWAS_locus_source_data.tsv"
)

if (!file.exists(ukb_file)) {
  stop("UKB Figure 13 input file not found: ", ukb_file)
}

if (!file.exists(finngen_file)) {
  stop("FinnGen Figure 13 input file not found: ", finngen_file)
}

# -------------------------------------------------------------------------
# Helper
# -------------------------------------------------------------------------

run_locus_panels <- function(
    study,
    input_file,
    index_snp,
    labels,
    label_x,
    pcutoff,
    locus_panel,
    eqtl_panel) {

  cat("\n===== Running Figure 13:", study, "=====\n")

  dat <- fread(input_file)

  loc <- locus(
    data = dat,
    gene = "SMAD3",
    fix_window = 1e6,
    ens_db = "EnsDb.Hsapiens.v75",
    LD = "r2",
    index_snp = index_snp
  )

  cat(study, "locus SNPs:", nrow(loc$data), "\n")

  # -----------------------------------------------------------------------
  # LD from LDlink / 1000 Genomes EUR
  # -----------------------------------------------------------------------

  loc <- link_LD(
    loc,
    token = ldlink_token,
    genome_build = "grch37"
  )

  cat(
    study,
    "SNPs with non-missing LD:",
    sum(!is.na(loc$data$ld)),
    "\n"
  )

  # -----------------------------------------------------------------------
  # Recombination data
  # -----------------------------------------------------------------------

  loc <- link_recomb(
    loc,
    genome = "hg19"
  )

  # Save numerical source data behind locus plot
  fwrite(
    as.data.table(loc$data),
    file.path(
      source_dir,
      paste0(
        "Figure13",
        locus_panel,
        "_",
        study,
        "_locusplot_source_data.tsv"
      )
    ),
    sep = "\t"
  )

  fwrite(
    as.data.table(loc$recomb),
    file.path(
      source_dir,
      paste0(
        "Figure13",
        locus_panel,
        "_",
        study,
        "_recombination_source_data.tsv"
      )
    ),
    sep = "\t"
  )

  # -----------------------------------------------------------------------
  # Locus plot
  # -----------------------------------------------------------------------

  png(
    file.path(
      figdir,
      paste0(
        "Figure13",
        locus_panel,
        "_",
        study,
        "_locusplot.png"
      )
    ),
    width = 3000,
    height = 3000,
    res = 300
  )

  print(
    locus_plot(
      loc,
      labels = labels,
      label_x = label_x,
      recomb_offset = 0.1,
      highlight = "SMAD3",
      pcutoff = pcutoff
    )
  )

  dev.off()

  # -----------------------------------------------------------------------
  # eQTL data from LDlink
  # -----------------------------------------------------------------------

  loc <- link_eqtl(
    loc,
    token = ldlink_token
  )

  eqtl_source <- as.data.table(loc$LDexp)[
    Gene_Symbol == "SMAD3" &
    Tissue == "Adipose - Subcutaneous"
  ]

  fwrite(
    eqtl_source,
    file.path(
      source_dir,
      paste0(
        "Figure13",
        eqtl_panel,
        "_",
        study,
        "_SMAD3_Adipose_Subcutaneous_eQTL_source_data.tsv"
      )
    ),
    sep = "\t"
  )

  cat(
    study,
    "SMAD3 subcutaneous adipose eQTL rows:",
    nrow(eqtl_source),
    "\n"
  )

  # -----------------------------------------------------------------------
  # eQTL overlay plot
  # -----------------------------------------------------------------------

  png(
    file.path(
      figdir,
      paste0(
        "Figure13",
        eqtl_panel,
        "_",
        study,
        "_eQTL_overlay.png"
      )
    ),
    width = 3000,
    height = 3000,
    res = 300
  )

  print(
    overlay_plot(
      loc,
      eqtl_gene = "SMAD3",
      tissue = "Adipose - Subcutaneous",
      labels = labels,
      label_x = label_x,
      highlight = "SMAD3",
      recomb_col = NA,
      pcutoff = NA
    )
  )

  dev.off()

  invisible(loc)
}

# =========================================================================
# UK Biobank — Figure 13A and 13C
# =========================================================================

loc_ukb <- run_locus_panels(
  study = "UKB",
  input_file = ukb_file,
  index_snp = "rs7181556",
  labels = c(
    "rs7181556",
    "rs10152595",
    "rs35874463",
    "rs12912045",
    "rs7181877"
  ),
  label_x = c(4, -5, 5),
  pcutoff = 1e-6,
  locus_panel = "A",
  eqtl_panel = "C"
)

# =========================================================================
# FinnGen — Figure 13B and 13D
# =========================================================================

loc_finngen <- run_locus_panels(
  study = "FinnGen",
  input_file = finngen_file,
  index_snp = "rs11315136",
  labels = c(
    "rs11315136",
    "rs1965269",
    "rs3784681"
  ),
  label_x = c(-4, -5, 5),
  pcutoff = 1e-4,
  locus_panel = "B",
  eqtl_panel = "D"
)

cat("\nFigure 13A-D source-data generation completed.\n")
