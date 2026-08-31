# ===== Reproducible SuSiE fine-mapping for UKB and FinnGen =====

rm(list = ls(all.names = TRUE))
gc()

library(data.table)
library(dplyr)
library(susieR)
library(coloc)
library(Matrix)
library(tibble)
library(purrr)
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

# -------------------------------------------------------------------------
# Input paths
#
# Third-party/local inputs are supplied through environment variables.
# Repository-local defaults follow DATA_SOURCES.md.
# -------------------------------------------------------------------------

ukb_gwas <- Sys.getenv(
  "SUSIE_UKB_GWAS",
  unset = file.path(repo_root, "data", "GWAS", "UKB", "UKB_SMAD3_250kb.txt")
)

fin_gwas <- Sys.getenv(
  "SUSIE_FINNGEN_GWAS",
  unset = file.path(repo_root, "data", "GWAS", "FinnGen", "Finn_SMAD3_250kb.txt")
)

ukb_bim <- Sys.getenv(
  "SUSIE_UKB_BIM",
  unset = file.path(repo_root, "data", "1000Genomes", "UKB_SMAD3_1000G_EUR.bim")
)

ukb_ld <- Sys.getenv(
  "SUSIE_UKB_LD",
  unset = file.path(repo_root, "data", "1000Genomes", "UKB_SMAD3_1000G_EUR_LD.ld")
)

fin_bim <- Sys.getenv(
  "SUSIE_FINNGEN_BIM",
  unset = file.path(repo_root, "data", "1000Genomes", "Finn_SMAD3_1000G_EUR.bim")
)

fin_ld <- Sys.getenv(
  "SUSIE_FINNGEN_LD",
  unset = file.path(repo_root, "data", "1000Genomes", "Finn_SMAD3_1000G_EUR_LD.ld")
)

required_files <- c(
  ukb_gwas,
  fin_gwas,
  ukb_bim,
  ukb_ld,
  fin_bim,
  fin_ld
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "Required SuSiE input file(s) not found:\n",
    paste(missing_files, collapse = "\n")
  )
}

# -------------------------------------------------------------------------
# Output directories
# -------------------------------------------------------------------------

processed_dir <- file.path(
  repo_root,
  "processed_results",
  "09_SuSiE"
)

source_dir <- file.path(
  repo_root,
  "source_data",
  "figure_source_data",
  "Figure_11"
)

figdir <- file.path(
  repo_root,
  "outputs",
  "SuSiE"
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
    outlier = NA_character_
  ),
  FinnGen = list(
    gwas = fin_gwas,
    bim = fin_bim,
    ld = fin_ld,
    ncase = 17096,
    ncontrol = 190557,
    outlier = "rs141195834"
  )
)

cat("SuSiE input preflight passed.\n")
cat("UKB total N:", cohorts$UKB$ncase + cohorts$UKB$ncontrol, "\n")
cat("FinnGen total N:", cohorts$FinnGen$ncase + cohorts$FinnGen$ncontrol, "\n")
cat("FinnGen diagnostic outlier:", cohorts$FinnGen$outlier, "\n")

# =========================================================================
# Main SuSiE workflow for one cohort
# =========================================================================

run_susie <- function(
    study,
    gwas_file,
    bim_file,
    ld_file,
    ncase,
    ncontrol,
    outlier = NA_character_) {

  cat("\n===== Running", study, "=====\n")

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

  # -----------------------------------------------------------------------
  # Basic LD QC
  # -----------------------------------------------------------------------

  diag_vals <- diag(ld_matrix)

  bad_snps <- names(diag_vals)[
    diag_vals <= 0 | is.na(diag_vals)
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
  # Allele alignment
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

  A1_g <- merged$A1
  A2_g <- merged$A2
  A1_r <- merged$A1_ld
  A2_r <- merged$A2_ld

  same <- (
    A1_g == A1_r &
      A2_g == A2_r
  )

  flipped <- (
    A1_g == A2_r &
      A2_g == A1_r
  )

  keep <- same | flipped

  merged <- merged[keep, ]

  merged$beta_aligned <- merged$b
  merged$freq_aligned <- merged$freq

  flip_idx <- which(flipped[keep])
  flip_snps <- merged$rsID[flip_idx]

  # -----------------------------------------------------------------------
  # Figure 11A/B source data: marginal z before allele flipping
  # -----------------------------------------------------------------------

  z_raw <- merged$beta_aligned / merged$se

  marginal_source <- data.table(
    SNP_index = seq_along(z_raw),
    SNP = merged$rsID,
    z_raw = z_raw,
    allele_flipped = merged$rsID %in% flip_snps
  )

  fwrite(
    marginal_source,
    file.path(
      source_dir,
      paste0(
        "Figure11_",
        study,
        "_marginal_z_source_data.csv"
      )
    )
  )

  # Plot marginal z
  png(
    file.path(
      figdir,
      paste0(study, "_marginal_associations.png")
    ),
    width = 3000,
    height = 1500,
    res = 300
  )

  plot(
    marginal_source$SNP_index,
    marginal_source$z_raw,
    pch = 16,
    col = "#767676",
    main = paste(study, "Marginal Associations"),
    xlab = "SNP index along locus",
    ylab = "z-scores"
  )

  if (any(marginal_source$allele_flipped)) {
    points(
      marginal_source$SNP_index[
        marginal_source$allele_flipped
      ],
      marginal_source$z_raw[
        marginal_source$allele_flipped
      ],
      col = "yellow",
      pch = 16,
      cex = 1.2
    )
  }

  dev.off()

  # Actually align effects
  merged$beta_aligned[flip_idx] <-
    -merged$beta_aligned[flip_idx]

  merged$freq_aligned[flip_idx] <-
    1 - merged$freq_aligned[flip_idx]

  # -----------------------------------------------------------------------
  # MAF filter
  # -----------------------------------------------------------------------

  merged$MAF <- pmin(
    merged$freq_aligned,
    1 - merged$freq_aligned
  )

  merged <- merged[
    merged$MAF >= 0.01,
  ]

  # -----------------------------------------------------------------------
  # Align GWAS and LD
  # -----------------------------------------------------------------------

  snps_gwas <- merged$rsID
  snps_ld <- rownames(ld_matrix_clean)

  final_snps <- intersect(
    snps_ld,
    snps_gwas
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
  # Conditional diagnostic BEFORE explicit outlier removal
  # -----------------------------------------------------------------------

  z_diag <- (
    merged_final$beta_aligned /
      merged_final$se
  )

  lambda_pre <- estimate_s_rss(
    z_diag,
    ld_pd,
    n = ntotal
  )

  cond_pre <- kriging_rss(
    z_diag,
    ld_pd,
    n = ntotal
  )

  cond_pre_df <- as.data.frame(
    cond_pre$conditional_dist
  )

  cond_pre_df$SNP <-
    rownames(ld_pd)

  cond_pre_df$study <- study
  cond_pre_df$stage <- "pre_outlier_removal"

  fwrite(
    cond_pre_df,
    file.path(
      source_dir,
      paste0(
        "Figure11_",
        study,
        "_conditional_pre_outlier_source_data.csv"
      )
    )
  )

  ggsave(
    file.path(
      figdir,
      paste0(
        study,
        "_conditional_pre_outlier.png"
      )
    ),
    cond_pre$plot,
    width = 10,
    height = 5,
    dpi = 300
  )

  cat(
    study,
    "lambda before outlier removal:",
    lambda_pre,
    "\n"
  )

  # -----------------------------------------------------------------------
  # Explicit FinnGen outlier removal
  # -----------------------------------------------------------------------

  if (
    !is.na(outlier) &&
      outlier %in% merged_final$rsID
  ) {

    cat(
      study,
      "removing diagnostic outlier:",
      outlier,
      "\n"
    )

    keep_post <- (
      merged_final$rsID != outlier
    )

    merged_final <- merged_final[
      keep_post,
    ]

    ld_post <- ld_pd[
      keep_post,
      keep_post,
      drop = FALSE
    ]

    ld_pd <- as.matrix(
      nearPD(
        ld_post,
        corr = TRUE
      )$mat
    )

    final_snps_post <-
      merged_final$rsID

    rownames(ld_pd) <- final_snps_post
    colnames(ld_pd) <- final_snps_post
  }

  # -----------------------------------------------------------------------
  # Conditional diagnostic AFTER explicit outlier removal
  # -----------------------------------------------------------------------

  z_post <- (
    merged_final$beta_aligned /
      merged_final$se
  )

  lambda_post <- estimate_s_rss(
    z_post,
    ld_pd,
    n = ntotal
  )

  cond_post <- kriging_rss(
    z_post,
    ld_pd,
    n = ntotal
  )

  cond_post_df <- as.data.frame(
    cond_post$conditional_dist
  )

  cond_post_df$SNP <-
    rownames(ld_pd)

  cond_post_df$study <- study
  cond_post_df$stage <- "post_outlier_removal"

  fwrite(
    cond_post_df,
    file.path(
      source_dir,
      paste0(
        "Figure11_",
        study,
        "_conditional_post_outlier_source_data.csv"
      )
    )
  )

  ggsave(
    file.path(
      figdir,
      paste0(
        study,
        "_conditional_post_outlier.png"
      )
    ),
    cond_post$plot,
    width = 10,
    height = 5,
    dpi = 300
  )

  cat(
    study,
    "lambda after outlier removal:",
    lambda_post,
    "\n"
  )

  # -----------------------------------------------------------------------
  # SuSiE
  # -----------------------------------------------------------------------

  beta_aligned <- merged_final$beta_aligned
  se_aligned <- merged_final$se
  varbeta <- se_aligned^2
  snp <- merged_final$rsID
  pos <- merged_final$pos

  named_beta <- setNames(
    beta_aligned,
    snp
  )

  named_varbeta <- setNames(
    varbeta,
    snp
  )

  dataset1 <- list(
    beta = named_beta,
    varbeta = named_varbeta,
    s = case_fraction,
    type = "cc",
    snp = snp,
    LD = ld_pd,
    position = pos,
    N = ntotal
  )

  coloc:::check_dataset(
    dataset1,
    1
  )

  susie_fit <- runsusie(
    dataset1,
    nref = 503,
    prior_variance = 0.2^2,
    estimate_prior_variance = FALSE,
    coverage = 0.90
  )

  # -----------------------------------------------------------------------
  # Figure 11F/G: PIP source data
  # -----------------------------------------------------------------------

  pip_source <- data.table(
    SNP_index = seq_along(susie_fit$pip),
    SNP = snp,
    PIP = susie_fit$pip
  )

  pip_source[, study := study]

  fwrite(
    pip_source,
    file.path(
      source_dir,
      paste0(
        "Figure11_",
        study,
        "_SuSiE_PIP_source_data.csv"
      )
    )
  )

  png(
    file.path(
      figdir,
      paste0(
        study,
        "_SuSiE_PIP.png"
      )
    ),
    width = 3000,
    height = 1500,
    res = 300
  )

  susie_plot(
    susie_fit,
    y = "PIP"
  )

  lead_pip <- which.max(
    susie_fit$pip
  )

  points(
    x = lead_pip,
    y = susie_fit$pip[lead_pip],
    col = 2,
    pch = 16
  )

  dev.off()

  # -----------------------------------------------------------------------
  # Credible-set summary
  # -----------------------------------------------------------------------

  if (
    !is.null(susie_fit$sets$cs) &&
      length(susie_fit$sets$cs) > 0
  ) {

    cs_list <- lapply(
      seq_along(susie_fit$sets$cs),
      function(i) {

        idx <- susie_fit$sets$cs[[i]]

        data.table(
          study = study,
          cs = i,
          variable = idx,
          SNP = snp[idx],
          PIP = susie_fit$pip[idx]
        )
      }
    )

    cs_result <- rbindlist(
      cs_list,
      use.names = TRUE,
      fill = TRUE
    )

  } else {

    cs_result <- data.table(
      study = character(),
      cs = integer(),
      variable = integer(),
      SNP = character(),
      PIP = numeric()
    )
  }

  fwrite(
    cs_result,
    file.path(
      processed_dir,
      paste0(
        study,
        "_susie_credible_sets.csv"
      )
    )
  )

  cat(
    study,
    "final SNP count:",
    length(snp),
    "\n"
  )

  cat(
    study,
    "lead PIP SNP:",
    snp[which.max(susie_fit$pip)],
    "\n"
  )

  cat(
    study,
    "lead PIP:",
    max(susie_fit$pip),
    "\n"
  )

  invisible(
    list(
      fit = susie_fit,
      pip = pip_source,
      credible_sets = cs_result,
      lambda_pre = lambda_pre,
      lambda_post = lambda_post
    )
  )
}

# =========================================================================
# Run UKB and FinnGen
# =========================================================================

ukb_result <- run_susie(
  study = "UKB",
  gwas_file = cohorts$UKB$gwas,
  bim_file = cohorts$UKB$bim,
  ld_file = cohorts$UKB$ld,
  ncase = cohorts$UKB$ncase,
  ncontrol = cohorts$UKB$ncontrol,
  outlier = cohorts$UKB$outlier
)

finngen_result <- run_susie(
  study = "FinnGen",
  gwas_file = cohorts$FinnGen$gwas,
  bim_file = cohorts$FinnGen$bim,
  ld_file = cohorts$FinnGen$ld,
  ncase = cohorts$FinnGen$ncase,
  ncontrol = cohorts$FinnGen$ncontrol,
  outlier = cohorts$FinnGen$outlier
)

cat("\nAll SuSiE analyses completed.\n")