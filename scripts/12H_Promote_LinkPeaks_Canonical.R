# ============================================================================
# Production update of canonical SMAD3 LinkPeaks downstream resources
#
# Purpose:
#   Materialize the reviewer-C7 multiple-testing correction in the canonical
#   TF-network tables without rerunning the computationally expensive upstream
#   TF-network / motif pipeline.
#
# The complete production workflow remains scripts/12A_TF_Network_Analysis_StC.R.
#
# Validated framework:
#   - targeted gene: SMAD3
#   - eligible cis LinkPeaks family: 200 peaks
#   - correction: Bonferroni across those 200 tests
#   - corrected chromatin criterion: adjusted P < 0.05
#
# This script:
#   1) updates the existing 22-row canonical direct table with multiplicity
#      metadata;
#   2) reconstructs the final integrated network using the current canonical
#      pseudobulk SMAD3 result;
#   3) validates the exact 5-key result before overwriting canonical files.
# ============================================================================

rm(list = ls(all.names = TRUE))

suppressPackageStartupMessages({
  library(dplyr)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (length(file_arg) == 1L) {
  script_dir <- dirname(
    normalizePath(sub("^--file=", "", file_arg))
  )
  repo_root <- normalizePath(
    file.path(script_dir, "..")
  )
} else {
  repo_root <- normalizePath(".")
}

processed_dir <- file.path(
  repo_root,
  "processed_results",
  "10_TF_network"
)

direct_file <- file.path(
  processed_dir,
  "TF_gene_peak_DIRECT_motifValidated_fap3.csv"
)

highconf_file <- file.path(
  processed_dir,
  "TF_gene_peak_DIRECT_motifValidated_fap3_highconf.csv"
)

bulk_de_file <- file.path(
  repo_root,
  "processed_results",
  "02_differential_expression",
  "bulk_de_sig_faps.csv"
)

for (f in c(
  direct_file,
  highconf_file,
  bulk_de_file
)) {
  if (!file.exists(f)) {
    stop("Required canonical file not found: ", f)
  }
}

# --------------------------------------------------------------------------
# Read canonical CSV while removing the legacy write.csv row-index column.
# --------------------------------------------------------------------------

read_canonical_csv <- function(path) {

  x <- read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  blank_names <- is.na(names(x)) | names(x) == ""

  if (any(blank_names)) {
    x <- x[
      ,
      !blank_names,
      drop = FALSE
    ]
  }

  x
}

direct <- read_canonical_csv(
  direct_file
)

existing_highconf <- read_canonical_csv(
  highconf_file
)

bulk_de <- read_canonical_csv(
  bulk_de_file
)

# --------------------------------------------------------------------------
# Validate canonical direct input
# --------------------------------------------------------------------------

required_direct_cols <- c(
  "TF",
  "target_gene",
  "reg_score",
  "Cor",
  "Frequency",
  "avg_log2FC_deg",
  "p_val_adj_deg",
  "peak_region",
  "peak_score",
  "peak_z",
  "peak_pval"
)

missing_direct_cols <- setdiff(
  required_direct_cols,
  names(direct)
)

if (length(missing_direct_cols) > 0L) {
  stop(
    "Canonical direct table is missing columns: ",
    paste(missing_direct_cols, collapse = ", ")
  )
}

if (nrow(direct) != 22L) {
  stop(
    "Expected 22 canonical direct TF-gene-peak rows; found ",
    nrow(direct),
    "."
  )
}

if (
  !identical(
    unique(direct$target_gene),
    "SMAD3"
  )
) {
  stop(
    "Canonical direct table is not restricted to SMAD3 as expected."
  )
}

# Make the operation idempotent if the script is rerun.
direct <- direct %>%
  select(
    -any_of(
      c(
        "n_linkpeaks_tests",
        "peak_pval_bonferroni",
        "peak_p_bonferroni_200"
      )
    )
  )

# --------------------------------------------------------------------------
# Add production LinkPeaks multiple-testing metadata
# --------------------------------------------------------------------------

n_tests <- 200L

direct <- direct %>%
  mutate(
    n_linkpeaks_tests = n_tests,
    peak_pval_bonferroni = pmin(
      peak_pval * n_tests,
      1
    )
  )

if (
  any(is.na(direct$peak_pval_bonferroni))
) {
  stop(
    "Missing Bonferroni-adjusted P value in canonical direct table."
  )
}

# --------------------------------------------------------------------------
# Current canonical SMAD3 pseudobulk result
# --------------------------------------------------------------------------

required_de_cols <- c(
  "gene",
  "avg_log2FC",
  "p_val",
  "p_val_adj"
)

missing_de_cols <- setdiff(
  required_de_cols,
  names(bulk_de)
)

if (length(missing_de_cols) > 0L) {
  stop(
    "Canonical FAP DEG table is missing columns: ",
    paste(missing_de_cols, collapse = ", ")
  )
}

smad3_de_df <- bulk_de %>%
  filter(gene == "SMAD3") %>%
  transmute(
    target_gene = gene,
    logFC_SMAD3 = avg_log2FC,
    p_val_SMAD3 = p_val,
    p_adj_SMAD3 = p_val_adj
  )

if (nrow(smad3_de_df) != 1L) {
  stop(
    "Expected exactly one canonical SMAD3 FAP DEG row; found ",
    nrow(smad3_de_df),
    "."
  )
}

# Numerical safeguards from the independently validated pseudobulk rerun.
if (
  abs(
    smad3_de_df$logFC_SMAD3 -
      0.993597836739163
  ) > 1e-10
) {
  stop(
    "Canonical SMAD3 logFC differs from the validated pseudobulk result."
  )
}

if (
  abs(
    smad3_de_df$p_adj_SMAD3 -
      0.044459358726055
  ) > 1e-10
) {
  stop(
    "Canonical SMAD3 adjusted P differs from the validated pseudobulk result."
  )
}

# --------------------------------------------------------------------------
# Reconstruct the exact downstream integration used by production 12A
# --------------------------------------------------------------------------

base <- direct %>%
  left_join(
    smad3_de_df,
    by = "target_gene"
  ) %>%
  mutate(
    s_reg = sign(avg_log2FC_deg),
    s_cor = sign(Cor),
    s_SMAD3 = sign(logFC_SMAD3),

    consistent_three_way = case_when(
      s_cor > 0 &
        s_reg == s_SMAD3 &
        s_reg != 0 ~ TRUE,

      s_cor < 0 &
        s_reg == -s_SMAD3 &
        s_reg != 0 ~ TRUE,

      TRUE ~ FALSE
    )
  )

corrected_highconf <- base %>%
  filter(
    consistent_three_way,
    !is.na(logFC_SMAD3),

    # TFNet
    abs(Cor) >= 0.25,
    reg_score >= 0.02,
    Frequency >= 0.02,

    # Regulon activity
    abs(avg_log2FC_deg) >= 0.25,
    p_val_adj_deg < 0.01,

    # Chromatin link
    peak_z >= 4,
    peak_pval_bonferroni < 0.05,
    peak_score >= 0.05
  )

# --------------------------------------------------------------------------
# Validate exact corrected result before modifying canonical files
# --------------------------------------------------------------------------

key_cols <- c(
  "TF",
  "target_gene",
  "peak_region"
)

make_keys <- function(x) {
  sort(
    paste(
      x$TF,
      x$target_gene,
      x$peak_region,
      sep = "|"
    )
  )
}

expected_corrected_keys <- sort(
  c(
    "ETS1|SMAD3|chr15-67198215-67198971",
    "ETS2|SMAD3|chr15-67109227-67110381",
    "ETS2|SMAD3|chr15-67198215-67198971",
    "FOS|SMAD3|chr15-67109227-67110381",
    "FOS|SMAD3|chr15-67106387-67107602"
  )
)

expected_historical_keys <- sort(
  c(
    "ETS1|SMAD3|chr15-67198215-67198971",
    "ETS2|SMAD3|chr15-67109227-67110381",
    "ETS2|SMAD3|chr15-67198215-67198971",
    "FOS|SMAD3|chr15-67109227-67110381"
  )
)

observed_corrected_keys <- make_keys(
  corrected_highconf
)

if (nrow(corrected_highconf) != 5L) {
  stop(
    "Expected 5 corrected high-confidence rows; found ",
    nrow(corrected_highconf),
    "."
  )
}

if (
  !identical(
    observed_corrected_keys,
    expected_corrected_keys
  )
) {
  stop(
    "Corrected high-confidence key set differs from the ",
    "independently validated C7 result."
  )
}

# Existing canonical highconf may be either the historical 4-row resource
# (first promotion) or the corrected 5-row resource (idempotent rerun).
existing_keys <- make_keys(
  existing_highconf
)

if (
  !identical(
    existing_keys,
    expected_historical_keys
  ) &&
  !identical(
    existing_keys,
    expected_corrected_keys
  )
) {
  stop(
    "Existing canonical high-confidence table has an unexpected key set."
  )
}

# Principal reported peak safeguard.
reported_peak <- "chr15-67109227-67110381"

reported_row <- direct %>%
  filter(
    peak_region == reported_peak
  ) %>%
  distinct(
    peak_region,
    peak_pval,
    n_linkpeaks_tests,
    peak_pval_bonferroni
  )

if (nrow(reported_row) != 1L) {
  stop(
    "Expected exactly one distinct LinkPeaks result for the reported peak."
  )
}

expected_reported_bonf <- 3.96983310662e-05

if (
  abs(
    reported_row$peak_pval_bonferroni -
      expected_reported_bonf
  ) > 1e-10
) {
  stop(
    "Reported SMAD3 Bonferroni P differs from the validated value."
  )
}

# Added fifth-link safeguard.
added_key <- "FOS|SMAD3|chr15-67106387-67107602"

if (
  !added_key %in% observed_corrected_keys
) {
  stop(
    "Expected corrected FOS-SMAD3 fifth link is absent."
  )
}

# --------------------------------------------------------------------------
# All checks passed: promote canonical resources
# --------------------------------------------------------------------------

write.csv(
  direct,
  direct_file,
  row.names = TRUE
)

write.csv(
  corrected_highconf,
  highconf_file,
  row.names = TRUE
)

cat("Canonical LinkPeaks production update complete.\n")
cat("Direct TF-gene-peak rows       :", nrow(direct), "\n")
cat("Corrected high-confidence rows :", nrow(corrected_highconf), "\n")
cat("LinkPeaks test family          :", n_tests, "\n")
cat(
  "Reported peak Bonferroni P     :",
  format(
    reported_row$peak_pval_bonferroni,
    scientific = TRUE,
    digits = 12
  ),
  "\n"
)
cat(
  "Canonical SMAD3 logFC          :",
  format(
    smad3_de_df$logFC_SMAD3,
    digits = 12
  ),
  "\n"
)
cat(
  "Canonical SMAD3 FDR            :",
  format(
    smad3_de_df$p_adj_SMAD3,
    scientific = TRUE,
    digits = 12
  ),
  "\n"
)
cat("Validated corrected keys:\n")
cat(
  paste0("  ", observed_corrected_keys),
  sep = "\n"
)
cat("\n")
