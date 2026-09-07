# ============================================================================
# Reviewer C7 audit:
# Effect of multiple-testing correction for SMAD3 LinkPeaks on the downstream
# high-confidence TF-gene-peak network.
#
# This script reproduces the historical downstream filtering from 12A and
# changes ONLY the chromatin-link significance criterion:
#
#   historical: raw LinkPeaks p < 1e-5
#   corrected : Bonferroni p < 0.05 across 200 eligible SMAD3 cis peaks
#
# No canonical files are modified.
# ============================================================================

rm(list = ls(all.names = TRUE))

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (length(file_arg) == 1L) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg)))
  repo_root <- normalizePath(file.path(script_dir, ".."))
} else {
  repo_root <- normalizePath(".")
}

processed_dir <- file.path(
  repo_root,
  "processed_results",
  "10_TF_network"
)

audit_dir <- file.path(
  processed_dir,
  "multiple_testing_audit"
)

dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

direct_file <- file.path(
  processed_dir,
  "TF_gene_peak_DIRECT_motifValidated_fap3.csv"
)

historical_highconf_file <- file.path(
  processed_dir,
  "TF_gene_peak_DIRECT_motifValidated_fap3_highconf.csv"
)

bulk_de_file <- file.path(
  repo_root,
  "processed_results",
  "02_differential_expression",
  "bulk_de_sig_faps.csv"
)

summary_file <- file.path(
  audit_dir,
  "SMAD3_LinkPeaks_multiple_testing_summary.tsv"
)

for (f in c(
  direct_file,
  historical_highconf_file,
  bulk_de_file,
  summary_file
)) {
  if (!file.exists(f)) {
    stop("Required file not found: ", f)
  }
}

# --------------------------------------------------------------------------
# Read inputs
# --------------------------------------------------------------------------

read_historical_csv <- function(path) {

  x <- read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  # Historical write.csv(..., row.names = TRUE) exports may contain an
  # unnamed first column. dplyr cannot operate on data frames with "" names.
  bad_names <- is.na(names(x)) | names(x) == ""

  if (any(bad_names)) {
    names(x)[bad_names] <- paste0(
      "historical_row_id_",
      seq_len(sum(bad_names))
    )
  }

  x
}

direct <- read_historical_csv(
  direct_file
)

historical_highconf <- read_historical_csv(
  historical_highconf_file
)

bulk_de <- read_historical_csv(
  bulk_de_file
)

link_summary <- read.delim(
  summary_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (nrow(link_summary) != 1L) {
  stop("Expected exactly one LinkPeaks summary row.")
}

n_tests <- link_summary$n_eligible_SMAD3_cis_peaks_500kb[1]

if (is.na(n_tests) || n_tests < 1L) {
  stop("Invalid number of reconstructed SMAD3 cis tests.")
}

cat("SMAD3 LinkPeaks test family:", n_tests, "cis peaks\n")

# --------------------------------------------------------------------------
# Reproduce the historical SMAD3 DEG join exactly as in 12A
# --------------------------------------------------------------------------

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
    "Expected exactly one historical SMAD3 DEG row; found ",
    nrow(smad3_de_df)
  )
}

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
      s_cor > 0 & s_reg ==  s_SMAD3 & s_reg != 0 ~ TRUE,
      s_cor < 0 & s_reg == -s_SMAD3 & s_reg != 0 ~ TRUE,
      TRUE ~ FALSE
    ),

    peak_p_bonferroni_200 = pmin(
      peak_pval * n_tests,
      1
    )
  )

# --------------------------------------------------------------------------
# Common filtering criteria from historical 12A
# --------------------------------------------------------------------------

common_pass <- base %>%
  filter(
    consistent_three_way,
    !is.na(logFC_SMAD3),

    abs(Cor) >= 0.25,
    reg_score >= 0.02,
    Frequency >= 0.02,

    abs(avg_log2FC_deg) >= 0.25,
    p_val_adj_deg < 0.01,

    peak_z >= 4,
    peak_score >= 0.05
  )

# Historical raw-P rule
old_pass <- common_pass %>%
  filter(
    peak_pval < 1e-5
  )

# Corrected multiple-testing rule
bonf_pass <- common_pass %>%
  filter(
    peak_p_bonferroni_200 < 0.05
  )

# --------------------------------------------------------------------------
# Validate that the historical reconstruction matches the saved highconf file
# --------------------------------------------------------------------------

key_cols <- c(
  "TF",
  "target_gene",
  "peak_region"
)

old_keys <- old_pass %>%
  select(all_of(key_cols)) %>%
  distinct() %>%
  arrange(across(everything()))

saved_keys <- historical_highconf %>%
  select(all_of(key_cols)) %>%
  distinct() %>%
  arrange(across(everything()))

validation_exact <- identical(
  old_keys,
  saved_keys
)

cat("\n============================================================\n")
cat("HISTORICAL RECONSTRUCTION VALIDATION\n")
cat("============================================================\n\n")

cat("Saved historical high-confidence rows :", nrow(historical_highconf), "\n")
cat("Reconstructed historical rows          :", nrow(old_pass), "\n")
cat("Saved unique TF-gene-peak keys          :", nrow(saved_keys), "\n")
cat("Reconstructed unique keys               :", nrow(old_keys), "\n")
cat("Exact key-set match                     :", validation_exact, "\n")

if (!validation_exact) {

  old_only <- anti_join(
    old_keys,
    saved_keys,
    by = key_cols
  )

  saved_only <- anti_join(
    saved_keys,
    old_keys,
    by = key_cols
  )

  cat("\nReconstructed-only keys:\n")
  print(old_only)

  cat("\nSaved-only keys:\n")
  print(saved_only)

  stop(
    "Historical high-confidence network was not exactly reproduced. ",
    "Do not interpret the corrected comparison yet."
  )
}

# --------------------------------------------------------------------------
# Compare old versus Bonferroni-corrected network
# --------------------------------------------------------------------------

bonf_keys <- bonf_pass %>%
  select(all_of(key_cols)) %>%
  distinct() %>%
  arrange(across(everything()))

added_keys <- anti_join(
  bonf_keys,
  old_keys,
  by = key_cols
)

removed_keys <- anti_join(
  old_keys,
  bonf_keys,
  by = key_cols
)

cat("\n============================================================\n")
cat("OLD RAW-P VS BONFERRONI-CORRECTED NETWORK\n")
cat("============================================================\n\n")

cat("Historical high-confidence rows       :", nrow(old_pass), "\n")
cat("Bonferroni-corrected rows              :", nrow(bonf_pass), "\n")
cat("Historical unique TF-gene-peak keys    :", nrow(old_keys), "\n")
cat("Bonferroni unique TF-gene-peak keys    :", nrow(bonf_keys), "\n")
cat("Added TF-gene-peak keys                :", nrow(added_keys), "\n")
cat("Removed TF-gene-peak keys              :", nrow(removed_keys), "\n")

cat("\nAdded keys:\n")
print(added_keys)

cat("\nRemoved keys:\n")
print(removed_keys)

# --------------------------------------------------------------------------
# Specifically inspect the fifth SMAD3 link that failed raw p < 1e-5 but
# passes Bonferroni across the actual 200-test family.
# --------------------------------------------------------------------------

fifth_peak <- "chr15-67106387-67107602"

fifth_direct <- base %>%
  filter(
    peak_region == fifth_peak
  ) %>%
  arrange(TF)

fifth_common <- common_pass %>%
  filter(
    peak_region == fifth_peak
  ) %>%
  arrange(TF)

fifth_final <- bonf_pass %>%
  filter(
    peak_region == fifth_peak
  ) %>%
  arrange(TF)

cat("\n============================================================\n")
cat("FIFTH SMAD3 LINK AUDIT\n")
cat("============================================================\n\n")

cat("Direct motif/network rows for fifth peak :", nrow(fifth_direct), "\n")
cat("Rows passing non-P common criteria        :", nrow(fifth_common), "\n")
cat("Rows entering corrected final network     :", nrow(fifth_final), "\n")

if (nrow(fifth_direct) > 0L) {
  print(
    fifth_direct %>%
      select(
        TF,
        target_gene,
        Cor,
        reg_score,
        Frequency,
        avg_log2FC_deg,
        p_val_adj_deg,
        peak_region,
        peak_score,
        peak_z,
        peak_pval,
        peak_p_bonferroni_200,
        consistent_three_way
      )
  )
}

# --------------------------------------------------------------------------
# Export audit outputs
# --------------------------------------------------------------------------

write.csv(
  old_pass,
  file.path(
    audit_dir,
    "TF_network_historical_rawP_reconstructed.csv"
  ),
  row.names = FALSE
)

write.csv(
  bonf_pass,
  file.path(
    audit_dir,
    "TF_network_Bonferroni200_corrected.csv"
  ),
  row.names = FALSE
)

write.csv(
  added_keys,
  file.path(
    audit_dir,
    "TF_network_Bonferroni200_added_keys.csv"
  ),
  row.names = FALSE
)

write.csv(
  fifth_direct,
  file.path(
    audit_dir,
    "SMAD3_fifth_peak_downstream_audit.csv"
  ),
  row.names = FALSE
)

cat("\nAudit outputs written to:\n")
cat(audit_dir, "\n")
cat("\nNo canonical file was modified.\n")
