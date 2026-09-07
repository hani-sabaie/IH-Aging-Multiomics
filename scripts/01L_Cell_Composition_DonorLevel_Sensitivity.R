# ============================================================================
# Donor-level cell-composition sensitivity analysis
#
# Compares three donor-level inferential frameworks:
#   1) exact permutation test of mean proportion difference
#   2) exact Wilcoxon rank-sum test
#   3) Welch t-test
#
# BH correction is applied separately across the same 14 cell-type hypotheses.
# Audit only: canonical files are not modified.
# ============================================================================

rm(list = ls(all.names = TRUE))

audit_dir <- file.path(
  "processed_results",
  "03_cell_composition",
  "donor_level_audit"
)

prop_file <- file.path(
  audit_dir,
  "donor_celltype_proportions.csv"
)

primary_file <- file.path(
  audit_dir,
  "donor_level_exact_permutation_results.csv"
)

if (!file.exists(prop_file)) {
  stop("Missing donor proportions: ", prop_file)
}

if (!file.exists(primary_file)) {
  stop("Missing exact permutation results: ", primary_file)
}

dat <- read.csv(
  prop_file,
  check.names = FALSE
)

primary <- read.csv(
  primary_file,
  check.names = FALSE
)

required <- c("sample", "condition")

if (!all(required %in% names(dat))) {
  stop("Required donor metadata columns missing.")
}

celltypes <- setdiff(
  names(dat),
  required
)

if (length(celltypes) != 14L) {
  stop(
    "Expected 14 cell types; found ",
    length(celltypes)
  )
}

if (
  sum(dat$condition == "Aged") != 5L ||
  sum(dat$condition == "Young") != 5L
) {
  stop("Expected 5 Aged and 5 Young donors.")
}

# --------------------------------------------------------------------------
# Exact two-sided permutation test of mean difference
# --------------------------------------------------------------------------

donors <- dat$sample
aged_sets <- combn(
  donors,
  5,
  simplify = FALSE
)

exact_perm <- function(x) {

  names(x) <- donors

  aged_obs <- donors[
    dat$condition == "Aged"
  ]

  young_obs <- donors[
    dat$condition == "Young"
  ]

  obs <- mean(x[aged_obs]) -
    mean(x[young_obs])

  null <- vapply(
    aged_sets,
    function(a) {
      y <- setdiff(donors, a)

      mean(x[a]) -
        mean(x[y])
    },
    numeric(1)
  )

  sum(
    abs(null) >= abs(obs) - 1e-15
  ) / length(null)
}

# --------------------------------------------------------------------------
# Run all three methods
# --------------------------------------------------------------------------

out <- lapply(
  celltypes,
  function(ct) {

    x <- dat[[ct]]

    aged <- x[
      dat$condition == "Aged"
    ]

    young <- x[
      dat$condition == "Young"
    ]

    perm_p <- exact_perm(x)

    wilcox_p <- suppressWarnings(
      wilcox.test(
        aged,
        young,
        alternative = "two.sided",
        exact = TRUE,
        correct = FALSE
      )$p.value
    )

    welch_p <- t.test(
      aged,
      young,
      alternative = "two.sided",
      var.equal = FALSE
    )$p.value

    data.frame(
      clusters = ct,
      Aged_mean = mean(aged),
      Young_mean = mean(young),
      mean_difference = mean(aged) - mean(young),
      exact_permutation_p = perm_p,
      exact_wilcoxon_p = wilcox_p,
      welch_t_p = welch_p,
      stringsAsFactors = FALSE
    )
  }
)

out <- do.call(
  rbind,
  out
)

out$BH_exact_permutation <- p.adjust(
  out$exact_permutation_p,
  method = "BH"
)

out$BH_exact_wilcoxon <- p.adjust(
  out$exact_wilcoxon_p,
  method = "BH"
)

out$BH_welch_t <- p.adjust(
  out$welch_t_p,
  method = "BH"
)

out$sig_exact_permutation <-
  out$BH_exact_permutation < 0.05

out$sig_exact_wilcoxon <-
  out$BH_exact_wilcoxon < 0.05

out$sig_welch_t <-
  out$BH_welch_t < 0.05

# --------------------------------------------------------------------------
# Confirm primary exact-permutation results reproduce 01K
# --------------------------------------------------------------------------

cmp <- merge(
  out[, c(
    "clusters",
    "exact_permutation_p",
    "BH_exact_permutation"
  )],
  primary[, c(
    "clusters",
    "exact_permutation_p",
    "BH_FDR_14"
  )],
  by = "clusters",
  suffixes = c("_01L", "_01K")
)

delta_p <- max(
  abs(
    cmp$exact_permutation_p_01L -
      cmp$exact_permutation_p_01K
  )
)

delta_fdr <- max(
  abs(
    cmp$BH_exact_permutation -
      cmp$BH_FDR_14
  )
)

if (
  delta_p > 1e-12 ||
  delta_fdr > 1e-12
) {
  stop("01L does not reproduce 01K exact permutation results.")
}

# --------------------------------------------------------------------------
# Save and report
# --------------------------------------------------------------------------

write.csv(
  out,
  file.path(
    audit_dir,
    "donor_level_method_sensitivity.csv"
  ),
  row.names = FALSE
)

cat("===== DONOR-LEVEL METHOD SENSITIVITY =====\n\n")

print(
  out[, c(
    "clusters",
    "exact_permutation_p",
    "BH_exact_permutation",
    "exact_wilcoxon_p",
    "BH_exact_wilcoxon",
    "welch_t_p",
    "BH_welch_t",
    "sig_exact_permutation",
    "sig_exact_wilcoxon",
    "sig_welch_t"
  )],
  row.names = FALSE,
  digits = 6
)

cat("\n===== SIGNIFICANCE COUNTS =====\n")
cat(
  "Exact permutation:",
  sum(out$sig_exact_permutation),
  "/ 14\n"
)
cat(
  "Exact Wilcoxon   :",
  sum(out$sig_exact_wilcoxon),
  "/ 14\n"
)
cat(
  "Welch t-test     :",
  sum(out$sig_welch_t),
  "/ 14\n"
)

cat("\n===== FAP-SPECIFIC RESULTS =====\n")

print(
  subset(
    out,
    clusters %in% c(
      "FAP1",
      "FAP2",
      "FAP3",
      "FAP4"
    )
  )[, c(
    "clusters",
    "mean_difference",
    "BH_exact_permutation",
    "BH_exact_wilcoxon",
    "BH_welch_t",
    "sig_exact_permutation",
    "sig_exact_wilcoxon",
    "sig_welch_t"
  )],
  row.names = FALSE,
  digits = 6
)

cat("\n===== METHOD AGREEMENT =====\n")

out$n_methods_significant <-
  rowSums(
    out[, c(
      "sig_exact_permutation",
      "sig_exact_wilcoxon",
      "sig_welch_t"
    )]
  )

print(
  out[, c(
    "clusters",
    "n_methods_significant"
  )],
  row.names = FALSE
)

cat(
  "\nExact-permutation reproduction max delta P:",
  delta_p,
  "\n"
)

cat(
  "Exact-permutation reproduction max delta FDR:",
  delta_fdr,
  "\n"
)

cat(
  "\nPASS: donor-level sensitivity analysis completed.\n"
)
