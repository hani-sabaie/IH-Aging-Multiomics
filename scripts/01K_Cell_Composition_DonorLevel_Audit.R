# ============================================================================
# Donor-level exact permutation audit of cell-type composition
#
# Purpose:
#   Re-evaluate Young-vs-Aged cell-type composition using biological donors
#   (5 Young, 5 Aged) as the independent units.
#
# Primary statistic:
#   mean donor proportion in Aged - mean donor proportion in Young
#
# Inference:
#   exact two-sided permutation test over all choose(10, 5) = 252 possible
#   allocations of five donors to each condition.
#
# Multiple testing:
#   BH correction across the 14 prespecified cell-type tests.
#
# This script writes audit outputs only and does NOT overwrite canonical results.
# ============================================================================

rm(list = ls(all.names = TRUE))
gc()

suppressPackageStartupMessages({
  library(Seurat)
})

set.seed(1234)

repo_root <- normalizePath(".")

obj_file <- Sys.getenv(
  "CELLPROP_SEURAT_RDS",
  unset = "F:/Hani's Files/Hernia/outputs/decont_merged_filt_nodoub_cc_sct_reduc_clust_integ_annot_obj.rds"
)

audit_dir <- file.path(
  repo_root,
  "processed_results",
  "03_cell_composition",
  "donor_level_audit"
)

dir.create(
  audit_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("===== DONOR-LEVEL CELL-COMPOSITION AUDIT =====\n")
cat("Object:", obj_file, "\n")
cat("Exists:", file.exists(obj_file), "\n\n")

if (!file.exists(obj_file)) {
  stop("Required Seurat object not found: ", obj_file)
}

obj <- readRDS(obj_file)
md <- obj@meta.data

required <- c(
  "sample",
  "condition",
  "skeletal_muscle"
)

missing_cols <- setdiff(required, colnames(md))

if (length(missing_cols) > 0) {
  stop(
    "Missing metadata columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

md <- md[, required, drop = FALSE]

md$sample <- as.character(md$sample)
md$condition <- as.character(md$condition)
md$skeletal_muscle <- as.character(md$skeletal_muscle)

# --------------------------------------------------------------------------
# Validate donor structure
# --------------------------------------------------------------------------

donors <- sort(unique(md$sample))

donor_condition <- sapply(
  donors,
  function(d) {
    z <- unique(md$condition[md$sample == d])

    if (length(z) != 1L) {
      stop(
        "Donor ",
        d,
        " maps to ",
        length(z),
        " conditions."
      )
    }

    z
  }
)

names(donor_condition) <- donors

cat("Donors:", length(donors), "\n")
print(table(donor_condition))

if (length(donors) != 10L) {
  stop(
    "Expected 10 donors; found ",
    length(donors)
  )
}

if (
  sum(donor_condition == "Aged") != 5L ||
  sum(donor_condition == "Young") != 5L
) {
  stop("Expected exactly 5 Aged and 5 Young donors.")
}

clusters <- sort(unique(md$skeletal_muscle))

cat("\nCell types:", length(clusters), "\n")
print(clusters)

if (length(clusters) != 14L) {
  stop(
    "Expected 14 cell types; found ",
    length(clusters)
  )
}

# --------------------------------------------------------------------------
# Donor x cell-type counts and proportions
# --------------------------------------------------------------------------

count_mat <- table(
  factor(md$sample, levels = donors),
  factor(md$skeletal_muscle, levels = clusters)
)

prop_mat <- sweep(
  count_mat,
  1,
  rowSums(count_mat),
  "/"
)

donor_prop_df <- data.frame(
  sample = donors,
  condition = unname(donor_condition[donors]),
  as.data.frame.matrix(prop_mat),
  check.names = FALSE
)

write.csv(
  donor_prop_df,
  file.path(
    audit_dir,
    "donor_celltype_proportions.csv"
  ),
  row.names = FALSE
)

donor_count_df <- data.frame(
  sample = donors,
  condition = unname(donor_condition[donors]),
  as.data.frame.matrix(count_mat),
  check.names = FALSE
)

write.csv(
  donor_count_df,
  file.path(
    audit_dir,
    "donor_celltype_counts.csv"
  ),
  row.names = FALSE
)

# --------------------------------------------------------------------------
# Exact permutation universe
# --------------------------------------------------------------------------

aged_obs <- donors[
  donor_condition[donors] == "Aged"
]

young_obs <- donors[
  donor_condition[donors] == "Young"
]

aged_combinations <- combn(
  donors,
  5,
  simplify = FALSE
)

if (length(aged_combinations) != 252L) {
  stop(
    "Expected 252 exact allocations; found ",
    length(aged_combinations)
  )
}

cat(
  "\nExact permutation allocations:",
  length(aged_combinations),
  "\n"
)

# --------------------------------------------------------------------------
# Per-cell-type exact tests
# --------------------------------------------------------------------------

result_list <- vector(
  "list",
  length(clusters)
)

names(result_list) <- clusters

for (ct in clusters) {

  x <- prop_mat[, ct]
  names(x) <- donors

  aged_mean <- mean(
    x[aged_obs]
  )

  young_mean <- mean(
    x[young_obs]
  )

  obs_diff <- aged_mean - young_mean

  obs_log2_ratio <- if (
    aged_mean > 0 &&
    young_mean > 0
  ) {
    log2(
      aged_mean /
        young_mean
    )
  } else {
    NA_real_
  }

  perm_diff <- vapply(
    aged_combinations,
    function(aged_set) {

      young_set <- setdiff(
        donors,
        aged_set
      )

      mean(x[aged_set]) -
        mean(x[young_set])
    },
    numeric(1)
  )

  # Exact two-sided permutation P.
  # All 252 allocations are enumerated, including the observed allocation.
  tol <- 1e-15

  n_extreme <- sum(
    abs(perm_diff) >=
      abs(obs_diff) - tol
  )

  exact_p <- n_extreme /
    length(perm_diff)

  result_list[[ct]] <- data.frame(
    clusters = ct,
    Aged_mean_donor_proportion = aged_mean,
    Young_mean_donor_proportion = young_mean,
    mean_difference_Aged_minus_Young = obs_diff,
    donor_mean_log2_ratio_Aged_vs_Young = obs_log2_ratio,
    n_Aged_donors = length(aged_obs),
    n_Young_donors = length(young_obs),
    n_exact_allocations = length(perm_diff),
    n_as_or_more_extreme = n_extreme,
    exact_permutation_p = exact_p,
    stringsAsFactors = FALSE
  )
}

res <- do.call(
  rbind,
  result_list
)

rownames(res) <- NULL

res$BH_FDR_14 <- p.adjust(
  res$exact_permutation_p,
  method = "BH"
)

res$significant_BH_0_05 <-
  res$BH_FDR_14 < 0.05

# --------------------------------------------------------------------------
# Compare with historical cell-level scProportionTest results
# --------------------------------------------------------------------------

historical_file <- file.path(
  repo_root,
  "processed_results",
  "03_cell_composition",
  "cell_proportion_test_results.csv"
)

if (file.exists(historical_file)) {

  hist <- read.csv(
    historical_file,
    check.names = FALSE
  )

  hist_keep <- hist[, c(
    "clusters",
    "Aged",
    "Young",
    "obs_log2FD",
    "pval",
    "FDR"
  )]

  names(hist_keep) <- c(
    "clusters",
    "historical_pooled_Aged",
    "historical_pooled_Young",
    "historical_pooled_log2FD",
    "historical_celllevel_p",
    "historical_celllevel_FDR"
  )

  res <- merge(
    res,
    hist_keep,
    by = "clusters",
    all.x = TRUE,
    sort = FALSE
  )
}

res <- res[
  match(
    clusters,
    res$clusters
  ),
]

write.csv(
  res,
  file.path(
    audit_dir,
    "donor_level_exact_permutation_results.csv"
  ),
  row.names = FALSE
)

# --------------------------------------------------------------------------
# Console summary
# --------------------------------------------------------------------------

cat("\n===== DONOR-LEVEL EXACT RESULTS =====\n\n")

print(
  res[, c(
    "clusters",
    "Aged_mean_donor_proportion",
    "Young_mean_donor_proportion",
    "mean_difference_Aged_minus_Young",
    "donor_mean_log2_ratio_Aged_vs_Young",
    "exact_permutation_p",
    "BH_FDR_14",
    "significant_BH_0_05"
  )],
  row.names = FALSE,
  digits = 6
)

cat(
  "\nSignificant after BH:",
  sum(
    res$significant_BH_0_05
  ),
  "/",
  nrow(res),
  "\n"
)

cat(
  "Minimum attainable two-sided exact P observed:",
  min(
    res$exact_permutation_p
  ),
  "\n"
)

cat("\n===== HISTORICAL VS DONOR-LEVEL SIGNIFICANCE =====\n")

if (
  "historical_celllevel_FDR" %in%
    names(res)
) {

  comparison <- data.frame(
    clusters = res$clusters,
    historical_celllevel_sig =
      res$historical_celllevel_FDR < 0.05,
    donor_level_BH_sig =
      res$BH_FDR_14 < 0.05
  )

  comparison$changed <-
    comparison$historical_celllevel_sig !=
    comparison$donor_level_BH_sig

  print(
    comparison,
    row.names = FALSE
  )

  cat(
    "\nSignificance-status changes:",
    sum(comparison$changed),
    "\n"
  )
}

cat(
  "\nPASS: donor-level exact permutation audit completed.\n"
)
cat(
  "No canonical result files were overwritten.\n"
)
