# ============================================================================
# Reviewer C7 audit:
# Full 18-pathway CellChat rerun with improved permutation resolution and
# family-wide multiple-testing correction.
#
# Usage:
#   Rscript 14C7_CellChat_Full18_nboot1000_BH_Audit.R Young
#   Rscript 14C7_CellChat_Full18_nboot1000_BH_Audit.R Aged
#
# Framework:
#   - use the canonical historical merged CellChat object
#   - preserve the exact historical LR family for each condition
#   - preserve the full historical data.signaling matrix
#   - computeCommunProb(type="triMean", nboot=1000, seed.use=1)
#   - convert CellChat b/B permutation P to plus-one empirical P:
#         (b + 1) / (B + 1)
#   - BH across all LR x source x target tests within each condition
#   - rebuild pathway probabilities, aggregate network, and centrality
#
# No canonical CellChat object or Figure 18 file is modified.
# ============================================================================

rm(list = ls(all.names = TRUE))

suppressPackageStartupMessages({
  library(CellChat)
  library(data.table)
  library(future)
})

args <- commandArgs(trailingOnly = TRUE)

if (
  length(args) != 1L ||
  !(args[1] %in% c("Young", "Aged"))
) {
  stop(
    "Provide exactly one condition: Young or Aged."
  )
}

ds <- args[1]

cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd, value = TRUE)

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

hist_file <- file.path(
  repo_root,
  "processed_results",
  "12_CellChat",
  "cellchat_merged.rds"
)

outdir <- file.path(
  repo_root,
  "processed_results",
  "12_CellChat",
  "multiple_testing_audit",
  "full18_nboot1000_BH"
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(hist_file)) {
  stop(
    "Canonical historical CellChat object not found: ",
    hist_file
  )
}

cat("============================================================\n")
cat("FULL 18-PATHWAY CELLCHAT C7 AUDIT\n")
cat("============================================================\n\n")

cat("Condition :", ds, "\n")
cat("Input     :", hist_file, "\n")
cat("Output    :", outdir, "\n")
cat("CellChat  :", as.character(packageVersion("CellChat")), "\n\n")

x <- readRDS(hist_file)

# --------------------------------------------------------------------------
# Historical condition-specific probability family
# --------------------------------------------------------------------------

if (!(ds %in% names(x@net))) {
  stop("Dataset absent from historical @net: ", ds)
}

if (!(ds %in% names(x@LR))) {
  stop("Dataset absent from historical @LR: ", ds)
}

historical_prob <- x@net[[ds]]$prob
historical_pval <- x@net[[ds]]$pval

if (!identical(dim(historical_prob), dim(historical_pval))) {
  stop("Historical probability and P-value arrays differ in dimensions.")
}

historical_lr <- dimnames(historical_prob)[[3]]
historical_sources <- dimnames(historical_prob)[[1]]
historical_targets <- dimnames(historical_prob)[[2]]

if (
  length(historical_sources) != 7L ||
  length(historical_targets) != 7L
) {
  stop(
    "Expected 7 x 7 historical CellChat groups; found ",
    length(historical_sources),
    " x ",
    length(historical_targets),
    "."
  )
}

# --------------------------------------------------------------------------
# Recover exact historical LR family
# --------------------------------------------------------------------------

lr_use <- x@LR[[ds]]$LRsig

if (
  is.null(lr_use) ||
  nrow(lr_use) == 0L
) {
  stop("Historical LRsig table is empty for ", ds)
}

if (!("interaction_name" %in% colnames(lr_use))) {
  stop(
    "Historical LRsig table lacks interaction_name for ",
    ds
  )
}

idx_lr <- match(
  historical_lr,
  as.character(lr_use$interaction_name)
)

if (anyNA(idx_lr)) {
  stop(
    "Some historical @net LR interactions are absent from @LR$LRsig."
  )
}

lr_use <- lr_use[
  idx_lr,
  ,
  drop = FALSE
]

rownames(lr_use) <- historical_lr

if (
  !identical(
    as.character(lr_use$interaction_name),
    historical_lr
  )
) {
  stop(
    "Failed to reproduce historical LR ordering."
  )
}

expected_n_lr <- if (ds == "Young") 488L else 468L
expected_n_tests <- if (ds == "Young") 23912L else 22932L

if (nrow(lr_use) != expected_n_lr) {
  stop(
    "Expected ",
    expected_n_lr,
    " historical LR interactions for ",
    ds,
    "; found ",
    nrow(lr_use),
    "."
  )
}

n_tests <- nrow(lr_use) *
  length(historical_sources) *
  length(historical_targets)

if (n_tests != expected_n_tests) {
  stop(
    "Expected ",
    expected_n_tests,
    " CellChat tests for ",
    ds,
    "; reconstructed ",
    n_tests,
    "."
  )
}

# --------------------------------------------------------------------------
# Confirm all LR interactions belong to the 18 prespecified pathways
# --------------------------------------------------------------------------

db <- x@DB$interaction

db_idx <- match(
  historical_lr,
  as.character(db$interaction_name)
)

if (anyNA(db_idx)) {
  stop(
    "Some historical LR interactions could not be matched to the DB."
  )
}

lr_pathways <- as.character(
  db$pathway_name[db_idx]
)

expected_pathways <- c(
  "TGFb",
  "COLLAGEN",
  "FN1",
  "LAMININ",
  "THBS",
  "TENASCIN",
  "VTN",
  "HSPG",
  "CDH",
  "CDH1",
  "CDH5",
  "PCDH",
  "ICAM",
  "VCAM",
  "JAM",
  "MMP",
  "SPP1",
  "PERIOSTIN"
)

if (
  !setequal(
    unique(lr_pathways),
    expected_pathways
  )
) {
  stop(
    "Historical LR family does not contain exactly the expected 18 pathways."
  )
}

cat("Historical LR interactions :", nrow(lr_use), "\n")
cat("Source groups              :", length(historical_sources), "\n")
cat("Target groups              :", length(historical_targets), "\n")
cat("Formal test family         :", n_tests, "\n")
cat("Prespecified pathways      :", length(unique(lr_pathways)), "\n\n")

# --------------------------------------------------------------------------
# Reconstruct condition directly from historical merged CellChat resources
# --------------------------------------------------------------------------

if (!("datasets" %in% colnames(x@meta))) {
  stop("Merged CellChat metadata lacks datasets column.")
}

if (!("group_cellchat" %in% colnames(x@meta))) {
  stop("Merged CellChat metadata lacks group_cellchat column.")
}

cells <- rownames(x@meta)[
  x@meta$datasets == ds
]

if (length(cells) == 0L) {
  stop("No cells recovered for ", ds)
}

meta <- x@meta[
  cells,
  ,
  drop = FALSE
]

full_data <- x@data.signaling[
  ,
  cells,
  drop = FALSE
]

if (
  !identical(
    colnames(full_data),
    rownames(meta)
  )
) {
  stop(
    "Historical signaling matrix and metadata are not aligned."
  )
}

meta$group_cellchat <- factor(
  meta$group_cellchat,
  levels = historical_sources
)

if (anyNA(meta$group_cellchat)) {
  stop(
    "Some cells could not be assigned to historical CellChat groups."
  )
}

group_sizes <- table(meta$group_cellchat)

if (any(group_sizes <= 10L)) {
  stop(
    "At least one historical CellChat group has <=10 cells."
  )
}

cat("Cells:", length(cells), "\n")
cat("Historical signaling genes:", nrow(full_data), "\n")
cat("Group sizes:\n")
print(group_sizes)
cat("\n")

# --------------------------------------------------------------------------
# Lightweight exact-probability reconstruction
# --------------------------------------------------------------------------

cc <- createCellChat(
  object = full_data,
  meta = meta,
  group.by = "group_cellchat"
)

cc@DB <- x@DB

# Use the exact historical signaling matrix.
cc@data.signaling <- full_data

cc@idents <- factor(
  meta$group_cellchat,
  levels = historical_sources
)

names(cc@idents) <- rownames(meta)

options(
  future.globals.maxSize = 8 * 1024^3
)

future::plan(
  "multisession",
  workers = 4
)

cat("Running computeCommunProb with nboot=1000...\n")

tm <- system.time({

  cc <- computeCommunProb(
    cc,
    type = "triMean",
    LR.use = lr_use,
    raw.use = TRUE,
    population.size = FALSE,
    nboot = 1000,
    seed.use = 1L
  )

})

future::plan("sequential")

cat(
  "computeCommunProb elapsed seconds:",
  unname(tm["elapsed"]),
  "\n\n"
)

cc <- filterCommunication(
  cc,
  min.cells = 10
)

# --------------------------------------------------------------------------
# Validate that observed communication probabilities are unchanged
# --------------------------------------------------------------------------

new_prob <- cc@net$prob
new_p_raw <- cc@net$pval

if (
  !identical(
    dimnames(new_prob),
    dimnames(historical_prob)
  )
) {
  stop(
    "New probability-array ordering does not exactly match historical array."
  )
}

if (
  !identical(
    dimnames(new_p_raw),
    dimnames(historical_pval)
  )
) {
  stop(
    "New P-value-array ordering does not exactly match historical array."
  )
}

prob_delta <- abs(
  as.numeric(new_prob) -
    as.numeric(historical_prob)
)

max_prob_delta <- max(
  prob_delta,
  na.rm = TRUE
)

mean_prob_delta <- mean(
  prob_delta,
  na.rm = TRUE
)

probability_match <- isTRUE(
  all.equal(
    as.numeric(new_prob),
    as.numeric(historical_prob),
    tolerance = 1e-12,
    check.attributes = FALSE
  )
)

cat("Probability validation:\n")
cat(
  "Max |delta|   :",
  format(max_prob_delta, scientific = TRUE, digits = 12),
  "\n"
)
cat(
  "Mean |delta|  :",
  format(mean_prob_delta, scientific = TRUE, digits = 12),
  "\n"
)
cat("Match <=1e-12:", probability_match, "\n\n")

if (!probability_match) {
  stop(
    "Observed probabilities were not reproduced to 1e-12; ",
    "do not interpret corrected P-values."
  )
}

# --------------------------------------------------------------------------
# Plus-one empirical P and BH across the complete condition-specific family
# --------------------------------------------------------------------------

if (anyNA(new_p_raw)) {
  stop(
    "New CellChat P-value array contains NA values."
  )
}

B <- 1000L

exceedance_count <- round(
  new_p_raw * B
)

exceedance_count <- pmin(
  B,
  pmax(
    0,
    exceedance_count
  )
)

p_plus1 <- (
  exceedance_count + 1
) / (
  B + 1
)

p_plus1_BH <- array(
  p.adjust(
    as.numeric(p_plus1),
    method = "BH"
  ),
  dim = dim(p_plus1),
  dimnames = dimnames(p_plus1)
)

# Historical nominal mask versus corrected family-wide mask.
historical_sig <- historical_pval < 0.05
corrected_sig <- p_plus1_BH < 0.05

n_historical_sig <- sum(
  historical_sig,
  na.rm = TRUE
)

n_raw_nboot1000_sig <- sum(
  new_p_raw < 0.05,
  na.rm = TRUE
)

n_corrected_sig <- sum(
  corrected_sig,
  na.rm = TRUE
)

n_sig_overlap <- sum(
  historical_sig &
    corrected_sig,
  na.rm = TRUE
)

n_removed <- sum(
  historical_sig &
    !corrected_sig,
  na.rm = TRUE
)

n_added <- sum(
  !historical_sig &
    corrected_sig,
  na.rm = TRUE
)

# --------------------------------------------------------------------------
# Long test table
# --------------------------------------------------------------------------

tests <- as.data.table(
  as.table(new_p_raw)
)

setnames(
  tests,
  c(
    "source",
    "target",
    "interaction_name",
    "p_raw_nboot1000"
  )
)

tests[
  ,
  probability :=
    as.numeric(new_prob)
]

tests[
  ,
  pathway :=
    lr_pathways[
      match(
        interaction_name,
        historical_lr
      )
    ]
]

tests[
  ,
  permutation_exceedance_count :=
    as.integer(
      as.numeric(exceedance_count)
    )
]

tests[
  ,
  p_plus1 :=
    as.numeric(p_plus1)
]

tests[
  ,
  p_plus1_BH_condition :=
    as.numeric(p_plus1_BH)
]

tests[
  ,
  historical_p_nboot100 :=
    as.numeric(historical_pval)
]

tests[
  ,
  historical_nominal_sig :=
    historical_p_nboot100 < 0.05
]

tests[
  ,
  corrected_BH_sig :=
    p_plus1_BH_condition < 0.05
]

tests[
  ,
  dataset := ds
]

setcolorder(
  tests,
  c(
    "dataset",
    "pathway",
    "interaction_name",
    "source",
    "target",
    "probability",
    "historical_p_nboot100",
    "historical_nominal_sig",
    "p_raw_nboot1000",
    "permutation_exceedance_count",
    "p_plus1",
    "p_plus1_BH_condition",
    "corrected_BH_sig"
  )
)

if (nrow(tests) != n_tests) {
  stop(
    "Long test table has ",
    nrow(tests),
    " rows; expected ",
    n_tests,
    "."
  )
}

# --------------------------------------------------------------------------
# Rebuild CellChat network using corrected P values
# --------------------------------------------------------------------------

cc@net$pval <- p_plus1_BH

cc <- computeCommunProbPathway(
  cc
)

cc <- aggregateNet(
  cc
)

cc <- netAnalysis_computeCentrality(
  cc,
  slot.name = "netP"
)

# --------------------------------------------------------------------------
# Compact comparison summary
# --------------------------------------------------------------------------

historical_pathways <- x@netP[[ds]]$pathways
corrected_pathways <- cc@netP$pathways

historical_count_sum <- sum(
  x@net[[ds]]$count,
  na.rm = TRUE
)

corrected_count_sum <- sum(
  cc@net$count,
  na.rm = TRUE
)

historical_weight_sum <- sum(
  x@net[[ds]]$weight,
  na.rm = TRUE
)

corrected_weight_sum <- sum(
  cc@net$weight,
  na.rm = TRUE
)

summary_dt <- data.table(
  dataset = ds,
  n_cells = length(cells),
  n_LR = nrow(lr_use),
  n_source_groups = length(historical_sources),
  n_target_groups = length(historical_targets),
  n_tests = n_tests,
  n_pathways_in_family = length(unique(lr_pathways)),
  probability_match_1e12 = probability_match,
  max_abs_probability_delta = max_prob_delta,
  mean_abs_probability_delta = mean_prob_delta,
  n_raw_p_zero_nboot1000 =
    sum(new_p_raw == 0, na.rm = TRUE),
  n_raw_p_lt_0_05_nboot1000 =
    n_raw_nboot1000_sig,
  n_historical_nominal_p_lt_0_05 =
    n_historical_sig,
  n_plus1_BH_lt_0_05 =
    n_corrected_sig,
  n_historical_corrected_overlap =
    n_sig_overlap,
  n_historical_removed_by_correction =
    n_removed,
  n_newly_added_vs_historical =
    n_added,
  n_historical_netP_pathways =
    length(historical_pathways),
  n_corrected_netP_pathways =
    length(corrected_pathways),
  historical_total_network_count =
    historical_count_sum,
  corrected_total_network_count =
    corrected_count_sum,
  historical_total_network_weight =
    historical_weight_sum,
  corrected_total_network_weight =
    corrected_weight_sum
)

# --------------------------------------------------------------------------
# Write audit resources only
# --------------------------------------------------------------------------

fwrite(
  tests,
  file.path(
    outdir,
    paste0(
      ds,
      "_full18_nboot1000_tests.tsv"
    )
  ),
  sep = "\t"
)

fwrite(
  summary_dt,
  file.path(
    outdir,
    paste0(
      ds,
      "_full18_nboot1000_BH_summary.tsv"
    )
  ),
  sep = "\t"
)

fwrite(
  data.table(
    interaction_name = historical_lr,
    pathway = lr_pathways
  ),
  file.path(
    outdir,
    paste0(
      ds,
      "_historical_LR_family.tsv"
    )
  ),
  sep = "\t"
)

saveRDS(
  cc,
  file.path(
    outdir,
    paste0(
      "cellchat_",
      ds,
      "_full18_nboot1000_BH.rds"
    )
  )
)

cat("\n============================================================\n")
cat("AUDIT SUMMARY\n")
cat("============================================================\n\n")

print(summary_dt)

cat("\nHistorical pathways:\n")
print(historical_pathways)

cat("\nCorrected pathways:\n")
print(corrected_pathways)

cat("\nAudit files written to:\n")
cat(outdir, "\n")

cat(
  "\nNo canonical CellChat object or Figure 18 file was modified.\n"
)
