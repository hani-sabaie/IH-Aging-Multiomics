# ============================================================================
# Reviewer C7 audit:
# Batched nboot=1000 rerun of the historical 18-pathway CellChat LR family.
#
# Usage:
#   Rscript scripts/14C7A_CellChat_Full18_Batch_Audit.R Young 1 25
#
# Each batch:
#   - uses the exact historical condition-specific LR family
#   - uses the full historical data.signaling matrix
#   - runs computeCommunProb(..., nboot=1000, seed.use=1)
#   - validates observed probabilities against the historical object
#   - exports raw permutation P and plus-one empirical P
#
# BH correction is NOT performed here. It will be applied once across all
# condition-specific batches after the complete family has been recovered.
#
# No canonical file is modified.
# ============================================================================

rm(list = ls(all.names = TRUE))

suppressPackageStartupMessages({
  library(CellChat)
  library(data.table)
  library(future)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3L) {
  stop("Usage: <Young|Aged> <start_LR_index> <end_LR_index>")
}

ds <- args[1]
start_idx <- as.integer(args[2])
end_idx <- as.integer(args[3])

if (!(ds %in% c("Young", "Aged"))) {
  stop("Condition must be Young or Aged.")
}

cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd, value = TRUE)

if (length(file_arg) == 1L) {
  script_dir <- dirname(
    normalizePath(sub("^--file=", "", file_arg))
  )
  repo_root <- normalizePath(file.path(script_dir, ".."))
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
  "full18_nboot1000_BH",
  "batches"
)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(hist_file)) {
  stop("Historical CellChat object not found: ", hist_file)
}

cat("Loading historical CellChat object...\n")
x <- readRDS(hist_file)

historical_prob <- x@net[[ds]]$prob
historical_lr <- dimnames(historical_prob)[[3]]
historical_sources <- dimnames(historical_prob)[[1]]

expected_n_lr <- if (ds == "Young") 488L else 468L

if (length(historical_lr) != expected_n_lr) {
  stop(
    "Expected ", expected_n_lr, " historical LR interactions for ",
    ds, "; found ", length(historical_lr), "."
  )
}

if (
  is.na(start_idx) ||
  is.na(end_idx) ||
  start_idx < 1L ||
  end_idx < start_idx ||
  end_idx > expected_n_lr
) {
  stop(
    "Invalid LR index range: ",
    start_idx, "-", end_idx,
    "; valid range is 1-", expected_n_lr, "."
  )
}

# --------------------------------------------------------------------------
# Recover exact historical LRsig family and ordering
# --------------------------------------------------------------------------

lr_full <- x@LR[[ds]]$LRsig

idx <- match(
  historical_lr,
  as.character(lr_full$interaction_name)
)

if (anyNA(idx)) {
  stop("Historical @net LR family does not match @LR$LRsig.")
}

lr_full <- lr_full[idx, , drop = FALSE]
rownames(lr_full) <- historical_lr

batch_indices <- seq.int(start_idx, end_idx)

lr_batch <- lr_full[
  batch_indices,
  ,
  drop = FALSE
]

batch_lr_names <- historical_lr[batch_indices]

if (
  !identical(
    as.character(lr_batch$interaction_name),
    batch_lr_names
  )
) {
  stop("Batch LR ordering mismatch.")
}

# --------------------------------------------------------------------------
# Historical cells / metadata / signaling matrix
# --------------------------------------------------------------------------

cells <- rownames(x@meta)[
  x@meta$datasets == ds
]

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
  stop("Historical data.signaling and metadata are misaligned.")
}

meta$group_cellchat <- factor(
  meta$group_cellchat,
  levels = historical_sources
)

if (anyNA(meta$group_cellchat)) {
  stop("CellChat grouping reconstruction failed.")
}

cat("\n============================================================\n")
cat("CELLCHAT FULL-18 BATCH AUDIT\n")
cat("============================================================\n\n")

cat("Condition        :", ds, "\n")
cat("LR batch         :", start_idx, "-", end_idx, "\n")
cat("LRs in batch     :", nrow(lr_batch), "\n")
cat("Total LR family  :", expected_n_lr, "\n")
cat("Cells            :", length(cells), "\n")
cat("Signaling genes  :", nrow(full_data), "\n\n")

# --------------------------------------------------------------------------
# Reconstruct CellChat object
# --------------------------------------------------------------------------

cc <- createCellChat(
  object = full_data,
  meta = meta,
  group.by = "group_cellchat"
)

cc@DB <- x@DB
cc@data.signaling <- full_data

cc@idents <- factor(
  meta$group_cellchat,
  levels = historical_sources
)

names(cc@idents) <- rownames(meta)

# Sequential avoids the large Windows multisession overhead encountered in
# the monolithic full-family run.
future::plan("sequential")

cat("Running computeCommunProb(nboot=1000)...\n")

tm <- system.time({

  cc <- computeCommunProb(
    cc,
    type = "triMean",
    LR.use = lr_batch,
    raw.use = TRUE,
    population.size = FALSE,
    nboot = 1000,
    seed.use = 1L
  )

})

cat(
  "Elapsed seconds:",
  unname(tm["elapsed"]),
  "\n\n"
)

cc <- filterCommunication(
  cc,
  min.cells = 10
)

# --------------------------------------------------------------------------
# Validate observed probability matrix
# --------------------------------------------------------------------------

new_prob <- cc@net$prob
new_p <- cc@net$pval

new_lr <- dimnames(new_prob)[[3]]

new_idx <- match(
  batch_lr_names,
  new_lr
)

if (anyNA(new_idx)) {
  stop("Some requested LR interactions are absent from batch result.")
}

new_prob <- new_prob[
  ,
  ,
  new_idx,
  drop = FALSE
]

new_p <- new_p[
  ,
  ,
  new_idx,
  drop = FALSE
]

dimnames(new_prob)[[3]] <- batch_lr_names
dimnames(new_p)[[3]] <- batch_lr_names

historical_batch_prob <- historical_prob[
  ,
  ,
  batch_indices,
  drop = FALSE
]

same_dim <- identical(
  dim(new_prob),
  dim(historical_batch_prob)
)

same_source <- identical(
  dimnames(new_prob)[[1]],
  dimnames(historical_batch_prob)[[1]]
)

same_target <- identical(
  dimnames(new_prob)[[2]],
  dimnames(historical_batch_prob)[[2]]
)

same_lr <- identical(
  dimnames(new_prob)[[3]],
  dimnames(historical_batch_prob)[[3]]
)

delta <- abs(
  as.numeric(new_prob) -
    as.numeric(historical_batch_prob)
)

max_delta <- max(delta, na.rm = TRUE)
mean_delta <- mean(delta, na.rm = TRUE)

probability_match <- (
  same_dim &&
  same_source &&
  same_target &&
  same_lr &&
  isTRUE(
    all.equal(
      as.numeric(new_prob),
      as.numeric(historical_batch_prob),
      tolerance = 1e-12,
      check.attributes = FALSE
    )
  )
)

cat("Probability validation:\n")
cat("Same dimensions   :", same_dim, "\n")
cat("Same source order :", same_source, "\n")
cat("Same target order :", same_target, "\n")
cat("Same LR order     :", same_lr, "\n")
cat(
  "Max |delta|      :",
  format(max_delta, scientific = TRUE, digits = 12),
  "\n"
)
cat(
  "Mean |delta|     :",
  format(mean_delta, scientific = TRUE, digits = 12),
  "\n"
)
cat("Match <= 1e-12   :", probability_match, "\n\n")

if (!probability_match) {
  stop(
    "Batch observed probabilities do not reproduce historical values."
  )
}

# --------------------------------------------------------------------------
# Export marginal permutation statistics
# --------------------------------------------------------------------------

B <- 1000L

tests <- as.data.table(
  as.table(new_p)
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
  probability := as.numeric(new_prob)
]

tests[
  ,
  permutation_exceedance_count :=
    as.integer(
      pmin(
        B,
        pmax(
          0,
          round(p_raw_nboot1000 * B)
        )
      )
    )
]

tests[
  ,
  p_plus1 :=
    (
      permutation_exceedance_count + 1
    ) / (
      B + 1
    )
]

tests[
  ,
  `:=`(
    dataset = ds,
    LR_family_index =
      match(
        interaction_name,
        historical_lr
      )
  )
]

setcolorder(
  tests,
  c(
    "dataset",
    "LR_family_index",
    "interaction_name",
    "source",
    "target",
    "probability",
    "p_raw_nboot1000",
    "permutation_exceedance_count",
    "p_plus1"
  )
)

expected_rows <-
  length(batch_indices) *
  length(historical_sources) *
  length(historical_sources)

if (nrow(tests) != expected_rows) {
  stop(
    "Expected ", expected_rows,
    " test rows; recovered ", nrow(tests), "."
  )
}

tag <- sprintf(
  "%s_LR%04d-%04d",
  ds,
  start_idx,
  end_idx
)

fwrite(
  tests,
  file.path(
    outdir,
    paste0(tag, "_tests.tsv")
  ),
  sep = "\t"
)

validation <- data.table(
  dataset = ds,
  LR_start = start_idx,
  LR_end = end_idx,
  n_LR = length(batch_indices),
  n_tests = nrow(tests),
  elapsed_seconds = unname(tm["elapsed"]),
  same_dimensions = same_dim,
  same_source_order = same_source,
  same_target_order = same_target,
  same_LR_order = same_lr,
  max_abs_probability_delta = max_delta,
  mean_abs_probability_delta = mean_delta,
  probability_match_1e12 = probability_match
)

fwrite(
  validation,
  file.path(
    outdir,
    paste0(tag, "_validation.tsv")
  ),
  sep = "\t"
)

cat("Batch written successfully:\n")
print(validation)

cat("\nNo canonical file was modified.\n")

future::plan("sequential")
