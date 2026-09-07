# ============================================================================
# Reviewer C7 audit
# CellChat rerun with increased permutation resolution.
#
# Historical analysis:
#   computeCommunProb(..., type = "triMean")
#   default nboot = 100
#
# Audit rerun:
#   computeCommunProb(..., type = "triMean",
#                     nboot = 1000, seed.use = 1L)
#
# All other analytical settings are preserved.
#
# No canonical file is modified.
# ============================================================================

rm(list = ls(all.names = TRUE))
gc()

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(data.table)
})

set.seed(1234)

# --------------------------------------------------------------------------
# Repository root
# --------------------------------------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (length(file_arg) == 1L) {
  script_path <- sub("^--file=", "", file_arg)
  script_dir <- dirname(normalizePath(script_path))
  repo_root <- normalizePath(file.path(script_dir, ".."))
} else {
  repo_root <- normalizePath(".")
}

# --------------------------------------------------------------------------
# Resolve original Seurat object
# --------------------------------------------------------------------------
input_candidates <- c(
  file.path(
    repo_root,
    "outputs",
    "decont_merged_filt_nodoub_cc_sct_reduc_clust_integ_annot_obj.rds"
  ),
  "F:/Hani's Files/Hernia/outputs/decont_merged_filt_nodoub_cc_sct_reduc_clust_integ_annot_obj.rds"
)

input_hits <- input_candidates[file.exists(input_candidates)]

if (length(input_hits) == 0L) {
  stop(
    "Original annotated Seurat object not found. Checked:\n  ",
    paste(input_candidates, collapse = "\n  ")
  )
}

input_obj_file <- normalizePath(
  input_hits[1],
  winslash = "/",
  mustWork = TRUE
)

historical_cellchat_file <- file.path(
  repo_root,
  "processed_results",
  "12_CellChat",
  "cellchat_merged.rds"
)

if (!file.exists(historical_cellchat_file)) {
  stop(
    "Historical CellChat merged object not found: ",
    historical_cellchat_file
  )
}

audit_dir <- file.path(
  repo_root,
  "processed_results",
  "12_CellChat",
  "multiple_testing_audit",
  "nboot1000"
)

dir.create(
  audit_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("Original Seurat input:\n")
cat(input_obj_file, "\n\n")

cat("Historical CellChat object:\n")
cat(historical_cellchat_file, "\n\n")

cat("Audit output directory:\n")
cat(audit_dir, "\n\n")

cat("CellChat version: ")
cat(as.character(packageVersion("CellChat")), "\n")

cat("Audit nboot: 1000\n")
cat("computeCommunProb seed.use: 1\n\n")

# --------------------------------------------------------------------------
# Historical CellChat object
# --------------------------------------------------------------------------
historical_merged <- readRDS(
  historical_cellchat_file
)

# --------------------------------------------------------------------------
# Verify the exact pre-specified CellChat database
# --------------------------------------------------------------------------
pathways_use <- c(
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

CellChatDB_use <- subsetDB(
  CellChatDB.human,
  search = pathways_use,
  key = "pathway_name"
)

historical_interactions <- sort(
  unique(
    as.character(
      historical_merged@DB$interaction$interaction_name
    )
  )
)

current_interactions <- sort(
  unique(
    as.character(
      CellChatDB_use$interaction$interaction_name
    )
  )
)

db_exact_match <- identical(
  historical_interactions,
  current_interactions
)

cat("============================================================\n")
cat("DATABASE VALIDATION\n")
cat("============================================================\n\n")

cat(
  "Historical DB interactions:",
  length(historical_interactions),
  "\n"
)

cat(
  "Current subset DB interactions:",
  length(current_interactions),
  "\n"
)

cat(
  "Exact interaction-set match:",
  db_exact_match,
  "\n\n"
)

if (!db_exact_match) {

  cat("Historical-only interactions:\n")
  print(
    setdiff(
      historical_interactions,
      current_interactions
    )
  )

  cat("\nCurrent-only interactions:\n")
  print(
    setdiff(
      current_interactions,
      historical_interactions
    )
  )

  stop(
    "Current CellChat DB subset does not exactly match the historical object."
  )
}

# --------------------------------------------------------------------------
# Load original Seurat object
# --------------------------------------------------------------------------
cat("Loading original annotated Seurat object...\n")

obj <- readRDS(
  input_obj_file
)

cat(
  "Cells in original object:",
  ncol(obj),
  "\n\n"
)

# --------------------------------------------------------------------------
# Reproduce historical CellChat grouping exactly
# --------------------------------------------------------------------------
obj$samples <- obj$sample

obj$cc_group <- NA_character_

obj$cc_group[
  obj$skeletal_muscle %in% c(
    "FAP1",
    "FAP2",
    "FAP4"
  )
] <- "Other FAPs"

obj$cc_group[
  obj$skeletal_muscle == "FAP3"
] <- "FAP3"

obj$cc_group[
  obj$skeletal_muscle == "MuSC"
] <- "MuSC"

obj$cc_group[
  obj$skeletal_muscle == "MSM"
] <- "Myogenic"

obj$cc_group[
  obj$skeletal_muscle %in% c(
    "Macrophage",
    "B/T/NK"
  )
] <- "Immune"

obj$cc_group[
  obj$skeletal_muscle %in% c(
    "EC1",
    "EC2",
    "EC3"
  )
] <- "Endothelial"

obj$cc_group[
  obj$skeletal_muscle %in% c(
    "Pericyte",
    "SMC1",
    "SMC2"
  )
] <- "Vascular stromal"

if (anyNA(obj$cc_group)) {
  stop(
    "Some cells were not assigned to a CellChat group: ",
    sum(is.na(obj$cc_group))
  )
}

cat("CellChat groups in full object:\n")
print(
  table(
    obj$condition,
    obj$cc_group
  )
)

# --------------------------------------------------------------------------
# Historical condition split
# --------------------------------------------------------------------------
cat("\nCreating Young and Aged subsets...\n")

obj_young <- subset(
  obj,
  subset = condition == "Young"
)

obj_aged <- subset(
  obj,
  subset = condition == "Aged"
)

rm(obj)
gc()

obj_young$group_cellchat <- obj_young$cc_group
obj_aged$group_cellchat <- obj_aged$cc_group

cat("\nYoung grouping:\n")
print(
  table(obj_young$group_cellchat)
)

cat("\nAged grouping:\n")
print(
  table(obj_aged$group_cellchat)
)

# --------------------------------------------------------------------------
# CellChat function
# Only intentional analytical change: nboot = 1000
# --------------------------------------------------------------------------
options(
  future.globals.maxSize = 8 * 1024^3
)

future::plan(
  "multisession",
  workers = 4
)

make_cellchat_nboot1000 <- function(
    seu,
    dataset_name
) {

  cat("\n============================================================\n")
  cat("RUNNING CELLCHAT:", dataset_name, "\n")
  cat("============================================================\n\n")

  data.input <- GetAssayData(
    seu,
    assay = "SCT",
    slot = "data"
  )

  meta <- seu@meta.data

  cellchat <- createCellChat(
    object = data.input,
    meta = meta,
    group.by = "group_cellchat",
    assay = "SCT"
  )

  # Use the exact pathway subset validated above.
  cellchat@DB <- CellChatDB_use

  cat(dataset_name, "- subsetData\n")
  cellchat <- subsetData(
    cellchat
  )

  cat(dataset_name, "- identifyOverExpressedGenes\n")
  cellchat <- identifyOverExpressedGenes(
    cellchat
  )

  cat(dataset_name, "- identifyOverExpressedInteractions\n")
  cellchat <- identifyOverExpressedInteractions(
    cellchat
  )

  cat(
    dataset_name,
    "- computeCommunProb: nboot = 1000, seed.use = 1\n"
  )

  tm <- system.time({

    cellchat <- computeCommunProb(
      cellchat,
      type = "triMean",
      nboot = 1000,
      seed.use = 1L
    )

  })

  cat(
    dataset_name,
    "- computeCommunProb elapsed seconds:",
    unname(tm["elapsed"]),
    "\n"
  )

  cat(dataset_name, "- filterCommunication\n")

  cellchat <- filterCommunication(
    cellchat,
    min.cells = 10
  )

  cat(dataset_name, "- computeCommunProbPathway\n")

  cellchat <- computeCommunProbPathway(
    cellchat
  )

  cat(dataset_name, "- aggregateNet\n")

  cellchat <- aggregateNet(
    cellchat
  )

  cat(dataset_name, "- netAnalysis_computeCentrality\n")

  cellchat <- netAnalysis_computeCentrality(
    cellchat,
    slot.name = "netP"
  )

  cellchat
}

# --------------------------------------------------------------------------
# Young
# --------------------------------------------------------------------------
cellchat_young <- make_cellchat_nboot1000(
  obj_young,
  "Young"
)

saveRDS(
  cellchat_young,
  file.path(
    audit_dir,
    "cellchat_Young_nboot1000.rds"
  )
)

rm(obj_young)
gc()

# --------------------------------------------------------------------------
# Aged
# --------------------------------------------------------------------------
cellchat_aged <- make_cellchat_nboot1000(
  obj_aged,
  "Aged"
)

saveRDS(
  cellchat_aged,
  file.path(
    audit_dir,
    "cellchat_Aged_nboot1000.rds"
  )
)

rm(obj_aged)
gc()

# --------------------------------------------------------------------------
# Validate observed communication probabilities against historical run.
#
# nboot controls permutation p-value resolution and should not alter the
# observed communication probability estimates.
# --------------------------------------------------------------------------
flatten_probability <- function(
    arr,
    value_name
) {

  lr_names <- dimnames(arr)[[3]]

  if (is.null(lr_names)) {
    stop(
      "Interaction names missing from probability array."
    )
  }

  dt <- as.data.table(
    as.table(arr)
  )

  setnames(
    dt,
    c(
      "source",
      "target",
      "interaction_name",
      value_name
    )
  )

  dt
}

compare_probabilities <- function(
    dataset,
    new_obj
) {

  old_arr <- historical_merged@net[[dataset]]$prob
  new_arr <- new_obj@net$prob

  old_dt <- flatten_probability(
    old_arr,
    "prob_historical"
  )

  new_dt <- flatten_probability(
    new_arr,
    "prob_nboot1000"
  )

  keys <- c(
    "source",
    "target",
    "interaction_name"
  )

  old_only <- fsetdiff(
    old_dt[, ..keys],
    new_dt[, ..keys]
  )

  new_only <- fsetdiff(
    new_dt[, ..keys],
    old_dt[, ..keys]
  )

  merged <- merge(
    old_dt,
    new_dt,
    by = keys,
    all = TRUE,
    sort = FALSE
  )

  shared <- merged[
    is.finite(prob_historical) &
      is.finite(prob_nboot1000)
  ]

  if (nrow(shared) > 0L) {

    shared[
      ,
      abs_delta :=
        abs(
          prob_historical -
            prob_nboot1000
        )
    ]

    max_abs_delta <- max(
      shared$abs_delta,
      na.rm = TRUE
    )

    mean_abs_delta <- mean(
      shared$abs_delta,
      na.rm = TRUE
    )

  } else {

    max_abs_delta <- NA_real_
    mean_abs_delta <- NA_real_
  }

  data.table(
    dataset = dataset,
    historical_probability_rows =
      nrow(old_dt),
    nboot1000_probability_rows =
      nrow(new_dt),
    shared_probability_rows =
      nrow(shared),
    historical_only_keys =
      nrow(old_only),
    nboot1000_only_keys =
      nrow(new_only),
    max_abs_probability_delta =
      max_abs_delta,
    mean_abs_probability_delta =
      mean_abs_delta,
    probability_match_1e12 =
      nrow(old_only) == 0L &&
      nrow(new_only) == 0L &&
      is.finite(max_abs_delta) &&
      max_abs_delta < 1e-12
  )
}

validation_young <- compare_probabilities(
  "Young",
  cellchat_young
)

validation_aged <- compare_probabilities(
  "Aged",
  cellchat_aged
)

validation <- rbindlist(
  list(
    validation_young,
    validation_aged
  )
)

# --------------------------------------------------------------------------
# Merge audit objects
# --------------------------------------------------------------------------
object.list <- list(
  Young = cellchat_young,
  Aged = cellchat_aged
)

cellchat_merged_nboot1000 <- mergeCellChat(
  object.list = object.list,
  add.names = names(object.list)
)

saveRDS(
  cellchat_merged_nboot1000,
  file.path(
    audit_dir,
    "cellchat_merged_nboot1000.rds"
  )
)

# --------------------------------------------------------------------------
# Basic p-value resolution summary
# --------------------------------------------------------------------------
summarize_pvalues <- function(
    dataset,
    obj_condition
) {

  p <- obj_condition@net$pval

  vals <- as.numeric(p)

  data.table(
    dataset = dataset,
    n_tests = length(vals),
    n_p_zero =
      sum(vals == 0, na.rm = TRUE),
    n_p_lt_0_05 =
      sum(vals < 0.05, na.rm = TRUE),
    n_unique_p =
      length(unique(vals)),
    minimum_nonzero_p =
      ifelse(
        any(vals > 0, na.rm = TRUE),
        min(
          vals[vals > 0],
          na.rm = TRUE
        ),
        NA_real_
      )
  )
}

p_summary <- rbindlist(
  list(
    summarize_pvalues(
      "Young",
      cellchat_young
    ),
    summarize_pvalues(
      "Aged",
      cellchat_aged
    )
  )
)

# --------------------------------------------------------------------------
# Write audit summaries
# --------------------------------------------------------------------------
fwrite(
  validation,
  file.path(
    audit_dir,
    "nboot1000_probability_validation.tsv"
  ),
  sep = "\t"
)

fwrite(
  p_summary,
  file.path(
    audit_dir,
    "nboot1000_pvalue_resolution_summary.tsv"
  ),
  sep = "\t"
)

# --------------------------------------------------------------------------
# Console output
# --------------------------------------------------------------------------
cat("\n============================================================\n")
cat("OBSERVED-PROBABILITY VALIDATION\n")
cat("============================================================\n\n")

print(validation)

cat("\n============================================================\n")
cat("NBOOT 1000 P-VALUE RESOLUTION\n")
cat("============================================================\n\n")

print(p_summary)

cat("\nAudit objects written to:\n")
cat(audit_dir, "\n")

cat("\nNo canonical CellChat object or Figure 18 file was modified.\n")

future::plan("sequential")
