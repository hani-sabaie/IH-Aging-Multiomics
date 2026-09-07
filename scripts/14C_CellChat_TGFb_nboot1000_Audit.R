# ============================================================================
# Reviewer C7 audit
# Targeted TGFb CellChat rerun with increased permutation resolution.
#
# Historical:
#   18 pre-specified pathways
#   computeCommunProb(type = "triMean", default nboot = 100)
#
# Audit:
#   TGFb only
#   computeCommunProb(type = "triMean", nboot = 1000, seed.use = 1L)
#
# The observed TGFb communication probabilities are compared against the
# historical 18-pathway run before the new p-values are interpreted.
#
# No canonical files are modified.
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
# Input paths
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
    "Original Seurat object not found. Checked:\n  ",
    paste(input_candidates, collapse = "\n  ")
  )
}

input_obj_file <- normalizePath(
  input_hits[1],
  winslash = "/",
  mustWork = TRUE
)

historical_file <- file.path(
  repo_root,
  "processed_results",
  "12_CellChat",
  "cellchat_merged.rds"
)

if (!file.exists(historical_file)) {
  stop(
    "Historical CellChat object not found: ",
    historical_file
  )
}

audit_dir <- file.path(
  repo_root,
  "processed_results",
  "12_CellChat",
  "multiple_testing_audit",
  "TGFb_nboot1000"
)

dir.create(
  audit_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("Original Seurat input:\n", input_obj_file, "\n\n", sep = "")
cat("Historical CellChat:\n", historical_file, "\n\n", sep = "")
cat("CellChat version: ", as.character(packageVersion("CellChat")), "\n", sep = "")
cat("Targeted pathway: TGFb\n")
cat("Audit nboot: 1000\n")
cat("seed.use: 1\n\n")

# --------------------------------------------------------------------------
# Historical object and TGFb DB validation
# --------------------------------------------------------------------------
historical <- readRDS(historical_file)

tgfb_db <- subsetDB(
  CellChatDB.human,
  search = "TGFb",
  key = "pathway_name"
)

current_tgfb_lr <- sort(
  unique(
    as.character(
      tgfb_db$interaction$interaction_name
    )
  )
)

historical_tgfb_lr <- sort(
  unique(
    as.character(
      historical@DB$interaction$interaction_name[
        historical@DB$interaction$pathway_name == "TGFb"
      ]
    )
  )
)

db_match <- identical(
  current_tgfb_lr,
  historical_tgfb_lr
)

cat("============================================================\n")
cat("TGFb DATABASE VALIDATION\n")
cat("============================================================\n\n")

cat("Historical TGFb LR pairs :", length(historical_tgfb_lr), "\n")
cat("Current TGFb LR pairs    :", length(current_tgfb_lr), "\n")
cat("Exact LR-set match       :", db_match, "\n\n")

if (!db_match) {

  cat("Historical-only:\n")
  print(setdiff(historical_tgfb_lr, current_tgfb_lr))

  cat("\nCurrent-only:\n")
  print(setdiff(current_tgfb_lr, historical_tgfb_lr))

  stop("TGFb DB mismatch. Audit stopped before rerun.")
}

print(current_tgfb_lr)

# --------------------------------------------------------------------------
# Load Seurat object and reconstruct historical grouping
# --------------------------------------------------------------------------
cat("\nLoading original Seurat object...\n")

obj <- readRDS(input_obj_file)

obj$samples <- obj$sample
obj$cc_group <- NA_character_

obj$cc_group[
  obj$skeletal_muscle %in% c("FAP1", "FAP2", "FAP4")
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
  obj$skeletal_muscle %in% c("Macrophage", "B/T/NK")
] <- "Immune"

obj$cc_group[
  obj$skeletal_muscle %in% c("EC1", "EC2", "EC3")
] <- "Endothelial"

obj$cc_group[
  obj$skeletal_muscle %in% c("Pericyte", "SMC1", "SMC2")
] <- "Vascular stromal"

if (anyNA(obj$cc_group)) {
  stop(
    "Unassigned CellChat cells: ",
    sum(is.na(obj$cc_group))
  )
}

cat("\nHistorical grouping reconstructed:\n")
print(table(obj$condition, obj$cc_group))

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

# --------------------------------------------------------------------------
# Parallel configuration
# --------------------------------------------------------------------------
options(
  future.globals.maxSize = 8 * 1024^3
)

future::plan(
  "multisession",
  workers = 4
)

# --------------------------------------------------------------------------
# Targeted TGFb CellChat rerun
# --------------------------------------------------------------------------
make_tgfb_cellchat <- function(seu, dataset_name) {

  cat("\n============================================================\n")
  cat("RUNNING TARGETED TGFb:", dataset_name, "\n")
  cat("============================================================\n\n")

  data.input <- tryCatch(
    GetAssayData(
      seu,
      assay = "SCT",
      layer = "data"
    ),
    error = function(e) {
      GetAssayData(
        seu,
        assay = "SCT",
        slot = "data"
      )
    }
  )

  meta <- seu@meta.data

  cellchat <- createCellChat(
    object = data.input,
    meta = meta,
    group.by = "group_cellchat",
    assay = "SCT"
  )

  cellchat@DB <- tgfb_db

  cat(dataset_name, ": subsetData\n")
  cellchat <- subsetData(cellchat)

  cat(dataset_name, ": identifyOverExpressedGenes\n")
  cellchat <- identifyOverExpressedGenes(cellchat)

  cat(dataset_name, ": identifyOverExpressedInteractions\n")
  cellchat <- identifyOverExpressedInteractions(cellchat)

  cat(dataset_name, ": computeCommunProb nboot=1000\n")

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
    ": computeCommunProb elapsed =",
    unname(tm["elapsed"]),
    "seconds\n"
  )

  cellchat <- filterCommunication(
    cellchat,
    min.cells = 10
  )

  cellchat <- computeCommunProbPathway(
    cellchat
  )

  cellchat <- aggregateNet(
    cellchat
  )

  cellchat
}

# --------------------------------------------------------------------------
# Run Young
# --------------------------------------------------------------------------
young_new <- make_tgfb_cellchat(
  obj_young,
  "Young"
)

saveRDS(
  young_new,
  file.path(
    audit_dir,
    "cellchat_Young_TGFb_nboot1000.rds"
  )
)

rm(obj_young)
gc()

# --------------------------------------------------------------------------
# Run Aged
# --------------------------------------------------------------------------
aged_new <- make_tgfb_cellchat(
  obj_aged,
  "Aged"
)

saveRDS(
  aged_new,
  file.path(
    audit_dir,
    "cellchat_Aged_TGFb_nboot1000.rds"
  )
)

rm(obj_aged)
gc()

# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------
flatten_array <- function(arr, value_name) {

  if (is.null(dimnames(arr)[[3]])) {
    stop("LR interaction names missing from CellChat array.")
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

# --------------------------------------------------------------------------
# Validate observed TGFb probabilities against historical analysis
# --------------------------------------------------------------------------
validate_condition <- function(
    ds,
    new_obj
) {

  old_prob <- historical@net[[ds]]$prob

  old_dt <- flatten_array(
    old_prob,
    "prob_historical"
  )[
    interaction_name %in% historical_tgfb_lr
  ]

  new_dt <- flatten_array(
    new_obj@net$prob,
    "prob_nboot1000"
  )

  keys <- c(
    "source",
    "target",
    "interaction_name"
  )

  old_keys <- unique(
    old_dt[, ..keys]
  )

  new_keys <- unique(
    new_dt[, ..keys]
  )

  old_only <- fsetdiff(
    old_keys,
    new_keys
  )

  new_only <- fsetdiff(
    new_keys,
    old_keys
  )

  m <- merge(
    old_dt,
    new_dt,
    by = keys,
    all = TRUE,
    sort = FALSE
  )

  shared <- m[
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

    max_delta <- max(
      shared$abs_delta,
      na.rm = TRUE
    )

    mean_delta <- mean(
      shared$abs_delta,
      na.rm = TRUE
    )

  } else {

    max_delta <- NA_real_
    mean_delta <- NA_real_
  }

  data.table(
    dataset = ds,

    historical_TGFb_rows =
      nrow(old_dt),

    nboot1000_TGFb_rows =
      nrow(new_dt),

    shared_rows =
      nrow(shared),

    historical_only_keys =
      nrow(old_only),

    nboot1000_only_keys =
      nrow(new_only),

    max_abs_probability_delta =
      max_delta,

    mean_abs_probability_delta =
      mean_delta,

    probability_match_1e12 =
      nrow(old_only) == 0L &&
      nrow(new_only) == 0L &&
      is.finite(max_delta) &&
      max_delta < 1e-12
  )
}

validation <- rbindlist(
  list(
    validate_condition(
      "Young",
      young_new
    ),
    validate_condition(
      "Aged",
      aged_new
    )
  )
)

# --------------------------------------------------------------------------
# Stop interpretation if targeted rerun changes observed probabilities
# --------------------------------------------------------------------------
cat("\n============================================================\n")
cat("OBSERVED TGFb PROBABILITY VALIDATION\n")
cat("============================================================\n\n")

print(validation)

fwrite(
  validation,
  file.path(
    audit_dir,
    "TGFb_nboot1000_probability_validation.tsv"
  ),
  sep = "\t"
)

if (!all(validation$probability_match_1e12)) {

  cat(
    "\nTargeted TGFb rerun did not reproduce the historical ",
    "observed probabilities exactly.\n"
  )

  cat(
    "The nboot=1000 p-values will NOT be interpreted as an ",
    "isolated permutation-resolution change.\n"
  )

  stop(
    "Probability validation failed."
  )
}

# --------------------------------------------------------------------------
# Multiple-testing audit using nboot=1000 p-values
# --------------------------------------------------------------------------
audit_condition <- function(
    ds,
    new_obj
) {

  prob_dt <- flatten_array(
    new_obj@net$prob,
    "probability"
  )

  p_dt <- flatten_array(
    new_obj@net$pval,
    "p_raw"
  )

  dt <- merge(
    prob_dt,
    p_dt,
    by = c(
      "source",
      "target",
      "interaction_name"
    ),
    all = TRUE,
    sort = FALSE
  )

  dt[
    ,
    dataset := ds
  ]

  dt[
    ,
    p_BH_TGFb :=
      p.adjust(
        p_raw,
        method = "BH"
      )
  ]

  # CellChat p = b / B, with B = 1000 here.
  dt[
    ,
    permutation_exceedance_count :=
      round(p_raw * 1000)
  ]

  # Plus-one empirical p-value sensitivity:
  # (b + 1) / (B + 1)
  dt[
    ,
    p_plus1 :=
      (
        permutation_exceedance_count + 1
      ) / 1001
  ]

  dt[
    ,
    p_plus1_BH_TGFb :=
      p.adjust(
        p_plus1,
        method = "BH"
      )
  ]

  dt
}

young_audit <- audit_condition(
  "Young",
  young_new
)

aged_audit <- audit_condition(
  "Aged",
  aged_new
)

audit_dt <- rbindlist(
  list(
    young_audit,
    aged_audit
  )
)

# Pooled-condition sensitivity
audit_dt[
  ,
  p_BH_TGFb_pooled :=
    p.adjust(
      p_raw,
      method = "BH"
    )
]

audit_dt[
  ,
  p_plus1_BH_TGFb_pooled :=
    p.adjust(
      p_plus1,
      method = "BH"
    )
]

# --------------------------------------------------------------------------
# Overall summary
# --------------------------------------------------------------------------
summary_dt <- audit_dt[
  ,
  .(
    n_tests = .N,

    n_LR_pairs =
      uniqueN(interaction_name),

    n_raw_p_zero =
      sum(p_raw == 0),

    n_raw_p_lt_0_05 =
      sum(p_raw < 0.05),

    n_BH_lt_0_05 =
      sum(p_BH_TGFb < 0.05),

    n_plus1_BH_lt_0_05 =
      sum(p_plus1_BH_TGFb < 0.05),

    n_pooled_BH_lt_0_05 =
      sum(p_BH_TGFb_pooled < 0.05),

    n_pooled_plus1_BH_lt_0_05 =
      sum(p_plus1_BH_TGFb_pooled < 0.05),

    n_unique_raw_p =
      uniqueN(p_raw),

    min_nonzero_raw_p =
      if (
        any(p_raw > 0)
      ) {
        min(p_raw[p_raw > 0])
      } else {
        NA_real_
      }
  ),
  by = dataset
]

# --------------------------------------------------------------------------
# Figure 18F hypothesis subset
# --------------------------------------------------------------------------
bubble_sources <- c(
  "FAP3",
  "Other FAPs"
)

bubble_targets <- c(
  "MuSC",
  "Myogenic",
  "Immune",
  "Endothelial",
  "Vascular stromal"
)

fig18f <- audit_dt[
  source %in% bubble_sources &
    target %in% bubble_targets
]

fig18f_summary <- fig18f[
  ,
  .(
    n_hypotheses = .N,

    n_probability_gt_0 =
      sum(probability > 0),

    n_raw_p_lt_0_05 =
      sum(p_raw < 0.05),

    n_BH_lt_0_05 =
      sum(p_BH_TGFb < 0.05),

    n_plus1_BH_lt_0_05 =
      sum(p_plus1_BH_TGFb < 0.05)
  ),
  by = dataset
]

# --------------------------------------------------------------------------
# TGFb source-target network sensitivity
# --------------------------------------------------------------------------
edge_dt <- audit_dt[
  ,
  .(
    raw_weight =
      sum(
        fifelse(
          p_raw < 0.05,
          probability,
          0
        )
      ),

    BH_weight =
      sum(
        fifelse(
          p_BH_TGFb < 0.05,
          probability,
          0
        )
      ),

    plus1_BH_weight =
      sum(
        fifelse(
          p_plus1_BH_TGFb < 0.05,
          probability,
          0
        )
      )
  ),
  by = .(
    dataset,
    source,
    target
  )
]

edge_summary <- edge_dt[
  ,
  .(
    n_raw_nonzero_edges =
      sum(raw_weight > 0),

    n_BH_nonzero_edges =
      sum(BH_weight > 0),

    n_plus1_BH_nonzero_edges =
      sum(plus1_BH_weight > 0),

    total_raw_weight =
      sum(raw_weight),

    total_BH_weight =
      sum(BH_weight),

    total_plus1_BH_weight =
      sum(plus1_BH_weight)
  ),
  by = dataset
]

# --------------------------------------------------------------------------
# Save targeted merged audit object
# --------------------------------------------------------------------------
targeted_list <- list(
  Young = young_new,
  Aged = aged_new
)

targeted_merged <- mergeCellChat(
  targeted_list,
  add.names = names(targeted_list)
)

saveRDS(
  targeted_merged,
  file.path(
    audit_dir,
    "cellchat_merged_TGFb_nboot1000.rds"
  )
)

# --------------------------------------------------------------------------
# Export
# --------------------------------------------------------------------------
fwrite(
  audit_dt,
  file.path(
    audit_dir,
    "TGFb_nboot1000_all_tests.tsv"
  ),
  sep = "\t"
)

fwrite(
  summary_dt,
  file.path(
    audit_dir,
    "TGFb_nboot1000_multiple_testing_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  fig18f,
  file.path(
    audit_dir,
    "Figure18F_TGFb_nboot1000_audit.tsv"
  ),
  sep = "\t"
)

fwrite(
  fig18f_summary,
  file.path(
    audit_dir,
    "Figure18F_TGFb_nboot1000_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  edge_dt,
  file.path(
    audit_dir,
    "TGFb_nboot1000_source_target_edges.tsv"
  ),
  sep = "\t"
)

fwrite(
  edge_summary,
  file.path(
    audit_dir,
    "TGFb_nboot1000_network_summary.tsv"
  ),
  sep = "\t"
)

# --------------------------------------------------------------------------
# Console summary
# --------------------------------------------------------------------------
cat("\n============================================================\n")
cat("TGFb NBOOT=1000 MULTIPLE-TESTING SUMMARY\n")
cat("============================================================\n\n")

print(summary_dt)

cat("\n============================================================\n")
cat("FIGURE 18F SUBSET\n")
cat("============================================================\n\n")

print(fig18f_summary)

cat("\n============================================================\n")
cat("TGFb NETWORK SUMMARY\n")
cat("============================================================\n\n")

print(edge_summary)

cat("\n============================================================\n")
cat("SMALLEST P-VALUES\n")
cat("============================================================\n\n")

for (ds in c("Young", "Aged")) {

  z <- audit_dt[
    dataset == ds
  ]

  cat(ds, "\n")

  print(
    head(
      z[
        order(p_raw),
        .(
          source,
          target,
          interaction_name,
          probability,
          p_raw,
          p_BH_TGFb,
          p_plus1,
          p_plus1_BH_TGFb
        )
      ],
      20
    )
  )

  cat("\n")
}

cat("\nAudit outputs written to:\n")
cat(audit_dir, "\n")

cat("\nNo canonical CellChat object or Figure 18 file was modified.\n")

future::plan("sequential")
