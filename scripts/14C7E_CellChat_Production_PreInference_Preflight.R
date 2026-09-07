# ============================================================================
# Reviewer C7 production preflight
#
# Validate that the production CellChat preprocessing stage, starting from
# the Seurat object, reproduces the exact historical signaling matrix and
# condition-specific LRsig families used in the validated nboot=1000 audit.
#
# No computeCommunProb() call is made.
# No canonical result is modified.
# ============================================================================

rm(list = ls(all.names = TRUE))
gc()

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(data.table)
})

# --------------------------------------------------------------------------
# Repository paths
# --------------------------------------------------------------------------

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (length(file_arg) == 1L) {
  script_dir <- dirname(
    normalizePath(
      sub("^--file=", "", file_arg)
    )
  )
  repo_root <- normalizePath(
    file.path(script_dir, "..")
  )
} else {
  repo_root <- normalizePath(".")
}

input_obj_file <- Sys.getenv(
  "CELLCHAT_INPUT_RDS",
  unset = file.path(
    repo_root,
    "outputs",
    "decont_merged_filt_nodoub_cc_sct_reduc_clust_integ_annot_obj.rds"
  )
)

historical_file <- file.path(
  repo_root,
  "processed_results",
  "12_CellChat",
  "cellchat_merged.rds"
)

db_order_file <- Sys.getenv(
  "CELLCHAT_DB_ORDER_FILE",
  unset = file.path(
    repo_root,
    "processed_results",
    "12_CellChat",
    "CellChatDB_human_18pathway_interaction_order.tsv"
  )
)

outdir <- file.path(
  repo_root,
  "processed_results",
  "12_CellChat",
  "multiple_testing_audit",
  "production_preflight"
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

for (f in c(input_obj_file, historical_file, db_order_file)) {
  if (!file.exists(f)) {
    stop("Required input not found: ", f)
  }
}

db_order <- fread(
  db_order_file
)

required_db_order_cols <- c(
  "db_order",
  "interaction_name",
  "pathway_name"
)

missing_db_order_cols <- setdiff(
  required_db_order_cols,
  names(db_order)
)

if (length(missing_db_order_cols) > 0L) {
  stop(
    "DB-order resource is missing columns: ",
    paste(missing_db_order_cols, collapse = ", ")
  )
}

if (
  nrow(db_order) != 518L ||
  uniqueN(db_order$interaction_name) != 518L ||
  uniqueN(db_order$pathway_name) != 18L ||
  !identical(
    as.integer(db_order$db_order),
    seq_len(518L)
  )
) {
  stop(
    "Standalone CellChat DB-order resource failed structural validation."
  )
}

cat("Loading Seurat production input...\n")
obj <- readRDS(input_obj_file)

cat("Loading historical CellChat object...\n")
historical <- readRDS(historical_file)

# --------------------------------------------------------------------------
# Exact production grouping logic
# --------------------------------------------------------------------------

obj$samples <- obj$sample

obj$cc_group <- NA_character_

obj$cc_group[
  obj$skeletal_muscle %in%
    c("FAP1", "FAP2", "FAP4")
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
  obj$skeletal_muscle %in%
    c("Macrophage", "B/T/NK")
] <- "Immune"

obj$cc_group[
  obj$skeletal_muscle %in%
    c("EC1", "EC2", "EC3")
] <- "Endothelial"

obj$cc_group[
  obj$skeletal_muscle %in%
    c("Pericyte", "SMC1", "SMC2")
] <- "Vascular stromal"

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

obj_young$group_cellchat <-
  obj_young$cc_group

obj_aged$group_cellchat <-
  obj_aged$cc_group

group_order <- c(
  "Endothelial",
  "FAP3",
  "Immune",
  "MuSC",
  "Myogenic",
  "Other FAPs",
  "Vascular stromal"
)

pathways.use <- c(
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

CellChatDB.use <- subsetDB(
  CellChatDB.human,
  search = pathways.use,
  key = "pathway_name"
)

# --------------------------------------------------------------------------
# Validate database subset against standalone reproducibility resource
# --------------------------------------------------------------------------

fresh_db_lr <- as.character(
  CellChatDB.use$interaction$interaction_name
)

resource_db_lr <- as.character(
  db_order$interaction_name
)

resource_db_pathway <- as.character(
  db_order$pathway_name
)

same_db_set <- setequal(
  fresh_db_lr,
  resource_db_lr
)

same_db_order_before <- identical(
  fresh_db_lr,
  resource_db_lr
)

cat(
  "\nDB interaction set match to standalone resource:",
  same_db_set,
  "\n"
)

cat(
  "DB interaction order match before reordering:",
  same_db_order_before,
  "\n"
)

if (!same_db_set) {
  stop(
    "Installed CellChat DB interaction content differs from ",
    "the standalone reproducibility resource."
  )
}

db_idx <- match(
  resource_db_lr,
  fresh_db_lr
)

if (anyNA(db_idx)) {
  stop(
    "Failed to map installed CellChat DB to standalone interaction order."
  )
}

CellChatDB.use$interaction <- CellChatDB.use$interaction[
  db_idx,
  ,
  drop = FALSE
]

same_db_order_after <- identical(
  as.character(
    CellChatDB.use$interaction$interaction_name
  ),
  resource_db_lr
)

same_db_pathway_after <- identical(
  as.character(
    CellChatDB.use$interaction$pathway_name
  ),
  resource_db_pathway
)

cat(
  "DB interaction order match after reordering:",
  same_db_order_after,
  "\n"
)

cat(
  "DB pathway annotations match resource:",
  same_db_pathway_after,
  "\n"
)

if (
  !same_db_order_after ||
  !same_db_pathway_after
) {
  stop(
    "Standalone DB-order resource did not reproduce the validated DB state."
  )
}

# Historical object is used only as an independent validation reference.
historical_db_lr <- as.character(
  historical@DB$interaction$interaction_name
)

historical_db_pathway <- as.character(
  historical@DB$interaction$pathway_name
)

historical_db_matches_resource <- (
  identical(
    historical_db_lr,
    resource_db_lr
  ) &&
  identical(
    historical_db_pathway,
    resource_db_pathway
  )
)

cat(
  "Historical validation DB matches standalone resource:",
  historical_db_matches_resource,
  "\n"
)

if (!historical_db_matches_resource) {
  stop(
    "Historical validation object does not match the standalone DB resource."
  )
}

# --------------------------------------------------------------------------
# Condition validation
# --------------------------------------------------------------------------

prepare_and_validate <- function(
    seu,
    ds
) {

  cat("\n============================================================\n")
  cat("PRODUCTION PREFLIGHT:", ds, "\n")
  cat("============================================================\n\n")

  data.input <- GetAssayData(
    seu,
    assay = "SCT",
    slot = "data"
  )

  meta <- seu@meta.data

  meta$group_cellchat <- factor(
    meta$group_cellchat,
    levels = group_order
  )

  if (anyNA(meta$group_cellchat)) {
    stop(
      ds,
      ": production group assignment contains NA."
    )
  }

  cc <- createCellChat(
    object = data.input,
    meta = meta,
    group.by = "group_cellchat",
    assay = "SCT"
  )

  cc@DB <- CellChatDB.use

  cc <- subsetData(cc)

  cc <- identifyOverExpressedGenes(cc)

  cc <- identifyOverExpressedInteractions(cc)

  # ------------------------------------------------------------------------
  # Exact fresh LR family
  # ------------------------------------------------------------------------

  fresh_lr <- as.character(
    cc@LR$LRsig$interaction_name
  )

  hist_lr <- dimnames(
    historical@net[[ds]]$prob
  )[[3]]

  if (is.null(hist_lr)) {
    stop(
      ds,
      ": historical LR dimnames are missing."
    )
  }

  same_lr_count <-
    length(fresh_lr) ==
      length(hist_lr)

  same_lr_set <-
    setequal(
      fresh_lr,
      hist_lr
    )

  same_lr_order <-
    identical(
      fresh_lr,
      hist_lr
    )

  # ------------------------------------------------------------------------
  # Cell identity/order
  # ------------------------------------------------------------------------

  hist_cells <- rownames(
    historical@meta
  )[
    historical@meta$datasets ==
      ds
  ]

  fresh_cells <- colnames(
    cc@data.signaling
  )

  same_cell_count <-
    length(fresh_cells) ==
      length(hist_cells)

  same_cell_set <-
    setequal(
      fresh_cells,
      hist_cells
    )

  same_cell_order <-
    identical(
      fresh_cells,
      hist_cells
    )

  # ------------------------------------------------------------------------
  # Signaling matrix
  # ------------------------------------------------------------------------

  hist_signaling <- historical@data.signaling[
    ,
    hist_cells,
    drop = FALSE
  ]

  fresh_signaling <- cc@data.signaling

  same_signal_dim <- identical(
    dim(fresh_signaling),
    dim(hist_signaling)
  )

  same_signal_gene_order <- identical(
    rownames(fresh_signaling),
    rownames(hist_signaling)
  )

  same_signal_cell_order <- identical(
    colnames(fresh_signaling),
    colnames(hist_signaling)
  )

  signal_max_delta <- NA_real_
  signal_mean_delta <- NA_real_
  signaling_numeric_match <- FALSE

  if (
    same_signal_dim &&
    same_signal_gene_order &&
    same_signal_cell_order
  ) {

    delta <- abs(
      as.numeric(fresh_signaling) -
        as.numeric(hist_signaling)
    )

    signal_max_delta <- max(
      delta,
      na.rm = TRUE
    )

    signal_mean_delta <- mean(
      delta,
      na.rm = TRUE
    )

    signaling_numeric_match <- isTRUE(
      all.equal(
        as.numeric(fresh_signaling),
        as.numeric(hist_signaling),
        tolerance = 1e-12,
        check.attributes = FALSE
      )
    )
  }

  # ------------------------------------------------------------------------
  # Group identity comparison
  # ------------------------------------------------------------------------

  fresh_group <- as.character(
    meta[
      fresh_cells,
      "group_cellchat"
    ]
  )

  hist_group <- as.character(
    historical@meta[
      hist_cells,
      "group_cellchat"
    ]
  )

  same_group_order <- identical(
    fresh_group,
    hist_group
  )

  validation <- data.table(
    condition = ds,

    fresh_cells =
      length(fresh_cells),

    historical_cells =
      length(hist_cells),

    same_cell_count =
      same_cell_count,

    same_cell_set =
      same_cell_set,

    same_cell_order =
      same_cell_order,

    fresh_signaling_genes =
      nrow(fresh_signaling),

    historical_signaling_genes =
      nrow(hist_signaling),

    same_signaling_dimensions =
      same_signal_dim,

    same_signaling_gene_order =
      same_signal_gene_order,

    same_signaling_cell_order =
      same_signal_cell_order,

    signaling_numeric_match_1e12 =
      signaling_numeric_match,

    max_abs_signaling_delta =
      signal_max_delta,

    mean_abs_signaling_delta =
      signal_mean_delta,

    fresh_LR =
      length(fresh_lr),

    historical_LR =
      length(hist_lr),

    same_LR_count =
      same_lr_count,

    same_LR_set =
      same_lr_set,

    same_LR_order =
      same_lr_order,

    same_group_order =
      same_group_order
  )

  print(validation)

  required_pass <- c(
    same_cell_count,
    same_cell_set,
    same_cell_order,
    same_signal_dim,
    same_signal_gene_order,
    same_signal_cell_order,
    signaling_numeric_match,
    same_lr_count,
    same_lr_set,
    same_lr_order,
    same_group_order
  )

  if (!all(required_pass)) {

    if (!same_lr_order) {

      lr_cmp <- data.table(
        index = seq_len(
          max(
            length(fresh_lr),
            length(hist_lr)
          )
        ),
        fresh = c(
          fresh_lr,
          rep(
            NA_character_,
            max(
              0,
              length(hist_lr) -
                length(fresh_lr)
            )
          )
        ),
        historical = c(
          hist_lr,
          rep(
            NA_character_,
            max(
              0,
              length(fresh_lr) -
                length(hist_lr)
            )
          )
        )
      )

      fwrite(
        lr_cmp,
        file.path(
          outdir,
          paste0(
            ds,
            "_LR_order_comparison.tsv"
          )
        ),
        sep = "\t"
      )
    }

    stop(
      ds,
      ": production pre-inference state does not exactly reproduce historical state."
    )
  }

  validation
}

young_validation <- prepare_and_validate(
  obj_young,
  "Young"
)

rm(obj_young)
gc()

aged_validation <- prepare_and_validate(
  obj_aged,
  "Aged"
)

rm(obj_aged)
gc()

validation <- rbindlist(
  list(
    young_validation,
    aged_validation
  ),
  use.names = TRUE
)

fwrite(
  validation,
  file.path(
    outdir,
    "CellChat_production_preinference_validation.tsv"
  ),
  sep = "\t"
)

cat("\n============================================================\n")
cat("PRODUCTION PREFLIGHT COMPLETE\n")
cat("============================================================\n\n")

print(validation)

cat(
  "\nPASS: production preprocessing exactly reproduces the ",
  "validated historical CellChat state.\n",
  sep = ""
)

cat(
  "\nNo computeCommunProb() call was made and no canonical file was modified.\n"
)
