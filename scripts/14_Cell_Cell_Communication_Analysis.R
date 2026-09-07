# ========================
# Setting up environment
# ========================
# ===== Clean environment =====
rm(list = ls(all.names = TRUE))
gc()

# ===== Loading relevant libraries =====
library(Seurat)
library(ggplot2)
library(CellChat)
library(data.table)

# ===== set seed =====
set.seed(1234)

# ===== Repository paths =====
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (length(file_arg) == 1) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg)))
  repo_root <- normalizePath(file.path(script_dir, ".."))
} else {
  repo_root <- normalizePath(".")
}

outdir <- file.path(repo_root, "outputs")

processed_dir <- file.path(
  repo_root,
  "processed_results",
  "12_CellChat"
)

figdir <- Sys.getenv(
  "CELLCHAT_FIGDIR",
  unset = file.path(
    outdir,
    "CellChat"
  )
)

input_obj_file <- Sys.getenv(
  "CELLCHAT_INPUT_RDS",
  unset = file.path(
    outdir,
    "decont_merged_filt_nodoub_cc_sct_reduc_clust_integ_annot_obj.rds"
  )
)

cellchat_db_order_file <- Sys.getenv(
  "CELLCHAT_DB_ORDER_FILE",
  unset = file.path(
    processed_dir,
    "CellChatDB_human_18pathway_interaction_order.tsv"
  )
)

cellchat_file <- Sys.getenv(
  "CELLCHAT_OUTPUT_RDS",
  unset = file.path(
    processed_dir,
    "cellchat_merged.rds"
  )
)

multiple_testing_dir <- Sys.getenv(
  "CELLCHAT_MT_DIR",
  unset = file.path(
    processed_dir,
    "multiple_testing"
  )
)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(
  multiple_testing_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

save_gg <- function(p, filename, w=7, h=5, dpi=300){
  if (inherits(p, "ggplot") || inherits(p, "patchwork")){
    ggsave(file.path(figdir, filename), p, width=w, height=h, units="in", dpi=dpi, bg="white")
  } else {
    png(file.path(figdir, filename), width=w, height=h, units="in", res=dpi)
    print(p); dev.off()
  }
}

# ============================================================================ #

# ========================
# Load reproducibility resources
# ========================

if (!file.exists(input_obj_file)) {
  stop(
    "CellChat Seurat input not found: ",
    input_obj_file,
    "\nSet CELLCHAT_INPUT_RDS to its local path."
  )
}

if (!file.exists(cellchat_db_order_file)) {
  stop(
    "CellChat DB-order resource not found: ",
    cellchat_db_order_file
  )
}

cellchat_db_order <- fread(
  cellchat_db_order_file
)

required_db_order_cols <- c(
  "db_order",
  "interaction_name",
  "pathway_name"
)

missing_db_order_cols <- setdiff(
  required_db_order_cols,
  names(cellchat_db_order)
)

if (length(missing_db_order_cols) > 0L) {
  stop(
    "CellChat DB-order resource is missing columns: ",
    paste(
      missing_db_order_cols,
      collapse = ", "
    )
  )
}

if (
  nrow(cellchat_db_order) != 518L ||
  uniqueN(cellchat_db_order$interaction_name) != 518L ||
  uniqueN(cellchat_db_order$pathway_name) != 18L
) {
  stop(
    "Unexpected CellChat DB-order resource dimensions/content."
  )
}

if (
  !identical(
    as.integer(cellchat_db_order$db_order),
    seq_len(518L)
  )
) {
  stop(
    "CellChat DB-order resource has invalid db_order indexing."
  )
}

# ========================
# Load the object
# ========================
obj <- readRDS(input_obj_file)

# ========================
# Cell-cell communication analysis
# ========================
options(future.globals.maxSize = 8 * 1024^3)

# Each LR batch is evaluated sequentially. This avoids the large Windows
# multisession overhead observed for the full LR family and reproduces the
# reviewer-C7 validation runs exactly.
future::plan("sequential")
obj$samples <- obj$sample
obj$cc_group <- NA_character_
obj$cc_group[obj$skeletal_muscle %in% c("FAP1","FAP2","FAP4")] <- "Other FAPs"
obj$cc_group[obj$skeletal_muscle == "FAP3"] <- "FAP3"
obj$cc_group[obj$skeletal_muscle == "MuSC"] <- "MuSC"
obj$cc_group[obj$skeletal_muscle == "MSM"] <- "Myogenic"
obj$cc_group[obj$skeletal_muscle %in% c("Macrophage","B/T/NK")] <- "Immune"
obj$cc_group[obj$skeletal_muscle %in% c("EC1","EC2","EC3")] <- "Endothelial"
obj$cc_group[obj$skeletal_muscle %in% c("Pericyte","SMC1","SMC2")] <- "Vascular stromal"

obj_young <- subset(obj, subset = condition == "Young")
obj_aged  <- subset(obj, subset = condition == "Aged")
rm(obj)

obj_young$group_cellchat <- obj_young$cc_group
obj_aged$group_cellchat  <- obj_aged$cc_group

# ============================================================================
# CellChat communication inference with reviewer-C7 multiple-testing control
#
# Primary inference:
#   - triMean
#   - nboot = 1000
#   - seed.use = 1
#   - plus-one empirical permutation P = (b + 1) / (B + 1)
#   - Benjamini-Hochberg correction across the COMPLETE condition-specific
#     LR x source x target family
#
# The LR family is evaluated in deterministic batches to avoid the large
# Windows multisession overhead of a monolithic nboot=1000 run.
# ============================================================================

cellchat_group_order <- c(
  "Endothelial",
  "FAP3",
  "Immune",
  "MuSC",
  "Myogenic",
  "Other FAPs",
  "Vascular stromal"
)

expected_lr_by_condition <- c(
  Young = 488L,
  Aged = 468L
)

B_cellchat <- 1000L

lr_batch_size <- suppressWarnings(
  as.integer(
    Sys.getenv(
      "CELLCHAT_LR_BATCH_SIZE",
      unset = "25"
    )
  )
)

if (
  length(lr_batch_size) != 1L ||
  is.na(lr_batch_size) ||
  lr_batch_size < 1L
) {
  stop(
    "CELLCHAT_LR_BATCH_SIZE must be a positive integer."
  )
}

make_cellchat <- function(
    seu,
    condition_name
) {

  if (
    !(condition_name %in% names(expected_lr_by_condition))
  ) {
    stop(
      "Unexpected CellChat condition: ",
      condition_name
    )
  }

  data.input <- GetAssayData(
    seu,
    assay = "SCT",
    slot = "data"
  )

  meta <- seu@meta.data

  meta$group_cellchat <- factor(
    meta$group_cellchat,
    levels = cellchat_group_order
  )

  if (anyNA(meta$group_cellchat)) {
    stop(
      condition_name,
      ": CellChat group assignment contains NA."
    )
  }

  # ------------------------------------------------------------------------
  # Establish the condition-specific signaling data and exact LR family.
  # ------------------------------------------------------------------------

  cellchat <- createCellChat(
    object = data.input,
    meta = meta,
    group.by = "group_cellchat",
    assay = "SCT"
  )

  CellChatDB <- CellChatDB.human

  pathways.use <- c(
    # SMAD3 / fibrosis axis
    "TGFb",

    # ECM structural components
    "COLLAGEN",
    "FN1",
    "LAMININ",
    "THBS",
    "TENASCIN",
    "VTN",
    "HSPG",

    # Cell-cell / cell-matrix adhesion
    "CDH",
    "CDH1",
    "CDH5",
    "PCDH",
    "ICAM",
    "VCAM",
    "JAM",

    # ECM remodeling
    "MMP",

    # FAP / fibro-inflammatory markers
    "SPP1",
    "PERIOSTIN"
  )

  CellChatDB.use <- subsetDB(
    CellChatDB,
    search = pathways.use,
    key = "pathway_name"
  )

  fresh_db_names <- as.character(
    CellChatDB.use$interaction$interaction_name
  )

  expected_db_names <- as.character(
    cellchat_db_order$interaction_name
  )

  if (
    length(fresh_db_names) != 518L ||
    !setequal(
      fresh_db_names,
      expected_db_names
    )
  ) {
    stop(
      condition_name,
      ": installed CellChatDB interaction content does not match ",
      "the validated 18-pathway reproducibility resource."
    )
  }

  db_idx <- match(
    expected_db_names,
    fresh_db_names
  )

  if (anyNA(db_idx)) {
    stop(
      condition_name,
      ": failed to map CellChatDB to the validated interaction order."
    )
  }

  CellChatDB.use$interaction <-
    CellChatDB.use$interaction[
      db_idx,
      ,
      drop = FALSE
    ]

  if (
    !identical(
      as.character(
        CellChatDB.use$interaction$interaction_name
      ),
      expected_db_names
    )
  ) {
    stop(
      condition_name,
      ": validated CellChatDB interaction order was not reproduced."
    )
  }

  if (
    !identical(
      as.character(
        CellChatDB.use$interaction$pathway_name
      ),
      as.character(
        cellchat_db_order$pathway_name
      )
    )
  ) {
    stop(
      condition_name,
      ": CellChatDB pathway annotations differ from the ",
      "validated reproducibility resource."
    )
  }

  cellchat@DB <- CellChatDB.use

  cellchat <- subsetData(cellchat)

  cellchat <- identifyOverExpressedGenes(
    cellchat
  )

  cellchat <- identifyOverExpressedInteractions(
    cellchat
  )

  lr_full <- cellchat@LR$LRsig

  if (
    is.null(lr_full) ||
    nrow(lr_full) == 0L ||
    !("interaction_name" %in% names(lr_full))
  ) {
    stop(
      condition_name,
      ": no valid LRsig family was generated."
    )
  }

  lr_names <- as.character(
    lr_full$interaction_name
  )

  if (anyDuplicated(lr_names) > 0L) {
    stop(
      condition_name,
      ": duplicated LR interaction names detected."
    )
  }

  expected_n_lr <-
    unname(
      expected_lr_by_condition[
        condition_name
      ]
    )

  if (length(lr_names) != expected_n_lr) {
    stop(
      condition_name,
      ": expected ",
      expected_n_lr,
      " LR interactions; found ",
      length(lr_names),
      "."
    )
  }

  signaling_data <- cellchat@data.signaling

  meta_batch <- meta[
    colnames(signaling_data),
    ,
    drop = FALSE
  ]

  if (
    !identical(
      colnames(signaling_data),
      rownames(meta_batch)
    )
  ) {
    stop(
      condition_name,
      ": signaling matrix and metadata are misaligned."
    )
  }

  meta_batch$group_cellchat <- factor(
    meta_batch$group_cellchat,
    levels = cellchat_group_order
  )

  if (anyNA(meta_batch$group_cellchat)) {
    stop(
      condition_name,
      ": reconstructed CellChat groups contain NA."
    )
  }

  groups <- cellchat_group_order

  # ------------------------------------------------------------------------
  # Allocate the complete formal LR x source x target family.
  # ------------------------------------------------------------------------

  prob_all <- array(
    0,
    dim = c(
      length(groups),
      length(groups),
      length(lr_names)
    ),
    dimnames = list(
      groups,
      groups,
      lr_names
    )
  )

  p_raw_all <- array(
    1,
    dim = dim(prob_all),
    dimnames = dimnames(prob_all)
  )

  batch_starts <- seq.int(
    1L,
    length(lr_names),
    by = lr_batch_size
  )

  batch_elapsed <- numeric(
    length(batch_starts)
  )

  cat("\n============================================================\n")
  cat("CELLCHAT PRODUCTION INFERENCE:", condition_name, "\n")
  cat("============================================================\n\n")

  cat("Cells          :", ncol(signaling_data), "\n")
  cat("Signaling genes:", nrow(signaling_data), "\n")
  cat("LR family      :", length(lr_names), "\n")
  cat("Batch size     :", lr_batch_size, "\n")
  cat("Batches        :", length(batch_starts), "\n")
  cat("nboot          :", B_cellchat, "\n\n")

  # ------------------------------------------------------------------------
  # Validated deterministic LR-batch inference.
  # ------------------------------------------------------------------------

  for (batch_number in seq_along(batch_starts)) {

    start_idx <- batch_starts[batch_number]

    end_idx <- min(
      start_idx + lr_batch_size - 1L,
      length(lr_names)
    )

    batch_indices <- seq.int(
      start_idx,
      end_idx
    )

    lr_batch <- lr_full[
      batch_indices,
      ,
      drop = FALSE
    ]

    batch_lr_names <- lr_names[
      batch_indices
    ]

    if (
      !identical(
        as.character(
          lr_batch$interaction_name
        ),
        batch_lr_names
      )
    ) {
      stop(
        condition_name,
        ": LR batch ordering mismatch."
      )
    }

    cat(
      "Running batch ",
      batch_number,
      "/",
      length(batch_starts),
      " [LR ",
      start_idx,
      "-",
      end_idx,
      "]...\n",
      sep = ""
    )

    cc <- createCellChat(
      object = signaling_data,
      meta = meta_batch,
      group.by = "group_cellchat"
    )

    cc@DB <- CellChatDB.use
    cc@data.signaling <- signaling_data

    cc@idents <- factor(
      meta_batch$group_cellchat,
      levels = groups
    )

    names(cc@idents) <- rownames(
      meta_batch
    )

    tm <- system.time({

      cc <- computeCommunProb(
        cc,
        type = "triMean",
        LR.use = lr_batch,
        raw.use = TRUE,
        population.size = FALSE,
        nboot = B_cellchat,
        seed.use = 1L
      )

    })

    batch_elapsed[batch_number] <-
      unname(
        tm["elapsed"]
      )

    cc <- filterCommunication(
      cc,
      min.cells = 10
    )

    new_prob <- cc@net$prob
    new_p <- cc@net$pval

    new_lr <- dimnames(
      new_prob
    )[[3]]

    new_idx <- match(
      batch_lr_names,
      new_lr
    )

    if (anyNA(new_idx)) {
      stop(
        condition_name,
        ": requested LR interaction absent from batch result."
      )
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

    dimnames(new_prob)[[3]] <-
      batch_lr_names

    dimnames(new_p)[[3]] <-
      batch_lr_names

    if (
      !identical(
        dimnames(new_prob)[[1]],
        groups
      ) ||
      !identical(
        dimnames(new_prob)[[2]],
        groups
      ) ||
      !identical(
        dimnames(new_prob)[[3]],
        batch_lr_names
      )
    ) {
      stop(
        condition_name,
        ": batch source/target/LR ordering mismatch."
      )
    }

    prob_all[
      ,
      ,
      batch_indices
    ] <- new_prob

    p_raw_all[
      ,
      ,
      batch_indices
    ] <- new_p
  }

  # ------------------------------------------------------------------------
  # Formal hypothesis table.
  # ------------------------------------------------------------------------

  tests <- as.data.table(
    as.table(
      p_raw_all
    )
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
      as.numeric(
        prob_all
      )
  ]

  tests[
    ,
    permutation_exceedance_count :=
      as.integer(
        pmin(
          B_cellchat,
          pmax(
            0,
            round(
              p_raw_nboot1000 *
                B_cellchat
            )
          )
        )
      )
  ]

  tests[
    ,
    p_plus1 :=
      (
        permutation_exceedance_count +
          1
      ) /
      (
        B_cellchat +
          1
      )
  ]

  tests[
    ,
    `:=`(
      dataset = condition_name,
      LR_family_index =
        match(
          interaction_name,
          lr_names
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

  expected_n_tests <-
    expected_n_lr *
      length(groups) *
      length(groups)

  if (nrow(tests) != expected_n_tests) {
    stop(
      condition_name,
      ": expected ",
      expected_n_tests,
      " formal tests; found ",
      nrow(tests),
      "."
    )
  }

  lr_indices <- sort(
    unique(
      tests$LR_family_index
    )
  )

  if (
    !identical(
      lr_indices,
      seq_len(expected_n_lr)
    )
  ) {
    stop(
      condition_name,
      ": incomplete LR family."
    )
  }

  lr_counts <- tests[
    ,
    .N,
    by = LR_family_index
  ]

  if (
    nrow(lr_counts) != expected_n_lr ||
    any(
      lr_counts$N !=
        length(groups)^2
    )
  ) {
    stop(
      condition_name,
      ": each LR must contribute exactly ",
      length(groups)^2,
      " source-target tests."
    )
  }

  key_cols <- c(
    "LR_family_index",
    "interaction_name",
    "source",
    "target"
  )

  if (
    anyDuplicated(
      tests[
        ,
        ..key_cols
      ]
    ) > 0L
  ) {
    stop(
      condition_name,
      ": duplicate LR-source-target hypothesis detected."
    )
  }

  if (
    anyNA(tests$p_plus1) ||
    any(
      tests$p_plus1 <= 0 |
        tests$p_plus1 > 1
    )
  ) {
    stop(
      condition_name,
      ": invalid plus-one empirical P values."
    )
  }

  # ------------------------------------------------------------------------
  # Multiple-testing correction.
  # ------------------------------------------------------------------------

  # Sensitivity analysis: BH applied directly to CellChat b/B values.
  tests[
    ,
    p_raw_BH_condition :=
      p.adjust(
        p_raw_nboot1000,
        method = "BH"
      )
  ]

  # Primary analysis: finite plus-one empirical P followed by BH across the
  # complete condition-specific LR x source x target family.
  tests[
    ,
    p_plus1_BH_condition :=
      p.adjust(
        p_plus1,
        method = "BH"
      )
  ]

  tests[
    ,
    `:=`(
      raw_nboot1000_sig =
        p_raw_nboot1000 < 0.05,

      raw_BH_sig =
        p_raw_BH_condition < 0.05,

      plus1_BH_sig =
        p_plus1_BH_condition < 0.05
    )
  ]

  # ------------------------------------------------------------------------
  # Reconstruct the multiplicity-corrected CellChat p-value array.
  # ------------------------------------------------------------------------

  p_bh_all <- array(
    1,
    dim = dim(prob_all),
    dimnames = dimnames(prob_all)
  )

  src_i <- match(
    as.character(tests$source),
    groups
  )

  tgt_i <- match(
    as.character(tests$target),
    groups
  )

  lr_i <- match(
    as.character(tests$interaction_name),
    lr_names
  )

  if (
    anyNA(src_i) ||
    anyNA(tgt_i) ||
    anyNA(lr_i)
  ) {
    stop(
      condition_name,
      ": failed to map corrected formal hypotheses back to the network array."
    )
  }

  p_bh_all[
    cbind(
      src_i,
      tgt_i,
      lr_i
    )
  ] <- tests$p_plus1_BH_condition

  if (anyNA(p_bh_all)) {
    stop(
      condition_name,
      ": corrected network P array contains NA."
    )
  }

  retained_lr <- lr_names[
    apply(
      prob_all *
        (
          p_bh_all <
            0.05
        ),
      3,
      sum,
      na.rm = TRUE
    ) > 0
  ]

  # The base CellChat object already contains the exact full LRsig family,
  # signaling matrix, metadata, identities, and database.
  cellchat@net$prob <- prob_all
  cellchat@net$pval <- p_bh_all
  cellchat@net$LRs <- retained_lr

  cellchat <- computeCommunProbPathway(
    cellchat
  )

  cellchat <- aggregateNet(
    cellchat
  )

  cellchat <- netAnalysis_computeCentrality(
    cellchat,
    slot.name = "netP"
  )

  # ------------------------------------------------------------------------
  # Reproducibility outputs.
  # ------------------------------------------------------------------------

  setorder(
    tests,
    LR_family_index,
    source,
    target
  )

  summary_dt <- data.table(
    dataset = condition_name,
    n_batches = length(batch_starts),
    n_LR = expected_n_lr,
    n_source_groups = length(groups),
    n_target_groups = length(groups),
    n_tests = nrow(tests),

    n_raw_p_zero =
      sum(
        tests$p_raw_nboot1000 ==
          0
      ),

    n_raw_p_lt_0_05 =
      sum(
        tests$p_raw_nboot1000 <
          0.05
      ),

    n_raw_BH_lt_0_05 =
      sum(
        tests$p_raw_BH_condition <
          0.05
      ),

    n_plus1_BH_lt_0_05 =
      sum(
        tests$p_plus1_BH_condition <
          0.05
      ),

    n_retained_LR =
      length(retained_lr),

    n_retained_pathways =
      length(
        cellchat@netP$pathways
      ),

    min_raw_p =
      min(
        tests$p_raw_nboot1000
      ),

    min_plus1_p =
      min(
        tests$p_plus1
      ),

    min_plus1_BH =
      min(
        tests$p_plus1_BH_condition
      ),

    total_batch_elapsed_seconds =
      sum(
        batch_elapsed
      )
  )

  fwrite(
    tests,
    file.path(
      multiple_testing_dir,
      paste0(
        condition_name,
        "_full18_nboot1000_conditionwide_BH_tests.tsv"
      )
    ),
    sep = "\t"
  )

  fwrite(
    summary_dt,
    file.path(
      multiple_testing_dir,
      paste0(
        condition_name,
        "_full18_nboot1000_conditionwide_BH_summary.tsv"
      )
    ),
    sep = "\t"
  )

  cat("\n")
  print(summary_dt)
  cat("\n")

  list(
    object = cellchat,
    tests = tests,
    summary = summary_dt
  )
}

young_result <- make_cellchat(
  obj_young,
  "Young"
)

aged_result <- make_cellchat(
  obj_aged,
  "Aged"
)

cellchat_young <- young_result$object
cellchat_aged <- aged_result$object

cellchat_multiple_testing_summary <- rbindlist(
  list(
    young_result$summary,
    aged_result$summary
  ),
  use.names = TRUE
)

fwrite(
  cellchat_multiple_testing_summary,
  file.path(
    multiple_testing_dir,
    "CellChat_full18_nboot1000_conditionwide_BH_summary.tsv"
  ),
  sep = "\t"
)

object.list <- list(
  Young = cellchat_young,
  Aged = cellchat_aged
)

cellchat_merged <- mergeCellChat(
  object.list = object.list,
  add.names = names(object.list)
)

saveRDS(
  cellchat_merged,
  file = cellchat_file
)

cat(
  "Corrected merged CellChat object written to:\n",
  cellchat_file,
  "\n"
)

# Visualization
# pathways
pathway_TGF <- "TGFb"

all_pw <- unique(CellChatDB.human$interaction$pathway_name)

# Circle plot (overall signalling)
plot_circle_png <- function(cellchat_obj, title_str) {
  mat <- cellchat_obj@net$weight
  group_order <- c("FAP3","Other FAPs","MuSC","Myogenic","Immune","Endothelial","Vascular stromal")
  go <- intersect(group_order, rownames(mat))
  mat_reorder <- mat[go, go, drop = FALSE]
  vertex.weight <- rowSums(mat_reorder)
  
  netVisual_circle(
    mat_reorder,
    vertex.weight = vertex.weight,
    weight.scale = TRUE,
    label.edge  = FALSE,
    edge.weight.max = max(mat_reorder),
    vertex.label.cex = 1.4,
    vertex.size.max  = 15
  )
  grid::grid.text(title_str, y = unit(0.97, "npc"),
                  gp = grid::gpar(fontsize = 18, fontface = "bold"))
}

png(file.path(figdir, "circle_overall_Young.png"),
    width = 1800, height = 1800, res = 300)
plot_circle_png(cellchat_young, "Young - overall signaling")
dev.off()

png(file.path(figdir, "circle_overall_Aged.png"),
    width = 1800, height = 1800, res = 300)
plot_circle_png(cellchat_aged,  "Aged - overall signaling")
dev.off()

# Circos plot (pathways)
png(file.path(figdir, "circle_TGFb_Young.png"),
    width = 2000, height = 2000, res = 300)

netVisual_aggregate(
  object = cellchat_young,
  signaling = pathway_TGF,
  layout = "circle",
  thresh = 0.05
)

dev.off()

png(file.path(figdir, "circle_TGFb_Aged.png"),
    width = 2000, height = 2000, res = 300)

netVisual_aggregate(
  object = cellchat_aged,
  signaling = pathway_TGF,
  layout = "circle",
  thresh = 0.05
)

dev.off()

# Bubble plot
signaling.use <- pathway_TGF

p_bubble_FAP <- netVisual_bubble(
  object  = cellchat_merged,
  sources.use = c("FAP3","Other FAPs"),
  targets.use = c("MuSC","Myogenic","Immune","Endothelial","Vascular stromal"),
  signaling = signaling.use,
  angle.x = 45,
  remove.isolate = TRUE,
  comparison = c(1,2)
)

save_gg(
  p_bubble_FAP, "bubble_FAP_TGFb_Young_vs_Aged.png",
  w = 8, h = 5
)

# Signaling role scatter
# Common range for point size (based on number of links)
num.link <- sapply(object.list, function(x) {
  rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)
})
weight.MinMax <- c(min(num.link), max(num.link))

p_scatter_list <- lapply(names(object.list), function(nm) {
  netAnalysis_signalingRole_scatter(
    object  = object.list[[nm]],
    slot.name = "netP",
    title = nm,
    weight.MinMax = weight.MinMax
  )
})

p_scatter_combined <- p_scatter_list[[1]] + p_scatter_list[[2]]

save_gg(
  p_scatter_combined, "signalingRole_scatter_Young_vs_Aged.png",
  w = 10, h = 5
)

# Signaling heatmap
# Pathway union so that Young/Aged are comparable
pathway.union <- union(cellchat_young@netP$pathways,
                       cellchat_aged@netP$pathways)

# Outgoing 
ht_out_young <- netAnalysis_signalingRole_heatmap(
  object = cellchat_young,
  pattern = "outgoing",
  slot.name = "netP",
  signaling = pathway.union,
  title = "Young"
)

ht_out_aged <- netAnalysis_signalingRole_heatmap(
  object = cellchat_aged,
  pattern = "outgoing",
  slot.name = "netP",
  signaling = pathway.union,
  title = "Aged"
)

png(file.path(figdir, "heatmap_outgoing_netP_Young_vs_Aged.png"),
    width = 3200, height = 1600, res = 300)
ComplexHeatmap::draw(ht_out_young + ht_out_aged, ht_gap = unit(0.5, "cm"))
dev.off()

# Incoming 
ht_in_young <- netAnalysis_signalingRole_heatmap(
  object = cellchat_young,
  pattern = "incoming",
  slot.name = "netP",
  signaling = pathway.union,
  title = "Young"
)

ht_in_aged <- netAnalysis_signalingRole_heatmap(
  object = cellchat_aged,
  pattern = "incoming",
  slot.name = "netP",
  signaling = pathway.union,
  title = "Aged"
)

png(file.path(figdir, "heatmap_incoming_netP_Young_vs_Aged.png"),
    width = 3200, height = 1600, res = 240)
ComplexHeatmap::draw(ht_in_young + ht_in_aged, ht_gap = unit(0.5, "cm"))
dev.off()

# Contribution plots
# TGFb in Aged
p_contrib_TGF_aged <- netAnalysis_contribution(
  object = cellchat_aged,
  signaling = pathway_TGF
)
save_gg(
  p_contrib_TGF_aged, "contribution_TGFb_Aged.png",
  w = 7, h = 5
)
