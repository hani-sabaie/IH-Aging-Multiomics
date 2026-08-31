# ===== Clean environment =====
rm(list = ls(all.names = TRUE))
gc()

# ===== Libraries =====
library(data.table)
library(dplyr)
library(Seurat)
library(Signac)
library(hdWGCNA)

# ===== Resolve repo root =====
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (length(file_arg) > 0) {
  script_path <- normalizePath(sub("^--file=", "", file_arg[1]))
  script_dir  <- dirname(script_path)
  repo_root   <- normalizePath(file.path(script_dir, ".."))
} else {
  repo_root <- normalizePath(".")
}

srcdir <- file.path(
  repo_root,
  "source_data",
  "figure_source_data",
  "Figure_16"
)
dir.create(srcdir, recursive = TRUE, showWarnings = FALSE)

obj_file <- Sys.getenv(
  "FIG16_TFNET_B_RDS",
  unset = file.path(repo_root, "outputs", "hdWGCNA_TFNet_B_obj.rds")
)

if (!file.exists(obj_file)) {
  stop("Required TF-network object not found: ", obj_file)
}

cat("Using TF-network object:\n", obj_file, "\n\n")

obj <- readRDS(obj_file)

pick_col <- function(nms, candidates, what) {
  hit <- candidates[candidates %in% nms]
  if (length(hit) == 0) {
    stop(
      "Could not find ", what, " column.\nAvailable columns:\n",
      paste(nms, collapse = ", ")
    )
  }
  hit[1]
}

# ============================================================================
# Figure 16A: RegulonBarPlot source data (SMAD3, Strategy B)
# ============================================================================

tf_regulons <- as.data.table(GetTFRegulons(obj))

tf_col     <- pick_col(names(tf_regulons), c("tf", "TF", "source"), "TF/source")
target_col <- pick_col(names(tf_regulons), c("gene", "target_gene", "target"), "target gene")

smad3_regulons <- copy(tf_regulons[get(tf_col) == "SMAD3"])

if ("reg_score" %in% names(smad3_regulons)) {
  smad3_regulons[, pass_reg_score_cutoff_0_10 := reg_score >= 0.10]
}
if ("reg_score_signed" %in% names(smad3_regulons)) {
  smad3_regulons[, pass_abs_reg_score_signed_cutoff_0_10 := abs(reg_score_signed) >= 0.10]
}
if ("Cor" %in% names(smad3_regulons)) {
  smad3_regulons[, pass_abs_Cor_cutoff_0_10 := abs(Cor) >= 0.10]
}

smad3_regulons[, figure_panel := "Figure16A"]

fwrite(
  smad3_regulons,
  file.path(srcdir, "Figure16A_SMAD3_regulon_barplot_source_data.csv")
)

# Optional filtered export matching displayed cutoff as closely as possible
if ("reg_score" %in% names(smad3_regulons)) {
  fwrite(
    smad3_regulons[reg_score >= 0.10],
    file.path(srcdir, "Figure16A_SMAD3_regulon_barplot_cutoff0.10.csv")
  )
}

# ============================================================================
# Figure 16B: FeaturePlot source data
# ============================================================================

reductions_available <- names(obj@reductions)
reduction_to_use <- c("wnn.umap", "umap", "harmony.sct.umap")[
  c("wnn.umap", "umap", "harmony.sct.umap") %in% reductions_available
][1]

if (is.na(reduction_to_use) || is.null(reduction_to_use)) {
  stop(
    "No suitable UMAP reduction found. Available reductions:\n",
    paste(reductions_available, collapse = ", ")
  )
}

cat("Using reduction for Figure 16B:", reduction_to_use, "\n")

emb <- Embeddings(obj, reduction = reduction_to_use)
umap_dt <- as.data.table(emb, keep.rownames = "cell")

if (ncol(umap_dt) < 3) {
  stop("Unexpected embedding matrix shape for reduction: ", reduction_to_use)
}

setnames(
  umap_dt,
  old = names(umap_dt)[2:3],
  new = c("UMAP_1", "UMAP_2")
)

meta_keep <- intersect(
  c("condition", "skeletal_muscle", "sample", "orig.ident"),
  colnames(obj@meta.data)
)

meta_dt <- as.data.table(obj@meta.data, keep.rownames = "cell")[, c("cell", meta_keep), with = FALSE]

expr_dt <- as.data.table(
  FetchData(obj, vars = c("SMAD3")),
  keep.rownames = "cell"
)
setnames(expr_dt, "SMAD3", "SMAD3_expression")

pos_regulon_scores <- GetRegulonScores(obj, target_type = "positive")
neg_regulon_scores <- GetRegulonScores(obj, target_type = "negative")

if (!"SMAD3" %in% colnames(pos_regulon_scores)) {
  stop("SMAD3 not found in positive regulon score matrix.")
}
if (!"SMAD3" %in% colnames(neg_regulon_scores)) {
  stop("SMAD3 not found in negative regulon score matrix.")
}

score_dt <- data.table(
  cell = rownames(pos_regulon_scores),
  positive_regulon_score = as.numeric(pos_regulon_scores[, "SMAD3"]),
  negative_regulon_score = as.numeric(neg_regulon_scores[, "SMAD3"])
)

fig16b_dt <- Reduce(
  function(x, y) merge(x, y, by = "cell", all = FALSE),
  list(umap_dt[, .(cell, UMAP_1, UMAP_2)], meta_dt, expr_dt, score_dt)
)

fig16b_dt[, reduction_used := reduction_to_use]

fwrite(
  fig16b_dt,
  file.path(srcdir, "Figure16B_SMAD3_regulon_featureplot_source_data.csv")
)

# Long format (often easier for journals)
fig16b_long <- rbindlist(list(
  copy(fig16b_dt)[, .(
    cell, UMAP_1, UMAP_2,
    condition = if ("condition" %in% names(fig16b_dt)) condition else NA,
    skeletal_muscle = if ("skeletal_muscle" %in% names(fig16b_dt)) skeletal_muscle else NA,
    feature = "SMAD3_expression",
    value = SMAD3_expression
  )],
  copy(fig16b_dt)[, .(
    cell, UMAP_1, UMAP_2,
    condition = if ("condition" %in% names(fig16b_dt)) condition else NA,
    skeletal_muscle = if ("skeletal_muscle" %in% names(fig16b_dt)) skeletal_muscle else NA,
    feature = "positive_regulon_score",
    value = positive_regulon_score
  )],
  copy(fig16b_dt)[, .(
    cell, UMAP_1, UMAP_2,
    condition = if ("condition" %in% names(fig16b_dt)) condition else NA,
    skeletal_muscle = if ("skeletal_muscle" %in% names(fig16b_dt)) skeletal_muscle else NA,
    feature = "negative_regulon_score",
    value = negative_regulon_score
  )]
))

fwrite(
  fig16b_long,
  file.path(srcdir, "Figure16B_SMAD3_regulon_featureplot_source_data_long.csv")
)

# ============================================================================
# Figure 16C: TFNetworkPlot source data
# ============================================================================

tf_network <- as.data.table(GetTFNetwork(obj))
net_tf_col <- pick_col(names(tf_network), c("tf", "TF", "source"), "TF/source in TF network")
net_target_col <- pick_col(names(tf_network), c("gene", "target_gene", "target"), "target gene in TF network")

smad3_network <- copy(tf_network[get(net_tf_col) == "SMAD3"])

hub_df <- as.data.table(GetHubGenes(obj, n_hubs = 20))
if (!all(c("module", "gene_name") %in% names(hub_df))) {
  stop(
    "GetHubGenes output does not have expected columns. Available columns:\n",
    paste(names(hub_df), collapse = ", ")
  )
}

brown_hubs <- copy(hub_df[module == "brown"])
brown_hub_genes <- unique(brown_hubs$gene_name)

smad3_network[, is_brown_hub_label := get(net_target_col) %in% brown_hub_genes]

fwrite(
  smad3_network,
  file.path(srcdir, "Figure16C_SMAD3_TFNetwork_source_data.csv")
)

fwrite(
  brown_hubs,
  file.path(srcdir, "Figure16C_brown_hub_genes_label_source_data.csv")
)

# Also export SMAD3 regulons annotated by brown hub labels
smad3_regulons_labeled <- copy(smad3_regulons)
smad3_regulons_labeled[, is_brown_hub_label := get(target_col) %in% brown_hub_genes]

fwrite(
  smad3_regulons_labeled,
  file.path(srcdir, "Figure16C_SMAD3_regulons_with_brown_hub_labels.csv")
)

# ============================================================================
# Summary
# ============================================================================

cat("\nFigure 16 source data generated.\n")
cat("Figure16A SMAD3 regulon rows:", nrow(smad3_regulons), "\n")
cat("Figure16B cell rows:", nrow(fig16b_dt), "\n")
cat("Figure16B long rows:", nrow(fig16b_long), "\n")
cat("Figure16C SMAD3 network rows:", nrow(smad3_network), "\n")
cat("Figure16C brown hub genes:", length(brown_hub_genes), "\n")
cat("Done.\n")
