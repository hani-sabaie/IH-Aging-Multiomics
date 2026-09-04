rm(list = ls(all.names = TRUE))
gc()

library(monocle3)
library(Seurat)
library(Matrix)
library(data.table)
library(dplyr)

set.seed(1234)

# -------------------------------------------------------------------------
# Repository root
# -------------------------------------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (length(file_arg) == 1) {
  script_path <- sub("^--file=", "", file_arg)
  script_dir <- dirname(normalizePath(script_path))
  repo_root <- normalizePath(file.path(script_dir, ".."))
} else {
  repo_root <- normalizePath(".")
}

# -------------------------------------------------------------------------
# Input / output
# -------------------------------------------------------------------------
obj_file <- Sys.getenv(
  "FIG19_HDWGCNA_RDS",
  unset = file.path(repo_root, "outputs", "hdWGCNA_obj.rds")
)

if (!file.exists(obj_file)) {
  stop(
    "hdWGCNA object not found: ",
    obj_file,
    "\nSet FIG19_HDWGCNA_RDS to its local path."
  )
}

srcdir <- file.path(
  repo_root,
  "source_data",
  "figure_source_data",
  "Figure_19"
)

dir.create(srcdir, recursive = TRUE, showWarnings = FALSE)

cat("Loading hdWGCNA object...\n")
obj <- readRDS(obj_file)

if (!"SMAD3" %in% rownames(obj)) {
  stop("SMAD3 not found in Seurat object.")
}

if (!"wnn.umap" %in% names(obj@reductions)) {
  stop("wnn.umap reduction not found.")
}

# -------------------------------------------------------------------------
# Reconstruct Monocle3 CDS exactly as in scripts/15_Monocle3_SMAD3.R
# -------------------------------------------------------------------------
cat("Building Monocle3 CDS...\n")

expr_mat <- GetAssayData(
  obj,
  assay = "SCT",
  slot = "data"
)

expr_mat <- as(expr_mat, "dgCMatrix")

cell_meta <- obj@meta.data

gene_anno <- data.frame(
  gene_short_name = rownames(expr_mat),
  row.names = rownames(expr_mat)
)

cds <- new_cell_data_set(
  expression_data = expr_mat,
  cell_metadata = cell_meta,
  gene_metadata = gene_anno
)

colData(cds)$cluster_fap <- cds$skeletal_muscle
colData(cds)$condition <- cds$condition

# -------------------------------------------------------------------------
# Restrict to FAP1-FAP3
# -------------------------------------------------------------------------
fap_levels <- c("FAP1", "FAP2", "FAP3")

keep_fap <- colnames(cds)[
  colData(cds)$cluster_fap %in% fap_levels
]

cds_fap <- cds[, keep_fap]

wnn_umap <- obj@reductions[["wnn.umap"]]@cell.embeddings
wnn_umap <- wnn_umap[colnames(cds_fap), , drop = FALSE]

reducedDims(cds_fap)$UMAP <- wnn_umap

# =========================================================================
# Figure 19A-B: Young-only trajectory
# =========================================================================
cat("\nReconstructing Young-only trajectory...\n")

keep_young <- colnames(cds_fap)[
  colData(cds_fap)$condition == "Young" &
    colData(cds_fap)$cluster_fap %in% fap_levels
]

cds_young <- cds_fap[, keep_young]

set.seed(1234)
cds_young <- cluster_cells(
  cds_young,
  reduction_method = "UMAP"
)

set.seed(1234)
cds_young <- learn_graph(
  cds_young,
  use_partition = FALSE
)

root_cells_young <- colnames(cds_young)[
  colData(cds_young)$cluster_fap == "FAP1"
]

if (length(root_cells_young) == 0) {
  stop("No FAP1 Young root cells found.")
}

cds_young <- order_cells(
  cds_young,
  root_cells = root_cells_young
)

pt_young <- monocle3::pseudotime(cds_young)
cds_young$pseudotime <- pt_young

smad3_expr_young <- as.numeric(
  exprs(cds_young)["SMAD3", ]
)

umap_young <- reducedDims(cds_young)$UMAP

fig19ab <- data.table(
  cell = colnames(cds_young),
  condition = as.character(colData(cds_young)$condition),
  FAP_cluster = as.character(colData(cds_young)$cluster_fap),
  pseudotime = as.numeric(pt_young),
  SMAD3_expression = smad3_expr_young,
  UMAP_1 = as.numeric(umap_young[, 1]),
  UMAP_2 = as.numeric(umap_young[, 2]),
  is_root_FAP1_cell = colnames(cds_young) %in% root_cells_young,
  finite_pseudotime = is.finite(pt_young)
)

setorder(fig19ab, pseudotime)

fwrite(
  fig19ab,
  file.path(
    srcdir,
    "Figure19AB_Young_SMAD3_pseudotime_source_data.csv"
  )
)

# -------------------------------------------------------------------------
# Moran's I for SMAD3 only
# -------------------------------------------------------------------------
cat("Running SMAD3 Moran's I test...\n")

smad3_gt <- graph_test(
  cds_young["SMAD3", ],
  neighbor_graph = "principal_graph",
  cores = 4
)

smad3_gt_dt <- as.data.table(
  as.data.frame(smad3_gt),
  keep.rownames = "gene_id"
)

fwrite(
  smad3_gt_dt,
  file.path(
    srcdir,
    "Figure19_SMAD3_Moran_test.csv"
  )
)

# =========================================================================
# Figure 19C-D: shared Young + Aged trajectory
# =========================================================================
cat("\nReconstructing Young + Aged trajectory...\n")

keep_ya <- colnames(cds_fap)[
  colData(cds_fap)$cluster_fap %in% fap_levels &
    colData(cds_fap)$condition %in% c("Young", "Aged")
]

cds_ya <- cds_fap[, keep_ya]

reducedDims(cds_ya)$UMAP <- wnn_umap[
  colnames(cds_ya),
  ,
  drop = FALSE
]

set.seed(1234)
cds_ya <- cluster_cells(
  cds_ya,
  reduction_method = "UMAP"
)

set.seed(1234)
cds_ya <- learn_graph(
  cds_ya,
  use_partition = FALSE
)

root_cells_ya <- colnames(cds_ya)[
  colData(cds_ya)$cluster_fap == "FAP1" &
    colData(cds_ya)$condition == "Young"
]

if (length(root_cells_ya) == 0) {
  stop("No FAP1 Young root cells found for combined trajectory.")
}

cds_ya <- order_cells(
  cds_ya,
  root_cells = root_cells_ya
)

pt_ya <- monocle3::pseudotime(cds_ya)
cds_ya$pseudotime <- pt_ya

smad3_expr_ya <- as.numeric(
  exprs(cds_ya)["SMAD3", ]
)

umap_ya <- reducedDims(cds_ya)$UMAP

fig19cd <- data.table(
  cell = colnames(cds_ya),
  condition = as.character(colData(cds_ya)$condition),
  FAP_cluster = as.character(colData(cds_ya)$cluster_fap),
  pseudotime = as.numeric(pt_ya),
  SMAD3_expression = smad3_expr_ya,
  UMAP_1 = as.numeric(umap_ya[, 1]),
  UMAP_2 = as.numeric(umap_ya[, 2]),
  is_root_FAP1_Young_cell = colnames(cds_ya) %in% root_cells_ya,
  finite_pseudotime = is.finite(pt_ya)
)

setorder(fig19cd, condition, pseudotime)

fwrite(
  fig19cd,
  file.path(
    srcdir,
    "Figure19CD_Young_Aged_SMAD3_pseudotime_source_data.csv"
  )
)

# -------------------------------------------------------------------------
# Condition effect model used in the original analysis
# -------------------------------------------------------------------------
cat("Running SMAD3 condition-effect model...\n")

cat(
  "Condition levels in reconstructed CDS:",
  paste(levels(factor(colData(cds_ya)$condition)), collapse = ", "),
  "\n"
)

# Explicitly use Aged as the reference level so that conditionYoung
# represents Young relative to Aged.
cds_ya$condition <- factor(
  cds_ya$condition,
  levels = c("Aged", "Young")
)

smad3_fit <- fit_models(
  cds_ya["SMAD3", ],
  model_formula_str =
    "~ condition + splines::ns(pseudotime, df = 3)"
)

coef_tab <- coefficient_table(smad3_fit)

cond_effect <- coef_tab[
  grepl("condition", coef_tab$term),
  ,
  drop = FALSE
]

# Remove heavy model columns if present
drop_cols <- intersect(
  c("model", "model_summary"),
  colnames(cond_effect)
)

if (length(drop_cols) > 0) {
  cond_effect <- cond_effect[
    ,
    setdiff(colnames(cond_effect), drop_cols),
    drop = FALSE
  ]
}

fwrite(
  as.data.table(cond_effect),
  file.path(
    srcdir,
    "Figure19_SMAD3_condition_effect.csv"
  )
)

# -------------------------------------------------------------------------
# Canonical published statistical results
#
# The original Moran analysis was performed genome-wide and SMAD3 was then
# extracted. Running graph_test() on SMAD3 alone reproduces the raw p-value
# and Moran's I, but not the genome-wide multiple-testing-adjusted q-value.
# Therefore, use the original processed result for the published statistic.
# -------------------------------------------------------------------------
orig_moran_file <- file.path(
  repo_root,
  "processed_results",
  "04_trajectory",
  "smad3_moran_test.csv"
)

orig_condition_file <- file.path(
  repo_root,
  "processed_results",
  "04_trajectory",
  "condition_effect_SMAD3_trajectory.csv"
)

if (!file.exists(orig_moran_file)) {
  stop("Original SMAD3 Moran result not found: ", orig_moran_file)
}

if (!file.exists(orig_condition_file)) {
  stop("Original SMAD3 condition-effect result not found: ", orig_condition_file)
}

orig_moran <- fread(orig_moran_file)
orig_condition <- fread(orig_condition_file)

# Remove row-name columns introduced by write.csv(), if present.
if (names(orig_moran)[1] %in% c("", "V1")) {
  orig_moran[, (names(orig_moran)[1]) := NULL]
}

if (names(orig_condition)[1] %in% c("", "V1")) {
  orig_condition[, (names(orig_condition)[1]) := NULL]
}

# Validate that the reconstructed single-gene test reproduces the raw
# Moran statistic and raw p-value from the original genome-wide analysis.
if (!isTRUE(all.equal(
  as.numeric(smad3_gt_dt$morans_I[1]),
  as.numeric(orig_moran$morans_I[1]),
  tolerance = 1e-10
))) {
  stop("Reconstructed Moran's I does not match the original result.")
}

if (!isTRUE(all.equal(
  as.numeric(smad3_gt_dt$p_value[1]),
  as.numeric(orig_moran$p_value[1]),
  tolerance = 1e-10
))) {
  stop("Reconstructed Moran raw p-value does not match the original result.")
}

fwrite(
  orig_moran,
  file.path(
    srcdir,
    "Figure19_SMAD3_Moran_test.csv"
  )
)

fwrite(
  orig_condition,
  file.path(
    srcdir,
    "Figure19_SMAD3_condition_effect.csv"
  )
)

cat(
  "\nPublished SMAD3 Moran q-value:",
  orig_moran$q_value[1],
  "\n"
)

cat(
  "Published SMAD3 condition-effect estimate:",
  orig_condition$estimate[1],
  "\n"
)

# -------------------------------------------------------------------------
# Validation summary
# -------------------------------------------------------------------------
cat("\n===== Figure 19 source-data summary =====\n")
cat("Figure19A-B Young cells:", nrow(fig19ab), "\n")
cat("Figure19C-D Young+Aged cells:", nrow(fig19cd), "\n")

cat("\nYoung-only cells by FAP:\n")
print(fig19ab[, .N, by = FAP_cluster])

cat("\nCombined cells by condition and FAP:\n")
print(fig19cd[, .N, by = .(condition, FAP_cluster)])

cat("\nReconstructed single-gene SMAD3 Moran test (validation only):\n")
print(smad3_gt_dt)

cat("\nPublished genome-wide-adjusted SMAD3 Moran result:\n")
print(orig_moran)

cat("\nPublished SMAD3 condition effect:\n")
print(orig_condition)

cat("\nDone.\n")
