# ========================
# Setting up environment
# ========================
# ===== Clean environment =====
rm(list = ls(all.names = TRUE))
gc()

# ===== Loading relevant libraries =====
library(monocle3)
library(Seurat)
library(ggplot2)
library(dplyr)

# ===== set seed =====
set.seed(1234)

# ===== Helpers =====
figdir <- "../outputs/sc/figs"
save_gg <- function(p, filename, w=7, h=5, dpi=300){
  if (inherits(p, "ggplot") || inherits(p, "patchwork")){
    ggsave(file.path(figdir, filename), p, width=w, height=h, units="in", dpi=dpi, bg="white")
  } else {
    png(file.path(figdir, filename), width=w, height=h, units="in", res=dpi)
    print(p); dev.off()
  }
}

wcsv <- function(x, path) write.csv(x, path, row.names = TRUE)

# ============================================================================ #

# ========================
# Load the object
# ========================
obj <- readRDS("../outputs/hdWGCNA_obj.rds")

# Build Monocle3 CDS from Seurat
# Expression matrix from SCT assay 
expr_mat <- GetAssayData(obj, assay = "SCT", slot = "data")
expr_mat <- as(expr_mat, "dgCMatrix")

# Cell metadata and gene metadata
cell_meta <- obj@meta.data
gene_anno <- data.frame(
  gene_short_name = rownames(expr_mat),
  row.names = rownames(expr_mat)
)

# Create cell_data_set
cds <- new_cell_data_set(
  expression_data = expr_mat,
  cell_metadata = cell_meta,
  gene_metadata = gene_anno
)

# Explicitly copy FAP cluster and condition
colData(cds)$cluster_fap <- cds$skeletal_muscle 
colData(cds)$condition <- cds$condition 


# Subset to FAP1-FAP3 and attach WNN UMAP
keep_fap <- colnames(cds)[
  colData(cds)$cluster_fap %in% c("FAP1", "FAP2", "FAP3")
]
cds_fap <- cds[, keep_fap]

# WNN UMAP from Seurat
wnn_umap <- obj@reductions[["wnn.umap"]]@cell.embeddings
wnn_umap <- wnn_umap[colnames(cds_fap), ]
reducedDims(cds_fap)$UMAP <- wnn_umap

# YOUNG-ONLY TRAJECTORY (FAP1–FAP3) FOR SMAD3
# Subset to Young FAP1–FAP3
keep_young <- colnames(cds_fap)[
  colData(cds_fap)$condition   == "Young" &
    colData(cds_fap)$cluster_fap %in% c("FAP1", "FAP2", "FAP3")
]
cds_young <- cds_fap[, keep_young]

# Cluster cells on UMAP and learn graph
cds_young <- cluster_cells(cds_young, reduction_method = "UMAP")
cds_young <- learn_graph(cds_young, use_partition = FALSE)

# Root cells = FAP1 (Young)
root_cells_young <- colnames(cds_young)[
  colData(cds_young)$cluster_fap == "FAP1"
]
cds_young <- order_cells(cds_young, root_cells = root_cells_young)

# Add pseudotime column to colData
cds_young$pseudotime <- monocle3::pseudotime(cds_young)

# SMAD3 expression along pseudotime (Young only)
# Check that SMAD3 exists
"SMAD3" %in% rownames(cds_young)

# Plot SMAD3 vs pseudotime, colored by FAP cluster
p_smad3_pt_young <- plot_genes_in_pseudotime(
  cds_young["SMAD3", ],
  color_cells_by = "cluster_fap"
)
save_gg(p_smad3_pt_young,
"monocle_FAP1_3_SMAD3_cluster.png",w=8,h=5)

# Same plot, colored by pseudotime (for visual smoothness)
p_smad3_pt_young_pt <- plot_genes_in_pseudotime(
  cds_young["SMAD3", ],
  color_cells_by = "pseudotime"
)

save_gg(p_smad3_pt_young_pt,
        "monocle_FAP1_3_SMAD3_pseudotime.png",w=8,h=5)

# Is SMAD3 trajectory-associated? (graph_test)
gt_young <- graph_test(
  cds_young,
  neighbor_graph = "principal_graph",
  cores = 4
)

smad3_gt_young <- gt_young[gt_young$gene_short_name == "SMAD3", ]
wcsv(smad3_gt_young, "../outputs/sc/smad3_moran_test.csv")
# Look at morans_I and q_value:
# q_value < 0.05 => SMAD3 significantly varies along the trajectory.

# Combined Young + Aged trajectory (FAP1–FAP3)
# Subset to FAP1-FAP3 in both conditions
keep_ya <- colnames(cds_fap)[
  colData(cds_fap)$cluster_fap %in% c("FAP1", "FAP2", "FAP3") &
    colData(cds_fap)$condition   %in% c("Young", "Aged")
]
cds_ya <- cds_fap[, keep_ya]

# Ensure UMAP is present for this subset
reducedDims(cds_ya)$UMAP <- wnn_umap[colnames(cds_ya), ]

# Cluster cells and learn shared graph
cds_ya <- cluster_cells(cds_ya, reduction_method = "UMAP")
cds_ya <- learn_graph(cds_ya, use_partition = FALSE)

# Root cells = FAP1 in Young
root_cells_ya <- colnames(cds_ya)[
  colData(cds_ya)$cluster_fap == "FAP1" &
    colData(cds_ya)$condition   == "Young"
]
cds_ya <- order_cells(cds_ya, root_cells = root_cells_ya)

# Add pseudotime column
cds_ya$pseudotime <- monocle3::pseudotime(cds_ya)

# SMAD3 along pseudotime, colored by condition (Young vs Aged)
# Pseudotime curves by condition
p_smad3_cond <- plot_genes_in_pseudotime(
  cds_ya["SMAD3", ],
  color_cells_by = "condition"
)

save_gg(p_smad3_cond,
        "monocle_FAP1_3_SMAD3_pseudotime_condition.png",w=8,h=5)

# Smooth expression vs pseudotime, faceted by condition
smad3_expr <- as.numeric(exprs(cds_ya)["SMAD3", ])
pt_df_ya <- as.data.frame(colData(cds_ya)) %>%
  mutate(SMAD3_expr = smad3_expr)

save_gg(ggplot(pt_df_ya, aes(x = pseudotime, y = SMAD3_expr, color = condition)) +
  geom_point(alpha = 0.4, size = 0.5) +
  geom_smooth(se = FALSE, method = "loess", span = 0.5) +
  theme_bw() +
  labs(
    title = "SMAD3 expression along pseudotime (FAP1–FAP3, Young vs Aged)",
    x = "Pseudotime",
    y = "SMAD3 expression"
  ),"monocle_FAP1_3_SMAD3_pseudotime_condition_alt.png",w=8,h=5)

# Test condition effect on SMAD3 along the trajectory (fit_models)
# Fit a model with condition + smooth pseudotime
# We use a spline on pseudotime to capture non-linear trends.
smad3_fit <- fit_models(
  cds_ya["SMAD3", ],
  model_formula_str = "~ condition + splines::ns(pseudotime, df = 3)"
)

coef_tab <- coefficient_table(smad3_fit)
coef_tab

# Extract the condition effect row (Aged vs Young)
# The exact term name depends on factor coding; here we search for "condition".
cond_effect <- coef_tab %>%
  dplyr::filter(grepl("condition", term))

cond_effect_clean <- cond_effect |>
  dplyr::select(-model, -model_summary)

wcsv(cond_effect_clean, "../outputs/sc/condition_effect_SMAD3_trajectory.csv")
# If q_value for condition term < 0.05:
# => SMAD3 is significantly different between Young and Aged along the trajectory.
