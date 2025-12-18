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

# ===== Build a Monocle3 cell_data_set =====
# Extract expression matrix from Seurat (SCT assay, log-normalized data)
expr_mat <- GetAssayData(obj, assay = "SCT", slot = "data")
expr_mat <- as(expr_mat, "dgCMatrix")  # monocle3 expects a sparse matrix

# Cell metadata from Seurat
cell_meta <- obj@meta.data

# Gene metadata
gene_anno <- data.frame(
  gene_short_name = rownames(expr_mat),
  row.names       = rownames(expr_mat)
)

# Create Monocle3 cell_data_set
cds <- new_cell_data_set(
  expression_data = expr_mat,
  cell_metadata   = cell_meta,
  gene_metadata   = gene_anno
)

# Add FAP cluster label and condition explicitly (for clarity)
colData(cds)$cluster_fap <- cds$skeletal_muscle   
colData(cds)$condition   <- cds$condition         

# Keep only FAP1–FAP3 
keep_fap <- colnames(cds)[
  colData(cds)$cluster_fap %in% c("FAP1", "FAP2", "FAP3")
]
cds_fap <- cds[, keep_fap]

# Extract WNN UMAP coordinates from Seurat
umap_mat <- obj@reductions[["wnn.umap"]]@cell.embeddings

# Match UMAP rows to the cells in cds_fap
umap_mat <- umap_mat[colnames(cds_fap), ]

# Insert Seurat UMAP into Monocle reducedDims
reducedDims(cds_fap)$UMAP <- umap_mat

# ===== TRAJECTORY on Young FAP1–FAP3: test if FAP3 is intermediate =====
# Subset Monocle CDS to Young + FAP1-3
keep_young <- colnames(cds_fap)[
  colData(cds_fap)$condition   == "Young" &
    colData(cds_fap)$cluster_fap %in% c("FAP1", "FAP2", "FAP3")
]

cds_young <- cds_fap[, keep_young]

# Quick sanity check
table(colData(cds_young)$cluster_fap, colData(cds_young)$condition)

# Cluster cells on the UMAP embedding
cds_young <- cluster_cells(
  cds_young,
  reduction_method = "UMAP"
)

# Learn principal graph on UMAP (single partition)
cds_young <- learn_graph(
  cds_young,
  use_partition = FALSE
)

# Define root cells: FAP1 in Young as early state
root_cells <- colnames(cds_young)[
  colData(cds_young)$cluster_fap == "FAP1"
]

length(root_cells)  # safety check

# Order cells in pseudotime
cds_young <- order_cells(
  cds_young,
  root_cells = root_cells
)

# Extract pseudotime and store it in colData
cds_young$pseudotime <- monocle3::pseudotime(cds_young)

# Visualization of trajectory and pseudotime
# Trajectory colored by FAP cluster (FAP1, FAP2, FAP3)
p_traj_clusters <- plot_cells(
  cds_young,
  reduction_method = "UMAP",
  color_cells_by = "cluster_fap",
  label_groups_by_cluster = FALSE,
  label_leaves = TRUE,
  label_branch_points = TRUE
)
save_gg(p_traj_clusters,"monocle_FAP1_3_clusters.png",w=8,h=5)

# Trajectory colored by pseudotime
p_traj_pt <- plot_cells(
  cds_young,
  reduction_method = "UMAP",
  color_cells_by = "pseudotime",
  label_groups_by_cluster = FALSE,
  label_leaves = TRUE,
  label_branch_points = TRUE
)
save_gg(p_traj_pt,"monocle_FAP1_3_pseudotime.png",w=8,h=5)

# Pseudotime distribution per FAP cluster (FAP1, FAP3, FAP2)
pt_df <- as.data.frame(colData(cds_young)) %>%
  dplyr::select(cluster_fap, pseudotime) %>%
  dplyr::filter(!is.na(pseudotime)) %>%
  dplyr::group_by(cluster_fap) %>%
  dplyr::summarise(
    n_cells = dplyr::n(),
    median_pt = median(pseudotime),
    q1_pt = quantile(pseudotime, 0.25),
    q3_pt = quantile(pseudotime, 0.75)
  ) %>%
  dplyr::arrange(median_pt)

wcsv(pt_df, "../outputs/sc/pseudotime_distribution.csv")
# Interpretation:
# median_pt(FAP1) < median_pt(FAP3) < median_pt(FAP2)
# -> FAP3 is a transcriptional intermediate between FAP1 and FAP2.

# Boxplot of pseudotime per FAP cluster
pt_cells_df <- as.data.frame(colData(cds_young))

save_gg(ggplot(pt_cells_df %>% dplyr::filter(!is.na(pseudotime)),
               aes(x = cluster_fap, y = pseudotime)) +
          geom_boxplot(outlier.size = 0.3) +
          theme_bw() +
          labs(
            title = "Pseudotime distribution across FAP1–FAP3 (Young)",
            x = "FAP cluster",
            y = "Pseudotime"
          ),"monocle_FAP1_3_pseudotime_distribution.png",w=8,h=5)

# Trajectory program between Young vs Aged
# Subset to FAP1–3 in both Young and Aged
keep_young_aged <- colnames(cds_fap)[
  colData(cds_fap)$cluster_fap %in% c("FAP1", "FAP2", "FAP3") &
    colData(cds_fap)$condition   %in% c("Young", "Aged")
]

cds_ya <- cds_fap[, keep_young_aged]

# Re-use same WNN UMAP 
# Cluster + learn graph
cds_ya <- cluster_cells(cds_ya, reduction_method = "UMAP")
cds_ya <- learn_graph(cds_ya, use_partition = FALSE)

# Root = FAP1-Young cells (early, baseline)
root_cells_ya <- colnames(cds_ya)[
  colData(cds_ya)$cluster_fap == "FAP1" &
    colData(cds_ya)$condition == "Young"
]

cds_ya <- order_cells(cds_ya, root_cells = root_cells_ya)
cds_ya$pseudotime <- monocle3::pseudotime(cds_ya)

# Visualize trajectory colored by condition (Young vs Aged)
save_gg(plot_cells(
  cds_ya,
  reduction_method = "UMAP",
  color_cells_by = "condition",
  label_groups_by_cluster = FALSE,
  label_leaves = TRUE,
  label_branch_points = TRUE,
  graph_label_size = 1
),"monocle_trajectory_condition.png",w=8,h=5)

# Subset to FAP1–FAP4 (Young only for clarity)
keep_all <- colnames(cds)[
  colData(cds)$cluster_fap %in% c("FAP1", "FAP2", "FAP3", "FAP4") &
    colData(cds)$condition == "Young"
]

cds_all <- cds[, keep_all]

# Attach full WNN UMAP from Seurat 
umap_full <- obj@reductions[["wnn.umap"]]@cell.embeddings
umap_full <- umap_full[colnames(cds_all), ]   # match cells

reducedDims(cds_all)$UMAP <- umap_full

# Cluster + learn graph on all FAPs
cds_all <- cluster_cells(cds_all, reduction_method = "UMAP")

cds_all <- learn_graph(
  cds_all,
  use_partition = FALSE
)

# Root = FAP1 (Young)
root_all <- colnames(cds_all)[
  colData(cds_all)$cluster_fap == "FAP1"
]

cds_all <- order_cells(
  cds_all,
  root_cells = root_all
)

cds_all$pseudotime <- monocle3::pseudotime(cds_all)

# Plot trajectory colored by FAP clusters
p_all <- plot_cells(
  cds_all,
  reduction_method = "UMAP",
  color_cells_by = "cluster_fap",
  label_groups_by_cluster = FALSE,
  label_leaves = TRUE,
  label_branch_points = TRUE,
  cell_size = 0.4,
  graph_label_size = 1.5
)

save_gg(p_all,"monocle_FAP1_4_clusters.png",w=8,h=5)

# Plot pseudotime for FAP1–4
p_pt_all <- plot_cells(
  cds_all,
  reduction_method = "UMAP",
  color_cells_by = "pseudotime",
  label_groups_by_cluster = FALSE,
  label_leaves = TRUE,
  label_branch_points = TRUE,
  cell_size = 0.4,
  graph_label_size = 1.5
)

save_gg(p_pt_all,"monocle_FAP1_4_pseudotime.png",w=8,h=5)
