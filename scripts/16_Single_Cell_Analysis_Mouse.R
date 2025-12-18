# ========================
# Setting up environment
# ========================
# ===== Clean environment =====
rm(list = ls(all.names = TRUE))
gc()

# ===== Loading relevant libraries =====
library(decontX)
library(patchwork)
library(scuttle) 
library(scDblFinder)
library(BiocParallel)
library(scater)
library(glmGamPoi) 
library(presto) 
library(Rsamtools)
library(data.table)
library(R.utils)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(readxl)
library(purrr)
library(stringr)
library(edgeR)
library(limma)
library(biovizBase)
library(Seurat)
library(tidyverse)
library(Signac)
library(clustree)
library(harmony)
library(cluster)

# ===== set seed =====
set.seed(1234)

# ===== Helpers =====
figdir <- "../outputs/mouse"
save_gg <- function(p, filename, w=7, h=5, dpi=300){
  if (inherits(p, "ggplot") || inherits(p, "patchwork")){
    ggsave(file.path(figdir, filename), p, width=w, height=h, units="in", dpi=dpi, bg="white")
  } else {
    png(file.path(figdir, filename), width=w, height=h, units="in", res=dpi)
    print(p); dev.off()
  }
}
save_dev <- function(filename, expr, w=7, h=5, dpi=300){
  png(file.path(figdir, filename), width=w, height=h, units="in", res=dpi)
  on.exit(dev.off(), add=TRUE); force(expr)
}
wcsv <- function(x, path) write.csv(x, path, row.names = TRUE)
wtxt <- function(v, path) write.table(v, path, quote = FALSE, row.names = FALSE, col.names = FALSE)

# ============================================================================ #

# ========================
# Load the object
# ========================
obj <- readRDS("../data/GSE288662/Processed_Seurat_Object/GSE288662_Processed_Seurat_Object.rds")

# ========================
# Quality Control
# ========================

# ===== Post-filter ATAC visualization =====
DefaultAssay(obj) <- "ATAC"
metadata <- obj@meta.data

cnames <- setNames(rep(c("cyan3","darkgoldenrod1"), each = 5), levels(factor(metadata$orig.ident)))

# Plot (ATAC metrics)
save_gg(
  VlnPlot(
    obj,
    features = c("nCount_ATAC","TSS.enrichment","nucleosome_signal"),
    group.by = "orig.ident", raster=FALSE, alpha=0.07, ncol = 3
  ) + scale_fill_manual(values = cnames),
  "VlnPlot_ATAC_metrics_post_filt.png", w = 12, h = 5)

# ===== Post-filter RNA visualization =====
DefaultAssay(obj) <- "RNA"

# set colors
cnames <-setNames(rep(c("cyan3","darkgoldenrod1"),each=5),levels(factor(metadata$orig.ident))) 

# plot (RNA metrics)
save_gg(
  VlnPlot(
    obj,
    features = c("nCount_RNA","nFeature_RNA","percent.mt"),
    layer="counts",
    group.by = "orig.ident", raster=FALSE, alpha=0.07, ncol = 3
  ) + scale_fill_manual(values = cnames),
  "VlnPlot_RNA_metrics_post_filt.png", w = 12, h = 5)

# Plot (RNA density plot)
p14 <- ggplot(metadata, aes(color=orig.ident, x=nCount_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() 
p15 <- ggplot(metadata, aes(color=orig.ident, x=nFeature_RNA,fill=orig.ident)) +
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() 
p16 <- ggplot(metadata, aes(color=orig.ident, x=percent.mt,fill=orig.ident)) +
  geom_density(alpha = 0.2) + 
  scale_x_log10()+
  theme_classic()
p17 <- wrap_plots(list(p14,p15,p16), ncol = 2, guides = "collect") & 
  theme(legend.position = "bottom")
save_gg(p17, "DensityPlot_RNA_Metrics_post_filt.png", w=8, h=5)

# Examine these metrics together
p18 <- FeatureScatter(obj, feature1 = "percent.mt", feature2 ="nFeature_RNA", 
                     group.by="orig.ident")
p19 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                     group.by="orig.ident",log=TRUE)
p20 <- wrap_plots(list(p18,p19), ncol = 2, guides = "collect") & 
  theme(legend.position = "right")
save_gg(p20, "FeatureScatter_RNA_Metrics_post_filt.png", w=8, h=5)

# ============================================================================ #

# ========================
# Load the object
# ========================
obj <- readRDS("../data/GSE288662/Processed_Seurat_Object/GSE288662_Processed_Seurat_Object.rds")

DefaultAssay(obj) <- "RNA"

obj_rna <- DietSeurat(
  obj,
  assays = "RNA",
  dimreducs = NULL,
  graphs = NULL
)

# ========================
# Explanatory Variables
# ========================
# ===== plotExplanatoryVariables =====
# Convert to SCE and create raw logcounts from RNA counts
obj_rna[["RNA"]] <- JoinLayers(obj_rna[["RNA"]])

sce <- as.SingleCellExperiment(obj_rna)
assay(sce, "logcounts_raw") <- log1p(counts(sce)) 

# Pull QC metrics directly from metadata
qc_vars <- c(
  "nFeature_RNA","nCount_RNA","percent.mt","orig.ident")
qc_vars <- qc_vars[ qc_vars %in% colnames(obj_rna@meta.data) ]

# Copy columns to colData
for (v in qc_vars) colData(sce)[[v]] <- obj_rna@meta.data[[v]]

# Choose which variables to visualize/assess
vars <- c("nFeature_RNA","nCount_RNA","percent.mt","orig.ident")
vars <- vars[ vars %in% colnames(colData(sce)) ]

# Density of gene-wise marginal R^2 for explanatory variables
# Curves further to the right => variable explains more variance across genes
save_gg(plotExplanatoryVariables(sce, exprs_values = "logcounts_raw", 
                                 variables    = vars), 
        "plotExplanatoryVariables.png", w=8, h=5)

# ============================================================================ #

# ========================
# Load the object
# ========================
obj <- readRDS("../data/GSE288662/Processed_Seurat_Object/GSE288662_Processed_Seurat_Object.rds")

# ========================
# Annotation
# ========================

# Plot with annotation
my_colors_21 <- c(
  "#FFD92F", # yellow
  "#377EB8", # blue
  "#E78AC3", # magenta
  "#7570B3", # muted violet
  "#984EA3", # purple
  "#E41A1C", # red
  "#F781BF", # pink
  "#A65628", # brown
  "#66C2A5", # teal
  "#FC8D62", # salmon / coral
  "#8DA0CB", # soft indigo
  "#4DAF4A", # green
  "#B3B3B3", # light gray
  "#A6D854", # lime green
  "#999999", # gray
  "#E5C494", # tan
  "#1B9E77", # deep teal/green
  "#D95F02", # burnt orange
  "#E7298A", # hot pink
  "#66A61E"  # olive/green
)
clust_col <- "wsnn_res.0.8"
clust_lvls <- levels(factor(obj[[clust_col]][,1]))
plot_colors <- my_colors_21[seq_along(clust_lvls)]
f1 <- DimPlot(
  obj,
  reduction = "wnn.umap",
  group.by  = clust_col,
  alpha = 0.4,
  label = T,
  label.size = 2
) + NoLegend() + scale_color_manual(values = plot_colors)
f2 <- DimPlot(
  obj,
  reduction = "wnn.umap",
  group.by  = "type",
  alpha = 0.4,
  label = T,
  label.size = 1.5
) + NoLegend() + scale_color_manual(values = plot_colors)
f <- (f1|f2)
save_gg(f, "Clusters_with_annotation.png", w=11, h=5)

# ============================================================================ #

# ========================
# Load the object
# ========================
obj <- readRDS("../data/GSE288662/Processed_Seurat_Object/GSE288662_Processed_Seurat_Object.rds")
# obj1 <- obj %>% subset(type == 'FAPs')

# ========================
# Differential expression analysis (pseudobulk with voom+limma)
# ========================

DefaultAssay(obj) <- "RNA"

# Drop ATAC and SCT assays, keep all embeddings/graphs
obj[["ATAC"]] <- NULL
obj[["SCT"]] <- NULL

# Join RNA layers to get a single-layer assay
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])

# AggregateExpression:
#   - genes x groups matrix
#   - group.by = sample_condition_celltype (orig.ident, condition, type)
pb_list <- AggregateExpression(
  obj,
  assays = "RNA",
  group.by = c("orig.ident", "condition", "type"),
  slot = "counts",
  return.seurat = FALSE
)

pb_counts <- pb_list$RNA   # genes x groups matrix

# ========================
# Build group metadata from column names
# ========================
# Column names assumed: sample_condition_celltype (with "_" separator)
sp <- stringr::str_split_fixed(colnames(pb_counts), "_", 3)

meta <- data.frame(
  sample = sp[, 1],
  condition = sp[, 2], # 3-group condition (e.g., Veh, EP, EPR)
  celltype = sp[, 3],
  stringsAsFactors = FALSE
)
rownames(meta) <- colnames(pb_counts)

# List of unique cell types
subtypes <- sort(unique(meta$celltype))

# ========================
# voom+limma helper: all pairwise contrasts
# ========================
run_voomlimma <- function(L) {
  message("voomlimma")
  
  # L$count = counts matrix (genes x samples)
  # L$condt = factor of condition (length = ncol(counts))
  stopifnot(ncol(L$count) == length(L$condt))
  
  dge <- DGEList(L$count, group = L$condt)
  dge <- calcNormFactors(dge)
  
  # Design without intercept; one column per condition level
  design <- model.matrix(~ 0 + L$condt)
  colnames(design) <- levels(L$condt)
  
  vm <- voom(dge, design = design, plot = FALSE)
  fit <- lmFit(vm, design = design)
  
  # Build all pairwise contrasts between condition levels
  levs <- colnames(design)  # e.g., c("Veh","EP","EPR")
  if (length(levs) < 2) {
    stop("Need at least two condition levels for contrasts.")
  }
  
  # Generate expressions like "EP - Veh", "EPR - Veh", ...
  combs <- combn(levs, 2)   # matrix with 2 x nContrasts
  contrast_expr <- apply(combs, 2, function(z) {
    paste0(z[2], " - ", z[1])
  })
  # Names like "EP_vs_Veh"
  contrast_names <- apply(combs, 2, function(z) {
    paste0(z[2], "_vs_", z[1])
  })
  
  # Build contrast matrix
  # do.call + makeContrasts allows dynamic number of contrasts
  cargs <- as.list(contrast_expr)
  names(cargs) <- contrast_names
  cargs$levels <- design
  
  contr_mat <- do.call(makeContrasts, cargs)
  
  fit2 <- contrasts.fit(fit, contr_mat)
  fit2 <- eBayes(fit2)
  
  # Extract topTable for each contrast
  res_list <- lapply(colnames(contr_mat), function(cn) {
    tt <- topTable(fit2, coef = cn, number = Inf, adjust.method = "BH")
    tt$contrast <- cn
    tt
  })
  names(res_list) <- colnames(contr_mat)
  
  return(res_list)
}

# ========================
# Run voom+limma per cell type (all pairwise contrasts)
# ========================
res_list <- lapply(subtypes, function(st) {
  message("Running voom+limma for celltype: ", st)
  
  # Subset columns for this cell type
  keep_cols <- meta$celltype == st
  counts_ct <- pb_counts[, keep_cols, drop = FALSE]
  meta_ct <- meta[keep_cols, , drop = FALSE]
  
  # Build condition factor
  cond_ct <- factor(meta_ct$condition)
  
  # Drop NAs and empty levels
  valid <- !is.na(cond_ct)
  counts_ct <- counts_ct[, valid, drop = FALSE]
  cond_ct <- droplevels(cond_ct[valid])
  
  # Quick sanity check on samples per condition
  tab_cond <- table(cond_ct)
  if (length(tab_cond) < 2) {
    warning("Skipping ", st, ": <2 condition levels.")
    return(NULL)
  }
  if (any(tab_cond < 2)) {
    warning("Skipping ", st, ": some conditions have <2 replicates.")
    return(NULL)
  }
  
  # Filter lowly expressed genes (edgeR recommended)
  dge_tmp <- DGEList(counts_ct)
  keep_genes <- filterByExpr(dge_tmp, group = cond_ct, min.count = 10)
  counts_ct <- counts_ct[keep_genes, , drop = FALSE]
  
  if (nrow(counts_ct) == 0) {
    warning("No genes left after filtering for ", st)
    return(NULL)
  }
  
  # Run voom+limma with all pairwise contrasts
  L <- list(count = counts_ct, condt = cond_ct)
  out <- run_voomlimma(L)  # list of topTable's, one per contrast
  
  # Convert each contrast result to a tidy data.frame
  res_ct <- lapply(names(out), function(cn) {
    tt <- out[[cn]]
    tt$gene <- rownames(tt)
    tt$celltype <- st
    
    tt %>%
      transmute(
        gene = gene,
        celltype = celltype,
        contrast = contrast,           # e.g. "EP_vs_Veh"
        avg_log2FC = round(logFC,     2),
        p_val = round(P.Value,   4),
        p_val_adj = round(adj.P.Val, 4)
      )
  })
  
  res_ct <- bind_rows(res_ct)
  return(res_ct)
})

# Drop NULL (celltypes with insufficient samples)
res_list <- Filter(Negate(is.null), res_list)

# ========================
# Combine across cell types and rank
# ========================
bulk_de <- bind_rows(res_list) %>%
  arrange(p_val_adj, p_val)

# Significant DE genes (change thresholds if you like)
bulk_de_sig <- bulk_de %>%
  filter(p_val_adj < 0.05 & abs(avg_log2FC) >= 1)

# ========================
# Save results
# ========================
# All DE (all contrasts, all cell types)
wcsv(bulk_de, "../outputs/bulk_de_all_contrasts_mouse.csv")

# Significant subset
wcsv(bulk_de_sig, "../outputs/bulk_de_sig_all_contrasts_mouse.csv")
