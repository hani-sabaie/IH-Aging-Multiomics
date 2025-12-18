# ========================
# Setting up environment
# ========================
# ===== Clean environment =====
rm(list = ls(all.names = TRUE))
gc()

# ===== Loading relevant libraries =====
library(WGCNA)
library(UCell)
library(hdWGCNA)
library(enrichR)
library(tidyverse)
library(cowplot)
library(patchwork)
library(writexl)
library(Seurat)
library(Signac)

# ===== set seed =====
set.seed(1234)

# ===== Helpers =====
resdir <- "../outputs/results_hdWGCNA"
save_gg <- function(p, filename, w=7, h=5, dpi=300){
  if (inherits(p, "ggplot") || inherits(p, "patchwork")){
    ggsave(file.path(resdir, filename), p, width=w, height=h, units="in", dpi=dpi, bg="white")
  } else {
    png(file.path(resdir, filename), width=w, height=h, units="in", res=dpi)
    print(p); dev.off()
  }
}
save_dev <- function(filename, expr, w=7, h=5, dpi=300){
  png(file.path(resdir, filename), width=w, height=h, units="in", res=dpi)
  on.exit(dev.off(), add=TRUE); force(expr)
}
wcsv <- function(x, path) write.csv(x, path, row.names = TRUE)
wtxt <- function(v, path) write.table(v, path, quote = FALSE, row.names = FALSE, col.names = FALSE)

# Using the cowplot theme for ggplot
theme_set(theme_cowplot())

# Optionally enable multithreading
enableWGCNAThreads(nThreads = 8)

# ============================================================================ #

# ========================
# Load the object
# ========================
obj <- readRDS("../outputs/decont_merged_filt_nodoub_cc_sct_reduc_clust_integ_annot_obj.rds")

# Subset just one cell type 
obj <- obj %>% subset(cell_type == 'FAPs')
table(obj$skeletal_muscle)
with(obj@meta.data, table(sample, skeletal_muscle))

# ========================
# SCTransform on single-cell data
# ========================
# ===== Set up Seurat object for WGCNA =====
obj <- SetupForWGCNA(
  obj,
  features = VariableFeatures(obj, nfeatures = 5000),
  wgcna_name = "SCT")

# Construct the metacell
obj <- MetacellsByGroups(
  seurat_obj = obj,
  group.by = c("sample"),
  k = 25,
  max_shared=12,
  min_cells = 50,
  reduction = 'harmony',
  ident.group = 'sample',
  slot = 'scale.data',
  assay = 'SCT')

# ========================
# Co-expression network analysis
# ========================
# Set expression matrix for hdWGCNA
obj <- SetDatExpr(obj)

# Test different soft power thresholds
obj <- TestSoftPowers(obj, powers = 1:20)
plot_list <- PlotSoftPowers(obj)

# Assemble with patchwork
save_gg(wrap_plots(plot_list, ncol=2), 
        "PlotSoftPowers.png", w=12, h=9)

# Construct the co-expression network and identify gene modules
obj <- ConstructNetwork(
  obj, 
  tom_name='SCT_cells', 
  overwrite_tom=TRUE)

# Compute module eigengenes and connectivity
obj <- ModuleEigengenes(obj)
obj <- ModuleConnectivity(obj)

# Plot the dendrogram
save_dev("PlotDendrogram.png", {
  PlotDendrogram(obj, main = "FAPs hdWGCNA SCT Dendrogram")
}, w = 8, h = 5, dpi = 300)

# ========================
# Other figures and tables
# ========================
# ===== UMAP =====
# Compute the co-expression network umap 
obj <- RunModuleUMAP(
  obj,
  n_hubs = 10,
  n_neighbors=30,
  min_dist=0.1)

# Get the hub gene UMAP table from the Seurat object
umap_df <- GetModuleUMAP(obj)

# UMAP
p <- ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
    color=umap_df$color, # color each point by WGCNA module
    size=umap_df$kME*2 # size of each point based on intramodular connectivity
  ) +
  umap_theme()

save_gg(p, "module_UMAP.png", w=8, h=5) 

# ===== Supervised UMAP =====
# Compute the co-expression network umap 
obj <- RunModuleUMAP(
  obj,
  n_hubs = 10,
  n_neighbors=30,
  min_dist=0.1,
  supervised=TRUE,
  target_weight=0.5)

# Get the hub gene UMAP table from the Seurat object
umap_df <- GetModuleUMAP(obj)

# Plot with ggplot
p <- ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
    color=umap_df$color,
    size=umap_df$kME*2
  ) + 
  umap_theme() 

# Add the module names to the plot by taking the mean coordinates
centroid_df <- umap_df %>% 
  dplyr::group_by(module) %>%
  dplyr::summarise(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2))

p <- p + geom_label(
  data = centroid_df, 
  label=as.character(centroid_df$module), vjust = -1.5,
  fontface='bold', size=2.5) + 
  theme(panel.background = element_rect(fill='white'))

save_gg(p, "module_UMAP_supervised.png", w=8, h=5) 

# ===== PlotKMEs =====
# Plot genes ranked by kME for each module
p <- PlotKMEs(obj, ncol=3)
save_gg(p, "PlotKMEs.png", w=8, h=5)

# ===== Tables =====
# Get the module assignment table
modules <- GetModules(obj) %>% subset(module != 'grey')
wcsv(modules, file.path("results_hdWGCNA", "module_assignment_table.csv"))

# Get hub genes
hub_df <- GetHubGenes(obj, n_hubs = 10)
wcsv(hub_df, file.path("results_hdWGCNA", "hub_genes.csv"))

# ===== Hub gene signature scores =====
# Compute gene scoring for the top 25 hub genes by kME for each module
# with UCell method
obj <- ModuleExprScore(
  obj,
  n_genes = 25,
  method='UCell')

# ===== Featureplots =====
# Make a featureplot of MEs for each module
plot_list <- ModuleFeaturePlot(
  obj,
  features='hMEs', # plot the hMEs
  order=TRUE, # order so the points with highest hMEs are on top
  reduction = "wnn.umap")

save_gg(wrap_plots(plot_list, ncol=4), "ModuleFeaturePlot_hMEs.png", w=8, h=5)

# Make a featureplot of hub scores for each module
plot_list <- ModuleFeaturePlot(
  obj,
  features='scores', # plot the hub gene scores
  order='shuffle', # order so cells are shuffled
  ucell = TRUE, # depending on Seurat vs UCell for gene scoring
  reduction = "wnn.umap")

save_gg(wrap_plots(plot_list, ncol=4), "ModuleFeaturePlot_hub_scores.png", w=8, h=5)

# ===== subclusters RadarPlot =====
# FAPs subclusters
plot_list <- ModuleRadarPlot(
  obj,
  group.by = "skeletal_muscle",
  barcodes = obj@meta.data %>% subset(cell_type == 'FAPs') %>% rownames(),
  axis.label.size=2,
  grid.label.size=2, 
  combine = F) 
save_gg(wrap_plots(plot_list, ncol=4), "ModuleRadarPlot.png", w=8, h=5)

# ========================
# Save the object
# ========================
saveRDS(obj,"../outputs/hdWGCNA_obj.rds")

# ============================================================================ #

# ========================
# Load the object
# ========================
obj <- readRDS("../outputs/hdWGCNA_obj.rds")

# ========================
# Module Trait Correlation
# ========================
# Convert condition to factor
obj$condition <- factor(obj$condition, levels = c("Young","Aged"))

# list of traits to correlate
cur_traits <- c('condition', 'nCount_RNA', 'nFeature_RNA', 
                'nCount_ATAC','nFeature_ATAC','percent_mt','percent_ribo')

obj <- ModuleTraitCorrelation(
  obj,
  traits = cur_traits,
  group.by='skeletal_muscle')

# Get the mt-correlation results
mt_cor <- GetModuleTraitCorrelation(obj)
names(mt_cor)
names(mt_cor$cor)

combine_ct <- function(ct){
  cor <- mt_cor$cor[[ct]]; pv <- mt_cor$pval[[ct]]; fd <- mt_cor$fdr[[ct]]
  colnames(cor) <- paste0(colnames(cor), "_cor")
  colnames(pv)  <- paste0(colnames(pv),  "_pval")
  colnames(fd)  <- paste0(colnames(fd),  "_fdr")
  out <- cbind(trait = rownames(cor), cor, pv, fd)
  rownames(out) <- NULL
  as.data.frame(out, check.names = FALSE)
}
cts <- names(mt_cor$cor)
write_xlsx(setNames(lapply(cts, combine_ct), paste0(cts, "_wide")), "mt_cor_wide_all.xlsx")

# Plot Correlation Heatmap
save_gg(PlotModuleTraitCorrelation(
  obj,
  label = 'fdr',
  label_symbol = 'stars',
  text_size = 4,
  text_digits = 3,
  text_color = 'black',
  high_color = 'red',
  mid_color = 'grey90',
  low_color = 'blue',
  plot_max = 1,
  combine=TRUE), 
  "PlotModuleTraitCorrelation.png", w=9, h=12)

# ============================================================================ #

# ========================
# Load the object
# ========================
obj <- readRDS("../outputs/hdWGCNA_obj.rds")

# ========================
# Differential module eigengene (DME) analysis
# ========================
# ===== DME analysis comparing two groups =====
group1 <- obj@meta.data %>% subset(condition == "Aged") %>% rownames
group2 <- obj@meta.data %>% subset(condition == "Young") %>% rownames
head(group1)

DMEs <- FindDMEs(
  obj,
  features = 'ModuleScores',
  barcodes1 = group1,
  barcodes2 = group2,
  test.use='wilcox',
  wgcna_name='SCT'
)
DMEs$Comparisons <- "Aged_vs_Young" 
head(DMEs)

save_gg(PlotDMEsLollipop(
  obj, 
  DMEs, 
  wgcna_name='SCT',
  group.by = "Comparisons",
  comparison = c("Aged_vs_Young"),
  pvalue = "p_val_adj"), 
  "PlotDMEsLollipop.png", w=8, h=5)

save_gg(PlotDMEsVolcano(
  obj,
  DMEs,
  wgcna_name = 'SCT'), 
  "PlotDMEsVolcano.png", w=8, h=5)

# ===== Looping through multiple clusters =====
# List of clusters to loop through
clusters <- c("FAP1","FAP2","FAP3","FAP4")

# Set up an empty dataframe for the DMEs
DMEs <- data.frame()

# Loop through the clusters
for(cur_cluster in clusters){
  
  # Identify barcodes for group1 and group2 in each cluster
  group1 <- obj@meta.data %>% subset(skeletal_muscle == cur_cluster & condition == "Aged") %>% rownames
  group2 <- obj@meta.data %>% subset(skeletal_muscle == cur_cluster & condition == "Young") %>% rownames
  
  # Run the DME test
  cur_DMEs <- FindDMEs(
    obj,
    barcodes1 = group1,
    barcodes2 = group2,
    test.use ='wilcox',
    pseudocount.use=0.01, # we can also change the pseudocount with this param
    wgcna_name = 'SCT'
  )
  
  # add the cluster info to the table
  cur_DMEs$cluster <- cur_cluster
  
  # append the table
  DMEs <- rbind(DMEs, cur_DMEs)
}

# Get the modules table
modules <- GetModules(obj)
mods <- levels(modules$module); mods <- mods[mods != 'grey']

# Make a copy of the DME table for plotting
plot_df <- DMEs

# Set the factor level for the modules so they plot in the right order:
plot_df$module <- factor(as.character(plot_df$module), levels=mods)

# Set a min/max threshold for plotting
maxval <- 0.5; minval <- -0.5
plot_df$avg_log2FC <- ifelse(plot_df$avg_log2FC > maxval, maxval, plot_df$avg_log2FC)
plot_df$avg_log2FC <- ifelse(plot_df$avg_log2FC < minval, minval, plot_df$avg_log2FC)

# Add significance levels
plot_df$Significance <- gtools::stars.pval(plot_df$p_val_adj)

# Change the text color to make it easier to see 
plot_df$textcolor <- ifelse(plot_df$avg_log2FC > 0.2, 'black', 'white')

# Make the heatmap with geom_tile
p <- plot_df %>% 
  ggplot(aes(y=cluster, x=module, fill=avg_log2FC)) +
  geom_tile() 

# Add the significance levels
p <- p + 
  geom_text(label=plot_df$Significance, color=plot_df$textcolor) 

# Customize the color and theme of the plot
p <- p + 
  scale_fill_gradient2(low='blue', mid='grey90', high='red') +
  RotatedAxis() +
  theme(
    panel.border = element_rect(fill=NA, color='black', size=1),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    plot.margin=margin(0,0,0,0)
  ) + xlab('') + ylab('') +
  coord_equal()

save_gg(p, 
  "heatmap_DME_effect_sizes.png", w=8, h=5)

# ========================
# Selected module gene network
# ========================
# Get modules and TOM from the seurat obj
col <- "brown"
modules <- GetModules(obj) %>% 
  subset(module != 'grey') %>% subset(module == col) %>%
  mutate(module = droplevels(module))
mods <- levels(modules$module)
TOM <- GetTOM(obj)

# Get genes in selected module
top_250 <- GetHubGenes(obj, n_hubs = 250)
cur_genes <- top_250[top_250$gene_name %in% modules$gene_name,]
cur_genes <- cur_genes[,c(-2,-3)]

# subset the TOM 
cur_TOM <- TOM[cur_genes,cur_genes] 

# Get module colors for plotting 
mod_colors <- dplyr::select(modules, c(module, color)) %>% distinct()
mod_cp <- mod_colors$color; names(mod_cp) <- as.character(mod_colors$module)

# Set up the graph object with igraph & tidygraph
graph <- cur_TOM %>% 
  igraph::graph_from_adjacency_matrix(mode='undirected', weighted=TRUE) %>% 
  tidygraph::as_tbl_graph(directed=FALSE) %>% 
  tidygraph::activate(nodes) 

# Make the plot with ggraph
p <- ggraph(graph) + 
  geom_edge_link(color= col, alpha=0.2) + 
  geom_node_point(color='black') +
  geom_node_label(aes(label=name), repel=TRUE, max.overlaps=Inf, fontface='italic') 

save_gg(p, "brown_module_gene_network.png", w=15, h=12)

# ============================================================================ #

# ========================
# Load the object
# ========================
obj <- readRDS("../outputs/hdWGCNA_obj.rds")

# ========================
# Enrichment analysis
# ========================
# ===== Enrichr =====
# Define the enrichr databases to test
dbs <- listEnrichrDbs()
dbs <- c('GO_Biological_Process_2025','GO_Cellular_Component_2025',
         'GO_Molecular_Function_2025','Reactome_Pathways_2024')

# Perform enrichment tests
obj <- RunEnrichr(
  obj,
  dbs=dbs,
  max_genes = Inf # use max_genes = Inf to choose all genes
  )

# Retrieve the output table
enrich_df <- GetEnrichrTable(obj)

# Look at the results
head(enrich_df)

# Enrichr dotplot
p1 <- EnrichrDotPlot(
  obj,
  mods = "all", # use all modules (default)
  database = "GO_Biological_Process_2025", # this must match one of the dbs used previously
  n_terms=5, # number of terms per module
  term_size=6, # font size for the terms
  p_adj = T # show the p-val or adjusted p-val?
)  + scale_color_stepsn(colors=rev(viridis::magma(256)))

p2 <- EnrichrDotPlot(
  obj,
  mods = "all", # use all modules (default)
  database = "GO_Cellular_Component_2025", # this must match one of the dbs used previously
  n_terms=5, # number of terms per module
  term_size=6, # font size for the terms
  p_adj = T # show the p-val or adjusted p-val?
)  + scale_color_stepsn(colors=rev(viridis::magma(256)))

p3 <- EnrichrDotPlot(
  obj,
  mods = "all", # use all modules (default)
  database = "GO_Molecular_Function_2025", # this must match one of the dbs used previously
  n_terms=5, # number of terms per module
  term_size=6, # font size for the terms
  p_adj = T # show the p-val or adjusted p-val?
)  + scale_color_stepsn(colors=rev(viridis::magma(256)))

p4 <- EnrichrDotPlot(
  obj,
  mods = "all", # use all modules (default)
  database = "Reactome_Pathways_2024", # this must match one of the dbs used previously
  n_terms=5, # number of terms per module
  term_size=6, # font size for the terms
  p_adj = T # show the p-val or adjusted p-val?
)  + scale_color_stepsn(colors=rev(viridis::magma(256)))

save_gg(wrap_plots(p1,p2,p3,p4, ncol = 2), 
        "EnrichrDotPlot.png", w=18, h=18)
