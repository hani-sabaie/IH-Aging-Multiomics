# ========================
# Setting up environment
# ========================
# ===== Clean environment =====
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage

# ===== Loading relevant libraries =====
library(WGCNA)
library(JASPAR2024)
library(motifmatchr)
library(TFBSTools)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(xgboost)
library(RSQLite)
library(enrichR)
library(igraph)
library(hdWGCNA)
library(tidyverse)
library(cowplot)
library(patchwork)
library(magrittr)
library(Seurat)
library(Signac)

# ===== set seed =====
set.seed(1234)

# ===== Helpers =====
figdir <- "../outputs/TF"
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
obj <- readRDS("../outputs/hdWGCNA_obj.rds")

# ========================
# Transcription Factor Network Analysis
# ========================
# ===== Identify TFs in promoter regions =====
# JASPAR 2024
JASPAR2024 <- JASPAR2024()
sq24 <- RSQLite::dbConnect(RSQLite::SQLite(), db(JASPAR2024))

pfm_core <- TFBSTools::getMatrixSet(
  x = sq24,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# Run the motif scan
obj <- MotifScan(
  obj,
  species_genome = 'hg38',
  pfm = pfm_core,
  EnsDb = EnsDb.Hsapiens.v86
)

# ===== Construct TF Regulatory Network =====
# Get the motif df
motif_df <- GetMotifs(obj)

# Keep all TFs, and then remove all genes from the grey module
tf_genes <- unique(motif_df$gene_name)
modules <- GetModules(obj)
nongrey_genes <- subset(modules, module != 'grey') %>% .$gene_name
genes_use <- c(tf_genes, nongrey_genes)

# Update the gene list and re-run SetDatExpr
obj <- SetWGCNAGenes(obj, genes_use)
obj <- SetDatExpr(obj)

# Define model params
model_params <- list(
  objective = 'reg:squarederror',
  max_depth = 1,
  eta = 0.1,
  nthread=16,
  alpha=0.5
)

# Construct the TF network
obj <- ConstructTFNetwork(obj, model_params)
tfnet <- GetTFNetwork(obj) # TF-gene with Gain/Cor/etc.

# Define TF Regulons
obj <- AssignTFRegulons(
  obj,
  strategy = "B", # Strategy “B” selects the top genes for each TF 
  reg_thresh = 0.01,
  n_genes = 50
)

tfrgs <- GetTFRegulons(obj) # filtered TF-gene pairs (regulons)

wcsv(tfnet, file.path("results_TF", "tf_network_full_B.csv"))
wcsv(tfrgs, file.path("results_TF", "tf_regulons_B.csv"))

# ===== Calculate regulon expression signatures =====
# Positive regulons
obj <- RegulonScores(
  obj,
  target_type = 'positive',
  ncores=8
)

# Negative regulons
obj <- RegulonScores(
  obj,
  target_type = 'negative',
  cor_thresh = -0.05,
  ncores=8
)

# Access the results
pos_regulon_scores <- GetRegulonScores(obj, target_type='positive')
neg_regulon_scores <- GetRegulonScores(obj, target_type='negative')

# ========================
# Save the object
# ========================
saveRDS(obj,"../outputs/hdWGCNA_TFNet_B_obj.rds")

# ============================================================================ #

# ========================
# Load the object
# ========================
obj <- readRDS("../outputs/hdWGCNA_TFNet_B_obj.rds")

# ========================
# RegulonBarPlot
# ========================
p1 <- RegulonBarPlot(obj, selected_tf='SMAD3', cutoff=0.10)
save_gg(p1, 
        "RegulonBarPlot_B.png", w=8, h=5) 

# ========================
# FeaturePlot
# ========================
# select a TF of interest
cur_tf <- 'SMAD3'

# add the regulon scores to the Seurat metadata
obj$pos_regulon_score <- pos_regulon_scores[,cur_tf]
obj$neg_regulon_score <- neg_regulon_scores[,cur_tf]

# plot using FeaturePlot
p1 <- FeaturePlot(obj, feature=cur_tf) + umap_theme()
p2 <- FeaturePlot(obj, feature='pos_regulon_score', cols=c('lightgrey', 'red')) + umap_theme()
p3 <- FeaturePlot(obj, feature='neg_regulon_score', cols=c('lightgrey', 'seagreen')) + umap_theme()

save_gg(p1 | p2 | p3, 
        "FeaturePlot_B.png", w=8, h=5) 

# ========================
# TFNetworkPlot
# ========================
cur_tf <- 'SMAD3'

# get a list of hub genes in the same module
tf_regulons <- GetTFRegulons(obj)
hub_df <- GetHubGenes(obj, n_hubs = 20)
cur_mod <- subset(hub_df, module == "brown") 
cur_mod_genes <- cur_mod$gene_name

# plot with custom gene labels
p1 <- TFNetworkPlot(
  obj, selected_tfs=cur_tf,
  label_TFs=0, label_genes=cur_mod_genes 
) + ggtitle('Custom gene labels')

save_gg(p1, 
        "TFNetworkPlot_B.png", w=8, h=5) 
