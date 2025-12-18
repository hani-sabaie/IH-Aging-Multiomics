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
  strategy = "C", # Strategy “C” retains all TF-gene pairs above a certain regulatory score 
  reg_thresh = 0.01,
)

tfrgs <- GetTFRegulons(obj) # filtered TF-gene pairs (regulons)

wcsv(tfnet, file.path("results_TF", "tf_network_full.csv"))
wcsv(tfrgs, file.path("results_TF", "tf_regulons.csv"))

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
saveRDS(obj,"../outputs/hdWGCNA_TFNet_obj.rds")

# ============================================================================ #

# ========================
# Load the object
# ========================
obj <- readRDS("../outputs/hdWGCNA_TFNet_obj.rds")

# ========================
# Differential regulon analysis
# ========================
# Get the cell barcodes for the groups of interest 
group1 <- obj@meta.data %>% subset(skeletal_muscle == "FAP3" & condition == "Aged") %>% rownames
group2 <- obj@meta.data %>% subset(skeletal_muscle == "FAP3" & condition == "Young") %>% rownames

# Calculate differential regulons
dregs <- FindDifferentialRegulons(
  obj,
  barcodes1 = group1,
  barcodes2 = group2,
  test.use='wilcox',
  wgcna_name = 'SCT',
  assay = 'SCT',
  pseudocount.use = 1
)

wcsv(dregs, file.path("results_TF", "differential_regulons.csv"))

# show the table
head(dregs)

# Use the dregs which we obtained above
p <- PlotDifferentialRegulons(obj, dregs)

# Show the plot
save_gg(p + theme(panel.background = element_rect(fill = 'white')), 
        "PlotDifferentialRegulons.png", w=8, h=5) 

# ========================
# Signac: Link peaks to genes
# ========================
# Provide the gene list
goi <- "SMAD3"

# Make sure we are on ATAC assay
DefaultAssay(obj) <- "ATAC"

# Add peak stats (GC etc.)
obj <- RegionStats(
  object = obj,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

# Link accessible regions to genes (cis peak-gene links)
obj <- LinkPeaks(
  object = obj,
  peak.assay = "ATAC",
  expression.assay = "SCT",
  genes.use = goi
)

# Extract the links as a data.frame
peak_links <- Links(obj)
peak_links_df <- as.data.frame(peak_links) %>%
  # Harmonize column names
  dplyr::rename(
    target_gene = gene, # Gene symbol linked to this peak
    peak_region = peak # Peak coordinates "chr-start-end"
  ) %>%
  dplyr::select(
    target_gene,
    peak_region,
    peak_score = score, # Correlation strength / link weight
    peak_z = zscore,   # z-scored strength
    peak_pval = pvalue  # p-value of link
  ) %>%
  # Keep only genes in the list of interest
  dplyr::filter(target_gene %in% goi) %>%
  dplyr::distinct()

# Clean TF-gene network (tfnet)
tfnet <- GetTFNetwork(obj)
tfnet_clean <- tfnet %>%
  dplyr::rename(
    TF = tf,    # TF name
    target_gene = gene   # Predicted target gene of that TF
  ) %>%
  dplyr::mutate(
    reg_score = Gain, # Feature importance from boosting model
    reg_score_signed = Gain * sign(Cor) # Sign the importance by correlation
  ) %>%
  dplyr::select(
    TF,
    target_gene,
    reg_score,
    reg_score_signed,
    Cor,
    Cover,
    Frequency
  )

# Clean differential regulon activity (Aged vs Young) from dregs
dregs_clean <- dregs %>%
  dplyr::mutate(
    direction = dplyr::case_when(
      avg_log2FC_deg > 0 ~ "Aged_up", # TF/regulon more active in Aged cells
      avg_log2FC_deg < 0 ~ "Young_up",  # TF/regulon more active in Young cells
      TRUE ~ "no_change"
    )
  ) %>%
  dplyr::rename(
    TF = tf
  ) %>%
  dplyr::select(
    TF,
    avg_log2FC_deg,
    p_val_adj_deg,
    avg_log2FC_positive,
    p_val_adj_positive,
    avg_log2FC_negative,
    p_val_adj_negative,
    module,
    kME,
    direction
  )

# Keep TFs with significant differential activity overall
sig_tfs <- dregs_clean %>%
  dplyr::filter(p_val_adj_deg < 0.05)

# Combine network evidence + TF differential status,
# and restrict to genes of interest
tf_targets_sig <- tfnet_clean %>%
  dplyr::inner_join(sig_tfs, by = "TF") %>%
  dplyr::filter(target_gene %in% goi) %>%
  dplyr::distinct()

# Make sure ATAC assay is active
DefaultAssay(obj) <- "ATAC"

# Add motifs 
JASPAR2024 <- JASPAR2024()
sq24 <- RSQLite::dbConnect(RSQLite::SQLite(), db(JASPAR2024))

pfm_core <- TFBSTools::getMatrixSet(
  x = sq24,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

obj <- AddMotifs(
  object = obj,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm_core
)

# Get motif object (motif PWMs and positions) from ATAC assay
motif_obj <- Motifs(obj[["ATAC"]])

# Collapse all motif hits into one GRanges with motif_id tag
motif_positions_list <- motif_obj@positions # named list of GRanges
motif_ids <- names(motif_positions_list) # motif IDs like "MA0004.1"

motif_hits_gr_list <- purrr::map(motif_ids, function(mtf) {
  gr <- motif_positions_list[[mtf]]
  if (length(gr) == 0) return(NULL)
  gr$motif_id <- mtf
  gr
})

motif_hits_gr <- do.call(c, motif_hits_gr_list)

# Make GRanges of peaks (peaks linked to our genes of interest)
peak_df_for_gr <- tidyr::separate(
  data = peak_links_df %>% dplyr::select(peak_region) %>% dplyr::distinct(),
  col = peak_region,
  into = c("seqnames", "start", "end"),
  sep = "-",
  remove = FALSE
)

peak_df_for_gr$start <- as.numeric(peak_df_for_gr$start)
peak_df_for_gr$end <- as.numeric(peak_df_for_gr$end)

peaks_gr <- GRanges(
  seqnames = peak_df_for_gr$seqnames,
  ranges = IRanges(start = peak_df_for_gr$start, end = peak_df_for_gr$end),
  peak_region = peak_df_for_gr$peak_region
)

# Overlap TF motif hits with those peaks
ov <- findOverlaps(motif_hits_gr, peaks_gr)

motif_hits_relevant <- motif_hits_gr[queryHits(ov)]
peak_matches <- peaks_gr[subjectHits(ov)]

tf_peak_direct_df <- data.frame(
  motif_id = motif_hits_relevant$motif_id,
  peak_region = peak_matches$peak_region,
  stringsAsFactors = FALSE
)

# motif_id -> TF name
motif_lookup <- data.frame(
  motif_id_raw = names(motif_obj@pwm),  # e.g. "MA0004.1"
  TF_raw = as.character(motif_obj@motif.names),  # e.g. "Arnt"
  stringsAsFactors = FALSE
)

motif_lookup <- motif_lookup %>%
  mutate(
    motif_id = motif_id_raw,
    TF = TF_raw
  ) %>%
  select(motif_id, TF) %>%
  distinct()

tf_peak_direct <- tf_peak_direct_df %>%
  left_join(motif_lookup, by = "motif_id")

tf_peak_direct <- tf_peak_direct %>%
  select(TF, peak_region) %>%
  distinct()

# peak -> gene link
tf_peak_gene_epigenetic <- tf_peak_direct %>%
  dplyr::inner_join(peak_links_df, by = "peak_region") %>%
  dplyr::distinct()

# Add TF -> target_gene regulatory evidence / differential aging evidence
tf_gene_peak_direct <- tf_targets_sig %>%
  dplyr::inner_join(
    tf_peak_gene_epigenetic,
    by = c("TF", "target_gene")
  ) %>%
  dplyr::distinct()


wcsv(
  tf_gene_peak_direct,
  file.path("results_TF", "TF_gene_peak_DIRECT_motifValidated_fap3.csv")
)

# High Confidence Full Integrated Filter
bulk_de <- read.csv(file = "bulk_de_sig_faps.csv", header = T, sep = ",")
smad3_de_df <- bulk_de %>%
  dplyr::filter(
    gene == goi) %>%
  dplyr::transmute(
    target_gene = gene,
    logFC_SMAD3 = avg_log2FC, # Aged vs Young
    p_val_SMAD3 = p_val,
    p_adj_SMAD3 = p_val_adj
  )

tf_gene_peak_final <- tf_gene_peak_direct %>%
  dplyr::left_join(smad3_de_df, by = "target_gene") %>%
  dplyr::mutate(
    s_reg = sign(avg_log2FC_deg),
    s_cor = sign(Cor),
    s_SMAD3 = sign(logFC_SMAD3),
    consistent_three_way = dplyr::case_when(
      s_cor > 0 & s_reg ==  s_SMAD3 & s_reg != 0 ~ TRUE,
      s_cor < 0 & s_reg == -s_SMAD3 & s_reg != 0 ~ TRUE,
      TRUE ~ FALSE
    )
  ) %>%
  dplyr::filter(
    # Sign Consistency 
    consistent_three_way,
    !is.na(logFC_SMAD3),
    
    # TFNet (expression network) 
    abs(Cor) >= 0.25, # strong TF-SMAD3 association
    reg_score >= 0.02, # moderate-strong boosting importance
    Frequency >= 0.02, # appears repeatedly across trees
    
    # Regulon activity 
    abs(avg_log2FC_deg) >= 0.25,
    p_val_adj_deg < 0.01,
    
    # Chromatin link
    peak_z >= 4,
    peak_pval < 1e-5,
    peak_score >= 0.05,
    
  )

wcsv(
  tf_gene_peak_final,
  file.path("results_TF", "TF_gene_peak_DIRECT_motifValidated_fap3_highconf.csv")
)

# ========================
# Save the object
# ========================
saveRDS(obj,"../outputs/hdWGCNA_TFNet_DEReg_L2G_obj.rds")

# ============================================================================ #

# ========================
# Load the object
# ========================
obj <- readRDS("../outputs/hdWGCNA_TFNet_obj.rds")

# ========================
# Module regulatory networks
# ========================
# ModuleRegulatoryHeatmap
p1 <- ModuleRegulatoryHeatmap(
  obj, feature='delta', dendrogram=FALSE
) + ggtitle('TFs only')
p2 <- ModuleRegulatoryHeatmap(
  obj, feature='delta',TFs_only=FALSE, 
  max_val=5, dendrogram=FALSE
) + ggtitle('All target genes')

save_gg(p1 | p2, 
        "ModuleRegulatoryHeatmap.png", w=8, h=5) 

# ModuleRegulatoryNetworkPlot
p1 <- ModuleRegulatoryNetworkPlot(
  obj, feature='positive',high_color='orange2', cutoff=0.5,max_val=50) + theme(panel.background = element_rect(fill = 'white'))
p2 <- ModuleRegulatoryNetworkPlot(
  obj, feature='negative',
  high_color='dodgerblue', cutoff=0.5,max_val=50) + 
  theme(panel.background = element_rect(fill = 'white'))

p3 <- ModuleRegulatoryNetworkPlot(
  obj, feature='positive',high_color='orange2',
  TFs_only=FALSE, cutoff=0.5,max_val=50) + theme(panel.background = element_rect(fill = 'white'))
p4 <- ModuleRegulatoryNetworkPlot(
  obj, feature='negative',
  high_color='dodgerblue',
  TFs_only=FALSE, cutoff=0.5,max_val=50) + 
  theme(panel.background = element_rect(fill = 'white'))

p_all <- (p1 + p2 + p3 + p4) + plot_layout(ncol = 2)
save_gg(p_all, 
        "ModuleRegulatoryNetworkPlot_TFs_and_all_target.png", w=12, h=9)

# Focus on specific module
# Source
modules <- GetModules(obj)
mods <- levels(modules$module)
p1 <- ModuleRegulatoryNetworkPlot(
  obj, feature='positive', 
  umap_background=TRUE,
  high_color='black',
  cutoff=0.1,
  loops=FALSE,
  focus_source = 'brown') + ggtitle('Focus on brown module as the source') +
  theme(panel.background = element_rect(fill = 'white'))

# Target
p2 <- ModuleRegulatoryNetworkPlot(
  obj, feature='positive', 
  umap_background=TRUE,
  high_color='black',
  loops=FALSE,
  cutoff=0.1,
  focus_target = 'brown') + ggtitle('Focus on brown module as the target') +
  theme(panel.background = element_rect(fill = 'white'))

save_gg(p1 | p2, 
        "ModuleRegulatoryNetworkPlot_focus_source_target.png", w=8, h=5)

# Visualize the regulon scores
# Select a TF of interest
cur_tf <- 'FOS'

# Add the regulon scores to the Seurat metadata
pos_regulon_scores <- GetRegulonScores(obj, target_type='positive')
neg_regulon_scores <- GetRegulonScores(obj, target_type='negative')
obj$pos_regulon_score <- pos_regulon_scores[,cur_tf]
obj$neg_regulon_score <- neg_regulon_scores[,cur_tf]

# Plot using FeaturePlot
p1 <- FeaturePlot(obj, feature=cur_tf) + umap_theme()
p2 <- FeaturePlot(obj, feature='pos_regulon_score', cols=c('lightgrey', 'red')) + umap_theme()
p3 <- FeaturePlot(obj, feature='neg_regulon_score', cols=c('lightgrey', 'seagreen')) + umap_theme()

p4 <- FeaturePlot(obj, feature=cur_tf) + umap_theme()
p5 <- FeaturePlot(obj, feature='pos_regulon_score', cols=c('lightgrey', 'red')) + umap_theme()
p6 <- FeaturePlot(obj, feature='neg_regulon_score', cols=c('lightgrey', 'seagreen')) + umap_theme()

p7 <- FeaturePlot(obj, feature=cur_tf) + umap_theme()
p8 <- FeaturePlot(obj, feature='pos_regulon_score', cols=c('lightgrey', 'red')) + umap_theme()
p9 <- FeaturePlot(obj, feature='neg_regulon_score', cols=c('lightgrey', 'seagreen')) + umap_theme()

p_all <- (p1+p2+p3+p4+p5+p6+p7+p8+p9) + plot_layout(ncol = 3)

save_gg(p_all, "FeaturePlot_Regulon_Scores.png", w=12, h=9)

