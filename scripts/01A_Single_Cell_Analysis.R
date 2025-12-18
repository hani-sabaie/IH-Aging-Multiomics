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
figdir <- "../outputs/sc/figs"
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
# Importing Data
# ========================
# List sample folders under ../data
files <- list.dirs(path = "../data/GSE268953/RNA_ATAC/", recursive = F, full.names = TRUE)

# Keep only folders that contain the full 10x trio
has_mtx  <- file.exists(file.path(files, "matrix.mtx.gz"))
has_feat <- file.exists(file.path(files, "features.tsv.gz"))
has_bc   <- file.exists(file.path(files, "barcodes.tsv.gz"))
files    <- files[ has_mtx & has_feat & has_bc ]

# Create a list of count matrices from 10x folders
input_list <- lapply(files, Read10X)

# Automated way to assign names
sample_names <- make.unique(basename(files))
names(input_list) <- sample_names

# Extract RNA and ATAC data
rna_counts_list  <- lapply(input_list, function(x) x$`Gene Expression`)
atac_counts_list <- lapply(input_list, function(x) x$Peaks)

# ========================
# Decontamination of ambient RNA (DecontX)
# ========================
sce_list <- setNames(vector("list", length(rna_counts_list)), names(rna_counts_list))
for (s in names(rna_counts_list)) {
  
  # Build a per-sample SCE
  sce_i <- SingleCellExperiment(list(counts = rna_counts_list[[s]]))
  
  # Run decontX
  sce_i <- decontX(sce_i)
  
  # UMAP
  umap <- reducedDim(sce_i, "decontX_UMAP")
  
  # Contamination plot
  p <- plotDecontXContamination(sce_i) + ggtitle(paste0("DecontX contamination — ", s))
  save_gg(p, paste0("Contam_", s, ".png"), w=8, h=5)
  
  # Store SCE for downstream steps
  sce_list[[s]] <- sce_i
}

# ========================
# Create Seurat object from SCE with decontX results
# ========================
dcx_list <- lapply(sce_list, function(sce_i) {
  m <- round(decontXcounts(sce_i))      
})

obj_list <- mapply(
    function(m, sname) {
    so <- CreateSeuratObject(counts = m, project = sname,
                             min.cells = 3, min.features = 200)
    so$sample <- sname
    so
  },
  m = dcx_list,
  sname = names(dcx_list),
  SIMPLIFY = FALSE
)

# Cells list
raw_cells_list <- lapply(obj_list, colnames)

# View obj_list 
obj_list

# Remove the original sparse matrices
rm(input_list)

# ========================
# Fragments indexing
# ========================
frag_dir   <- "../data/GSE268953/Fragments"
frag_files <- setNames(file.path(frag_dir, names(obj_list), "atac_fragments.tsv.gz"),
                       names(obj_list))

for (s in names(frag_files)) {
  fp    <- frag_files[[s]]                 
  plain <- sub("\\.gz$", "", fp)           
  
  # gunzip 
  R.utils::gunzip(fp, destname = plain, overwrite = TRUE, remove = FALSE)
  
  # sort by chr (col1) and start (col2)
  dt <- fread(plain, header = FALSE, sep = "\t", showProgress = FALSE)
  setorder(dt, V1, V2)
  fwrite(dt, plain, sep = "\t", col.names = FALSE)
  
  # bgzip back to the same .gz path and build tabix index
  Rsamtools::bgzip(plain, dest = fp, overwrite = TRUE)
  Rsamtools::indexTabix(fp, format = "bed")
  
  # remove the plain .tsv 
  file.remove(plain)
}

# ========================
# Build a unified peak set across samples
# ========================
# Collect GRanges of peaks from each sample's raw ATAC matrix
peak_gr_list <- lapply(atac_counts_list, function(mat){
  grange_counts <- StringToGRanges(rownames(mat), sep = c(":", "-"))
  
  # Keep standard chromosomes only
  grange_use <- seqnames(grange_counts) %in% standardChromosomes(grange_counts)
  grange_counts[grange_use]
})

# Union peaks
grl <- GenomicRanges::GRangesList(peak_gr_list)           
all_peaks <- unlist(grl, use.names = FALSE)               
peaks_union <- GenomicRanges::reduce(all_peaks)
peaks_union <- keepStandardChromosomes(peaks_union, pruning.mode = "coarse")

# Filter out bad peaks based on length 
pw <- width(peaks_union)
peaks_union <- peaks_union[pw > 20 & pw < 10000]

# Tidy seqlevels style / genome
seqlevelsStyle(peaks_union) <- "UCSC"
genome(peaks_union) <- "hg38"

# ========================
# Attach ATAC to the existing Seurat objects
# ========================
# Shared genome annotations for ATAC
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "hg38"  # UCSC name for GRCh38

for (s in names(obj_list)) {
  so <- obj_list[[s]]
  
  # Raw cells
  raw_cells <- intersect(colnames(so), raw_cells_list[[s]])
  
  # Build Fragment object; restrict to the cells of this sample
  frags_obj <- CreateFragmentObject(
    path = frag_files[[s]],
    cells = raw_cells,
    validate.fragments = TRUE)
    
  # Recount over the unified peak set
  mat <- FeatureMatrix(
    fragments = frags_obj,
    features  = peaks_union,
    cells     = raw_cells)
  common <- colnames(so)[colnames(so) %in% colnames(mat)]
  mat <- mat[, common, drop = FALSE]
  
  # Create ChromatinAssay (keep min.cells small to preserve union)
  chrom_assay <- CreateChromatinAssay(
    counts     = mat,
    genome     = "hg38",
    fragments  = frags_obj,
    min.cells  = 1,
    annotation = annotations
    )
  
    so[["ATAC"]] <- chrom_assay
    
    # Add prefix
    so <- RenameCells(so, add.cell.id = s)
    
    obj_list[[s]] <- so
  }

# ========================
# Merging Seurat objects
# ========================
keep <- c("obj_list")
rm(list = setdiff(ls(), keep))
gc()

merged <- obj_list[[1]]
for (i in 2:length(obj_list)) {
  merged <- merge(merged, y = obj_list[[i]], add.cell.ids = NULL, project = "GSE268953")
  gc()
}
obj <- merged

# ========================
# Inspect and work with the Seurat object
# ========================
# Add to metadata
obj$condition <- ifelse(str_detect(obj@meta.data$orig.ident, "^A"),
                        "Aged","Young")

# Cell names of the entire object
head(Cells(obj,layer="counts.Young_5"))
head(colnames(obj))

# Feature/gene names
head(Features(obj))
head(rownames(obj))

# List the assays
Assays(obj)

# Return number of cells across all layers
num_cells <- ncol(obj)

# Return number of features across all layers
num_features <- nrow(obj)

# To return row by column information
dim(obj)

# Visualize the number of cell counts per sample
save_gg(obj@meta.data %>% 
          ggplot(aes(x=orig.ident, fill=orig.ident)) + 
          geom_bar(color="black") +
          stat_count(geom = "text", colour = "black", size = 3.5, 
                     aes(label = ..count..),
                     position=position_stack(vjust=0.5), vjust = -0.25)+
          theme_classic() +
          theme(plot.title = element_text(hjust=0.5, face="bold")) +
          ggtitle("Number of Cells per Sample"), 
        "cell_counts_per_sample.png", w=8, h=5)

# ========================
# Save the object
# ========================
saveRDS(obj,"../outputs/decont_merged_obj.rds")

# ============================================================================ #

# ========================
# Load the object
# ========================
obj <- readRDS("../outputs/decont_merged_obj.rds")

# Examine the object
glimpse(obj)

# ========================
# Quality Control
# ========================
# ===== RNA QC metrics =====
DefaultAssay(obj) <- "RNA"

# Add percent mitochondrial
obj[["percent_mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

# Add percent ribosomal
obj[["percent_ribo"]] <- PercentageFeatureSet(obj, pattern="^RPS|^RPL")

# Add number of genes per UMI for each cell (complexity)
obj$log10GenesPerUMI <- log10(obj$nFeature_RNA) / log10(obj$nCount_RNA)

# ===== ATAC QC metrics =====
DefaultAssay(obj) <- "ATAC"

# Add 'nucleosome_signal'
obj <- NucleosomeSignal(object = obj)  

# Add 'TSS.enrichment'
obj <- TSSEnrichment(object = obj, fast = FALSE)   

# FRiP (fraction of reads in peaks)
sample_names <- c(paste0("Aged_", 1:5), paste0("Young_", 1:5))
fragments_files <- c(
  "C:/Users/Hani/Desktop/Hernia/data/GSE268953/Fragments/Aged_1/atac_fragments.tsv.gz",
  "C:/Users/Hani/Desktop/Hernia/data/GSE268953/Fragments/Aged_2/atac_fragments.tsv.gz",
  "C:/Users/Hani/Desktop/Hernia/data/GSE268953/Fragments/Aged_3/atac_fragments.tsv.gz",
  "C:/Users/Hani/Desktop/Hernia/data/GSE268953/Fragments/Aged_4/atac_fragments.tsv.gz",
  "C:/Users/Hani/Desktop/Hernia/data/GSE268953/Fragments/Aged_5/atac_fragments.tsv.gz",
  "C:/Users/Hani/Desktop/Hernia/data/GSE268953/Fragments/Young_1/atac_fragments.tsv.gz",
  "C:/Users/Hani/Desktop/Hernia/data/GSE268953/Fragments/Young_2/atac_fragments.tsv.gz",
  "C:/Users/Hani/Desktop/Hernia/data/GSE268953/Fragments/Young_3/atac_fragments.tsv.gz",
  "C:/Users/Hani/Desktop/Hernia/data/GSE268953/Fragments/Young_4/atac_fragments.tsv.gz",
  "C:/Users/Hani/Desktop/Hernia/data/GSE268953/Fragments/Young_5/atac_fragments.tsv.gz"
)
raw_from_merged <- function(x) sub("^.*_", "", x)
obj@meta.data$fragments_total <- NA_real_
cells_all  <- colnames(obj)
sample_vec <- as.character(obj@meta.data$sample)

for (i in seq_along(sample_names)) {
  cells_i <- cells_all[ sample_vec == sample_names[i] ]
  if (length(cells_i) == 0) next
  rb_i <- raw_from_merged(cells_i)
  
  cf <- CountFragments(fragments_files[i])   
  CB  <- as.character(cf[["CB"]])
  TOT <- as.numeric(cf[["frequency_count"]]) 
  
  m <- match(rb_i, CB)
  obj@meta.data[cells_i, "fragments_total"] <- TOT[m]
}

obj <- Signac::FRiP(
  object          = obj,
  assay           = "ATAC",
  total.fragments = "fragments_total",
  col.name        = "FRiP",
  verbose         = TRUE
)

obj$pct_reads_in_peaks <- 100 * obj$FRiP

# Blacklist ratio
data("blacklist_hg38", package = "Signac")
blacklist_hg38_unified <- blacklist_hg38
obj$blacklist_ratio <- FractionCountsInRegion(
  object = obj,
  assay  = "ATAC",
  regions = blacklist_hg38_unified)

# ===== Refresh metadata =====
metadata <- obj@meta.data

# ===== Pre-filter ATAC visualizations =====
# set colors
cnames <- setNames(rep(c("cyan3","darkgoldenrod1"),each=5),levels(factor(metadata$orig.ident))) 

# Plot (ATAC metrics)
save_gg(
  VlnPlot(
    obj,
    features = c("nCount_ATAC","TSS.enrichment","blacklist_ratio","nucleosome_signal","pct_reads_in_peaks"),
    group.by = "orig.ident", raster=FALSE, alpha=0.07, ncol = 5
  ) + scale_fill_manual(values = cnames),
  "VlnPlot_ATAC_metrics_pre_filt.png", w = 12, h = 5)

# Plot (ATAC density scatter)
samples <- unique(obj$orig.ident)
for (s in samples) {
  obj_s <- subset(obj, idents = s)
  p <- DensityScatter(obj_s, x = "nCount_ATAC", y = "TSS.enrichment",
                      log_x = TRUE, quantiles = TRUE) +
    ggplot2::ggtitle(s)
  save_gg(p, sprintf("DensityScatter_nCount_ATAC_vs_TSS_pre_filt_%s.png", s), w = 8, h = 5)
}

# Plot (ATAC fragment histogram)
obj$nucleosome_group <- ifelse(high_ns, "NS high (3×MAD)", "NS ok")
for (s in samples) {
  obj_s <- subset(obj, idents = s)
  p <- FragmentHistogram(object = obj_s, group.by = "nucleosome_group") +
    ggplot2::ggtitle(s)
  save_gg(p, sprintf("FragmentHistogram_nucleosome_group_pre_filt_%s.png", s), w = 8, h = 5)
}

# ===== Pre-filter RNA visualizations =====
DefaultAssay(obj) <- "RNA"

# Plot (RNA metrics)
save_gg(
  VlnPlot(
    obj,
    features = c("nCount_RNA","nFeature_RNA","percent_mt","percent_ribo","log10GenesPerUMI"),
    layer="counts",
    group.by = "orig.ident", raster=FALSE, alpha=0.07, ncol = 5
  ) + scale_fill_manual(values = cnames),
  "VlnPlot_RNA_metrics_pre_filt.png", w = 12, h = 5)

# Plot (RNA density plot)
p1 <- ggplot(metadata, aes(color=orig.ident, x=nCount_RNA, fill= orig.ident)) + 
          geom_density(alpha = 0.2) + 
          theme_classic() +
          scale_x_log10() 
p2 <- ggplot(metadata, aes(color=orig.ident, x=nFeature_RNA,fill=orig.ident)) +
          geom_density(alpha = 0.2) + 
          theme_classic() +
          scale_x_log10() 
p3 <- ggplot(metadata, aes(color=orig.ident, x=percent_mt,fill=orig.ident)) +
          geom_density(alpha = 0.2) + 
          scale_x_log10()+
          theme_classic()
p4 <- ggplot(metadata, aes(color=orig.ident, x=percent_ribo,fill=orig.ident)) +
          geom_density(alpha = 0.2) + 
          scale_x_log10()+
          theme_classic()
p5 <- ggplot(metadata, aes(color=orig.ident, x=log10GenesPerUMI,fill=orig.ident)) +
          geom_density(alpha = 0.2) + 
          scale_x_log10()+
          theme_classic()
p6 <- wrap_plots(list(p1,p2,p3,p4,p5), ncol = 2, guides = "collect") & 
  theme(legend.position = "bottom")
save_gg(p6, "DensityPlot_RNA_Metrics_pre_filt.png", w=8, h=5)

# Examine these metrics together
p7 <- FeatureScatter(obj, feature1 = "percent_mt", feature2 ="nFeature_RNA", 
                       group.by="orig.ident")
p8 <- FeatureScatter(obj, feature1 = "percent_ribo", feature2 ="nFeature_RNA", 
                       group.by="orig.ident")
p9 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                       group.by="orig.ident",log=TRUE)
p10 <- wrap_plots(list(p7,p8,p9), ncol = 2, guides = "collect") & 
  theme(legend.position = "right")
save_gg(p10, "FeatureScatter_RNA_Metrics_pre_filt.png", w=8, h=5)

p11 <- ggplot(metadata) +
          geom_point(aes(x=nCount_RNA,y=nFeature_RNA,fill=high_mito),shape=21,alpha=0.4) + 
          theme_classic() +
          scale_x_log10()+
          scale_y_log10()+
          facet_grid(.~sample)
        
p12 <- ggplot(metadata) +
          geom_point(aes(x=nCount_RNA,y=nFeature_RNA,fill=high_ribo),shape=21,alpha=0.4) + 
          theme_classic() +
          scale_x_log10()+
          scale_y_log10()+
          facet_grid(.~sample)
p13 <- wrap_plots(list(p11,p12), ncol = 1, guides = "collect") & 
  theme(legend.position = "right")
save_gg(p13, "FeatureScatter_RNA_Mito_Ribo_pre_filt.png", w=12, h=5)

# ===== Decide on thresholds =====
batch <- obj$sample

# RNA flags
extreme_lib <- isOutlier(obj$nCount_RNA, type = "both", log = TRUE, 
                         batch = batch, nmads = 3)
extreme_ngene <- isOutlier(obj$nFeature_RNA, type = "both", log = TRUE, 
                           batch = batch, nmads = 3)
high_mito <- isOutlier(obj$percent_mt, type = "higher", 
                       batch = batch, nmads = 3)
high_ribo <- isOutlier(obj$percent_ribo, type = "higher", 
                       batch = batch, nmads = 3)
low_complex <- isOutlier(obj$log10GenesPerUMI, type = "lower", 
                         batch = batch, nmads = 3)

# ATAC flags
extreme_npeaks <- isOutlier(obj$nCount_ATAC, type = "both", log = TRUE,
                            batch = batch, nmads = 3)
high_blacklist <- isOutlier(obj$blacklist_ratio, type = "higher", log = FALSE,
                            batch = batch, nmads = 3)
low_tss <- isOutlier(obj$TSS.enrichment, type = "lower", log = FALSE,
                            batch = batch, nmads = 3)
high_ns <- isOutlier(obj$nucleosome_signal, type = "higher", log = FALSE,
                            batch = batch, nmads = 3)
low_frip <- isOutlier(obj$pct_reads_in_peaks, type = "lower", log = FALSE,
                            batch = batch, nmads = 3)

# ===== Apply filtering =====
metadata <- obj@meta.data
metadata$extreme_lib <- extreme_lib
metadata$extreme_ngene <- extreme_ngene
metadata$high_mito <- high_mito
metadata$high_ribo <- high_ribo
metadata$low_complex <- low_complex

metadata$extreme_npeaks <- extreme_npeaks
metadata$high_blacklist <- high_blacklist
metadata$low_tss <- low_tss
metadata$high_ns <- high_ns
metadata$low_frip <- low_frip

keep <- metadata |> 
  filter(!extreme_lib, !extreme_ngene, !high_mito, !high_ribo, !low_complex,
         !extreme_npeaks, !high_blacklist, !low_tss, !high_ns, !low_frip) |> 
  rownames_to_column("Cell") |> pull(Cell)
obj_filt <- subset(obj, cells = keep)
metadata_filt <- obj_filt@meta.data

# ===== Post-filter ATAC visualization =====
DefaultAssay(obj_filt) <- "ATAC"

cnames <- setNames(rep(c("cyan3","darkgoldenrod1"), each = 5), levels(factor(metadata_filt$orig.ident)))

# Plot (ATAC metrics)
save_gg(
  VlnPlot(
    obj_filt,
    features = c("nCount_ATAC","TSS.enrichment","blacklist_ratio","nucleosome_signal","pct_reads_in_peaks"),
    group.by = "orig.ident", raster=FALSE, alpha=0.07, ncol = 5
  ) + scale_fill_manual(values = cnames),
  "VlnPlot_ATAC_metrics_post_filt.png", w = 12, h = 5)

# Plot (ATAC density scatter)
samples <- unique(obj_filt$orig.ident)
for (s in samples) {
  obj_s <- subset(obj_filt, idents = s)
  p <- DensityScatter(obj_s, x = "nCount_ATAC", y = "TSS.enrichment",
                      log_x = TRUE, quantiles = TRUE) +
    ggplot2::ggtitle(s)
  save_gg(p, sprintf("DensityScatter_nCount_ATAC_vs_TSS_post_filt_%s.png", s), w = 8, h = 5)
}

# ===== Post-filter RNA visualization =====
DefaultAssay(obj_filt) <- "RNA"

# set colors
cnames<-setNames(rep(c("cyan3","darkgoldenrod1"),each=5),levels(factor(metadata_filt$orig.ident))) 

# plot (RNA metrics)
save_gg(
  VlnPlot(
    obj_filt,
    features = c("nCount_RNA","nFeature_RNA","percent_mt","percent_ribo","log10GenesPerUMI"),
    layer="counts",
    group.by = "orig.ident", raster=FALSE, alpha=0.07, ncol = 5
  ) + scale_fill_manual(values = cnames),
  "VlnPlot_RNA_metrics_post_filt.png", w = 12, h = 5)

# Plot (RNA density plot)
p14 <- ggplot(metadata_filt, aes(color=orig.ident, x=nCount_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() 
p15 <- ggplot(metadata_filt, aes(color=orig.ident, x=nFeature_RNA,fill=orig.ident)) +
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() 
p16 <- ggplot(metadata_filt, aes(color=orig.ident, x=percent_mt,fill=orig.ident)) +
  geom_density(alpha = 0.2) + 
  scale_x_log10()+
  theme_classic()
p17 <- ggplot(metadata_filt, aes(color=orig.ident, x=percent_ribo,fill=orig.ident)) +
  geom_density(alpha = 0.2) + 
  scale_x_log10()+
  theme_classic()
p18 <- ggplot(metadata_filt, aes(color=orig.ident, x=log10GenesPerUMI,fill=orig.ident)) +
  geom_density(alpha = 0.2) + 
  scale_x_log10()+
  theme_classic()
p19 <- wrap_plots(list(p14,p15,p16,p17,p18), ncol = 2, guides = "collect") & 
  theme(legend.position = "bottom")
save_gg(p19, "DensityPlot_RNA_Metrics_post_filt.png", w=8, h=5)

# Examine these metrics together
p20 <- FeatureScatter(obj_filt, feature1 = "percent_mt", feature2 ="nFeature_RNA", 
                     group.by="orig.ident")
p21 <- FeatureScatter(obj_filt, feature1 = "percent_ribo", feature2 ="nFeature_RNA", 
                     group.by="orig.ident")
p22 <- FeatureScatter(obj_filt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                     group.by="orig.ident",log=TRUE)
p23 <- wrap_plots(list(p20,p21,p22), ncol = 2, guides = "collect") & 
  theme(legend.position = "right")
save_gg(p23, "FeatureScatter_RNA_Metrics_post_filt.png", w=8, h=5)

high_mito <- { x <- high_mito[match(colnames(obj_filt), names(high_mito))]; x[is.na(x)] <- FALSE; x }
brks <- c(300, 1000, 3000, 10000, 30000)
p24 <- ggplot(metadata_filt) +
  geom_point(aes(x=nCount_RNA,y=nFeature_RNA,fill=high_mito),shape=21,alpha=0.4) + 
  theme_classic() +
  scale_x_log10(breaks = brks,
                labels = function(v) {
                  labs <- formatC(v, format = "f", digits = 0)
                  labs[seq_along(labs) %% 2 == 0] <- ""   
                  labs
                })+
  scale_y_log10()+
  facet_grid(.~sample)
high_ribo <- { x <- high_ribo[match(colnames(obj_filt), names(high_ribo))]; x[is.na(x)] <- FALSE; x }
p25 <- ggplot(metadata_filt) +
  geom_point(aes(x=nCount_RNA,y=nFeature_RNA,fill=high_ribo),shape=21,alpha=0.4) + 
  theme_classic() +
  scale_x_log10(breaks = brks,
                labels = function(v) {
                  labs <- formatC(v, format = "f", digits = 0)
                  labs[seq_along(labs) %% 2 == 0] <- ""   
                  labs
                })+
  scale_y_log10()+
  facet_grid(.~sample)
p26 <- wrap_plots(list(p24,p25), ncol = 1, guides = "collect") & 
  theme(legend.position = "right")
save_gg(p26, "FeatureScatter_RNA_Mito_Ribo_post_filt.png", w=12, h=5)

# ========================
# Save the object
# ========================
saveRDS(obj_filt,"../outputs/decont_merged_filt_obj.rds")

# ============================================================================ #

# ========================
# Load the object
# ========================
obj <- readRDS("../outputs/decont_merged_filt_obj.rds")

# ========================
# Per-sample scDblFinder 
# ========================
DefaultAssay(obj) <- "RNA"
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])

# ===== Run scDblFinder =====
sce <- as.SingleCellExperiment(obj)
bp <- SnowParam(3, RNGseed=1234)
sce <- scDblFinder(sce, samples = "sample", BPPARAM = bp) # dbr = multiplet_rate
table(sce$scDblFinder.class)
sce@colData@listData %>% as.data.frame() %>% head()

# ===== Explore results =====
meta_scdblfinder <- sce@colData@listData %>% as.data.frame() %>% 
  dplyr::select(starts_with('scDblFinder')) # 'scDblFinder.class'
head(meta_scdblfinder)

# ===== Add the results to Seurat object =====
rownames(meta_scdblfinder) <- sce@colData@rownames
head(meta_scdblfinder)
obj <- AddMetaData(object = obj, metadata = meta_scdblfinder %>% dplyr::select('scDblFinder.class'))
head(obj@meta.data)
table(obj$scDblFinder.class)
rm(list = c('meta_scdblfinder', 'sce'))

# ===== QC stats and summary table =====
# Doublet stats
# Check how doublets singlets differ in QC measures per sample
save_gg(VlnPlot(obj, group.by = 'sample', split.by = "scDblFinder.class",
                features = c("nCount_RNA","nFeature_RNA","percent_mt","percent_ribo","log10GenesPerUMI"), 
                ncol = 3, pt.size = 0) + theme(legend.position = "bottom"),
        "VlnPlot_doublets_singlets_QC.png", w=12, h=5)

# Total count of cells in each group and the percentage of doublets and singlets for each sample
doublets_summary <- obj@meta.data %>% 
  group_by(sample, scDblFinder.class) %>% 
  summarise(total_count = n(),.groups = 'drop') %>% as.data.frame() %>% ungroup() %>%
  group_by(sample) %>%
  mutate(countT = sum(total_count)) %>%
  group_by(scDblFinder.class, .add = TRUE) %>%
  mutate(percent = paste0(round(100 * total_count/countT, 2),'%')) %>%
  dplyr::select(-countT)
doublets_summary
write.table(doublets_summary, file = file.path("../outputs", paste0('scDblFinder_doublets_summary.txt')), quote = FALSE, row.names = FALSE, sep = '\t')

# ===== Remove doublets =====
obj_nodoub <- subset(obj, scDblFinder.class == "singlet")

# ========================
# Save the object
# ========================
obj_nodoub[["RNA"]] <- split(obj_nodoub[["RNA"]], f = obj_nodoub$sample)
saveRDS(obj_nodoub,"../outputs/decont_merged_filt_nodoub_obj.rds")

# ============================================================================ #

# ========================
# Load the object
# ========================
obj <- readRDS("../outputs/decont_merged_filt_nodoub_obj.rds")

# ========================
# Normalization, find variable features, and scale
# ========================
# ===== Normalize =====
DefaultAssay(obj) <- "RNA"
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])

obj <- NormalizeData(obj,
                     normalization.method = "LogNormalize",
                     scale.factor = 10000)

# ===== Find variable features =====
obj <- FindVariableFeatures(obj, selection.method = "vst",
                            nfeatures = 2000)

# ===== Scaling the data =====
# all_genes <- rownames(x = obj)
# obj <- ScaleData(object = obj, features = all_genes)
obj <- ScaleData(obj)

# ========================
# Scaling and regression of sources of unwanted variation
# ========================
# ===== Cell Cycle Scoring =====
s_genes <- cc.genes$s.genes
g2m_genes <- cc.genes$g2m.genes
obj <- CellCycleScoring(obj,g2m.features = g2m_genes, s.features = s_genes)

# Visualize the distribution of cell cycle markers across
allm <- FindAllMarkers(obj, group.by = "Phase", only.pos = TRUE,
                       min.pct = 0.0, logfc.threshold = 0.0)
allm <- transform(allm, dpct = pct.1 - pct.2)

top_by_phase <- allm |>
  dplyr::filter(p_val_adj < 0.05,
         avg_log2FC >= 0.25,
         pct.1 >= 0.2,
         dpct >= 0.1) |>
  group_by(cluster) |>
  arrange(desc(avg_log2FC), desc(dpct), p_val_adj, .by_group = TRUE) |>
  slice_head(n = 4) |>
  ungroup()

top_by_phase 

save_gg(RidgePlot(obj, features = c("POLA1", "NASP", "SMC4", "COL1A2"), group.by = "Phase", ncol = 2), 
                  "RidgePlot_cell_cycle.png", w=8, h=5)

# ===== Perform PCA =====
obj <- RunPCA(obj)

# ===== Visualize the PCA, grouping by cell cycle phase =====
save_gg(DimPlot(obj, reduction = "pca", group.by = "Phase"), 
        "PCA_cell_cycle_before_regress.png", w=8, h=5)

# ===== Apply regression variables =====
# Define variables in metadata to regress
# vars_to_regress <- c("S.Score", "G2M.Score")

# Regress out the uninteresting sources of variation in the data
# obj <- ScaleData(object = obj,
                 # vars.to.regress = vars_to_regress, 
                 # verbose = FALSE)

# Re-run the PCA
# obj <- RunPCA(obj)
# save_gg(DimPlot(obj, reduction = "pca", group.by = "Phase"), 
        # "PCA_cell_cycle_after_regress.png", w=8, h=5)

# ========================
# Save the metadata
# ========================
cc_meta <- obj@meta.data %>%
  dplyr::select(S.Score, G2M.Score, Phase) %>%
  tibble::rownames_to_column("cell")
saveRDS(cc_meta, "../outputs/cellcycle_meta.rds")

# ============================================================================ #

# ========================
# Load the object
# ========================
obj <- readRDS("../outputs/decont_merged_filt_nodoub_obj.rds")

# ========================
# Explanatory Variables
# ========================
# ===== plotExplanatoryVariables =====
DefaultAssay(obj) <- "RNA"

# Convert to SCE and create raw logcounts from RNA counts
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
sce <- as.SingleCellExperiment(obj)
assay(sce, "logcounts_raw") <- log1p(counts(sce)) 

# Pull QC metrics directly from metadata
qc_vars <- c(
  "nFeature_RNA","nCount_RNA","percent_mt","percent_ribo",
  "log10GenesPerUMI","S.Score","G2M.Score","sample")
qc_vars <- qc_vars[ qc_vars %in% colnames(obj@meta.data) ]

# Copy columns to colData
for (v in qc_vars) colData(sce)[[v]] <- obj@meta.data[[v]]

# Choose which variables to visualize/assess
vars <- c("nFeature_RNA","nCount_RNA","percent_mt","percent_ribo",
          "log10GenesPerUMI","S.Score","G2M.Score","sample")
vars <- vars[ vars %in% colnames(colData(sce)) ]

# Density of gene-wise marginal R^2 for explanatory variables
# Curves further to the right => variable explains more variance across genes
save_gg(plotExplanatoryVariables(sce, exprs_values = "logcounts_raw", 
                                 variables    = vars), 
        "plotExplanatoryVariables.png", w=8, h=5)

# ===== PC/QC R^2 =====
sce <- runPCA(sce, exprs_values = "logcounts_raw", ncomponents = 20)
pcs_to_check <- 1:5
pc_mat <- reducedDim(sce, "PCA")

lin_r2 <- function(y, x) {
  # Fit a simple linear model and return R^2
  df <- data.frame(y = y, x = x)
  m  <- if (is.numeric(x)) lm(y ~ x, df) else lm(y ~ 0 + x, df)
  summary(m)$r.squared
}

pcqc_r2 <- sapply(vars, function(v) {
  x <- colData(sce)[[v]]
  sapply(pcs_to_check, function(k) lin_r2(pc_mat[,k], x))
})
rownames(pcqc_r2) <- paste0("PC", pcs_to_check)

# Print numeric table
print(round(pcqc_r2, 3))

# Heatmap
dt <- as.data.table(pcqc_r2, keep.rownames = "PC")  
df_pcqc <- reshape2::melt(dt, id.vars = "PC",
                variable.name = "Variable",
                value.name = "R2")
df_pcqc <- reshape2::melt(pcqc_r2, varnames = c("PC","Variable"), value.name = "R2")
save_gg(ggplot(df_pcqc, aes(Variable, PC, fill = R2)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)),
  "QC_PC_Heatmap.png", w=8, h=5)

# ============================================================================ #

# ========================
# Load the object and metadata
# ========================
obj <- readRDS("../outputs/decont_merged_filt_nodoub_obj.rds")
cc_meta <- readRDS("../outputs/cellcycle_meta.rds")
cc_meta <- cc_meta[cc_meta$cell %in% colnames(obj), ]
rownames(cc_meta) <- cc_meta$cell
cc_meta$cell <- NULL
obj <- AddMetaData(obj, metadata = cc_meta)

# ========================
# SCTransform and Regression of Sources of Unwanted Variation
# ========================
DefaultAssay(obj) <- "RNA"

vars_to_regress <- unique("percent_mt")
obj <- SCTransform(obj, 
                   vars.to.regress = vars_to_regress, 
                   verbose = FALSE)

# Check default assay
DefaultAssay(object = obj) 

# ========================
# Save the object
# ========================
saveRDS(obj,"../outputs/decont_merged_filt_nodoub_cc_sct_obj.rds")

# ============================================================================ #

# ========================
# Load the object
# ========================
obj <- readRDS("../outputs/decont_merged_filt_nodoub_cc_sct_obj.rds")

# ========================
# Linear dimension reduction (RNA analysis)
# ========================
DefaultAssay(obj) <- "SCT"

# ===== Run PCA =====
obj <- RunPCA(obj, verbose = FALSE, assay = "SCT")

# ===== Exploring the PCA results =====
# DimHeatmap 1 to 30
png(file.path(figdir, "DimHeatmap_PCA_1_to_9.png"), width = 12, height = 12, units = "in", res = 300, bg = "white")
DimHeatmap(obj, dims = 1:9, cells = 500, balanced = TRUE, ncol=3, combine = F)
dev.off()

png(file.path(figdir, "DimHeatmap_PCA_10_to_18.png"), width = 12, height = 12, units = "in", res = 300, bg = "white")
DimHeatmap(obj, dims = 10:18, cells = 500, balanced = TRUE, ncol=3, combine = F)
dev.off()

png(file.path(figdir, "DimHeatmap_PCA_19_to_30.png"), width = 12, height = 12, units = "in", res = 300, bg = "white")
DimHeatmap(obj, dims = 19:30, cells = 500, balanced = TRUE, ncol=3, combine = F)
dev.off()

# Elbow plot
save_gg(ElbowPlot(obj, ndims = 40), 
        "ElbowPlot_PCA_1_to_40.png", w=8, h=5)

# ========================
# Linear dimension reduction (ATAC analysis)
# ========================
DefaultAssay(obj) <- "ATAC"

# ===== Run LSI =====
obj <- FindTopFeatures(obj, min.cutoff = 'q0')
obj <- RunTFIDF(obj)
obj <- RunSVD(obj, reduction.key = "LSI_", reduction.name = "lsi")

# ===== Exploring the LSI results =====
# Elbow plot
save_gg(
  ElbowPlot(obj, reduction = "lsi", ndims = 40),
  "ElbowPlot_LSI_1_to_40.png", w = 8, h = 5)

# Depth correlation (ATAC-specific sanity check)
# LSI_1 is often highly correlated with library size; if so, skip it.
dc <- DepthCor(obj, reduction = "lsi")  # correlation of LSI components with depth (nCount_ATAC)
save_gg(dc, "DepthCor_LSI.png", w = 8, h = 5)

# ========================
# Save the object
# ========================
saveRDS(obj,"../outputs/decont_merged_filt_nodoub_cc_sct_reduc_obj.rds")

# ============================================================================ #

# ========================
# Load the object
# ========================
obj <- readRDS("../outputs/decont_merged_filt_nodoub_cc_sct_reduc_obj.rds")

# ========================
# SCT integration
# ========================
# ===== Perform integration =====
DefaultAssay(obj) <- "SCT"
obj <- IntegrateLayers(
  object = obj,
  method = HarmonyIntegration, 
  orig.reduction = "pca",
  new.reduction = "harmony",
  normalization.method = "SCT",
  verbose = F)

# ========================
# ATAC integration
# ========================
# ===== Perform integration =====
DefaultAssay(obj) <- "ATAC"
obj <- RunHarmony(
  object = obj,
  group.by.vars = "sample",
  reduction = "lsi",
  reduction.save = "harmony.atac",
  assay.use = "ATAC",
  project.dim = F)

# ========================
# Silhouette
# ========================
pcs <- 40
pcs_lsi <- 40

# Candidate k values (simple & practical)
k.candidates <- sort(unique(c(20, 30, 40, 50, 60, 70, 80, 90, 100)))

# Silhouette on a subsample to avoid O(N^2)
sil_max_n <- 5000
set.seed(1234)
idx_sil <- if (N > sil_max_n) sample(seq_len(N), sil_max_n) else seq_len(N)

# "Small cluster" threshold = 0.2% of total 
thresh_small <- 0.002 * N

score.df <- data.frame(
  k = numeric(),
  n_clusters = numeric(),
  mean_silhouette = numeric(),
  p_small = numeric(),
  min_size = numeric(),
  stringsAsFactors = FALSE
)

for (kval in k.candidates) {
  obj.tmp <- obj
  
  # Build WNN
  obj.tmp <- FindMultiModalNeighbors(
    obj.tmp,
    reduction.list = list("pca", "lsi"),
    dims.list      = list(1:pcs, 2:pcs_lsi),  # LSI from component 2
    k.nn = kval,
    modality.weight.name = paste0("RNA.weight.", kval)
  )
  
  # Cluster on WNN graph
  obj.tmp <- FindClusters(
    obj.tmp,
    graph.name = "wsnn",
    resolution = 0.8
  )
  
  # UMAP for visualization & silhouette
  obj.tmp <- RunUMAP(
    obj.tmp,
    nn.name = "weighted.nn",
    reduction.name = paste0("wnn.umap.k", kval),
    reduction.key  = paste0("wnnUMAPk", kval, "_"),
    verbose = FALSE
  )
  
  # Cluster size stats
  cl <- Idents(obj.tmp)
  sz <- as.numeric(table(cl))
  p_small <- mean(sz < thresh_small)
  min_size <- min(sz)
  
  # Silhouette on subsample
  emb <- Embeddings(obj.tmp, paste0("wnn.umap.k", kval))
  emb_sub <- emb[idx_sil, , drop = FALSE]
  cl_sub  <- cl[idx_sil]
  
  mean_sil <- if (length(unique(cl_sub)) < 2) {
    NA_real_
  } else {
    sil.obj <- silhouette(as.integer(factor(cl_sub)), dist(emb_sub))
    mean(sil.obj[, "sil_width"])
  }
  
  score.df <- rbind(
    score.df,
    data.frame(
      k = kval,
      n_clusters = length(unique(cl)),
      mean_silhouette = mean_sil,
      p_small = p_small,
      min_size = min_size
    )
  )
}

# See best ks by silhouette first, then by fewer tiny clusters
score.df[order(-score.df$mean_silhouette, score.df$p_small, -score.df$min_size), ]

# ========================
# WNN graph
# ========================
# ===== WNN analysis =====
obj <- FindMultiModalNeighbors(obj, reduction.list = list("harmony", "harmony.atac"), dims.list = list(1:pcs, 2:pcs_lsi), k.nn = 100)
obj <- FindClusters(obj, graph.name = "wsnn", algorithm = 4, random.seed = 1234,
                    resolution = c(0.4,0.6,0.8,1,1.2,1.4,1.6))

# ===== Dimplots =====
# WNN
obj <- RunUMAP(obj, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
groups <- c("orig.ident",
            "wsnn_res.0.4","wsnn_res.0.6",
            "wsnn_res.0.8","wsnn_res.1","wsnn_res.1.2","wsnn_res.1.4","wsnn_res.1.6")
p30 <- lapply(groups, function(g){
  pl <- DimPlot(
    obj,
    reduction = "wnn.umap",
    group.by  = g,
    alpha     = 0.4,
    label     = (g != "orig.ident"),
    label.size = 2.5
  )
  if (g != "orig.ident") {
    pl <- pl + NoLegend() +
      guides(color = "none", fill = "none", shape = "none",
             alpha = "none", size = "none") +
      theme(legend.position = "none")
  }
  pl
})

combo <- wrap_plots(p30, ncol = 2) & 
  theme(legend.position = "right")

save_gg(combo, "WNN_clusters.png", w = 12, h = 15)

# Clusteree
save_gg(clustree(obj@meta.data, prefix = "wsnn_res.", from = 0.4, to = 1.6), 
        "Clusteree_WNN.png", w=12, h=8)

# SCT
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:pcs, assay = "SCT", k.param = 100)
obj <- FindClusters(obj, graph.name = "SCT_snn", algorithm = 4, random.seed = 1234,
                    resolution = 0.8)
obj <- RunUMAP(obj, dims = 1:pcs, reduction.name = 'harmony.sct.umap', reduction.key = 'sctUMAP_')
groups <- c("orig.ident","SCT_snn_res.0.8")
p30 <- lapply(groups, function(g){
  pl <- DimPlot(
    obj,
    reduction = "harmony.sct.umap",
    group.by  = g,
    alpha     = 0.4,
    label     = (g != "orig.ident"),
    label.size = 2.5
  )
  if (g != "orig.ident") {
    pl <- pl + NoLegend() +
      guides(color = "none", fill = "none", shape = "none",
             alpha = "none", size = "none") +
      theme(legend.position = "none")
  }
  pl
})

combo <- wrap_plots(p30, ncol = 2) & 
  theme(legend.position = "right")

save_gg(combo, "SCT_clusters.png", w = 8, h = 5)

# ATAC
obj <- FindNeighbors(obj, reduction = "harmony.atac", dims = 2:pcs_lsi, assay = "ATAC", k.param = 100)
obj <- FindClusters(obj, graph.name = "ATAC_snn", algorithm = 4, random.seed = 1234,
                    resolution = 0.8)
obj <- RunUMAP(obj, dims = 2:pcs_lsi, reduction.name = 'harmony.atac.umap', reduction.key = 'atacUMAP_')
groups <- c("orig.ident","ATAC_snn_res.0.8")
p30 <- lapply(groups, function(g){
  pl <- DimPlot(
    obj,
    reduction = "harmony.atac.umap",
    group.by  = g,
    alpha     = 0.4,
    label     = (g != "orig.ident"),
    label.size = 2.5
  )
  if (g != "orig.ident") {
    pl <- pl + NoLegend() +
      guides(color = "none", fill = "none", shape = "none",
             alpha = "none", size = "none") +
      theme(legend.position = "none")
  }
  pl
})

combo <- wrap_plots(p30, ncol = 2) & 
  theme(legend.position = "right")

save_gg(combo, "ATAC_clusters.png", w = 8, h = 5)

# ===== Compare SCT, ATAC, WNN =====
p31 <- DimPlot(obj, reduction = "harmony.sct.umap", group.by = "SCT_snn_res.0.8", label = TRUE, label.size = 2.5) + ggtitle("SCT")
p32 <- DimPlot(obj, reduction = "harmony.atac.umap", group.by = "ATAC_snn_res.0.8", label = TRUE, label.size = 2.5) + ggtitle("ATAC")
p33 <- DimPlot(obj, reduction = "wnn.umap", group.by = "wsnn_res.0.8", label = TRUE, label.size = 2.5) + ggtitle("WNN")
p34 <- p31 + p32 + p33 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
save_gg(p34, "SCT_ATAC_WNN_Clusters.png", w=15, h=5)

# ===== Join layers =====
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])

# ========================
# Save the object
# ========================
saveRDS(obj,"../outputs/decont_merged_filt_nodoub_cc_sct_reduc_clust_integ_obj.rds")

# ============================================================================ #

# ========================
# Load the object
# ========================
obj <- readRDS("../outputs/decont_merged_filt_nodoub_cc_sct_reduc_clust_integ_obj.rds")

# Reviewing the data
table(obj$sample)

# ========================
# Annotation
# ========================
# ===== Preparation =====
# PrepSCTFindMarkers
obj <- PrepSCTFindMarkers(obj, verbose=T)

# Define identity
# Idents(obj) <- "wsnn_res.0.8"
# lev_num <- sort(as.integer(levels(Idents(obj))))
# Idents(obj) <- factor(Idents(obj), levels = as.character(lev_num))

Idents(obj) <- "skeletal_muscle"
Idents(obj) <- factor(Idents(obj),
                      levels = c("FAP1","FAP2","FAP3","FAP4","SMC1","SMC2",
                                 "Pericyte","EC1","EC2","EC3","Macrophage",
                                 "B/T/NK","MuSC","MSM"))
table(Idents(obj))

# Default assay
DefaultAssay(obj) <- "SCT"

# ===== Finding Cluster Biomarkers =====
# Find all markers
obj_markers <- FindAllMarkers(obj, test.use="wilcox", min.pct=0.25, only.pos=TRUE)

# Ordering the results
obj_markers <- obj_markers %>% 
  arrange(cluster,desc(avg_log2FC), desc(p_val_adj))

# Heatmap of top potential marker genes
top <- obj_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1 & p_val_adj < 0.05) %>%
  slice_head(n = 10) %>%
  ungroup() 
save_gg(DoHeatmap(obj, features = top$gene) + NoLegend(),
        "Heatmap_top_10_markers.png", w=15, h=15)

# Top 5 potential marker genes per cluster
# top5PerCluster <- matrix(ncol=7)
# colnames(top5PerCluster) = colnames(obj_markers)
# for (i in obj_markers$cluster){
#   top5PerCluster = rbind(top5PerCluster, head(obj_markers[which(obj_markers$cluster==i),], 5))
# }
# top5PerCluster <- top5PerCluster[-1,]
# top5PerCluster

# Heatmap of top 5 potential marker per cluster
# save_gg(DoHeatmap(obj,features = top5PerCluster$gene, slot="scale.data"),
#        "Heatmap_top5_per_cluster.png", w=15, h=12) #  In order to ensure that all genes are represented, the data slot can be directly called as the data reference, but it will not be z-scaled.

# Compare marker expression across clusters
mk <- c("MYLPF","TNNC2","ACTA1","CKM","TTN",
        "CHRDL2","PAX7","MYF5","APOC1","DLK1",
        "CD3D","IL7R","GNLY","NKG7","CCL5","CXCR4","SRGN","PTPRC","CCL4",
        "LYZ","TYROBP","CCL4","C1QA","LYVE1","SRGN","CD74",
        "CAV1","PECAM1","CLDN5",
        "CLDN5","PECAM1","CDH5","FABP4","CXCL2","VWF","ICAM1",
        "VCAM1","PLVAP","IL6","ACKR1","SELE","TM4SF1","CD74","PECAM1",
        "NDUFA4L2","CD36","RGS5",
        "TAGLN","ACTA2","MYL9",
        "MYH11","ACTA2","MYL9","TAGLN","NDUFA4L2","LYVE1","CD74","LYZ","TYROBP","IL7R",
        "COL1A1","PDGFRA","DCN","GSN","FBN1","COL1A2","FSTL1","NCAM1",
        "MME","PTGDS","MYOC","IGF1","CXCL14","PLAC9","ABCA8","APOD","SMOC2","CFD","COL15A1",
        "CD248","PCOLCE2","CD55","PRG4","TNXB","MFAP5","IGFBP6",
        "SFRP2","GPC3",
        "PLA2G2A","SEMA3C"
        )
mkv <- c("PDGFRA","NCAM1",
        "MME","PTGDS",
        "CD55","PCOLCE2",
        "GPC3","SFRP2",
        "PRG4","PLA2G2A"
        )
mkf <- c("PDGFRA","DCN","GSN","NCAM1",
         "MME","PTGDS",
         "CD55","PCOLCE2",
         "GPC3","SFRP2",
         "PRG4","PLA2G2A")

save_gg(VlnPlot(obj, features = mkv, pt.size = 0,stack = T,flip = T,
                idents = c("FAP1","FAP2","FAP3","FAP4")) + NoLegend() ,
        "VlnPlot_markers_clust.png", w=9, h=12)
save_gg(FeaturePlot(obj,features = mkf, reduction = "wnn.umap",ncol = 4),
        "FeaturePlot_markers_clust.png", w=15, h=12)
save_gg(DotPlot(obj, features= unique(mk),dot.scale = 9,dot.min = 0.25, col.min = 0,col.max = 2),
        "DotPlot_markers_clust.png", w=49, h=12)

# =====  Cell annotation by differential expressed gene markers =====
# Idents(obj) <- "wsnn_res.0.8"
# avgExp <- AverageExpression(obj, assay="SCT", 
#                             features = c()$SCT
# avgExp
# save_gg(DimPlot(obj,label=T),
#         "Umap_annotation.png", w=8, h=5)
# save_gg(FeaturePlot(obj, features=c(FAP1, FAP2, FAP3, FAP4), ncol=3, order=T, reduction = "wnn.umap"),
#         "FeaturePlot_annotation_markers.png", w=8, h=5)

# =====  Assign cell identities =====
skeletal_muscle <- vector(length=ncol(obj))
skeletal_muscle[which(obj$wsnn_res.0.8 %in% c(6))]="EC2"
skeletal_muscle[which(obj$wsnn_res.0.8 %in% c(2,8,11))]="EC1"
skeletal_muscle[which(obj$wsnn_res.0.8 %in% c(20))]="EC3"
skeletal_muscle[which(obj$wsnn_res.0.8 %in% c(7,19))]="MSM"
skeletal_muscle[which(obj$wsnn_res.0.8 %in% c(5))]="MuSC"
skeletal_muscle[which(obj$wsnn_res.0.8 %in% c(9))]="Macrophage"
skeletal_muscle[which(obj$wsnn_res.0.8 %in% c(3))]="Pericyte"
skeletal_muscle[which(obj$wsnn_res.0.8 %in% c(14))]="SMC1"
skeletal_muscle[which(obj$wsnn_res.0.8 %in% c(10))]="SMC2"
skeletal_muscle[which(obj$wsnn_res.0.8 %in% c(15,16))]="B/T/NK"
skeletal_muscle[which(obj$wsnn_res.0.8 %in% c(1,4,18))]="FAP1"
skeletal_muscle[which(obj$wsnn_res.0.8 %in% c(12))]="FAP2"
skeletal_muscle[which(obj$wsnn_res.0.8 %in% c(13))]="FAP3"
skeletal_muscle[which(obj$wsnn_res.0.8 %in% c(17))]="FAP4"
obj$skeletal_muscle <- skeletal_muscle

cell_type <- vector(length=ncol(obj))
cell_type[which(obj$wsnn_res.0.8 %in% c(2,6,8,11,20))]="ECs"
cell_type[which(obj$wsnn_res.0.8 %in% c(7,19))]="MSM"
cell_type[which(obj$wsnn_res.0.8 %in% c(5))]="MuSC"
cell_type[which(obj$wsnn_res.0.8 %in% c(9))]="Macrophage"
cell_type[which(obj$wsnn_res.0.8 %in% c(3))]="Pericyte"
cell_type[which(obj$wsnn_res.0.8 %in% c(10,14))]="SMCs"
cell_type[which(obj$wsnn_res.0.8 %in% c(15,16))]="B/T/NK"
cell_type[which(obj$wsnn_res.0.8 %in% c(1,4,12,13,17,18))]="FAPs"
obj$cell_type <- cell_type

# Plot with annotation
my_colors_21 <- c(
  "#FFD92F", # yellow
  "#377EB8", # blue
  "#E78AC3", # magenta
  "#7570B3", # muted violet
  "#984EA3", # purple
  "#FF7F00", # orange
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
  label.size = 2.5
) + NoLegend() + scale_color_manual(values = plot_colors)
f2 <- DimPlot(
  obj,
  reduction = "wnn.umap",
  group.by  = "skeletal_muscle",
  alpha = 0.4,
  label = T,
  label.size = 2.5
) + NoLegend() + scale_color_manual(values = plot_colors)
f <- (f1|f2)
save_gg(f, "Clusters_with_annotation.png", w=11, h=5)

# ========================
# Save the object
# ========================
saveRDS(obj,"../outputs/decont_merged_filt_nodoub_cc_sct_reduc_clust_integ_annot_obj.rds")

# ============================================================================ #

# ========================
# Load the object
# ========================
obj <- readRDS("../outputs/decont_merged_filt_nodoub_cc_sct_reduc_clust_integ_annot_obj.rds")
# obj <- obj %>% subset(cell_type == 'FAPs')

# ========================
# Differential expression analysis (pseudobulk with voom+limma)
# ========================

DefaultAssay(obj) <- "RNA"

run_voomlimma <- function(L) {
  message("voomlimma")
    dge <- DGEList(L$count, group = L$condt)
    dge <- calcNormFactors(dge)
    design <- model.matrix(~ L$condt) # coefficient is L$condtAged if levels = c("Young","Aged")
    vm <- voom(dge, design = design, plot = FALSE) # set to TRUE if you want the plot
    fit <- lmFit(vm, design = design)
    fit <- eBayes(fit)
    tt <- topTable(fit, n = Inf, adjust.method = "BH")
  
  list(tt = tt,
       df = data.frame(pval = tt$P.Value,
                       padj = tt$adj.P.Val,
                       row.names = rownames(tt)))
}

# =====  Pseudobulk aggregation (counts) =====
# Using return.seurat = FALSE to get a plain matrix for limma/edgeR
pb_list <- AggregateExpression(
  obj,
  assays = "RNA",
  group.by = c("sample", "condition", "skeletal_muscle"),
  slot = "counts",
  return.seurat = FALSE
)

# pb_counts is a genes x groups matrix (RNA assay)
pb_counts <- pb_list$RNA

# =====  Build group metadata from column names =====
sp <- str_split_fixed(colnames(pb_counts), "_", 3)  # adjust separator if needed
meta <- data.frame(
  sample = sp[,1],
  condition = factor(sp[,2], levels = c("Young","Aged")), # Young as reference
  celltype = sp[,3],
  stringsAsFactors = FALSE
)
rownames(meta) <- colnames(pb_counts)

# =====  Run voom+limma per cell type (Aged vs Young) =====
subtypes <- sort(unique(meta$celltype))

res_list <- lapply(subtypes, function(st) {
  message("Running voom+limma for ", st)
  
  # Subset columns for this cell type
  keep_cols <- meta$celltype == st
  counts_ct <- pb_counts[, keep_cols, drop = FALSE]
  cond_ct <- droplevels(meta$condition[keep_cols])
  
  # Ensure reference level is Young (design uses ~ cond; coef will be Aged vs Young)
  cond_ct <- relevel(cond_ct, ref = "Young")
  
  # (Optional) filter lowly expressed features
  dge_tmp <- DGEList(counts_ct)
  keep_genes <- filterByExpr(dge_tmp, group = cond_ct, min.count = 10) # conservative
  counts_ct <- counts_ct[keep_genes, , drop = FALSE]
  
  # Skip if not enough samples per group
  if (length(unique(cond_ct)) < 2 || any(table(cond_ct) < 2)) {
    warning("Skipping ", st, " due to insufficient replicates per condition.")
    return(NULL)
  }
  
  # Run voom+limma
  L <- list(count = counts_ct, condt = cond_ct)
  out <- run_voomlimma(L)
  
  # Prepare a Seurat-like marker table layout
  tt <- out$tt
  tt$gene <- rownames(tt)
  tt$celltype <- st
  
  # Map names to the previous pipeline for downstream compatibility
  # limma's logFC is log2(Aged/Young) because coef is condtAged vs Young
  tt <- tt %>%
    transmute(
      gene = gene,
      p_val = round(P.Value, 4),
      p_val_adj = round(adj.P.Val, 4),
      avg_log2FC = round(logFC, 2),
      celltype = celltype
    )
  
  tt
})

# Drop NULLs (cell types with insufficient replicates)
res_list <- Filter(Negate(is.null), res_list)

# =====  Combine and rank results =====
bulk_de <- bind_rows(res_list) %>%
  arrange(p_val_adj, p_val)

# Significant DE genes 
bulk_de_sig <- bulk_de %>% filter(p_val_adj < 0.05 & avg_log2FC >= 1)

wcsv(bulk_de_sig, file.path("results_sc", "bulk_de_sig.csv"))
wcsv(bulk_de_sig, file.path("results_sc", "bulk_de_sig_faps.csv"))


