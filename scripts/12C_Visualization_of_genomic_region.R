# ========================
# Setting up environment
# ========================
# ===== Loading relevant libraries =====
library(Signac)
library(ggplot2)
library(Seurat)
library(BSgenome.Hsapiens.UCSC.hg38)

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
obj <- readRDS("../outputs/hdWGCNA_TFNet_DEReg_L2G_obj.rds")

# ========================
# Motif footprinting
# ========================
# Gather the footprinting information for sets of motifs
obj <- Footprint(
  object = obj,
  motif.name = c("ETS2", "ETS1", "FOS"),
  genome = BSgenome.Hsapiens.UCSC.hg38,
  in.peaks = T
)

# Plot the footprint data for each group of cells
p1 <- PlotFootprint(obj, features = c("ETS2", "ETS1", "FOS"),label.top = 4)
p2 <- p1 + patchwork::plot_layout(ncol = 1)
save_gg(p2, "PlotFootprint.png", w=5, h=8)

# ========================
# Plotting aggregated signal
# ========================
# CoveragePlot
roi = "chr15-67065602-67195169"

# DA peaks overlapping gene of interest
open_fap3 <- c("chr15-67109227-67110381","chr15-67198215-67198971")
regions_highlight <- subsetByOverlaps(StringToGRanges(open_fap3), LookupGeneCoords(obj, "SMAD3"))

p3 <- CoveragePlot(
  object = obj,
  region = "SMAD3",
  features = "SMAD3",
  annotation = TRUE,
  peaks = TRUE,
  tile = TRUE,
  links = TRUE,
  expression.assay = "SCT",
  region.highlight = regions_highlight,
  extend.upstream = 1000,
  extend.downstream = 1000
)

save_gg(p3, "genomic_tracks.png", w=8, h=6)

# ========================
# Save the object
# ========================

saveRDS(obj,"../outputs/hdWGCNA_TFNet_DEReg_L2G_Foot_obj.rds")










