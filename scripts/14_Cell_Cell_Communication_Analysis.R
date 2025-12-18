# ========================
# Setting up environment
# ========================
# ===== Clean environment =====
rm(list = ls(all.names = TRUE))
gc()

# ===== Loading relevant libraries =====
library(Seurat)
library(ggplot2)
library(CellChat)

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

# ============================================================================ #

# ========================
# Load the object
# ========================
obj <- readRDS("../outputs/decont_merged_filt_nodoub_cc_sct_reduc_clust_integ_annot_obj.rds")

# ========================
# Cell-cell communication analysis
# ========================
options(future.globals.maxSize = 8 * 1024^3)
future::plan("multisession", workers = 4) # do parallel
obj$samples <- obj$sample
obj$cc_group <- NA_character_
obj$cc_group[obj$skeletal_muscle %in% c("FAP1","FAP2","FAP4")] <- "Other FAPs"
obj$cc_group[obj$skeletal_muscle == "FAP3"] <- "FAP3"
obj$cc_group[obj$skeletal_muscle == "MuSC"] <- "MuSC"
obj$cc_group[obj$skeletal_muscle == "MSM"] <- "Myogenic"
obj$cc_group[obj$skeletal_muscle %in% c("Macrophage","B/T/NK")] <- "Immune"
obj$cc_group[obj$skeletal_muscle %in% c("EC1","EC2","EC3")] <- "Endothelial"
obj$cc_group[obj$skeletal_muscle %in% c("Pericyte","SMC1","SMC2")] <- "Vascular stromal"

obj_young <- subset(obj, subset = condition == "Young")
obj_aged  <- subset(obj, subset = condition == "Aged")
rm(obj)

obj_young$group_cellchat <- obj_young$cc_group
obj_aged$group_cellchat  <- obj_aged$cc_group

make_cellchat <- function(seu) {
  data.input <- GetAssayData(seu, assay = "SCT", slot = "data")
  meta <- seu@meta.data
  
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group_cellchat", assay = "SCT")
  
  CellChatDB <- CellChatDB.human
  
  pathways.use <- c(
    # SMAD3 / fibrosis axis
    "TGFb",
    
    # ECM structural components
    "COLLAGEN",
    "FN1",
    "LAMININ",
    "THBS",   # thrombospondin
    "TENASCIN",
    "VTN",  # vitronectin
    "HSPG", # heparan sulfate proteoglycan
    
    # Cell–cell / cell–matrix adhesion
    "CDH", "CDH1", "CDH5", "PCDH",  # cadherins
    "ICAM", "VCAM", "JAM",
    
    # ECM remodeling
    "MMP", # metalloproteinases
    
    # FAP / fibro-inflammatory markers
    "SPP1",        
    "PERIOSTIN"   
  )
  
  CellChatDB.use <- subsetDB(CellChatDB, search = pathways.use, key = "pathway_name")
  
  cellchat@DB <- CellChatDB.use
  
  cellchat <- subsetData(cellchat)
  
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  cellchat <- computeCommunProb(
    cellchat,
    type  = "triMean")
  
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  
  return(cellchat)
}

cellchat_young <- make_cellchat(obj_young)
cellchat_aged <- make_cellchat(obj_aged)

object.list <- list(Young = cellchat_young, Aged = cellchat_aged)
cellchat_merged <- mergeCellChat(object.list = object.list, add.names = names(object.list))

# Visualization
# pathways
pathway_TGF <- "TGFb"

all_pw <- unique(CellChatDB.human$interaction$pathway_name)

# Circle plot (overall signalling)
plot_circle_png <- function(cellchat_obj, title_str) {
  mat <- cellchat_obj@net$weight
  group_order <- c("FAP3","Other FAPs","MuSC","Myogenic","Immune","Endothelial","Vascular stromal")
  go <- intersect(group_order, rownames(mat))
  mat_reorder <- mat[go, go, drop = FALSE]
  vertex.weight <- rowSums(mat_reorder)
  
  netVisual_circle(
    mat_reorder,
    vertex.weight = vertex.weight,
    weight.scale = TRUE,
    label.edge  = FALSE,
    edge.weight.max = max(mat_reorder),
    vertex.label.cex = 1.4,
    vertex.size.max  = 15
  )
  grid::grid.text(title_str, y = unit(0.97, "npc"),
                  gp = grid::gpar(fontsize = 18, fontface = "bold"))
}

png(file.path(figdir, "circle_overall_Young.png"),
    width = 1800, height = 1800, res = 300)
plot_circle_png(cellchat_young, "Young - overall signaling")
dev.off()

png(file.path(figdir, "circle_overall_Aged.png"),
    width = 1800, height = 1800, res = 300)
plot_circle_png(cellchat_aged,  "Aged - overall signaling")
dev.off()

# Circos plot (pathways)
png(file.path(figdir, "circle_TGFb_Young.png"),
    width = 2000, height = 2000, res = 300)

netVisual_aggregate(
  object = cellchat_young,
  signaling = pathway_TGF,
  layout = "circle",
  thresh = 0.05
)

dev.off()

png(file.path(figdir, "circle_TGFb_Aged.png"),
    width = 2000, height = 2000, res = 300)

netVisual_aggregate(
  object = cellchat_aged,
  signaling = pathway_TGF,
  layout = "circle",
  thresh = 0.05
)

dev.off()

# Bubble plot
signaling.use <- pathway_TGF

p_bubble_FAP <- netVisual_bubble(
  object  = cellchat_merged,
  sources.use = c("FAP3","Other FAPs"),
  targets.use = c("MuSC","Myogenic","Immune","Endothelial","Vascular stromal"),
  signaling = signaling.use,
  angle.x = 45,
  remove.isolate = TRUE,
  comparison = c(1,2)
)

save_gg(
  p_bubble_FAP, "bubble_FAP_TGFb_Young_vs_Aged.png",
  w = 8, h = 5
)

# Signaling role scatter
# Common range for point size (based on number of links)
num.link <- sapply(object.list, function(x) {
  rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)
})
weight.MinMax <- c(min(num.link), max(num.link))

p_scatter_list <- lapply(names(object.list), function(nm) {
  netAnalysis_signalingRole_scatter(
    object  = object.list[[nm]],
    slot.name = "netP",
    title = nm,
    weight.MinMax = weight.MinMax
  )
})

p_scatter_combined <- p_scatter_list[[1]] + p_scatter_list[[2]]

save_gg(
  p_scatter_combined, "signalingRole_scatter_Young_vs_Aged.png",
  w = 10, h = 5
)

# Signaling heatmap
# Pathway union so that Young/Aged are comparable
pathway.union <- union(cellchat_young@netP$pathways,
                       cellchat_aged@netP$pathways)

# Outgoing 
ht_out_young <- netAnalysis_signalingRole_heatmap(
  object = cellchat_young,
  pattern = "outgoing",
  slot.name = "netP",
  signaling = pathway.union,
  title = "Young"
)

ht_out_aged <- netAnalysis_signalingRole_heatmap(
  object = cellchat_aged,
  pattern = "outgoing",
  slot.name = "netP",
  signaling = pathway.union,
  title = "Aged"
)

png(file.path(figdir, "heatmap_outgoing_netP_Young_vs_Aged.png"),
    width = 3200, height = 1600, res = 300)
ComplexHeatmap::draw(ht_out_young + ht_out_aged, ht_gap = unit(0.5, "cm"))
dev.off()

# Incoming 
ht_in_young <- netAnalysis_signalingRole_heatmap(
  object = cellchat_young,
  pattern = "incoming",
  slot.name = "netP",
  signaling = pathway.union,
  title = "Young"
)

ht_in_aged <- netAnalysis_signalingRole_heatmap(
  object = cellchat_aged,
  pattern = "incoming",
  slot.name = "netP",
  signaling = pathway.union,
  title = "Aged"
)

png(file.path(figdir, "heatmap_incoming_netP_Young_vs_Aged.png"),
    width = 3200, height = 1600, res = 240)
ComplexHeatmap::draw(ht_in_young + ht_in_aged, ht_gap = unit(0.5, "cm"))
dev.off()

# Contribution plots
# TGFb in Aged
p_contrib_TGF_aged <- netAnalysis_contribution(
  object = cellchat_aged,
  signaling = pathway_TGF
)
save_gg(
  p_contrib_TGF_aged, "contribution_TGFb_Aged.png",
  w = 7, h = 5
)

saveRDS(cellchat_merged, file = "../outputs/cellchat_merged.rds")
