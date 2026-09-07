# ============================================================================
# Reviewer C7:
# Regenerate Figure 18A-C candidate panels using the corrected full-18
# CellChat networks (nboot=1000, plus-one empirical P, condition-wide BH).
#
# Historical plotting conventions are preserved.
# No canonical figure is modified.
# ============================================================================

rm(list = ls(all.names = TRUE))

suppressPackageStartupMessages({
  library(CellChat)
  library(ggplot2)
  library(patchwork)
  library(ComplexHeatmap)
  library(grid)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (length(file_arg) == 1L) {
  script_dir <- dirname(
    normalizePath(sub("^--file=", "", file_arg))
  )
  repo_root <- normalizePath(file.path(script_dir, ".."))
} else {
  repo_root <- normalizePath(".")
}

candidate_dir <- file.path(
  repo_root,
  "processed_results",
  "12_CellChat",
  "multiple_testing_audit",
  "full18_nboot1000_BH",
  "figure18ABC_corrected_candidate"
)

figdir <- file.path(
  candidate_dir,
  "plots"
)

dir.create(
  figdir,
  recursive = TRUE,
  showWarnings = FALSE
)

young_file <- file.path(
  candidate_dir,
  "cellchat_Young_full18_nboot1000_plus1BH_corrected_candidate.rds"
)

aged_file <- file.path(
  candidate_dir,
  "cellchat_Aged_full18_nboot1000_plus1BH_corrected_candidate.rds"
)

for (f in c(young_file, aged_file)) {
  if (!file.exists(f)) {
    stop("Required corrected CellChat object not found: ", f)
  }
}

cellchat_young <- readRDS(young_file)
cellchat_aged  <- readRDS(aged_file)

object.list <- list(
  Young = cellchat_young,
  Aged = cellchat_aged
)

# ============================================================================
# Figure 18A: overall signaling circles
# ============================================================================

plot_circle_png <- function(
    cellchat_obj,
    title_str,
    outfile
) {

  mat <- cellchat_obj@net$weight

  group_order <- c(
    "FAP3",
    "Other FAPs",
    "MuSC",
    "Myogenic",
    "Immune",
    "Endothelial",
    "Vascular stromal"
  )

  go <- intersect(
    group_order,
    rownames(mat)
  )

  mat_reorder <- mat[
    go,
    go,
    drop = FALSE
  ]

  vertex.weight <- rowSums(
    mat_reorder
  )

  png(
    outfile,
    width = 1800,
    height = 1800,
    res = 300
  )

  netVisual_circle(
    mat_reorder,
    vertex.weight = vertex.weight,
    weight.scale = TRUE,
    label.edge = FALSE,
    edge.weight.max = max(mat_reorder),
    vertex.label.cex = 1.4,
    vertex.size.max = 15
  )

  grid::grid.text(
    title_str,
    y = unit(0.97, "npc"),
    gp = grid::gpar(
      fontsize = 18,
      fontface = "bold"
    )
  )

  dev.off()
}

plot_circle_png(
  cellchat_young,
  "Young - overall signaling",
  file.path(
    figdir,
    "Figure18A_circle_overall_Young_corrected_candidate.png"
  )
)

plot_circle_png(
  cellchat_aged,
  "Aged - overall signaling",
  file.path(
    figdir,
    "Figure18A_circle_overall_Aged_corrected_candidate.png"
  )
)

# ============================================================================
# Figure 18B: signaling-role scatter
# ============================================================================

num.link <- sapply(
  object.list,
  function(z) {
    rowSums(z@net$count) +
      colSums(z@net$count) -
      diag(z@net$count)
  }
)

weight.MinMax <- c(
  min(num.link),
  max(num.link)
)

p_scatter_list <- lapply(
  names(object.list),
  function(nm) {

    netAnalysis_signalingRole_scatter(
      object = object.list[[nm]],
      slot.name = "netP",
      title = nm,
      weight.MinMax = weight.MinMax
    )
  }
)

p_scatter_combined <-
  p_scatter_list[[1]] +
  p_scatter_list[[2]]

ggsave(
  filename = file.path(
    figdir,
    "Figure18B_signalingRole_scatter_Young_vs_Aged_corrected_candidate.png"
  ),
  plot = p_scatter_combined,
  width = 10,
  height = 5,
  dpi = 300
)

ggsave(
  filename = file.path(
    figdir,
    "Figure18B_signalingRole_scatter_Young_vs_Aged_corrected_candidate.pdf"
  ),
  plot = p_scatter_combined,
  width = 10,
  height = 5
)

# ============================================================================
# Figure 18C: signaling-role heatmaps
# ============================================================================

pathway.union <- union(
  cellchat_young@netP$pathways,
  cellchat_aged@netP$pathways
)

cat(
  "Corrected pathway union:",
  paste(pathway.union, collapse = ", "),
  "\n"
)

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

png(
  file.path(
    figdir,
    "Figure18C_incoming_netP_Young_vs_Aged_corrected_candidate.png"
  ),
  width = 3200,
  height = 1600,
  res = 300
)

ComplexHeatmap::draw(
  ht_in_young + ht_in_aged,
  ht_gap = unit(0.5, "cm")
)

dev.off()

pdf(
  file.path(
    figdir,
    "Figure18C_incoming_netP_Young_vs_Aged_corrected_candidate.pdf"
  ),
  width = 3200 / 300,
  height = 1600 / 300
)

ComplexHeatmap::draw(
  ht_in_young + ht_in_aged,
  ht_gap = unit(0.5, "cm")
)

dev.off()

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

png(
  file.path(
    figdir,
    "Figure18C_outgoing_netP_Young_vs_Aged_corrected_candidate.png"
  ),
  width = 3200,
  height = 1600,
  res = 300
)

ComplexHeatmap::draw(
  ht_out_young + ht_out_aged,
  ht_gap = unit(0.5, "cm")
)

dev.off()

pdf(
  file.path(
    figdir,
    "Figure18C_outgoing_netP_Young_vs_Aged_corrected_candidate.pdf"
  ),
  width = 3200 / 300,
  height = 1600 / 300
)

ComplexHeatmap::draw(
  ht_out_young + ht_out_aged,
  ht_gap = unit(0.5, "cm")
)

dev.off()

cat("\n============================================================\n")
cat("FIGURE 18A-C CORRECTED CANDIDATES GENERATED\n")
cat("============================================================\n\n")

cat("Output directory:\n")
cat(figdir, "\n")

cat(
  "\nNo canonical Figure 18 file was modified.\n"
)
