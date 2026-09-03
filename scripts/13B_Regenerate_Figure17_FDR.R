library(data.table)
library(dplyr)
library(ggplot2)
library(gghalves)
library(ggpubr)

# -------------------------------------------------------------------------
# Repository paths
# -------------------------------------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (length(file_arg) == 1) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg)))
  repo_root <- normalizePath(file.path(script_dir, ".."))
} else {
  repo_root <- normalizePath(".")
}

source_dir <- file.path(
  repo_root,
  "source_data",
  "figure_source_data",
  "Figure_17"
)

figdir <- file.path(
  repo_root,
  "outputs",
  "chromVAR",
  "Raincloud"
)

dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# Load canonical Figure 17 source data
# -------------------------------------------------------------------------
df_plot <- fread(
  file.path(
    source_dir,
    "Figure17_chromVAR_cell_level_source_data.csv"
  )
)

stats_fdr <- fread(
  file.path(
    source_dir,
    "Figure17_Wilcoxon_Young_vs_Aged.csv"
  )
)

df_plot[, condition := factor(condition, levels = c("Young", "Aged"))]
df_plot[, TF := factor(TF, levels = c("ETS2", "ETS1", "FOS", "SMAD3"))]
df_plot[, FAP := factor(FAP, levels = c("FAP1", "FAP2", "FAP3", "FAP4"))]

stats_fdr[, TF := factor(TF, levels = c("ETS2", "ETS1", "FOS", "SMAD3"))]
stats_fdr[, FAP := factor(FAP, levels = c("FAP1", "FAP2", "FAP3", "FAP4"))]

# -------------------------------------------------------------------------
# FDR-based annotations
# -------------------------------------------------------------------------
yrange <- df_plot[
  ,
  .(
    y_max = max(value, na.rm = TRUE),
    y_min = min(value, na.rm = TRUE)
  ),
  by = .(FAP, TF)
]

stats_fdr <- merge(
  stats_fdr,
  yrange,
  by = c("FAP", "TF"),
  all.x = TRUE
)

stats_fdr[
  ,
  `:=`(
    y.position = y_max + 0.08 * (y_max - y_min),
    group1 = "Young",
    group2 = "Aged"
  )
]

stats_fdr_sig <- stats_fdr[
  p_adj_significance != "ns"
]

# -------------------------------------------------------------------------
# Figure 17
# -------------------------------------------------------------------------
p_facet <- ggplot(
  df_plot,
  aes(x = condition, y = value, fill = condition)
) +
  geom_half_violin(
    side = "l",
    alpha = 0.6,
    trim = FALSE
  ) +
  geom_half_boxplot(
    side = "r",
    width = 0.25,
    alpha = 1,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.1,
    alpha = 0.25,
    size = 0.6
  ) +
  ggpubr::stat_pvalue_manual(
    stats_fdr_sig,
    label = "p_adj_significance",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01,
    size = 3,
    inherit.aes = FALSE
  ) +
  facet_grid(
    TF ~ FAP,
    scales = "free_y"
  ) +
  scale_fill_manual(
    values = c("#1f78b4", "#e31a1c")
  ) +
  theme_bw() +
  labs(
    y = "chromVAR deviation",
    x = "",
    title = "chromVAR Motif Deviation Scores: Raincloud (FAP1-4 x ETS2/ETS1/FOS/SMAD3)"
  )

ggsave(
  file.path(
    figdir,
    "Raincloud_FACET_FAP1_4_all_TFs_FDR.png"
  ),
  p_facet,
  width = 15,
  height = 12,
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(
    figdir,
    "Raincloud_FACET_FAP1_4_all_TFs_FDR.pdf"
  ),
  p_facet,
  width = 15,
  height = 12,
  bg = "white"
)

cat("\nFigure 17 regenerated with BH-FDR annotations.\n")
cat("PNG:", file.path(figdir, "Raincloud_FACET_FAP1_4_all_TFs_FDR.png"), "\n")
cat("PDF:", file.path(figdir, "Raincloud_FACET_FAP1_4_all_TFs_FDR.pdf"), "\n")
