# ===== Clean environment =====
rm(list = ls(all.names = TRUE))
gc()

# ===== Loading relevant libraries =====
library(arrow)  # read_parquet
library(data.table)
library(dplyr)
library(ggplot2)
library(viridis)
library(magick) # morphology image
library(patchwork)
library(spatstat.geom) # point patterns (ppp)
library(spatstat.explore) # K-functions

# ============================================================================ #

# ===== Paths and gene sets =====
# Base directory that contains GSM* folders
base_dir <- "C:/Users/Hani/Desktop/Hernia/data/GSE288662/Spatial"

# 4 samples: 2 Veh + 2 EP
sample_info <- data.frame(
  sample = c("1450Veh", "1603Veh", "1460EP", "1461EP"),
  folder = c("GSM8772705_1450Veh",
             "GSM8772706_1603Veh",
             "GSM8772707_1460EP",
             "GSM8772708_1461EP"),
  stringsAsFactors = FALSE
)

# Paths to morphology PNGs
morph_paths <- c(
  "1450Veh" = file.path(base_dir, "GSM8772705_1450Veh", "morphology_simple.png"),
  "1603Veh" = file.path(base_dir, "GSM8772706_1603Veh", "morphology_simple.png"),
  "1460EP"  = file.path(base_dir, "GSM8772707_1460EP",  "morphology_simple.png"),
  "1461EP"  = file.path(base_dir, "GSM8772708_1461EP",  "morphology_simple.png")
)

# Smad3 TF-network targets (mouse gene symbols)
smad3_gene <- "Smad3"   # main TF
smad3_pos_genes <- c(
  "Gprc5a","Hivep1","Eml4","Agtr1a",
  "Medag","Galnt2","Hivep2","Spag17",
  "Clstn2","Glis3","Cxcl2","Trio"
)
smad3_neg_genes <- c("Wdfy1","Il34","Mtm1","Emb")

# ===== Helper: load transcripts.parquet for one sample =====
load_transcripts <- function(sample_row,
                             qv_min = 20,
                             keep_only_genes = NULL) {
  # sample_row is one row from sample_info (sample, folder)
  
  samp   <- sample_row$sample
  folder <- file.path(base_dir, sample_row$folder)
  parquet_path <- file.path(folder, "transcripts.parquet")
  
  message("Loading transcripts for sample: ", samp)
  tx <- as.data.frame(arrow::read_parquet(parquet_path))
  
  # Keep only true gene molecules if column exists
  if ("is_gene" %in% colnames(tx)) {
    tx <- tx[tx$is_gene == TRUE, , drop = FALSE]
  }
  
  # Quality filter
  if ("qv" %in% colnames(tx)) {
    tx <- tx[tx$qv >= qv_min, , drop = FALSE]
  }
  
  # Optionally restrict to selected genes
  if (!is.null(keep_only_genes)) {
    tx <- tx[tx$feature_name %in% keep_only_genes, , drop = FALSE]
  }
  
  tx$sample <- samp
  
  # Rename to generic columns
  stopifnot(all(c("feature_name","x_location","y_location") %in% colnames(tx)))
  tx <- tx %>%
    dplyr::rename(
      gene = feature_name,
      x    = x_location,
      y    = y_location
    )
  
  return(tx)
}

# ===== Load only Smad3 + its TF-network targets (all 4 samples) =====
genes_of_interest <- unique(c(smad3_gene,
                              smad3_pos_genes,
                              smad3_neg_genes))

tx_list <- list()
for (i in seq_len(nrow(sample_info))) {
  tx_list[[ sample_info$sample[i] ]] <-
    load_transcripts(sample_info[i, ],
                     qv_min = 20,
                     keep_only_genes = genes_of_interest)
}

# Combine all samples in one data.frame
tx_all <- data.table::rbindlist(tx_list, use.names = TRUE, fill = TRUE)

# Define groups: Smad3, Pos targets, Neg targets
tx_all$group <- "Other"
tx_all$group[ tx_all$gene == smad3_gene ] <- "Smad3"
tx_all$group[ tx_all$gene %in% smad3_pos_genes ] <- "SMAD3_Pos"
tx_all$group[ tx_all$gene %in% smad3_neg_genes ] <- "SMAD3_Neg"
tx_all$group <- factor(tx_all$group,
                       levels = c("Smad3","SMAD3_Pos","SMAD3_Neg","Other"))

# Basic counts per sample × group
count_table <- tx_all %>%
  count(sample, group) %>%
  tidyr::pivot_wider(names_from = group, values_from = n, values_fill = 0)
write.csv(count_table, "Transcript_counts_Smad3_network_by_sample.csv",
          row.names = FALSE)

# Interpretation:
# This table shows how many Smad3 molecules and target transcripts
# are detected in each sample (2 Veh vs 2 EP). Strong differences
# in counts suggest global up/down regulation.

# ===== Point plots + 2D density for each sample =====
plot_transcripts_on_image <- function(df_sample, sample_name) {
  # df_sample: transcripts for one sample (Smad3 + targets)
  
  # Load morphology image
  morph_path <- morph_paths[[sample_name]]
  img <- image_read(morph_path)
  rst <- as.raster(img)
  
  xmin <- min(df_sample$x); xmax <- max(df_sample$x)
  ymin <- min(df_sample$y); ymax <- max(df_sample$y)
  
  # (A) Point plot: global + Smad3/Pos/Neg
  p_points <- ggplot() +
    annotation_raster(rst, xmin = xmin, xmax = xmax,
                      ymin = ymin, ymax = ymax) +
    geom_point(
      data = df_sample,
      aes(x = x, y = y),
      color = "grey80",
      size  = 0.1,
      alpha = 0.4
    ) +
    geom_point(
      data = df_sample %>%
        filter(group %in% c("Smad3","SMAD3_Pos","SMAD3_Neg")),
      aes(x = x, y = y, color = group),
      size  = 0.35,
      alpha = 0.8
    ) +
    scale_color_manual(
      values = c("Smad3" = "gold",
                 "SMAD3_Pos" = "red",
                 "SMAD3_Neg" = "dodgerblue",
                 "Other" = "grey80")
    ) +
    scale_y_reverse() +
    coord_equal() +
    theme_void(base_size = 11) +
    ggtitle(paste0(sample_name,
                   " — Smad3 network transcripts on morphology")) +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5, size = 11))
  # Interpretation:
  # Yellow: Smad3, red: positive targets, blue: negative targets.
  # Co-local clusters of yellow+red suggest local network activation.
  
  # (B) 2D density only for SMAD3_Pos
  df_pos <- df_sample %>% filter(group == "SMAD3_Pos")
  
  p_density_pos <- ggplot() +
    annotation_raster(rst, xmin = xmin, xmax = xmax,
                      ymin = ymin, ymax = ymax) +
    stat_density_2d(
      data = df_pos,
      aes(x = x, y = y, fill = after_stat(level)),
      geom = "polygon",
      alpha = 0.5,
      contour = TRUE
    ) +
    scale_fill_viridis(option = "magma") +
    geom_point(
      data = df_pos,
      aes(x = x, y = y),
      color = "white",
      size = 0.2,
      alpha = 0.8
    ) +
    scale_y_reverse() +
    coord_equal() +
    theme_void(base_size = 11) +
    ggtitle(paste0(sample_name,
                   " — density of SMAD3_Pos transcripts")) +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5, size = 11))
  # Interpretation:
  # Filled contours highlight spatial hotspots of positive targets only.
  # Comparing 2 Veh vs 2 EP reveals where Smad3-driven activation emerges.
  
  list(points = p_points, density_pos = p_density_pos)
}

plots_per_sample <- list()
for (s in sample_info$sample) {
  df_s <- tx_all %>% filter(sample == s)
  plots_per_sample[[s]] <- plot_transcripts_on_image(df_s, s)
}

# Arrange point plots in 2×2
p_points_grid <- (plots_per_sample[["1450Veh"]]$points |
                    plots_per_sample[["1603Veh"]]$points) /
  (plots_per_sample[["1460EP"]]$points  |
     plots_per_sample[["1461EP"]]$points)

ggsave("Xenium_transcripts_Smad3_points_4samples.png",
       p_points_grid, width = 12, height = 8, dpi = 300)

# Arrange density maps in 2×2
p_density_grid <- (plots_per_sample[["1450Veh"]]$density_pos |
                     plots_per_sample[["1603Veh"]]$density_pos) /
  (plots_per_sample[["1460EP"]]$density_pos  |
     plots_per_sample[["1461EP"]]$density_pos)

ggsave("Xenium_transcripts_Smad3Pos_density_4samples.png",
       p_density_grid, width = 12, height = 8, dpi = 300)

# ===== Cross-K function: Smad3 vs Pos targets (each sample) =====
# We compute Kcross(Smad3, PosTarget) in each sample separately,
# to test spatial attraction vs independence between TF and targets.

for (s in sample_info$sample) {
  df_s <- tx_all %>%
    filter(sample == s,
           group %in% c("Smad3","SMAD3_Pos"))
  
  # If one of the marks is missing, skip
  if (sum(df_s$group == "Smad3") == 0 ||
      sum(df_s$group == "SMAD3_Pos") == 0) {
    message("Skipping Kcross for sample ", s,
            " (Smad3 or Pos targets not detected).")
    next
  }
  
  win_s <- owin(
    xrange = range(df_s$x),
    yrange = range(df_s$y)
  )
  
  pp_s <- ppp(
    x = df_s$x,
    y = df_s$y,
    window = win_s,
    marks = factor(ifelse(df_s$group == "Smad3",
                           "Smad3", "PosTarget"))
  )
  
  K_cross_s <- Kcross(pp_s, i = "Smad3", j = "PosTarget")
  
  png(paste0("Kcross_Smad3_PosTargets_", s, ".png"),
      width = 900, height = 700, res = 300)
  plot(K_cross_s,
       main = paste0("Kcross(Smad3, Pos targets) — ", s))
  dev.off()
  
  saveRDS(K_cross_s,
          file = paste0("Kcross_Smad3_PosTargets_", s, ".rds"))
  
  # Interpretation:
  # If Kcross(r) > theoretical pi*r^2 curve across distances r,
  # Smad3 and its positive targets are spatially attracted at scale r.
  # EP samples showing stronger attraction vs Veh would support
  # treatment-specific co-localization of TF and its targets.
}
