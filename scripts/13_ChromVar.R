# ===== Clean environment =====
rm(list = ls(all.names = TRUE))
gc()

# ===== Loading relevant libraries =====
library(future)
plan("sequential")
library(BiocParallel)
register(SerialParam())
library(Signac)
library(JASPAR2024)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)
library(ggplot2)
library(Seurat)
library(chromVAR)
library(gghalves)
library(ggpubr)
library(tibble)
library(ggplot2)

# ===== Helpers =====
figdir <- "results_TF/Raincloud"
save_gg <- function(p, filename, w=7, h=5, dpi=300){
  if (inherits(p, "ggplot") || inherits(p, "patchwork")){
    ggsave(file.path(figdir, filename), p, width=w, height=h, units="in", dpi=dpi, bg="white")
  } else {
    png(file.path(figdir, filename), width=w, height=h, units="in", res=dpi)
    print(p); dev.off()
  }
}

# ============================================================================ #

# ===== chromVAR motif activity =====
# Load object with SCT + ATAC + TFNet
obj <- readRDS("../outputs/hdWGCNA_TFNet_DEReg_L2G_obj.rds")

# Use ATAC assay
DefaultAssay(obj) <- "ATAC"
obj <- FindTopFeatures(obj, min.cutoff = "q0.7")

goi <- "SMAD3"
peak_links <- Links(obj)
peak_links_df <- as.data.frame(peak_links) %>%
  # Harmonize column names
  dplyr::rename(
    target_gene = gene, # Gene symbol linked to this peak
    peak_region = peak  # Peak coordinates "chr-start-end"
  ) %>%
  dplyr::select(
    target_gene,
    peak_region,
    peak_score = score, # Correlation strength / link weight
    peak_z = zscore,  # z-scored strength
    peak_pval = pvalue  # p-value of link
  ) %>%
  # Keep only genes in the list of interest
  dplyr::filter(target_gene %in% goi) %>%
  dplyr::distinct()

# Run chromVAR
obj <- RunChromVAR(
  object = obj,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

saveRDS(obj, "../outputs/hdWGCNA_TFNet_DEReg_L2G_chromVAR_obj.rds")

DefaultAssay(obj) <- "chromvar"

# Find motif IDs for TFs of interest
motif_obj <- Motifs(obj[["ATAC"]])
motif_names <- motif_obj@motif.names
motif_ids <- names(motif_obj@pwm)

get_motif_ids <- function(tf, motif_names, motif_ids) {
  hits <- grep(tf, motif_names, ignore.case = TRUE)
  motif_ids[hits]
}

motifs_smad3 <- get_motif_ids("SMAD3", motif_names, motif_ids)
motifs_ets1 <- get_motif_ids("ETS1",  motif_names, motif_ids)
motifs_ets2 <- get_motif_ids("ETS2",  motif_names, motif_ids)
motifs_fos <- get_motif_ids("FOS",   motif_names, motif_ids)

motifs_smad3; motifs_ets1; motifs_ets2; motifs_fos

# ===== Helper: finding motif in chromvar =====
get_chromvar_feature <- function(seu, motif_id) {
  rn <- rownames(seu[["chromvar"]])
  cand <- c(motif_id, paste0("motif_", motif_id))
  hit  <- cand[cand %in% rn]
  if (length(hit) == 0) {
    warning("No chromVAR feature found for motif_id = ", motif_id)
    return(NA_character_)
  }
  hit[1]
}

# ===== Raincloud function =====
plot_raincloud <- function(df, feat, cf, tf_name, save_dir="results_TF/Raincloud") {
  
  if (is.na(feat)) {
    warning("Skipping ", tf_name, " in ", cf, " : feature is NA")
    return(NULL)
  }
  
  if (!feat %in% colnames(df)) {
    warning("Feature ", feat, " not found in df; skipping ", tf_name, " in ", cf)
    return(NULL)
  }
  
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  
  p <- ggplot(df, aes(
    x = condition,
    y = .data[[feat]],
    fill = condition
  )) +
    geom_half_violin(side = "l", alpha = 0.6, trim = FALSE) +
    geom_half_boxplot(side = "r", width = 0.25, alpha = 1, outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.3, size = 0.7) +
    stat_compare_means(
      method = "wilcox.test",
      label = "p.format",
      size = 3,
      comparisons = list(c("Young","Aged"))
    ) +
    scale_fill_manual(values=c("#1f78b4", "#e31a1c")) +
    theme_bw() +
    labs(
      title = paste0(cf, " — ", tf_name),
      y = "chromVAR deviation",
      x = ""
    )
  
  ggsave(
    filename = file.path(save_dir, paste0("Raincloud_", cf, "_", tf_name, ".png")),
    plot = p,
    width = 5, height = 4, dpi = 300, bg = "white"
  )
  
  return(p)
}

# ===== Loop for FAP1–4: Raincloud plots + DA peaks =====
# Peaks_use for SMAD3
peaks_use <- unique(peak_links_df$peak_region)

# Loop over FAP1–4
fap_levels <- c("FAP1", "FAP2", "FAP3", "FAP4")
results_list <- list()
all_plots_for_facet <- list()

for (cf in fap_levels) {
  
  message("==== Analyzing ", cf, " ====")
  
  # subset 
  obj_cf <- subset(
    obj,
    subset = skeletal_muscle == cf & condition %in% c("Young","Aged")
  )
  if (ncol(obj_cf) == 0) {
    warning("No cells for ", cf, " with Young/Aged; skipping.")
    next
  }
  
  DefaultAssay(obj_cf) <- "chromvar"
  Idents(obj_cf) <- obj_cf$condition
  
  feat_ets2  <- get_chromvar_feature(obj_cf, motifs_ets2[1])
  feat_ets1  <- get_chromvar_feature(obj_cf, motifs_ets1[1])
  feat_fos   <- get_chromvar_feature(obj_cf, motifs_fos[1])
  feat_smad3 <- get_chromvar_feature(obj_cf, motifs_smad3[1])
  
  vars_use <- c("condition", feat_ets2, feat_ets1, feat_fos, feat_smad3)
  vars_use <- vars_use[!is.na(vars_use)]
  
  df_cf <- FetchData(
    obj_cf,
    vars = vars_use
  )
  df_cf$condition <- factor(df_cf$condition, levels = c("Young","Aged"))
  
  # Raincloud ETS2
  p_ets2 <- plot_raincloud(df_cf, feat_ets2, cf, "ETS2")
  # Raincloud ETS1
  p_ets1 <- plot_raincloud(df_cf, feat_ets1, cf, "ETS1")
  # Raincloud FOS
  p_fos  <- plot_raincloud(df_cf, feat_fos,  cf, "FOS")
  # Raincloud SMAD3 motif
  p_smad3 <- plot_raincloud(df_cf, feat_smad3, cf, "SMAD3")
  
  if (!is.na(feat_ets2))
    all_plots_for_facet[[paste0(cf,"_ETS2")]] <- df_cf %>% mutate(TF="ETS2",  feature=feat_ets2,  FAP=cf)
  if (!is.na(feat_ets1))
    all_plots_for_facet[[paste0(cf,"_ETS1")]] <- df_cf %>% mutate(TF="ETS1",  feature=feat_ets1,  FAP=cf)
  if (!is.na(feat_fos))
    all_plots_for_facet[[paste0(cf,"_FOS")]] <- df_cf %>% mutate(TF="FOS",   feature=feat_fos,   FAP=cf)
  if (!is.na(feat_smad3))
    all_plots_for_facet[[paste0(cf,"_SMAD3")]] <- df_cf %>% mutate(TF="SMAD3", feature=feat_smad3, FAP=cf)
  
  # DA peaks for SMAD3
  DefaultAssay(obj_cf) <- "ATAC"
  
  da_peaks_smad3 <- FindMarkers(
    obj_cf,
    ident.1 = "Aged",
    ident.2 = "Young",
    features = peaks_use,
    test.use = "LR",
    latent.vars = "nCount_ATAC"
  )
  
  da_peaks_smad3 <- da_peaks_smad3 %>%
    tibble::rownames_to_column("peak_region") %>%
    left_join(peak_links_df, by = "peak_region")
  
  out_csv <- paste0("results_TF/DA_peaks_linked_to_SMAD3_", cf, "_Aged_vs_Young.csv")
  write.csv(da_peaks_smad3, out_csv, row.names = FALSE)
  
  results_list[[cf]] <- list(
    raincloud = list(
      ETS2 = p_ets2,
      ETS1 = p_ets1,
      FOS = p_fos,
      SMAD3 = p_smad3
    ),
    da_peaks = da_peaks_smad3
  )
}

saveRDS(results_list, "results_TF/FAP1_4_chromVAR_DA_results.rds")

# ===== Master facet figure (FAP1–4 × 4 motif) =====
df_all <- bind_rows(all_plots_for_facet)

df_plot <- df_all %>%
  tidyr::pivot_longer(
    cols = starts_with("MA"),
    names_to = "feature_real",
    values_to = "value"
  ) %>%
  dplyr::filter(feature_real == feature) %>%
  dplyr::select(condition, FAP, TF, feature, value)

df_plot$condition <- factor(df_plot$condition, levels=c("Young","Aged"))
df_plot$TF <- factor(df_plot$TF,  levels=c("ETS2","ETS1","FOS","SMAD3"))
df_plot$FAP <- factor(df_plot$FAP, levels=c("FAP1","FAP2","FAP3","FAP4"))

p_facet <- ggplot(df_plot, aes(x=condition, y=value, fill=condition)) +
  geom_half_violin(side="l", alpha=0.6, trim=FALSE) +
  geom_half_boxplot(side="r", width=0.25, alpha=1, outlier.shape = NA) +
  geom_jitter(width=0.1, alpha=0.25, size=0.6) +
  stat_compare_means(
    method="wilcox.test",
    label="p.signif",
    hide.ns=TRUE,
    size=3,
    comparisons=list(c("Young","Aged"))
  ) +
  facet_grid(TF ~ FAP, scales="free_y") +
  scale_fill_manual(values=c("#1f78b4","#e31a1c")) +
  theme_bw() +
  labs(
    y="chromVAR deviation",
    x="",
    title="chromVAR Motif Activity — Raincloud (FAP1–4 × ETS2/ETS1/FOS/SMAD3)"
  )

ggsave(
  "results_TF/Raincloud/Raincloud_FACET_FAP1_4_all_TFs.png",
  p_facet,
  width=15, height=12, dpi=300, bg="white"
)
