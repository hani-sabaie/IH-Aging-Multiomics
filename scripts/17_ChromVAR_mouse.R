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
library(BSgenome.Mmusculus.UCSC.mm10)
library(dplyr)
library(ggplot2)
library(Seurat)
library(chromVAR)
library(gghalves)
library(ggpubr)
library(tibble)
library(tidyr)

# ===== Helpers =====
figdir <- "results_TF_mouse/Raincloud"
save_gg <- function(p, filename, w=7, h=5, dpi=300){
  if (!dir.exists(figdir)) dir.create(figdir, recursive = TRUE)
  if (inherits(p, "ggplot") || inherits(p, "patchwork")){
    ggsave(file.path(figdir, filename), p, width=w, height=h, units="in", dpi=dpi, bg="white")
  } else {
    png(file.path(figdir, filename), width=w, height=h, units="in", res=dpi)
    print(p); dev.off()
  }
}

# ============================================================================ #

# ============================
# Load mouse object
# ============================
obj <- readRDS("../data/GSE288662/Processed_Seurat_Object/GSE288662_Processed_Seurat_Object.rds")

# Harmonize metadata to reuse logic
obj$skeletal_muscle <- obj$type   

cond1 <- "EP"
cond2 <- "Veh"

# FAP-like fibroblast clusters in mouse
fap_levels <- c("Pgr- Fibroblast", "Pgr+ Fibroblast", "Adipocytes")

# ============================
# Prepare ATAC assay (no Links)
# ============================
DefaultAssay(obj) <- "ATAC"
obj <- FindTopFeatures(obj, min.cutoff = "q0.7")

# ============================
# Add motifs Run chromVAR
# ============================
motif_obj <- tryCatch(Motifs(obj_mouse[["ATAC"]]), error=function(e) NULL)
if (is.null(motif_obj)) {
  db <- JASPAR2024()
  sq <- RSQLite::dbConnect(RSQLite::SQLite(), db(db))
  pfm <- getMatrixSet(sq, opts=list(collection="CORE", tax_group="vertebrates"))
  obj <- AddMotifs(obj, genome=BSgenome.Mmusculus.UCSC.mm10, pfm=pfm)
}

# Run chromVAR
obj <- RunChromVAR(
  object = obj,
  genome = BSgenome.Mmusculus.UCSC.mm10
)

saveRDS(obj, "../outputs/mouse_chromVAR_obj.rds")

# ============================
# Motif IDs for TFs of interest (mouse)
# ============================
obj <- readRDS("../data/GSE288662/Processed_Seurat_Object/GSE288662_Processed_Seurat_Object.rds")
DefaultAssay(obj) <- "chromvar"

motif_obj <- Motifs(obj[["ATAC"]])
motif_names <- motif_obj@motif.names
motif_ids <- names(motif_obj@pwm)

get_motif_ids <- function(tf, motif_names, motif_ids) {
  hits <- grep(tf, motif_names, ignore.case = TRUE)
  motif_ids[hits]
}

# JASPAR motif 
motifs_smad3 <- get_motif_ids("SMAD3", motif_names, motif_ids)
motifs_ets1 <- get_motif_ids("ETS1",  motif_names, motif_ids)
motifs_ets2 <- get_motif_ids("ETS2",  motif_names, motif_ids)
motifs_fos <- get_motif_ids("FOS",   motif_names, motif_ids)

motifs_smad3; motifs_ets1; motifs_ets2; motifs_fos

# ============================
# Helper: featur in chromVAR
# ============================
get_chromvar_feature <- function(seu, motif_id) {
  rn <- rownames(seu[["chromvar"]])
  cand <- c(motif_id, paste0("motif_", motif_id))
  hit <- cand[cand %in% rn]
  if (length(hit) == 0) {
    warning("No chromVAR feature found for motif_id = ", motif_id)
    return(NA_character_)
  }
  hit[1]
}

# ============================
# Raincloud
# ============================
plot_raincloud <- function(df, feat, cf, tf_name, save_dir="results_TF_mouse/Raincloud") {
  
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
      comparisons = list(c(cond2, cond1))   # Veh vs EP
    ) +
    scale_fill_manual(values=c("#1f78b4", "#e31a1c")) +
    theme_bw() +
    labs(
      title = paste0(cf, " — ", tf_name, " (mouse)"),
      y = "chromVAR deviation",
      x = ""
    )
  
  ggsave(
    filename = file.path(save_dir, paste0("Raincloud_mouse_", cf, "_", tf_name, ".png")),
    plot = p,
    width = 5, height = 4, dpi = 300, bg = "white"
  )
  
  return(p)
}

# ============================
# Loop over fibroblast-like clusters
# ============================
results_list <- list()
all_plots_for_facet <- list()

for (cf in fap_levels) {
  
  message("==== Analyzing ", cf, " (mouse) ====")
  
  obj_cf <- subset(
    obj,
    subset = skeletal_muscle == cf & condition %in% c(cond1, cond2)
  )
  if (ncol(obj_cf) == 0) {
    warning("No cells for ", cf, " with ", cond1, "/", cond2, "; skipping.")
    next
  }
  
  DefaultAssay(obj_cf) <- "chromvar"
  Idents(obj_cf) <- obj_cf$condition
  
  feat_ets2 <- get_chromvar_feature(obj_cf, motifs_ets2[1])
  feat_ets1 <- get_chromvar_feature(obj_cf, motifs_ets1[1])
  feat_fos <- get_chromvar_feature(obj_cf, motifs_fos[1])
  feat_smad3 <- get_chromvar_feature(obj_cf, motifs_smad3[1])
  
  vars_use <- c("condition", feat_ets2, feat_ets1, feat_fos, feat_smad3)
  vars_use <- vars_use[!is.na(vars_use)]
  
  df_cf <- FetchData(
    obj_cf,
    vars = vars_use
  )
  df_cf$condition <- factor(df_cf$condition, levels = c(cond2, cond1))  # Veh, EP
  
  # Raincloud 
  p_ets2 <- plot_raincloud(df_cf, feat_ets2,  cf, "ETS2")
  p_ets1 <- plot_raincloud(df_cf, feat_ets1,  cf, "ETS1")
  p_fos <- plot_raincloud(df_cf, feat_fos,   cf, "FOS")
  p_smad3 <- plot_raincloud(df_cf, feat_smad3, cf, "SMAD3")
  
  # Facet:
  if (!is.na(feat_ets2))
    all_plots_for_facet[[paste0(cf,"_ETS2")]]  <- df_cf %>% mutate(TF="ETS2",  feature=feat_ets2,  FAP=cf)
  if (!is.na(feat_ets1))
    all_plots_for_facet[[paste0(cf,"_ETS1")]]  <- df_cf %>% mutate(TF="ETS1",  feature=feat_ets1,  FAP=cf)
  if (!is.na(feat_fos))
    all_plots_for_facet[[paste0(cf,"_FOS")]]   <- df_cf %>% mutate(TF="FOS",   feature=feat_fos,   FAP=cf)
  if (!is.na(feat_smad3))
    all_plots_for_facet[[paste0(cf,"_SMAD3")]] <- df_cf %>% mutate(TF="SMAD3", feature=feat_smad3, FAP=cf)
  
  results_list[[cf]] <- list(
    raincloud = list(
      ETS2  = p_ets2,
      ETS1  = p_ets1,
      FOS   = p_fos,
      SMAD3 = p_smad3
    )
  )
}

saveRDS(results_list, "results_TF_mouse/Fibro_mouse_chromVAR_results.rds")

# ============================
# Master facet figure (mouse: FAP-like × 4 motif)
# ============================
df_all <- bind_rows(all_plots_for_facet)

df_plot <- df_all %>%
  tidyr::pivot_longer(
    cols = starts_with("MA"),  
    names_to = "feature_real",
    values_to = "value"
  ) %>%
  dplyr::filter(feature_real == feature) %>%
  dplyr::select(condition, FAP, TF, feature, value)

df_plot$condition <- factor(df_plot$condition, levels=c(cond2,cond1))  # Veh, EP
df_plot$TF <- factor(df_plot$TF,  levels=c("ETS2","ETS1","FOS","SMAD3"))
df_plot$FAP <- factor(df_plot$FAP, levels=fap_levels)

p_facet <- ggplot(df_plot, aes(x=condition, y=value, fill=condition)) +
  geom_half_violin(side="l", alpha=0.6, trim=FALSE) +
  geom_half_boxplot(side="r", width=0.25, alpha=1, outlier.shape = NA) +
  geom_jitter(width=0.1, alpha=0.25, size=0.6) +
  stat_compare_means(
    method="wilcox.test",
    label="p.signif",
    hide.ns=TRUE,
    size=3,
    comparisons=list(c(cond2,cond1))
  ) +
  facet_grid(TF ~ FAP, scales="free_y") +
  scale_fill_manual(values=c("#1f78b4","#e31a1c")) +
  theme_bw() +
  labs(
    y="chromVAR deviation",
    x="",
    title=paste0("Mouse chromVAR Motif Activity — Raincloud (", cond1, " vs ", cond2, ")")
  )

ggsave(
  "results_TF_mouse/Raincloud/Raincloud_FACET_mouse_fibro_all_TFs.png",
  p_facet,
  width=15, height=12, dpi=300, bg="white"
)

df_summary <- df_plot %>%
  group_by(FAP, TF, condition) %>%
  summarise(
    n = n(),
    mean = mean(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    .groups = "drop"
  )

df_change <- df_summary %>%
  pivot_wider(
    id_cols = c(FAP, TF),
    names_from = condition,
    values_from = mean
  ) %>%
  mutate(
    change = EP - Veh,
    direction = case_when(
      change > 0  ~ "Up (in EP)",
      change < 0  ~ "Down (in EP)",
      TRUE        ~ "No change"
    )
  )

View(df_change)
