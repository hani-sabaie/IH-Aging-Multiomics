library(dplyr)
library(ggplot2)
library(patchwork)  # for combining two plots

# 1) Read the CSV file
smr_all <- read.csv(
  "C:\\Users\\Hani\\Desktop\\smr-1.3.1-win-x86_64\\SMR_IHR_Local\\all_results\\SMR_all_studies_all_tissues.csv",  # Adjust your path here
  header = TRUE
)

# 2) Compute −log10(p)
smr_all <- smr_all %>%
  mutate(
    log_pSMR   = -log10(p_SMR),
    log_pHEIDI = -log10(p_HEIDI)
  )

gene_of_interest <- "SMAD3"   # Set your gene here

# Colors for highlighted bars (top 2 by pSMR) and others
fill_cols <- c(
  top1  = "#6EC7D4",   # teal
  top2  = "#F5D4A4",   # beige
  other = "grey80"     # grey
)

# Thresholds based on your method:
thr_smr   <- -log10(0.05)   # pSMR = 0.05
thr_heidi <- -log10(0.01)   # pHEIDI = 0.01

# ------------- LOOP OVER STUDIES ------------- #

for (this_study in unique(smr_all$study)) {
  
  # 3) Subset by gene AND study, keep best SMR per tissue
  dat_gene <- smr_all %>%
    filter(Gene == gene_of_interest,
           study == this_study) %>%
    group_by(tissue) %>%
    slice_min(p_SMR, with_ties = FALSE) %>%  # Keep best SMR per tissue
    ungroup() %>%
    arrange(desc(log_pSMR)) %>%
    mutate(
      tissue_snp = paste(tissue, topSNP, sep = "\n"),
      tissue_snp = factor(tissue_snp, levels = tissue_snp),
      rank = dplyr::row_number(),
      highlight = dplyr::case_when(
        rank == 1 ~ "top1",
        rank == 2 ~ "top2",
        TRUE ~ "other"
      )
    )
  
  if (nrow(dat_gene) == 0) next
  
  # Calculate y-limits based on the data for SMR and HEIDI plots
  ymax_smr   <- max(dat_gene$log_pSMR, na.rm = TRUE) + 0.8
  ymax_heidi <- max(dat_gene$log_pHEIDI, na.rm = TRUE) + 1
  ymin_heidi <- min(0, min(dat_gene$log_pHEIDI, na.rm = TRUE)) + 0.12  # Set ymin_heidi slightly above zero to fix the issue
  
  # Small data frames for nice threshold labels
  thr_lab_smr <- data.frame(
    tissue_snp = tail(levels(dat_gene$tissue_snp), 1),
    y = thr_smr
  )
  
  thr_lab_heidi <- data.frame(
    tissue_snp = tail(levels(dat_gene$tissue_snp), 1),
    y = thr_heidi
  )
  
  # 4) Upper plot: −log10(p_SMR)
  p1 <- ggplot(dat_gene,
               aes(x = tissue_snp, y = log_pSMR, fill = highlight)) +
    geom_col(width = 0.8) +
    geom_label(aes(label = round(log_pSMR, 2)),
               fill = "white",
               color = "black",
               linewidth = 0.25,
               size = 3,
               vjust = -0.2) +
    geom_hline(yintercept = thr_smr, linetype = "dashed",
               color = "#d73027", linewidth = 0.6) +
    annotate("label",
             x = length(levels(dat_gene$tissue_snp)) - 0.1,
             y = thr_smr + 0.4,
             label = "pSMR = 0.05",
             hjust = 1, vjust = 1.1,
             size = 3,
             linewidth = 0.25,
             fill = "white", color = "#5b0000") +
    scale_fill_manual(values = fill_cols, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(
      title = paste0("SMR Results for", " (", this_study, ")"),
      x = NULL, y = NULL
    ) +
    scale_x_discrete(position = "bottom") +
    coord_cartesian(ylim = c(0, ymax_smr), clip = "off") +
    theme_classic() +
    theme(
      plot.title = element_text(face = "italic", hjust = 0),
      axis.text.x.bottom = element_text(angle = 30, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 9),
      axis.ticks.y = element_line(color = "black"),
      plot.margin = margin(t = 5, r = 40, b = 0, l = 5)
    )
  
  # 5) Lower plot: −log10(p_HEIDI) (Flip HEIDI axis)
  p2 <- ggplot(dat_gene,
               aes(x = tissue_snp, y = log_pHEIDI, fill = highlight)) +
    geom_col(width = 0.8) +
    geom_label(aes(label = round(log_pHEIDI, 2)),
               fill = "white",
               color = "black",
               linewidth = 0.25,
               size = 3,
               vjust = 1.2) +
    geom_hline(yintercept = thr_heidi, linetype = "dashed",
               color = "#4575b4", linewidth = 0.6) +
    geom_label(data = thr_lab_heidi,
               aes(x = tissue_snp, y = y, label = "pHEIDI = 0.01"),
               inherit.aes = FALSE,
               hjust = 1.1, vjust = 1.3,
               size = 3, linewidth = 0.25,
               fill = "white", color = "#2b3f73") +
    scale_fill_manual(values = fill_cols, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = NULL, y = NULL) +
    scale_x_discrete(position = "top") +  # Keep x-axis but remove labels
    coord_cartesian(ylim = c(ymin_heidi, ymax_heidi), clip = "off") +
    scale_y_reverse() +  # Flip HEIDI axis here
    theme_classic() +
    theme(
      axis.text.x.top = element_blank(),  # Remove x-axis labels from p2
      axis.text.y = element_text(size = 9),
      axis.ticks.y = element_line(color = "black"),
      plot.margin = margin(t = 0, r = 40, b = 5, l = 5)
    )
  
  # 6) Combine HEIDI (bottom) and SMR (top)
  final_plot <- p1 / p2 + plot_layout(heights = c(2, 1))  # Adjust the height ratios
  
  # 7) Save per study
  out_name <- paste0(gene_of_interest, "_", this_study, "_SMR_HEIDI.png")
  ggsave(out_name, final_plot, width = 6, height = 5, dpi = 300)
}

