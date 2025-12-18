library(dplyr)
library(ggplot2)
library(patchwork)  # for combining two plots

# 1) Read the CSV file for GCTA-COJO results
gcta_cojo_data <- data.frame(
  Chr = c(15, 15, 15, 15),
  SNP = c("rs35874463", "rs12912045", "rs1965269", "rs3784681"),  # Two SNPs for each study
  bp = c(67457698, 67467297, 67218945, 67472185),
  refA = c("G", "T", "G", "C"),
  freq = c(0.0576908, 0.209548, 0.09438, 0.2432),
  b = c(0.0937634, -0.0611749, -0.088, -0.0638),
  se = c(0.0196655, 0.01139, 0.0218, 0.0149),
  p = c(1.86136e-06, 7.83337e-08, 5.4209e-05, 1.85322e-05),
  n = c(375534, 367402, 208053, 206819),
  freq_geno = c(0.0526839, 0.229622, 0.106362, 0.282306),
  bJ = c(0.082597, -0.0556112, -0.0961371, -0.0690894),
  bJ_se = c(0.0198004, 0.0114683, 0.0218718, 0.0149491),
  pJ = c(3.02623e-05, 1.23998e-06, 1.10526e-05, 3.80729e-06),
  LD_r = c(-0.117578, 0, -0.0807382, 0)
)

# Adding study column
gcta_cojo_data$study <- c("UKB", "UKB", "FinnGen", "FinnGen")

# 2) Transform the data for visualization
gcta_cojo_data <- gcta_cojo_data %>%
  mutate(
    log_p = -log10(p),    # Calculate log10(p) for visualization
    rank = rank(-log_p),  # Rank based on log(p)
    highlight = case_when(
      rank == 1 ~ "top1",
      rank == 2 ~ "top2",
      TRUE ~ "other"
    )
  )

# 3) Plot the GCTA-COJO results
p_gcta_cojo <- ggplot(gcta_cojo_data, aes(x = reorder(SNP, log_p), y = log_p, fill = study)) +
  geom_col(width = 0.8) +
  geom_label(aes(label = round(log_p, 2)),
             fill = "white", color = "black", 
             linewidth = 0.25, size = 3, vjust = -0.2) +
  scale_fill_manual(values = c("UKB" = "#6EC7D4", "FinnGen" = "#F5D4A4")) +  # Different colors for studies
  geom_hline(yintercept = -log10(5e-5), linetype = "dashed", color = "#d73027", linewidth = 0.6) +  # p = 5e-5 threshold
  annotate("label", x = 2, y = -log10(5e-5) + 0.1, label = "p = 5e-5", hjust = 0.5, vjust = 0, size = 3, color = "#5b0000") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(title = "GCTA-COJO Results for SMAD3", x = "SNP", y = "-log10(p)") +
  scale_x_discrete(position = "bottom") +
  coord_cartesian(ylim = c(0, max(gcta_cojo_data$log_p) + 0.2), clip = "off") +
  theme_classic() +
  theme(
    plot.title = element_text(face = "italic", hjust = 0),
    axis.text.x.bottom = element_text(angle = 30, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 9),
    axis.ticks.y = element_line(color = "black"),
    plot.margin = margin(t = 5, r = 40, b = 0, l = 5)
  )

# 4) Save the plot
ggsave("SMAD3_GCTA_COJO_results_with_tissue_and_study_threshold_5e-5.png", p_gcta_cojo, width = 6, height = 4, dpi = 300)
