library(dplyr)
library(ggplot2)
library(patchwork)  # for combining two plots

# Locate repository root from this script
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (length(file_arg) == 1) {
  script_path <- sub("^--file=", "", file_arg)
  script_dir <- dirname(normalizePath(script_path))
  repo_root <- normalizePath(file.path(script_dir, ".."))

figdir <- file.path(repo_root, "outputs", "colocalization")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
} else {
  repo_root <- normalizePath(".")
}

# 1) Read the CSV file for coloc results
coloc_data <- read.csv(
  file.path(
    repo_root,
    "processed_results",
    "08_colocalization",
    "COLOC_all_studies_PPH4_gt_0p75.csv"
  )
) %>%
  filter(gene == "SMAD3") %>%
  rename(
    PP.H4 = pph4
  )

# 2) Transform the data for visualization
coloc_data <- coloc_data %>%
  mutate(
    rank = rank(-PP.H4),    # Rank the data based on PP.H4 for ordering
    highlight = case_when(
      rank == 1 ~ "top1",
      rank == 2 ~ "top2",
      TRUE ~ "other"
    )
  )

# 3) Plot the coloc results with tissue and study color differentiation
p_coloc <- ggplot(coloc_data, aes(x = reorder(top_snp, PP.H4), y = PP.H4, fill = study)) +
  geom_col(width = 0.8) +
  geom_label(aes(label = round(PP.H4, 2)),
             fill = "white", color = "black", 
             linewidth = 0.25, size = 3, vjust = -0.2) +
  scale_fill_manual(values = c("UKB" = "#6EC7D4", "FinnGen" = "#F5D4A4")) +  # Different colors for studies
  geom_hline(yintercept = 0.75, linetype = "dashed", color = "#d73027", linewidth = 0.6) +  # Changed threshold to 0.75
  annotate("label", x = 2, y = 0.75 + 0.02, label = "PP.H4 = 0.75", hjust = 2, vjust = 0, size = 3, color = "#5b0000") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(title = "Coloc Results for SMAD3", x = "Top SNP", y = "PP.H4", subtitle = "Adipose Subcutaneous Tissue") +
  scale_x_discrete(position = "bottom") +
  coord_cartesian(ylim = c(0, max(coloc_data$PP.H4) + 0.2), clip = "off") +
  theme_classic() +
  theme(
    plot.title = element_text(face = "italic", hjust = 0),
    axis.text.x.bottom = element_text(angle = 30, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 9),
    axis.ticks.y = element_line(color = "black"),
    plot.margin = margin(t = 5, r = 40, b = 0, l = 5)
  )

# 4) Save the plot
ggsave(
  file.path(figdir, "SMAD3_coloc_results_with_tissue_and_study_threshold_0.75.png"),
  p_coloc,
  width = 6, height = 4, dpi = 300
)
