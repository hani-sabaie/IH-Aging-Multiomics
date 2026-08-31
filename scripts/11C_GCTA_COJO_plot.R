library(dplyr)
library(ggplot2)
library(patchwork)  # for combining two plots

# 1) Read the CSV file for GCTA-COJO results
# Locate repository root from this script
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (length(file_arg) == 1) {
  script_path <- sub("^--file=", "", file_arg)
  script_dir <- dirname(normalizePath(script_path))
  repo_root <- normalizePath(file.path(script_dir, ".."))

figdir <- file.path(repo_root, "outputs", "GCTA_COJO")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
} else {
  repo_root <- normalizePath(".")
}

cojo_dir <- file.path(
  repo_root,
  "processed_results",
  "07_GCTA_COJO"
)

ukb_file <- file.path(
  cojo_dir,
  "UKB_SMAD3_GCTA_corrected.jma.cojo"
)

finngen_file <- file.path(
  cojo_dir,
  "Finn_SMAD3_GCTA_corrected.jma.cojo"
)

if (!file.exists(ukb_file) || !file.exists(finngen_file)) {
  stop("Corrected GCTA-COJO .jma.cojo files not found.")
}

ukb_cojo <- read.table(
  ukb_file,
  header = TRUE,
  stringsAsFactors = FALSE
)
ukb_cojo$study <- "UKB"

finngen_cojo <- read.table(
  finngen_file,
  header = TRUE,
  stringsAsFactors = FALSE
)
finngen_cojo$study <- "FinnGen"

gcta_cojo_data <- bind_rows(
  ukb_cojo,
  finngen_cojo
)

# 2) Transform the data for visualization
gcta_cojo_data <- gcta_cojo_data %>%
  mutate(
    log_p = -log10(pJ),    # Calculate log10(p) for visualization
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
  labs(title = "GCTA-COJO Results for SMAD3", x = "SNP", y = expression(-log[10](P[J]))) +
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
ggsave(
  file.path(figdir, "SMAD3_GCTA_COJO_results_with_tissue_and_study_threshold_5e-5.png"),
  p_gcta_cojo,
  width = 6, height = 4, dpi = 300
)
