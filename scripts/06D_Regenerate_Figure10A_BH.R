# Reviewer C7
# Regenerate Figure 10A using BH-corrected UKB SMR discovery results.
#
# This is an audit/revision script.
# It does NOT overwrite the current canonical Figure 10A.

rm(list = ls(all.names = TRUE))
gc()

suppressPackageStartupMessages({
  library(ggVennDiagram)
  library(ggplot2)
  library(data.table)
})

# -------------------------------------------------------------------------
# Locate repository root
# -------------------------------------------------------------------------

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (length(file_arg) == 1) {
  script_path <- sub("^--file=", "", file_arg)
  script_dir <- dirname(normalizePath(script_path))
  repo_root <- normalizePath(file.path(script_dir, ".."))
} else {
  repo_root <- normalizePath(".")
}

# -------------------------------------------------------------------------
# Paths
# -------------------------------------------------------------------------

module_file <- file.path(
  repo_root,
  "processed_results",
  "05_hdWGCNA",
  "module_assignment_table.csv"
)

deg_file <- file.path(
  repo_root,
  "processed_results",
  "02_differential_expression",
  "bulk_de_sig_faps.csv"
)

ukb_bh_file <- file.path(
  repo_root,
  "processed_results",
  "06_SMR_HEIDI",
  "multiple_testing_audit",
  "discovery_replication",
  "UKB_all_tissues_with_BH.tsv"
)

outdir <- file.path(repo_root, "outputs", "Venn")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

audit_dir <- file.path(
  repo_root,
  "processed_results",
  "06_SMR_HEIDI",
  "multiple_testing_audit",
  "discovery_replication"
)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# Load data
# -------------------------------------------------------------------------

modules <- fread(module_file)
deg <- fread(deg_file)
ukb <- fread(ukb_bh_file)

# -------------------------------------------------------------------------
# Define gene sets
# -------------------------------------------------------------------------

yellow_genes <- unique(
  na.omit(modules[tolower(module) == "yellow", gene_name])
)

brown_genes <- unique(
  na.omit(modules[tolower(module) == "brown", gene_name])
)

blue_genes <- unique(
  na.omit(modules[tolower(module) == "blue", gene_name])
)

deg_genes <- unique(
  na.omit(deg$gene)
)

# Primary corrected UKB discovery definition:
# BH across all UKB gene-tissue tests + HEIDI criterion
smr_ukb_bh_genes <- unique(
  na.omit(
    ukb[
      !is.na(FDR_cohortwide) &
      FDR_cohortwide < 0.05 &
      !is.na(p_HEIDI) &
      p_HEIDI > 0.01,
      Gene
    ]
  )
)

# -------------------------------------------------------------------------
# Report set sizes
# -------------------------------------------------------------------------

cat("\n============================================================\n")
cat("FIGURE 10A: BH-CORRECTED UKB DISCOVERY SETS\n")
cat("============================================================\n\n")

cat("Yellow module:", length(yellow_genes), "genes\n")
cat("Brown module :", length(brown_genes), "genes\n")
cat("Blue module  :", length(blue_genes), "genes\n")
cat("FAP DEGs     :", length(deg_genes), "unique genes\n")
cat("UKB SMR BH   :", length(smr_ukb_bh_genes), "unique genes\n")

# -------------------------------------------------------------------------
# Key convergence
# -------------------------------------------------------------------------

deg_brown <- intersect(deg_genes, brown_genes)
convergent <- Reduce(
  intersect,
  list(
    deg_genes,
    brown_genes,
    smr_ukb_bh_genes
  )
)

cat("\nFAP DEG ∩ Brown module:", length(deg_brown), "\n")
cat(
  "FAP DEG ∩ Brown module ∩ UKB SMR BH:",
  length(convergent),
  "\n"
)

cat("\nConvergent gene(s):\n")
print(convergent)

# -------------------------------------------------------------------------
# UpSet plot
# -------------------------------------------------------------------------

sets <- list(
  Yellow_Module = yellow_genes,
  Brown_Module = brown_genes,
  Blue_Module = blue_genes,
  DEGs_FAPs = deg_genes,
  SMR_UKB_BH = smr_ukb_bh_genes
)

venn <- Venn(sets)

p <- plot_upset(
  venn,
  nintersects = 26,
  order.intersect.by = "size",
  relative_height = 2,
  relative_width = 0.3
)

ggsave(
  filename = file.path(
    outdir,
    "plot_upset_BH_UKB_discovery.png"
  ),
  plot = p,
  width = 8,
  height = 5,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = file.path(
    outdir,
    "plot_upset_BH_UKB_discovery.pdf"
  ),
  plot = p,
  width = 8,
  height = 5,
  units = "in",
  bg = "white"
)

# -------------------------------------------------------------------------
# Revised Figure 10A source-data table
# -------------------------------------------------------------------------

all_genes <- sort(unique(c(
  yellow_genes,
  brown_genes,
  blue_genes,
  deg_genes,
  smr_ukb_bh_genes
)))

fig10a <- data.table(
  gene = all_genes,
  Yellow_Module = as.integer(all_genes %in% yellow_genes),
  Brown_Module = as.integer(all_genes %in% brown_genes),
  Blue_Module = as.integer(all_genes %in% blue_genes),
  DEGs_FAPs = as.integer(all_genes %in% deg_genes),
  SMR_UKB_BH = as.integer(all_genes %in% smr_ukb_bh_genes)
)

fwrite(
  fig10a,
  file.path(
    audit_dir,
    "Figure10A_Upset_BH_source_data.csv"
  )
)

# Save exact convergent genes
fwrite(
  data.table(gene = convergent),
  file.path(
    audit_dir,
    "Figure10A_convergent_genes_BH.csv"
  )
)

# Save UKB SMR rows contributing to corrected discovery set
fwrite(
  ukb[
    Gene %in% smr_ukb_bh_genes &
    !is.na(FDR_cohortwide) &
    FDR_cohortwide < 0.05 &
    !is.na(p_HEIDI) &
    p_HEIDI > 0.01
  ],
  file.path(
    audit_dir,
    "Figure10A_UKB_SMR_BH_discovery_rows.tsv"
  ),
  sep = "\t"
)

cat("\nFiles written:\n")
cat(
  file.path(outdir, "plot_upset_BH_UKB_discovery.png"),
  "\n"
)
cat(
  file.path(outdir, "plot_upset_BH_UKB_discovery.pdf"),
  "\n"
)
cat(
  file.path(audit_dir, "Figure10A_Upset_BH_source_data.csv"),
  "\n"
)
cat(
  file.path(audit_dir, "Figure10A_convergent_genes_BH.csv"),
  "\n"
)
