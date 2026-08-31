# ===== Prepare source data for Figure 10 =====

library(data.table)

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

source_dir <- file.path(
  repo_root,
  "source_data",
  "figure_source_data",
  "Figure_10"
)

dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

# =========================================================================
# Figure 10A: UpSet plot
# =========================================================================

modules <- fread(
  file.path(
    repo_root,
    "processed_results",
    "05_hdWGCNA",
    "module_assignment_table.csv"
  )
)

deg <- fread(
  file.path(
    repo_root,
    "processed_results",
    "02_differential_expression",
    "bulk_de_sig_faps.csv"
  )
)

smr_sig <- fread(
  file.path(
    repo_root,
    "processed_results",
    "06_SMR_HEIDI",
    "SMR_all_studies_all_tissues_sig.csv"
  )
)

yellow_genes <- unique(
  modules[tolower(module) == "yellow", gene_name]
)

brown_genes <- unique(
  modules[tolower(module) == "brown", gene_name]
)

blue_genes <- unique(
  modules[tolower(module) == "blue", gene_name]
)

deg_genes <- unique(deg$gene)

smr_ukb_genes <- unique(
  smr_sig[study == "UKB", Gene]
)

smr_finngen_genes <- unique(
  smr_sig[study == "FinnGen", Gene]
)

all_genes <- sort(unique(c(
  yellow_genes,
  brown_genes,
  blue_genes,
  deg_genes,
  smr_ukb_genes,
  smr_finngen_genes
)))

fig10a <- data.table(
  gene = all_genes,
  Yellow_Module = as.integer(all_genes %in% yellow_genes),
  Brown_Module = as.integer(all_genes %in% brown_genes),
  Blue_Module = as.integer(all_genes %in% blue_genes),
  DEGs_FAPs = as.integer(all_genes %in% deg_genes),
  SMR_UKB = as.integer(all_genes %in% smr_ukb_genes),
  SMR_FinnGen = as.integer(all_genes %in% smr_finngen_genes)
)

fwrite(
  fig10a,
  file.path(source_dir, "Figure10A_Upset_source_data.csv")
)

# =========================================================================
# Figure 10B-C: SMR/HEIDI results for SMAD3
# =========================================================================

smr_all <- fread(
  file.path(
    repo_root,
    "processed_results",
    "06_SMR_HEIDI",
    "SMR_all_studies_all_tissues.csv"
  )
)

fig10bc <- smr_all[
  Gene == "SMAD3"
][
  order(study, tissue, p_SMR)
][
  ,
  .SD[1],
  by = .(study, tissue)
]

fig10bc[, `:=`(
  neg_log10_pSMR = -log10(p_SMR),
  neg_log10_pHEIDI = -log10(p_HEIDI)
)]

fwrite(
  fig10bc,
  file.path(source_dir, "Figure10BC_SMR_HEIDI_source_data.csv")
)

# =========================================================================
# Figure 10D: corrected GCTA-COJO independent signals
# =========================================================================

ukb_cojo <- fread(
  file.path(
    repo_root,
    "processed_results",
    "07_GCTA_COJO",
    "UKB_SMAD3_GCTA_corrected.jma.cojo"
  )
)
ukb_cojo[, study := "UKB"]

fin_cojo <- fread(
  file.path(
    repo_root,
    "processed_results",
    "07_GCTA_COJO",
    "Finn_SMAD3_GCTA_corrected.jma.cojo"
  )
)
fin_cojo[, study := "FinnGen"]

fig10d <- rbindlist(
  list(ukb_cojo, fin_cojo),
  use.names = TRUE,
  fill = TRUE
)

fig10d[, neg_log10_pJ := -log10(pJ)]

setcolorder(
  fig10d,
  c(
    "study",
    setdiff(names(fig10d), "study")
  )
)

fwrite(
  fig10d,
  file.path(source_dir, "Figure10D_GCTA_COJO_source_data.csv")
)

# =========================================================================
# Figure 10E: standard COLOC results for SMAD3
# =========================================================================

coloc <- fread(
  file.path(
    repo_root,
    "processed_results",
    "08_colocalization",
    "COLOC_all_studies_sig.csv"
  )
)

fig10e <- coloc[gene == "SMAD3"]

fwrite(
  fig10e,
  file.path(source_dir, "Figure10E_COLOC_source_data.csv")
)

# =========================================================================
# Summary
# =========================================================================

cat("Figure 10 source-data directory:", source_dir, "\n")
cat("Figure 10A rows:", nrow(fig10a), "\n")
cat("Figure 10B-C rows:", nrow(fig10bc), "\n")
cat("Figure 10D rows:", nrow(fig10d), "\n")
cat("Figure 10E rows:", nrow(fig10e), "\n")