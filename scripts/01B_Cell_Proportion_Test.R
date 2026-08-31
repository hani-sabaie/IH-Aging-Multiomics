# ========================
# Setting up environment
# ========================
# ===== Clean environment =====
rm(list = ls(all.names = TRUE))
gc()

# ===== Loading relevant libraries =====
library(Seurat)
library(ggplot2)
library(scProportionTest)

# Locate repository root from this script
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (length(file_arg) == 1) {
  script_path <- sub("^--file=", "", file_arg)
  script_dir <- dirname(normalizePath(script_path))
  repo_root <- normalizePath(file.path(script_dir, ".."))
} else {
  repo_root <- normalizePath(".")
}

# ===== set seed =====
set.seed(1234)

# ===== Helpers =====
outdir <- file.path(repo_root, "outputs")
figdir <- file.path(outdir, "sc", "figs")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

save_gg <- function(p, filename, w=7, h=5, dpi=300){
  if (inherits(p, "ggplot") || inherits(p, "patchwork")){
    ggsave(file.path(figdir, filename), p, width=w, height=h, units="in", dpi=dpi, bg="white")
  } else {
    png(file.path(figdir, filename), width=w, height=h, units="in", res=dpi)
    print(p); dev.off()
  }
}

# ============================================================================ #

# ========================
# Load the object
# ========================
obj_file <- Sys.getenv(
  "CELLPROP_SEURAT_RDS",
  unset = file.path(
    repo_root,
    "outputs",
    "decont_merged_filt_nodoub_cc_sct_reduc_clust_integ_annot_obj.rds"
  )
)

if (!file.exists(obj_file)) {
  stop("Required Seurat object not found: ", obj_file)
}

obj <- readRDS(obj_file)

# ========================
# Single Cell Proportion Test 
# ========================
prop_test <- sc_utils(obj)

prop_test <- permutation_test(
  prop_test, cluster_identity = "skeletal_muscle",
  sample_1 = "Young", sample_2 = "Aged",
  sample_identity = "condition"
)

# Save numerical permutation-test results
prop_results <- as.data.frame(prop_test@results$permutation)

processed_dir <- file.path(
  repo_root,
  "processed_results",
  "03_cell_composition"
)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(
  prop_results,
  file.path(processed_dir, "cell_proportion_test_results.csv"),
  row.names = FALSE
)

p1 <- permutation_plot(prop_test)

save_gg(p1, "Cell_Proportion.png", w= 8, h=5)
