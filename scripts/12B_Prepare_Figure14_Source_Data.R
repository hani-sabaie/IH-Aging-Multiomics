rm(list = ls(all.names = TRUE))
gc()

library(hdWGCNA)
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

# -------------------------------------------------------------------------
# Input object
# -------------------------------------------------------------------------

obj_file <- Sys.getenv(
  "TFNET_OBJECT_RDS",
  unset = file.path(
    repo_root,
    "outputs",
    "hdWGCNA_TFNet_obj.rds"
  )
)

if (!file.exists(obj_file)) {
  stop(
    "TF-network Seurat object not found: ",
    obj_file,
    "\nSet TFNET_OBJECT_RDS to its local path."
  )
}

# -------------------------------------------------------------------------
# Output directory
# -------------------------------------------------------------------------

source_dir <- file.path(
  repo_root,
  "source_data",
  "figure_source_data",
  "Figure_14"
)

dir.create(
  source_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# -------------------------------------------------------------------------
# Load object
# -------------------------------------------------------------------------

cat("Loading TF-network object...\n")
obj <- readRDS(obj_file)

cat(
  "Active WGCNA:",
  obj@misc$active_wgcna,
  "\n"
)

# -------------------------------------------------------------------------
# Core module-regulatory tables used by hdWGCNA plotting functions
# -------------------------------------------------------------------------

net_tf <- as.data.table(
  ModuleRegulatoryNetwork(
    obj,
    TFs_only = TRUE
  )
)

net_all <- as.data.table(
  ModuleRegulatoryNetwork(
    obj,
    TFs_only = FALSE
  )
)

# Heatmap feature in Figure 14A
net_tf[, delta := score_pos - score_neg]
net_all[, delta := score_pos - score_neg]

fwrite(
  net_tf,
  file.path(
    source_dir,
    "Figure14A_module_regulatory_scores_TFs_only.csv"
  )
)

fwrite(
  net_all,
  file.path(
    source_dir,
    "Figure14A_module_regulatory_scores_all_targets.csv"
  )
)

# -------------------------------------------------------------------------
# Figure 14B
# Original plotting cutoff = 0.5
# -------------------------------------------------------------------------

fwrite(
  net_tf[score_pos >= 0.5],
  file.path(
    source_dir,
    "Figure14B_positive_network_TFs_only_cutoff0.5.csv"
  )
)

fwrite(
  net_tf[score_neg >= 0.5],
  file.path(
    source_dir,
    "Figure14B_negative_network_TFs_only_cutoff0.5.csv"
  )
)

fwrite(
  net_all[score_pos >= 0.5],
  file.path(
    source_dir,
    "Figure14B_positive_network_all_targets_cutoff0.5.csv"
  )
)

fwrite(
  net_all[score_neg >= 0.5],
  file.path(
    source_dir,
    "Figure14B_negative_network_all_targets_cutoff0.5.csv"
  )
)

# -------------------------------------------------------------------------
# Figure 14C
# Original plotting calls use TFs_only = TRUE, positive score,
# cutoff = 0.1, with brown as source or target.
# These subsets follow the hdWGCNA plotting output directly.
# -------------------------------------------------------------------------

fwrite(
  net_tf[
    score_pos >= 0.1 &
    as.character(source) == "brown"
  ],
  file.path(
    source_dir,
    "Figure14C_brown_source_positive_cutoff0.1.csv"
  )
)

fwrite(
  net_tf[
    score_pos >= 0.1 &
    as.character(target) == "brown"
  ],
  file.path(
    source_dir,
    "Figure14C_brown_target_positive_cutoff0.1.csv"
  )
)

# -------------------------------------------------------------------------
# Module metadata / UMAP coordinates used for network layout
# -------------------------------------------------------------------------

modules <- as.data.table(
  GetModules(obj)
)

module_umap <- as.data.table(
  GetModuleUMAP(obj)
)

fwrite(
  modules,
  file.path(
    source_dir,
    "Figure14_module_assignments_source_data.csv"
  )
)

fwrite(
  module_umap,
  file.path(
    source_dir,
    "Figure14_module_UMAP_source_data.csv"
  )
)

# -------------------------------------------------------------------------
# Summary
# -------------------------------------------------------------------------

cat("\nFigure 14 source data generated.\n")
cat("TF-only module pairs:", nrow(net_tf), "\n")
cat("All-target module pairs:", nrow(net_all), "\n")

cat(
  "Figure14B TF-only positive edges:",
  nrow(net_tf[score_pos >= 0.5]),
  "\n"
)

cat(
  "Figure14B TF-only negative edges:",
  nrow(net_tf[score_neg >= 0.5]),
  "\n"
)

cat(
  "Figure14C brown-source edges:",
  nrow(
    net_tf[
      score_pos >= 0.1 &
      as.character(source) == "brown"
    ]
  ),
  "\n"
)

cat(
  "Figure14C brown-target edges:",
  nrow(
    net_tf[
      score_pos >= 0.1 &
      as.character(target) == "brown"
    ]
  ),
  "\n"
)

cat("Done.\n")
