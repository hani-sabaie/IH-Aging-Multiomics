# ===== Clean environment =====
rm(list = ls(all.names = TRUE))
gc()

# ===== Libraries =====
library(data.table)
library(dplyr)
library(Seurat)
library(Signac)
library(hdWGCNA)

# ===== Repo root =====
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (length(file_arg) == 1) {
  script_path <- sub("^--file=", "", file_arg)
  script_dir <- dirname(normalizePath(script_path))
  repo_root <- normalizePath(file.path(script_dir, ".."))
} else {
  repo_root <- normalizePath(".")
}

pick_first_existing <- function(paths) {
  paths <- paths[nzchar(paths)]
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) {
    stop("None of the candidate files exist:\n", paste(paths, collapse = "\n"))
  }
  normalizePath(hit[1], winslash = "/", mustWork = TRUE)
}

srcdir <- file.path(repo_root, "source_data", "figure_source_data", "Figure_15")
dir.create(srcdir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Figure 15A: differential regulons
# =============================================================================
dregs_file <- pick_first_existing(c(
  Sys.getenv("FIG15_DREGS_FILE"),
  file.path(
    repo_root,
    "processed_results",
    "10_TF_network",
    "differential_regulons.csv"
  )
))

cat("Using differential regulons file:\n", dregs_file, "\n\n")

dregs <- fread(dregs_file)

if (!"tf" %in% names(dregs)) {
  first_col <- names(dregs)[1]
  setnames(dregs, first_col, "tf")
}

dregs_clean <- dregs %>%
  dplyr::mutate(
    direction = dplyr::case_when(
      avg_log2FC_deg > 0 ~ "Aged_up",
      avg_log2FC_deg < 0 ~ "Young_up",
      TRUE ~ "no_change"
    )
  ) %>%
  dplyr::rename(TF = tf) %>%
  dplyr::select(
    TF,
    avg_log2FC_deg,
    p_val_adj_deg,
    avg_log2FC_positive,
    p_val_adj_positive,
    avg_log2FC_negative,
    p_val_adj_negative,
    module,
    kME,
    direction
  )

fwrite(
  as.data.table(dregs_clean),
  file.path(srcdir, "Figure15A_differential_regulons_source_data.csv")
)

# =============================================================================
# Figure 15D: FOS expression + regulon scores + UMAP
# =============================================================================
tfnet_obj_file <- pick_first_existing(c(
  Sys.getenv("FIG15_TFNET_OBJ"),
  file.path(repo_root, "outputs", "hdWGCNA_TFNet_obj.rds")
))

cat("Using TF-network object:\n", tfnet_obj_file, "\n\n")
obj_tf <- readRDS(tfnet_obj_file)

# =============================================================================
# Figure 15D: ETS2 / ETS1 / FOS expression + regulon scores + UMAP
# =============================================================================

tfs_fig15d <- c("ETS2", "ETS1", "FOS")

pos_regulon_scores <- GetRegulonScores(
  obj_tf,
  target_type = "positive"
)

neg_regulon_scores <- GetRegulonScores(
  obj_tf,
  target_type = "negative"
)

red_candidates <- c(
  "wnn.umap",
  "umap",
  "harmony.sct.umap",
  "harmony.atac.umap"
)

use_red <- red_candidates[
  red_candidates %in% names(obj_tf@reductions)
][1]

if (is.na(use_red) || !nzchar(use_red)) {
  stop("No suitable UMAP reduction found in object.")
}

emb <- Embeddings(
  obj_tf,
  reduction = use_red
)

cells <- rownames(emb)

fig15d_list <- lapply(
  tfs_fig15d,
  function(cur_tf) {

    expr <- FetchData(
      obj_tf,
      vars = cur_tf,
      cells = cells
    )[, 1]

    data.table(
      cell = cells,
      TF = cur_tf,
      UMAP_1 = emb[cells, 1],
      UMAP_2 = emb[cells, 2],
      TF_expression = expr,
      positive_regulon_score =
        pos_regulon_scores[
          match(cells, rownames(pos_regulon_scores)),
          cur_tf
        ],
      negative_regulon_score =
        neg_regulon_scores[
          match(cells, rownames(neg_regulon_scores)),
          cur_tf
        ],
      condition = obj_tf@meta.data[cells, "condition"],
      skeletal_muscle =
        obj_tf@meta.data[cells, "skeletal_muscle"]
    )
  }
)

fig15d <- rbindlist(fig15d_list)

fwrite(
  fig15d,
  file.path(
    srcdir,
    "Figure15D_ETS2_ETS1_FOS_regulon_UMAP_source_data.csv"
  )
)

# =============================================================================
# Figure 15B and 15C: Signac object
# =============================================================================
foot_obj_file <- pick_first_existing(c(
  Sys.getenv("FIG15_FOOT_OBJ"),
  file.path(repo_root, "outputs", "hdWGCNA_TFNet_DEReg_L2G_Foot_obj.rds"),
  file.path(repo_root, "outputs", "hdWGCNA_TFNet_DEReg_L2G_obj.rds")
))

cat("Using Signac footprint object:\n", foot_obj_file, "\n\n")
obj_foot <- readRDS(foot_obj_file)

# ----- Figure 15B: SMAD3 region metadata -----
smad3_coords <- LookupGeneCoords(obj_foot, "SMAD3")
smad3_coords_df <- as.data.frame(smad3_coords)
fwrite(
  as.data.table(smad3_coords_df),
  file.path(srcdir, "Figure15B_SMAD3_gene_coords_source_data.tsv"),
  sep = "\t"
)

open_fap3 <- c("chr15-67109227-67110381", "chr15-67198215-67198971")
regions_highlight <- subsetByOverlaps(
  StringToGRanges(open_fap3),
  LookupGeneCoords(obj_foot, "SMAD3")
)
regions_highlight_df <- as.data.frame(regions_highlight)

fwrite(
  as.data.table(regions_highlight_df),
  file.path(srcdir, "Figure15B_SMAD3_highlight_regions_source_data.tsv"),
  sep = "\t"
)

smad3_links <- Links(obj_foot)
smad3_links <- subsetByOverlaps(smad3_links, LookupGeneCoords(obj_foot, "SMAD3"))
smad3_links_df <- as.data.frame(smad3_links)

fwrite(
  as.data.table(smad3_links_df),
  file.path(srcdir, "Figure15B_SMAD3_peak_gene_links_source_data.tsv"),
  sep = "\t"
)

# ----- Figure 15C: footprint matrices -----
pe_list <- obj_foot[["ATAC"]]@positionEnrichment

if (length(pe_list) == 0) {
  warning("No positionEnrichment entries found in ATAC assay.")
} else {
  pe_names <- names(pe_list)
  cat("positionEnrichment entries:\n")
  print(pe_names)

  keep_entries <- pe_names[grepl("ETS2|ETS1|FOS", pe_names, ignore.case = TRUE)]

  if (length(keep_entries) == 0) {
    keep_entries <- pe_names
  }

  for (nm in keep_entries) {
    cur <- pe_list[[nm]]
    cur_df <- as.data.frame(cur)
    cur_df <- cbind(position = rownames(cur_df), cur_df)
    out_file <- file.path(
      srcdir,
      paste0("Figure15C_", gsub("[^A-Za-z0-9_]+", "_", nm), "_footprint_source_data.tsv")
    )
    fwrite(as.data.table(cur_df), out_file, sep = "\t")
  }
}

# =============================================================================
# Summary
# =============================================================================
cat("\nFigure 15 source data generated.\n")
cat("Figure15A rows:", nrow(dregs_clean), "\n")
cat("Figure15D rows:", nrow(fig15d), "\n")
cat("Figure15B gene-coord rows:", nrow(smad3_coords_df), "\n")
cat("Figure15B highlighted-region rows:", nrow(regions_highlight_df), "\n")
cat("Figure15B peak-gene-link rows:", nrow(smad3_links_df), "\n")
cat("Done.\n")
