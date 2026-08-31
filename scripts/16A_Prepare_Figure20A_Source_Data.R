rm(list = ls(all.names = TRUE))
gc()

library(Seurat)
library(data.table)

# -------------------------------------------------------------------------
# Repository root
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
# Input / output
# -------------------------------------------------------------------------
obj_file <- Sys.getenv("FIG20_MOUSE_RDS")

if (obj_file == "" || !file.exists(obj_file)) {
  stop(
    "Mouse Seurat object not found. Set FIG20_MOUSE_RDS to ",
    "GSE288662_Processed_Seurat_Object.rds"
  )
}

srcdir <- file.path(
  repo_root,
  "source_data",
  "figure_source_data",
  "Figure_20"
)

dir.create(srcdir, recursive = TRUE, showWarnings = FALSE)

cat("Loading mouse Seurat object...\n")
cat("Input:", obj_file, "\n")

obj <- readRDS(obj_file)

# -------------------------------------------------------------------------
# Validate objects used in original Figure 20A
# -------------------------------------------------------------------------
required_meta <- c(
  "wsnn_res.0.8",
  "type"
)

missing_meta <- setdiff(
  required_meta,
  colnames(obj@meta.data)
)

if (length(missing_meta) > 0) {
  stop(
    "Missing metadata columns: ",
    paste(missing_meta, collapse = ", ")
  )
}

if (!"wnn.umap" %in% names(obj@reductions)) {
  stop("wnn.umap reduction not found.")
}

emb <- Embeddings(
  obj,
  reduction = "wnn.umap"
)

if (ncol(emb) < 2) {
  stop("wnn.umap does not contain at least two dimensions.")
}

if (!identical(rownames(emb), rownames(obj@meta.data))) {
  meta <- obj@meta.data[rownames(emb), , drop = FALSE]
} else {
  meta <- obj@meta.data
}

# -------------------------------------------------------------------------
# Figure 20A source data
#
# Original plot:
#   left  = wsnn_res.0.8
#   right = type
# -------------------------------------------------------------------------
fig20a <- data.table(
  cell = rownames(emb),
  WNN_UMAP_1 = as.numeric(emb[, 1]),
  WNN_UMAP_2 = as.numeric(emb[, 2]),
  wsnn_res_0_8 = as.character(meta[["wsnn_res.0.8"]]),
  cell_type = as.character(meta[["type"]])
)

# Include useful provenance metadata if present.
optional_meta <- c(
  "condition",
  "orig.ident",
  "sample",
  "Sample"
)

for (nm in optional_meta) {
  if (nm %in% colnames(meta)) {
    fig20a[[nm]] <- as.character(meta[[nm]])
  }
}

fwrite(
  fig20a,
  file.path(
    srcdir,
    "Figure20A_mouse_WNN_UMAP_source_data.csv"
  )
)

# -------------------------------------------------------------------------
# Annotation counts
# -------------------------------------------------------------------------
cluster_counts <- fig20a[
  ,
  .N,
  by = wsnn_res_0_8
][order(wsnn_res_0_8)]

type_counts <- fig20a[
  ,
  .N,
  by = cell_type
][order(cell_type)]

fwrite(
  cluster_counts,
  file.path(
    srcdir,
    "Figure20A_mouse_cluster_counts.csv"
  )
)

fwrite(
  type_counts,
  file.path(
    srcdir,
    "Figure20A_mouse_cell_type_counts.csv"
  )
)

# -------------------------------------------------------------------------
# Summary
# -------------------------------------------------------------------------
cat("\n===== Figure 20A source-data summary =====\n")
cat("Cells:", nrow(fig20a), "\n")
cat("WNN dimensions:", ncol(emb), "\n")
cat("Number of wsnn_res.0.8 clusters:", uniqueN(fig20a$wsnn_res_0_8), "\n")
cat("Number of annotated cell types:", uniqueN(fig20a$cell_type), "\n")

cat("\nCluster counts:\n")
print(cluster_counts)

cat("\nCell-type counts:\n")
print(type_counts)

cat("\nDone.\n")

rm(obj, emb, meta)
gc()
