# ============================================================================
# GSE288662 mouse pseudobulk DE — original-design full-precision rebuild
#
# Reproduces the original analysis design:
#   - biological replicate = individual mouse (orig.ident)
#   - conditions = EP, EPR, Veh
#   - all cell types with >=2 pseudobulk replicates in every represented
#     condition are retained, matching the historical production logic
#   - no new per-sample minimum-nuclei threshold
#   - raw RNA counts -> pseudobulk -> edgeR filterByExpr -> TMM ->
#     voom -> limma -> eBayes
#   - all three pairwise contrasts:
#       EPR_vs_EP
#       Veh_vs_EP
#       Veh_vs_EPR
#
# Multiple testing:
#   BH separately across all tested genes within each
#   cell-type x contrast analysis.
#
# Full precision is preserved.
# Primary DEG definition:
#   BH FDR < 0.05
#
# Sensitivity:
#   BH FDR < 0.05 AND |log2FC| >= 1
#
# Audit/candidate only: canonical files are NOT overwritten.
# ============================================================================

rm(list = ls(all.names = TRUE))
gc()

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(edgeR)
  library(limma)
  library(data.table)
})

set.seed(1234)

repo_root <- normalizePath(".")

obj_file <- Sys.getenv(
  "MOUSE_SEURAT_RDS",
  unset = paste0(
    "C:/Users/Hani/Desktop/Hernia/data/GSE288662/",
    "Processed_Seurat_Object/",
    "GSE288662_Processed_Seurat_Object.rds"
  )
)

out_dir <- file.path(
  repo_root,
  "processed_results",
  "13_mouse_validation",
  "original_design_fullprecision"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("===== MOUSE ORIGINAL-DESIGN FULL-PRECISION REBUILD =====\n")
cat("Object:", obj_file, "\n")
cat("Exists:", file.exists(obj_file), "\n\n")

if (!file.exists(obj_file)) {
  stop("Mouse Seurat object not found.")
}

obj <- readRDS(obj_file)

md <- obj@meta.data

required <- c(
  "orig.ident",
  "condition",
  "type"
)

if (!all(required %in% names(md))) {
  stop("Required metadata columns missing.")
}

md <- md[, required, drop = FALSE]

names(md) <- c(
  "sample",
  "condition",
  "celltype"
)

md$sample <- as.character(md$sample)
md$condition <- as.character(md$condition)
md$celltype <- as.character(md$celltype)

# ============================================================================
# Validate experimental design
# ============================================================================

sample_design <- unique(
  md[, c(
    "sample",
    "condition"
  )]
)

if (any(table(sample_design$sample) != 1L)) {
  stop("At least one mouse maps to multiple conditions.")
}

sample_design <- sample_design[
  order(
    sample_design$condition,
    sample_design$sample
  ),
]

cat("===== SAMPLE DESIGN =====\n")
print(sample_design, row.names = FALSE)

cat("\nMice per condition:\n")
print(table(sample_design$condition))

if (
  sum(sample_design$condition == "EP") != 4L ||
  sum(sample_design$condition == "EPR") != 3L ||
  sum(sample_design$condition == "Veh") != 3L
) {
  stop("Expected EP=4, EPR=3, Veh=3.")
}

# ============================================================================
# Observed sample x cell-type groups
# ============================================================================

group_counts <- as.data.frame(
  table(
    sample = md$sample,
    celltype = md$celltype
  ),
  stringsAsFactors = FALSE
)

group_counts <- group_counts[
  group_counts$Freq > 0,
  ,
  drop = FALSE
]

names(group_counts)[3] <- "n_cells"

pb_meta <- merge(
  group_counts,
  sample_design,
  by = "sample",
  all.x = TRUE,
  sort = FALSE
)

pb_meta <- pb_meta[
  order(
    pb_meta$celltype,
    pb_meta$condition,
    pb_meta$sample
  ),
]

pb_meta$pb_key <- paste(
  pb_meta$sample,
  pb_meta$celltype,
  sep = "|||"
)

pb_meta$pb_group <- sprintf(
  "PB%03d",
  seq_len(nrow(pb_meta))
)

if (anyDuplicated(pb_meta$pb_key)) {
  stop("Duplicate pseudobulk key.")
}

fwrite(
  pb_meta,
  file.path(
    out_dir,
    "pseudobulk_metadata.csv"
  )
)

cat(
  "\nObserved pseudobulk groups:",
  nrow(pb_meta),
  "\n"
)

# ============================================================================
# Replicate structure per cell type
# ============================================================================

rep_table <- aggregate(
  sample ~ celltype + condition,
  data = pb_meta,
  FUN = function(x) length(unique(x))
)

names(rep_table)[3] <- "n_samples"

cat("\n===== REPLICATES PER CELL TYPE x CONDITION =====\n")

print(
  rep_table[
    order(
      rep_table$celltype,
      rep_table$condition
    ),
  ],
  row.names = FALSE
)

celltypes <- sort(
  unique(pb_meta$celltype)
)

eligible <- vapply(
  celltypes,
  function(ct) {

    x <- rep_table[
      rep_table$celltype == ct,
      ,
      drop = FALSE
    ]

    # Historical production requirement:
    # every condition represented for the cell type must have >=2 samples.
    length(unique(x$condition)) >= 2L &&
      all(x$n_samples >= 2L)
  },
  logical(1)
)

eligible_celltypes <- celltypes[
  eligible
]

excluded_celltypes <- celltypes[
  !eligible
]

cat(
  "\nEligible cell types:",
  length(eligible_celltypes),
  "\n"
)

cat(
  paste(
    eligible_celltypes,
    collapse = "\n"
  ),
  "\n"
)

cat(
  "\nExcluded cell types:",
  length(excluded_celltypes),
  "\n"
)

if (length(excluded_celltypes) > 0L) {
  cat(
    paste(
      excluded_celltypes,
      collapse = "\n"
    ),
    "\n"
  )
}

if (length(eligible_celltypes) != 16L) {
  stop(
    "Expected 16 historical cell types; found ",
    length(eligible_celltypes)
  )
}

# ============================================================================
# Map cells to pseudobulk groups
# ============================================================================

cell_map <- md

cell_map$pb_key <- paste(
  cell_map$sample,
  cell_map$celltype,
  sep = "|||"
)

cell_map$pb_group <- pb_meta$pb_group[
  match(
    cell_map$pb_key,
    pb_meta$pb_key
  )
]

if (anyNA(cell_map$pb_group)) {
  stop("Some cells could not be mapped to pseudobulk groups.")
}

# ============================================================================
# Layer-wise raw-count pseudobulk aggregation
# ============================================================================

rna <- obj[["RNA"]]

count_layers <- Layers(
  rna,
  search = "^counts"
)

cat("\n===== RNA COUNT LAYERS =====\n")
print(count_layers)

if (length(count_layers) < 1L) {
  stop("No RNA count layer found.")
}

first <- LayerData(
  rna,
  layer = count_layers[1]
)

genes <- rownames(first)

if (is.null(genes)) {
  stop("RNA count layer has no row names.")
}

pb_counts <- matrix(
  0,
  nrow = length(genes),
  ncol = nrow(pb_meta),
  dimnames = list(
    genes,
    pb_meta$pb_group
  )
)

rm(first)
gc()

cells_seen <- 0L
cells_used <- 0L

for (lay in count_layers) {

  cat(
    "Aggregating layer:",
    lay,
    "\n"
  )

  mat <- LayerData(
    rna,
    layer = lay
  )

  if (!identical(
    rownames(mat),
    genes
  )) {
    stop(
      "Gene order differs in layer ",
      lay
    )
  }

  layer_cells <- colnames(mat)

  cells_seen <- cells_seen +
    length(layer_cells)

  idx <- match(
    layer_cells,
    rownames(cell_map)
  )

  use <- !is.na(idx)

  if (!any(use)) {
    rm(mat)
    gc()
    next
  }

  mat_use <- mat[
    ,
    use,
    drop = FALSE
  ]

  groups <- cell_map$pb_group[
    idx[use]
  ]

  cells_used <- cells_used +
    ncol(mat_use)

  for (g in unique(groups)) {

    ii <- which(
      groups == g
    )

    pb_counts[, g] <-
      pb_counts[, g] +
      Matrix::rowSums(
        mat_use[
          ,
          ii,
          drop = FALSE
        ]
      )
  }

  rm(
    mat,
    mat_use
  )

  gc()
}

cat("\nCells seen:", cells_seen, "\n")
cat("Cells used:", cells_used, "\n")
cat("Expected cells:", nrow(md), "\n")

if (
  cells_seen != nrow(md) ||
  cells_used != nrow(md)
) {
  stop("Cell-count validation failed.")
}

if (any(colSums(pb_counts) <= 0)) {
  stop("At least one pseudobulk has zero library size.")
}

rm(
  obj,
  rna,
  md,
  cell_map
)

gc()

cat(
  "PASS: raw-count pseudobulk aggregation completed.\n"
)

# ============================================================================
# Original three-group voom-limma analysis
# ============================================================================

all_results <- list()
summary_results <- list()

for (ct in eligible_celltypes) {

  cat(
    "\n===== ",
    ct,
    " =====\n",
    sep = ""
  )

  meta_ct <- pb_meta[
    pb_meta$celltype == ct,
    ,
    drop = FALSE
  ]

  counts_ct <- pb_counts[
    ,
    meta_ct$pb_group,
    drop = FALSE
  ]

  # Match historical factor ordering:
  # alphabetical => EP, EPR, Veh.
  condition <- factor(
    meta_ct$condition
  )

  tab_cond <- table(condition)

  cat("Replicates:\n")
  print(tab_cond)

  if (
    length(tab_cond) < 2L ||
    any(tab_cond < 2L)
  ) {
    stop(
      "Historical replicate rule failed for ",
      ct
    )
  }

  dge <- DGEList(
    counts = counts_ct,
    group = condition
  )

  keep <- filterByExpr(
    dge,
    group = condition,
    min.count = 10
  )

  dge <- dge[
    keep,
    ,
    keep.lib.sizes = FALSE
  ]

  dge <- calcNormFactors(
    dge
  )

  design <- model.matrix(
    ~ 0 + condition
  )

  colnames(design) <- levels(
    condition
  )

  v <- voom(
    dge,
    design,
    plot = FALSE
  )

  fit <- lmFit(
    v,
    design
  )

  levs <- colnames(design)

  combs <- combn(
    levs,
    2
  )

  contrast_expr <- apply(
    combs,
    2,
    function(z) {
      paste0(
        z[2],
        " - ",
        z[1]
      )
    }
  )

  contrast_names <- apply(
    combs,
    2,
    function(z) {
      paste0(
        z[2],
        "_vs_",
        z[1]
      )
    }
  )

  cargs <- as.list(
    contrast_expr
  )

  names(cargs) <- contrast_names
  cargs$levels <- design

  cm <- do.call(
    makeContrasts,
    cargs
  )

  fit2 <- contrasts.fit(
    fit,
    cm
  )

  fit2 <- eBayes(
    fit2
  )

  for (cn in colnames(cm)) {

    tt <- topTable(
      fit2,
      coef = cn,
      number = Inf,
      adjust.method = "BH",
      sort.by = "none"
    )

    # Independent BH validation.
    bh <- p.adjust(
      tt$P.Value,
      method = "BH"
    )

    max_delta <- max(
      abs(
        bh -
          tt$adj.P.Val
      ),
      na.rm = TRUE
    )

    if (
      !is.finite(max_delta) ||
      max_delta > 1e-12
    ) {
      stop(
        "BH validation failed: ",
        ct,
        " / ",
        cn
      )
    }

    res <- data.frame(
      gene = rownames(tt),
      celltype = ct,
      contrast = cn,
      n_tests_within_family = nrow(tt),
      avg_log2FC = tt$logFC,
      AveExpr = tt$AveExpr,
      t = tt$t,
      p_val = tt$P.Value,
      p_val_adj = tt$adj.P.Val,
      B = tt$B,
      stringsAsFactors = FALSE
    )

    res$primary_DEG <-
      res$p_val_adj < 0.05

    res$absLog2FCge1_sensitivity <-
      res$p_val_adj < 0.05 &
      abs(res$avg_log2FC) >= 1

    key <- paste(
      ct,
      cn,
      sep = "|||"
    )

    all_results[[key]] <- res

    summary_results[[key]] <- data.frame(
      celltype = ct,
      contrast = cn,
      n_genes_tested = nrow(tt),
      n_BH_FDR_lt_005 = sum(
        tt$adj.P.Val < 0.05
      ),
      n_BH_FDR_lt_005_absLog2FCge1 = sum(
        tt$adj.P.Val < 0.05 &
        abs(tt$logFC) >= 1
      ),
      max_BH_validation_delta = max_delta,
      stringsAsFactors = FALSE
    )
  }
}

bulk_de <- rbindlist(
  all_results,
  use.names = TRUE,
  fill = TRUE
)

summary_dt <- rbindlist(
  summary_results,
  use.names = TRUE,
  fill = TRUE
)

setorder(
  bulk_de,
  p_val_adj,
  p_val
)

primary_sig <- bulk_de[
  p_val_adj < 0.05
]

fc_sensitivity <- bulk_de[
  p_val_adj < 0.05 &
    abs(avg_log2FC) >= 1
]

fwrite(
  bulk_de,
  file.path(
    out_dir,
    "bulk_de_all_contrasts_mouse_fullprecision.csv"
  )
)

fwrite(
  primary_sig,
  file.path(
    out_dir,
    "bulk_de_sig_BH_FDR005.csv"
  )
)

fwrite(
  fc_sensitivity,
  file.path(
    out_dir,
    "bulk_de_sig_BH_FDR005_absLog2FCge1_sensitivity.csv"
  )
)

fwrite(
  summary_dt,
  file.path(
    out_dir,
    "analysis_summary.csv"
  )
)

# ============================================================================
# Smad3
# ============================================================================

smad3 <- bulk_de[
  tolower(gene) == "smad3"
]

fwrite(
  smad3,
  file.path(
    out_dir,
    "Smad3_all_contrasts_fullprecision.csv"
  )
)

cat("\n===== ANALYSIS SUMMARY =====\n")

print(
  summary_dt[
    order(
      celltype,
      contrast
    )
  ],
  row.names = FALSE
)

cat(
  "\nTotal tests:",
  nrow(bulk_de),
  "\n"
)

cat(
  "Primary FDR-only significant rows:",
  nrow(primary_sig),
  "\n"
)

cat(
  "FDR + |log2FC|>=1 sensitivity rows:",
  nrow(fc_sensitivity),
  "\n"
)

cat("\n===== SMAD3 RESULTS =====\n")

print(
  smad3[, .(
    gene,
    celltype,
    contrast,
    n_tests_within_family,
    avg_log2FC,
    p_val,
    p_val_adj,
    primary_DEG
  )],
  row.names = FALSE,
  digits = 10
)

cat(
  "\nSmad3 BH-significant rows:",
  sum(smad3$p_val_adj < 0.05),
  "/",
  nrow(smad3),
  "\n"
)

cat(
  "\nPASS: original-design full-precision rebuild completed.\n"
)

cat(
  "Canonical mouse DE files were NOT overwritten.\n"
)
