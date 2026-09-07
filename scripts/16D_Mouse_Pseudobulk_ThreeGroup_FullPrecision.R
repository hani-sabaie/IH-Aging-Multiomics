# ============================================================================
# GSE288662 mouse pseudobulk DE — full three-group rebuild
#
# Biological replicate:
#   individual mouse (orig.ident)
#
# Primary biological contrast:
#   EP vs Veh
#
# Model:
#   common three-group model (Veh + EP + EPR), matching the original analysis.
#
# Eligibility for the primary EP-vs-Veh comparison:
#   >= 3 mice per EP/Veh condition with >= 5 nuclei for the cell type.
#
# Samples with <5 nuclei for a given cell type are not used for that
# cell-type-specific pseudobulk model.
#
# DE:
#   raw RNA count layers -> pseudobulk -> edgeR filterByExpr ->
#   TMM -> voom -> limma -> eBayes
#
# Primary multiple testing:
#   BH across all genes tested within each cell-type-specific EP-vs-Veh
#   contrast.
#
# Sensitivity:
#   BH across all eligible gene x cell-type EP-vs-Veh tests.
#
# Full precision is preserved throughout.
# No canonical files are overwritten.
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

audit_dir <- file.path(
  repo_root,
  "processed_results",
  "13_mouse_validation",
  "three_group_fullprecision"
)

dir.create(
  audit_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("===== MOUSE THREE-GROUP PSEUDOBULK REBUILD =====\n")
cat("Object:", obj_file, "\n")
cat("Exists:", file.exists(obj_file), "\n\n")

if (!file.exists(obj_file)) {
  stop("Mouse object not found.")
}

obj <- readRDS(obj_file)

md <- obj@meta.data

required <- c(
  "orig.ident",
  "condition",
  "type"
)

if (!all(required %in% names(md))) {
  stop("Required mouse metadata columns are missing.")
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

if (!all(
  unique(md$condition) %in%
    c("Veh", "EP", "EPR")
)) {
  stop("Unexpected condition detected.")
}

# ============================================================================
# Validate mouse-level experimental design
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
print(
  sample_design,
  row.names = FALSE
)

cat("\nMice per condition:\n")
print(
  table(sample_design$condition)
)

if (
  sum(sample_design$condition == "Veh") != 3L ||
  sum(sample_design$condition == "EP") != 4L ||
  sum(sample_design$condition == "EPR") != 3L
) {
  stop("Expected Veh=3, EP=4, EPR=3.")
}

# ============================================================================
# Mouse x cell-type nuclei counts
# ============================================================================

all_samples <- sort(
  unique(md$sample)
)

all_celltypes <- sort(
  unique(md$celltype)
)

grid <- expand.grid(
  sample = all_samples,
  celltype = all_celltypes,
  stringsAsFactors = FALSE
)

grid <- merge(
  grid,
  sample_design,
  by = "sample",
  all.x = TRUE,
  sort = FALSE
)

obs <- as.data.frame(
  table(
    sample = md$sample,
    celltype = md$celltype
  ),
  stringsAsFactors = FALSE
)

names(obs)[3] <- "n_cells"

qc <- merge(
  grid,
  obs,
  by = c(
    "sample",
    "celltype"
  ),
  all.x = TRUE,
  sort = FALSE
)

qc$n_cells[
  is.na(qc$n_cells)
] <- 0L

qc$adequate_nuclei <-
  qc$n_cells >= 5L

qc <- qc[
  order(
    qc$celltype,
    qc$condition,
    qc$sample
  ),
]

fwrite(
  qc,
  file.path(
    audit_dir,
    "sample_celltype_nuclei_counts.csv"
  )
)

# ============================================================================
# Eligibility for primary EP-vs-Veh comparison
# ============================================================================

epveh_qc <- qc[
  qc$condition %in%
    c("EP", "Veh"),
]

adequacy <- aggregate(
  adequate_nuclei ~ celltype + condition,
  data = epveh_qc,
  FUN = sum
)

names(adequacy)[3] <-
  "n_adequate_samples"

wide <- reshape(
  adequacy,
  idvar = "celltype",
  timevar = "condition",
  direction = "wide"
)

for (nm in c(
  "n_adequate_samples.EP",
  "n_adequate_samples.Veh"
)) {
  if (!nm %in% names(wide)) {
    wide[[nm]] <- 0L
  }
}

wide$eligible <-
  wide$n_adequate_samples.EP >= 3L &
  wide$n_adequate_samples.Veh >= 3L

wide <- wide[
  order(wide$celltype),
]

fwrite(
  wide,
  file.path(
    audit_dir,
    "celltype_eligibility_EP_vs_Veh.csv"
  )
)

eligible_celltypes <-
  wide$celltype[
    wide$eligible
  ]

cat("\n===== CELL-TYPE ELIGIBILITY =====\n")
print(
  wide,
  row.names = FALSE
)

cat(
  "\nEligible cell types:",
  length(eligible_celltypes),
  "\n"
)

if (length(eligible_celltypes) != 12L) {
  stop(
    "Expected 12 eligible EP-vs-Veh cell types; found ",
    length(eligible_celltypes)
  )
}

# ============================================================================
# Define valid pseudobulk groups
# ============================================================================

pb_meta <- qc[
  qc$celltype %in% eligible_celltypes &
    qc$adequate_nuclei,
  c(
    "sample",
    "condition",
    "celltype",
    "n_cells"
  )
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
    audit_dir,
    "pseudobulk_metadata.csv"
  )
)

cat(
  "\nPseudobulk groups retained:",
  nrow(pb_meta),
  "\n"
)

# Map every retained single cell to one safe pseudobulk group.
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

cell_map <- cell_map[
  !is.na(cell_map$pb_group),
  ,
  drop = FALSE
]

cat(
  "Single cells contributing to valid pseudobulks:",
  nrow(cell_map),
  "\n"
)

# ============================================================================
# Memory-efficient raw-count aggregation layer by layer
# ============================================================================

rna <- obj[["RNA"]]

count_layers <- Layers(
  rna,
  search = "^counts"
)

cat("\n===== RNA COUNT LAYERS =====\n")
print(count_layers)

if (length(count_layers) < 1L) {
  stop("No RNA count layers found.")
}

# Establish common gene universe from first count layer.
m1 <- LayerData(
  rna,
  layer = count_layers[1]
)

genes <- rownames(m1)

if (is.null(genes)) {
  stop("RNA count layer has no gene names.")
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

rm(m1)
gc()

total_cells_seen <- 0L
total_cells_used <- 0L

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

  total_cells_seen <-
    total_cells_seen +
    length(layer_cells)

  map_idx <- match(
    layer_cells,
    rownames(cell_map)
  )

  use <- !is.na(map_idx)

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
    map_idx[use]
  ]

  total_cells_used <-
    total_cells_used +
    ncol(mat_use)

  for (g in unique(groups)) {

    cols <- which(
      groups == g
    )

    pb_counts[, g] <-
      pb_counts[, g] +
      Matrix::rowSums(
        mat_use[
          ,
          cols,
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

cat(
  "\nCells seen across RNA count layers:",
  total_cells_seen,
  "\n"
)

cat(
  "Cells used for eligible pseudobulks:",
  total_cells_used,
  "\n"
)

expected_cells <- sum(
  pb_meta$n_cells
)

cat(
  "Expected contributing cells from metadata:",
  expected_cells,
  "\n"
)

if (total_cells_used != expected_cells) {
  stop(
    "RNA-layer cell count does not match metadata expectation."
  )
}

# Library sizes must match sum of raw counts.
lib_sizes <- colSums(
  pb_counts
)

if (any(lib_sizes <= 0)) {
  stop("At least one pseudobulk has zero library size.")
}

fwrite(
  data.frame(
    pb_group = names(lib_sizes),
    library_size = as.numeric(lib_sizes)
  ),
  file.path(
    audit_dir,
    "pseudobulk_library_sizes.csv"
  )
)

# We no longer need the single-cell object.
rm(
  obj,
  rna,
  md,
  cell_map
)

gc()

cat(
  "PASS: layer-wise raw-count pseudobulk aggregation completed.\n"
)

# ============================================================================
# Three-group voom-limma per eligible cell type
# ============================================================================

result_list <- list()
summary_list <- list()

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

  condition <- factor(
    meta_ct$condition,
    levels = c(
      "Veh",
      "EP",
      "EPR"
    )
  )

  n_veh <- sum(
    condition == "Veh"
  )

  n_ep <- sum(
    condition == "EP"
  )

  n_epr <- sum(
    condition == "EPR"
  )

  cat(
    "Adequate mice: Veh=",
    n_veh,
    " EP=",
    n_ep,
    " EPR=",
    n_epr,
    "\n",
    sep = ""
  )

  if (
    n_veh < 3L ||
    n_ep < 3L
  ) {
    stop(
      "Primary EP/Veh adequacy failed for ",
      ct
    )
  }

  if (n_epr < 2L) {
    stop(
      "Too few adequate EPR mice to retain the three-group model for ",
      ct
    )
  }

  dge <- DGEList(
    counts = counts_ct,
    group = condition
  )

  keep_gene <- filterByExpr(
    dge,
    group = condition,
    min.count = 10
  )

  dge <- dge[
    keep_gene,
    ,
    keep.lib.sizes = FALSE
  ]

  dge <- calcNormFactors(
    dge
  )

  design <- model.matrix(
    ~ 0 + condition
  )

  colnames(design) <- c(
    "Veh",
    "EP",
    "EPR"
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

  cm <- makeContrasts(
    EP_vs_Veh = EP - Veh,
    levels = design
  )

  fit2 <- contrasts.fit(
    fit,
    cm
  )

  fit2 <- eBayes(
    fit2
  )

  tt <- topTable(
    fit2,
    coef = "EP_vs_Veh",
    number = Inf,
    adjust.method = "BH",
    sort.by = "none"
  )

  # Independent BH validation.
  bh_check <- p.adjust(
    tt$P.Value,
    method = "BH"
  )

  max_delta <- max(
    abs(
      bh_check -
      tt$adj.P.Val
    ),
    na.rm = TRUE
  )

  if (
    !is.finite(max_delta) ||
    max_delta > 1e-12
  ) {
    stop(
      "BH validation failed for ",
      ct
    )
  }

  res <- data.frame(
    gene = rownames(tt),
    celltype = ct,
    contrast = "EP_vs_Veh",
    n_Veh = n_veh,
    n_EP = n_ep,
    n_EPR = n_epr,
    n_tests_within_celltype = nrow(tt),
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

  result_list[[ct]] <- res

  summary_list[[ct]] <- data.frame(
    celltype = ct,
    n_Veh = n_veh,
    n_EP = n_ep,
    n_EPR = n_epr,
    n_genes_tested = nrow(res),
    n_BH_FDR_lt_005 = sum(
      res$p_val_adj < 0.05
    ),
    max_BH_validation_delta = max_delta,
    stringsAsFactors = FALSE
  )
}

bulk_de <- rbindlist(
  result_list,
  use.names = TRUE,
  fill = TRUE
)

summary_dt <- rbindlist(
  summary_list,
  use.names = TRUE,
  fill = TRUE
)

# ============================================================================
# Global BH sensitivity
# ============================================================================

bulk_de[
  ,
  FDR_global_eligible :=
    p.adjust(
      p_val,
      method = "BH"
    )
]

bulk_de[
  ,
  global_BH_sig :=
    FDR_global_eligible < 0.05
]

setorder(
  bulk_de,
  celltype,
  p_val_adj,
  p_val
)

fwrite(
  bulk_de,
  file.path(
    audit_dir,
    "mouse_EP_vs_Veh_threegroup_fullprecision.csv"
  )
)

fwrite(
  summary_dt,
  file.path(
    audit_dir,
    "analysis_summary.csv"
  )
)

# ============================================================================
# Smad3 extraction
# ============================================================================

smad3 <- bulk_de[
  tolower(gene) == "smad3"
]

fwrite(
  smad3,
  file.path(
    audit_dir,
    "Smad3_EP_vs_Veh_threegroup.csv"
  )
)

cat(
  "\n===== ANALYSIS SUMMARY =====\n"
)

print(
  summary_dt,
  row.names = FALSE
)

cat(
  "\nTotal eligible gene x cell-type tests:",
  nrow(bulk_de),
  "\n"
)

cat(
  "Within-celltype BH significant tests:",
  sum(bulk_de$p_val_adj < 0.05),
  "\n"
)

cat(
  "Global-BH sensitivity significant tests:",
  sum(bulk_de$FDR_global_eligible < 0.05),
  "\n"
)

cat(
  "\n===== SMAD3 RESULTS =====\n"
)

print(
  smad3[, .(
    gene,
    celltype,
    n_Veh,
    n_EP,
    n_EPR,
    n_tests_within_celltype,
    avg_log2FC,
    p_val,
    p_val_adj,
    FDR_global_eligible,
    primary_DEG
  )],
  row.names = FALSE,
  digits = 10
)

cat(
  "\nSmad3 significant within-celltype BH:",
  sum(smad3$p_val_adj < 0.05),
  "/",
  nrow(smad3),
  "\n"
)

cat(
  "Smad3 significant global-BH sensitivity:",
  sum(smad3$FDR_global_eligible < 0.05),
  "/",
  nrow(smad3),
  "\n"
)

cat(
  "\nPASS: three-group full-precision mouse pseudobulk rebuild completed.\n"
)

cat(
  "Canonical files were NOT overwritten.\n"
)
