# ============================================================================
# Mouse GSE288662 pseudobulk DE rebuild — audit only
#
# Primary comparison:
#   EP vs Veh
#
# Biological replicate:
#   individual mouse (orig.ident)
#
# Eligibility:
#   >= 3 mice per condition with >= 5 nuclei for the cell type
#
# DE:
#   raw RNA counts -> pseudobulk -> edgeR filterByExpr -> voom-limma
#
# Primary multiple testing:
#   BH within each eligible cell-type-specific transcriptome-wide analysis
#
# Sensitivity:
#   BH across all eligible gene x cell-type tests
#
# Primary DEG definition:
#   within-cell-type BH FDR < 0.05
#   no hard fold-change threshold
#
# Audit only: canonical files are NOT overwritten.
# ============================================================================

rm(list = ls(all.names = TRUE))
gc()

suppressPackageStartupMessages({
  library(Seurat)
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
  "pseudobulk_rebuild_audit"
)

dir.create(
  audit_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("===== MOUSE PSEUDOBULK REBUILD =====\n")
cat("Object:", obj_file, "\n")
cat("Exists:", file.exists(obj_file), "\n\n")

if (!file.exists(obj_file)) {
  stop("Mouse Seurat object not found: ", obj_file)
}

# ============================================================================
# Load object and validate metadata
# ============================================================================

obj <- readRDS(obj_file)

required_meta <- c(
  "orig.ident",
  "condition",
  "type"
)

missing_meta <- setdiff(
  required_meta,
  colnames(obj@meta.data)
)

if (length(missing_meta) > 0L) {
  stop(
    "Missing metadata columns: ",
    paste(missing_meta, collapse = ", ")
  )
}

md <- obj@meta.data[, required_meta, drop = FALSE]

names(md) <- c(
  "sample",
  "condition",
  "celltype"
)

md$sample <- as.character(md$sample)
md$condition <- as.character(md$condition)
md$celltype <- as.character(md$celltype)

# Primary analysis is EP vs Veh only.
keep_cells <- rownames(md)[
  md$condition %in% c("Veh", "EP")
]

md <- md[
  keep_cells,
  ,
  drop = FALSE
]

cat("Cells in EP/Veh analysis:", nrow(md), "\n")

# --------------------------------------------------------------------------
# Confirm each mouse belongs to one condition
# --------------------------------------------------------------------------

sample_condition <- unique(
  md[, c("sample", "condition")]
)

if (
  any(
    table(sample_condition$sample) != 1L
  )
) {
  stop(
    "At least one sample maps to multiple conditions."
  )
}

cat("\n===== SAMPLE DESIGN =====\n")
print(
  sample_condition[
    order(
      sample_condition$condition,
      sample_condition$sample
    ),
  ],
  row.names = FALSE
)

cat("\nSamples per condition:\n")
print(
  table(sample_condition$condition)
)

if (
  sum(sample_condition$condition == "Veh") != 3L ||
  sum(sample_condition$condition == "EP") != 4L
) {
  stop(
    "Expected 3 Veh and 4 EP mice."
  )
}

# ============================================================================
# Sample x cell-type nuclei counts
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
  sample_condition,
  by = "sample",
  all.x = TRUE,
  sort = FALSE
)

observed_counts <- as.data.frame(
  table(
    sample = md$sample,
    celltype = md$celltype
  ),
  stringsAsFactors = FALSE
)

names(observed_counts)[3] <- "n_cells"

cell_qc <- merge(
  grid,
  observed_counts,
  by = c("sample", "celltype"),
  all.x = TRUE,
  sort = FALSE
)

cell_qc$n_cells[
  is.na(cell_qc$n_cells)
] <- 0L

cell_qc$adequate_nuclei <-
  cell_qc$n_cells >= 5L

cell_qc <- cell_qc[
  order(
    cell_qc$celltype,
    cell_qc$condition,
    cell_qc$sample
  ),
]

fwrite(
  cell_qc,
  file.path(
    audit_dir,
    "sample_celltype_nuclei_counts.csv"
  )
)

# --------------------------------------------------------------------------
# Eligibility: >=3 adequate mice in BOTH EP and Veh
# --------------------------------------------------------------------------

adequacy <- aggregate(
  adequate_nuclei ~ celltype + condition,
  data = cell_qc,
  FUN = sum
)

names(adequacy)[3] <- "n_adequate_samples"

adequacy_wide <- reshape(
  adequacy,
  idvar = "celltype",
  timevar = "condition",
  direction = "wide"
)

for (cn in c(
  "n_adequate_samples.EP",
  "n_adequate_samples.Veh"
)) {
  if (!cn %in% names(adequacy_wide)) {
    adequacy_wide[[cn]] <- 0L
  }
}

adequacy_wide$eligible <-
  adequacy_wide$n_adequate_samples.EP >= 3L &
  adequacy_wide$n_adequate_samples.Veh >= 3L

adequacy_wide <- adequacy_wide[
  order(adequacy_wide$celltype),
]

fwrite(
  adequacy_wide,
  file.path(
    audit_dir,
    "celltype_eligibility.csv"
  )
)

eligible_celltypes <-
  adequacy_wide$celltype[
    adequacy_wide$eligible
  ]

excluded_celltypes <-
  adequacy_wide$celltype[
    !adequacy_wide$eligible
  ]

cat("\n===== CELL-TYPE ELIGIBILITY =====\n")
print(
  adequacy_wide,
  row.names = FALSE
)

cat(
  "\nEligible:",
  length(eligible_celltypes),
  "\n"
)

cat(
  "Excluded:",
  length(excluded_celltypes),
  "\n"
)

cat(
  "Eligible cell types:\n",
  paste(
    eligible_celltypes,
    collapse = "\n"
  ),
  "\n",
  sep = ""
)

cat(
  "\nExcluded cell types:\n",
  paste(
    excluded_celltypes,
    collapse = "\n"
  ),
  "\n",
  sep = ""
)

if (length(eligible_celltypes) != 12L) {
  warning(
    "Expected 12 eligible cell types; found ",
    length(eligible_celltypes)
  )
}

# ============================================================================
# Memory-efficient pseudobulk aggregation
# ============================================================================

cat(
  "\n===== MEMORY-EFFICIENT RNA PSEUDOBULK AGGREGATION =====\n"
)

# Make RNA the default assay before reducing the object.
# The original object may have ATAC as DefaultAssay, and Seurat does not
# allow deletion of the current default assay.
DefaultAssay(obj) <- "RNA"

# Restrict first to the primary EP-vs-Veh comparison.
obj <- subset(
  obj,
  cells = keep_cells
)

gc()

# Retain RNA only. ATAC/SCT, reductions, and graphs are unnecessary for
# pseudobulk DE and consume substantial memory.
obj <- DietSeurat(
  obj,
  assays = "RNA",
  dimreducs = NULL,
  graphs = NULL
)

DefaultAssay(obj) <- "RNA"

gc()

# Align metadata with the actual cell order retained by Seurat.
md <- md[
  colnames(obj),
  ,
  drop = FALSE
]

if (!identical(
  rownames(md),
  colnames(obj)
)) {
  stop("Metadata/Seurat cell alignment failed.")
}

# --------------------------------------------------------------------------
# Create safe pseudobulk group IDs.
#
# IMPORTANT:
# We do NOT reconstruct metadata by splitting pseudobulk column names.
# This avoids problems caused by sample IDs containing both "-" and "_",
# e.g. 1583-EP versus 1696_EP.
# --------------------------------------------------------------------------

pb_key <- paste(
  md$sample,
  md$condition,
  md$celltype,
  sep = "|||"
)

unique_keys <- unique(pb_key)

key_map <- data.frame(
  pb_group = sprintf(
    "PB%03d",
    seq_along(unique_keys)
  ),
  pb_key = unique_keys,
  stringsAsFactors = FALSE
)

md$pb_group <- key_map$pb_group[
  match(
    pb_key,
    key_map$pb_key
  )
]

if (anyNA(md$pb_group)) {
  stop("Failed to assign pseudobulk group IDs.")
}

obj$pb_group <- md$pb_group

pb_meta <- unique(
  md[, c(
    "pb_group",
    "sample",
    "condition",
    "celltype"
  )]
)

if (anyDuplicated(pb_meta$pb_group)) {
  stop(
    "A pseudobulk ID maps to multiple metadata combinations."
  )
}

cat(
  "Cells retained:",
  ncol(obj),
  "\n"
)

cat(
  "Pseudobulk groups expected:",
  nrow(pb_meta),
  "\n"
)

# Join RNA count layers only after reducing the object.
obj[["RNA"]] <- JoinLayers(
  obj[["RNA"]]
)

gc()

# --------------------------------------------------------------------------
# Aggregate raw RNA counts using Seurat's pseudobulk implementation.
# --------------------------------------------------------------------------

pb_list <- AggregateExpression(
  obj,
  assays = "RNA",
  group.by = "pb_group",
  slot = "counts",
  return.seurat = FALSE,
  verbose = FALSE
)

pb_counts <- pb_list$RNA

cat(
  "Pseudobulk matrix:",
  nrow(pb_counts),
  "genes x",
  ncol(pb_counts),
  "groups\n"
)

cat(
  "First pseudobulk column names:\n"
)

print(
  head(
    colnames(pb_counts),
    20
  )
)

# Validate that Seurat preserved the safe IDs.
if (!setequal(
  colnames(pb_counts),
  pb_meta$pb_group
)) {

  cat(
    "\nExpected IDs:\n"
  )
  print(pb_meta$pb_group)

  cat(
    "\nObserved IDs:\n"
  )
  print(colnames(pb_counts))

  stop(
    "AggregateExpression pseudobulk IDs do not match metadata."
  )
}

pb_meta <- pb_meta[
  match(
    colnames(pb_counts),
    pb_meta$pb_group
  ),
  ,
  drop = FALSE
]

if (!identical(
  pb_meta$pb_group,
  colnames(pb_counts)
)) {
  stop(
    "Pseudobulk metadata alignment failed."
  )
}

# Attach nuclei counts / adequacy status.
pb_meta <- merge(
  pb_meta,
  cell_qc[, c(
    "sample",
    "condition",
    "celltype",
    "n_cells",
    "adequate_nuclei"
  )],
  by = c(
    "sample",
    "condition",
    "celltype"
  ),
  all.x = TRUE,
  sort = FALSE
)

pb_meta <- pb_meta[
  match(
    colnames(pb_counts),
    pb_meta$pb_group
  ),
  ,
  drop = FALSE
]

if (
  anyNA(pb_meta$n_cells) ||
  anyNA(pb_meta$adequate_nuclei)
) {
  stop(
    "Failed to attach sample x cell-type QC metadata."
  )
}

fwrite(
  pb_meta,
  file.path(
    audit_dir,
    "pseudobulk_metadata.csv"
  )
)

# We no longer need the full single-cell Seurat object.
rm(
  obj,
  pb_list
)

gc()

cat(
  "PASS: memory-efficient pseudobulk aggregation completed.\n"
)

# ============================================================================
# Differential expression per eligible cell type
# ============================================================================

result_list <- vector(
  "list",
  length(eligible_celltypes)
)

names(result_list) <-
  eligible_celltypes

analysis_summary <- vector(
  "list",
  length(eligible_celltypes)
)

names(analysis_summary) <-
  eligible_celltypes

for (ct in eligible_celltypes) {

  cat(
    "\n===== ",
    ct,
    " =====\n",
    sep = ""
  )

  use <- (
    pb_meta$celltype == ct &
    pb_meta$adequate_nuclei
  )

  meta_ct <- pb_meta[
    use,
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
      "EP"
    )
  )

  n_veh <- sum(
    condition == "Veh"
  )

  n_ep <- sum(
    condition == "EP"
  )

  cat(
    "Adequate mice: Veh =",
    n_veh,
    ", EP =",
    n_ep,
    "\n"
  )

  if (
    n_veh < 3L ||
    n_ep < 3L
  ) {
    stop(
      "Eligibility inconsistency for ",
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

  n_tested <- sum(
    keep_gene
  )

  cat(
    "Genes retained by filterByExpr:",
    n_tested,
    "\n"
  )

  if (n_tested < 1L) {
    stop(
      "No genes retained for ",
      ct
    )
  }

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
    "EP"
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

  contrast_matrix <- makeContrasts(
    EP_vs_Veh = EP - Veh,
    levels = design
  )

  fit2 <- contrasts.fit(
    fit,
    contrast_matrix
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

  # Independently validate limma's BH values.
  bh_check <- p.adjust(
    tt$P.Value,
    method = "BH"
  )

  max_bh_delta <- max(
    abs(
      bh_check -
      tt$adj.P.Val
    ),
    na.rm = TRUE
  )

  if (
    !is.finite(max_bh_delta) ||
    max_bh_delta > 1e-12
  ) {
    stop(
      "BH validation failed for ",
      ct,
      "; max delta = ",
      max_bh_delta
    )
  }

  res <- data.frame(
    gene = rownames(tt),
    celltype = ct,
    contrast = "EP_vs_Veh",
    n_Veh = n_veh,
    n_EP = n_ep,
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

  res$absLog2FCge1_sensitivity <-
    res$p_val_adj < 0.05 &
    abs(res$avg_log2FC) >= 1

  result_list[[ct]] <- res

  analysis_summary[[ct]] <- data.frame(
    celltype = ct,
    n_Veh = n_veh,
    n_EP = n_ep,
    n_genes_tested = nrow(res),
    n_BH_FDR_lt_005 = sum(
      res$p_val_adj < 0.05
    ),
    n_BH_FDR_lt_005_absLog2FCge1 = sum(
      res$p_val_adj < 0.05 &
      abs(res$avg_log2FC) >= 1
    ),
    max_BH_validation_delta = max_bh_delta,
    stringsAsFactors = FALSE
  )
}

bulk_de <- rbindlist(
  result_list,
  use.names = TRUE,
  fill = TRUE
)

analysis_summary <- rbindlist(
  analysis_summary,
  use.names = TRUE,
  fill = TRUE
)

# ============================================================================
# Global eligible-test BH sensitivity
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

# Sort only AFTER all full-precision calculations.
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
    "mouse_EP_vs_Veh_pseudobulk_fullprecision.csv"
  )
)

fwrite(
  analysis_summary,
  file.path(
    audit_dir,
    "mouse_EP_vs_Veh_analysis_summary.csv"
  )
)

# ============================================================================
# SMAD3 targeted extraction
# ============================================================================

smad3 <- bulk_de[
  tolower(gene) == "smad3"
]

fwrite(
  smad3,
  file.path(
    audit_dir,
    "Smad3_EP_vs_Veh_all_eligible_celltypes.csv"
  )
)

# ============================================================================
# Console report
# ============================================================================

cat(
  "\n===== ANALYSIS SUMMARY =====\n"
)

print(
  analysis_summary,
  row.names = FALSE
)

cat(
  "\nTotal eligible gene x cell-type tests:",
  nrow(bulk_de),
  "\n"
)

cat(
  "Primary DEGs (within-celltype BH FDR < 0.05):",
  sum(bulk_de$primary_DEG),
  "\n"
)

cat(
  "Global-BH sensitivity significant tests:",
  sum(bulk_de$global_BH_sig),
  "\n"
)

cat(
  "\n===== SMAD3 RESULTS =====\n"
)

if (nrow(smad3) == 0L) {

  cat(
    "Smad3 was not retained by filterByExpr in any eligible cell type.\n"
  )

} else {

  print(
    smad3[, .(
      gene,
      celltype,
      n_Veh,
      n_EP,
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
}

cat(
  "\nSMAD3 significant within-celltype BH:",
  sum(smad3$p_val_adj < 0.05),
  "/",
  nrow(smad3),
  "\n"
)

cat(
  "SMAD3 significant global-BH sensitivity:",
  sum(smad3$FDR_global_eligible < 0.05),
  "/",
  nrow(smad3),
  "\n"
)

cat(
  "\nPASS: mouse EP-vs-Veh pseudobulk rebuild completed.\n"
)

cat(
  "Canonical mouse DE files were NOT overwritten.\n"
)
