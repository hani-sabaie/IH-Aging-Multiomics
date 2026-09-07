# ============================================================================
# Revised human pseudobulk differential-expression candidate pipeline
#
# Purpose:
#   - Re-run pseudobulk DE directly from the annotated human Seurat object
#   - Apply explicit donor/cell-type adequacy QC
#   - Preserve full numerical precision
#   - Apply BH correction within each cell-type-specific transcriptome-wide test
#   - Define primary DEGs by BH-FDR < 0.05 (no hard fold-change threshold)
#
# IMPORTANT:
#   This script writes only to a revision-candidate directory.
#   It does NOT overwrite current canonical outputs.
# ============================================================================

rm(list = ls(all.names = TRUE))
gc()

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(edgeR)
  library(limma)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript script.R <annotated_human_object.rds>")
}

rds_file <- args[1]

outdir <- file.path(
  "processed_results",
  "02_differential_expression",
  "revised_candidate"
)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cat("\n============================================================\n")
cat("LOAD OBJECT\n")
cat("============================================================\n")

obj <- readRDS(rds_file)

required_meta <- c("sample", "condition", "skeletal_muscle")

missing_meta <- setdiff(
  required_meta,
  colnames(obj@meta.data)
)

if (length(missing_meta)) {
  stop(
    "Missing required metadata columns: ",
    paste(missing_meta, collapse = ", ")
  )
}

DefaultAssay(obj) <- "RNA"

cat("Cells:", ncol(obj), "\n")
cat("RNA features:", nrow(GetAssayData(obj, assay = "RNA", layer = "counts")), "\n")

# ============================================================================
# 1. Donor x cell-type adequacy
# ============================================================================

cat("\n============================================================\n")
cat("DONOR x CELL-TYPE ADEQUACY\n")
cat("============================================================\n")

md <- as.data.table(obj@meta.data)

cell_n <- md[
  ,
  .N,
  by = .(
    sample,
    condition,
    celltype = skeletal_muscle
  )
]

MIN_NUCLEI <- 5L
MIN_DONORS_PER_CONDITION <- 3L

all_celltypes <- sort(unique(cell_n$celltype))

eligibility <- rbindlist(
  lapply(all_celltypes, function(ct) {

    z <- cell_n[celltype == ct]

    ny <- uniqueN(
      z[
        condition == "Young" &
        N >= MIN_NUCLEI,
        sample
      ]
    )

    na <- uniqueN(
      z[
        condition == "Aged" &
        N >= MIN_NUCLEI,
        sample
      ]
    )

    data.table(
      celltype = ct,
      min_nuclei_per_donor = MIN_NUCLEI,
      Young_donors_meeting_threshold = ny,
      Aged_donors_meeting_threshold = na,
      eligible =
        ny >= MIN_DONORS_PER_CONDITION &
        na >= MIN_DONORS_PER_CONDITION
    )
  })
)

print(eligibility)

eligible_celltypes <- eligibility[
  eligible == TRUE,
  celltype
]

excluded_celltypes <- eligibility[
  eligible == FALSE,
  celltype
]

cat("\nEligible cell types:\n")
print(eligible_celltypes)

cat("\nExcluded from formal pseudobulk DE:\n")
print(excluded_celltypes)

fwrite(
  eligibility,
  file.path(outdir, "pseudobulk_celltype_eligibility.csv")
)

fwrite(
  cell_n,
  file.path(outdir, "donor_celltype_nucleus_counts.csv")
)

# ============================================================================
# 2. Pseudobulk aggregation
# ============================================================================

cat("\n============================================================\n")
cat("PSEUDOBULK AGGREGATION\n")
cat("============================================================\n")

pb_counts <- AggregateExpression(
  obj,
  assays = "RNA",
  group.by = c(
    "sample",
    "condition",
    "skeletal_muscle"
  ),
  slot = "counts",
  return.seurat = FALSE
)$RNA

cat(
  "Matrix:",
  nrow(pb_counts),
  "genes x",
  ncol(pb_counts),
  "pseudobulk groups\n"
)

# Seurat replaces underscores in sample identity values with dashes.
# Parse from the right using the known condition labels.

parsed <- str_match(
  colnames(pb_counts),
  "^(.+?)_(Young|Aged)_(.+)$"
)

if (any(is.na(parsed[, 1]))) {
  stop("Failed to parse one or more pseudobulk column names.")
}

pb_meta <- data.table(
  column = colnames(pb_counts),
  sample = parsed[, 2],
  condition = parsed[, 3],
  celltype = parsed[, 4]
)

pb_meta[, condition := factor(
  condition,
  levels = c("Young", "Aged")
)]

# ============================================================================
# 3. Cell-type-specific voom + limma
# ============================================================================

cat("\n============================================================\n")
cat("VOOM + LIMMA\n")
cat("============================================================\n")

res_list <- list()

for (ct in eligible_celltypes) {

  cat("Running:", ct, "\n")

  ix <- pb_meta$celltype == ct

  counts_ct <- pb_counts[
    ,
    ix,
    drop = FALSE
  ]

  meta_ct <- pb_meta[ix]

  cond <- droplevels(
    factor(
      meta_ct$condition,
      levels = c("Young", "Aged")
    )
  )

  if (length(levels(cond)) != 2) {
    stop("Condition factor problem for ", ct)
  }

  n_young <- sum(cond == "Young")
  n_aged  <- sum(cond == "Aged")

  if (
    n_young < MIN_DONORS_PER_CONDITION ||
    n_aged  < MIN_DONORS_PER_CONDITION
  ) {
    stop(
      "Unexpected insufficient pseudobulk replicates for eligible cell type: ",
      ct
    )
  }

  # Low-expression filtering
  dge <- DGEList(
    counts = counts_ct,
    group = cond
  )

  keep <- filterByExpr(
    dge,
    group = cond,
    min.count = 10
  )

  dge <- dge[
    keep,
    ,
    keep.lib.sizes = FALSE
  ]

  # TMM normalization
  dge <- calcNormFactors(dge)

  # Explicit Aged vs Young design
  design <- model.matrix(~ cond)

  colnames(design) <- c(
    "Intercept",
    "Aged_vs_Young"
  )

  v <- voom(
    dge,
    design = design,
    plot = FALSE
  )

  fit <- lmFit(
    v,
    design = design
  )

  fit <- eBayes(fit)

  tt <- topTable(
    fit,
    coef = "Aged_vs_Young",
    number = Inf,
    adjust.method = "BH",
    sort.by = "none"
  )

  tt[, "gene"] <- rownames(tt)

  z <- data.table(
    gene = tt$gene,
    celltype = ct,
    avg_log2FC = tt$logFC,
    AveExpr = tt$AveExpr,
    t = tt$t,
    p_val = tt$P.Value,
    p_val_adj = tt$adj.P.Val,
    B = tt$B,
    n_Young = n_young,
    n_Aged = n_aged
  )

  z[, direction := fifelse(
    avg_log2FC > 0,
    "Up",
    fifelse(
      avg_log2FC < 0,
      "Down",
      "Zero"
    )
  )]

  res_list[[ct]] <- z
}

all_de <- rbindlist(
  res_list,
  use.names = TRUE
)

setorder(
  all_de,
  celltype,
  p_val_adj,
  p_val
)

# ============================================================================
# 4. Primary DEG definition
# ============================================================================

primary_deg <- all_de[
  p_val_adj < 0.05
]

# FAP4 is excluded by the adequacy rule.
# Formal FAP DE therefore consists of FAP1-FAP3.

formal_fap_types <- intersect(
  c("FAP1", "FAP2", "FAP3", "FAP4"),
  eligible_celltypes
)

fap_deg <- primary_deg[
  celltype %in% formal_fap_types
]

# Effect-size-filtered sensitivity only; NOT primary DEG definition

effect_size_sensitivity <- primary_deg[
  abs(avg_log2FC) >= 1
]

fap_effect_size_sensitivity <- fap_deg[
  abs(avg_log2FC) >= 1
]

# ============================================================================
# 5. Summaries
# ============================================================================

summary_by_celltype <- primary_deg[
  ,
  .(
    DEGs = .N,
    Up = sum(avg_log2FC > 0),
    Down = sum(avg_log2FC < 0)
  ),
  by = celltype
][order(celltype)]

summary_total <- primary_deg[
  ,
  .(
    DEGs = .N,
    Up = sum(avg_log2FC > 0),
    Down = sum(avg_log2FC < 0)
  )
]

fap_summary <- fap_deg[
  ,
  .(
    DEGs = .N,
    Up = sum(avg_log2FC > 0),
    Down = sum(avg_log2FC < 0)
  ),
  by = celltype
][order(celltype)]

candidate_audit <- all_de[
  gene %in% c("SMAD3", "PLEKHA6") &
  celltype %in% c("FAP1", "FAP2", "FAP3"),
  .(
    gene,
    celltype,
    avg_log2FC,
    p_val,
    p_val_adj,
    direction,
    primary_DEG = p_val_adj < 0.05,
    passes_abs_log2FC_1 =
      p_val_adj < 0.05 &
      abs(avg_log2FC) >= 1
  )
]

# ============================================================================
# 6. Validation against 01D audit
# ============================================================================

audit_file <- file.path(
  "processed_results",
  "02_differential_expression",
  "full_precision_audit",
  "all_celltypes_full_precision.tsv"
)

validation <- NULL

if (file.exists(audit_file)) {

  old <- fread(audit_file)

  old <- old[
    celltype %in% eligible_celltypes,
    .(
      gene,
      celltype,
      old_log2FC = avg_log2FC,
      old_p = p_val,
      old_FDR = p_val_adj
    )
  ]

  new <- all_de[
    ,
    .(
      gene,
      celltype,
      new_log2FC = avg_log2FC,
      new_p = p_val,
      new_FDR = p_val_adj
    )
  ]

  cmp <- merge(
    old,
    new,
    by = c("gene", "celltype"),
    all = TRUE
  )

  validation <- data.table(
    old_rows = nrow(old),
    new_rows = nrow(new),
    shared_rows = sum(
      !is.na(cmp$old_p) &
      !is.na(cmp$new_p)
    ),
    old_only = sum(
      !is.na(cmp$old_p) &
      is.na(cmp$new_p)
    ),
    new_only = sum(
      is.na(cmp$old_p) &
      !is.na(cmp$new_p)
    ),
    max_abs_log2FC_difference = max(
      abs(cmp$old_log2FC - cmp$new_log2FC),
      na.rm = TRUE
    ),
    max_abs_P_difference = max(
      abs(cmp$old_p - cmp$new_p),
      na.rm = TRUE
    ),
    max_abs_FDR_difference = max(
      abs(cmp$old_FDR - cmp$new_FDR),
      na.rm = TRUE
    )
  )

  fwrite(
    cmp,
    file.path(
      outdir,
      "validation_against_01D_full_precision.tsv"
    ),
    sep = "\t"
  )

  fwrite(
    validation,
    file.path(
      outdir,
      "validation_summary.csv"
    )
  )
}

# ============================================================================
# 7. Write revision-candidate outputs
# ============================================================================

fwrite(
  all_de,
  file.path(
    outdir,
    "all_gene_results_full_precision.tsv"
  ),
  sep = "\t"
)

fwrite(
  primary_deg,
  file.path(
    outdir,
    "DEGs_BH_FDR005.csv"
  )
)

fwrite(
  fap_deg,
  file.path(
    outdir,
    "DEGs_FAP_BH_FDR005.csv"
  )
)

fwrite(
  effect_size_sensitivity,
  file.path(
    outdir,
    "DEGs_BH_FDR005_absLog2FCge1_sensitivity.csv"
  )
)

fwrite(
  fap_effect_size_sensitivity,
  file.path(
    outdir,
    "DEGs_FAP_BH_FDR005_absLog2FCge1_sensitivity.csv"
  )
)

fwrite(
  summary_by_celltype,
  file.path(
    outdir,
    "DEG_summary_by_celltype.csv"
  )
)

fwrite(
  fap_summary,
  file.path(
    outdir,
    "FAP_DEG_summary.csv"
  )
)

fwrite(
  candidate_audit,
  file.path(
    outdir,
    "SMAD3_PLEKHA6_full_precision.csv"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    outdir,
    "R_session_info.txt"
  )
)

# ============================================================================
# 8. Console report
# ============================================================================

cat("\n============================================================\n")
cat("PRIMARY DEG SUMMARY\n")
cat("============================================================\n\n")

print(summary_total)

cat("\nBY CELL TYPE:\n")
print(summary_by_celltype)

cat("\n============================================================\n")
cat("FORMAL FAP DEG SUMMARY\n")
cat("============================================================\n\n")

print(fap_summary)

cat("\nFAP total rows =", nrow(fap_deg), "\n")
cat(
  "FAP unique genes =",
  uniqueN(fap_deg$gene),
  "\n"
)

cat("\n============================================================\n")
cat("SMAD3 / PLEKHA6\n")
cat("============================================================\n\n")

print(candidate_audit)

cat("\n============================================================\n")
cat("VALIDATION AGAINST 01D\n")
cat("============================================================\n\n")

if (is.null(validation)) {
  cat("01D audit file not found.\n")
} else {
  print(validation)
}

cat("\n============================================================\n")
cat("OUTPUT DIRECTORY\n")
cat("============================================================\n\n")

cat(outdir, "\n")
