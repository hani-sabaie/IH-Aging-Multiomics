# ============================================================================
# Reviewer C7 / reproducibility audit
# Human pseudobulk DE full-precision reconstruction
#
# Goals:
#   1) Reconstruct the historical rounded analysis.
#   2) Re-run the same model using full-precision thresholds.
#   3) Apply the intended two-sided effect-size criterion:
#        BH FDR < 0.05 AND |log2FC| >= 1
#   4) Audit SMAD3 and PLEKHA6 explicitly.
#
# IMPORTANT:
#   This script does NOT overwrite canonical processed results.
# ============================================================================

rm(list = ls(all.names = TRUE))
gc()

suppressPackageStartupMessages({
  library(Seurat)
  library(stringr)
  library(dplyr)
  library(data.table)
  library(edgeR)
  library(limma)
})

# -------------------------------------------------------------------------
# Arguments / repository root
# -------------------------------------------------------------------------

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)

if (length(file_arg) == 1) {
  script_path <- sub("^--file=", "", file_arg)
  script_dir <- dirname(normalizePath(script_path))
  repo_root <- normalizePath(file.path(script_dir, ".."))
} else {
  repo_root <- normalizePath(".")
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop(
    "Usage: Rscript 01D_Pseudobulk_DE_FullPrecision_Audit.R ",
    "<annotated_human_rds>"
  )
}

rds_file <- args[1]

if (!file.exists(rds_file)) {
  stop("RDS not found: ", rds_file)
}

audit_dir <- file.path(
  repo_root,
  "processed_results",
  "02_differential_expression",
  "full_precision_audit"
)

dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

current_all_file <- file.path(
  repo_root,
  "processed_results",
  "02_differential_expression",
  "bulk_de_sig.csv"
)

current_fap_file <- file.path(
  repo_root,
  "processed_results",
  "02_differential_expression",
  "bulk_de_sig_faps.csv"
)

# -------------------------------------------------------------------------
# Load object
# -------------------------------------------------------------------------

cat("\n============================================================\n")
cat("LOAD ANNOTATED HUMAN OBJECT\n")
cat("============================================================\n")

cat("RDS:", rds_file, "\n")

obj <- readRDS(rds_file)

cat("Genes in object:", nrow(obj), "\n")
cat("Cells:", ncol(obj), "\n")

DefaultAssay(obj) <- "RNA"

# -------------------------------------------------------------------------
# Aggregate counts exactly as in historical analysis
# -------------------------------------------------------------------------

cat("\n============================================================\n")
cat("PSEUDOBULK AGGREGATION\n")
cat("============================================================\n")

pb_list <- AggregateExpression(
  obj,
  assays = "RNA",
  group.by = c("sample", "condition", "skeletal_muscle"),
  slot = "counts",
  return.seurat = FALSE
)

pb_counts <- pb_list$RNA

cat(
  "Pseudobulk matrix:",
  nrow(pb_counts), "genes x",
  ncol(pb_counts), "sample-condition-celltype groups\n"
)

# Seurat sanitizes sample underscores to dashes, e.g.
# Aged_1 -> Aged-1, producing Aged-1_Aged_FAP3.
#
# Parse with an explicit condition field rather than relying on unrestricted
# splitting.
parsed <- str_match(
  colnames(pb_counts),
  "^(.+?)_(Young|Aged)_(.+)$"
)

if (any(is.na(parsed[, 1]))) {
  bad <- colnames(pb_counts)[is.na(parsed[, 1])]
  stop(
    "Could not parse pseudobulk column names:\n",
    paste(bad, collapse = "\n")
  )
}

meta <- data.frame(
  sample = parsed[, 2],
  condition = factor(parsed[, 3], levels = c("Young", "Aged")),
  celltype = parsed[, 4],
  stringsAsFactors = FALSE
)

rownames(meta) <- colnames(pb_counts)

cat("\nReplicate counts per cell type and condition:\n")
rep_tab <- with(meta, table(celltype, condition))
print(rep_tab)

fwrite(
  as.data.table(as.data.frame.matrix(rep_tab), keep.rownames = "celltype"),
  file.path(audit_dir, "pseudobulk_replicates_by_celltype.csv")
)

# -------------------------------------------------------------------------
# Historical voom + limma model
# -------------------------------------------------------------------------

run_voomlimma <- function(count_mat, cond) {

  dge <- DGEList(count_mat, group = cond)
  dge <- calcNormFactors(dge)

  design <- model.matrix(~ cond)

  vm <- voom(
    dge,
    design = design,
    plot = FALSE
  )

  fit <- lmFit(vm, design = design)
  fit <- eBayes(fit)

  tt <- topTable(
    fit,
    n = Inf,
    adjust.method = "BH",
    sort.by = "P"
  )

  tt
}

subtypes <- sort(unique(meta$celltype))

cat("\n============================================================\n")
cat("VOOM + LIMMA BY CELL TYPE\n")
cat("============================================================\n")

res_list <- lapply(subtypes, function(st) {

  cat("Running:", st, "\n")

  keep_cols <- meta$celltype == st

  counts_ct <- pb_counts[, keep_cols, drop = FALSE]

  cond_ct <- droplevels(meta$condition[keep_cols])
  cond_ct <- relevel(cond_ct, ref = "Young")

  n_young <- sum(cond_ct == "Young")
  n_aged  <- sum(cond_ct == "Aged")

  if (
    length(unique(cond_ct)) < 2 ||
    any(table(cond_ct) < 2)
  ) {
    warning(
      "Skipping ", st,
      " because of insufficient biological replicates."
    )
    return(NULL)
  }

  dge_tmp <- DGEList(counts_ct)

  keep_genes <- filterByExpr(
    dge_tmp,
    group = cond_ct,
    min.count = 10
  )

  counts_ct <- counts_ct[keep_genes, , drop = FALSE]

  tt <- run_voomlimma(
    count_mat = counts_ct,
    cond = cond_ct
  )

  out <- data.table(
    gene = rownames(tt),
    p_val = tt$P.Value,
    p_val_adj = tt$adj.P.Val,
    avg_log2FC = tt$logFC,
    celltype = st,
    n_Young = n_young,
    n_Aged = n_aged,
    n_genes_tested_celltype = nrow(tt)
  )

  # Historical values after rounding, exactly matching the old script logic.
  out[, p_val_round4 := round(p_val, 4)]
  out[, p_val_adj_round4 := round(p_val_adj, 4)]
  out[, avg_log2FC_round2 := round(avg_log2FC, 2)]

  # Historical committed selection logic.
  out[, pass_historical_rounded :=
        p_val_adj_round4 < 0.05 &
        avg_log2FC_round2 >= 1]

  # Same one-direction analysis but threshold BEFORE rounding.
  out[, pass_positive_fullprecision :=
        p_val_adj < 0.05 &
        avg_log2FC >= 1]

  # Intended Methods criterion:
  # two-sided effect-size threshold and full-precision statistics.
  out[, pass_abs_fullprecision :=
        p_val_adj < 0.05 &
        abs(avg_log2FC) >= 1]

  out
})

res_list <- Filter(Negate(is.null), res_list)

bulk_full <- rbindlist(res_list, use.names = TRUE)

setorder(
  bulk_full,
  p_val_adj,
  p_val,
  celltype,
  gene
)

# -------------------------------------------------------------------------
# Load current committed result sets
# -------------------------------------------------------------------------

current_all <- fread(current_all_file)
current_fap <- fread(current_fap_file)

current_all_key <- paste(
  current_all$gene,
  current_all$celltype,
  sep = "||"
)

current_fap_key <- paste(
  current_fap$gene,
  current_fap$celltype,
  sep = "||"
)

bulk_full[, key := paste(gene, celltype, sep = "||")]

bulk_full[, current_all_member := key %in% current_all_key]
bulk_full[, current_fap_member := key %in% current_fap_key]

# -------------------------------------------------------------------------
# Create result sets
# -------------------------------------------------------------------------

hist_rounded <- bulk_full[
  pass_historical_rounded == TRUE
]

positive_full <- bulk_full[
  pass_positive_fullprecision == TRUE
]

abs_full <- bulk_full[
  pass_abs_fullprecision == TRUE
]

fap_names <- c("FAP1", "FAP2", "FAP3", "FAP4")

hist_rounded_fap <- hist_rounded[celltype %in% fap_names]
positive_full_fap <- positive_full[celltype %in% fap_names]
abs_full_fap <- abs_full[celltype %in% fap_names]

# -------------------------------------------------------------------------
# Summary
# -------------------------------------------------------------------------

summary_dt <- rbindlist(
  lapply(subtypes, function(st) {

    x <- bulk_full[celltype == st]

    data.table(
      celltype = st,
      genes_tested = nrow(x),
      historical_rounded = sum(x$pass_historical_rounded),
      positive_fullprecision = sum(x$pass_positive_fullprecision),
      abs_fullprecision = sum(x$pass_abs_fullprecision),
      abs_fullprecision_up = sum(
        x$pass_abs_fullprecision &
        x$avg_log2FC >= 1
      ),
      abs_fullprecision_down = sum(
        x$pass_abs_fullprecision &
        x$avg_log2FC <= -1
      ),
      rounded_only = sum(
        x$pass_historical_rounded &
        !x$pass_positive_fullprecision
      ),
      fullprecision_positive_only = sum(
        !x$pass_historical_rounded &
        x$pass_positive_fullprecision
      )
    )
  }),
  use.names = TRUE
)

# Add TOTAL
summary_dt <- rbind(
  summary_dt,
  data.table(
    celltype = "TOTAL",
    genes_tested = nrow(bulk_full),
    historical_rounded = nrow(hist_rounded),
    positive_fullprecision = nrow(positive_full),
    abs_fullprecision = nrow(abs_full),
    abs_fullprecision_up = sum(abs_full$avg_log2FC >= 1),
    abs_fullprecision_down = sum(abs_full$avg_log2FC <= -1),
    rounded_only = sum(
      bulk_full$pass_historical_rounded &
      !bulk_full$pass_positive_fullprecision
    ),
    fullprecision_positive_only = sum(
      !bulk_full$pass_historical_rounded &
      bulk_full$pass_positive_fullprecision
    )
  ),
  fill = TRUE
)

# -------------------------------------------------------------------------
# Compare historical reconstruction with committed files
# -------------------------------------------------------------------------

hist_key <- hist_rounded$key
hist_fap_key <- hist_rounded_fap$key

comparison <- data.table(
  comparison = c(
    "Current all DEG rows",
    "Historical-rounded reconstructed all",
    "Overlap: current all vs reconstructed",
    "Current-only all",
    "Reconstructed-only all",
    "Current FAP DEG rows",
    "Historical-rounded reconstructed FAP",
    "Overlap: current FAP vs reconstructed",
    "Current-only FAP",
    "Reconstructed-only FAP"
  ),
  n = c(
    length(current_all_key),
    length(hist_key),
    length(intersect(current_all_key, hist_key)),
    length(setdiff(current_all_key, hist_key)),
    length(setdiff(hist_key, current_all_key)),
    length(current_fap_key),
    length(hist_fap_key),
    length(intersect(current_fap_key, hist_fap_key)),
    length(setdiff(current_fap_key, hist_fap_key)),
    length(setdiff(hist_fap_key, current_fap_key))
  )
)

# -------------------------------------------------------------------------
# Borderline genes affected by pre-filter rounding
# -------------------------------------------------------------------------

rounding_changes <- bulk_full[
  pass_historical_rounded != pass_positive_fullprecision |
  (
    p_val_adj < 0.055 &
    p_val_adj > 0.045
  ) |
  (
    abs(avg_log2FC) >= 0.95 &
    abs(avg_log2FC) <= 1.05
  )
][
  order(celltype, gene)
]

# -------------------------------------------------------------------------
# Candidate audit: SMAD3 and PLEKHA6
# -------------------------------------------------------------------------

candidate_audit <- bulk_full[
  gene %in% c("SMAD3", "PLEKHA6")
][
  order(gene, celltype)
]

# -------------------------------------------------------------------------
# Write audit outputs
# -------------------------------------------------------------------------

fwrite(
  bulk_full,
  file.path(audit_dir, "all_celltypes_full_precision.tsv"),
  sep = "\t"
)

fwrite(
  hist_rounded,
  file.path(audit_dir, "historical_rounded_reconstruction.csv")
)

fwrite(
  positive_full,
  file.path(audit_dir, "positive_fullprecision_FDR005_logFCge1.csv")
)

fwrite(
  abs_full,
  file.path(audit_dir, "intended_fullprecision_FDR005_absLogFCge1.csv")
)

fwrite(
  hist_rounded_fap,
  file.path(audit_dir, "historical_rounded_FAP.csv")
)

fwrite(
  positive_full_fap,
  file.path(audit_dir, "positive_fullprecision_FAP.csv")
)

fwrite(
  abs_full_fap,
  file.path(audit_dir, "intended_fullprecision_FAP.csv")
)

fwrite(
  summary_dt,
  file.path(audit_dir, "DEG_count_summary.csv")
)

fwrite(
  comparison,
  file.path(audit_dir, "historical_vs_committed_comparison.csv")
)

fwrite(
  rounding_changes,
  file.path(audit_dir, "borderline_and_rounding_audit.csv")
)

fwrite(
  candidate_audit,
  file.path(audit_dir, "SMAD3_PLEKHA6_full_precision_audit.csv")
)

# -------------------------------------------------------------------------
# Console report
# -------------------------------------------------------------------------

cat("\n============================================================\n")
cat("DEG COUNT SUMMARY\n")
cat("============================================================\n\n")
print(summary_dt)

cat("\n============================================================\n")
cat("HISTORICAL RECONSTRUCTION VS COMMITTED FILES\n")
cat("============================================================\n\n")
print(comparison)

cat("\n============================================================\n")
cat("FAP COUNTS\n")
cat("============================================================\n\n")

cat("Historical rounded FAP:\n")
print(table(hist_rounded_fap$celltype))

cat("\nPositive full-precision FAP:\n")
print(table(positive_full_fap$celltype))

cat("\nIntended |logFC| full-precision FAP:\n")
print(table(abs_full_fap$celltype))

cat("\n============================================================\n")
cat("SMAD3 / PLEKHA6 FULL-PRECISION RESULTS\n")
cat("============================================================\n\n")

print(
  candidate_audit[
    celltype %in% fap_names,
    .(
      gene,
      celltype,
      p_val,
      p_val_adj,
      avg_log2FC,
      p_val_round4,
      p_val_adj_round4,
      avg_log2FC_round2,
      pass_historical_rounded,
      pass_positive_fullprecision,
      pass_abs_fullprecision,
      current_fap_member
    )
  ]
)

cat("\n============================================================\n")
cat("DOWNREGULATED DEGs ADDED BY |logFC| CRITERION\n")
cat("============================================================\n\n")

down <- abs_full[avg_log2FC <= -1]

cat("N =", nrow(down), "\n")

if (nrow(down) > 0) {
  print(
    down[
      order(p_val_adj),
      .(
        gene,
        celltype,
        avg_log2FC,
        p_val_adj
      )
    ][1:min(.N, 30)]
  )
}

cat("\nAudit outputs written to:\n")
cat(audit_dir, "\n")

cat("\nSession info:\n")
print(sessionInfo())
