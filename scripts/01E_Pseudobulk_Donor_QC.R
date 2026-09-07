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
if (length(args) < 1) stop("Provide annotated human RDS path.")
rds_file <- args[1]

obj <- readRDS(rds_file)
DefaultAssay(obj) <- "RNA"

# ------------------------------------------------------------
# Pseudobulk
# ------------------------------------------------------------

pb <- AggregateExpression(
  obj,
  assays = "RNA",
  group.by = c("sample", "condition", "skeletal_muscle"),
  slot = "counts",
  return.seurat = FALSE
)$RNA

parsed <- str_match(
  colnames(pb),
  "^(.+?)_(Young|Aged)_(.+)$"
)

meta <- data.table(
  column = colnames(pb),
  sample = parsed[,2],
  condition = parsed[,3],
  celltype = parsed[,4]
)

# Cell numbers per donor/cell type
md <- as.data.table(obj@meta.data)

cell_n <- md[
  ,
  .N,
  by = .(sample, condition, skeletal_muscle)
]

cell_n[, sample_clean := gsub("_", "-", sample, fixed = TRUE)]

# ------------------------------------------------------------
# FAP3 donor-level expression
# ------------------------------------------------------------

ix <- meta$celltype == "FAP3"

cnt <- pb[, ix, drop = FALSE]
m <- copy(meta[ix])

cond <- factor(m$condition, levels = c("Young", "Aged"))

dge <- DGEList(cnt, group = cond)

keep <- filterByExpr(
  dge,
  group = cond,
  min.count = 10
)

dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)

cpm_raw <- cpm(
  dge,
  normalized.lib.sizes = TRUE,
  log = FALSE
)

logcpm <- cpm(
  dge,
  normalized.lib.sizes = TRUE,
  log = TRUE,
  prior.count = 0.25
)

genes <- c("SMAD3", "PLEKHA6")

candidate_rows <- rbindlist(lapply(genes, function(g) {

  if (!g %in% rownames(dge)) return(NULL)

  x <- data.table(
    gene = g,
    sample = m$sample,
    condition = m$condition,
    raw_count = as.numeric(dge$counts[g, ]),
    CPM = as.numeric(cpm_raw[g, ]),
    logCPM = as.numeric(logcpm[g, ]),
    library_size = dge$samples$lib.size,
    norm_factor = dge$samples$norm.factors
  )

  x
}))

candidate_rows <- merge(
  candidate_rows,
  cell_n[
    skeletal_muscle == "FAP3",
    .(
      sample = sample_clean,
      condition,
      nuclei = N
    )
  ],
  by = c("sample", "condition"),
  all.x = TRUE
)

setorder(candidate_rows, gene, condition, sample)

# ------------------------------------------------------------
# Leave-one-donor-out FAP3 DE
# ------------------------------------------------------------

run_de <- function(count_mat, cond_vec, candidate) {

  d <- DGEList(count_mat, group = cond_vec)

  keep <- filterByExpr(
    d,
    group = cond_vec,
    min.count = 10
  )

  d <- d[keep, , keep.lib.sizes = FALSE]

  if (!candidate %in% rownames(d)) {
    return(data.table(
      gene = candidate,
      logFC = NA_real_,
      P = NA_real_,
      FDR = NA_real_
    ))
  }

  d <- calcNormFactors(d)

  design <- model.matrix(~ cond_vec)

  v <- voom(d, design, plot = FALSE)
  fit <- eBayes(lmFit(v, design))

  tt <- topTable(
    fit,
    coef = 2,
    number = Inf,
    adjust.method = "BH",
    sort.by = "none"
  )

  data.table(
    gene = candidate,
    logFC = tt[candidate, "logFC"],
    P = tt[candidate, "P.Value"],
    FDR = tt[candidate, "adj.P.Val"]
  )
}

loo <- list()

# Full model
for (g in genes) {
  z <- run_de(cnt, cond, g)
  z[, excluded_sample := "NONE_FULL_MODEL"]
  z[, n_Young := sum(cond == "Young")]
  z[, n_Aged := sum(cond == "Aged")]
  loo[[length(loo) + 1]] <- z
}

# Leave one donor out
for (j in seq_len(ncol(cnt))) {

  cnt_j <- cnt[, -j, drop = FALSE]
  cond_j <- droplevels(cond[-j])

  for (g in genes) {

    z <- run_de(cnt_j, cond_j, g)

    z[, excluded_sample := m$sample[j]]
    z[, n_Young := sum(cond_j == "Young")]
    z[, n_Aged := sum(cond_j == "Aged")]

    loo[[length(loo) + 1]] <- z
  }
}

loo <- rbindlist(loo)

# ------------------------------------------------------------
# FAP4 donor QC
# ------------------------------------------------------------

fap4_meta <- meta[celltype == "FAP4"]

fap4_qc <- data.table(
  sample = fap4_meta$sample,
  condition = fap4_meta$condition,
  library_size = colSums(
    pb[, fap4_meta$column, drop = FALSE]
  )
)

fap4_qc <- merge(
  fap4_qc,
  cell_n[
    skeletal_muscle == "FAP4",
    .(
      sample = sample_clean,
      condition,
      nuclei = N
    )
  ],
  by = c("sample", "condition"),
  all.x = TRUE
)

# ------------------------------------------------------------
# Output
# ------------------------------------------------------------

outdir <- file.path(
  "processed_results",
  "02_differential_expression",
  "full_precision_audit"
)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

fwrite(
  candidate_rows,
  file.path(outdir, "FAP3_SMAD3_PLEKHA6_donor_expression.csv")
)

fwrite(
  loo,
  file.path(outdir, "FAP3_SMAD3_PLEKHA6_leave_one_out.csv")
)

fwrite(
  fap4_qc,
  file.path(outdir, "FAP4_donor_QC.csv")
)

cat("\n============================================================\n")
cat("FAP3 DONOR-LEVEL EXPRESSION\n")
cat("============================================================\n")
print(candidate_rows)

cat("\n============================================================\n")
cat("LEAVE-ONE-DONOR-OUT\n")
cat("============================================================\n")
print(loo)

cat("\n============================================================\n")
cat("LEAVE-ONE-OUT SUMMARY\n")
cat("============================================================\n")

print(
  loo[
    excluded_sample != "NONE_FULL_MODEL",
    .(
      min_logFC = min(logFC, na.rm = TRUE),
      max_logFC = max(logFC, na.rm = TRUE),
      median_logFC = median(logFC, na.rm = TRUE),
      min_FDR = min(FDR, na.rm = TRUE),
      max_FDR = max(FDR, na.rm = TRUE),
      runs_FDR_lt_005 = sum(FDR < 0.05, na.rm = TRUE),
      runs_absFC_ge_1 = sum(abs(logFC) >= 1, na.rm = TRUE),
      runs_both = sum(
        FDR < 0.05 & abs(logFC) >= 1,
        na.rm = TRUE
      )
    ),
    by = gene
  ]
)

cat("\n============================================================\n")
cat("FAP4 DONOR QC\n")
cat("============================================================\n")
print(fap4_qc)
