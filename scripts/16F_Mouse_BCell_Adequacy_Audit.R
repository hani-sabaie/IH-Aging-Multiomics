# ============================================================================
# Mouse B-cell targeted audit
#
# Purpose:
#   Determine whether the B-cell differences between the historical canonical
#   analysis and 16D are caused by excluding EPR mouse 1693_EPR, which has
#   only 4 B cells.
#
# Model:
#   Veh + EP + EPR, contrast EP - Veh
#
# A) all original pseudobulk samples
# B) exclude 1693_EPR
#
# Full transcriptome + BH within B Cell.
# Audit only.
# ============================================================================

rm(list = ls(all.names = TRUE))
gc()

suppressPackageStartupMessages({
  library(Seurat)
  library(edgeR)
  library(limma)
  library(data.table)
})

obj_file <- paste0(
  "C:/Users/Hani/Desktop/Hernia/data/GSE288662/",
  "Processed_Seurat_Object/",
  "GSE288662_Processed_Seurat_Object.rds"
)

if (!file.exists(obj_file)) {
  stop("Mouse object not found.")
}

obj <- readRDS(obj_file)

md <- obj@meta.data

cells <- rownames(md)[
  md$type == "B Cell"
]

cat("===== B-CELL ADEQUACY AUDIT =====\n")
cat("B cells:", length(cells), "\n\n")

DefaultAssay(obj) <- "RNA"

obj <- subset(
  obj,
  cells = cells
)

obj <- DietSeurat(
  obj,
  assays = "RNA",
  dimreducs = NULL,
  graphs = NULL
)

DefaultAssay(obj) <- "RNA"

md <- obj@meta.data[
  ,
  c(
    "orig.ident",
    "condition"
  ),
  drop = FALSE
]

names(md) <- c(
  "sample",
  "condition"
)

md$sample <- as.character(md$sample)
md$condition <- as.character(md$condition)

cat("===== B CELLS PER MOUSE =====\n")

tab <- as.data.frame(
  table(
    sample = md$sample,
    condition = md$condition
  )
)

tab <- tab[
  tab$Freq > 0,
]

print(
  tab,
  row.names = FALSE
)

# Safe sample IDs
samples <- unique(md$sample)

sample_map <- data.frame(
  sample = samples,
  pb_group = sprintf(
    "PB%02d",
    seq_along(samples)
  ),
  stringsAsFactors = FALSE
)

md$pb_group <- sample_map$pb_group[
  match(
    md$sample,
    sample_map$sample
  )
]

obj$pb_group <- md$pb_group

obj[["RNA"]] <- JoinLayers(
  obj[["RNA"]]
)

pb_counts <- AggregateExpression(
  obj,
  assays = "RNA",
  group.by = "pb_group",
  slot = "counts",
  return.seurat = FALSE,
  verbose = FALSE
)$RNA

pb_meta <- unique(
  md[, c(
    "pb_group",
    "sample",
    "condition"
  )]
)

pb_meta <- pb_meta[
  match(
    colnames(pb_counts),
    pb_meta$pb_group
  ),
]

rm(obj)
gc()

run_model <- function(label, exclude = character()) {

  meta_x <- pb_meta[
    !pb_meta$sample %in% exclude,
    ,
    drop = FALSE
  ]

  counts_x <- pb_counts[
    ,
    meta_x$pb_group,
    drop = FALSE
  ]

  cond <- factor(
    meta_x$condition,
    levels = c(
      "Veh",
      "EP",
      "EPR"
    )
  )

  dge <- DGEList(
    counts = counts_x,
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

  dge <- calcNormFactors(dge)

  design <- model.matrix(
    ~ 0 + cond
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

  fit <- contrasts.fit(
    fit,
    cm
  )

  fit <- eBayes(fit)

  tt <- topTable(
    fit,
    coef = "EP_vs_Veh",
    number = Inf,
    adjust.method = "BH",
    sort.by = "none"
  )

  data.table(
    gene = rownames(tt),
    model = label,
    n_Veh = sum(cond == "Veh"),
    n_EP = sum(cond == "EP"),
    n_EPR = sum(cond == "EPR"),
    n_genes_tested = nrow(tt),
    log2FC = tt$logFC,
    p_val = tt$P.Value,
    FDR = tt$adj.P.Val
  )
}

all_model <- run_model(
  "A_all_original_samples"
)

filtered_model <- run_model(
  "B_exclude_1693_EPR",
  exclude = "1693_EPR"
)

cat("\n===== MODEL SUMMARY =====\n")

print(
  rbind(
    all_model[, .(
      model = unique(model),
      n_Veh = unique(n_Veh),
      n_EP = unique(n_EP),
      n_EPR = unique(n_EPR),
      n_genes_tested = unique(n_genes_tested),
      n_FDR_lt_005 = sum(FDR < 0.05)
    )],
    filtered_model[, .(
      model = unique(model),
      n_Veh = unique(n_Veh),
      n_EP = unique(n_EP),
      n_EPR = unique(n_EPR),
      n_genes_tested = unique(n_genes_tested),
      n_FDR_lt_005 = sum(FDR < 0.05)
    )]
  )
)

# --------------------------------------------------------------------------
# Compare all-original model with historical canonical B Cell
# --------------------------------------------------------------------------

old <- fread(
  file.path(
    "processed_results",
    "13_mouse_validation",
    "bulk_de_all_contrasts_mouse.csv"
  )
)

old <- old[
  celltype == "B Cell" &
    contrast == "Veh_vs_EP"
]

old[, log2FC_EP_vs_Veh := -avg_log2FC]

cmp <- merge(
  old[, .(
    gene,
    old_log2FC = log2FC_EP_vs_Veh,
    old_p = p_val,
    old_FDR = p_val_adj
  )],
  all_model[, .(
    gene,
    new_log2FC = log2FC,
    new_p = p_val,
    new_FDR = FDR
  )],
  by = "gene"
)

cmp[
  ,
  old_sig := old_FDR < 0.05
]

cmp[
  ,
  new_sig := new_FDR < 0.05
]

cat("\n===== ORIGINAL-DESIGN VS CANONICAL =====\n")

cat(
  "Shared genes:",
  nrow(cmp),
  "\n"
)

cat(
  "Median |delta logFC|:",
  median(
    abs(
      cmp$old_log2FC -
        cmp$new_log2FC
    )
  ),
  "\n"
)

cat(
  "Median |delta P|:",
  median(
    abs(
      cmp$old_p -
        cmp$new_p
    )
  ),
  "\n"
)

cat(
  "Median |delta FDR|:",
  median(
    abs(
      cmp$old_FDR -
        cmp$new_FDR
    )
  ),
  "\n"
)

cat(
  "Significance-status changes:",
  sum(
    cmp$old_sig !=
      cmp$new_sig
  ),
  "\n"
)

cat("\n===== SMAD3 =====\n")

print(
  cmp[
    tolower(gene) == "smad3"
  ],
  digits = 10
)

outdir <- file.path(
  "processed_results",
  "13_mouse_validation",
  "three_group_fullprecision"
)

fwrite(
  all_model,
  file.path(
    outdir,
    "BCell_all_original_samples_fullprecision.csv"
  )
)

fwrite(
  filtered_model,
  file.path(
    outdir,
    "BCell_exclude_1693_EPR_fullprecision.csv"
  )
)

cat(
  "\nPASS: B-cell adequacy audit completed.\n"
)
