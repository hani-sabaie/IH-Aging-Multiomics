# ============================================================================
# Targeted audit: why did Smad3 Type IIa FDR change?
#
# Compares:
#   A) correct metadata, EP + Veh + EPR
#   B) correct metadata, EP + Veh only
#   C) historical-parser equivalent: drop underscore-named mice, 3 groups
#   D) EP + Veh only, but drop 1696_EP
#
# Full transcriptome is tested in Type IIa so BH denominator is preserved.
# Audit only. No canonical files are modified.
# ============================================================================

rm(list = ls(all.names = TRUE))
gc()

suppressPackageStartupMessages({
  library(Seurat)
  library(edgeR)
  library(limma)
})

obj_file <- paste0(
  "C:/Users/Hani/Desktop/Hernia/data/GSE288662/",
  "Processed_Seurat_Object/",
  "GSE288662_Processed_Seurat_Object.rds"
)

if (!file.exists(obj_file)) {
  stop("Mouse object not found.")
}

cat("===== TYPE IIA SMAD3 MODEL AUDIT =====\n")

obj <- readRDS(obj_file)

md0 <- obj@meta.data

needed <- c(
  "orig.ident",
  "condition",
  "type"
)

if (!all(needed %in% names(md0))) {
  stop("Required metadata missing.")
}

# --------------------------------------------------------------------------
# Keep only Type IIa cells immediately to minimize memory use.
# --------------------------------------------------------------------------

cells_use <- rownames(md0)[
  md0$type == "Type IIa Myofiber"
]

cat("Type IIa cells:", length(cells_use), "\n")

DefaultAssay(obj) <- "RNA"

obj <- subset(
  obj,
  cells = cells_use
)

gc()

obj <- DietSeurat(
  obj,
  assays = "RNA",
  dimreducs = NULL,
  graphs = NULL
)

DefaultAssay(obj) <- "RNA"

md <- obj@meta.data[, needed, drop = FALSE]

names(md) <- c(
  "sample",
  "condition",
  "celltype"
)

md$sample <- as.character(md$sample)
md$condition <- as.character(md$condition)

# --------------------------------------------------------------------------
# Show historical parser behavior.
# --------------------------------------------------------------------------

sample_design <- unique(
  md[, c(
    "sample",
    "condition"
  )]
)

sample_design <- sample_design[
  order(
    sample_design$condition,
    sample_design$sample
  ),
]

cat("\n===== TRUE SAMPLE DESIGN =====\n")
print(
  sample_design,
  row.names = FALSE
)

cat("\n===== HISTORICAL PARSER CHECK =====\n")

parser_check <- do.call(
  rbind,
  lapply(
    seq_len(nrow(sample_design)),
    function(i) {

      sample <- sample_design$sample[i]
      cond <- sample_design$condition[i]

      old_name <- paste(
        sample,
        cond,
        "Type IIa Myofiber",
        sep = "_"
      )

      sp <- strsplit(
        old_name,
        "_",
        fixed = TRUE
      )[[1]]

      # Equivalent to str_split_fixed(..., "_", 3):
      parsed_sample <- sp[1]
      parsed_condition <- sp[2]
      parsed_celltype <- paste(
        sp[3:length(sp)],
        collapse = "_"
      )

      data.frame(
        sample = sample,
        condition = cond,
        historical_name = old_name,
        parsed_sample = parsed_sample,
        parsed_condition = parsed_condition,
        parsed_celltype = parsed_celltype,
        parser_correct =
          parsed_sample == sample &&
          parsed_condition == cond &&
          parsed_celltype == "Type IIa Myofiber",
        stringsAsFactors = FALSE
      )
    }
  )
)

print(
  parser_check,
  row.names = FALSE
)

bad_samples <- parser_check$sample[
  !parser_check$parser_correct
]

cat(
  "\nSamples misparsed historically:",
  paste(bad_samples, collapse = ", "),
  "\n"
)

# --------------------------------------------------------------------------
# Safe pseudobulk IDs
# --------------------------------------------------------------------------

unique_samples <- unique(md$sample)

sample_map <- data.frame(
  sample = unique_samples,
  pb_group = sprintf(
    "PB%02d",
    seq_along(unique_samples)
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

gc()

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

if (!identical(
  pb_meta$pb_group,
  colnames(pb_counts)
)) {
  stop("Pseudobulk metadata alignment failed.")
}

cat(
  "\nPseudobulk matrix:",
  nrow(pb_counts),
  "genes x",
  ncol(pb_counts),
  "mice\n"
)

rm(obj)
gc()

# ============================================================================
# DE helper
# ============================================================================

run_model <- function(
  label,
  keep_samples,
  condition_levels
) {

  meta_x <- pb_meta[
    pb_meta$sample %in% keep_samples,
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
    levels = condition_levels
  )

  if (
    !all(c("EP", "Veh") %in% levels(cond))
  ) {
    stop("EP/Veh missing in ", label)
  }

  dge <- DGEList(
    counts = counts_x,
    group = cond
  )

  keep_gene <- filterByExpr(
    dge,
    group = cond,
    min.count = 10
  )

  dge <- dge[
    keep_gene,
    ,
    keep.lib.sizes = FALSE
  ]

  dge <- calcNormFactors(dge)

  design <- model.matrix(
    ~ 0 + cond
  )

  colnames(design) <- levels(cond)

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

  fit2 <- eBayes(fit2)

  tt <- topTable(
    fit2,
    coef = "EP_vs_Veh",
    number = Inf,
    adjust.method = "BH",
    sort.by = "none"
  )

  sm <- tt[
    tolower(rownames(tt)) == "smad3",
    ,
    drop = FALSE
  ]

  if (nrow(sm) != 1L) {
    stop(
      "Expected exactly one Smad3 row in ",
      label
    )
  }

  data.frame(
    model = label,
    n_EP = sum(cond == "EP"),
    n_Veh = sum(cond == "Veh"),
    n_EPR = sum(cond == "EPR"),
    n_genes_tested = nrow(tt),
    Smad3_log2FC_EP_vs_Veh = sm$logFC,
    Smad3_P = sm$P.Value,
    Smad3_BH_FDR = sm$adj.P.Val,
    stringsAsFactors = FALSE
  )
}

all_samples <- pb_meta$sample

epveh_samples <- pb_meta$sample[
  pb_meta$condition %in% c(
    "EP",
    "Veh"
  )
]

# Historical parser removes the two underscore-named samples
# from the correctly named Type IIa group.
historical_true_type_samples <- setdiff(
  all_samples,
  bad_samples
)

epveh_without_1696 <- setdiff(
  epveh_samples,
  "1696_EP"
)

res <- rbind(
  run_model(
    "A_correct_3group",
    all_samples,
    c("Veh", "EP", "EPR")
  ),

  run_model(
    "B_correct_EP_vs_Veh",
    epveh_samples,
    c("Veh", "EP")
  ),

  run_model(
    "C_historical_parser_3group",
    historical_true_type_samples,
    c("Veh", "EP", "EPR")
  ),

  run_model(
    "D_EP_vs_Veh_without_1696_EP",
    epveh_without_1696,
    c("Veh", "EP")
  )
)

cat("\n===== CONTROLLED MODEL COMPARISON =====\n")

print(
  res,
  row.names = FALSE,
  digits = 10
)

# --------------------------------------------------------------------------
# Validate model B against 16B rebuild
# --------------------------------------------------------------------------

rebuild_file <- file.path(
  "processed_results",
  "13_mouse_validation",
  "pseudobulk_rebuild_audit",
  "mouse_EP_vs_Veh_pseudobulk_fullprecision.csv"
)

if (file.exists(rebuild_file)) {

  x <- read.csv(
    rebuild_file,
    check.names = FALSE
  )

  x <- x[
    x$celltype == "Type IIa Myofiber" &
      tolower(x$gene) == "smad3",
    ,
    drop = FALSE
  ]

  if (nrow(x) != 1L) {
    stop(
      "Could not identify unique 16B Type IIa Smad3 row."
    )
  }

  b <- res[
    res$model == "B_correct_EP_vs_Veh",
    ,
    drop = FALSE
  ]

  delta_fc <- abs(
    b$Smad3_log2FC_EP_vs_Veh -
      x$avg_log2FC
  )

  delta_p <- abs(
    b$Smad3_P -
      x$p_val
  )

  delta_fdr <- abs(
    b$Smad3_BH_FDR -
      x$p_val_adj
  )

  cat("\n===== 16B REPRODUCTION =====\n")
  cat("delta logFC:", delta_fc, "\n")
  cat("delta P    :", delta_p, "\n")
  cat("delta FDR  :", delta_fdr, "\n")

  if (
    delta_fc > 1e-10 ||
    delta_p > 1e-10 ||
    delta_fdr > 1e-10
  ) {
    stop(
      "Type IIa targeted model does not reproduce 16B."
    )
  }

  cat("PASS: model B exactly reproduces 16B.\n")
}

# --------------------------------------------------------------------------
# Historical canonical comparison
# --------------------------------------------------------------------------

old_file <- file.path(
  "processed_results",
  "13_mouse_validation",
  "bulk_de_all_contrasts_mouse.csv"
)

if (file.exists(old_file)) {

  old <- read.csv(
    old_file,
    check.names = FALSE
  )

  old_sm <- old[
    old$celltype == "Type IIa Myofiber" &
      old$contrast == "Veh_vs_EP" &
      tolower(old$gene) == "smad3",
    ,
    drop = FALSE
  ]

  cat("\n===== HISTORICAL CANONICAL =====\n")

  if (nrow(old_sm) == 1L) {

    cat(
      "Converted to EP_vs_Veh:\n"
    )

    cat(
      "logFC:",
      -old_sm$avg_log2FC,
      "\n"
    )

    cat(
      "P:",
      old_sm$p_val,
      "\n"
    )

    cat(
      "FDR:",
      old_sm$p_val_adj,
      "\n"
    )

  } else {

    cat(
      "Historical Type IIa Smad3 row not uniquely identified.\n"
    )
  }
}

write.csv(
  res,
  file.path(
    "processed_results",
    "13_mouse_validation",
    "pseudobulk_rebuild_audit",
    "Smad3_TypeIIa_model_sensitivity.csv"
  ),
  row.names = FALSE
)

cat(
  "\nPASS: Type IIa Smad3 model audit completed.\n"
)
