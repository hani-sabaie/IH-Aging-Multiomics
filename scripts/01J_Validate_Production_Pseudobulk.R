# ============================================================================
# Validate the revised production pseudobulk logic before canonicalization.
#
# This independently reproduces the pseudobulk block now implemented in 01A
# and compares it against the previously validated 01I revised-candidate
# results.
#
# No canonical processed result is modified.
# ============================================================================

rm(list = ls(all.names = TRUE))
gc()

suppressPackageStartupMessages({
  library(Seurat)
  library(edgeR)
  library(limma)
  library(dplyr)
  library(data.table)
})

# --------------------------------------------------------------------------
# Repository root
# --------------------------------------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (length(file_arg) == 1L) {
  script_dir <- dirname(
    normalizePath(
      sub("^--file=", "", file_arg)
    )
  )
  repo_root <- normalizePath(
    file.path(script_dir, "..")
  )
} else {
  repo_root <- normalizePath(".")
}

# --------------------------------------------------------------------------
# Input object
# --------------------------------------------------------------------------
obj_candidates <- c(
  file.path(
    repo_root,
    "outputs",
    "decont_merged_filt_nodoub_cc_sct_reduc_clust_integ_annot_obj.rds"
  ),
  "F:/Hani's Files/Hernia/outputs/decont_merged_filt_nodoub_cc_sct_reduc_clust_integ_annot_obj.rds"
)

obj_hits <- obj_candidates[
  file.exists(obj_candidates)
]

if (length(obj_hits) == 0L) {
  stop(
    "Annotated Seurat object not found. Checked:\n  ",
    paste(obj_candidates, collapse = "\n  ")
  )
}

obj_file <- normalizePath(
  obj_hits[1],
  winslash = "/",
  mustWork = TRUE
)

candidate_dir <- file.path(
  repo_root,
  "processed_results",
  "02_differential_expression",
  "revised_candidate"
)

candidate_all_file <- file.path(
  candidate_dir,
  "all_gene_results_full_precision.tsv"
)

candidate_deg_file <- file.path(
  candidate_dir,
  "DEGs_BH_FDR005.csv"
)

candidate_fap_file <- file.path(
  candidate_dir,
  "DEGs_FAP_BH_FDR005.csv"
)

for (f in c(
  candidate_all_file,
  candidate_deg_file,
  candidate_fap_file
)) {
  if (!file.exists(f)) {
    stop("Required 01I candidate file missing: ", f)
  }
}

outdir <- file.path(
  repo_root,
  "processed_results",
  "02_differential_expression",
  "production_validation"
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("Annotated object:\n", obj_file, "\n\n", sep = "")

# --------------------------------------------------------------------------
# Load object
# --------------------------------------------------------------------------
cat("Loading annotated Seurat object...\n")

obj <- readRDS(obj_file)

cat("Cells:", ncol(obj), "\n")
cat("RNA genes:", nrow(obj[["RNA"]]), "\n\n")

DefaultAssay(obj) <- "RNA"

# --------------------------------------------------------------------------
# 1. Production eligibility rule
# --------------------------------------------------------------------------
donor_celltype_nuclei <- obj@meta.data %>%
  transmute(
    sample = as.character(sample),
    condition = as.character(condition),
    celltype = as.character(skeletal_muscle)
  ) %>%
  count(
    sample,
    condition,
    celltype,
    name = "n_nuclei"
  ) %>%
  arrange(
    celltype,
    condition,
    sample
  )

adequacy <- donor_celltype_nuclei %>%
  mutate(
    donor_ge5_nuclei = n_nuclei >= 5
  ) %>%
  group_by(
    celltype,
    condition
  ) %>%
  summarise(
    n_donors_present = n_distinct(sample),
    n_donors_ge5_nuclei = sum(donor_ge5_nuclei),
    total_nuclei = sum(n_nuclei),
    .groups = "drop"
  )

eligibility <- adequacy %>%
  group_by(celltype) %>%
  summarise(
    n_Young_donors_ge5 =
      sum(
        n_donors_ge5_nuclei[
          condition == "Young"
        ],
        na.rm = TRUE
      ),

    n_Aged_donors_ge5 =
      sum(
        n_donors_ge5_nuclei[
          condition == "Aged"
        ],
        na.rm = TRUE
      ),

    eligible_for_formal_DE =
      n_Young_donors_ge5 >= 3 &
      n_Aged_donors_ge5 >= 3,

    .groups = "drop"
  ) %>%
  arrange(celltype)

eligible_celltypes <- eligibility %>%
  filter(eligible_for_formal_DE) %>%
  pull(celltype)

excluded_celltypes <- eligibility %>%
  filter(!eligible_for_formal_DE) %>%
  pull(celltype)

cat("============================================================\n")
cat("CELL-TYPE ELIGIBILITY\n")
cat("============================================================\n\n")

print(eligibility)

cat(
  "\nEligible:",
  paste(eligible_celltypes, collapse = ", "),
  "\n"
)

cat(
  "Excluded:",
  paste(excluded_celltypes, collapse = ", "),
  "\n\n"
)

# --------------------------------------------------------------------------
# 2. Aggregate raw RNA counts
# --------------------------------------------------------------------------
cat("Aggregating pseudobulk RNA counts...\n")

pb_list <- AggregateExpression(
  obj,
  assays = "RNA",
  group.by = c(
    "sample",
    "condition",
    "skeletal_muscle"
  ),
  slot = "counts",
  return.seurat = FALSE
)

pb_counts <- pb_list$RNA

cat(
  "Pseudobulk matrix:",
  nrow(pb_counts),
  "genes x",
  ncol(pb_counts),
  "groups\n\n"
)

# Object no longer needed after aggregation/metadata extraction.
rm(obj)
gc()

# --------------------------------------------------------------------------
# 3. Parse pseudobulk column metadata exactly as revised 01A
# --------------------------------------------------------------------------
parse_pb_column <- function(x) {

  m <- regexec(
    "^(.*)_([^_]+)_([^_]+)$",
    x
  )

  z <- regmatches(
    x,
    m
  )[[1]]

  if (length(z) != 4L) {
    stop(
      "Unable to parse pseudobulk column name: ",
      x
    )
  }

  data.frame(
    sample = z[2],
    condition = z[3],
    celltype = z[4],
    stringsAsFactors = FALSE
  )
}

meta <- bind_rows(
  lapply(
    colnames(pb_counts),
    parse_pb_column
  )
)

rownames(meta) <- colnames(pb_counts)

meta$condition <- factor(
  meta$condition,
  levels = c(
    "Young",
    "Aged"
  )
)

if (anyNA(meta$condition)) {
  stop("Unexpected condition in pseudobulk metadata.")
}

# --------------------------------------------------------------------------
# 4. Production voom/limma implementation
# --------------------------------------------------------------------------
run_voomlimma <- function(
    counts,
    condition
) {

  condition <- droplevels(condition)

  condition <- relevel(
    condition,
    ref = "Young"
  )

  dge <- DGEList(
    counts = counts,
    group = condition
  )

  keep_genes <- filterByExpr(
    dge,
    group = condition,
    min.count = 10
  )

  dge <- dge[
    keep_genes,
    ,
    keep.lib.sizes = FALSE
  ]

  dge <- calcNormFactors(dge)

  design <- model.matrix(
    ~ condition
  )

  vm <- voom(
    dge,
    design = design,
    plot = FALSE
  )

  fit <- lmFit(
    vm,
    design = design
  )

  fit <- eBayes(fit)

  tt <- topTable(
    fit,
    coef = "conditionAged",
    n = Inf,
    adjust.method = "BH",
    sort.by = "none"
  )

  tt
}

res_list <- lapply(
  sort(eligible_celltypes),
  function(st) {

    cat("Running:", st, "\n")

    keep_cols <- meta$celltype == st

    counts_ct <- pb_counts[
      ,
      keep_cols,
      drop = FALSE
    ]

    cond_ct <- droplevels(
      meta$condition[
        keep_cols
      ]
    )

    if (
      length(unique(cond_ct)) < 2L ||
      any(table(cond_ct) < 2L)
    ) {
      stop(
        "Unexpected model-level replicate failure for ",
        st
      )
    }

    tt <- run_voomlimma(
      counts = counts_ct,
      condition = cond_ct
    )

    data.frame(
      gene = rownames(tt),
      p_val = tt$P.Value,
      p_val_adj = tt$adj.P.Val,
      avg_log2FC = tt$logFC,
      celltype = st,
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }
)

production_all <- bind_rows(
  res_list
) %>%
  arrange(
    celltype,
    p_val_adj,
    p_val,
    gene
  )

production_deg <- production_all %>%
  filter(
    p_val_adj < 0.05
  ) %>%
  arrange(
    celltype,
    p_val_adj,
    p_val,
    gene
  )

production_fap <- production_deg %>%
  filter(
    celltype %in% c(
      "FAP1",
      "FAP2",
      "FAP3"
    )
  )

# --------------------------------------------------------------------------
# 5. Compare with validated 01I all-gene result
# --------------------------------------------------------------------------
candidate_all <- fread(
  candidate_all_file
)

required_cols <- c(
  "gene",
  "celltype",
  "p_val",
  "p_val_adj",
  "avg_log2FC"
)

missing_candidate <- setdiff(
  required_cols,
  names(candidate_all)
)

if (length(missing_candidate) > 0L) {
  stop(
    "01I all-gene table missing columns: ",
    paste(missing_candidate, collapse = ", ")
  )
}

candidate_core <- candidate_all[
  ,
  ..required_cols
]

production_dt <- as.data.table(
  production_all
)

setkey(
  production_dt,
  gene,
  celltype
)

setkey(
  candidate_core,
  gene,
  celltype
)

shared <- merge(
  production_dt,
  candidate_core,
  by = c(
    "gene",
    "celltype"
  ),
  suffixes = c(
    "_production",
    "_01I"
  )
)

production_only <- fsetdiff(
  production_dt[
    ,
    .(
      gene,
      celltype
    )
  ],
  candidate_core[
    ,
    .(
      gene,
      celltype
    )
  ]
)

candidate_only <- fsetdiff(
  candidate_core[
    ,
    .(
      gene,
      celltype
    )
  ],
  production_dt[
    ,
    .(
      gene,
      celltype
    )
  ]
)

shared[
  ,
  delta_logFC :=
    abs(
      avg_log2FC_production -
        avg_log2FC_01I
    )
]

shared[
  ,
  delta_P :=
    abs(
      p_val_production -
        p_val_01I
    )
]

shared[
  ,
  delta_FDR :=
    abs(
      p_val_adj_production -
        p_val_adj_01I
    )
]

# --------------------------------------------------------------------------
# 6. Compare significant DEG key sets
# --------------------------------------------------------------------------
candidate_deg <- fread(
  candidate_deg_file
)

candidate_fap <- fread(
  candidate_fap_file
)

key_cols <- c(
  "gene",
  "celltype"
)

for (nm in key_cols) {
  if (!(nm %in% names(candidate_deg))) {
    stop("01I DEG file missing column: ", nm)
  }

  if (!(nm %in% names(candidate_fap))) {
    stop("01I FAP DEG file missing column: ", nm)
  }
}

production_deg_keys <- unique(
  as.data.table(production_deg)[
    ,
    ..key_cols
  ]
)

candidate_deg_keys <- unique(
  candidate_deg[
    ,
    ..key_cols
  ]
)

production_fap_keys <- unique(
  as.data.table(production_fap)[
    ,
    ..key_cols
  ]
)

candidate_fap_keys <- unique(
  candidate_fap[
    ,
    ..key_cols
  ]
)

deg_production_only <- fsetdiff(
  production_deg_keys,
  candidate_deg_keys
)

deg_candidate_only <- fsetdiff(
  candidate_deg_keys,
  production_deg_keys
)

fap_production_only <- fsetdiff(
  production_fap_keys,
  candidate_fap_keys
)

fap_candidate_only <- fsetdiff(
  candidate_fap_keys,
  production_fap_keys
)

# --------------------------------------------------------------------------
# 7. Critical candidates
# --------------------------------------------------------------------------
critical <- production_dt[
  celltype == "FAP3" &
    gene %in% c(
      "SMAD3",
      "PLEKHA6"
    )
][
  order(gene)
]

# --------------------------------------------------------------------------
# 8. Summary
# --------------------------------------------------------------------------
summary_dt <- data.table(
  metric = c(
    "production_all_tests",
    "candidate_01I_all_tests",
    "shared_tests",
    "production_only_tests",
    "candidate_only_tests",
    "max_abs_delta_log2FC",
    "max_abs_delta_P",
    "max_abs_delta_FDR",
    "production_DEG_rows",
    "candidate_01I_DEG_rows",
    "production_only_DEG_keys",
    "candidate_only_DEG_keys",
    "production_FAP_rows",
    "candidate_01I_FAP_rows",
    "production_FAP_unique_genes",
    "production_only_FAP_keys",
    "candidate_only_FAP_keys"
  ),

  value = c(
    nrow(production_dt),
    nrow(candidate_core),
    nrow(shared),
    nrow(production_only),
    nrow(candidate_only),

    max(
      shared$delta_logFC,
      na.rm = TRUE
    ),

    max(
      shared$delta_P,
      na.rm = TRUE
    ),

    max(
      shared$delta_FDR,
      na.rm = TRUE
    ),

    nrow(production_deg),
    nrow(candidate_deg),

    nrow(deg_production_only),
    nrow(deg_candidate_only),

    nrow(production_fap),
    nrow(candidate_fap),

    n_distinct(
      production_fap$gene
    ),

    nrow(fap_production_only),
    nrow(fap_candidate_only)
  )
)

# --------------------------------------------------------------------------
# 9. Expected-result safeguards
# --------------------------------------------------------------------------
validation_pass <-
  nrow(production_dt) == 106597L &&
  nrow(candidate_core) == 106597L &&
  nrow(shared) == 106597L &&
  nrow(production_only) == 0L &&
  nrow(candidate_only) == 0L &&
  max(shared$delta_logFC, na.rm = TRUE) < 1e-12 &&
  max(shared$delta_P, na.rm = TRUE) < 1e-12 &&
  max(shared$delta_FDR, na.rm = TRUE) < 1e-12 &&
  nrow(production_deg) == 1154L &&
  nrow(production_fap) == 210L &&
  n_distinct(production_fap$gene) == 194L &&
  nrow(deg_production_only) == 0L &&
  nrow(deg_candidate_only) == 0L &&
  nrow(fap_production_only) == 0L &&
  nrow(fap_candidate_only) == 0L &&
  setequal(
    excluded_celltypes,
    c(
      "EC3",
      "FAP4"
    )
  )

# --------------------------------------------------------------------------
# 10. Export validation-only resources
# --------------------------------------------------------------------------
fwrite(
  production_dt,
  file.path(
    outdir,
    "production_logic_all_gene_results.tsv"
  ),
  sep = "\t"
)

fwrite(
  summary_dt,
  file.path(
    outdir,
    "production_vs_01I_validation_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  eligibility,
  file.path(
    outdir,
    "production_celltype_eligibility.tsv"
  ),
  sep = "\t"
)

fwrite(
  critical,
  file.path(
    outdir,
    "production_SMAD3_PLEKHA6.tsv"
  ),
  sep = "\t"
)

# --------------------------------------------------------------------------
# Console output
# --------------------------------------------------------------------------
cat("\n============================================================\n")
cat("PRODUCTION VS 01I VALIDATION\n")
cat("============================================================\n\n")

print(summary_dt)

cat("\n============================================================\n")
cat("CRITICAL FAP3 CANDIDATES\n")
cat("============================================================\n\n")

print(critical)

cat("\n============================================================\n")
cat("FINAL STATUS\n")
cat("============================================================\n\n")

cat(
  "Production validation pass:",
  validation_pass,
  "\n"
)

cat("\nValidation outputs:\n")
cat(outdir, "\n")

cat("\nNo canonical DEG or supplementary-table file was modified.\n")

if (!validation_pass) {
  stop(
    "Production pseudobulk validation failed. ",
    "Do not canonicalize outputs."
  )
}
