rm(list = ls(all.names = TRUE))
gc()

suppressPackageStartupMessages({
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

audit_file <- file.path(
  repo_root,
  "processed_results",
  "12_CellChat",
  "multiple_testing_audit",
  "TGFb_lightweight",
  "TGFb_nboot1000_all_tests_combined.tsv"
)

canonical_source_file <- file.path(
  repo_root,
  "source_data",
  "figure_source_data",
  "Figure_18",
  "Figure18D_TGFb_pathway_network_source_data.csv"
)

candidate_dir <- file.path(
  repo_root,
  "processed_results",
  "12_CellChat",
  "multiple_testing_audit",
  "figure18_revised_candidate"
)

dir.create(
  candidate_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(audit_file)) {
  stop("TGFb nboot1000 audit table not found.")
}

if (!file.exists(canonical_source_file)) {
  stop("Canonical Figure18D source data not found.")
}

dt <- fread(audit_file)

required <- c(
  "dataset",
  "source",
  "target",
  "interaction_name",
  "probability",
  "p_raw",
  "p_plus1",
  "p_plus1_BH_TGFb"
)

missing_cols <- setdiff(
  required,
  names(dt)
)

if (length(missing_cols) > 0L) {
  stop(
    "Missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

# --------------------------------------------------------------------------
# Corrected TGFb pathway aggregation:
# retain LR-level communications passing plus-one empirical P + BH.
# --------------------------------------------------------------------------
corrected <- dt[
  p_plus1_BH_TGFb < 0.05,
  .(
    TGFb_pathway_probability =
      sum(
        probability,
        na.rm = TRUE
      ),

    n_significant_LR_links =
      .N,

    minimum_plus1_BH =
      min(
        p_plus1_BH_TGFb,
        na.rm = TRUE
      )
  ),
  by = .(
    condition = dataset,
    source,
    target
  )
]

# Preserve all source-target combinations, including zeros, as in Figure18D.
groups <- sort(
  unique(
    c(
      dt$source,
      dt$target
    )
  )
)

grid <- CJ(
  condition = c("Young", "Aged"),
  source = groups,
  target = groups,
  unique = TRUE
)

corrected_full <- merge(
  grid,
  corrected,
  by = c(
    "condition",
    "source",
    "target"
  ),
  all.x = TRUE,
  sort = FALSE
)

corrected_full[
  is.na(TGFb_pathway_probability),
  `:=`(
    TGFb_pathway_probability = 0,
    n_significant_LR_links = 0L,
    minimum_plus1_BH = NA_real_
  )
]

# Match canonical source-data ordering as closely as possible.
canonical <- fread(
  canonical_source_file
)

cat("Canonical columns:\n")
print(names(canonical))

# Detect canonical source/target field names.
source_col <- intersect(
  c("source", "source_group", "source_name"),
  names(canonical)
)

target_col <- intersect(
  c("target", "target_group", "target_name"),
  names(canonical)
)

condition_col <- intersect(
  c("condition", "dataset"),
  names(canonical)
)

value_col <- intersect(
  c(
    "TGFb_pathway_probability",
    "value",
    "probability"
  ),
  names(canonical)
)

if (
  length(source_col) != 1L ||
  length(target_col) != 1L ||
  length(condition_col) != 1L ||
  length(value_col) != 1L
) {
  stop(
    "Could not uniquely identify canonical Figure18D columns."
  )
}

# Compact corrected table matching canonical core schema.
candidate_core <- corrected_full[
  ,
  .(
    condition,
    source,
    target,
    TGFb_pathway_probability
  )
]

# --------------------------------------------------------------------------
# Compare against canonical source-data.
# --------------------------------------------------------------------------
canonical_cmp <- canonical[
  ,
  .(
    condition =
      as.character(get(condition_col)),

    source =
      as.character(get(source_col)),

    target =
      as.character(get(target_col)),

    canonical_probability =
      as.numeric(get(value_col))
  )
]

candidate_cmp <- candidate_core[
  ,
  .(
    condition,
    source,
    target,
    corrected_probability =
      TGFb_pathway_probability
  )
]

cmp <- merge(
  canonical_cmp,
  candidate_cmp,
  by = c(
    "condition",
    "source",
    "target"
  ),
  all = TRUE,
  sort = FALSE
)

cmp[
  ,
  delta :=
    corrected_probability -
      canonical_probability
]

cmp[
  ,
  abs_delta :=
    abs(delta)
]

changed <- cmp[
  is.na(canonical_probability) |
    is.na(corrected_probability) |
    abs_delta > 1e-12
]

setorder(
  changed,
  condition,
  -abs_delta
)

# --------------------------------------------------------------------------
# Save candidate-only outputs.
# --------------------------------------------------------------------------
fwrite(
  candidate_core,
  file.path(
    candidate_dir,
    "Figure18D_TGFb_pathway_network_corrected_candidate.csv"
  )
)

fwrite(
  corrected_full,
  file.path(
    candidate_dir,
    "Figure18D_TGFb_pathway_network_corrected_with_audit.csv"
  )
)

fwrite(
  cmp,
  file.path(
    candidate_dir,
    "Figure18D_canonical_vs_corrected_comparison.csv"
  )
)

fwrite(
  changed,
  file.path(
    candidate_dir,
    "Figure18D_changed_edges.csv"
  )
)

# --------------------------------------------------------------------------
# Console summary
# --------------------------------------------------------------------------
cat("\n============================================================\n")
cat("FIGURE 18D CORRECTED SOURCE-DATA AUDIT\n")
cat("============================================================\n\n")

cat(
  "Canonical rows :",
  nrow(canonical_cmp),
  "\n"
)

cat(
  "Candidate rows :",
  nrow(candidate_cmp),
  "\n"
)

cat(
  "Changed rows (>1e-12):",
  nrow(changed),
  "\n\n"
)

if (nrow(changed) > 0L) {
  print(changed)
}

cat("\nExpected affected edge:\n")

print(
  candidate_core[
    condition == "Young" &
      source == "Vascular stromal" &
      target == "FAP3"
  ]
)

cat("\nCandidate files written to:\n")
cat(candidate_dir, "\n")

cat("\nNo canonical source-data or figure file was modified.\n")
