# ============================================================================
# Reviewer C7
# Promote the fully validated corrected CellChat object and Figure 18
# source-data resources to canonical repository locations.
#
# No inference is performed here.
# ============================================================================

rm(list = ls(all.names = TRUE))

suppressPackageStartupMessages({
  library(data.table)
})

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
# Validated candidate inputs
# --------------------------------------------------------------------------

audit_root <- file.path(
  repo_root,
  "processed_results",
  "12_CellChat",
  "multiple_testing_audit",
  "full18_nboot1000_BH"
)

candidate_root <- file.path(
  audit_root,
  "figure18ABC_corrected_candidate"
)

candidate_obj <- file.path(
  candidate_root,
  "cellchat_merged_full18_nboot1000_plus1BH_corrected_candidate.rds"
)

candidate_source_dir <- file.path(
  audit_root,
  "figure18_full_corrected_source_candidate"
)

# --------------------------------------------------------------------------
# Canonical outputs
# --------------------------------------------------------------------------

processed_dir <- file.path(
  repo_root,
  "processed_results",
  "12_CellChat"
)

canonical_obj <- file.path(
  processed_dir,
  "cellchat_merged.rds"
)

canonical_source_dir <- file.path(
  repo_root,
  "source_data",
  "figure_source_data",
  "Figure_18"
)

canonical_mt_dir <- file.path(
  processed_dir,
  "multiple_testing"
)

dir.create(
  canonical_source_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  canonical_mt_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

source_files <- c(
  "Figure18A_overall_signaling_network_source_data.csv",
  "Figure18B_signaling_role_scatter_source_data.csv",
  "Figure18C_signaling_role_heatmap_source_data.csv",
  "Figure18D_TGFb_pathway_network_source_data.csv",
  "Figure18E_TGFb_LR_contribution_Aged_source_data.csv",
  "Figure18F_TGFb_bubble_source_data.csv"
)

source_expected_rows <- c(
  98L,
  14L,
  308L,
  98L,
  4L,
  4L
)

datasets <- c(
  "Young",
  "Aged"
)

expected_lr <- c(
  Young = 488L,
  Aged = 468L
)

expected_tests <- c(
  Young = 23912L,
  Aged = 22932L
)

# These are validated outcomes, used only to safeguard promotion.
# They are NOT analysis-selection criteria in the production inference code.
expected_sig <- c(
  Young = 680L,
  Aged = 702L
)

expected_retained_lr <- c(
  Young = 104L,
  Aged = 111L
)

# --------------------------------------------------------------------------
# Input existence
# --------------------------------------------------------------------------

required <- c(
  candidate_obj,
  file.path(
    candidate_source_dir,
    source_files
  )
)

for (ds in datasets) {

  required <- c(
    required,

    file.path(
      audit_root,
      paste0(
        ds,
        "_full18_nboot1000_conditionwide_BH_tests.tsv"
      )
    ),

    file.path(
      audit_root,
      paste0(
        ds,
        "_full18_nboot1000_conditionwide_BH_summary.tsv"
      )
    )
  )
}

missing <- required[
  !file.exists(required)
]

if (length(missing) > 0L) {
  stop(
    "Required validated input(s) missing:\n",
    paste(
      missing,
      collapse = "\n"
    )
  )
}

# --------------------------------------------------------------------------
# Validate candidate Figure18 source files
# --------------------------------------------------------------------------

cat("\nValidating Figure 18 source-data candidates...\n")

for (i in seq_along(source_files)) {

  f <- file.path(
    candidate_source_dir,
    source_files[i]
  )

  z <- fread(f)

  if (
    nrow(z) !=
      source_expected_rows[i]
  ) {
    stop(
      source_files[i],
      ": expected ",
      source_expected_rows[i],
      " rows; found ",
      nrow(z),
      "."
    )
  }

  cat(
    "  ",
    source_files[i],
    ": ",
    nrow(z),
    " rows\n",
    sep = ""
  )
}

# Figure18F revised inferential semantics must be explicit.
fig18f <- fread(
  file.path(
    candidate_source_dir,
    "Figure18F_TGFb_bubble_source_data.csv"
  )
)

required_f_cols <- c(
  "CellChat_probability",
  "plus1_BH_adjusted_p_value",
  "significant_plus1_BH_p_lt_0_05"
)

missing_f_cols <- setdiff(
  required_f_cols,
  names(fig18f)
)

if (length(missing_f_cols) > 0L) {
  stop(
    "Corrected Figure18F source data is missing revised columns: ",
    paste(
      missing_f_cols,
      collapse = ", "
    )
  )
}

# --------------------------------------------------------------------------
# Validate corrected CellChat object
# --------------------------------------------------------------------------

cat("\nValidating corrected CellChat object...\n")

x <- readRDS(
  candidate_obj
)

if (
  !all(
    datasets %in%
      names(x@net)
  ) ||
  !all(
    datasets %in%
      names(x@netP)
  )
) {
  stop(
    "Corrected merged CellChat object is missing Young/Aged networks."
  )
}

production_summaries <- list()

for (ds in datasets) {

  prob <- x@net[[ds]]$prob
  pval <- x@net[[ds]]$pval

  if (
    length(dim(prob)) != 3L ||
    dim(prob)[1] != 7L ||
    dim(prob)[2] != 7L ||
    dim(prob)[3] != expected_lr[[ds]]
  ) {
    stop(
      ds,
      ": unexpected corrected probability-array dimensions."
    )
  }

  if (
    !identical(
      dim(prob),
      dim(pval)
    )
  ) {
    stop(
      ds,
      ": probability/P arrays have different dimensions."
    )
  }

  if (
    anyNA(pval) ||
    any(
      pval < 0 |
        pval > 1
    )
  ) {
    stop(
      ds,
      ": invalid corrected network P values."
    )
  }

  n_sig <- sum(
    pval < 0.05
  )

  if (
    n_sig !=
      expected_sig[[ds]]
  ) {
    stop(
      ds,
      ": expected ",
      expected_sig[[ds]],
      " validated significant tests; found ",
      n_sig,
      "."
    )
  }

  if (
    sum(
      x@net[[ds]]$count,
      na.rm = TRUE
    ) !=
      expected_sig[[ds]]
  ) {
    stop(
      ds,
      ": aggregate interaction count does not match validated result."
    )
  }

  retained_lr <- length(
    x@net[[ds]]$LRs
  )

  if (
    retained_lr !=
      expected_retained_lr[[ds]]
  ) {
    stop(
      ds,
      ": expected ",
      expected_retained_lr[[ds]],
      " retained LR interactions; found ",
      retained_lr,
      "."
    )
  }

  if (
    length(
      x@netP[[ds]]$pathways
    ) != 11L
  ) {
    stop(
      ds,
      ": expected 11 retained pathways."
    )
  }

  # ------------------------------------------------------------------------
  # Validate full multiple-testing table against corrected object.
  # ------------------------------------------------------------------------

  tests_file <- file.path(
    audit_root,
    paste0(
      ds,
      "_full18_nboot1000_conditionwide_BH_tests.tsv"
    )
  )

  tests <- fread(
    tests_file
  )

  if (
    nrow(tests) !=
      expected_tests[[ds]]
  ) {
    stop(
      ds,
      ": expected ",
      expected_tests[[ds]],
      " tests; found ",
      nrow(tests),
      "."
    )
  }

  required_test_cols <- c(
    "LR_family_index",
    "interaction_name",
    "source",
    "target",
    "probability",
    "p_raw_nboot1000",
    "permutation_exceedance_count",
    "p_plus1",
    "p_raw_BH_condition",
    "p_plus1_BH_condition"
  )

  missing_test_cols <- setdiff(
    required_test_cols,
    names(tests)
  )

  if (length(missing_test_cols) > 0L) {
    stop(
      ds,
      ": multiple-testing table missing columns: ",
      paste(
        missing_test_cols,
        collapse = ", "
      )
    )
  }

  key_cols <- c(
    "LR_family_index",
    "interaction_name",
    "source",
    "target"
  )

  if (
    anyDuplicated(
      tests[
        ,
        ..key_cols
      ]
    ) > 0L
  ) {
    stop(
      ds,
      ": duplicate formal hypotheses detected."
    )
  }

  if (
    uniqueN(
      tests$LR_family_index
    ) !=
      expected_lr[[ds]]
  ) {
    stop(
      ds,
      ": incomplete LR family in multiple-testing table."
    )
  }

  if (
    sum(
      tests$p_plus1_BH_condition <
        0.05
    ) !=
      expected_sig[[ds]]
  ) {
    stop(
      ds,
      ": multiple-testing table does not reproduce validated significant count."
    )
  }

  groups <- dimnames(
    prob
  )[[1]]

  lr_names <- dimnames(
    prob
  )[[3]]

  src_i <- match(
    tests$source,
    groups
  )

  tgt_i <- match(
    tests$target,
    groups
  )

  lr_i <- match(
    tests$interaction_name,
    lr_names
  )

  if (
    anyNA(src_i) ||
    anyNA(tgt_i) ||
    anyNA(lr_i)
  ) {
    stop(
      ds,
      ": failed to map formal hypotheses to corrected object."
    )
  }

  obj_prob <- prob[
    cbind(
      src_i,
      tgt_i,
      lr_i
    )
  ]

  obj_p <- pval[
    cbind(
      src_i,
      tgt_i,
      lr_i
    )
  ]

  max_prob_delta <- max(
    abs(
      obj_prob -
        tests$probability
    )
  )

  max_p_delta <- max(
    abs(
      obj_p -
        tests$p_plus1_BH_condition
    )
  )

  if (
    max_prob_delta > 1e-12 ||
    max_p_delta > 1e-12
  ) {
    stop(
      ds,
      ": corrected object and multiple-testing table disagree."
    )
  }

  audit_summary <- fread(
    file.path(
      audit_root,
      paste0(
        ds,
        "_full18_nboot1000_conditionwide_BH_summary.tsv"
      )
    )
  )

  if (nrow(audit_summary) != 1L) {
    stop(
      ds,
      ": unexpected audit summary structure."
    )
  }

  production_summary <- data.table(
    dataset = ds,
    n_batches =
      audit_summary$n_batches,
    n_LR =
      expected_lr[[ds]],
    n_source_groups = 7L,
    n_target_groups = 7L,
    n_tests =
      nrow(tests),

    n_raw_p_zero =
      sum(
        tests$p_raw_nboot1000 ==
          0
      ),

    n_raw_p_lt_0_05 =
      sum(
        tests$p_raw_nboot1000 <
          0.05
      ),

    n_raw_BH_lt_0_05 =
      sum(
        tests$p_raw_BH_condition <
          0.05
      ),

    n_plus1_BH_lt_0_05 =
      sum(
        tests$p_plus1_BH_condition <
          0.05
      ),

    n_retained_LR =
      retained_lr,

    n_retained_pathways =
      length(
        x@netP[[ds]]$pathways
      ),

    min_raw_p =
      min(
        tests$p_raw_nboot1000
      ),

    min_plus1_p =
      min(
        tests$p_plus1
      ),

    min_plus1_BH =
      min(
        tests$p_plus1_BH_condition
      ),

    total_batch_elapsed_seconds =
      audit_summary$total_batch_elapsed_seconds
  )

  production_summaries[[ds]] <-
    production_summary

  cat(
    "  ",
    ds,
    ": ",
    nrow(tests),
    " tests; ",
    n_sig,
    " significant; ",
    retained_lr,
    " retained LR; ",
    length(
      x@netP[[ds]]$pathways
    ),
    " pathways\n",
    sep = ""
  )

  cat(
    "    object/table max probability delta = ",
    format(
      max_prob_delta,
      scientific = TRUE,
      digits = 12
    ),
    "\n",
    sep = ""
  )

  cat(
    "    object/table max corrected-P delta = ",
    format(
      max_p_delta,
      scientific = TRUE,
      digits = 12
    ),
    "\n",
    sep = ""
  )
}

# --------------------------------------------------------------------------
# Promotion
# --------------------------------------------------------------------------

cat("\nPromoting corrected canonical resources...\n")

# Binary corrected CellChat object
if (
  !file.copy(
    candidate_obj,
    canonical_obj,
    overwrite = TRUE
  )
) {
  stop(
    "Failed to promote corrected CellChat object."
  )
}

# Figure18 source data
for (f in source_files) {

  src <- file.path(
    candidate_source_dir,
    f
  )

  dst <- file.path(
    canonical_source_dir,
    f
  )

  if (
    !file.copy(
      src,
      dst,
      overwrite = TRUE
    )
  ) {
    stop(
      "Failed to promote Figure18 source file: ",
      f
    )
  }
}

# Full formal multiple-testing tables
for (ds in datasets) {

  src_tests <- file.path(
    audit_root,
    paste0(
      ds,
      "_full18_nboot1000_conditionwide_BH_tests.tsv"
    )
  )

  dst_tests <- file.path(
    canonical_mt_dir,
    basename(
      src_tests
    )
  )

  if (
    !file.copy(
      src_tests,
      dst_tests,
      overwrite = TRUE
    )
  ) {
    stop(
      "Failed to promote ",
      ds,
      " multiple-testing table."
    )
  }

  fwrite(
    production_summaries[[ds]],
    file.path(
      canonical_mt_dir,
      paste0(
        ds,
        "_full18_nboot1000_conditionwide_BH_summary.tsv"
      )
    ),
    sep = "\t"
  )
}

combined_summary <- rbindlist(
  production_summaries,
  use.names = TRUE
)

fwrite(
  combined_summary,
  file.path(
    canonical_mt_dir,
    "CellChat_full18_nboot1000_conditionwide_BH_summary.tsv"
  ),
  sep = "\t"
)

# --------------------------------------------------------------------------
# Byte-level promotion validation
# --------------------------------------------------------------------------

if (
  unname(
    tools::md5sum(candidate_obj)
  ) !=
    unname(
      tools::md5sum(canonical_obj)
    )
) {
  stop(
    "Canonical CellChat object MD5 does not match validated candidate."
  )
}

for (f in source_files) {

  src_md5 <- unname(
    tools::md5sum(
      file.path(
        candidate_source_dir,
        f
      )
    )
  )

  dst_md5 <- unname(
    tools::md5sum(
      file.path(
        canonical_source_dir,
        f
      )
    )
  )

  if (src_md5 != dst_md5) {
    stop(
      "Canonical source-data MD5 mismatch: ",
      f
    )
  }
}

cat("\n============================================================\n")
cat("CELLCHAT CANONICAL PROMOTION COMPLETE\n")
cat("============================================================\n\n")

print(
  combined_summary
)

cat(
  "\nCanonical CellChat object:\n",
  canonical_obj,
  "\n",
  sep = ""
)

cat(
  "\nCanonical Figure18 source-data directory:\n",
  canonical_source_dir,
  "\n",
  sep = ""
)

cat(
  "\nCanonical multiple-testing directory:\n",
  canonical_mt_dir,
  "\n",
  sep = ""
)

cat(
  "\nPASS: canonical resources exactly match the validated corrected candidates.\n"
)
