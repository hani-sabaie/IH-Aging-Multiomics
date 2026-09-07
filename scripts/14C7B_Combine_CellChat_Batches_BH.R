# ============================================================================
# Reviewer C7 audit:
# Combine completed CellChat nboot=1000 batches and apply the final
# condition-wide multiple-testing correction.
#
# Usage:
#   Rscript scripts/14C7B_Combine_CellChat_Batches_BH.R Young
#   Rscript scripts/14C7B_Combine_CellChat_Batches_BH.R Aged
#
# Primary correction:
#   plus-one empirical permutation P = (b + 1) / (B + 1)
#   BH across the COMPLETE condition-specific LR x source x target family.
#
# No canonical CellChat object or Figure 18 file is modified.
# ============================================================================

rm(list = ls(all.names = TRUE))

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

if (
  length(args) != 1L ||
  !(args[1] %in% c("Young", "Aged"))
) {
  stop("Provide exactly one condition: Young or Aged.")
}

ds <- args[1]

cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd, value = TRUE)

if (length(file_arg) == 1L) {
  script_dir <- dirname(
    normalizePath(sub("^--file=", "", file_arg))
  )
  repo_root <- normalizePath(
    file.path(script_dir, "..")
  )
} else {
  repo_root <- normalizePath(".")
}

base_dir <- file.path(
  repo_root,
  "processed_results",
  "12_CellChat",
  "multiple_testing_audit",
  "full18_nboot1000_BH"
)

batch_dir <- file.path(
  base_dir,
  "batches"
)

dir.create(
  base_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

expected_n_lr <- if (ds == "Young") 488L else 468L
expected_n_tests <- expected_n_lr * 7L * 7L
batch_size <- 25L

starts <- seq.int(
  1L,
  expected_n_lr,
  by = batch_size
)

ends <- pmin(
  starts + batch_size - 1L,
  expected_n_lr
)

tags <- sprintf(
  "%s_LR%04d-%04d",
  ds,
  starts,
  ends
)

test_files <- file.path(
  batch_dir,
  paste0(
    tags,
    "_tests.tsv"
  )
)

validation_files <- file.path(
  batch_dir,
  paste0(
    tags,
    "_validation.tsv"
  )
)

missing_test <- test_files[
  !file.exists(test_files)
]

missing_validation <- validation_files[
  !file.exists(validation_files)
]

if (length(missing_test) > 0L) {
  stop(
    "Missing batch test file(s):\n",
    paste(
      missing_test,
      collapse = "\n"
    )
  )
}

if (length(missing_validation) > 0L) {
  stop(
    "Missing batch validation file(s):\n",
    paste(
      missing_validation,
      collapse = "\n"
    )
  )
}

cat("============================================================\n")
cat("COMBINE CELLCHAT BATCHES + CONDITION-WIDE BH\n")
cat("============================================================\n\n")

cat("Condition             :", ds, "\n")
cat("Expected LR family    :", expected_n_lr, "\n")
cat("Expected formal tests :", expected_n_tests, "\n")
cat("Expected batches      :", length(tags), "\n\n")

# --------------------------------------------------------------------------
# Validate every batch
# --------------------------------------------------------------------------

validation <- rbindlist(
  lapply(
    validation_files,
    fread
  ),
  use.names = TRUE,
  fill = TRUE
)

if (nrow(validation) != length(tags)) {
  stop(
    "Expected ",
    length(tags),
    " validation rows; found ",
    nrow(validation),
    "."
  )
}

if (
  !"probability_match_1e12" %in%
    names(validation)
) {
  stop(
    "Batch validation tables lack probability_match_1e12."
  )
}

if (
  any(
    is.na(validation$probability_match_1e12)
  ) ||
  !all(
    validation$probability_match_1e12
  )
) {
  stop(
    "At least one batch failed exact historical probability validation."
  )
}

if (
  sum(validation$n_LR) != expected_n_lr
) {
  stop(
    "Validated batch LR total is ",
    sum(validation$n_LR),
    "; expected ",
    expected_n_lr,
    "."
  )
}

if (
  sum(validation$n_tests) != expected_n_tests
) {
  stop(
    "Validated batch test total is ",
    sum(validation$n_tests),
    "; expected ",
    expected_n_tests,
    "."
  )
}

cat(
  "All batch probability validations:",
  all(validation$probability_match_1e12),
  "\n"
)

cat(
  "Total batch runtime (seconds)       :",
  sum(validation$elapsed_seconds),
  "\n\n"
)

# --------------------------------------------------------------------------
# Combine all tests
# --------------------------------------------------------------------------

tests <- rbindlist(
  lapply(
    test_files,
    fread
  ),
  use.names = TRUE,
  fill = TRUE
)

required_cols <- c(
  "dataset",
  "LR_family_index",
  "interaction_name",
  "source",
  "target",
  "probability",
  "p_raw_nboot1000",
  "permutation_exceedance_count",
  "p_plus1"
)

missing_cols <- setdiff(
  required_cols,
  names(tests)
)

if (length(missing_cols) > 0L) {
  stop(
    "Combined test table is missing columns: ",
    paste(
      missing_cols,
      collapse = ", "
    )
  )
}

if (nrow(tests) != expected_n_tests) {
  stop(
    "Combined table contains ",
    nrow(tests),
    " tests; expected ",
    expected_n_tests,
    "."
  )
}

if (
  !all(tests$dataset == ds)
) {
  stop(
    "Combined batches contain an unexpected dataset label."
  )
}

# --------------------------------------------------------------------------
# Validate complete LR family
# --------------------------------------------------------------------------

lr_indices <- sort(
  unique(
    tests$LR_family_index
  )
)

if (
  !identical(
    lr_indices,
    seq_len(expected_n_lr)
  )
) {
  stop(
    "Combined batches do not contain exactly LR indices 1-",
    expected_n_lr,
    "."
  )
}

lr_counts <- tests[
  ,
  .N,
  by = LR_family_index
]

if (
  nrow(lr_counts) != expected_n_lr ||
  any(lr_counts$N != 49L)
) {
  stop(
    "Each LR must contribute exactly 49 source-target tests."
  )
}

# Every LR index must map to exactly one interaction name.
lr_name_check <- tests[
  ,
  .(
    n_interaction_names =
      uniqueN(interaction_name)
  ),
  by = LR_family_index
]

if (
  any(
    lr_name_check$n_interaction_names != 1L
  )
) {
  stop(
    "At least one LR index maps to multiple interaction names."
  )
}

# Check duplicate formal hypotheses.
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
    "Duplicate LR-source-target hypothesis detected."
  )
}

cat(
  "Combined tests validated           :",
  nrow(tests),
  "\n"
)

cat(
  "Unique LR interactions             :",
  uniqueN(tests$interaction_name),
  "\n"
)

cat(
  "Tests per LR                       : 49\n\n"
)

# --------------------------------------------------------------------------
# Condition-wide correction
# --------------------------------------------------------------------------

if (
  anyNA(tests$p_plus1) ||
  any(
    tests$p_plus1 <= 0 |
      tests$p_plus1 > 1
  )
) {
  stop(
    "Invalid plus-one empirical P values."
  )
}

# Sensitivity: BH applied directly to CellChat b/B P values.
tests[
  ,
  p_raw_BH_condition :=
    p.adjust(
      p_raw_nboot1000,
      method = "BH"
    )
]

# Primary analysis:
# finite plus-one empirical P, followed by BH across the COMPLETE
# condition-specific family.
tests[
  ,
  p_plus1_BH_condition :=
    p.adjust(
      p_plus1,
      method = "BH"
    )
]

tests[
  ,
  raw_nboot1000_sig :=
    p_raw_nboot1000 < 0.05
]

tests[
  ,
  raw_BH_sig :=
    p_raw_BH_condition < 0.05
]

tests[
  ,
  plus1_BH_sig :=
    p_plus1_BH_condition < 0.05
]

# Sort back into canonical LR/source/target family order.
setorder(
  tests,
  LR_family_index,
  source,
  target
)

# --------------------------------------------------------------------------
# Summary
# --------------------------------------------------------------------------

summary_dt <- data.table(
  dataset = ds,
  n_batches = length(tags),
  n_LR = expected_n_lr,
  n_source_groups = 7L,
  n_target_groups = 7L,
  n_tests = nrow(tests),

  n_raw_p_zero =
    sum(
      tests$p_raw_nboot1000 == 0
    ),

  n_raw_p_lt_0_05 =
    sum(
      tests$p_raw_nboot1000 < 0.05
    ),

  n_raw_BH_lt_0_05 =
    sum(
      tests$p_raw_BH_condition < 0.05
    ),

  n_plus1_BH_lt_0_05 =
    sum(
      tests$p_plus1_BH_condition < 0.05
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

  max_probability_validation_delta =
    max(
      validation$max_abs_probability_delta
    ),

  all_probability_batches_match =
    all(
      validation$probability_match_1e12
    ),

  total_batch_elapsed_seconds =
    sum(
      validation$elapsed_seconds
    )
)

# --------------------------------------------------------------------------
# Write audit outputs
# --------------------------------------------------------------------------

fwrite(
  tests,
  file.path(
    base_dir,
    paste0(
      ds,
      "_full18_nboot1000_conditionwide_BH_tests.tsv"
    )
  ),
  sep = "\t"
)

fwrite(
  summary_dt,
  file.path(
    base_dir,
    paste0(
      ds,
      "_full18_nboot1000_conditionwide_BH_summary.tsv"
    )
  ),
  sep = "\t"
)

fwrite(
  validation,
  file.path(
    base_dir,
    paste0(
      ds,
      "_full18_batch_validation_summary.tsv"
    )
  ),
  sep = "\t"
)

cat("\n============================================================\n")
cat("FINAL CONDITION-WIDE SUMMARY\n")
cat("============================================================\n\n")

print(summary_dt)

cat("\nOutput directory:\n")
cat(base_dir, "\n")

cat(
  "\nNo canonical CellChat object or Figure 18 file was modified.\n"
)
