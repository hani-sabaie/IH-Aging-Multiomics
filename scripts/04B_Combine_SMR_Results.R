# ===== Combine SMR/HEIDI results across studies and GTEx tissues =====

library(data.table)

# Locate repository root from this script
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (length(file_arg) == 1) {
  script_path <- sub("^--file=", "", file_arg)
  script_dir <- dirname(normalizePath(script_path))
  repo_root <- normalizePath(file.path(script_dir, ".."))
} else {
  repo_root <- normalizePath(".")
}

smr_dir <- file.path(
  repo_root,
  "processed_results",
  "06_SMR_HEIDI"
)

# Expected 2 GWAS x 4 GTEx tissues
expected_files <- c(
  "UKB_Adipose_Subcutaneous.smr",
  "UKB_Adipose_Visceral_Omentum.smr",
  "UKB_Muscle_Skeletal.smr",
  "UKB_Cells_Cultured_fibroblasts.smr",
  "FinnGen_Adipose_Subcutaneous.smr",
  "FinnGen_Adipose_Visceral_Omentum.smr",
  "FinnGen_Muscle_Skeletal.smr",
  "FinnGen_Cells_Cultured_fibroblasts.smr"
)

missing_files <- expected_files[
  !file.exists(file.path(smr_dir, expected_files))
]

if (length(missing_files) > 0) {
  stop(
    "Missing SMR result file(s): ",
    paste(missing_files, collapse = ", ")
  )
}

# Read and annotate each result file
smr_list <- lapply(expected_files, function(filename) {

  dat <- fread(file.path(smr_dir, filename))

  # Verify expected SMR output columns
  required_cols <- c(
    "probeID", "ProbeChr", "Gene", "Probe_bp",
    "topSNP", "topSNP_chr", "topSNP_bp",
    "A1", "A2", "Freq",
    "b_GWAS", "se_GWAS", "p_GWAS",
    "b_eQTL", "se_eQTL", "p_eQTL",
    "b_SMR", "se_SMR", "p_SMR",
    "p_HEIDI", "nsnp_HEIDI"
  )

  missing_cols <- setdiff(required_cols, names(dat))

  if (length(missing_cols) > 0) {
    stop(
      "Unexpected schema in ", filename,
      "; missing column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }

  stem <- sub("\\.smr$", "", filename)

  study <- if (grepl("^UKB_", stem)) {
    "UKB"
  } else if (grepl("^FinnGen_", stem)) {
    "FinnGen"
  } else {
    stop("Cannot determine study from filename: ", filename)
  }

  tissue <- sub("^(UKB|FinnGen)_", "", stem)

  # Human-readable tissue label
  tissue <- gsub("_", " ", tissue)

  dat[, `:=`(
    study = study,
    tissue = tissue
  )]

  setcolorder(
    dat,
    c("study", "tissue",
      setdiff(names(dat), c("study", "tissue")))
  )

  dat
})

smr_all <- rbindlist(smr_list, use.names = TRUE, fill = FALSE)

# Save combined author-generated processed result
out_file <- file.path(
  smr_dir,
  "SMR_all_studies_all_tissues.csv"
)

fwrite(smr_all, out_file)

# This script only combines the complete SMR/HEIDI results.
# Statistical significance is defined downstream using the revised
# multiple-testing framework (cohort-wide BH correction for UKB discovery,
# followed by targeted BH correction for the carried-forward FinnGen
# replication hypotheses). No nominal-P "significant" subset is generated
# here.

cat("Combined SMR rows:", nrow(smr_all), "\n")
cat("Studies:", paste(unique(smr_all$study), collapse = ", "), "\n")
cat("Tissues:", paste(unique(smr_all$tissue), collapse = " | "), "\n")
cat("Output:", out_file, "\n")
