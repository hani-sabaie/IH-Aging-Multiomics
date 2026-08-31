# ===== Clean environment =====
rm(list = ls(all.names = TRUE))
gc()

# ===== Loading relevant libraries =====
library(data.table)

# ===== Repository paths =====
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (length(file_arg) == 1) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg)))
  repo_root <- normalizePath(file.path(script_dir, ".."))
} else {
  repo_root <- normalizePath(".")
}

ukb_dir <- file.path(repo_root, "data", "GWAS", "UKB")
fin_dir <- file.path(repo_root, "data", "GWAS", "FinnGen")

# ===== Input locus files =====
ukb_file <- Sys.getenv(
  "UKB_SMAD3_LOCUS",
  unset = file.path(ukb_dir, "UKB_SMAD3_250kb.txt")
)

fin_file <- Sys.getenv(
  "FINNGEN_SMAD3_LOCUS",
  unset = file.path(fin_dir, "Finn_SMAD3_250kb.txt")
)

if (!file.exists(ukb_file)) {
  stop("UKB locus file not found: ", ukb_file)
}

if (!file.exists(fin_file)) {
  stop("FinnGen locus file not found: ", fin_file)
}

# ============================================================================ #
# ===== UK Biobank =====
ukb_loc <- fread(ukb_file)

gcta_ukb <- ukb_loc[, .(
  SNP = SNP,
  A1 = A1,
  A2 = A2,
  freq = freq,
  b = b,
  se = se,
  p = p,
  N = n
)]

fwrite(
  gcta_ukb,
  file.path(ukb_dir, "UKB_SMAD3_for_GCTA.tsv"),
  sep = "\t"
)

# ============================================================================ #
# ===== FinnGen =====
fin_loc <- fread(fin_file)

gcta_fin <- fin_loc[, .(
  SNP = SNP,
  A1 = A1,
  A2 = A2,
  freq = freq,
  b = b,
  se = se,
  p = p,
  N = n
)]

fwrite(
  gcta_fin,
  file.path(fin_dir, "Finn_SMAD3_for_GCTA.tsv"),
  sep = "\t"
)

message("Prepared GCTA-COJO summary files for UKB and FinnGen.")
