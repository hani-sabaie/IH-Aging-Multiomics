# ===== Clean environment =====
rm(list = ls(all.names = TRUE))
gc()

# ===== Loading relevant libraries =====
library(data.table)
library(dplyr)

# ============================================================================ #
# ===== SMAD3 locus (hg19) ± 250kb =====
target_chr <- 15L
locus_start <- 67107940L
locus_end <- 67737507L

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

dir.create(ukb_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fin_dir, recursive = TRUE, showWarnings = FALSE)

# ===== Input files =====
ukb_file <- Sys.getenv(
  "UKB_GWAS_XQTL",
  unset = file.path(ukb_dir, "UKB_gwas_for_xqtlbiolinks.txt")
)

fin_file <- Sys.getenv(
  "FINNGEN_GWAS_XQTL",
  unset = file.path(fin_dir, "Finn_gwas_for_xqtlbiolinks.txt")
)

# ===== Load UKB GWAS =====
ukb <- fread(ukb_file)

ukb <- ukb %>%
  transmute(SNP, chr, pos, A1, A2, freq, b, se, p, n)

ukb_loc <- ukb[
  chr == target_chr &
    pos >= locus_start &
    pos <= locus_end
]

fwrite(ukb_loc, file.path(ukb_dir, "UKB_SMAD3_250kb.txt"), sep = "\t")

# ===== Load FinnGen GWAS =====
fin <- fread(fin_file)

fin <- fin %>%
  transmute(SNP, chr, pos, A1, A2, freq, b, se, p, n)

fin_loc <- fin[
  chr == target_chr &
    pos >= locus_start &
    pos <= locus_end
]

fwrite(fin_loc, file.path(fin_dir, "Finn_SMAD3_250kb.txt"), sep = "\t")

# SNP list in SMAD3 locus (for LD extraction)
all_snps_ukb <- sort(unique(ukb_loc$SNP))
rs_file_ukb <- file.path(ukb_dir, "SMAD3_UKB_rs.tsv")
fwrite(data.table(SNP = all_snps_ukb), rs_file_ukb,
       col.names = FALSE, sep = "\t")

all_snps_finn <- sort(unique(fin_loc$SNP))
rs_file_finn <- file.path(fin_dir, "SMAD3_Finn_rs.tsv")
fwrite(data.table(SNP = all_snps_finn), rs_file_finn,
       col.names = FALSE, sep = "\t")
