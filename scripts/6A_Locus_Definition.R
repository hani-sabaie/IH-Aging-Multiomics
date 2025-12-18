# ===== Clean environment =====
rm(list = ls(all.names = TRUE))
gc()

# ===== Loading relevant libraries =====
library(data.table)
library(dplyr)

# ============================================================================ #
# ===== SMAD3 locus (hg19) ± 250kb =====
chr <- 15L
start <- 67107940L
end <- 67737507L

# Load UKB GWAS 
ukb <- fread("C:\\Users\\Hani\\Desktop\\Hernia\\script\\UKB_gwas_for_xqtlbiolinks.txt")
ukb <- ukb %>% 
  transmute(SNP, chr, pos, A1, A2, freq, b, se, p, n)

ukb_loc <- ukb[chr == chr & pos >= start & pos <= end]
fwrite(ukb_loc, "UKB_SMAD3_250kb.txt", sep = "\t")

# Same for FinnGen
fin <- fread("C:\\Users\\Hani\\Desktop\\Hernia\\script\\Finn_gwas_for_xqtlbiolinks.txt")
fin <- fin %>% 
  transmute(SNP, chr, pos, A1, A2, freq, b, se, p, n)
fin_loc <- fin[chr == chr & pos >= start & pos <= end]
fwrite(fin_loc, "Finn_SMAD3_250kb.txt", sep = "\t")

# SNP list in SMAD3 locus (for LD extraction)
all_snps_ukb <- sort(unique(ukb_loc$SNP))
rs_file_ukb  <- "SMAD3_UKB_rs.tsv"
fwrite(data.table(SNP = all_snps_ukb), rs_file_ukb,
       col.names = FALSE, sep = "\t")

all_snps_finn <- sort(unique(fin_loc$SNP))
rs_file_finn  <- "SMAD3_Finn_rs.tsv"
fwrite(data.table(SNP = all_snps_finn), rs_file_finn,
       col.names = FALSE, sep = "\t")
