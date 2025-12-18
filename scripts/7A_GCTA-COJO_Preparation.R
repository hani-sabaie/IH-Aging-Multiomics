# ===== Clean environment =====
rm(list = ls(all.names = TRUE))
gc()

# ===== Loading relevant libraries =====
library(data.table)

# ============================================================================ #

# ===== Preparation =====
# Locus file
ukb_loc <- fread("UKB_SMAD3_250kb.txt")

# GCTA expects columns 
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

fwrite(gcta_ukb,
       "UKB_SMAD3_for_GCTA.tsv",
       sep = "\t")
