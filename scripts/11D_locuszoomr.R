# ===== Clean environment =====
rm(list = ls(all.names = TRUE))
gc()

# ===== Loading relevant libraries =====
library(locuszoomr)
library(EnsDb.Hsapiens.v75)
library(data.table)
library(ggplot2)

# ============================================================================ #
# ===== locusplot =====
gwas_ukb <- fread("UKB_gwas_for_xqtlbiolinks.txt")
loc_ukb <- locus(data = gwas_ukb, gene = 'SMAD3', fix_window = 1e6,
             ens_db = "EnsDb.Hsapiens.v75", LD = "r2",index_snp = "rs7181556")
loc_ukb <- link_LD(loc_ukb, token = "",genome_build = "grch37")
loc_ukb <- link_recomb(loc_ukb, genome = "hg19")
png("locusplot_UKB.png", width = 3000, height = 3000, res = 300)
locus_plot(loc_ukb, labels = c("rs7181556", "rs10152595",
                           "rs35874463", "rs12912045",
                           "rs7181877"), 
           label_x = c(4, -5, 5),
           recomb_offset = 0.1,
           highlight = "SMAD3",
           pcutoff = 1e-6)
dev.off()

gwas_finn <- fread("Finn_gwas_for_xqtlbiolinks.txt")
loc_finn <- locus(data = gwas_finn, gene = 'SMAD3', fix_window = 1e6,
                  ens_db = "EnsDb.Hsapiens.v75", LD = "r2",index_snp = "rs11315136")
loc_finn <- link_LD(loc_finn, token = "",genome_build = "grch37")
loc_finn <- link_recomb(loc_finn, genome = "hg19")
png("locusplot_Finn.png", width = 3000, height = 3000, res = 300)
locus_plot(loc_finn, labels = c("rs11315136", "rs1965269", "rs3784681"), 
            label_x = c(-4, -5, 5),
            recomb_offset = 0.1,
            highlight = "SMAD3",
            pcutoff = 1e-4)
dev.off()

# ===== eqtlplot =====
loc_ukb <- link_eqtl(loc_ukb, token = "")
table(loc_ukb$LDexp$Gene_Symbol)
table(loc_ukb$LDexp$Tissue)
png("eqtlplot_UKB.png", width = 3000, height = 3000, res = 300)
overlay_plot(loc_ukb, eqtl_gene = "SMAD3",tissue = "Adipose - Subcutaneous",
             labels = c("rs7181556", "rs10152595",
                        "rs35874463", "rs12912045",
                        "rs7181877"),
             label_x = c(4, -5, 5), highlight = "SMAD3",
             recomb_col = NA, pcutoff = NA)
dev.off()

loc_finn <- link_eqtl(loc_finn, token = "")
png("eqtlplot_Finn.png", width = 3000, height = 3000, res = 300)
overlay_plot(loc_finn, eqtl_gene = "SMAD3",tissue = "Adipose - Subcutaneous",
             labels = c("rs11315136", "rs1965269", "rs3784681"), 
             label_x = c(-4, -5, 5), highlight = "SMAD3",
             recomb_col = NA,pcutoff = NA)
dev.off()

