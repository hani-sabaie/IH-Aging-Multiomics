# ===== Clean environment =====
rm(list = ls(all.names = TRUE))
gc()

# ===== Loading relevant libraries =====library(data.table)
library(dplyr)
library(stringr)

# ============================================================================ #
# ===== UKB dataset =====
gwas <- fread("shortfinalINGallMAFge05padjplusxxy.txt")

lambda_gc <- 1.084
N_gwas <- 371810   

gwas_xqtl <- gwas %>%
  transmute(
    SNP = SNPID,
    A1 = Allele2,
    A2 = Allele1,
    freq  = AFAllele2,
    b = BETA,
    # back-calc raw p and SE from padj and lambdaGC
    chisq_gc = qchisq(1 - padj, df = 1),
    chisq_raw = chisq_gc * lambda_gc,
    z_raw = sqrt(chisq_raw),
    se = abs(b) / pmax(z_raw, .Machine$double.eps),
    p = 2 * pnorm(-abs(z_raw)),
    n = N_gwas,
    pos = POS,
    chr = CHR
  ) %>%
  select(SNP, A1, A2, freq, b, se, p, n, pos, chr)

head(gwas)
head(gwas_xqtl)

# write.table(gwas_xqtl, "UKB_gwas_for_xqtlbiolinks.txt", quote = FALSE, sep = "\t", row.names = FALSE)

gwas_smr <- gwas_xqtl %>%
  select(SNP, A1, A2, freq, b, se, p, n)

write.table(gwas_smr, "UKB_gwas_for_smr.txt", quote = FALSE, sep = "\t", row.names = FALSE)

# ============================================================================ #
# ===== FinnGen dataset =====
gwas <- fread("finn-b-K11_HERING.vcf.gz", skip = "#CHROM")
dt <- as.data.table(gwas)
ph <- setdiff(colnames(gwas),
                     c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"))
fields <- c("ES","SE","LP","AF","IDx")
dt[, c("ES","SE","LP","AF","IDx") :=
     tstrsplit(get(ph), ":", fixed = TRUE, fill = NA)]
numcols <- c("ES","SE","LP","AF")
dt[, (numcols) := lapply(.SD, as.numeric), .SDcols = numcols]
dt[, P := 10^(-LP)]
dt[, AF_info := as.numeric(sub(".*AF=([^;]+).*", "\\1", INFO))]
dt[, Freq := fifelse(is.na(AF), AF_info, AF)]

N_gwas <- 207653 
colnames(dt)[colnames(dt) == "#CHROM"] <- "CHROM"

gwas_xqtl <- dt[!is.na(ES) & !is.na(SE) & !is.na(P),
               .(SNP = ifelse(!is.na(IDx) & IDx != "", IDx, ID),
                 A1 = ALT,
                 A2 = REF,
                 freq = Freq,
                 b = ES,
                 se = SE,
                 p = P,
                 n = N_gwas,
                 pos = POS,
                 chr = CHROM)
]

head(dt)
head(gwas_xqtl)

# fwrite(gwas_xqtl, "Finn_gwas_for_xqtlbiolinks.txt", sep = "\t")

gwas_smr <- gwas_xqtl %>%
  select(SNP, A1, A2, freq, b, se, p, n)

write.table(gwas_smr, "Finn_gwas_for_smr.txt", quote = FALSE, sep = "\t", row.names = FALSE)
