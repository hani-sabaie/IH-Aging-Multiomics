# ===== Clean environment =====
rm(list = ls(all.names = TRUE))
gc()

# ===== Loading relevant libraries =====
library(data.table)
library(dplyr)
library(susieR)
library(coloc)
library(Matrix)
library(tibble)
library(purrr)
library(magrittr)
library(locuscomparer)
library(ggplot2)
library(cowplot)
library(ggrepel)

# ============================================================================ #

# ===== Inputs =====
gwas_file <- "C:/Users/Hani/Desktop/FM/UKB_SMAD3_250kb.txt"
bim_file <- "C:/Users/Hani/Desktop/plink2_win64_20251111/INHR/UKB_SMAD3_1000G_EUR.bim"
ld_file <- "C:/Users/Hani/Desktop/plink2_win64_20251111/INHR/UKB_SMAD3_1000G_EUR_LD.ld"
eqtl_file <- "C:/Users/Hani/Desktop/smr-1.3.1-win-x86_64/For_COLOC/myLiteCisEqtl.txt"

# GWAS case–control numbers
ncase <- 28707
ncontrol <- 343103
ntotal <- ncase + ncontrol
s <- ncase / ntotal # case fraction

# eQTL sample size 
N_eqtl <- 663              
sdY_eqtl <- 1 # assume standardized expression (sd ~ 1)

# ===== GWAS: load and basic QC =====
gwas <- fread(gwas_file)
gwas$rsID <- gwas$SNP

# Keep one row per SNP (smallest p)
gwas <- gwas %>%
  arrange(rsID, p) %>%
  group_by(rsID) %>%
  slice(1) %>%
  ungroup()

# Remove NA and zero betas
gwas <- gwas[complete.cases(gwas), ]
gwas <- gwas[gwas$b != 0, ]

# Add varbeta and case fraction
gwas <- gwas %>%
  mutate(
    s       = s,
    varbeta = se^2
  ) %>%
  dplyr::select(rsID, A1, A2, freq, b, se, p, pos, s, varbeta)

# ===== LD panel: BIM + LD matrix + basic QC =====
bim <- fread(
  bim_file,
  col.names = c("CHR","rsID","CM","BP_ld","A1_ld","A2_ld")
)

# Remove duplicated SNP IDs in BIM
dup_bim <- duplicated(bim$rsID)
bim2 <- bim[!dup_bim, ]

# Load PLINK LD matrix (no header, numeric only)
ld_matrix <- as.matrix(read.table(ld_file, header = FALSE, check.names = FALSE))
ld_matrix <- ld_matrix[!dup_bim, !dup_bim]

rownames(ld_matrix) <- bim2$rsID
colnames(ld_matrix) <- bim2$rsID

# Remove SNPs with invalid diagonal entries
diag_vals <- diag(ld_matrix)
bad_snps <- names(diag_vals)[diag_vals <= 0 | is.na(diag_vals)]
# bad_snps <- c(bad_snps, "rs141195834")

good_snps <- setdiff(rownames(ld_matrix), bad_snps)

ld_matrix_clean <- ld_matrix[good_snps, good_snps]
bim_clean <- bim2[bim2$rsID %in% good_snps, ]

# ===== Allele alignment between GWAS and LD =====
merged <- gwas %>%
  inner_join(
    bim_clean %>% dplyr::select(rsID, A1_ld, A2_ld),
    by = "rsID"
  )

A1_g <- merged$A1 # GWAS effect allele
A2_g <- merged$A2 # GWAS other allele
A1_r <- merged$A1_ld  # LD panel allele 1
A2_r <- merged$A2_ld  # LD panel allele 2

same <- (A1_g == A1_r & A2_g == A2_r)
flipped <- (A1_g == A2_r & A2_g == A1_r)

# Keep only SNPs with consistent alleles (allow flipped orientation)
keep  <- (same | flipped)
merged <- merged[keep, ]

# Build aligned beta and frequency
merged$beta_aligned <- merged$b
merged$freq_aligned <- merged$freq

# Flip beta and freq where alleles are reversed
flip_idx <- which(flipped[keep])
merged$beta_aligned[flip_idx] <- -merged$beta_aligned[flip_idx]
merged$freq_aligned[flip_idx] <- 1 - merged$freq_aligned[flip_idx]

# ===== MAF filter and intersect GWAS with LD =====
merged$MAF <- pmin(merged$freq_aligned, 1 - merged$freq_aligned)
merged <- merged[merged$MAF >= 0.01, ]

snps_gwas <- merged$rsID
snps_ld <- rownames(ld_matrix_clean)

final_snps <- intersect(snps_ld, snps_gwas)

# Subset LD to final SNP set
ld_matrix_final <- ld_matrix_clean[final_snps, final_snps]

# Make LD numerically positive-definite
ld_pd <- as.matrix(nearPD(ld_matrix_final, corr = TRUE)$mat)

# Subset merged GWAS to final SNP set and reorder to match LD
merged_final <- merged[merged$rsID %in% final_snps, ]
merged_final <- merged_final[match(final_snps, merged_final$rsID), ]

# Sanity checks
stopifnot(all(rownames(ld_pd) == final_snps))
stopifnot(all(merged_final$rsID == final_snps))

# ===== Build GWAS dataset1 on full GWAS+LD SNP set =====
snp <- merged_final$rsID
beta_aligned <- merged_final$beta_aligned
se_aligned <- merged_final$se
varbeta <- se_aligned^2
pos <- merged_final$pos

named_beta <- beta_aligned %>% set_names(snp)
named_varbeta <- varbeta %>% set_names(snp)

R <- ld_pd

snps_gwas_ld <- snp   # SNPs present in GWAS and LD (after all QC)

# ===== eQTL: load SMR-lite file and prepare SMAD3 eQTL =====
eqtl <- fread(eqtl_file)
eqtl <- eqtl %>%
  filter(
    BP > (67357940 - 1e6) &
      BP < (67487507 + 1e6)
  )

# Rename columns we need
setnames(eqtl,
         old = c("Chr","se"),
         new = c("CHR","SE"),
         skip_absent = TRUE)

# Keep only necessary columns
eqtl <- eqtl[, .(SNP, CHR, BP, A1, A2, Gene, b, SE, p)]

# Ensure numeric types
eqtl$b <- as.numeric(eqtl$b)
eqtl$SE <- as.numeric(eqtl$SE)
eqtl$p <- as.numeric(eqtl$p)
eqtl$BP <- as.numeric(eqtl$BP)

# Remove rows with missing essential values
eqtl <- eqtl[complete.cases(eqtl)]

# Keep SMAD3 only
eqtl <- eqtl[Gene == "SMAD3", ]

# Compute varbeta and z
eqtl[, `:=`(
  varbeta = SE^2,
  z       = b / SE
)]

# ===== Build common SNP set: GWAS+LD ∩ eQTL =====
snps_eqtl <- eqtl$SNP
common_snps <- intersect(snps_gwas_ld, snps_eqtl)

# Subset LD to common SNPs
R_common <- R[common_snps, common_snps]

# Subset GWAS to common SNPs and match order
gwas_sub <- data.frame(
  SNP = merged_final$rsID,
  p = merged_final$p, # add p-value here
  beta = merged_final$beta_aligned,
  se = merged_final$se,
  varbeta = merged_final$se^2,
  pos = merged_final$pos,
  s = s,
  stringsAsFactors = FALSE
)
gwas_sub <- gwas_sub[gwas_sub$SNP %in% common_snps, ]
gwas_sub <- gwas_sub[match(common_snps, gwas_sub$SNP), ]

# Subset eQTL to common SNPs and match order
eqtl_sub <- eqtl[eqtl$SNP %in% common_snps, ]
eqtl_sub <- eqtl_sub[match(common_snps, eqtl_sub$SNP), ]

# Sanity checks on alignment
stopifnot(all(rownames(R_common) == common_snps))
stopifnot(all(gwas_sub$SNP  == common_snps))
stopifnot(all(eqtl_sub$SNP  == common_snps))

# ===== Build dataset1 (GWAS) for coloc.susie and run SuSiE =====
named_beta_gwas <- set_names(gwas_sub$beta, gwas_sub$SNP)
named_varbeta_gwas <- set_names(gwas_sub$varbeta, gwas_sub$SNP)

dataset1 <- list(
  beta = named_beta_gwas,
  varbeta = named_varbeta_gwas,
  s = gwas_sub$s[1],
  type = "cc",
  snp = gwas_sub$SNP,
  LD  = R_common,
  position = gwas_sub$pos,
  N = ntotal
)

coloc:::check_dataset(dataset1, 1)

S3 <- runsusie(
  dataset1,
  nref = 503, coverage = 0.90
)
summary(S3)

# ===== Build dataset2 (eQTL) for coloc.susie and run SuSiE =====
named_beta_eqtl <- set_names(eqtl_sub$b, eqtl_sub$SNP)
named_varbeta_eqtl <- set_names(eqtl_sub$varbeta, eqtl_sub$SNP)

dataset2 <- list(
  beta = named_beta_eqtl,
  varbeta = named_varbeta_eqtl,
  N = N_eqtl,
  sdY = sdY_eqtl,
  type = "quant",
  snp = eqtl_sub$SNP,
  LD = R_common,
  position = eqtl_sub$BP
)

coloc:::check_dataset(dataset2, 2)

S4 <- runsusie(
  dataset2,
  nref = 503, coverage = 0.90
)
summary(S4)

# ===== Colocalization with coloc.susie =====
res <- coloc.susie(S3, S4)

res$summary[res$summary$PP.H4.abf > 0.75, ]
colsus <- res$summary[res$summary$PP.H4.abf > 0.75, ]

write.table(x = colsus, file = "UKB_coloc_susie_res.txt", sep = "\t")

# =====  LocusCompare + multiple SNP labels =====
# add_label
add_label <- function(merged, snps) {
  # Add a "label" column: show rsid only for SNPs in 'snps'
  merged$label <- ifelse(merged$rsid %in% snps, merged$rsid, "")
  return(merged)
}

# make_combined_plot  (supports multiple highlighted SNPs)
make_combined_plot <- function(merged, title1, title2, ld, chr, snp = NULL,
                               combine = TRUE, legend = TRUE,
                               legend_position = c("bottomright","topright","topleft"),
                               lz_ylab_linebreak = FALSE,
                               highlight_snps = NULL) {
  
  # Choose lead SNP (if snp is NULL, it will use the smallest p-value)
  lead_snp <- get_lead_snp(merged, snp)
  
  # SNPs to label: lead SNP + highlighted SNPs
  snps_to_label <- unique(c(lead_snp, highlight_snps))
  snps_to_label <- snps_to_label[snps_to_label %in% merged$rsid]
  
  # Assign LD-based colors using the lead SNP
  color <- assign_color(merged$rsid, lead_snp, ld)
  
  # Shapes:
  #  - lead SNP: diamond (23)
  #  - highlighted SNPs: triangle (24)
  #  - others: circle (21)
  shape <- ifelse(merged$rsid == lead_snp, 23,
                  ifelse(merged$rsid %in% snps_to_label, 24, 21))
  names(shape) <- merged$rsid
  
  # Sizes:
  #  - lead SNP: biggest
  #  - highlighted SNPs: medium
  #  - others: small
  size <- ifelse(merged$rsid == lead_snp, 2,
                 ifelse(merged$rsid %in% snps_to_label, 2, 2))
  names(size) <- merged$rsid
  
  # Add text labels for selected SNPs
  merged <- add_label(merged, snps_to_label)
  
  # Scatter plot (locuscompare)
  p1 <- make_scatterplot(merged, title1, title2, color,
                         shape, size, legend, legend_position)
  
  # Locuszoom for study 1
  metal1 <- merged[, c("rsid", "logp1", "chr", "pos", "label")]
  colnames(metal1)[colnames(metal1) == "logp1"] <- "logp"
  p2 <- make_locuszoom(metal1, title1, chr, color, shape, size, lz_ylab_linebreak)
  
  # Locuszoom for study 2
  metal2 <- merged[, c("rsid", "logp2", "chr", "pos", "label")]
  colnames(metal2)[colnames(metal2) == "logp2"] <- "logp"
  p3 <- make_locuszoom(metal2, title2, chr, color, shape, size, lz_ylab_linebreak)
  
  # Combine panels or return individually
  if (combine) {
    p2 <- p2 + theme(axis.text.x = element_blank(),
                     axis.title.x = element_blank())
    p4 <- cowplot::plot_grid(p2, p3, align = "v", nrow = 2,
                             rel_heights = c(0.8, 1))
    p5 <- cowplot::plot_grid(p1, p4)
    return(p5)
  } else {
    return(list(locuscompare = p1, locuszoom1 = p2, locuszoom2 = p3))
  }
}

# locuscompare  (new argument: highlight_snps)
locuscompare <- function(in_fn1, in_fn2,
                         marker_col1 = "rsid", pval_col1 = "pval",
                         title1 = "eQTL",
                         marker_col2 = "rsid", pval_col2 = "pval",
                         title2 = "GWAS",
                         snp = NULL, # lead SNP (optional)
                         highlight_snps = NULL, # extra SNPs to label
                         population = "EUR",
                         combine = TRUE, legend = TRUE,
                         legend_position = c("bottomright","topright","topleft"),
                         lz_ylab_linebreak = FALSE,
                         genome = c("hg19","hg38")) {
  
  # Read and merge two studies
  d1 <- read_metal(in_fn1, marker_col1, pval_col1)
  d2 <- read_metal(in_fn2, marker_col2, pval_col2)
  merged <- merge(d1, d2, by = "rsid", suffixes = c("1", "2"), all = FALSE)
  
  # Add chr/pos
  genome <- match.arg(genome)
  merged <- get_position(merged, genome)
  
  chr <- unique(merged$chr)
  if (length(chr) != 1)
    stop("There must be one and only one chromosome.")
  
  # Choose lead SNP (if snp is NULL, use minimum p-value across studies)
  lead_snp <- get_lead_snp(merged, snp)
  
  # Retrieve LD for the lead SNP
  ld <- retrieve_LD(chr, lead_snp, population)
  
  # Build combined plot (scatter + 2 locuszoom)
  p <- make_combined_plot(
    merged, title1, title2, ld, chr,
    snp = lead_snp,
    combine = combine,
    legend = legend,
    legend_position = legend_position,
    lz_ylab_linebreak = lz_ylab_linebreak,
    highlight_snps = highlight_snps
  )
  
  return(p)
}

df_lc_gwas <- gwas_sub %>%
  transmute(
    rsid = SNP,
    pval = p
  )

df_lc_eqtl <- eqtl_sub %>%
  transmute(
    rsid = SNP,
    pval = p
  )

df_lc_gwas <- df_lc_gwas[df_lc_gwas$rsid %in% common_snps, ]
df_lc_eqtl <- df_lc_eqtl[df_lc_eqtl$rsid %in% common_snps, ]

plot <- locuscompare(
  in_fn1 = df_lc_gwas,
  in_fn2 = df_lc_eqtl,
  title1 = "GWAS",
  title2 = "eQTL",
  genome = "hg19",
  population = "EUR",
  snp = "rs7181556",   # lead SNP
  highlight_snps = c("rs10152595", "rs7181877", "rs35874463", "rs12912045")  # extra SNPs
)

png("UKB_locuscompare.png",
    width = 3600, height = 1800, res = 300)
plot
dev.off()
