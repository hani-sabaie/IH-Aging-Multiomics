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

# ============================================================================ #
# ===== Inputs =====
# Paths
gwas_file <- "C:/Users/Hani/Desktop/FM/Finn_SMAD3_250kb.txt"
bim_file <- "C:/Users/Hani/Desktop/plink2_win64_20251111/INHR/Finn_SMAD3_1000G_EUR.bim"
ld_file <- "C:/Users/Hani/Desktop/plink2_win64_20251111/INHR/Finn_SMAD3_1000G_EUR_LD.ld"

# Case/control numbers
ncase <- 17096
ncontrol <- 190557
ntotal <- ncase + ncontrol
s <- ncase / ntotal # case fraction

# ===== Load GWAS and keep best row per rsID =====
gwas <- fread(gwas_file)
gwas$rsID <- gwas$SNP

# Keep only one row per rsID: the one with the smallest p-value
gwas <- gwas %>%
  arrange(rsID, p) %>%
  group_by(rsID) %>%
  slice(1) %>%
  ungroup()

# Basic QC: remove incomplete rows and zero betas
gwas <- gwas[complete.cases(gwas), ]
gwas <- gwas[gwas$b != 0, ]

# Add case fraction and varbeta
gwas <- gwas %>%
  mutate(
    s       = s,
    varbeta = se^2
  ) %>%
  dplyr::select(rsID, A1, A2, freq, b, se, p, pos, s, varbeta)

# ===== Load BIM and remove duplicated rsIDs =====
bim <- fread(
  bim_file,
  col.names = c("CHR","rsID","CM","BP_ld","A1_ld","A2_ld")
)

# Identify duplicated SNP IDs
dup_bim <- duplicated(bim$rsID)

# Keep only unique SNPs
bim2 <- bim[!dup_bim, ]

# ===== Load LD matrix and align with BIM ordering =====
# Load PLINK LD matrix
ld_matrix <- as.matrix(
  read.table(ld_file, header = FALSE, check.names = FALSE)
)

# Apply the same duplicated filter to LD matrix rows/columns
ld_matrix <- ld_matrix[!dup_bim, !dup_bim]

# Assign SNP names to LD matrix
rownames(ld_matrix) <- bim2$rsID
colnames(ld_matrix) <- bim2$rsID

# ===== Remove SNPs with LD diagonal <= 0 (bad variance/no LD) =====
diag_vals <- diag(ld_matrix)

# SNPs with zero or negative diagonal entries are invalid
bad_snps <- names(diag_vals)[diag_vals <= 0 | is.na(diag_vals)]
# bad_snps <- c(bad_snps, "rs141195834")

# Keep only good SNPs
good_snps <- setdiff(rownames(ld_matrix), bad_snps)

ld_matrix_clean <- ld_matrix[good_snps, good_snps]

# Also trim BIM to good SNPs
bim_clean <- bim2[bim2$rsID %in% good_snps, ]

# ===== Allele alignment between GWAS and LD (BIM) =====
# Merge GWAS with BIM alleles, restricted to SNPs that exist in LD
merged <- gwas %>%
  inner_join(
    bim_clean %>% dplyr::select(rsID, A1_ld, A2_ld),
    by = "rsID"
  )

# Rename for clarity
A1_g <- merged$A1 # GWAS effect allele
A2_g <- merged$A2 # GWAS non-effect allele
A1_r <- merged$A1_ld  # LD panel allele 1
A2_r <- merged$A2_ld  # LD panel allele 2

# Allele configurations
same <- (A1_g == A1_r & A2_g == A2_r)
flipped <- (A1_g == A2_r & A2_g == A1_r)

# Keep only SNPs with consistent, non-ambiguous alleles
keep <- (same | flipped)
merged <- merged[keep, ]

# Build aligned beta and freq
merged$beta_aligned <- merged$b
merged$freq_aligned <- merged$freq

# Indices (within 'merged') of flipped SNPs
flip_idx <- which(flipped[keep])
flip_snps <- merged$rsID[flip_idx]

# ===== Diagnostic plot 1: marginal z, highlight flipped SNPs =====
z_raw <- merged$beta_aligned / merged$se

plot_z_with_highlight <- function(z, snp_ids, highlight_snps = NULL,
                                  main = "Marginal Associations") {
  snp_index <- seq_along(z)
  
  plot(
    snp_index, z,
    pch = 16,
    col = "#767676",
    main = main,
    xlab = "SNP index along locus",
    ylab = "z-scores"
  )
  
  # highlight arbitrary SNPs by rsID (e.g. flipped SNPs)
  if (!is.null(highlight_snps)) {
    idx_h <- which(snp_ids %in% highlight_snps)
    if (length(idx_h) > 0) {
      points(
        snp_index[idx_h],
        z[idx_h],
        col = "yellow",
        pch = 16,
        cex = 1.2
      )
    }
  }
}

png(filename = "Finn_marg_asso.png", width = 3000, height = 1500, res = 300)
plot_z_with_highlight(z_raw, snp_ids = merged$rsID, highlight_snps = flip_snps)
dev.off()

# Now actually flip beta and freq for flipped SNPs
merged$beta_aligned[flip_idx] <- -merged$beta_aligned[flip_idx]
merged$freq_aligned[flip_idx] <- 1 - merged$freq_aligned[flip_idx]

# ===== MAF filter =====
merged$MAF <- pmin(merged$freq_aligned, 1 - merged$freq_aligned)
merged <- merged[merged$MAF >= 0.01, ]

# ===== Define final SNP set and subset/align LD and GWAS =====
# SNPs present after allele QC & MAF filter
snps_gwas <- merged$rsID

# SNPs present in cleaned LD
snps_ld <- rownames(ld_matrix_clean)

# Intersection
final_snps <- intersect(snps_ld, snps_gwas)

length(snps_ld) # number of SNPs in LD
length(snps_gwas) # number of SNPs in GWAS after MAF filter
length(final_snps) # number of SNPs in final overlap

# Subset LD to final SNPs
ld_matrix_final <- ld_matrix_clean[final_snps, final_snps]

# Check eigenvalues (optional)
eig_vals <- eigen(ld_matrix_final, symmetric = TRUE, only.values = TRUE)$values
min(eig_vals)
summary(eig_vals)

# Make LD matrix numerically positive-definite (nearPD)
ld_pd <- as.matrix(nearPD(ld_matrix_final, corr = TRUE)$mat)

# Check eigenvalues again (optional)
eig_vals_pd <- eigen(ld_pd, symmetric = TRUE, only.values = TRUE)$values
min(eig_vals_pd)
summary(eig_vals_pd)

# Subset merged GWAS to final SNPs and reorder to match LD ordering
merged_final <- merged[merged$rsID %in% final_snps, ]
merged_final <- merged_final[match(final_snps, merged_final$rsID), ]

# Sanity checks
stopifnot(all(rownames(ld_pd) == final_snps))
stopifnot(all(merged_final$rsID == final_snps))

# ===== Diagnostic plot 2: lambda_hat + kriging_rss =====
z_diag <- merged_final$beta_aligned / merged_final$se

# Estimate lambda (residual variance scaling)
lambda_hat <- estimate_s_rss(z_diag, ld_pd, n = ntotal)
lambda_hat 

# Kriging-based conditional z-scores and diagnostic plot
condz <- kriging_rss(z_diag, ld_pd, n = ntotal)

png(filename = "Finn_cndzplot.png", width = 3000, height = 1500, res = 300)
print(condz$plot)   # condz$plot is a ggplot object
dev.off()

# z_obs <- condz[["conditional_dist"]][["z"]]
# z_exp <- condz[["conditional_dist"]][["condmean"]]
# diff_val <- abs(z_obs - z_exp)
# bad_idx <- which(diff_val > 2.5) # 3 or 2.5
# bad_idx
# snp_names <- rownames(ld_pd)
# bad_snps_condz <- snp_names[bad_idx]
# bad_snps_condz
# Add to bad_snps and run again

# Keep only those SNPs in LD as well
keep_snps <- merged_final$rsID

# Subset LD to the SNP set
R <- ld_pd[keep_snps, keep_snps]

# Reorder merged_final to match R ordering (for safety)
merged_final <- merged_final[match(keep_snps, merged_final$rsID), ]

# Sanity checks
stopifnot(all(rownames(R) == merged_final$rsID))

# ===== Build z, beta, varbeta and dataset1 for coloc/SuSiE =====
beta_aligned <- merged_final$beta_aligned
se_aligned <- merged_final$se
varbeta <- se_aligned^2
z <- beta_aligned / se_aligned
pos <- merged_final$pos
snp <- merged_final$rsID

# Named vectors
named_beta <- beta_aligned %>% set_names(snp)
named_varbeta <- varbeta %>% set_names(snp)

# Build dataset1 for coloc (summary + LD)
dataset1 <- list(
  beta = named_beta,
  varbeta = named_varbeta,
  s = s,
  type = "cc",
  snp = snp,
  LD = R,
  position = pos,
  N = ntotal
)

# Sanity check (coloc)
coloc:::check_dataset(dataset1, 1)

# ===== Run SuSiE with LD via coloc::runsusie =====
S3 <- runsusie(dataset1,
               nref = 503,
               prior_variance = 0.2^2,
               estimate_prior_variance = FALSE,
               coverage = 0.90)

# ===== SuSiE PIP plot =====
png(filename = "UKB_SuSiEPlot.png", width = 3000, height = 1500, res = 300)
susie_plot(S3, y = "PIP")

# Optionally, highlight the lead SNP (by PIP index)
pips <- S3$pip
lead_pip <- which.max(pips)

points(
  x = lead_pip,
  y = pips[lead_pip],
  col = 2,
  pch = 16
)
dev.off()

# ===== Summarize SuSiE credible sets =====
SNPs <- lapply(S3$sets$cs, function(x) data.frame(SNP = names(x), variable = x))
SNPs_number <- do.call(rbind, SNPs)

probab <- summary(S3)$vars
susie_result <- inner_join(probab, SNPs_number, by = "variable") %>% 
  dplyr::select(SNP, variable_prob, cs) %>% 
  dplyr::mutate(LD_reference = "1000G", .before = 1)

print(susie_result)
write.table(x = susie_result, file = "UKB_susie_res.txt", sep = "\t")
