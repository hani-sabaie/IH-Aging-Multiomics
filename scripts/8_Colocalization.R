# ===== Clean environment =====
rm(list = ls(all.names = TRUE))
gc()

# ===== Loading relevant libraries =====
library(data.table)
library(dplyr)
library(coloc)
library(locuscomparer)
library(purrr)

# ============================================================================ #

# ===== COLOC =====
# eQTL data process
# SNPs in SMAD3 gene region ± 1000kb
eqtl <- fread("myLiteCisEqtl.txt", data.table = F)
eqtl_bin1 <- eqtl %>%
  filter(
    BP > (67357940 - 1e6) &
      BP < (67487507 + 1e6)
  )

# GWAS summary
gwas <- fread("UKB_gwas_for_smr.txt", data.table = F)

# Prepare input data for coloc
gwas <- gwas %>%
  mutate(
    MAF_gwas = ifelse(freq <= 0.5, freq, 1 - freq)
  )

eqtl_sub <- eqtl_bin1 %>%
  transmute(
    SNP,
    A1_eqtl = A1,
    A2_eqtl = A2,
    beta_eqtl = b,
    se_eqtl = SE,
    p_eqtl = p,
    freq_eqtl = Freq,
    gene = Gene
  )

gwas_sub <- gwas %>%
  transmute(
    SNP,
    A1_gwas = A1,
    A2_gwas = A2,
    beta_gwas = b,
    se_gwas = se,
    p_gwas = p,
    MAF_gwas
  )

df_coloc <- inner_join(eqtl_sub, gwas_sub, by = "SNP")
same_alleles <- (df_coloc$A1_eqtl == df_coloc$A1_gwas & df_coloc$A2_eqtl == df_coloc$A2_gwas) |
  (df_coloc$A1_eqtl == df_coloc$A2_gwas & df_coloc$A2_eqtl == df_coloc$A1_gwas)

df_coloc <- df_coloc[same_alleles, ]

flip <- (df_coloc$A1_eqtl == df_coloc$A2_gwas & df_coloc$A2_eqtl == df_coloc$A1_gwas)
df_coloc$beta_gwas[flip] <- -df_coloc$beta_gwas[flip]

# Run coloc 
N_all_gwas <- 371810   
N_cases <- 28707    
N_all_eqtl <- 663       

# Case/control GWAS design
genes <- unique(df_coloc$gene)

results <- lapply(genes, function(g) {
  tmp <- df_coloc %>%
    filter(gene == g) %>%     
    arrange(SNP) %>%
    distinct(SNP, .keep_all = TRUE) %>%
    filter(p_eqtl > 0, p_gwas > 0)
  
  if (nrow(tmp) < 1) return(NULL)
  
  D1 <- list(
    snp = tmp$SNP,
    beta = tmp$beta_gwas,
    varbeta = (tmp$se_gwas)^2,
    type = "cc",
    s = N_cases / N_all_gwas, 
    N = N_all_gwas
  )
  
  D2 <- list(
    snp = tmp$SNP,
    beta = tmp$beta_eqtl,
    varbeta = (tmp$se_eqtl)^2,
    type = "quant",
    N = N_all_eqtl
  )
  
  coloc.abf(dataset1 = D1,
            dataset2 = D2,
            MAF = tmp$MAF_gwas)
})
names(results) <- genes
coloc_sum <- map_dfr(results, function(res) {
  if (is.null(res)) return(NULL)
  
  snp_tab <- res$results
  snp_tab <- snp_tab[order(-snp_tab$SNP.PP.H4), ]
  
  data.frame(
    nsnps = res$summary["nsnps"],
    pph0 = res$summary["PP.H0.abf"],
    pph1 = res$summary["PP.H1.abf"],
    pph2 = res$summary["PP.H2.abf"],
    pph3 = res$summary["PP.H3.abf"],
    pph4 = res$summary["PP.H4.abf"],
    top_snp = snp_tab$snp[1],
    top_snp_pph4 = snp_tab$SNP.PP.H4[1]
  )
}, .id = "gene")

hits_pph4 <- coloc_sum %>%
  filter(pph4 > 0.5)

write.csv(x = hits_pph4, file = "ColocResults.csv", sep = ",")

