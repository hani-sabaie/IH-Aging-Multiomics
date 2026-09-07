# Reviewer C7
# SMR discovery-replication multiple-testing audit
#
# Evaluates:
# 1) BH separately within each UKB tissue
# 2) BH across all UKB gene-tissue tests jointly
# 3) Integration with FAP DEGs and brown hdWGCNA module
# 4) Targeted replication of discovery hypotheses in FinnGen
#
# No original result files are modified.

suppressPackageStartupMessages({
  library(data.table)
})

smr_dir <- "processed_results/06_SMR_HEIDI"
audit_dir <- file.path(smr_dir, "multiple_testing_audit")

deg_file <- "processed_results/02_differential_expression/bulk_de_sig_faps.csv"
module_file <- "processed_results/05_hdWGCNA/module_assignment_table.csv"

out_dir <- file.path(audit_dir, "discovery_replication")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Load FAP DEG and brown-module genes
# ------------------------------------------------------------

deg <- fread(deg_file)
modules <- fread(module_file)

deg_genes <- unique(deg$gene)
brown_genes <- unique(
  modules[tolower(module) == "brown", gene_name]
)

deg_brown <- intersect(deg_genes, brown_genes)

cat("\n============================================================\n")
cat("FAP DEG / BROWN MODULE INTEGRATION\n")
cat("============================================================\n")
cat("FAP DEG genes:", length(deg_genes), "\n")
cat("Brown-module genes:", length(brown_genes), "\n")
cat("FAP DEG ∩ Brown:", length(deg_brown), "\n\n")


# ------------------------------------------------------------
# Function to read all four files for one cohort
# ------------------------------------------------------------

read_cohort <- function(cohort) {

  files <- list.files(
    smr_dir,
    pattern = paste0("^", cohort, "_.*\\.smr$"),
    full.names = TRUE
  )

  if (length(files) != 4) {
    stop("Expected 4 SMR files for ", cohort, ", found ", length(files))
  }

  ans <- lapply(files, function(f) {

    x <- fread(f)

    tissue <- sub(
      paste0("^", cohort, "_"),
      "",
      sub("\\.smr$", "", basename(f))
    )

    x[, cohort := cohort]
    x[, tissue := tissue]
    x[, source_file := basename(f)]

    x
  })

  rbindlist(ans, fill = TRUE)
}


ukb <- read_cohort("UKB")
fin <- read_cohort("FinnGen")


# ------------------------------------------------------------
# Validate p values
# ------------------------------------------------------------

ukb[, valid_p := is.finite(p_SMR) & p_SMR >= 0 & p_SMR <= 1]
fin[, valid_p := is.finite(p_SMR) & p_SMR >= 0 & p_SMR <= 1]


# ------------------------------------------------------------
# A. BH within each tissue
# ------------------------------------------------------------

ukb[, FDR_within_tissue := NA_real_]
ukb[
  valid_p == TRUE,
  FDR_within_tissue := p.adjust(p_SMR, method = "BH"),
  by = tissue
]

fin[, FDR_within_tissue := NA_real_]
fin[
  valid_p == TRUE,
  FDR_within_tissue := p.adjust(p_SMR, method = "BH"),
  by = tissue
]


# ------------------------------------------------------------
# B. BH across all four tissues within cohort
# ------------------------------------------------------------

ukb[, FDR_cohortwide := NA_real_]
ukb[
  valid_p == TRUE,
  FDR_cohortwide := p.adjust(p_SMR, method = "BH")
]

fin[, FDR_cohortwide := NA_real_]
fin[
  valid_p == TRUE,
  FDR_cohortwide := p.adjust(p_SMR, method = "BH")
]


# ------------------------------------------------------------
# HEIDI-compatible shared-signal definitions
# ------------------------------------------------------------

ukb[, pass_per_tissue :=
      !is.na(FDR_within_tissue) &
      FDR_within_tissue < 0.05 &
      !is.na(p_HEIDI) &
      p_HEIDI > 0.01]

ukb[, pass_cohortwide :=
      !is.na(FDR_cohortwide) &
      FDR_cohortwide < 0.05 &
      !is.na(p_HEIDI) &
      p_HEIDI > 0.01]

fin[, pass_per_tissue :=
      !is.na(FDR_within_tissue) &
      FDR_within_tissue < 0.05 &
      !is.na(p_HEIDI) &
      p_HEIDI > 0.01]

fin[, pass_cohortwide :=
      !is.na(FDR_cohortwide) &
      FDR_cohortwide < 0.05 &
      !is.na(p_HEIDI) &
      p_HEIDI > 0.01]


# ------------------------------------------------------------
# SMAD3 audit
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SMAD3: UKB DISCOVERY\n")
cat("============================================================\n\n")

print(
  ukb[
    Gene == "SMAD3",
    .(
      tissue,
      probeID,
      topSNP,
      b_SMR,
      se_SMR,
      p_SMR,
      FDR_within_tissue,
      FDR_cohortwide,
      p_HEIDI,
      pass_per_tissue,
      pass_cohortwide
    )
  ]
)


cat("\n============================================================\n")
cat("SMAD3: FINNGEN\n")
cat("============================================================\n\n")

print(
  fin[
    Gene == "SMAD3",
    .(
      tissue,
      probeID,
      topSNP,
      b_SMR,
      se_SMR,
      p_SMR,
      FDR_within_tissue,
      FDR_cohortwide,
      p_HEIDI,
      pass_per_tissue,
      pass_cohortwide
    )
  ]
)


# ------------------------------------------------------------
# Discovery candidate genes:
# DEG ∩ Brown ∩ corrected UKB SMR
# ------------------------------------------------------------

ukb_per_tissue_genes <- unique(
  ukb[pass_per_tissue == TRUE, Gene]
)

ukb_cohortwide_genes <- unique(
  ukb[pass_cohortwide == TRUE, Gene]
)

cand_per_tissue <- intersect(
  deg_brown,
  ukb_per_tissue_genes
)

cand_cohortwide <- intersect(
  deg_brown,
  ukb_cohortwide_genes
)


cat("\n============================================================\n")
cat("UKB DISCOVERY CANDIDATES\n")
cat("============================================================\n")

cat("\nUsing BH separately within each tissue:\n")
cat("N =", length(cand_per_tissue), "\n")
print(cand_per_tissue)

cat("\nUsing BH across all four UKB tissues:\n")
cat("N =", length(cand_cohortwide), "\n")
print(cand_cohortwide)


# ------------------------------------------------------------
# Discovery gene-tissue hypotheses
# ------------------------------------------------------------

disc_per_tissue <- ukb[
  pass_per_tissue == TRUE &
  Gene %in% deg_brown,
  .(
    Gene,
    tissue,
    probeID,
    topSNP,
    b_SMR_UKB = b_SMR,
    se_SMR_UKB = se_SMR,
    p_SMR_UKB = p_SMR,
    FDR_UKB = FDR_within_tissue,
    p_HEIDI_UKB = p_HEIDI
  )
]

disc_cohortwide <- ukb[
  pass_cohortwide == TRUE &
  Gene %in% deg_brown,
  .(
    Gene,
    tissue,
    probeID,
    topSNP,
    b_SMR_UKB = b_SMR,
    se_SMR_UKB = se_SMR,
    p_SMR_UKB = p_SMR,
    FDR_UKB = FDR_cohortwide,
    p_HEIDI_UKB = p_HEIDI
  )
]


# ------------------------------------------------------------
# Targeted FinnGen replication
#
# Match exact discovery Gene + tissue hypothesis.
# Correct replication P values only across discovery hypotheses.
# ------------------------------------------------------------

replicate_candidates <- function(discovery, label) {

  if (nrow(discovery) == 0) {
    cat("\n", label, ": no discovery hypotheses to replicate.\n", sep = "")
    return(data.table())
  }

  rep <- merge(
    discovery,
    fin[
      ,
      .(
        Gene,
        tissue,
        topSNP_FinnGen = topSNP,
        b_SMR_FinnGen = b_SMR,
        se_SMR_FinnGen = se_SMR,
        p_SMR_FinnGen = p_SMR,
        p_HEIDI_FinnGen = p_HEIDI
      )
    ],
    by = c("Gene", "tissue"),
    all.x = TRUE
  )

  rep[, replication_BH := p.adjust(
    p_SMR_FinnGen,
    method = "BH"
  )]

  rep[, replication_nominal :=
        !is.na(p_SMR_FinnGen) &
        p_SMR_FinnGen < 0.05 &
        !is.na(p_HEIDI_FinnGen) &
        p_HEIDI_FinnGen > 0.01]

  rep[, replication_BH_pass :=
        !is.na(replication_BH) &
        replication_BH < 0.05 &
        !is.na(p_HEIDI_FinnGen) &
        p_HEIDI_FinnGen > 0.01]

  cat("\n============================================================\n")
  cat(label, "\n")
  cat("============================================================\n\n")

  print(rep)

  rep
}


rep_per_tissue <- replicate_candidates(
  disc_per_tissue,
  "FINNGEN TARGETED REPLICATION: PER-TISSUE UKB BH DISCOVERY"
)

rep_cohortwide <- replicate_candidates(
  disc_cohortwide,
  "FINNGEN TARGETED REPLICATION: COHORT-WIDE UKB BH DISCOVERY"
)


# ------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------

fwrite(
  ukb,
  file.path(out_dir, "UKB_all_tissues_with_BH.tsv"),
  sep = "\t"
)

fwrite(
  fin,
  file.path(out_dir, "FinnGen_all_tissues_with_BH.tsv"),
  sep = "\t"
)

fwrite(
  disc_per_tissue,
  file.path(out_dir, "UKB_discovery_per_tissue_BH.csv")
)

fwrite(
  disc_cohortwide,
  file.path(out_dir, "UKB_discovery_cohortwide_BH.csv")
)

fwrite(
  rep_per_tissue,
  file.path(out_dir, "FinnGen_replication_per_tissue_discovery.csv")
)

fwrite(
  rep_cohortwide,
  file.path(out_dir, "FinnGen_replication_cohortwide_discovery.csv")
)

cat("\nOutputs written to:\n")
cat(out_dir, "\n")
