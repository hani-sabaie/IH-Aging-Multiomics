# ============================================================================
# Reviewer C7
# Human pseudobulk DE multiple-testing and integration audit
#
# No canonical files are overwritten.
# ============================================================================

rm(list = ls(all.names = TRUE))
gc()

suppressPackageStartupMessages({
  library(data.table)
})

# -------------------------------------------------------------------------
# Paths
# -------------------------------------------------------------------------

de_file <- file.path(
  "processed_results",
  "02_differential_expression",
  "full_precision_audit",
  "all_celltypes_full_precision.tsv"
)

module_file <- file.path(
  "processed_results",
  "05_hdWGCNA",
  "module_assignment_table.csv"
)

ukb_file <- file.path(
  "processed_results",
  "06_SMR_HEIDI",
  "multiple_testing_audit",
  "discovery_replication",
  "UKB_all_tissues_with_BH.tsv"
)

fin_file <- file.path(
  "processed_results",
  "06_SMR_HEIDI",
  "multiple_testing_audit",
  "discovery_replication",
  "FinnGen_all_tissues_with_BH.tsv"
)

outdir <- file.path(
  "processed_results",
  "02_differential_expression",
  "full_precision_audit",
  "integration_audit"
)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# Load
# -------------------------------------------------------------------------

de <- fread(de_file)
modules <- fread(module_file)
ukb <- fread(ukb_file)
fin <- fread(fin_file)

# -------------------------------------------------------------------------
# Cell-type adequacy
#
# EC3 and FAP4 showed inadequate donor representation across a broad range
# of nuclei-count thresholds (5-20 nuclei).
#
# IMPORTANT:
# We exclude these cell TYPES from formal DE rather than selectively
# excluding individual donor pseudobulks.
# -------------------------------------------------------------------------

excluded_celltypes <- c("EC3", "FAP4")

de_eligible <- de[
  !celltype %in% excluded_celltypes
]

fap_types <- c("FAP1", "FAP2", "FAP3")

fap <- de_eligible[
  celltype %in% fap_types
]

# -------------------------------------------------------------------------
# Multiple-testing frameworks
# -------------------------------------------------------------------------

# Existing BH correction performed separately within each cell type
de_eligible[, FDR_within_celltype := p_val_adj]

# BH across all eligible gene x cell-type tests
de_eligible[, FDR_global_eligible := p.adjust(p_val, method = "BH")]

# FAP-specific table with BH across all FAP1-FAP3 tests
fap <- copy(de_eligible[celltype %in% fap_types])

fap[, FDR_FAPwide := p.adjust(p_val, method = "BH")]

# bring global FDR already calculated above
# FDR_within_celltype and FDR_global_eligible are already present

# -------------------------------------------------------------------------
# Significance definitions
# -------------------------------------------------------------------------

de_eligible[, pass_within_abs :=
  FDR_within_celltype < 0.05 &
  abs(avg_log2FC) >= 1
]

de_eligible[, pass_global_abs :=
  FDR_global_eligible < 0.05 &
  abs(avg_log2FC) >= 1
]

fap[, pass_within_abs :=
  FDR_within_celltype < 0.05 &
  abs(avg_log2FC) >= 1
]

fap[, pass_FAPwide_abs :=
  FDR_FAPwide < 0.05 &
  abs(avg_log2FC) >= 1
]

fap[, pass_global_abs :=
  FDR_global_eligible < 0.05 &
  abs(avg_log2FC) >= 1
]

# FDR-only sensitivity: no hard fold-change threshold
fap[, pass_within_FDRonly :=
  FDR_within_celltype < 0.05
]

fap[, pass_FAPwide_FDRonly :=
  FDR_FAPwide < 0.05
]

# -------------------------------------------------------------------------
# Count summaries
# -------------------------------------------------------------------------

cat("\n============================================================\n")
cat("ELIGIBLE CELL TYPES\n")
cat("============================================================\n\n")

cat("Excluded from formal pseudobulk DE:\n")
print(excluded_celltypes)

cat("\nIncluded:\n")
print(sort(unique(de_eligible$celltype)))

cat("\n============================================================\n")
cat("DEG COUNTS: ELIGIBLE CELL TYPES\n")
cat("============================================================\n\n")

all_summary <- de_eligible[
  ,
  .(
    genes_tested = .N,
    within_BH_absFC = sum(pass_within_abs),
    global_BH_absFC = sum(pass_global_abs),
    within_up = sum(pass_within_abs & avg_log2FC >= 1),
    within_down = sum(pass_within_abs & avg_log2FC <= -1)
  ),
  by = celltype
]

print(all_summary)

cat("\nTOTAL:\n")
print(
  de_eligible[
    ,
    .(
      within_BH_absFC = sum(pass_within_abs),
      global_BH_absFC = sum(pass_global_abs),
      within_up = sum(pass_within_abs & avg_log2FC >= 1),
      within_down = sum(pass_within_abs & avg_log2FC <= -1)
    )
  ]
)

cat("\n============================================================\n")
cat("FAP1-FAP3 DEG COUNTS\n")
cat("============================================================\n\n")

fap_summary <- fap[
  ,
  .(
    genes_tested = .N,
    within_BH_absFC = sum(pass_within_abs),
    FAPwide_BH_absFC = sum(pass_FAPwide_abs),
    global_BH_absFC = sum(pass_global_abs),
    within_FDRonly = sum(pass_within_FDRonly),
    FAPwide_FDRonly = sum(pass_FAPwide_FDRonly)
  ),
  by = celltype
]

print(fap_summary)

cat("\nTOTAL FAP1-FAP3:\n")
print(
  fap[
    ,
    .(
      within_BH_absFC = sum(pass_within_abs),
      FAPwide_BH_absFC = sum(pass_FAPwide_abs),
      global_BH_absFC = sum(pass_global_abs),
      within_FDRonly = sum(pass_within_FDRonly),
      FAPwide_FDRonly = sum(pass_FAPwide_FDRonly)
    )
  ]
)

# -------------------------------------------------------------------------
# Candidate audit
# -------------------------------------------------------------------------

cat("\n============================================================\n")
cat("SMAD3 / PLEKHA6 UNDER ALL CORRECTION FRAMEWORKS\n")
cat("============================================================\n\n")

candidate <- fap[
  gene %in% c("SMAD3", "PLEKHA6"),
  .(
    gene,
    celltype,
    avg_log2FC,
    p_val,
    FDR_within_celltype,
    FDR_FAPwide,
    FDR_global_eligible,
    pass_within_abs,
    pass_FAPwide_abs,
    pass_global_abs,
    pass_within_FDRonly,
    pass_FAPwide_FDRonly
  )
]

print(candidate)

# -------------------------------------------------------------------------
# Module sets
# -------------------------------------------------------------------------

yellow <- unique(
  modules[tolower(module) == "yellow", gene_name]
)

brown <- unique(
  modules[tolower(module) == "brown", gene_name]
)

blue <- unique(
  modules[tolower(module) == "blue", gene_name]
)

module_sets <- list(
  Yellow = yellow,
  Brown = brown,
  Blue = blue
)

# -------------------------------------------------------------------------
# Corrected UKB discovery set
# -------------------------------------------------------------------------

ukb_sig <- ukb[
  !is.na(FDR_cohortwide) &
  FDR_cohortwide < 0.05 &
  !is.na(p_HEIDI) &
  p_HEIDI > 0.01
]

ukb_genes <- unique(ukb_sig$Gene)

# -------------------------------------------------------------------------
# Integration function
# -------------------------------------------------------------------------

audit_framework <- function(label, deg_genes) {

  cat("\n------------------------------------------------------------\n")
  cat(label, "\n")
  cat("------------------------------------------------------------\n")

  cat("Unique FAP DEG genes:", length(deg_genes), "\n")

  candidates <- character()

  for (m in names(module_sets)) {

    z <- sort(
      intersect(
        intersect(
          deg_genes,
          module_sets[[m]]
        ),
        ukb_genes
      )
    )

    cat(
      m,
      "∩ FAP-DEG ∩ UKB-SMR-BH:",
      length(z),
      "\n"
    )

    if (length(z)) print(z)

    candidates <- union(candidates, z)
  }

  candidates <- sort(unique(candidates))

  cat("\nAll candidate genes:\n")
  print(candidates)

  if (length(candidates) == 0) {
    return(NULL)
  }

  # Exact discovery gene-tissue hypotheses
  disc <- ukb_sig[
    Gene %in% candidates,
    .(
      Gene,
      tissue,
      probeID,
      topSNP_UKB = topSNP,
      b_SMR_UKB = b_SMR,
      p_SMR_UKB = p_SMR,
      FDR_UKB = FDR_cohortwide,
      p_HEIDI_UKB = p_HEIDI
    )
  ]

  # FinnGen exact gene+tissue replication
  fin_sub <- fin[
    ,
    .(
      Gene,
      tissue,
      topSNP_FinnGen = topSNP,
      b_SMR_FinnGen = b_SMR,
      p_SMR_FinnGen = p_SMR,
      p_HEIDI_FinnGen = p_HEIDI
    )
  ]

  rep <- merge(
    disc,
    fin_sub,
    by = c("Gene", "tissue"),
    all.x = TRUE
  )

  rep[, replication_BH :=
    p.adjust(p_SMR_FinnGen, method = "BH")
  ]

  rep[, replication_pass :=
    !is.na(replication_BH) &
    replication_BH < 0.05 &
    !is.na(p_HEIDI_FinnGen) &
    p_HEIDI_FinnGen > 0.01
  ]

  cat("\nExact FinnGen targeted replication:\n")
  print(rep)

  rep[, framework := label]

  rep
}

# -------------------------------------------------------------------------
# Run frameworks
# -------------------------------------------------------------------------

deg_within_abs <- unique(
  fap[pass_within_abs == TRUE, gene]
)

deg_fapwide_abs <- unique(
  fap[pass_FAPwide_abs == TRUE, gene]
)

deg_global_abs <- unique(
  fap[pass_global_abs == TRUE, gene]
)

deg_within_fdr <- unique(
  fap[pass_within_FDRonly == TRUE, gene]
)

deg_fapwide_fdr <- unique(
  fap[pass_FAPwide_FDRonly == TRUE, gene]
)

cat("\n============================================================\n")
cat("MODULE + UKB SMR INTEGRATION\n")
cat("============================================================\n")

rep1 <- audit_framework(
  "A: Within-celltype BH + |log2FC| >= 1",
  deg_within_abs
)

rep2 <- audit_framework(
  "B: FAP1-FAP3-wide BH + |log2FC| >= 1",
  deg_fapwide_abs
)

rep3 <- audit_framework(
  "C: Global eligible-test BH + |log2FC| >= 1",
  deg_global_abs
)

rep4 <- audit_framework(
  "D: Within-celltype BH, FDR-only sensitivity",
  deg_within_fdr
)

rep5 <- audit_framework(
  "E: FAP1-FAP3-wide BH, FDR-only sensitivity",
  deg_fapwide_fdr
)

# -------------------------------------------------------------------------
# Outputs
# -------------------------------------------------------------------------

fwrite(
  fap,
  file.path(outdir, "FAP1_FAP3_multiple_testing_audit.tsv"),
  sep = "\t"
)

fwrite(
  candidate,
  file.path(outdir, "SMAD3_PLEKHA6_correction_frameworks.csv")
)

reps <- Filter(Negate(is.null), list(rep1, rep2, rep3, rep4, rep5))

if (length(reps)) {
  fwrite(
    rbindlist(reps, fill = TRUE),
    file.path(outdir, "integration_FinnGen_replication_all_frameworks.csv")
  )
}

cat("\nOutputs written to:\n")
cat(outdir, "\n")
