# ============================================================================
# Reviewer C7 pre-production validation
#
# Validate that the canonical combined SMR resource reproduces the corrected
# UKB discovery set and revised FAP/module integration.
#
# No canonical result or figure is modified.
# ============================================================================

rm(list = ls(all.names = TRUE))
gc()

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (length(file_arg) == 1L) {
  script_dir <- dirname(
    normalizePath(
      sub("^--file=", "", file_arg)
    )
  )
  repo_root <- normalizePath(
    file.path(script_dir, "..")
  )
} else {
  repo_root <- normalizePath(".")
}

smr_file <- file.path(
  repo_root,
  "processed_results",
  "06_SMR_HEIDI",
  "SMR_all_studies_all_tissues.csv"
)

deg_file <- file.path(
  repo_root,
  "processed_results",
  "02_differential_expression",
  "bulk_de_sig_faps.csv"
)

module_file <- file.path(
  repo_root,
  "processed_results",
  "05_hdWGCNA",
  "module_assignment_table.csv"
)

audit_file <- file.path(
  repo_root,
  "processed_results",
  "06_SMR_HEIDI",
  "multiple_testing_audit",
  "figure10_revised_candidate",
  "Figure10A_UKB_corrected_discovery_rows.tsv"
)

for (f in c(
  smr_file,
  deg_file,
  module_file,
  audit_file
)) {
  if (!file.exists(f)) {
    stop("Required file missing: ", f)
  }
}

# --------------------------------------------------------------------------
# Canonical resources
# --------------------------------------------------------------------------
smr <- fread(smr_file)
deg <- fread(deg_file)
modules <- fread(module_file)

required_smr <- c(
  "study",
  "tissue",
  "probeID",
  "Gene",
  "topSNP",
  "p_SMR",
  "p_HEIDI"
)

miss <- setdiff(
  required_smr,
  names(smr)
)

if (length(miss)) {
  stop(
    "Canonical SMR table missing: ",
    paste(miss, collapse = ", ")
  )
}

# --------------------------------------------------------------------------
# UKB cohort-wide BH directly from canonical combined SMR table
# --------------------------------------------------------------------------
ukb <- copy(
  smr[study == "UKB"]
)

ukb[
  ,
  tissue_normalized :=
    gsub(
      " ",
      "_",
      tissue,
      fixed = TRUE
    )
]

ukb[
  ,
  valid_p :=
    is.finite(p_SMR) &
    p_SMR >= 0 &
    p_SMR <= 1
]

ukb[
  ,
  FDR_cohortwide :=
    NA_real_
]

ukb[
  valid_p == TRUE,
  FDR_cohortwide :=
    p.adjust(
      p_SMR,
      method = "BH"
    )
]

ukb[
  ,
  pass_cohortwide :=
    !is.na(FDR_cohortwide) &
    FDR_cohortwide < 0.05 &
    !is.na(p_HEIDI) &
    p_HEIDI > 0.01
]

discovery <- ukb[
  pass_cohortwide == TRUE
]

discovery_genes <- sort(
  unique(
    na.omit(
      discovery$Gene
    )
  )
)

# --------------------------------------------------------------------------
# Compare with previously validated audit set
# --------------------------------------------------------------------------
audit <- fread(audit_file)

canonical_cmp <- discovery[
  ,
  .(
    probeID,
    Gene,
    topSNP,
    tissue = tissue_normalized,
    p_SMR,
    FDR_cohortwide,
    p_HEIDI
  )
]

audit_cmp <- audit[
  ,
  .(
    probeID,
    Gene,
    topSNP,
    tissue,
    p_SMR,
    FDR_cohortwide,
    p_HEIDI
  )
]

key_cols <- c(
  "probeID",
  "Gene",
  "topSNP",
  "tissue"
)

canonical_keys <- unique(
  canonical_cmp[
    ,
    ..key_cols
  ]
)

audit_keys <- unique(
  audit_cmp[
    ,
    ..key_cols
  ]
)

canonical_only <- fsetdiff(
  canonical_keys,
  audit_keys
)

audit_only <- fsetdiff(
  audit_keys,
  canonical_keys
)

merged <- merge(
  canonical_cmp,
  audit_cmp,
  by = key_cols,
  suffixes = c(
    "_canonical",
    "_audit"
  )
)

max_delta_p <- max(
  abs(
    merged$p_SMR_canonical -
      merged$p_SMR_audit
  ),
  na.rm = TRUE
)

max_delta_fdr <- max(
  abs(
    merged$FDR_cohortwide_canonical -
      merged$FDR_cohortwide_audit
  ),
  na.rm = TRUE
)

max_delta_heidi <- max(
  abs(
    merged$p_HEIDI_canonical -
      merged$p_HEIDI_audit
  ),
  na.rm = TRUE
)

# --------------------------------------------------------------------------
# Revised FAP/module integration
# --------------------------------------------------------------------------
deg_genes <- sort(
  unique(
    na.omit(
      deg$gene
    )
  )
)

module_gene_set <- function(m) {
  sort(
    unique(
      na.omit(
        modules[
          tolower(module) == tolower(m),
          gene_name
        ]
      )
    )
  )
}

yellow <- module_gene_set("Yellow")
brown  <- module_gene_set("Brown")
blue   <- module_gene_set("Blue")

yellow_hits <- Reduce(
  intersect,
  list(
    yellow,
    deg_genes,
    discovery_genes
  )
)

brown_hits <- Reduce(
  intersect,
  list(
    brown,
    deg_genes,
    discovery_genes
  )
)

blue_hits <- Reduce(
  intersect,
  list(
    blue,
    deg_genes,
    discovery_genes
  )
)

# --------------------------------------------------------------------------
# Critical genes
# --------------------------------------------------------------------------
critical <- discovery[
  Gene %in% c(
    "SMAD3",
    "PLEKHA6"
  ),
  .(
    Gene,
    tissue,
    probeID,
    topSNP,
    p_SMR,
    FDR_cohortwide,
    p_HEIDI
  )
][
  order(Gene, tissue)
]

# --------------------------------------------------------------------------
# Validation
# --------------------------------------------------------------------------
pass <-
  nrow(discovery) == 138L &&
  length(discovery_genes) == 84L &&
  length(deg_genes) == 194L &&
  nrow(canonical_only) == 0L &&
  nrow(audit_only) == 0L &&
  max_delta_p < 1e-12 &&
  max_delta_fdr < 1e-12 &&
  max_delta_heidi < 1e-12 &&
  identical(
    sort(brown_hits),
    "SMAD3"
  ) &&
  identical(
    sort(blue_hits),
    "PLEKHA6"
  ) &&
  length(yellow_hits) == 0L

cat("\n============================================================\n")
cat("CANONICAL UKB DISCOVERY VALIDATION\n")
cat("============================================================\n\n")

cat("UKB rows:", nrow(ukb), "\n")
cat("Corrected discovery rows:", nrow(discovery), "\n")
cat("Corrected discovery unique genes:", length(discovery_genes), "\n")
cat("FAP unique genes:", length(deg_genes), "\n\n")

cat("Canonical-only discovery keys:", nrow(canonical_only), "\n")
cat("Audit-only discovery keys:", nrow(audit_only), "\n")

cat(
  "Max |delta p_SMR|:",
  format(max_delta_p, scientific = TRUE),
  "\n"
)

cat(
  "Max |delta FDR|:",
  format(max_delta_fdr, scientific = TRUE),
  "\n"
)

cat(
  "Max |delta HEIDI|:",
  format(max_delta_heidi, scientific = TRUE),
  "\n\n"
)

cat("Yellow ∩ FAP ∩ UKB:\n")
print(yellow_hits)

cat("\nBrown ∩ FAP ∩ UKB:\n")
print(brown_hits)

cat("\nBlue ∩ FAP ∩ UKB:\n")
print(blue_hits)

cat("\nCritical corrected SMR rows:\n")
print(critical)

cat("\n============================================================\n")
cat("FINAL STATUS\n")
cat("============================================================\n\n")

cat(
  "Canonical SMR integration validation pass:",
  pass,
  "\n"
)

cat("\nNo canonical file was modified.\n")

if (!pass) {
  stop(
    "Canonical SMR integration validation failed. ",
    "Do not patch production integration."
  )
}
