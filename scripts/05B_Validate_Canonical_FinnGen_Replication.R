# ============================================================================
# Reviewer C7
# Validate targeted FinnGen replication directly from canonical resources.
#
# Discovery candidates:
#   Brown -> SMAD3
#   Blue  -> PLEKHA6
#
# FinnGen replication:
#   exact gene + tissue + probe carried forward from UKB discovery
#   BH correction across the two carried-forward hypotheses
#   replication requires BH < 0.05 and HEIDI P > 0.01
#
# No canonical file is modified.
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

candidate_file <- file.path(
  repo_root,
  "processed_results",
  "06_SMR_HEIDI",
  "integration_FAP_hdWGCNA_UKB_candidate_details.csv"
)

for (f in c(smr_file, candidate_file)) {
  if (!file.exists(f)) {
    stop("Required canonical file missing: ", f)
  }
}

smr <- fread(smr_file)
cand <- fread(candidate_file)

# --------------------------------------------------------------------------
# Validate candidate set
# --------------------------------------------------------------------------
if (
  nrow(cand) != 2L ||
  !setequal(cand$Gene, c("SMAD3", "PLEKHA6"))
) {
  stop(
    "Expected exactly two UKB discovery candidates: SMAD3 and PLEKHA6."
  )
}

# Normalize tissue naming for exact cross-cohort matching.
normalize_tissue <- function(x) {
  gsub(
    " ",
    "_",
    x,
    fixed = TRUE
  )
}

cand[
  ,
  tissue_norm :=
    normalize_tissue(tissue)
]

fin <- copy(
  smr[study == "FinnGen"]
)

fin[
  ,
  tissue_norm :=
    normalize_tissue(tissue)
]

# --------------------------------------------------------------------------
# Match the exact UKB-carried hypothesis in FinnGen:
# gene + tissue + probeID
# --------------------------------------------------------------------------
rep <- merge(
  cand[
    ,
    .(
      Gene,
      Module,
      tissue_UKB = tissue,
      tissue_norm,
      probeID,
      topSNP_UKB = topSNP,
      p_SMR_UKB = p_SMR,
      FDR_UKB = FDR_cohortwide,
      p_HEIDI_UKB = p_HEIDI
    )
  ],
  fin[
    ,
    .(
      Gene,
      tissue_FinnGen = tissue,
      tissue_norm,
      probeID,
      topSNP_FinnGen = topSNP,
      b_SMR_FinnGen = b_SMR,
      se_SMR_FinnGen = se_SMR,
      p_SMR_FinnGen = p_SMR,
      p_HEIDI_FinnGen = p_HEIDI,
      nsnp_HEIDI_FinnGen = nsnp_HEIDI
    )
  ],
  by = c(
    "Gene",
    "tissue_norm",
    "probeID"
  ),
  all.x = TRUE,
  allow.cartesian = FALSE
)

if (nrow(rep) != 2L) {
  stop(
    "Expected two exact FinnGen replication rows; found ",
    nrow(rep)
  )
}

if (
  anyNA(rep$p_SMR_FinnGen) ||
  anyNA(rep$p_HEIDI_FinnGen)
) {
  stop(
    "Missing FinnGen SMR/HEIDI value for a carried-forward hypothesis."
  )
}

# --------------------------------------------------------------------------
# Targeted replication BH across the two hypotheses
# --------------------------------------------------------------------------
rep[
  ,
  replication_BH :=
    p.adjust(
      p_SMR_FinnGen,
      method = "BH"
    )
]

rep[
  ,
  replication_pass :=
    replication_BH < 0.05 &
    p_HEIDI_FinnGen > 0.01
]

setorder(
  rep,
  Gene
)

# --------------------------------------------------------------------------
# Expected values from the independent C7 audit
# --------------------------------------------------------------------------
smad3 <- rep[Gene == "SMAD3"]
plekha6 <- rep[Gene == "PLEKHA6"]

pass <-
  nrow(smad3) == 1L &&
  nrow(plekha6) == 1L &&
  abs(smad3$p_SMR_FinnGen - 0.001071558) < 1e-12 &&
  abs(smad3$replication_BH - 0.002143116) < 1e-12 &&
  smad3$p_HEIDI_FinnGen > 0.01 &&
  isTRUE(smad3$replication_pass) &&
  abs(plekha6$replication_BH - 0.373608) < 1e-12 &&
  !isTRUE(plekha6$replication_pass)

cat("\n============================================================\n")
cat("CANONICAL FINNGEN TARGETED REPLICATION\n")
cat("============================================================\n\n")

print(
  rep[
    ,
    .(
      Gene,
      Module,
      tissue_norm,
      probeID,
      p_SMR_UKB,
      FDR_UKB,
      p_HEIDI_UKB,
      p_SMR_FinnGen,
      replication_BH,
      p_HEIDI_FinnGen,
      replication_pass
    )
  ]
)

cat("\n============================================================\n")
cat("FINAL STATUS\n")
cat("============================================================\n\n")

cat(
  "Canonical FinnGen replication validation pass:",
  pass,
  "\n"
)

cat("\nNo canonical file was modified.\n")

if (!pass) {
  stop(
    "Canonical FinnGen replication validation failed. ",
    "Do not productionize Figure 10B-C yet."
  )
}
