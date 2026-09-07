rm(list = ls(all.names = TRUE))

suppressPackageStartupMessages({
  library(data.table)
})

base_dir <- file.path(
  "processed_results",
  "13_mouse_validation"
)

candidate_file <- file.path(
  base_dir,
  "original_design_fullprecision",
  "bulk_de_all_contrasts_mouse_fullprecision.csv"
)

canonical_all <- file.path(
  base_dir,
  "bulk_de_all_contrasts_mouse.csv"
)

canonical_sig <- file.path(
  base_dir,
  "bulk_de_sig_all_contrasts_mouse.csv"
)

canonical_sensitivity <- file.path(
  base_dir,
  "bulk_de_sig_all_contrasts_mouse_absLog2FCge1_sensitivity.csv"
)

table_s6 <- file.path(
  "processed_results",
  "supplementary_tables",
  "Table_S6_mouse_Smad3_pseudobulk_results.csv"
)

cat("===== MOUSE PSEUDOBULK CANONICAL PROMOTION =====\n")

old_sig <- fread(canonical_sig)
candidate <- fread(candidate_file)

required <- c(
  "gene",
  "celltype",
  "contrast",
  "n_tests_within_family",
  "avg_log2FC",
  "p_val",
  "p_val_adj"
)

if (!all(required %in% names(candidate))) {
  stop("Candidate is missing required production columns.")
}

# Match the production-script output schema exactly.
all_new <- candidate[
  ,
  ..required
]

setorder(
  all_new,
  p_val_adj,
  p_val
)

keyfun <- function(x) {
  paste(
    x$gene,
    x$celltype,
    x$contrast,
    sep = "|||"
  )
}

all_new[, key := keyfun(all_new)]

if (anyDuplicated(all_new$key)) {
  stop("Duplicate gene/celltype/contrast keys.")
}

if (nrow(all_new) != 213960L) {
  stop(
    "Expected 213,960 total tests; found ",
    nrow(all_new)
  )
}

# Validate family sizes and independent BH adjustment.
family_check <- all_new[
  ,
  {
    expected_n <- .N

    declared_n <- unique(
      n_tests_within_family
    )

    if (
      length(declared_n) != 1L ||
      declared_n != expected_n
    ) {
      stop(
        "Family-size mismatch: ",
        celltype[1],
        " / ",
        contrast[1]
      )
    }

    bh <- p.adjust(
      p_val,
      method = "BH"
    )

    .(
      n_tests = .N,
      max_BH_delta = max(
        abs(
          bh -
          p_val_adj
        )
      )
    )
  },
  by = .(
    celltype,
    contrast
  )
]

if (
  nrow(family_check) != 48L ||
  any(family_check$max_BH_delta > 1e-12)
) {
  stop("BH family validation failed.")
}

primary <- all_new[
  p_val_adj < 0.05
]

sensitivity <- all_new[
  p_val_adj < 0.05 &
    abs(avg_log2FC) >= 1
]

if (nrow(primary) != 1576L) {
  stop(
    "Expected 1,576 primary FDR-only rows; found ",
    nrow(primary)
  )
}

if (nrow(sensitivity) != 1479L) {
  stop(
    "Expected 1,479 sensitivity rows; found ",
    nrow(sensitivity)
  )
}

# Historical 1,479-row membership must be exactly preserved.
old_sig[, key := keyfun(old_sig)]

if (
  !setequal(
    old_sig$key,
    sensitivity$key
  )
) {
  stop(
    "Historical sensitivity membership is not identical."
  )
}

# Remove temporary keys before writing production files.
all_new[, key := NULL]
primary[, key := NULL]
sensitivity[, key := NULL]

# Smad3 supplementary table.
smad3 <- all_new[
  tolower(gene) == "smad3"
]

if (nrow(smad3) != 30L) {
  stop(
    "Expected 30 Smad3 rows; found ",
    nrow(smad3)
  )
}

if (any(smad3$p_val_adj < 0.05)) {
  stop(
    "Unexpected BH-significant Smad3 result."
  )
}

# Explicitly validate the manuscript-relevant EP/Veh values.
type2a <- smad3[
  celltype == "Type IIa Myofiber" &
    contrast == "Veh_vs_EP"
]

pgrneg <- smad3[
  celltype == "Pgr- Fibroblast" &
    contrast == "Veh_vs_EP"
]

if (
  nrow(type2a) != 1L ||
  abs(
    type2a$p_val_adj -
      0.0598454692221
  ) > 1e-10
) {
  stop("Type IIa Smad3 validation failed.")
}

if (
  nrow(pgrneg) != 1L ||
  abs(
    pgrneg$p_val_adj -
      0.7823168402853
  ) > 1e-10
) {
  stop("Pgr- fibroblast Smad3 validation failed.")
}

# --------------------------------------------------------------------------
# Promote
# --------------------------------------------------------------------------

fwrite(
  all_new,
  canonical_all
)

fwrite(
  primary,
  canonical_sig
)

fwrite(
  sensitivity,
  canonical_sensitivity
)

fwrite(
  smad3,
  table_s6
)

cat("\n===== PROMOTED OUTPUTS =====\n")
cat("All tests:", nrow(all_new), "\n")
cat("Primary FDR-only:", nrow(primary), "\n")
cat("FC-filtered sensitivity:", nrow(sensitivity), "\n")
cat("Table S6 Smad3 rows:", nrow(smad3), "\n")
cat("BH families:", nrow(family_check), "\n")
cat(
  "Maximum independent BH delta:",
  max(family_check$max_BH_delta),
  "\n"
)

cat("\n===== SMAD3 EP/Veh =====\n")

print(
  smad3[
    contrast == "Veh_vs_EP",
    .(
      gene,
      celltype,
      avg_log2FC,
      p_val,
      p_val_adj
    )
  ],
  digits = 12
)

cat(
  "\nPASS: canonical mouse pseudobulk promotion completed.\n"
)
