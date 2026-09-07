# ============================================================================
# Compare historical canonical mouse pseudobulk results with the new
# full-precision three-group rebuild.
#
# Focus:
#   EP vs Veh
#   overlapping eligible cell types / genes
#   effect sizes, raw P, BH-FDR
#   impact of historical 4-decimal rounding
#
# Audit only.
# ============================================================================

rm(list = ls(all.names = TRUE))

suppressPackageStartupMessages({
  library(data.table)
})

old_file <- file.path(
  "processed_results",
  "13_mouse_validation",
  "bulk_de_all_contrasts_mouse.csv"
)

new_file <- file.path(
  "processed_results",
  "13_mouse_validation",
  "three_group_fullprecision",
  "mouse_EP_vs_Veh_threegroup_fullprecision.csv"
)

old <- fread(old_file)
new <- fread(new_file)

cat("===== FILE DIMENSIONS =====\n")
cat("Old:", nrow(old), "rows\n")
cat("New:", nrow(new), "rows\n\n")

cat("Old columns:\n")
print(names(old))

cat("\nNew columns:\n")
print(names(new))

# --------------------------------------------------------------------------
# Historical table uses Veh_vs_EP, so convert sign to EP_vs_Veh.
# --------------------------------------------------------------------------

old_epveh <- old[
  contrast == "Veh_vs_EP"
]

cat("\nHistorical Veh_vs_EP rows:", nrow(old_epveh), "\n")

old_epveh[
  ,
  avg_log2FC_EP_vs_Veh :=
    -avg_log2FC
]

old_epveh[
  ,
  key :=
    paste(
      celltype,
      gene,
      sep = "|||"
    )
]

new[
  ,
  key :=
    paste(
      celltype,
      gene,
      sep = "|||"
    )
]

shared <- merge(
  old_epveh[, .(
    key,
    celltype_old = celltype,
    gene_old = gene,
    old_log2FC_EP_vs_Veh = avg_log2FC_EP_vs_Veh,
    old_p = p_val,
    old_FDR = p_val_adj
  )],
  new[, .(
    key,
    celltype_new = celltype,
    gene_new = gene,
    new_log2FC_EP_vs_Veh = avg_log2FC,
    new_p = p_val,
    new_FDR = p_val_adj
  )],
  by = "key"
)

cat("Shared rows:", nrow(shared), "\n")

if (nrow(shared) < 1L) {
  stop("No shared rows.")
}

# --------------------------------------------------------------------------
# Numerical differences
# --------------------------------------------------------------------------

shared[
  ,
  delta_logFC :=
    abs(
      old_log2FC_EP_vs_Veh -
      new_log2FC_EP_vs_Veh
    )
]

shared[
  ,
  delta_p :=
    abs(
      old_p -
      new_p
    )
]

shared[
  ,
  delta_FDR :=
    abs(
      old_FDR -
      new_FDR
    )
]

cat("\n===== NUMERICAL AGREEMENT =====\n")

cat(
  "Max |delta log2FC|:",
  max(shared$delta_logFC, na.rm = TRUE),
  "\n"
)

cat(
  "Median |delta log2FC|:",
  median(shared$delta_logFC, na.rm = TRUE),
  "\n"
)

cat(
  "Max |delta P|:",
  max(shared$delta_p, na.rm = TRUE),
  "\n"
)

cat(
  "Median |delta P|:",
  median(shared$delta_p, na.rm = TRUE),
  "\n"
)

cat(
  "Max |delta FDR|:",
  max(shared$delta_FDR, na.rm = TRUE),
  "\n"
)

cat(
  "Median |delta FDR|:",
  median(shared$delta_FDR, na.rm = TRUE),
  "\n"
)

# --------------------------------------------------------------------------
# Significance comparison
# --------------------------------------------------------------------------

shared[
  ,
  old_sig :=
    old_FDR < 0.05
]

shared[
  ,
  new_sig :=
    new_FDR < 0.05
]

shared[
  ,
  sig_changed :=
    old_sig != new_sig
]

cat("\n===== SIGNIFICANCE AGREEMENT =====\n")

cat(
  "Old significant among shared:",
  sum(shared$old_sig),
  "\n"
)

cat(
  "New significant among shared:",
  sum(shared$new_sig),
  "\n"
)

cat(
  "Status changes:",
  sum(shared$sig_changed),
  "\n"
)

if (any(shared$sig_changed)) {

  cat("\nRows with significance-status changes:\n")

  print(
    shared[
      sig_changed == TRUE,
      .(
        celltype = celltype_new,
        gene = gene_new,
        old_log2FC_EP_vs_Veh,
        new_log2FC_EP_vs_Veh,
        old_p,
        new_p,
        old_FDR,
        new_FDR
      )
    ][
      order(
        pmin(old_FDR, new_FDR)
      )
    ],
    nrows = 100
  )
}

# --------------------------------------------------------------------------
# Historical rounding audit
# --------------------------------------------------------------------------

cat("\n===== HISTORICAL ROUNDING AUDIT =====\n")

# Canonical old p/FDR values were written at 4 decimals.
# Identify boundary values where rounding could affect <0.05 membership.

boundary <- old_epveh[
  p_val_adj >= 0.045 &
    p_val_adj <= 0.055
]

cat(
  "Historical rows with FDR in [0.045, 0.055]:",
  nrow(boundary),
  "\n"
)

if (nrow(boundary) > 0) {
  print(
    boundary[
      order(p_val_adj),
      .(
        gene,
        celltype,
        contrast,
        avg_log2FC,
        p_val,
        p_val_adj
      )
    ],
    nrows = 100
  )
}

# --------------------------------------------------------------------------
# Smad3 explicit check
# --------------------------------------------------------------------------

cat("\n===== SMAD3 COMPARISON =====\n")

print(
  shared[
    tolower(gene_new) == "smad3",
    .(
      celltype = celltype_new,
      old_log2FC_EP_vs_Veh,
      new_log2FC_EP_vs_Veh,
      old_p,
      new_p,
      old_FDR,
      new_FDR,
      old_sig,
      new_sig
    )
  ],
  nrows = 100
)

# --------------------------------------------------------------------------
# Save comparison
# --------------------------------------------------------------------------

outdir <- file.path(
  "processed_results",
  "13_mouse_validation",
  "three_group_fullprecision"
)

fwrite(
  shared,
  file.path(
    outdir,
    "canonical_vs_fullprecision_EP_vs_Veh.csv"
  )
)

cat(
  "\nPASS: canonical-vs-fullprecision comparison completed.\n"
)
