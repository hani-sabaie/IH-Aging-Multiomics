rm(list = ls(all.names = TRUE))

suppressPackageStartupMessages({
  library(data.table)
})

base_dir <- file.path(
  "processed_results",
  "13_mouse_validation"
)

primary_file <- file.path(
  base_dir,
  "bulk_de_all_contrasts_mouse.csv"
)

sensitivity_source <- file.path(
  base_dir,
  "pseudobulk_rebuild_audit",
  "Smad3_TypeIIa_model_sensitivity.csv"
)

out_file <- file.path(
  base_dir,
  "Smad3_TypeIIa_primary_vs_EP_Veh_sensitivity.csv"
)

primary <- fread(primary_file)
sens <- fread(sensitivity_source)

# --------------------------------------------------------------------------
# Primary three-condition model
# Historical contrast is Veh_vs_EP, so reverse sign for EP_vs_Veh reporting.
# --------------------------------------------------------------------------

p <- primary[
  tolower(gene) == "smad3" &
    celltype == "Type IIa Myofiber" &
    contrast == "Veh_vs_EP"
]

if (nrow(p) != 1L) {
  stop("Primary Type IIa Smad3 row not uniquely identified.")
}

# --------------------------------------------------------------------------
# Direct EP-vs-Veh-only sensitivity model.
#
# This model uses the same original-design mouse-level inclusion rule as the
# primary analysis and the same four EP and three Veh mice. No additional
# minimum-nuclei threshold is applied; the only design change is omission of
# EPR from the fitted model.
# --------------------------------------------------------------------------

s <- sens[
  model == "B_correct_EP_vs_Veh"
]

if (nrow(s) != 1L) {
  stop("EP-vs-Veh sensitivity row not uniquely identified.")
}

# Validate exact values established by the independent audits.
if (
  abs((-p$avg_log2FC) - 1.3517707722949) > 1e-10 ||
  abs(p$p_val - 0.000449950771698) > 1e-10 ||
  abs(p$p_val_adj - 0.0598454692221) > 1e-10
) {
  stop("Primary result validation failed.")
}

if (
  abs(s$Smad3_log2FC_EP_vs_Veh - 1.35946755488249) > 1e-10 ||
  abs(s$Smad3_P - 0.000105832074432827) > 1e-10 ||
  abs(s$Smad3_BH_FDR - 0.0331016201968554) > 1e-10
) {
  stop("Sensitivity result validation failed.")
}

out <- rbind(
  data.table(
    analysis = "Primary_three_condition_model",
    celltype = "Type IIa Myofiber",
    contrast = "EP_vs_Veh",
    n_EP = 4L,
    n_Veh = 3L,
    n_EPR = 3L,
    n_genes_tested = p$n_tests_within_family,
    Smad3_log2FC_EP_vs_Veh = -p$avg_log2FC,
    p_val = p$p_val,
    BH_FDR = p$p_val_adj,
    significant = p$p_val_adj < 0.05
  ),
  data.table(
    analysis = "Sensitivity_EP_vs_Veh_only",
    celltype = "Type IIa Myofiber",
    contrast = "EP_vs_Veh",
    n_EP = s$n_EP,
    n_Veh = s$n_Veh,
    n_EPR = 0L,
    n_genes_tested = s$n_genes_tested,
    Smad3_log2FC_EP_vs_Veh = s$Smad3_log2FC_EP_vs_Veh,
    p_val = s$Smad3_P,
    BH_FDR = s$Smad3_BH_FDR,
    significant = s$Smad3_BH_FDR < 0.05
  )
)

fwrite(
  out,
  out_file
)

cat("===== SMAD3 PRIMARY VS SENSITIVITY =====\n")
print(out, digits = 12)

cat(
  "\nEffect-size difference:",
  abs(
    out$Smad3_log2FC_EP_vs_Veh[1] -
      out$Smad3_log2FC_EP_vs_Veh[2]
  ),
  "\n"
)

cat(
  "\nPASS: Smad3 mouse sensitivity result promoted.\n"
)
