rm(list = ls(all.names = TRUE))

suppressPackageStartupMessages({
  library(data.table)
})

old_all <- fread(
  "processed_results/13_mouse_validation/bulk_de_all_contrasts_mouse.csv"
)

old_sig <- fread(
  "processed_results/13_mouse_validation/bulk_de_sig_all_contrasts_mouse.csv"
)

new_all <- fread(
  paste0(
    "processed_results/13_mouse_validation/",
    "original_design_fullprecision/",
    "bulk_de_all_contrasts_mouse_fullprecision.csv"
  )
)

new_sens <- fread(
  paste0(
    "processed_results/13_mouse_validation/",
    "original_design_fullprecision/",
    "bulk_de_sig_BH_FDR005_absLog2FCge1_sensitivity.csv"
  )
)

new_primary <- fread(
  paste0(
    "processed_results/13_mouse_validation/",
    "original_design_fullprecision/",
    "bulk_de_sig_BH_FDR005.csv"
  )
)

cat("===== FINAL MOUSE PSEUDOBULK VALIDATION =====\n")

cat("Old all rows:", nrow(old_all), "\n")
cat("New all rows:", nrow(new_all), "\n")
cat("Old significant rows:", nrow(old_sig), "\n")
cat("New FC-sensitivity rows:", nrow(new_sens), "\n")
cat("New FDR-only rows:", nrow(new_primary), "\n\n")

keyfun <- function(x) {
  paste(
    x$gene,
    x$celltype,
    x$contrast,
    sep = "|||"
  )
}

old_all[, key := keyfun(old_all)]
new_all[, key := keyfun(new_all)]

if (
  anyDuplicated(old_all$key) ||
  anyDuplicated(new_all$key)
) {
  stop("Duplicate test keys detected.")
}

cat(
  "All-test key sets identical:",
  setequal(old_all$key, new_all$key),
  "\n"
)

shared <- merge(
  old_all[, .(
    key,
    old_logFC = avg_log2FC,
    old_p = p_val,
    old_FDR = p_val_adj
  )],
  new_all[, .(
    key,
    new_logFC = avg_log2FC,
    new_p = p_val,
    new_FDR = p_val_adj
  )],
  by = "key"
)

cat("Shared tests:", nrow(shared), "\n")

cat(
  "Median |delta logFC|:",
  median(abs(shared$old_logFC - shared$new_logFC)),
  "\n"
)

cat(
  "Median |delta P|:",
  median(abs(shared$old_p - shared$new_p)),
  "\n"
)

cat(
  "Median |delta FDR|:",
  median(abs(shared$old_FDR - shared$new_FDR)),
  "\n"
)

# --------------------------------------------------------------------------
# Historical significant table vs full-precision FC sensitivity
# --------------------------------------------------------------------------

old_sig[, key := keyfun(old_sig)]
new_sens[, key := keyfun(new_sens)]
new_primary[, key := keyfun(new_primary)]

only_old <- setdiff(
  old_sig$key,
  new_sens$key
)

only_new <- setdiff(
  new_sens$key,
  old_sig$key
)

cat("\n===== FC-SENSITIVITY MEMBERSHIP =====\n")
cat("Old sig:", length(unique(old_sig$key)), "\n")
cat("New sensitivity:", length(unique(new_sens$key)), "\n")
cat("Only historical:", length(only_old), "\n")
cat("Only full precision:", length(only_new), "\n")

if (length(only_old) > 0) {
  cat("\nOnly historical:\n")
  print(
    old_sig[
      key %in% only_old,
      .(
        gene,
        celltype,
        contrast,
        avg_log2FC,
        p_val,
        p_val_adj
      )
    ]
  )
}

if (length(only_new) > 0) {
  cat("\nOnly full precision:\n")
  print(
    new_sens[
      key %in% only_new,
      .(
        gene,
        celltype,
        contrast,
        avg_log2FC,
        p_val,
        p_val_adj
      )
    ]
  )
}

# --------------------------------------------------------------------------
# New FDR-only rows added when FC threshold is removed
# --------------------------------------------------------------------------

added_fdr_only <- setdiff(
  new_primary$key,
  new_sens$key
)

cat("\n===== FDR-ONLY PRIMARY =====\n")
cat(
  "Additional significant rows after removing hard FC threshold:",
  length(added_fdr_only),
  "\n"
)

cat(
  "Expected difference:",
  nrow(new_primary) - nrow(new_sens),
  "\n"
)

if (
  length(added_fdr_only) !=
    nrow(new_primary) - nrow(new_sens)
) {
  stop("FDR-only membership inconsistency.")
}

# --------------------------------------------------------------------------
# Smad3 exact values
# --------------------------------------------------------------------------

cat("\n===== SMAD3 EP/Veh =====\n")

print(
  new_all[
    tolower(gene) == "smad3" &
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

cat("\nPASS: final mouse pseudobulk validation completed.\n")
