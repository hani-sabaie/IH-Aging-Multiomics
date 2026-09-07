# ============================================================================
# Reviewer C7 audit
# Lightweight reconstruction of historical TGFb CellChat probabilities.
#
# Purpose:
#   Validate a memory-efficient TGFb-only representation BEFORE running
#   nboot = 1000.
#
# Strategy:
#   - Use historical merged CellChat data.signaling (161 genes), not Seurat.
#   - Retain only genes needed by the 12 TGFb LR interactions.
#   - Add an unused scale-anchor row when necessary so max(data) is identical
#     to the historical condition-specific signaling matrix.
#   - Restrict inference using LR.use = TGFb interactions.
#   - nboot = 1 because observed probability does not depend on nboot.
#
# No canonical file is modified.
# ============================================================================

rm(list = ls(all.names = TRUE))
gc()

suppressPackageStartupMessages({
  library(CellChat)
  library(data.table)
})

future::plan("sequential")

# --------------------------------------------------------------------------
# Repository root
# --------------------------------------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (length(file_arg) == 1L) {
  script_path <- sub("^--file=", "", file_arg)
  script_dir <- dirname(normalizePath(script_path))
  repo_root <- normalizePath(file.path(script_dir, ".."))
} else {
  repo_root <- normalizePath(".")
}

hist_file <- file.path(
  repo_root,
  "processed_results",
  "12_CellChat",
  "cellchat_merged.rds"
)

if (!file.exists(hist_file)) {
  stop("Historical CellChat object not found: ", hist_file)
}

outdir <- file.path(
  repo_root,
  "processed_results",
  "12_CellChat",
  "multiple_testing_audit",
  "TGFb_lightweight"
)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cat("Loading historical merged CellChat object...\n")
x <- readRDS(hist_file)

# --------------------------------------------------------------------------
# TGFb LR set
# --------------------------------------------------------------------------
db <- x@DB$interaction

tgfb_lr <- db[
  db$pathway_name == "TGFb",
  ,
  drop = FALSE
]

if (nrow(tgfb_lr) != 12L) {
  stop(
    "Expected 12 historical TGFb LR pairs; found ",
    nrow(tgfb_lr)
  )
}

rownames(tgfb_lr) <- tgfb_lr$interaction_name

cat("Historical TGFb LR pairs:", nrow(tgfb_lr), "\n")

# --------------------------------------------------------------------------
# Determine all genes needed by the selected LR interactions
# --------------------------------------------------------------------------
complex_db <- x@DB$complex
cofactor_db <- x@DB$cofactor

expand_complex <- function(items) {

  items <- unique(
    items[
      !is.na(items) &
        items != ""
    ]
  )

  out <- character()

  for (z in items) {

    if (z %in% rownames(x@data.signaling)) {

      out <- c(out, z)

    } else if (z %in% rownames(complex_db)) {

      vals <- unlist(
        complex_db[z, , drop = FALSE],
        use.names = FALSE
      )

      vals <- as.character(vals)

      vals <- vals[
        !is.na(vals) &
          vals != ""
      ]

      out <- c(out, vals)

    } else {

      stop(
        "Unable to resolve ligand/receptor term: ",
        z
      )
    }
  }

  unique(out)
}

lr_genes <- expand_complex(
  c(
    as.character(tgfb_lr$ligand),
    as.character(tgfb_lr$receptor)
  )
)

# --------------------------------------------------------------------------
# Add genes referenced through agonist / antagonist / coreceptor definitions,
# if TGFb uses any of them.
# --------------------------------------------------------------------------
cofactor_ref_cols <- intersect(
  c(
    "agonist",
    "antagonist",
    "co_A_receptor",
    "co_I_receptor"
  ),
  colnames(tgfb_lr)
)

cofactor_refs <- character()

for (cc in cofactor_ref_cols) {

  vals <- as.character(
    tgfb_lr[[cc]]
  )

  vals <- vals[
    !is.na(vals) &
      vals != ""
  ]

  cofactor_refs <- c(
    cofactor_refs,
    vals
  )
}

cofactor_refs <- unique(cofactor_refs)

cofactor_genes <- character()

if (length(cofactor_refs) > 0L) {

  for (z in cofactor_refs) {

    if (!(z %in% rownames(cofactor_db))) {
      stop(
        "TGFb cofactor reference not found in cofactor DB: ",
        z
      )
    }

    vals <- unlist(
      cofactor_db[z, , drop = FALSE],
      use.names = FALSE
    )

    vals <- as.character(vals)

    vals <- vals[
      !is.na(vals) &
        vals != ""
    ]

    cofactor_genes <- c(
      cofactor_genes,
      vals
    )
  }
}

needed_genes <- unique(
  c(
    lr_genes,
    cofactor_genes
  )
)

needed_genes <- intersect(
  needed_genes,
  rownames(x@data.signaling)
)

if (length(needed_genes) == 0L) {
  stop("No TGFb genes recovered from historical data.signaling.")
}

cat("\nTGFb genes required for inference:\n")
print(needed_genes)

cat(
  "\nNumber of required signaling genes:",
  length(needed_genes),
  "\n"
)

if (length(cofactor_refs) > 0L) {
  cat("\nTGFb cofactor references:\n")
  print(cofactor_refs)
}

# --------------------------------------------------------------------------
# Lightweight condition reconstruction
# --------------------------------------------------------------------------
run_validation <- function(ds) {

  cat("\n============================================================\n")
  cat("VALIDATING:", ds, "\n")
  cat("============================================================\n\n")

  cells <- rownames(x@meta)[
    x@meta$datasets == ds
  ]

  if (length(cells) == 0L) {
    stop("No cells found for dataset ", ds)
  }

  meta <- x@meta[
    cells,
    ,
    drop = FALSE
  ]

  full_data <- x@data.signaling[
    ,
    cells,
    drop = FALSE
  ]

  historical_max <- max(
    full_data,
    na.rm = TRUE
  )

  light_data <- full_data[
    needed_genes,
    ,
    drop = FALSE
  ]

  light_max_before_anchor <- max(
    light_data,
    na.rm = TRUE
  )

  anchor_needed <-
    abs(
      historical_max -
        light_max_before_anchor
    ) > 1e-15

  if (anchor_needed) {

    anchor <- matrix(
      0,
      nrow = 1L,
      ncol = ncol(light_data),
      dimnames = list(
        "CELLCHAT_SCALE_ANCHOR",
        colnames(light_data)
      )
    )

    anchor[1, 1] <- historical_max

    light_data <- rbind(
      light_data,
      anchor
    )
  }

  cat("Cells                      :", length(cells), "\n")
  cat("Historical signaling genes :", nrow(full_data), "\n")
  cat("Lightweight genes          :", nrow(light_data), "\n")
  cat("Historical max(data)       :", format(historical_max, digits = 16), "\n")
  cat(
    "TGFb-only max before anchor:",
    format(light_max_before_anchor, digits = 16),
    "\n"
  )
  cat("Scale anchor required      :", anchor_needed, "\n")
  cat(
    "Lightweight max after anchor:",
    format(max(light_data), digits = 16),
    "\n\n"
  )

  # createCellChat establishes identities/options/meta.
  cc <- createCellChat(
    object = light_data,
    meta = meta,
    group.by = "group_cellchat"
  )

  # Use the exact historical DB structures.
  cc@DB <- x@DB

  # We deliberately supply the compact signaling matrix directly.
  # computeCommunProb(raw.use=TRUE) reads object@data.signaling.
  cc@data.signaling <- light_data

  # Historical groups are all represented in each condition.
  cc@idents <- droplevels(
    factor(
      meta$group_cellchat
    )
  )

  names(cc@idents) <- rownames(meta)

  tm <- system.time({

    cc <- computeCommunProb(
      cc,
      type = "triMean",
      LR.use = tgfb_lr,
      raw.use = TRUE,
      population.size = FALSE,
      nboot = 1,
      seed.use = 1L
    )

  })

  cat(
    "Validation computeCommunProb elapsed:",
    unname(tm["elapsed"]),
    "seconds\n"
  )

  # ------------------------------------------------------------------------
  # Compare new observed probabilities with the historical TGFb probabilities
  # ------------------------------------------------------------------------
  old_prob <- x@net[[ds]]$prob

  old_lr <- dimnames(old_prob)[[3]]

  idx <- match(
    rownames(tgfb_lr),
    old_lr
  )

  if (anyNA(idx)) {
    stop(
      "Some TGFb LR pairs are absent from historical probability array."
    )
  }

  old_tgfb <- old_prob[
    ,
    ,
    idx,
    drop = FALSE
  ]

  new_tgfb <- cc@net$prob

  # Force identical third-dimension ordering for comparison.
  new_idx <- match(
    rownames(tgfb_lr),
    dimnames(new_tgfb)[[3]]
  )

  if (anyNA(new_idx)) {
    stop(
      "Some TGFb LR pairs are absent from reconstructed probability array."
    )
  }

  new_tgfb <- new_tgfb[
    ,
    ,
    new_idx,
    drop = FALSE
  ]

  same_dim <- identical(
    dim(old_tgfb),
    dim(new_tgfb)
  )

  same_source <- identical(
    dimnames(old_tgfb)[[1]],
    dimnames(new_tgfb)[[1]]
  )

  same_target <- identical(
    dimnames(old_tgfb)[[2]],
    dimnames(new_tgfb)[[2]]
  )

  delta <- abs(
    as.numeric(old_tgfb) -
      as.numeric(new_tgfb)
  )

  max_delta <- max(
    delta,
    na.rm = TRUE
  )

  mean_delta <- mean(
    delta,
    na.rm = TRUE
  )

  exact_numeric_match <-
    same_dim &&
    same_source &&
    same_target &&
    isTRUE(
      all.equal(
        as.numeric(old_tgfb),
        as.numeric(new_tgfb),
        tolerance = 1e-12,
        check.attributes = FALSE
      )
    )

  cat("\nProbability comparison:\n")
  cat("Same dimensions      :", same_dim, "\n")
  cat("Same source order    :", same_source, "\n")
  cat("Same target order    :", same_target, "\n")
  cat(
    "Max |delta prob|    :",
    format(max_delta, scientific = TRUE, digits = 12),
    "\n"
  )
  cat(
    "Mean |delta prob|   :",
    format(mean_delta, scientific = TRUE, digits = 12),
    "\n"
  )
  cat(
    "Match within 1e-12 :",
    exact_numeric_match,
    "\n"
  )

  data.table(
    dataset = ds,
    n_cells = length(cells),
    n_historical_signaling_genes = nrow(full_data),
    n_lightweight_rows = nrow(light_data),
    n_required_TGFb_genes = length(needed_genes),
    scale_anchor_required = anchor_needed,
    historical_max = historical_max,
    TGFb_max_before_anchor = light_max_before_anchor,
    same_dimensions = same_dim,
    same_source_order = same_source,
    same_target_order = same_target,
    max_abs_probability_delta = max_delta,
    mean_abs_probability_delta = mean_delta,
    probability_match_1e12 = exact_numeric_match,
    elapsed_seconds = unname(tm["elapsed"])
  )
}

# --------------------------------------------------------------------------
# Young + Aged validation
# --------------------------------------------------------------------------
validation <- rbindlist(
  list(
    run_validation("Young"),
    run_validation("Aged")
  )
)

cat("\n============================================================\n")
cat("FINAL LIGHTWEIGHT VALIDATION SUMMARY\n")
cat("============================================================\n\n")

print(validation)

fwrite(
  validation,
  file.path(
    outdir,
    "TGFb_lightweight_probability_validation.tsv"
  ),
  sep = "\t"
)

cat("\nOutput written to:\n")
cat(outdir, "\n")

cat("\nNo canonical CellChat or Figure 18 file was modified.\n")

future::plan("sequential")
