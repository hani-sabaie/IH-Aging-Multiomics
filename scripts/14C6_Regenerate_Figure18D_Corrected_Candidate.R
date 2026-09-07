# ============================================================================
# Reviewer C7
# Regenerate Figure 18D candidate using multiplicity-corrected TGFb
# communications.
#
# Inference is NOT rerun.
#
# Observed communication probabilities come from the validated historical
# CellChat analysis and are identical to the nboot=1000 lightweight rerun.
#
# Significance criterion used for plotting:
#   plus-one empirical permutation P = (b + 1)/(1000 + 1)
#   followed by BH correction across 588 pre-specified TGFb hypotheses
#   separately within Young and Aged.
#
# The original CellChat plotting function and circle layout are preserved.
#
# No canonical figure/source-data file is modified.
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

historical_file <- file.path(
  repo_root,
  "processed_results",
  "12_CellChat",
  "cellchat_merged.rds"
)

audit_file <- file.path(
  repo_root,
  "processed_results",
  "12_CellChat",
  "multiple_testing_audit",
  "TGFb_lightweight",
  "TGFb_nboot1000_all_tests_combined.tsv"
)

candidate_source_file <- file.path(
  repo_root,
  "processed_results",
  "12_CellChat",
  "multiple_testing_audit",
  "figure18_revised_candidate",
  "Figure18D_TGFb_pathway_network_corrected_candidate.csv"
)

figdir <- file.path(
  repo_root,
  "outputs",
  "reviewer_c7",
  "Figure18_revised_candidate"
)

dir.create(
  figdir,
  recursive = TRUE,
  showWarnings = FALSE
)

for (f in c(
  historical_file,
  audit_file,
  candidate_source_file
)) {
  if (!file.exists(f)) {
    stop("Required file not found: ", f)
  }
}

x <- readRDS(historical_file)
audit <- fread(audit_file)
candidate_source <- fread(candidate_source_file)

# --------------------------------------------------------------------------
# Exact historical TGFb LR definitions
# --------------------------------------------------------------------------
tgfb_lr <- x@DB$interaction[
  x@DB$interaction$pathway_name == "TGFb",
  ,
  drop = FALSE
]

rownames(tgfb_lr) <- tgfb_lr$interaction_name

if (nrow(tgfb_lr) != 12L) {
  stop(
    "Expected 12 historical TGFb LR interactions; found ",
    nrow(tgfb_lr)
  )
}

datasets <- c("Young", "Aged")

# --------------------------------------------------------------------------
# Construct a plotting-only CellChat object for one condition.
# --------------------------------------------------------------------------
make_corrected_plot_object <- function(ds) {

  cat("\nPreparing corrected plotting object:", ds, "\n")

  cells <- rownames(x@meta)[
    x@meta$datasets == ds
  ]

  meta <- x@meta[
    cells,
    ,
    drop = FALSE
  ]

  # The merged object retains the exact signaling matrix used historically.
  data_input <- x@data.signaling[
    ,
    cells,
    drop = FALSE
  ]

  cc <- createCellChat(
    object = data_input,
    meta = meta,
    group.by = "group_cellchat"
  )

  cc@DB <- x@DB

  # netVisual_aggregate() searches object@LR$LRsig for the requested pathway.
  cc@LR$LRsig <- tgfb_lr

  z <- copy(
    audit[
      dataset == ds
    ]
  )

  expected_n <- 7L * 7L * 12L

  if (nrow(z) != expected_n) {
    stop(
      ds,
      ": expected ",
      expected_n,
      " TGFb hypotheses; found ",
      nrow(z)
    )
  }

  groups <- levels(cc@idents)

  if (length(groups) != 7L) {
    stop(
      ds,
      ": expected seven CellChat groups; found ",
      length(groups)
    )
  }

  lr_order <- rownames(tgfb_lr)

  # ------------------------------------------------------------------------
  # Reconstruct exact 7 x 7 x 12 probability and corrected-P arrays.
  # ------------------------------------------------------------------------
  prob_arr <- array(
    0,
    dim = c(
      length(groups),
      length(groups),
      length(lr_order)
    ),
    dimnames = list(
      groups,
      groups,
      lr_order
    )
  )

  pval_arr <- array(
    1,
    dim = dim(prob_arr),
    dimnames = dimnames(prob_arr)
  )

  for (i in seq_len(nrow(z))) {

    src <- as.character(z$source[i])
    tgt <- as.character(z$target[i])
    lr  <- as.character(z$interaction_name[i])

    if (
      !(src %in% groups) ||
      !(tgt %in% groups) ||
      !(lr %in% lr_order)
    ) {
      stop(
        "Unexpected source/target/LR in ",
        ds,
        ": ",
        src,
        " / ",
        tgt,
        " / ",
        lr
      )
    }

    prob_arr[src, tgt, lr] <-
      z$probability[i]

    # This is the inferential quantity used by the revised plot.
    pval_arr[src, tgt, lr] <-
      z$p_plus1_BH_TGFb[i]
  }

  cc@net$prob <- prob_arr
  cc@net$pval <- pval_arr

  # Store retained LR names for consistency with CellChat objects.
  retained_lr <- lr_order[
    apply(
      prob_arr * (pval_arr < 0.05),
      3,
      sum,
      na.rm = TRUE
    ) > 0
  ]

  cc@net$LRs <- retained_lr

  cc
}

# --------------------------------------------------------------------------
# Validate aggregation of each corrected plotting object against the
# candidate source-data produced in 14C5.
# --------------------------------------------------------------------------
validate_corrected_object <- function(
    ds,
    cc
) {

  prob <- cc@net$prob
  pval <- cc@net$pval

  prob[
    pval >= 0.05
  ] <- 0

  pathway_mat <- apply(
    prob,
    c(1, 2),
    sum,
    na.rm = TRUE
  )

  candidate <- candidate_source[
    condition == ds
  ]

  expected <- matrix(
    0,
    nrow = length(rownames(pathway_mat)),
    ncol = length(colnames(pathway_mat)),
    dimnames = dimnames(pathway_mat)
  )

  for (i in seq_len(nrow(candidate))) {

    expected[
      candidate$source[i],
      candidate$target[i]
    ] <- candidate$TGFb_pathway_probability[i]
  }

  delta <- abs(
    pathway_mat -
      expected
  )

  max_delta <- max(
    delta,
    na.rm = TRUE
  )

  match <- isTRUE(
    all.equal(
      as.numeric(pathway_mat),
      as.numeric(expected),
      tolerance = 1e-12,
      check.attributes = FALSE
    )
  )

  cat(
    ds,
    " corrected aggregation max |delta| = ",
    format(
      max_delta,
      scientific = TRUE,
      digits = 12
    ),
    "\n",
    sep = ""
  )

  cat(
    ds,
    " corrected aggregation match = ",
    match,
    "\n",
    sep = ""
  )

  if (!match) {
    stop(
      ds,
      ": corrected plotting object does not match candidate source-data."
    )
  }

  invisible(pathway_mat)
}

# --------------------------------------------------------------------------
# Build corrected objects
# --------------------------------------------------------------------------
young_cc <- make_corrected_plot_object("Young")
aged_cc  <- make_corrected_plot_object("Aged")

young_mat <- validate_corrected_object(
  "Young",
  young_cc
)

aged_mat <- validate_corrected_object(
  "Aged",
  aged_cc
)

# --------------------------------------------------------------------------
# Plot with the same CellChat circle-network function used historically.
# --------------------------------------------------------------------------
plot_one <- function(
    cc,
    ds,
    filename_stem
) {

  png(
    file.path(
      figdir,
      paste0(
        filename_stem,
        ".png"
      )
    ),
    width = 2000,
    height = 2000,
    res = 300
  )

  netVisual_aggregate(
    object = cc,
    signaling = "TGFb",
    layout = "circle",
    thresh = 0.05
  )

  dev.off()

  pdf(
    file.path(
      figdir,
      paste0(
        filename_stem,
        ".pdf"
      )
    ),
    width = 2000 / 300,
    height = 2000 / 300,
    useDingbats = FALSE
  )

  netVisual_aggregate(
    object = cc,
    signaling = "TGFb",
    layout = "circle",
    thresh = 0.05
  )

  dev.off()
}

plot_one(
  young_cc,
  "Young",
  "Figure18D_TGFb_Young_corrected_candidate"
)

plot_one(
  aged_cc,
  "Aged",
  "Figure18D_TGFb_Aged_corrected_candidate"
)

# --------------------------------------------------------------------------
# Compact validation record
# --------------------------------------------------------------------------
validation <- data.table(
  condition = c(
    "Young",
    "Aged"
  ),

  correction =
    "plus-one empirical permutation P; BH within 588 TGFb tests",

  nboot = 1000L,

  aggregation_match_candidate_source =
    c(
      TRUE,
      TRUE
    ),

  affected_edge =
    c(
      "Vascular stromal -> FAP3",
      "None"
    ),

  historical_weight =
    c(
      0.00451896545096,
      NA_real_
    ),

  corrected_weight =
    c(
      0.00298370327737,
      NA_real_
    )
)

fwrite(
  validation,
  file.path(
    figdir,
    "Figure18D_corrected_candidate_validation.tsv"
  ),
  sep = "\t"
)

cat("\n============================================================\n")
cat("FIGURE 18D CANDIDATE COMPLETE\n")
cat("============================================================\n\n")

cat("Files written to:\n")
cat(figdir, "\n\n")

print(
  list.files(
    figdir,
    full.names = FALSE
  )
)

cat("\nNo canonical Figure 18 or source-data file was modified.\n")
