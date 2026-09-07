# ============================================================================
# Reviewer C7 audit:
# Multiple-testing assessment for CellChat communication probabilities
# underlying Figure 18.
#
# Primary targeted family:
#   Within each condition, all source x target x ligand-receptor tests
#   belonging to the pre-specified TGFb pathway.
#
# Broader sensitivity family:
#   Within each condition, all tests across the full pre-specified
#   18-pathway CellChat database used in the original analysis.
#
# Additional sensitivity:
#   Pool Young and Aged hypotheses into a single family.
#
# Permutation-resolution sensitivity:
#   If stored CellChat p-values are compatible with the default nboot grid,
#   reconstruct a plus-one empirical p-value:
#       (b + 1) / (B + 1)
#
# No canonical files are modified.
# ============================================================================

rm(list = ls(all.names = TRUE))
gc()

suppressPackageStartupMessages({
  library(CellChat)
  library(data.table)
})

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

# --------------------------------------------------------------------------
# Input / output
# --------------------------------------------------------------------------
obj_file <- file.path(
  repo_root,
  "processed_results",
  "12_CellChat",
  "cellchat_merged.rds"
)

outdir <- file.path(
  repo_root,
  "processed_results",
  "12_CellChat",
  "multiple_testing_audit"
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(obj_file)) {
  stop("CellChat merged object not found: ", obj_file)
}

cat("CellChat input:\n", obj_file, "\n\n", sep = "")

x <- readRDS(obj_file)

datasets <- c("Young", "Aged")

if (!all(datasets %in% names(x@net))) {
  stop("Expected Young and Aged networks were not found.")
}

# --------------------------------------------------------------------------
# Original pre-specified pathway set
# --------------------------------------------------------------------------
pathways_prespecified <- c(
  "TGFb",
  "COLLAGEN",
  "FN1",
  "LAMININ",
  "THBS",
  "TENASCIN",
  "VTN",
  "HSPG",
  "CDH",
  "CDH1",
  "CDH5",
  "PCDH",
  "ICAM",
  "VCAM",
  "JAM",
  "MMP",
  "SPP1",
  "PERIOSTIN"
)

cat("Number of pre-specified pathways:", length(pathways_prespecified), "\n")

# --------------------------------------------------------------------------
# CellChat / permutation information
# --------------------------------------------------------------------------
cellchat_version <- as.character(
  packageVersion("CellChat")
)

nboot_formal <- formals(CellChat::computeCommunProb)$nboot

nboot_default <- tryCatch(
  eval(nboot_formal),
  error = function(e) NA_real_
)

cat("CellChat version:", cellchat_version, "\n")
cat("computeCommunProb() default nboot:", nboot_default, "\n\n")

# --------------------------------------------------------------------------
# Interaction -> pathway lookup
# --------------------------------------------------------------------------
db <- x@DB$interaction

required_db_cols <- c(
  "interaction_name",
  "pathway_name"
)

if (!all(required_db_cols %in% names(db))) {
  stop(
    "CellChat DB does not contain required interaction/pathway columns."
  )
}

lr_map <- unique(
  data.table(
    interaction_name = as.character(db$interaction_name),
    pathway = as.character(db$pathway_name)
  )
)

# --------------------------------------------------------------------------
# Flatten one condition
# --------------------------------------------------------------------------
flatten_condition <- function(ds) {

  prob <- x@net[[ds]]$prob
  pval <- x@net[[ds]]$pval

  if (!identical(dim(prob), dim(pval))) {
    stop("prob and pval dimensions differ for ", ds)
  }

  lr_names <- dimnames(prob)[[3]]

  if (
    is.null(lr_names) ||
    length(lr_names) != dim(prob)[3]
  ) {
    lr_names <- x@net[[ds]]$LRs

    if (is.data.frame(lr_names)) {
      if ("interaction_name" %in% names(lr_names)) {
        lr_names <- lr_names$interaction_name
      } else {
        lr_names <- lr_names[[1]]
      }
    }

    lr_names <- as.character(lr_names)

    dimnames(prob)[[3]] <- lr_names
    dimnames(pval)[[3]] <- lr_names
  }

  pdt <- as.data.table(
    as.table(pval)
  )

  setnames(
    pdt,
    c(
      "source",
      "target",
      "interaction_name",
      "p_raw"
    )
  )

  prob_dt <- as.data.table(
    as.table(prob)
  )

  setnames(
    prob_dt,
    c(
      "source",
      "target",
      "interaction_name",
      "probability"
    )
  )

  dt <- merge(
    pdt,
    prob_dt,
    by = c(
      "source",
      "target",
      "interaction_name"
    ),
    all = TRUE,
    sort = FALSE
  )

  dt[
    ,
    dataset := ds
  ]

  dt <- merge(
    dt,
    lr_map,
    by = "interaction_name",
    all.x = TRUE,
    sort = FALSE
  )

  dt[
    ,
    in_prespecified_pathway :=
      pathway %in% pathways_prespecified
  ]

  dt[
    ,
    is_TGFb := pathway == "TGFb"
  ]

  dt
}

all_dt <- rbindlist(
  lapply(
    datasets,
    flatten_condition
  ),
  fill = TRUE
)

# Keep only hypotheses belonging to the DB used in this analysis.
analysis_dt <- all_dt[
  in_prespecified_pathway == TRUE &
    is.finite(p_raw)
]

if (nrow(analysis_dt) == 0L) {
  stop("No valid CellChat hypotheses recovered.")
}

# --------------------------------------------------------------------------
# Verify actual pathways present
# --------------------------------------------------------------------------
actual_pathways <- sort(
  unique(
    analysis_dt$pathway
  )
)

cat("Pathways represented in stored network:", length(actual_pathways), "\n")
cat(
  paste(actual_pathways, collapse = ", "),
  "\n\n"
)

missing_pathways <- setdiff(
  pathways_prespecified,
  actual_pathways
)

if (length(missing_pathways) > 0L) {
  cat(
    "Pre-specified pathways absent from stored test array:\n",
    paste(missing_pathways, collapse = ", "),
    "\n\n",
    sep = ""
  )
}

# --------------------------------------------------------------------------
# BH correction within each condition
# --------------------------------------------------------------------------
analysis_dt[
  ,
  p_BH_all18_within_condition :=
    p.adjust(
      p_raw,
      method = "BH"
    ),
  by = dataset
]

analysis_dt[
  ,
  p_BH_TGFb_within_condition := NA_real_
]

analysis_dt[
  is_TGFb == TRUE,
  p_BH_TGFb_within_condition :=
    p.adjust(
      p_raw,
      method = "BH"
    ),
  by = dataset
]

# --------------------------------------------------------------------------
# Pooled Young + Aged sensitivities
# --------------------------------------------------------------------------
analysis_dt[
  ,
  p_BH_all18_pooled_conditions :=
    p.adjust(
      p_raw,
      method = "BH"
    )
]

analysis_dt[
  ,
  p_BH_TGFb_pooled_conditions := NA_real_
]

analysis_dt[
  is_TGFb == TRUE,
  p_BH_TGFb_pooled_conditions :=
    p.adjust(
      p_raw,
      method = "BH"
    )
]

# --------------------------------------------------------------------------
# Check empirical permutation grid and build plus-one sensitivity
# --------------------------------------------------------------------------
grid_compatible <- FALSE

if (
  is.finite(nboot_default) &&
    nboot_default > 0
) {

  grid_error <- abs(
    analysis_dt$p_raw * nboot_default -
      round(analysis_dt$p_raw * nboot_default)
  )

  grid_compatible <- all(
    grid_error < 1e-8,
    na.rm = TRUE
  )
}

cat(
  "Stored p-values compatible with default nboot grid:",
  grid_compatible,
  "\n"
)

analysis_dt[
  ,
  p_plus1 := NA_real_
]

if (grid_compatible) {

  analysis_dt[
    ,
    permutation_exceedance_count :=
      round(p_raw * nboot_default)
  ]

  analysis_dt[
    ,
    p_plus1 :=
      (
        permutation_exceedance_count + 1
      ) /
      (
        nboot_default + 1
      )
  ]

  analysis_dt[
    ,
    p_plus1_BH_all18_within_condition :=
      p.adjust(
        p_plus1,
        method = "BH"
      ),
    by = dataset
  ]

  analysis_dt[
    ,
    p_plus1_BH_TGFb_within_condition :=
      NA_real_
  ]

  analysis_dt[
    is_TGFb == TRUE,
    p_plus1_BH_TGFb_within_condition :=
      p.adjust(
        p_plus1,
        method = "BH"
      ),
    by = dataset
  ]

} else {

  analysis_dt[
    ,
    permutation_exceedance_count := NA_real_
  ]

  analysis_dt[
    ,
    p_plus1_BH_all18_within_condition := NA_real_
  ]

  analysis_dt[
    ,
    p_plus1_BH_TGFb_within_condition := NA_real_
  ]
}

# --------------------------------------------------------------------------
# Condition-level summary
# --------------------------------------------------------------------------
summary_list <- lapply(
  datasets,
  function(ds) {

    z <- analysis_dt[
      dataset == ds
    ]

    t <- z[
      is_TGFb == TRUE
    ]

    data.table(
      dataset = ds,

      n_all18_tests =
        nrow(z),

      n_TGFb_tests =
        nrow(t),

      n_all18_LR_pairs =
        uniqueN(z$interaction_name),

      n_TGFb_LR_pairs =
        uniqueN(t$interaction_name),

      n_raw_p_zero_all18 =
        sum(z$p_raw == 0, na.rm = TRUE),

      n_raw_p_zero_TGFb =
        sum(t$p_raw == 0, na.rm = TRUE),

      n_raw_p_lt_0_05_all18 =
        sum(z$p_raw < 0.05, na.rm = TRUE),

      n_raw_p_lt_0_05_TGFb =
        sum(t$p_raw < 0.05, na.rm = TRUE),

      n_BH_all18_lt_0_05 =
        sum(
          z$p_BH_all18_within_condition < 0.05,
          na.rm = TRUE
        ),

      n_BH_TGFb_lt_0_05 =
        sum(
          t$p_BH_TGFb_within_condition < 0.05,
          na.rm = TRUE
        ),

      n_BH_TGFb_pooled_lt_0_05 =
        sum(
          t$p_BH_TGFb_pooled_conditions < 0.05,
          na.rm = TRUE
        ),

      n_plus1_BH_all18_lt_0_05 =
        sum(
          z$p_plus1_BH_all18_within_condition < 0.05,
          na.rm = TRUE
        ),

      n_plus1_BH_TGFb_lt_0_05 =
        sum(
          t$p_plus1_BH_TGFb_within_condition < 0.05,
          na.rm = TRUE
        )
    )
  }
)

summary_dt <- rbindlist(
  summary_list
)

# --------------------------------------------------------------------------
# Figure 18F targeted bubble-plot hypothesis subset
# --------------------------------------------------------------------------
bubble_sources <- c(
  "FAP3",
  "Other FAPs"
)

bubble_targets <- c(
  "MuSC",
  "Myogenic",
  "Immune",
  "Endothelial",
  "Vascular stromal"
)

fig18f_audit <- analysis_dt[
  is_TGFb == TRUE &
    source %in% bubble_sources &
    target %in% bubble_targets
]

fig18f_audit[
  ,
  raw_significant :=
    p_raw < 0.05
]

fig18f_audit[
  ,
  BH_TGFb_significant :=
    p_BH_TGFb_within_condition < 0.05
]

fig18f_audit[
  ,
  BH_all18_significant :=
    p_BH_all18_within_condition < 0.05
]

fig18f_audit[
  ,
  plus1_BH_TGFb_significant :=
    p_plus1_BH_TGFb_within_condition < 0.05
]

fig18f_summary <- fig18f_audit[
  ,
  .(
    n_hypotheses = .N,

    n_probability_gt_0 =
      sum(probability > 0, na.rm = TRUE),

    n_raw_p_lt_0_05 =
      sum(raw_significant, na.rm = TRUE),

    n_BH_TGFb_lt_0_05 =
      sum(BH_TGFb_significant, na.rm = TRUE),

    n_BH_all18_lt_0_05 =
      sum(BH_all18_significant, na.rm = TRUE),

    n_plus1_BH_TGFb_lt_0_05 =
      sum(
        plus1_BH_TGFb_significant,
        na.rm = TRUE
      )
  ),
  by = dataset
]

# --------------------------------------------------------------------------
# Figure 18E Aged TGFb ligand-receptor contribution comparison
# --------------------------------------------------------------------------
contribution_table <- function(
    dt,
    criterion_name,
    pass_col
) {

  z <- copy(
    dt[
      dataset == "Aged" &
        is_TGFb == TRUE
    ]
  )

  z[
    ,
    pass := get(pass_col)
  ]

  z[
    is.na(pass),
    pass := FALSE
  ]

  z[
    ,
    probability_used :=
      fifelse(
        pass,
        probability,
        0
      )
  ]

  out <- z[
    ,
    .(
      raw_significant_probability_sum =
        sum(
          probability_used,
          na.rm = TRUE
        ),

      number_of_significant_cellpair_edges =
        sum(
          probability_used > 0,
          na.rm = TRUE
        )
    ),
    by = .(
      interaction_name,
      pathway
    )
  ]

  out <- out[
    raw_significant_probability_sum > 0
  ]

  if (nrow(out) == 0L) {
    return(
      data.table(
        criterion = character(),
        interaction_name = character(),
        pathway = character(),
        raw_significant_probability_sum = numeric(),
        number_of_significant_cellpair_edges = integer(),
        relative_raw_contribution = numeric()
      )
    )
  }

  total <- sum(
    out$raw_significant_probability_sum,
    na.rm = TRUE
  )

  out[
    ,
    relative_raw_contribution :=
      raw_significant_probability_sum / total
  ]

  out[
    ,
    criterion := criterion_name
  ]

  setcolorder(
    out,
    c(
      "criterion",
      "interaction_name",
      "pathway",
      "raw_significant_probability_sum",
      "number_of_significant_cellpair_edges",
      "relative_raw_contribution"
    )
  )

  setorder(
    out,
    -relative_raw_contribution
  )

  out
}

analysis_dt[
  ,
  raw_significant :=
    p_raw < 0.05
]

analysis_dt[
  ,
  BH_TGFb_significant :=
    p_BH_TGFb_within_condition < 0.05
]

analysis_dt[
  ,
  BH_all18_significant :=
    p_BH_all18_within_condition < 0.05
]

analysis_dt[
  ,
  plus1_BH_TGFb_significant :=
    p_plus1_BH_TGFb_within_condition < 0.05
]

fig18e_comparison <- rbindlist(
  list(
    contribution_table(
      analysis_dt,
      "raw_p_lt_0.05",
      "raw_significant"
    ),

    contribution_table(
      analysis_dt,
      "BH_TGFb_within_condition",
      "BH_TGFb_significant"
    ),

    contribution_table(
      analysis_dt,
      "BH_all18_within_condition",
      "BH_all18_significant"
    ),

    contribution_table(
      analysis_dt,
      "plus1_BH_TGFb_within_condition",
      "plus1_BH_TGFb_significant"
    )
  ),
  fill = TRUE
)

# --------------------------------------------------------------------------
# TGFb source-target network edge sensitivity
# Relevant to aggregate TGFb network panels
# --------------------------------------------------------------------------
tgfb_edge_summary <- analysis_dt[
  is_TGFb == TRUE,
  .(
    raw_weight =
      sum(
        fifelse(
          raw_significant,
          probability,
          0
        ),
        na.rm = TRUE
      ),

    BH_TGFb_weight =
      sum(
        fifelse(
          BH_TGFb_significant,
          probability,
          0
        ),
        na.rm = TRUE
      ),

    BH_all18_weight =
      sum(
        fifelse(
          BH_all18_significant,
          probability,
          0
        ),
        na.rm = TRUE
      ),

    plus1_BH_TGFb_weight =
      sum(
        fifelse(
          plus1_BH_TGFb_significant,
          probability,
          0
        ),
        na.rm = TRUE
      )
  ),
  by = .(
    dataset,
    source,
    target
  )
]

tgfb_network_summary <- tgfb_edge_summary[
  ,
  .(
    n_raw_nonzero_edges =
      sum(raw_weight > 0),

    n_BH_TGFb_nonzero_edges =
      sum(BH_TGFb_weight > 0),

    n_BH_all18_nonzero_edges =
      sum(BH_all18_weight > 0),

    n_plus1_BH_TGFb_nonzero_edges =
      sum(plus1_BH_TGFb_weight > 0),

    total_raw_weight =
      sum(raw_weight),

    total_BH_TGFb_weight =
      sum(BH_TGFb_weight),

    total_BH_all18_weight =
      sum(BH_all18_weight),

    total_plus1_BH_TGFb_weight =
      sum(plus1_BH_TGFb_weight)
  ),
  by = dataset
]

# --------------------------------------------------------------------------
# Export
# --------------------------------------------------------------------------
fwrite(
  analysis_dt,
  file.path(
    outdir,
    "CellChat_all_hypotheses_with_multiple_testing.tsv"
  ),
  sep = "\t"
)

fwrite(
  summary_dt,
  file.path(
    outdir,
    "CellChat_multiple_testing_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  fig18f_audit,
  file.path(
    outdir,
    "Figure18F_TGFb_bubble_multiple_testing_audit.tsv"
  ),
  sep = "\t"
)

fwrite(
  fig18f_summary,
  file.path(
    outdir,
    "Figure18F_TGFb_bubble_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  fig18e_comparison,
  file.path(
    outdir,
    "Figure18E_TGFb_contribution_multiple_testing_audit.tsv"
  ),
  sep = "\t"
)

fwrite(
  tgfb_edge_summary,
  file.path(
    outdir,
    "TGFb_source_target_edge_multiple_testing_audit.tsv"
  ),
  sep = "\t"
)

fwrite(
  tgfb_network_summary,
  file.path(
    outdir,
    "TGFb_network_multiple_testing_summary.tsv"
  ),
  sep = "\t"
)

# --------------------------------------------------------------------------
# Console summary
# --------------------------------------------------------------------------
cat("\n============================================================\n")
cat("CELLCHAT MULTIPLE-TESTING SUMMARY\n")
cat("============================================================\n\n")
print(summary_dt)

cat("\n============================================================\n")
cat("FIGURE 18F TARGETED BUBBLE SUBSET\n")
cat("============================================================\n\n")
print(fig18f_summary)

cat("\n============================================================\n")
cat("TGFb NETWORK EDGE SUMMARY\n")
cat("============================================================\n\n")
print(tgfb_network_summary)

cat("\n============================================================\n")
cat("AGED TGFb LR CONTRIBUTIONS\n")
cat("============================================================\n\n")

if (nrow(fig18e_comparison) == 0L) {
  cat("No LR pair survived one or more corrected criteria.\n")
} else {
  print(fig18e_comparison)
}

cat("\n============================================================\n")
cat("P-VALUE DISTRIBUTION CHECK\n")
cat("============================================================\n\n")

for (ds in datasets) {

  z <- analysis_dt[
    dataset == ds &
      is_TGFb == TRUE
  ]

  cat(ds, "\n")
  cat("  TGFb tests :", nrow(z), "\n")
  cat("  raw p = 0  :", sum(z$p_raw == 0), "\n")
  cat(
    "  unique raw p-values:",
    length(unique(z$p_raw)),
    "\n"
  )

  cat("  smallest raw p-values:\n")
  print(
    head(
      sort(
        unique(z$p_raw)
      ),
      10
    )
  )

  cat("\n")
}

cat("\nAudit outputs written to:\n")
cat(outdir, "\n")
cat("\nNo canonical file was modified.\n")
