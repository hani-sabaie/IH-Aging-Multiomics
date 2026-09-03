rm(list = ls(all.names = TRUE))
gc()

library(ggplot2)
library(data.table)

# -------------------------------------------------------------------------
# Repository root
# -------------------------------------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (length(file_arg) == 1) {
  script_path <- sub("^--file=", "", file_arg)
  script_dir <- dirname(normalizePath(script_path))
  repo_root <- normalizePath(file.path(script_dir, ".."))
} else {
  repo_root <- normalizePath(".")
}

# -------------------------------------------------------------------------
# Input / output
# -------------------------------------------------------------------------
input_file <- file.path(
  repo_root,
  "processed_results",
  "13_mouse_validation",
  "Fibro_mouse_chromVAR_results.rds"
)

if (!file.exists(input_file)) {
  stop("Mouse chromVAR plot object not found: ", input_file)
}

srcdir <- file.path(
  repo_root,
  "source_data",
  "figure_source_data",
  "Figure_20"
)

dir.create(srcdir, recursive = TRUE, showWarnings = FALSE)

cat("Loading original Figure 20B plot objects...\n")
x <- readRDS(input_file)

expected_fap <- c(
  "Pgr- Fibroblast",
  "Pgr+ Fibroblast",
  "Adipocytes"
)

tf_to_motif <- c(
  ETS2  = "MA1484.2",
  ETS1  = "MA0098.4",
  FOS   = "MA1127.1",
  SMAD3 = "MA0795.1"
)

if (!all(expected_fap %in% names(x))) {
  stop("Expected mouse cell types are missing from the plot object.")
}

# =========================================================================
# Cell-level source data
# =========================================================================
long_list <- list()
k <- 1L

for (cf in expected_fap) {

  if (
    !"raincloud" %in% names(x[[cf]]) ||
    !all(names(tf_to_motif) %in% names(x[[cf]]$raincloud))
  ) {
    stop("Incomplete raincloud objects for: ", cf)
  }

  # Each TF-specific ggplot contains the same underlying 5-column
  # cell-level data frame. Use one plot as the canonical source.
  dat <- x[[cf]]$raincloud[["ETS2"]]$data

  if (is.null(dat) || nrow(dat) == 0) {
    stop("No plot data found for: ", cf)
  }

  required_cols <- c("condition", unname(tf_to_motif))

  if (!all(required_cols %in% colnames(dat))) {
    stop(
      "Required chromVAR columns missing for ",
      cf,
      ": ",
      paste(setdiff(required_cols, colnames(dat)), collapse = ", ")
    )
  }

  cell_ids <- rownames(dat)

  if (is.null(cell_ids) || any(cell_ids == "")) {
    stop("Cell identifiers are missing for: ", cf)
  }

  for (tf in names(tf_to_motif)) {

    motif <- tf_to_motif[[tf]]

    long_list[[k]] <- data.table(
      cell = cell_ids,
      cell_type = cf,
      condition = as.character(dat$condition),
      TF = tf,
      motif_id = motif,
      chromVAR_deviation = as.numeric(dat[[motif]])
    )

    k <- k + 1L
  }
}

fig20b <- rbindlist(long_list, use.names = TRUE)

fig20b[
  ,
  cell_type_order := match(cell_type, expected_fap)
]

fig20b[
  ,
  TF_order := match(TF, names(tf_to_motif))
]

fig20b[
  ,
  condition_order := match(condition, c("Veh", "EP"))
]

setorder(
  fig20b,
  cell_type_order,
  TF_order,
  condition_order,
  cell
)

fwrite(
  fig20b,
  file.path(
    srcdir,
    "Figure20B_mouse_chromVAR_cell_level_source_data.csv"
  )
)

# =========================================================================
# Wilcoxon tests: Veh vs EP
# =========================================================================
stats_list <- list()
k <- 1L

for (cf in expected_fap) {
  for (tf in names(tf_to_motif)) {

    d <- fig20b[
      cell_type == cf &
        TF == tf &
        condition %in% c("Veh", "EP")
    ]

    veh <- d[condition == "Veh", chromVAR_deviation]
    ep  <- d[condition == "EP", chromVAR_deviation]

    n_veh <- length(veh)
    n_ep  <- length(ep)

    if (n_veh > 0 && n_ep > 0) {

      wt <- suppressWarnings(
        stats::wilcox.test(
          veh,
          ep,
          alternative = "two.sided"
        )
      )

      pval <- unname(wt$p.value)

    } else {

      pval <- NA_real_
    }

    stats_list[[k]] <- data.table(
      cell_type = cf,
      TF = tf,
      motif_id = tf_to_motif[[tf]],
      n_Veh = n_veh,
      n_EP = n_ep,
      mean_Veh = mean(veh, na.rm = TRUE),
      mean_EP = mean(ep, na.rm = TRUE),
      mean_diff_EP_minus_Veh =
        mean(ep, na.rm = TRUE) - mean(veh, na.rm = TRUE),
      median_Veh = median(veh, na.rm = TRUE),
      median_EP = median(ep, na.rm = TRUE),
      median_diff_EP_minus_Veh =
        median(ep, na.rm = TRUE) - median(veh, na.rm = TRUE),
      wilcoxon_p_value = pval
    )

    k <- k + 1L
  }
}

fig20b_stats <- rbindlist(stats_list)

# Benjamini-Hochberg correction across the 12 prespecified
# cell-type x transcription-factor comparisons.
fig20b_stats[
  ,
  p_adj_BH := p.adjust(wilcoxon_p_value, method = "BH")
]

fig20b_stats[
  ,
  p_significance := fifelse(
    is.na(wilcoxon_p_value), NA_character_,
    fifelse(
      wilcoxon_p_value <= 0.0001, "****",
      fifelse(
        wilcoxon_p_value <= 0.001, "***",
        fifelse(
          wilcoxon_p_value <= 0.01, "**",
          fifelse(
            wilcoxon_p_value <= 0.05, "*", "ns"
          )
        )
      )
    )
  )
]

fig20b_stats[
  ,
  p_adj_significance := fifelse(
    is.na(p_adj_BH), NA_character_,
    fifelse(
      p_adj_BH <= 0.0001, "****",
      fifelse(
        p_adj_BH <= 0.001, "***",
        fifelse(
          p_adj_BH <= 0.01, "**",
          fifelse(
            p_adj_BH <= 0.05, "*", "ns"
          )
        )
      )
    )
  )
]

fwrite(
  fig20b_stats,
  file.path(
    srcdir,
    "Figure20B_mouse_chromVAR_Wilcoxon_Veh_vs_EP.csv"
  )
)

# =========================================================================
# Descriptive statistics
# =========================================================================
fig20b_summary <- fig20b[
  ,
  .(
    n = .N,
    mean_chromVAR = mean(chromVAR_deviation, na.rm = TRUE),
    median_chromVAR = median(chromVAR_deviation, na.rm = TRUE),
    sd_chromVAR = sd(chromVAR_deviation, na.rm = TRUE)
  ),
  by = .(
    cell_type,
    TF,
    motif_id,
    condition
  )
]

fig20b_summary[
  ,
  cell_type_order := match(cell_type, expected_fap)
]

fig20b_summary[
  ,
  TF_order := match(TF, names(tf_to_motif))
]

fig20b_summary[
  ,
  condition_order := match(condition, c("Veh", "EP"))
]

setorder(
  fig20b_summary,
  cell_type_order,
  TF_order,
  condition_order
)

fwrite(
  fig20b_summary,
  file.path(
    srcdir,
    "Figure20B_mouse_chromVAR_summary_statistics.csv"
  )
)

# =========================================================================
# Motif mapping
# =========================================================================
motif_mapping <- data.table(
  TF = names(tf_to_motif),
  motif_id = unname(tf_to_motif)
)

fwrite(
  motif_mapping,
  file.path(
    srcdir,
    "Figure20B_mouse_chromVAR_motif_mapping.csv"
  )
)

# =========================================================================
# Validation
# =========================================================================
cat("\n===== Figure 20B source-data summary =====\n")

cat("Cell-level rows:", nrow(fig20b), "\n")

cat("\nUnique cells by cell type and condition:\n")
print(
  unique(
    fig20b[
      ,
      .(cell, cell_type, condition)
    ]
  )[
    ,
    .N,
    by = .(cell_type, condition)
  ]
)

cat("\nRows by TF and cell type:\n")
print(
  fig20b[
    ,
    .N,
    by = .(cell_type, TF)
  ]
)

cat("\nWilcoxon results:\n")
print(fig20b_stats)

cat("\nDone.\n")
