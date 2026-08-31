rm(list = ls(all.names = TRUE))
gc()

library(arrow)
library(data.table)
library(dplyr)
library(ggplot2)

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
spatial_dir <- Sys.getenv(
  "FIG21_SPATIAL_DIR",
  unset = file.path(repo_root, "data", "GSE288663", "Spatial")
)

if (!dir.exists(spatial_dir)) {
  stop(
    "Spatial input directory not found: ",
    spatial_dir,
    "\nSet FIG21_SPATIAL_DIR to the local GSE288663 Spatial directory."
  )
}

srcdir <- file.path(
  repo_root,
  "source_data",
  "figure_source_data",
  "Figure_21"
)

dir.create(srcdir, recursive = TRUE, showWarnings = FALSE)

sample_info <- data.frame(
  sample = c("1450Veh", "1603Veh", "1460EP", "1461EP"),
  folder = c(
    "GSM8772705_1450Veh",
    "GSM8772706_1603Veh",
    "GSM8772707_1460EP",
    "GSM8772708_1461EP"
  ),
  stringsAsFactors = FALSE
)

smad3_gene <- "Smad3"

smad3_pos_genes <- c(
  "Gprc5a", "Hivep1", "Eml4", "Agtr1a",
  "Medag", "Galnt2", "Hivep2", "Spag17",
  "Clstn2", "Glis3", "Cxcl2", "Trio"
)

smad3_neg_genes <- c(
  "Wdfy1", "Il34", "Mtm1", "Emb"
)

genes_of_interest <- unique(
  c(
    smad3_gene,
    smad3_pos_genes,
    smad3_neg_genes
  )
)

# =========================================================================
# Figure 21E-F: transcript coordinates
# Same filtering as scripts/18_Spatial_Data_Analysis.R
# =========================================================================
tx_list <- list()

for (i in seq_len(nrow(sample_info))) {

  samp <- sample_info$sample[i]

  parquet_file <- file.path(
    spatial_dir,
    sample_info$folder[i],
    "transcripts.parquet"
  )

  if (!file.exists(parquet_file)) {
    stop("Missing transcript parquet: ", parquet_file)
  }

  cat("Loading:", samp, "\n")

  tx <- as.data.frame(
    arrow::read_parquet(parquet_file)
  )

  required <- c(
    "feature_name",
    "x_location",
    "y_location"
  )

  if (!all(required %in% colnames(tx))) {
    stop(
      "Required columns missing from ",
      samp,
      ": ",
      paste(
        setdiff(required, colnames(tx)),
        collapse = ", "
      )
    )
  }

  # Original analysis filters
  if ("is_gene" %in% colnames(tx)) {
    tx <- tx[
      tx$is_gene == TRUE,
      ,
      drop = FALSE
    ]
  }

  if ("qv" %in% colnames(tx)) {
    tx <- tx[
      tx$qv >= 20,
      ,
      drop = FALSE
    ]
  }

  tx <- tx[
    tx$feature_name %in% genes_of_interest,
    ,
    drop = FALSE
  ]

  dt <- as.data.table(tx)

  setnames(
    dt,
    c("feature_name", "x_location", "y_location"),
    c("gene", "x", "y")
  )

  dt[, sample := samp]

  dt[
    ,
    group := fifelse(
      gene == smad3_gene,
      "Smad3",
      fifelse(
        gene %in% smad3_pos_genes,
        "SMAD3_Pos",
        fifelse(
          gene %in% smad3_neg_genes,
          "SMAD3_Neg",
          "Other"
        )
      )
    )
  ]

  # Keep the variables required to reproduce the figure,
  # plus useful Xenium identifiers/QV when available.
  keep_cols <- intersect(
    c(
      "sample",
      "gene",
      "group",
      "x",
      "y",
      "qv",
      "is_gene",
      "transcript_id",
      "cell_id"
    ),
    names(dt)
  )

  dt <- dt[, ..keep_cols]

  tx_list[[samp]] <- dt

  cat(
    "  retained transcripts:",
    nrow(dt),
    "\n"
  )

  rm(tx, dt)
  gc()
}

tx_all <- rbindlist(
  tx_list,
  use.names = TRUE,
  fill = TRUE
)

sample_order <- sample_info$sample
group_order <- c(
  "Smad3",
  "SMAD3_Pos",
  "SMAD3_Neg"
)

tx_all[
  ,
  sample_order := match(sample, sample_order)
]

tx_all[
  ,
  group_order := match(group, group_order)
]

setorder(
  tx_all,
  sample_order,
  group_order,
  gene,
  x,
  y
)

fwrite(
  tx_all,
  file.path(
    srcdir,
    "Figure21E_spatial_transcript_coordinates_source_data.csv"
  )
)

# -------------------------------------------------------------------------
# Positive-target coordinates underlying Figure 21F
# -------------------------------------------------------------------------
fig21f_pos <- tx_all[
  group == "SMAD3_Pos"
]

fwrite(
  fig21f_pos,
  file.path(
    srcdir,
    "Figure21F_SMAD3_Pos_transcript_coordinates_source_data.csv"
  )
)

# -------------------------------------------------------------------------
# Transcript counts by sample and group
# -------------------------------------------------------------------------
count_table <- tx_all[
  ,
  .N,
  by = .(sample, group)
]

count_wide <- dcast(
  count_table,
  sample ~ group,
  value.var = "N",
  fill = 0
)

count_wide[
  ,
  sample_order_tmp := match(sample, sample_order)
]

setorder(
  count_wide,
  sample_order_tmp
)

count_wide[
  ,
  sample_order_tmp := NULL
]

fwrite(
  count_wide,
  file.path(
    srcdir,
    "Figure21_transcript_counts_by_sample.csv"
  )
)

# Validate against original processed count table
orig_count_file <- file.path(
  repo_root,
  "processed_results",
  "14_spatial",
  "Transcript_counts_Smad3_network_by_sample.csv"
)

if (file.exists(orig_count_file)) {

  orig_counts <- fread(orig_count_file)

  common_cols <- intersect(
    names(orig_counts),
    names(count_wide)
  )

  orig_cmp <- copy(orig_counts)[
    ,
    ..common_cols
  ]

  new_cmp <- copy(count_wide)[
    ,
    ..common_cols
  ]

  setorder(orig_cmp, sample)
  setorder(new_cmp, sample)

  if (!identical(
    as.data.frame(orig_cmp),
    as.data.frame(new_cmp)
  )) {
    stop(
      "Reconstructed transcript counts do not match ",
      "the original processed count table."
    )
  }

  cat("Transcript counts match original processed results.\n")
}

# =========================================================================
# Figure 21F: exact density-contour source data
#
# Uses the same stat_density_2d() settings as the original plot.
# Morphology raster is visual background only and is not needed here.
# =========================================================================
density_list <- list()

for (samp in sample_order) {

  d <- as.data.frame(
    fig21f_pos[
      sample == samp,
      .(x, y)
    ]
  )

  if (nrow(d) == 0) {
    warning("No SMAD3_Pos transcripts for ", samp)
    next
  }

  p <- ggplot(
    d,
    aes(x = x, y = y)
  ) +
    stat_density_2d(
      aes(fill = after_stat(level)),
      geom = "polygon",
      alpha = 0.5,
      contour = TRUE
    ) +
    scale_y_reverse() +
    coord_equal()

  built <- ggplot_build(p)$data[[1]]
  bd <- as.data.table(built)

  bd[, sample := samp]

  density_list[[samp]] <- bd
}

fig21f_density <- rbindlist(
  density_list,
  use.names = TRUE,
  fill = TRUE
)

fwrite(
  fig21f_density,
  file.path(
    srcdir,
    "Figure21F_SMAD3_Pos_density_contour_source_data.csv"
  )
)

# =========================================================================
# Figure 21A-D: original Kcross objects
# =========================================================================
kcross_list <- list()

for (samp in sample_order) {

  k_file <- file.path(
    repo_root,
    "processed_results",
    "14_spatial",
    paste0(
      "Kcross_Smad3_PosTargets_",
      samp,
      ".rds"
    )
  )

  if (!file.exists(k_file)) {
    stop("Missing original Kcross object: ", k_file)
  }

  k_obj <- readRDS(k_file)

  kd <- as.data.table(
    as.data.frame(k_obj)
  )

  kd[, sample := samp]

  setcolorder(
    kd,
    c(
      "sample",
      setdiff(names(kd), "sample")
    )
  )

  kcross_list[[samp]] <- kd
}

fig21ad <- rbindlist(
  kcross_list,
  use.names = TRUE,
  fill = TRUE
)

fwrite(
  fig21ad,
  file.path(
    srcdir,
    "Figure21AD_Kcross_Smad3_PosTargets_source_data.csv"
  )
)

# =========================================================================
# Gene-set definitions
# =========================================================================
gene_sets <- rbindlist(
  list(
    data.table(
      group = "Smad3",
      gene = smad3_gene
    ),
    data.table(
      group = "SMAD3_Pos",
      gene = smad3_pos_genes
    ),
    data.table(
      group = "SMAD3_Neg",
      gene = smad3_neg_genes
    )
  )
)

fwrite(
  gene_sets,
  file.path(
    srcdir,
    "Figure21_SMAD3_gene_sets.csv"
  )
)

# =========================================================================
# Summary
# =========================================================================
cat("\n===== Figure 21 source-data summary =====\n")

cat(
  "Transcript-level rows:",
  nrow(tx_all),
  "\n"
)

cat(
  "SMAD3_Pos rows:",
  nrow(fig21f_pos),
  "\n"
)

cat(
  "Density-contour rows:",
  nrow(fig21f_density),
  "\n"
)

cat(
  "Kcross rows:",
  nrow(fig21ad),
  "\n"
)

cat("\nCounts by sample/group:\n")
print(
  tx_all[
    ,
    .N,
    by = .(sample, group)
  ][
    order(
      match(sample, sample_order),
      match(group, group_order)
    )
  ]
)

cat("\nKcross columns:\n")
print(names(fig21ad))

cat("\nDone.\n")
