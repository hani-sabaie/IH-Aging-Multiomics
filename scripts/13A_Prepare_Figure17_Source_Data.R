rm(list = ls(all.names = TRUE))
gc()

library(data.table)
library(dplyr)
library(Seurat)
library(Signac)

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
obj_file <- Sys.getenv(
  "FIG17_CHROMVAR_RDS",
  unset = file.path(
    repo_root,
    "outputs",
    "hdWGCNA_TFNet_DEReg_L2G_chromVAR_obj.rds"
  )
)

if (!file.exists(obj_file)) {
  stop(
    "chromVAR object not found: ",
    obj_file,
    "\nSet FIG17_CHROMVAR_RDS to its local path."
  )
}

source_dir <- file.path(
  repo_root,
  "source_data",
  "figure_source_data",
  "Figure_17"
)

dir.create(
  source_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("Loading chromVAR object...\n")
obj <- readRDS(obj_file)

# -------------------------------------------------------------------------
# Identify the same motif features used in scripts/13_ChromVar.R
# -------------------------------------------------------------------------
motif_obj <- Motifs(obj[["ATAC"]])

motif_names <- motif_obj@motif.names
motif_ids <- names(motif_obj@pwm)

get_motif_ids <- function(tf, motif_names, motif_ids) {
  hits <- grep(tf, motif_names, ignore.case = TRUE)
  motif_ids[hits]
}

get_chromvar_feature <- function(seu, motif_id) {
  rn <- rownames(seu[["chromvar"]])

  cand <- c(
    motif_id,
    paste0("motif_", motif_id)
  )

  hit <- cand[cand %in% rn]

  if (length(hit) == 0) {
    return(NA_character_)
  }

  hit[1]
}

tf_names <- c(
  "ETS2",
  "ETS1",
  "FOS",
  "SMAD3"
)

motif_map <- rbindlist(
  lapply(
    tf_names,
    function(tf) {

      ids <- get_motif_ids(
        tf,
        motif_names,
        motif_ids
      )

      if (length(ids) == 0) {
        stop("No motif ID found for ", tf)
      }

      selected_id <- ids[1]

      feature <- get_chromvar_feature(
        obj,
        selected_id
      )

      if (is.na(feature)) {
        stop(
          "No chromVAR feature found for ",
          tf,
          " / ",
          selected_id
        )
      }

      data.table(
        TF = tf,
        motif_id = selected_id,
        chromVAR_feature = feature
      )
    }
  )
)

fwrite(
  motif_map,
  file.path(
    source_dir,
    "Figure17_motif_feature_mapping.csv"
  )
)

print(motif_map)

# -------------------------------------------------------------------------
# Reconstruct the cell-level data underlying the raincloud plots
# -------------------------------------------------------------------------
fap_levels <- c(
  "FAP1",
  "FAP2",
  "FAP3",
  "FAP4"
)

all_data <- list()
k <- 1L

for (cf in fap_levels) {

  cat("Processing", cf, "...\n")

  cells <- rownames(
    obj@meta.data[
      obj@meta.data$skeletal_muscle == cf &
        obj@meta.data$condition %in% c("Young", "Aged"),
      ,
      drop = FALSE
    ]
  )

  if (length(cells) == 0) {
    warning("No Young/Aged cells found for ", cf)
    next
  }

  obj_cf <- subset(
    obj,
    cells = cells
  )

  DefaultAssay(obj_cf) <- "chromvar"

  for (tf in tf_names) {

    feat <- motif_map[
      TF == tf,
      chromVAR_feature
    ][1]

    df <- FetchData(
      obj_cf,
      vars = c("condition", feat)
    )

    cur <- data.table(
      cell = rownames(df),
      condition = as.character(df$condition),
      FAP = cf,
      TF = tf,
      feature = feat,
      value = as.numeric(df[[feat]])
    )

    all_data[[k]] <- cur
    k <- k + 1L
  }
}

df_plot <- rbindlist(all_data)

df_plot[
  ,
  condition := factor(
    condition,
    levels = c("Young", "Aged")
  )
]

df_plot[
  ,
  TF := factor(
    TF,
    levels = c("ETS2", "ETS1", "FOS", "SMAD3")
  )
]

df_plot[
  ,
  FAP := factor(
    FAP,
    levels = c("FAP1", "FAP2", "FAP3", "FAP4")
  )
]

fwrite(
  df_plot,
  file.path(
    source_dir,
    "Figure17_chromVAR_cell_level_source_data.csv"
  )
)

# -------------------------------------------------------------------------
# Young vs Aged Wilcoxon statistics used in the figure
# -------------------------------------------------------------------------
stats_list <- list()
k <- 1L

for (cf in fap_levels) {

  for (tf in tf_names) {

    cur <- df_plot[
      as.character(FAP) == cf &
        as.character(TF) == tf
    ]

    young <- cur[
      condition == "Young",
      value
    ]

    aged <- cur[
      condition == "Aged",
      value
    ]

    wt <- suppressWarnings(
      wilcox.test(
        young,
        aged
      )
    )

    stats_list[[k]] <- data.table(
      FAP = cf,
      TF = tf,
      n_Young = length(young),
      n_Aged = length(aged),
      mean_Young = mean(young, na.rm = TRUE),
      mean_Aged = mean(aged, na.rm = TRUE),
      median_Young = median(young, na.rm = TRUE),
      median_Aged = median(aged, na.rm = TRUE),
      W = unname(wt$statistic),
      p_value = wt$p.value
    )

    k <- k + 1L
  }
}

wilcox_stats <- rbindlist(stats_list)

wilcox_stats[
  ,
  p_significance := fifelse(
    p_value <= 0.0001, "****",
    fifelse(
      p_value <= 0.001, "***",
      fifelse(
        p_value <= 0.01, "**",
        fifelse(
          p_value <= 0.05, "*",
          "ns"
        )
      )
    )
  )
]

fwrite(
  wilcox_stats,
  file.path(
    source_dir,
    "Figure17_Wilcoxon_Young_vs_Aged.csv"
  )
)

# -------------------------------------------------------------------------
# Summary
# -------------------------------------------------------------------------
cat("\nFigure 17 source data generated.\n")
cat("Cell-level rows:", nrow(df_plot), "\n")
cat("Expected FAP x TF comparisons: 16\n")
cat("Wilcoxon rows:", nrow(wilcox_stats), "\n")

cat("\nRows by FAP x TF:\n")
print(
  df_plot[
    ,
    .N,
    by = .(FAP, TF)
  ]
)

cat("\nDone.\n")
