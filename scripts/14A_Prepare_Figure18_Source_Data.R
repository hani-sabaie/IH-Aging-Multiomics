rm(list = ls(all.names = TRUE))
gc()

library(CellChat)
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
obj_file <- Sys.getenv(
  "FIG18_CELLCHAT_RDS",
  unset = file.path(
    repo_root,
    "processed_results",
    "12_CellChat",
    "cellchat_merged.rds"
  )
)

if (!file.exists(obj_file)) {
  stop(
    "CellChat merged object not found: ",
    obj_file,
    "\nSet FIG18_CELLCHAT_RDS to its local path."
  )
}

srcdir <- Sys.getenv(
  "FIG18_SOURCE_DIR",
  unset = file.path(
    repo_root,
    "source_data",
    "figure_source_data",
    "Figure_18"
  )
)

dir.create(
  srcdir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("Loading CellChat object...\n")
x <- readRDS(obj_file)

datasets <- c("Young", "Aged")

if (!all(datasets %in% names(x@net))) {
  stop("Expected Young and Aged datasets were not found in @net.")
}

if (!all(datasets %in% names(x@netP))) {
  stop("Expected Young and Aged datasets were not found in @netP.")
}

group_order <- c(
  "FAP3",
  "Other FAPs",
  "MuSC",
  "Myogenic",
  "Immune",
  "Endothelial",
  "Vascular stromal"
)

# -------------------------------------------------------------------------
# Helper: matrix -> long format
# -------------------------------------------------------------------------
matrix_long <- function(mat, condition, value_name) {

  dt <- as.data.table(as.table(mat))

  if (ncol(dt) != 3L) {
    stop(
      "Unexpected matrix-to-long structure: ",
      paste(names(dt), collapse = ", ")
    )
  }

  setnames(
    dt,
    names(dt),
    c("source", "target", value_name)
  )

  dt[, condition := condition]

  dt[
    ,
    source_order := match(source, group_order)
  ]

  dt[
    ,
    target_order := match(target, group_order)
  ]

  dt[]
}

# =========================================================================
# Figure 18A
# Overall signaling networks
# =========================================================================
fig18a_list <- list()

for (ds in datasets) {

  weight <- x@net[[ds]]$weight
  count  <- x@net[[ds]]$count

  go <- intersect(group_order, rownames(weight))

  weight <- weight[go, go, drop = FALSE]
  count  <- count[go, go, drop = FALSE]

  dt_w <- matrix_long(
    weight,
    ds,
    "interaction_weight"
  )

  dt_c <- matrix_long(
    count,
    ds,
    "interaction_count"
  )

  dt <- merge(
    dt_w,
    dt_c[
      ,
      .(
        condition,
        source,
        target,
        interaction_count
      )
    ],
    by = c("condition", "source", "target"),
    all.x = TRUE
  )

  vertex_weight <- rowSums(weight)

  dt[
    ,
    source_vertex_weight :=
      vertex_weight[as.character(source)]
  ]

  dt[
    ,
    source_order := match(source, group_order)
  ]

  dt[
    ,
    target_order := match(target, group_order)
  ]

  fig18a_list[[ds]] <- dt
}

fig18a <- rbindlist(fig18a_list, use.names = TRUE)

setorder(
  fig18a,
  condition,
  source_order,
  target_order
)

fwrite(
  fig18a,
  file.path(
    srcdir,
    "Figure18A_overall_signaling_network_source_data.csv"
  )
)

# =========================================================================
# Figure 18B
# Signaling-role scatter
# x = summed outgoing strength
# y = summed incoming strength
# dot size = total number of inferred links
# =========================================================================
fig18b_list <- list()

for (ds in datasets) {

  centr <- x@netP[[ds]]$centr

  if (length(centr) == 0) {
    stop("Centrality data missing for ", ds)
  }

  groups <- names(centr[[1]]$outdeg)

  outgoing_mat <- sapply(
    centr,
    function(z) z$outdeg[groups]
  )

  incoming_mat <- sapply(
    centr,
    function(z) z$indeg[groups]
  )

  if (is.null(dim(outgoing_mat))) {
    outgoing_mat <- matrix(
      outgoing_mat,
      ncol = 1,
      dimnames = list(groups, names(centr))
    )
  }

  if (is.null(dim(incoming_mat))) {
    incoming_mat <- matrix(
      incoming_mat,
      ncol = 1,
      dimnames = list(groups, names(centr))
    )
  }

  count_mat <- x@net[[ds]]$count

  num_link <- (
    rowSums(count_mat) +
      colSums(count_mat) -
      diag(count_mat)
  )

  dt <- data.table(
    condition = ds,
    cell_group = groups,
    outgoing_interaction_strength =
      rowSums(outgoing_mat, na.rm = TRUE),
    incoming_interaction_strength =
      rowSums(incoming_mat, na.rm = TRUE),
    number_of_links =
      as.numeric(num_link[groups])
  )

  dt[
    ,
    group_order := match(cell_group, group_order)
  ]

  fig18b_list[[ds]] <- dt
}

fig18b <- rbindlist(fig18b_list)

fig18b[
  ,
  dot_size_global_min := min(number_of_links, na.rm = TRUE)
]

fig18b[
  ,
  dot_size_global_max := max(number_of_links, na.rm = TRUE)
]

setorder(fig18b, condition, group_order)

fwrite(
  fig18b,
  file.path(
    srcdir,
    "Figure18B_signaling_role_scatter_source_data.csv"
  )
)

# =========================================================================
# Figure 18C
# Incoming / outgoing signaling-role heatmaps
#
# CellChat row-scaling:
# each pathway divided by its maximum score across cell groups
# =========================================================================
pathway_union <- union(
  x@netP$Young$pathways,
  x@netP$Aged$pathways
)

fig18c_list <- list()
k <- 1L

for (ds in datasets) {

  centr <- x@netP[[ds]]$centr

  groups <- names(centr[[1]]$outdeg)

  for (pattern in c("outgoing", "incoming")) {

    measure <- if (
      pattern == "outgoing"
    ) {
      "outdeg"
    } else {
      "indeg"
    }

    for (pw in pathway_union) {

      if (pw %in% names(centr)) {

        score <- centr[[pw]][[measure]][groups]

      } else {

        score <- setNames(
          rep(0, length(groups)),
          groups
        )
      }

      mx <- max(score, na.rm = TRUE)

      if (is.finite(mx) && mx > 0) {
        scaled <- score / mx
      } else {
        scaled <- rep(0, length(score))
      }

      fig18c_list[[k]] <- data.table(
        condition = ds,
        pattern = pattern,
        pathway = pw,
        cell_group = groups,
        raw_score = as.numeric(score),
        row_scaled_score = as.numeric(scaled)
      )

      k <- k + 1L
    }
  }
}

fig18c <- rbindlist(fig18c_list)

fig18c[
  ,
  cell_total_raw := sum(raw_score, na.rm = TRUE),
  by = .(condition, pattern, cell_group)
]

fig18c[
  ,
  pathway_total_raw := sum(raw_score, na.rm = TRUE),
  by = .(condition, pattern, pathway)
]

fig18c[
  ,
  pathway_order := match(pathway, pathway_union)
]

fig18c[
  ,
  group_order := match(cell_group, group_order)
]

setorder(
  fig18c,
  condition,
  pattern,
  pathway_order,
  group_order
)

fwrite(
  fig18c,
  file.path(
    srcdir,
    "Figure18C_signaling_role_heatmap_source_data.csv"
  )
)

# =========================================================================
# Figure 18D
# TGFb pathway circle network
# =========================================================================
fig18d_list <- list()

for (ds in datasets) {

  pw <- x@netP[[ds]]$pathways
  idx <- match("TGFb", pw)

  if (is.na(idx)) {
    stop("TGFb pathway not found for ", ds)
  }

  prob_mat <- x@netP[[ds]]$prob[, , idx]

  go <- intersect(group_order, rownames(prob_mat))

  prob_mat <- prob_mat[
    go,
    go,
    drop = FALSE
  ]

  dt <- matrix_long(
    prob_mat,
    ds,
    "TGFb_pathway_probability"
  )

  fig18d_list[[ds]] <- dt
}

fig18d <- rbindlist(fig18d_list)

setorder(
  fig18d,
  condition,
  source_order,
  target_order
)

fwrite(
  fig18d,
  file.path(
    srcdir,
    "Figure18D_TGFb_pathway_network_source_data.csv"
  )
)

# =========================================================================
# Figure 18E
# TGFb ligand-receptor contribution in Aged
#
# Reproduces the normalization used by netAnalysis_contribution()
# =========================================================================
ds <- "Aged"

prob <- x@net[[ds]]$prob
pval <- x@net[[ds]]$pval

lr_names <- dimnames(prob)[[3]]

if (
  is.null(lr_names) ||
    length(lr_names) != dim(prob)[3]
) {
  lr_names <- x@net[[ds]]$LRs
  dimnames(prob)[[3]] <- lr_names
  dimnames(pval)[[3]] <- lr_names
}

db <- x@DB$interaction

tgfb_lr <- db[
  db$pathway_name == "TGFb",
  "interaction_name"
]

tgfb_lr <- intersect(
  as.character(tgfb_lr),
  lr_names
)

if (length(tgfb_lr) == 0) {
  stop("No TGFb ligand-receptor pairs found in Aged network.")
}

# CellChat contribution threshold
prob_sig <- prob
prob_sig[pval > 0.05] <- 0

lr_sum <- sapply(
  tgfb_lr,
  function(lr) {
    sum(prob_sig[, , lr], na.rm = TRUE)
  }
)

tgfb_lr_use <- names(lr_sum)[lr_sum != 0]

if (length(tgfb_lr_use) == 0) {
  stop("No significant TGFb communications found in Aged.")
}

prob_sel <- prob_sig[
  ,
  ,
  tgfb_lr_use,
  drop = FALSE
]

pmin_all <- min(prob_sel, na.rm = TRUE)
pmax_all <- max(prob_sel, na.rm = TRUE)

if (
  is.finite(pmax_all) &&
    is.finite(pmin_all) &&
    pmax_all > pmin_all
) {

  prob_scaled <- (
    prob_sel - pmin_all
  ) / (
    pmax_all - pmin_all
  )

} else {

  prob_scaled <- prob_sel * 0
}

scaled_sum <- apply(
  prob_scaled,
  3,
  sum,
  na.rm = TRUE
)

total_scaled <- sum(
  prob_scaled,
  na.rm = TRUE
)

if (total_scaled > 0) {
  contribution <- scaled_sum / total_scaled
} else {
  contribution <- rep(0, length(scaled_sum))
}

raw_sum <- apply(
  prob_sel,
  3,
  sum,
  na.rm = TRUE
)

n_sig_edges <- apply(
  prob_sel > 0,
  3,
  sum,
  na.rm = TRUE
)

db_idx <- match(
  tgfb_lr_use,
  db$interaction_name
)

fig18e <- data.table(
  condition = ds,
  interaction_name = tgfb_lr_use,
  interaction_name_2 =
    db$interaction_name_2[db_idx],
  ligand =
    db$ligand[db_idx],
  receptor =
    db$receptor[db_idx],
  raw_significant_probability_sum =
    as.numeric(raw_sum[tgfb_lr_use]),
  scaled_probability_sum =
    as.numeric(scaled_sum[tgfb_lr_use]),
  relative_contribution =
    as.numeric(contribution[tgfb_lr_use]),
  number_of_significant_cellpair_edges =
    as.numeric(n_sig_edges[tgfb_lr_use])
)

setorder(
  fig18e,
  -relative_contribution
)

fwrite(
  fig18e,
  file.path(
    srcdir,
    "Figure18E_TGFb_LR_contribution_Aged_source_data.csv"
  )
)

# =========================================================================
# Figure 18F
# Exact data used by the TGFb comparison bubble plot
# =========================================================================
cat("Extracting Figure 18F bubble-plot data...\n")

p_bubble <- netVisual_bubble(
  object = x,
  sources.use = c(
    "FAP3",
    "Other FAPs"
  ),
  targets.use = c(
    "MuSC",
    "Myogenic",
    "Immune",
    "Endothelial",
    "Vascular stromal"
  ),
  signaling = "TGFb",
  angle.x = 45,
  remove.isolate = TRUE,
  comparison = c(1, 2)
)

fig18f <- as.data.table(
  as.data.frame(p_bubble$data)
)

fig18f[
  ,
  figure_panel := "Figure18F"
]

# Add the communication probability and the revised inferential P value
# directly from @net. In the revised canonical CellChat object, @net$pval
# stores plus-one empirical P values after condition-wide BH correction.
fig18f[
  ,
  CellChat_probability := mapply(
    function(ds, src, tgt, lr) {
      x@net[[ds]]$prob[
        as.character(src),
        as.character(tgt),
        as.character(lr)
      ]
    },
    dataset,
    source,
    target,
    interaction_name
  )
]

fig18f[
  ,
  plus1_BH_adjusted_p_value := mapply(
    function(ds, src, tgt, lr) {
      x@net[[ds]]$pval[
        as.character(src),
        as.character(tgt),
        as.character(lr)
      ]
    },
    dataset,
    source,
    target,
    interaction_name
  )
]

fig18f[
  ,
  significant_plus1_BH_p_lt_0_05 :=
    plus1_BH_adjusted_p_value < 0.05
]

fwrite(
  fig18f,
  file.path(
    srcdir,
    "Figure18F_TGFb_bubble_source_data.csv"
  )
)

# =========================================================================
# Summary
# =========================================================================
cat("\nFigure 18 source data generated.\n")
cat("Figure18A rows:", nrow(fig18a), "\n")
cat("Figure18B rows:", nrow(fig18b), "\n")
cat("Figure18C rows:", nrow(fig18c), "\n")
cat("Figure18D rows:", nrow(fig18d), "\n")
cat("Figure18E rows:", nrow(fig18e), "\n")
cat("Figure18F rows:", nrow(fig18f), "\n")

cat("\nPathways used in Figure18C:\n")
print(pathway_union)

cat("\nTGFb contribution table:\n")
print(fig18e)

cat("\nDone.\n")
