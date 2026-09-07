# ============================================================================
# Reviewer C7 audit:
# Rebuild corrected full-18 CellChat networks from the completed
# condition-wide nboot=1000 + plus-one + BH results, then quantify effects
# on Figure 18A-C source data.
#
# No bootstrap inference is performed here.
# No canonical CellChat object or Figure 18 source-data file is modified.
# ============================================================================

rm(list = ls(all.names = TRUE))

suppressPackageStartupMessages({
  library(CellChat)
  library(data.table)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (length(file_arg) == 1L) {
  script_dir <- dirname(
    normalizePath(sub("^--file=", "", file_arg))
  )
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

audit_root <- file.path(
  repo_root,
  "processed_results",
  "12_CellChat",
  "multiple_testing_audit",
  "full18_nboot1000_BH"
)

outdir <- file.path(
  audit_root,
  "figure18ABC_corrected_candidate"
)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

datasets <- c("Young", "Aged")

test_files <- setNames(
  file.path(
    audit_root,
    paste0(
      datasets,
      "_full18_nboot1000_conditionwide_BH_tests.tsv"
    )
  ),
  datasets
)

required_files <- c(
  hist_file,
  unname(test_files)
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0L) {
  stop(
    "Required file(s) missing:\n",
    paste(missing_files, collapse = "\n")
  )
}

cat("Loading canonical historical CellChat object...\n")
x <- readRDS(hist_file)

# ============================================================================
# Construct one corrected condition-specific CellChat object
# ============================================================================

build_corrected_condition <- function(ds) {

  cat("\n============================================================\n")
  cat("REBUILDING CORRECTED NETWORK:", ds, "\n")
  cat("============================================================\n\n")

  z <- fread(test_files[[ds]])

  expected_n_lr <- if (ds == "Young") 488L else 468L
  expected_n_tests <- if (ds == "Young") 23912L else 22932L
  expected_n_sig <- if (ds == "Young") 680L else 702L

  required_cols <- c(
    "dataset",
    "LR_family_index",
    "interaction_name",
    "source",
    "target",
    "probability",
    "p_plus1_BH_condition"
  )

  missing_cols <- setdiff(
    required_cols,
    names(z)
  )

  if (length(missing_cols) > 0L) {
    stop(
      ds,
      " corrected test table missing columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  if (nrow(z) != expected_n_tests) {
    stop(
      ds,
      ": expected ",
      expected_n_tests,
      " corrected tests; found ",
      nrow(z),
      "."
    )
  }

  if (!all(z$dataset == ds)) {
    stop(ds, ": unexpected dataset label in corrected test table.")
  }

  historical_prob <- x@net[[ds]]$prob

  groups <- dimnames(historical_prob)[[1]]
  targets <- dimnames(historical_prob)[[2]]
  lr_order <- dimnames(historical_prob)[[3]]

  if (
    length(groups) != 7L ||
    length(targets) != 7L ||
    length(lr_order) != expected_n_lr
  ) {
    stop(
      ds,
      ": unexpected historical probability dimensions."
    )
  }

  if (!identical(groups, targets)) {
    stop(ds, ": source and target group order differs.")
  }

  # ------------------------------------------------------------------------
  # Validate full corrected hypothesis family
  # ------------------------------------------------------------------------

  src_i <- match(z$source, groups)
  tgt_i <- match(z$target, groups)
  lr_i <- match(z$interaction_name, lr_order)

  if (
    anyNA(src_i) ||
    anyNA(tgt_i) ||
    anyNA(lr_i)
  ) {
    stop(
      ds,
      ": corrected tests contain an unknown source/target/LR."
    )
  }

  formal_keys <- paste(
    z$interaction_name,
    z$source,
    z$target,
    sep = "|"
  )

  if (anyDuplicated(formal_keys) > 0L) {
    stop(ds, ": duplicate corrected formal hypotheses detected.")
  }

  # Probability values in the completed batch audit reproduced historical
  # probabilities exactly. Confirm this again before reconstruction.
  hist_prob_from_keys <- historical_prob[
    cbind(
      src_i,
      tgt_i,
      lr_i
    )
  ]

  max_prob_delta <- max(
    abs(
      hist_prob_from_keys -
        z$probability
    ),
    na.rm = TRUE
  )

  if (
    !isTRUE(
      all.equal(
        hist_prob_from_keys,
        z$probability,
        tolerance = 1e-12,
        check.attributes = FALSE
      )
    )
  ) {
    stop(
      ds,
      ": corrected test probabilities no longer match historical values."
    )
  }

  cat(
    "Max historical/test probability delta:",
    format(max_prob_delta, scientific = TRUE, digits = 12),
    "\n"
  )

  # ------------------------------------------------------------------------
  # Reconstruct corrected adjusted-P array
  # ------------------------------------------------------------------------

  corrected_p <- array(
    NA_real_,
    dim = dim(historical_prob),
    dimnames = dimnames(historical_prob)
  )

  corrected_p[
    cbind(
      src_i,
      tgt_i,
      lr_i
    )
  ] <- z$p_plus1_BH_condition

  if (anyNA(corrected_p)) {
    stop(
      ds,
      ": corrected P array was not completely populated."
    )
  }

  n_corrected_sig <- sum(
    corrected_p < 0.05,
    na.rm = TRUE
  )

  if (n_corrected_sig != expected_n_sig) {
    stop(
      ds,
      ": expected ",
      expected_n_sig,
      " corrected significant tests; found ",
      n_corrected_sig,
      "."
    )
  }

  cat(
    "Corrected significant tests:",
    n_corrected_sig,
    "\n"
  )

  # ------------------------------------------------------------------------
  # Recover condition-specific metadata and historical signaling matrix
  # ------------------------------------------------------------------------

  cells <- rownames(x@meta)[
    x@meta$datasets == ds
  ]

  if (length(cells) == 0L) {
    stop("No cells found for ", ds)
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

  if (
    !identical(
      colnames(full_data),
      rownames(meta)
    )
  ) {
    stop(
      ds,
      ": signaling matrix and metadata are not aligned."
    )
  }

  meta$group_cellchat <- factor(
    meta$group_cellchat,
    levels = groups
  )

  if (anyNA(meta$group_cellchat)) {
    stop(ds, ": reconstructed group identities contain NA.")
  }

  # ------------------------------------------------------------------------
  # Reconstruct a CellChat object without rerunning inference
  # ------------------------------------------------------------------------

  cc <- createCellChat(
    object = full_data,
    meta = meta,
    group.by = "group_cellchat"
  )

  cc@DB <- x@DB
  cc@data.signaling <- full_data

  cc@idents <- factor(
    meta$group_cellchat,
    levels = groups
  )

  names(cc@idents) <- rownames(meta)

  # Exact historical LR family/order.
  lr_sig <- x@LR[[ds]]$LRsig

  lr_idx <- match(
    lr_order,
    as.character(lr_sig$interaction_name)
  )

  if (anyNA(lr_idx)) {
    stop(
      ds,
      ": historical LR family cannot be reconstructed from @LR$LRsig."
    )
  }

  lr_sig <- lr_sig[
    lr_idx,
    ,
    drop = FALSE
  ]

  rownames(lr_sig) <- lr_order

  cc@LR <- list(
    LRsig = lr_sig
  )

  cc@net$prob <- historical_prob
  cc@net$pval <- corrected_p

  retained_lr <- lr_order[
    apply(
      historical_prob * (corrected_p < 0.05),
      3,
      sum,
      na.rm = TRUE
    ) > 0
  ]

  cc@net$LRs <- retained_lr

  cat(
    "LRs with >=1 corrected significant edge:",
    length(retained_lr),
    "\n"
  )

  # ------------------------------------------------------------------------
  # Recompute pathway aggregation / total network / centrality
  # ------------------------------------------------------------------------

  cc <- computeCommunProbPathway(
    cc
  )

  cc <- aggregateNet(
    cc
  )

  cc <- netAnalysis_computeCentrality(
    cc,
    slot.name = "netP"
  )

  cat(
    "Corrected retained pathways:",
    length(cc@netP$pathways),
    "\n"
  )

  cat(
    "Corrected total aggregate count:",
    sum(cc@net$count, na.rm = TRUE),
    "\n"
  )

  list(
    object = cc,
    summary = data.table(
      condition = ds,
      n_tests = nrow(z),
      n_significant_corrected = n_corrected_sig,
      n_retained_LR = length(retained_lr),
      n_historical_pathways =
        length(x@netP[[ds]]$pathways),
      n_corrected_pathways =
        length(cc@netP$pathways),
      historical_total_count =
        sum(x@net[[ds]]$count, na.rm = TRUE),
      corrected_total_count =
        sum(cc@net$count, na.rm = TRUE),
      historical_total_weight =
        sum(x@net[[ds]]$weight, na.rm = TRUE),
      corrected_total_weight =
        sum(cc@net$weight, na.rm = TRUE),
      max_probability_delta =
        max_prob_delta
    )
  )
}

young_result <- build_corrected_condition(
  "Young"
)

aged_result <- build_corrected_condition(
  "Aged"
)

object.list <- list(
  Young = young_result$object,
  Aged = aged_result$object
)

corrected_merged <- mergeCellChat(
  object.list = object.list,
  add.names = names(object.list)
)

network_summary <- rbindlist(
  list(
    young_result$summary,
    aged_result$summary
  )
)

saveRDS(
  young_result$object,
  file.path(
    outdir,
    "cellchat_Young_full18_nboot1000_plus1BH_corrected_candidate.rds"
  )
)

saveRDS(
  aged_result$object,
  file.path(
    outdir,
    "cellchat_Aged_full18_nboot1000_plus1BH_corrected_candidate.rds"
  )
)

saveRDS(
  corrected_merged,
  file.path(
    outdir,
    "cellchat_merged_full18_nboot1000_plus1BH_corrected_candidate.rds"
  )
)

fwrite(
  network_summary,
  file.path(
    outdir,
    "corrected_network_summary.tsv"
  ),
  sep = "\t"
)

# ============================================================================
# Shared source-data extraction helpers
# ============================================================================

matrix_long <- function(
    mat,
    condition,
    value_name
) {

  dt <- as.data.table(
    as.table(mat)
  )

  setnames(
    dt,
    c(
      "source",
      "target",
      value_name
    )
  )

  dt[
    ,
    condition := condition
  ]

  setcolorder(
    dt,
    c(
      "condition",
      "source",
      "target",
      value_name
    )
  )

  dt[]
}

# Use a single shared pathway universe for historical/corrected comparison.
pathway_union <- unique(
  c(
    x@netP$Young$pathways,
    x@netP$Aged$pathways,
    corrected_merged@netP$Young$pathways,
    corrected_merged@netP$Aged$pathways
  )
)

# ============================================================================
# Figure 18A extraction
# ============================================================================

extract_A <- function(obj) {

  out <- list()

  for (ds in datasets) {

    weight <- obj@net[[ds]]$weight
    count <- obj@net[[ds]]$count

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
      dt_c,
      by = c(
        "condition",
        "source",
        "target"
      ),
      all = TRUE
    )

    vertex_weight <- rowSums(
      weight,
      na.rm = TRUE
    )

    dt[
      ,
      source_vertex_weight :=
        vertex_weight[
          as.character(source)
        ]
    ]

    out[[ds]] <- dt
  }

  rbindlist(
    out,
    use.names = TRUE
  )
}

# ============================================================================
# Figure 18B extraction
# ============================================================================

extract_B <- function(obj) {

  out <- list()

  for (ds in datasets) {

    centr <- obj@netP[[ds]]$centr

    if (length(centr) == 0L) {
      stop(
        "Centrality data missing for ",
        ds
      )
    }

    groups <- names(
      centr[[1]]$outdeg
    )

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
        ncol = 1L,
        dimnames = list(
          groups,
          names(centr)
        )
      )
    }

    if (is.null(dim(incoming_mat))) {
      incoming_mat <- matrix(
        incoming_mat,
        ncol = 1L,
        dimnames = list(
          groups,
          names(centr)
        )
      )
    }

    count_mat <- obj@net[[ds]]$count

    num_link <- (
      rowSums(count_mat) +
        colSums(count_mat) -
        diag(count_mat)
    )

    out[[ds]] <- data.table(
      condition = ds,
      cell_group = groups,
      outgoing_interaction_strength =
        rowSums(
          outgoing_mat,
          na.rm = TRUE
        ),
      incoming_interaction_strength =
        rowSums(
          incoming_mat,
          na.rm = TRUE
        ),
      number_of_links =
        as.numeric(
          num_link[groups]
        )
    )
  }

  rbindlist(out)
}

# ============================================================================
# Figure 18C extraction
# ============================================================================

extract_C <- function(obj) {

  out <- list()
  k <- 1L

  for (ds in datasets) {

    centr <- obj@netP[[ds]]$centr

    groups <- names(
      centr[[1]]$outdeg
    )

    for (
      pattern in c(
        "outgoing",
        "incoming"
      )
    ) {

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
            rep(
              0,
              length(groups)
            ),
            groups
          )
        }

        mx <- max(
          score,
          na.rm = TRUE
        )

        if (
          is.finite(mx) &&
          mx > 0
        ) {
          scaled <- score / mx
        } else {
          scaled <- rep(
            0,
            length(score)
          )
        }

        out[[k]] <- data.table(
          condition = ds,
          pattern = pattern,
          pathway = pw,
          cell_group = groups,
          raw_score =
            as.numeric(score),
          row_scaled_score =
            as.numeric(scaled)
        )

        k <- k + 1L
      }
    }
  }

  rbindlist(out)
}

# ============================================================================
# Historical vs corrected source data
# ============================================================================

hist_A <- extract_A(x)
corr_A <- extract_A(corrected_merged)

hist_B <- extract_B(x)
corr_B <- extract_B(corrected_merged)

hist_C <- extract_C(x)
corr_C <- extract_C(corrected_merged)

# --------------------------------------------------------------------------
# Figure 18A comparison
# --------------------------------------------------------------------------

cmp_A <- merge(
  hist_A,
  corr_A,
  by = c(
    "condition",
    "source",
    "target"
  ),
  suffixes = c(
    "_historical",
    "_corrected"
  ),
  all = TRUE
)

cmp_A[
  ,
  `:=`(
    delta_interaction_weight =
      interaction_weight_corrected -
      interaction_weight_historical,

    delta_interaction_count =
      interaction_count_corrected -
      interaction_count_historical,

    delta_source_vertex_weight =
      source_vertex_weight_corrected -
      source_vertex_weight_historical
  )
]

cmp_A[
  ,
  changed :=
    abs(delta_interaction_weight) > 1e-12 |
    abs(delta_interaction_count) > 0 |
    abs(delta_source_vertex_weight) > 1e-12
]

# --------------------------------------------------------------------------
# Figure 18B comparison
# --------------------------------------------------------------------------

cmp_B <- merge(
  hist_B,
  corr_B,
  by = c(
    "condition",
    "cell_group"
  ),
  suffixes = c(
    "_historical",
    "_corrected"
  ),
  all = TRUE
)

cmp_B[
  ,
  `:=`(
    delta_outgoing =
      outgoing_interaction_strength_corrected -
      outgoing_interaction_strength_historical,

    delta_incoming =
      incoming_interaction_strength_corrected -
      incoming_interaction_strength_historical,

    delta_number_of_links =
      number_of_links_corrected -
      number_of_links_historical
  )
]

cmp_B[
  ,
  changed :=
    abs(delta_outgoing) > 1e-12 |
    abs(delta_incoming) > 1e-12 |
    abs(delta_number_of_links) > 0
]

# --------------------------------------------------------------------------
# Figure 18C comparison
# --------------------------------------------------------------------------

cmp_C <- merge(
  hist_C,
  corr_C,
  by = c(
    "condition",
    "pattern",
    "pathway",
    "cell_group"
  ),
  suffixes = c(
    "_historical",
    "_corrected"
  ),
  all = TRUE
)

cmp_C[
  ,
  `:=`(
    delta_raw_score =
      raw_score_corrected -
      raw_score_historical,

    delta_row_scaled_score =
      row_scaled_score_corrected -
      row_scaled_score_historical
  )
]

cmp_C[
  ,
  changed :=
    abs(delta_raw_score) > 1e-12 |
    abs(delta_row_scaled_score) > 1e-12
]

# ============================================================================
# Panel-level summary
# ============================================================================

panel_summary <- rbindlist(
  list(
    data.table(
      panel = "Figure18A",
      n_rows = nrow(cmp_A),
      n_changed_rows =
        sum(cmp_A$changed, na.rm = TRUE),
      max_abs_primary_delta =
        max(
          abs(
            cmp_A$delta_interaction_weight
          ),
          na.rm = TRUE
        )
    ),

    data.table(
      panel = "Figure18B",
      n_rows = nrow(cmp_B),
      n_changed_rows =
        sum(cmp_B$changed, na.rm = TRUE),
      max_abs_primary_delta =
        max(
          c(
            abs(cmp_B$delta_outgoing),
            abs(cmp_B$delta_incoming)
          ),
          na.rm = TRUE
        )
    ),

    data.table(
      panel = "Figure18C",
      n_rows = nrow(cmp_C),
      n_changed_rows =
        sum(cmp_C$changed, na.rm = TRUE),
      max_abs_primary_delta =
        max(
          abs(
            cmp_C$delta_row_scaled_score
          ),
          na.rm = TRUE
        )
    )
  )
)

# ============================================================================
# Pathway-level corrected-significance summary
# ============================================================================

db <- x@DB$interaction

sig_by_pathway_list <- list()

for (ds in datasets) {

  z <- fread(
    test_files[[ds]]
  )

  idx <- match(
    z$interaction_name,
    as.character(
      db$interaction_name
    )
  )

  if (anyNA(idx)) {
    stop(
      ds,
      ": unable to map some LR interactions to CellChat DB."
    )
  }

  z[
    ,
    pathway :=
      as.character(
        db$pathway_name[idx]
      )
  ]

  sig_by_pathway_list[[ds]] <- z[
    ,
    .(
      n_tests = .N,
      n_corrected_significant =
        sum(
          p_plus1_BH_condition < 0.05
        )
    ),
    by = .(
      condition = dataset,
      pathway
    )
  ]
}

sig_by_pathway <- rbindlist(
  sig_by_pathway_list
)

setorder(
  sig_by_pathway,
  condition,
  pathway
)

# ============================================================================
# Write audit outputs
# ============================================================================

fwrite(
  corr_A,
  file.path(
    outdir,
    "Figure18A_corrected_candidate_source_data.csv"
  )
)

fwrite(
  corr_B,
  file.path(
    outdir,
    "Figure18B_corrected_candidate_source_data.csv"
  )
)

fwrite(
  corr_C,
  file.path(
    outdir,
    "Figure18C_corrected_candidate_source_data.csv"
  )
)

fwrite(
  cmp_A,
  file.path(
    outdir,
    "Figure18A_historical_vs_corrected.csv"
  )
)

fwrite(
  cmp_B,
  file.path(
    outdir,
    "Figure18B_historical_vs_corrected.csv"
  )
)

fwrite(
  cmp_C,
  file.path(
    outdir,
    "Figure18C_historical_vs_corrected.csv"
  )
)

fwrite(
  cmp_A[changed == TRUE],
  file.path(
    outdir,
    "Figure18A_changed_rows.csv"
  )
)

fwrite(
  cmp_B[changed == TRUE],
  file.path(
    outdir,
    "Figure18B_changed_rows.csv"
  )
)

fwrite(
  cmp_C[changed == TRUE],
  file.path(
    outdir,
    "Figure18C_changed_rows.csv"
  )
)

fwrite(
  panel_summary,
  file.path(
    outdir,
    "Figure18ABC_change_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  sig_by_pathway,
  file.path(
    outdir,
    "corrected_significant_tests_by_pathway.tsv"
  ),
  sep = "\t"
)

cat("\n============================================================\n")
cat("CORRECTED NETWORK SUMMARY\n")
cat("============================================================\n\n")

print(network_summary)

cat("\n============================================================\n")
cat("FIGURE 18A-C CHANGE SUMMARY\n")
cat("============================================================\n\n")

print(panel_summary)

cat("\n============================================================\n")
cat("CORRECTED SIGNIFICANT TESTS BY PATHWAY\n")
cat("============================================================\n\n")

print(sig_by_pathway)

cat("\nAudit candidate directory:\n")
cat(outdir, "\n")

cat(
  "\nNo canonical CellChat object or Figure 18 source-data file was modified.\n"
)
