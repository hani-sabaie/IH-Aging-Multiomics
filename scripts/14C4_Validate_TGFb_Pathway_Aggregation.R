rm(list = ls(all.names = TRUE))
gc()

suppressPackageStartupMessages({
  library(CellChat)
  library(data.table)
})

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

hist_file <- file.path(
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

if (!file.exists(hist_file)) {
  stop("Historical CellChat object not found.")
}

if (!file.exists(audit_file)) {
  stop("nboot1000 TGFb audit table not found.")
}

x <- readRDS(hist_file)
dt <- fread(audit_file)

datasets <- c("Young", "Aged")

results <- list()

for (ds in datasets) {

  cat("\n============================================================\n")
  cat(ds, "\n")
  cat("============================================================\n\n")

  z <- copy(
    dt[
      dataset == ds
    ]
  )

  # Historical CellChat rule:
  # probability retained when raw CellChat p < 0.05.
  raw_mat <- z[
    p_raw < 0.05,
    .(
      reconstructed_raw_pathway_probability =
        sum(
          probability,
          na.rm = TRUE
        )
    ),
    by = .(
      source,
      target
    )
  ]

  groups <- dimnames(
    x@netP[[ds]]$prob
  )[[1]]

  reconstructed <- matrix(
    0,
    nrow = length(groups),
    ncol = length(groups),
    dimnames = list(
      groups,
      groups
    )
  )

  for (i in seq_len(nrow(raw_mat))) {

    reconstructed[
      raw_mat$source[i],
      raw_mat$target[i]
    ] <- raw_mat$reconstructed_raw_pathway_probability[i]
  }

  pw <- x@netP[[ds]]$pathways
  idx <- match("TGFb", pw)

  if (is.na(idx)) {
    stop("TGFb pathway absent from historical netP: ", ds)
  }

  historical <- x@netP[[ds]]$prob[
    ,
    ,
    idx
  ]

  # Match group order explicitly.
  reconstructed <- reconstructed[
    rownames(historical),
    colnames(historical),
    drop = FALSE
  ]

  delta <- abs(
    historical -
      reconstructed
  )

  max_delta <- max(
    delta,
    na.rm = TRUE
  )

  mean_delta <- mean(
    delta,
    na.rm = TRUE
  )

  exact_match <- isTRUE(
    all.equal(
      as.numeric(historical),
      as.numeric(reconstructed),
      tolerance = 1e-12,
      check.attributes = FALSE
    )
  )

  cat(
    "Max |historical netP - reconstructed| : ",
    format(
      max_delta,
      scientific = TRUE,
      digits = 12
    ),
    "\n",
    sep = ""
  )

  cat(
    "Mean |delta|                          : ",
    format(
      mean_delta,
      scientific = TRUE,
      digits = 12
    ),
    "\n",
    sep = ""
  )

  cat(
    "Exact match within 1e-12              : ",
    exact_match,
    "\n",
    sep = ""
  )

  # Specifically report the C7-affected Young edge.
  if (
    "Vascular stromal" %in% rownames(historical) &&
    "FAP3" %in% colnames(historical)
  ) {

    cat(
      "Historical Vascular stromal -> FAP3 : ",
      historical[
        "Vascular stromal",
        "FAP3"
      ],
      "\n",
      sep = ""
    )

    cat(
      "Reconstructed raw value              : ",
      reconstructed[
        "Vascular stromal",
        "FAP3"
      ],
      "\n",
      sep = ""
    )
  }

  results[[ds]] <- data.table(
    dataset = ds,
    max_abs_delta = max_delta,
    mean_abs_delta = mean_delta,
    exact_match_1e12 = exact_match
  )
}

summary_dt <- rbindlist(results)

cat("\n============================================================\n")
cat("FINAL VALIDATION\n")
cat("============================================================\n\n")

print(summary_dt)

out_file <- file.path(
  repo_root,
  "processed_results",
  "12_CellChat",
  "multiple_testing_audit",
  "TGFb_lightweight",
  "TGFb_pathway_aggregation_validation.tsv"
)

fwrite(
  summary_dt,
  out_file,
  sep = "\t"
)

cat("\nWritten:\n")
cat(out_file, "\n")

cat("\nNo canonical file was modified.\n")
