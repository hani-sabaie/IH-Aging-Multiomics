# ===== SMAD3 LinkPeaks multiple-testing audit =====
#
# Purpose:
#   Audit the reported SMAD3 peak-to-gene linkage for reviewer C7 without
#   rerunning the stochastic LinkPeaks background-sampling procedure.
#
# Strategy:
#   1) Recover the original SMAD3 links from the saved post-LinkPeaks object.
#   2) Identify the reported peak chr15:67109227-67110381.
#   3) Reconstruct the full number of cis peaks eligible for testing under
#      the original LinkPeaks settings:
#         distance = 500 kb
#         min.cells = 10
#   4) Apply a conservative Bonferroni correction across ALL eligible
#      SMAD3 cis peaks.
#   5) Also calculate an intentionally ultra-conservative correction across
#      every ATAC peak in the object.
#
# No canonical analysis object or manuscript file is modified.

rm(list = ls(all.names = TRUE))
gc()

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(Matrix)
  library(GenomicRanges)
  library(data.table)
})

# -------------------------------------------------------------------------
# Repository paths
# -------------------------------------------------------------------------

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (length(file_arg) == 1) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg)))
  repo_root <- normalizePath(file.path(script_dir, ".."))
} else {
  repo_root <- normalizePath(".")
}

# Resolve large intermediate RDS files.
#
# Prefer repository-local outputs when available. For the original analysis,
# the large intermediate objects were retained on the original F: drive and
# were not copied into the reproducibility repository.

original_outputs_dir <- "F:/Hani's Files/Hernia/outputs"

resolve_rds <- function(filename) {

  candidates <- c(
    file.path(repo_root, "outputs", filename),
    file.path(original_outputs_dir, filename)
  )

  hit <- candidates[file.exists(candidates)]

  if (length(hit) == 0L) {
    stop(
      "Required RDS object not found. Checked:\n  ",
      paste(candidates, collapse = "\n  ")
    )
  }

  normalizePath(hit[1], winslash = "/", mustWork = TRUE)
}

pre_link_file <- resolve_rds(
  "hdWGCNA_TFNet_obj.rds"
)

post_link_file <- resolve_rds(
  "hdWGCNA_TFNet_DEReg_L2G_obj.rds"
)

outdir <- file.path(
  repo_root,
  "processed_results",
  "10_TF_network",
  "multiple_testing_audit"
)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cat("Pre-LinkPeaks object : ", pre_link_file, "\n", sep = "")
cat("Post-LinkPeaks object: ", post_link_file, "\n", sep = "")

# -------------------------------------------------------------------------
# Original analysis settings
# -------------------------------------------------------------------------

goi <- "SMAD3"
distance_bp <- 5e5
min_cells <- 10L

reported_chr <- "chr15"
reported_start <- 67109227L
reported_end <- 67110381L
reported_peak_pattern <- "67109227.*67110381"

# -------------------------------------------------------------------------
# 1. Recover historical LinkPeaks result
# -------------------------------------------------------------------------

cat("\n============================================================\n")
cat("RECOVERING ORIGINAL SMAD3 LINKPEAKS RESULT\n")
cat("============================================================\n\n")

post_obj <- readRDS(post_link_file)
DefaultAssay(post_obj) <- "ATAC"

historical_links <- Links(post_obj)

if (length(historical_links) == 0L) {
  stop("No LinkPeaks results found in the post-LinkPeaks object.")
}

historical_df <- as.data.table(
  as.data.frame(historical_links)
)

required_link_cols <- c(
  "gene",
  "peak",
  "score",
  "zscore",
  "pvalue"
)

missing_link_cols <- setdiff(
  required_link_cols,
  names(historical_df)
)

if (length(missing_link_cols) > 0L) {
  stop(
    "Historical LinkPeaks table is missing columns: ",
    paste(missing_link_cols, collapse = ", ")
  )
}

smad3_links <- historical_df[
  gene == goi
]

if (nrow(smad3_links) == 0L) {
  stop("No historical SMAD3 links found.")
}

reported_link <- smad3_links[
  grepl(reported_peak_pattern, peak)
]

if (nrow(reported_link) != 1L) {
  cat("\nAll historical SMAD3 links:\n")
  print(
    smad3_links[
      ,
      .(gene, peak, score, zscore, pvalue)
    ]
  )

  stop(
    "Expected exactly one historical SMAD3 link matching ",
    reported_chr, ":",
    reported_start, "-", reported_end,
    "; found ", nrow(reported_link), "."
  )
}

raw_p <- reported_link$pvalue[1]
link_score <- reported_link$score[1]
link_z <- reported_link$zscore[1]

cat("Reported peak :", reported_link$peak[1], "\n")
cat("Link score    :", format(link_score, digits = 12), "\n")
cat("Z score       :", format(link_z, digits = 12), "\n")
cat("Raw P value   :", format(raw_p, scientific = TRUE, digits = 12), "\n")
cat("Historical SMAD3 links retained by LinkPeaks:", nrow(smad3_links), "\n")

# Keep only the compact link tables before loading the second multi-GB object.
rm(post_obj, historical_links)
gc()

# -------------------------------------------------------------------------
# 2. Reconstruct all cis peaks eligible under original LinkPeaks settings
# -------------------------------------------------------------------------

cat("\n============================================================\n")
cat("RECONSTRUCTING SMAD3 CIS TEST FAMILY\n")
cat("============================================================\n\n")

pre_obj <- readRDS(pre_link_file)
DefaultAssay(pre_obj) <- "ATAC"

if (!goi %in% rownames(pre_obj[["SCT"]])) {
  stop("SMAD3 not found in the SCT assay.")
}

peak_counts <- GetAssayData(
  object = pre_obj,
  assay = "ATAC",
  slot = "counts"
)

gene_data <- GetAssayData(
  object = pre_obj,
  assay = "SCT",
  slot = "data"
)

n_all_atac_peaks <- nrow(peak_counts)

peak_positive_cells <- Matrix::rowSums(peak_counts > 0)
gene_positive_cells <- Matrix::rowSums(gene_data[goi, , drop = FALSE] > 0)

smad3_positive_cells <- as.numeric(gene_positive_cells[1])

if (smad3_positive_cells <= min_cells) {
  stop(
    "SMAD3 is expressed in only ",
    smad3_positive_cells,
    " cells; it would not have been eligible for LinkPeaks."
  )
}

# Signac uses strictly > min.cells, not >= min.cells.
peaks_keep <- peak_positive_cells > min_cells

n_min_cells_peaks <- sum(peaks_keep)

peak_ranges <- granges(pre_obj[["ATAC"]])
peak_ranges <- peak_ranges[peaks_keep]

annot <- Annotation(pre_obj[["ATAC"]])

if (is.null(annot) || length(annot) == 0L) {
  stop("ATAC gene annotation is unavailable.")
}

CollapseToLongestTranscript <- getFromNamespace(
  "CollapseToLongestTranscript",
  "Signac"
)

DistanceToTSS <- getFromNamespace(
  "DistanceToTSS",
  "Signac"
)

gene_coords <- CollapseToLongestTranscript(
  ranges = annot
)

gene_coords_use <- gene_coords[
  gene_coords$gene_name == goi
]

if (length(gene_coords_use) != 1L) {
  stop(
    "Expected one collapsed coordinate entry for SMAD3; found ",
    length(gene_coords_use), "."
  )
}

peak_distance_matrix <- DistanceToTSS(
  peaks = peak_ranges,
  genes = gene_coords_use,
  distance = distance_bp
)

if (ncol(peak_distance_matrix) != 1L) {
  stop(
    "Unexpected SMAD3 distance-matrix dimensions: ",
    paste(dim(peak_distance_matrix), collapse = " x ")
  )
}

cis_keep <- as.logical(peak_distance_matrix[, 1])

n_cis_eligible <- sum(cis_keep)

if (n_cis_eligible < 1L) {
  stop("No SMAD3 cis peaks were reconstructed.")
}

cis_peak_names <- rownames(peak_distance_matrix)[cis_keep]

cat("All ATAC peaks                         :", n_all_atac_peaks, "\n")
cat("Peaks present in >10 cells             :", n_min_cells_peaks, "\n")
cat("SMAD3-positive cells                   :", smad3_positive_cells, "\n")
cat("Eligible SMAD3 cis peaks within 500 kb :", n_cis_eligible, "\n")

reported_peak_in_family <- any(
  grepl(reported_peak_pattern, cis_peak_names)
)

cat(
  "Reported peak present in reconstructed family:",
  reported_peak_in_family,
  "\n"
)

if (!reported_peak_in_family) {
  stop(
    "The reported peak was not found among reconstructed eligible cis peaks."
  )
}

# -------------------------------------------------------------------------
# 3. Conservative multiplicity corrections
# -------------------------------------------------------------------------

cat("\n============================================================\n")
cat("MULTIPLE-TESTING CORRECTION\n")
cat("============================================================\n\n")

bonf_cis <- min(
  raw_p * n_cis_eligible,
  1
)

# Deliberately extreme sensitivity analysis:
# pretend every peak in the entire ATAC assay belonged to the test family.
bonf_all_atac <- min(
  raw_p * n_all_atac_peaks,
  1
)

max_tests_for_bonf_005 <- floor(
  0.05 / raw_p
)

cat(
  "Raw P                               :",
  format(raw_p, scientific = TRUE, digits = 12),
  "\n"
)

cat(
  "Bonferroni P across SMAD3 cis peaks :",
  format(bonf_cis, scientific = TRUE, digits = 12),
  "\n"
)

cat(
  "Pass cis Bonferroni < 0.05           :",
  bonf_cis < 0.05,
  "\n"
)

cat(
  "Bonferroni P across ALL ATAC peaks   :",
  format(bonf_all_atac, scientific = TRUE, digits = 12),
  "\n"
)

cat(
  "Pass all-ATAC Bonferroni < 0.05      :",
  bonf_all_atac < 0.05,
  "\n"
)

cat(
  "Maximum number of tests compatible with\n",
  "Bonferroni-adjusted P < 0.05         :",
  max_tests_for_bonf_005,
  "\n"
)

# -------------------------------------------------------------------------
# 4. Export audit tables
# -------------------------------------------------------------------------

summary_dt <- data.table(
  gene = goi,
  reported_peak = reported_link$peak[1],
  link_score = link_score,
  link_zscore = link_z,
  raw_pvalue = raw_p,
  n_historical_links_retained = nrow(smad3_links),
  n_all_ATAC_peaks = n_all_atac_peaks,
  n_peaks_positive_gt_10_cells = n_min_cells_peaks,
  SMAD3_positive_cells = smad3_positive_cells,
  n_eligible_SMAD3_cis_peaks_500kb = n_cis_eligible,
  reported_peak_in_reconstructed_family = reported_peak_in_family,
  bonferroni_p_cis_family = bonf_cis,
  bonferroni_pass_cis_0_05 = bonf_cis < 0.05,
  bonferroni_p_all_ATAC_peaks = bonf_all_atac,
  bonferroni_pass_all_ATAC_0_05 = bonf_all_atac < 0.05,
  max_tests_for_bonferroni_p_lt_0_05 = max_tests_for_bonf_005
)

fwrite(
  summary_dt,
  file.path(
    outdir,
    "SMAD3_LinkPeaks_multiple_testing_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  smad3_links[
    order(pvalue),
    .(
      gene,
      peak,
      score,
      zscore,
      pvalue
    )
  ],
  file.path(
    outdir,
    "SMAD3_historical_LinkPeaks_retained_links.tsv"
  ),
  sep = "\t"
)

fwrite(
  data.table(
    peak = cis_peak_names
  ),
  file.path(
    outdir,
    "SMAD3_reconstructed_eligible_cis_peaks.tsv"
  ),
  sep = "\t"
)

cat("\n============================================================\n")
cat("FILES WRITTEN\n")
cat("============================================================\n\n")

cat(outdir, "\n")
cat("\nNo canonical analysis file was modified.\n")
