# Reviewer C7: multiple-testing audit for SMR results
# Applies Benjamini-Hochberg correction separately within each
# GWAS x GTEx tissue SMR analysis.
#
# Original .smr files are not modified.

input_dir <- "processed_results/06_SMR_HEIDI"
output_dir <- file.path(input_dir, "multiple_testing_audit")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(
  input_dir,
  pattern = "^(UKB|FinnGen)_.*\\.smr$",
  full.names = TRUE
)

if (length(files) != 8) {
  stop("Expected exactly 8 .smr files, but found ", length(files))
}

summary_list <- list()
smad3_list <- list()

for (f in sort(files)) {

  dat <- read.table(
    f,
    header = TRUE,
    sep = "",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = ""
  )

  required_cols <- c(
    "probeID", "Gene", "b_SMR", "se_SMR",
    "p_SMR", "p_HEIDI", "nsnp_HEIDI"
  )

  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0) {
    stop(
      "Missing columns in ", basename(f), ": ",
      paste(missing_cols, collapse = ", ")
    )
  }

  dat$p_SMR <- as.numeric(dat$p_SMR)
  dat$p_HEIDI <- as.numeric(dat$p_HEIDI)

  valid <- is.finite(dat$p_SMR) &
           dat$p_SMR >= 0 &
           dat$p_SMR <= 1

  # BH correction performed separately within each GWAS x tissue analysis
  dat$FDR_SMR_BH <- NA_real_
  dat$FDR_SMR_BH[valid] <- p.adjust(
    dat$p_SMR[valid],
    method = "BH"
  )

  # Rank among valid SMR P values
  dat$SMR_P_rank <- NA_integer_
  dat$SMR_P_rank[valid] <- rank(
    dat$p_SMR[valid],
    ties.method = "min"
  )

  dat$pass_nominal_SMR <- valid & dat$p_SMR < 0.05
  dat$pass_BH_SMR <- !is.na(dat$FDR_SMR_BH) &
                     dat$FDR_SMR_BH < 0.05

  # Mirrors the existing downstream shared-signal HEIDI criterion,
  # but replaces nominal P_SMR < 0.05 with BH-FDR < 0.05
  dat$pass_BH_SMR_HEIDI <- dat$pass_BH_SMR &
                           !is.na(dat$p_HEIDI) &
                           dat$p_HEIDI > 0.01

  filename <- basename(f)
  stem <- sub("\\.smr$", "", filename)

  cohort <- if (grepl("^UKB_", stem)) "UKB" else "FinnGen"
  tissue <- sub("^(UKB|FinnGen)_", "", stem)

  n_valid <- sum(valid)

  summary_list[[filename]] <- data.frame(
    file = filename,
    cohort = cohort,
    tissue = tissue,
    n_rows = nrow(dat),
    n_valid_p_SMR = n_valid,
    n_nominal_p_SMR_lt_0.05 = sum(dat$pass_nominal_SMR, na.rm = TRUE),
    n_BH_FDR_lt_0.05 = sum(dat$pass_BH_SMR, na.rm = TRUE),
    n_BH_FDR_lt_0.05_and_p_HEIDI_gt_0.01 =
      sum(dat$pass_BH_SMR_HEIDI, na.rm = TRUE),
    min_p_SMR = min(dat$p_SMR[valid], na.rm = TRUE),
    min_FDR_SMR_BH = min(dat$FDR_SMR_BH, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  smad3 <- dat[
    !is.na(dat$Gene) & toupper(dat$Gene) == "SMAD3",
    ,
    drop = FALSE
  ]

  if (nrow(smad3) == 0) {

    smad3_list[[filename]] <- data.frame(
      file = filename,
      cohort = cohort,
      tissue = tissue,
      SMAD3_found = FALSE,
      probeID = NA,
      topSNP = NA,
      b_SMR = NA,
      se_SMR = NA,
      p_SMR = NA,
      SMR_P_rank = NA,
      n_valid_p_SMR = n_valid,
      FDR_SMR_BH = NA,
      p_HEIDI = NA,
      nsnp_HEIDI = NA,
      pass_nominal_SMR = FALSE,
      pass_BH_SMR = FALSE,
      pass_BH_SMR_HEIDI = FALSE,
      stringsAsFactors = FALSE
    )

  } else {

    for (i in seq_len(nrow(smad3))) {

      x <- smad3[i, , drop = FALSE]

      smad3_list[[paste0(filename, "_", i)]] <- data.frame(
        file = filename,
        cohort = cohort,
        tissue = tissue,
        SMAD3_found = TRUE,
        probeID = x$probeID,
        topSNP = x$topSNP,
        b_SMR = x$b_SMR,
        se_SMR = x$se_SMR,
        p_SMR = x$p_SMR,
        SMR_P_rank = x$SMR_P_rank,
        n_valid_p_SMR = n_valid,
        FDR_SMR_BH = x$FDR_SMR_BH,
        p_HEIDI = x$p_HEIDI,
        nsnp_HEIDI = x$nsnp_HEIDI,
        pass_nominal_SMR = x$pass_nominal_SMR,
        pass_BH_SMR = x$pass_BH_SMR,
        pass_BH_SMR_HEIDI = x$pass_BH_SMR_HEIDI,
        stringsAsFactors = FALSE
      )
    }
  }

  # Save augmented audit table; original .smr remains untouched
  out_file <- file.path(
    output_dir,
    paste0(stem, "_BH_audit.tsv")
  )

  write.table(
    dat,
    out_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

summary_df <- do.call(rbind, summary_list)
smad3_df <- do.call(rbind, smad3_list)

row.names(summary_df) <- NULL
row.names(smad3_df) <- NULL

write.csv(
  summary_df,
  file.path(output_dir, "SMR_BH_file_summary.csv"),
  row.names = FALSE
)

write.csv(
  smad3_df,
  file.path(output_dir, "SMAD3_BH_summary.csv"),
  row.names = FALSE
)

cat("\n============================================================\n")
cat("SMR BH MULTIPLE-TESTING AUDIT\n")
cat("BH correction performed separately within each GWAS x tissue file\n")
cat("============================================================\n\n")

print(summary_df, row.names = FALSE)

cat("\n============================================================\n")
cat("SMAD3 RESULTS\n")
cat("============================================================\n\n")

print(
  smad3_df[
    ,
    c(
      "cohort", "tissue", "probeID", "topSNP",
      "p_SMR", "SMR_P_rank", "n_valid_p_SMR",
      "FDR_SMR_BH", "p_HEIDI",
      "pass_nominal_SMR", "pass_BH_SMR",
      "pass_BH_SMR_HEIDI"
    )
  ],
  row.names = FALSE
)

cat("\nAudit files written to:\n")
cat(output_dir, "\n")
