# ============================================================================
# Production Figure 10A-C
#
# Figure 10A:
#   Discovery-stage integration using:
#     - Yellow / Brown / Blue hdWGCNA modules
#     - canonical revised FAP pseudobulk DEGs
#       (BH-FDR < 0.05; FAP1-FAP3)
#     - UKB SMR significant after cohort-wide BH + HEIDI
#
# Figure 10B:
#   SMAD3 SMR/HEIDI across UKB tissues.
#   Raw SMR P values are displayed; highlighting requires:
#     cohort-wide BH FDR < 0.05 AND HEIDI P > 0.01.
#
# Figure 10C:
#   Targeted FinnGen replication of the exact hypotheses carried forward
#   from corrected UKB discovery.
#   BH correction is applied across the carried-forward hypotheses and
#   HEIDI P > 0.01 is additionally required.
#
# All inferential inputs are reconstructed directly from canonical
# processed-result resources. No reviewer-audit result is required.
# ============================================================================

rm(list = ls(all.names = TRUE))
gc()

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggVennDiagram)
  library(patchwork)
})

# ============================================================================
# Repository root
# ============================================================================

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (length(file_arg) == 1) {
  script_path <- sub("^--file=", "", file_arg)
  script_dir <- dirname(normalizePath(script_path))
  repo_root <- normalizePath(file.path(script_dir, ".."))
} else {
  repo_root <- normalizePath(".")
}

# ============================================================================
# Canonical inputs
# ============================================================================

module_file <- file.path(
  repo_root,
  "processed_results",
  "05_hdWGCNA",
  "module_assignment_table.csv"
)

deg_file <- file.path(
  repo_root,
  "processed_results",
  "02_differential_expression",
  "bulk_de_sig_faps.csv"
)

smr_file <- file.path(
  repo_root,
  "processed_results",
  "06_SMR_HEIDI",
  "SMR_all_studies_all_tissues.csv"
)

input_files <- c(
  module_file,
  deg_file,
  smr_file
)

missing_files <- input_files[
  !file.exists(input_files)
]

if (length(missing_files)) {
  stop(
    "Missing required canonical input file(s):\n",
    paste(missing_files, collapse = "\n")
  )
}

# ============================================================================
# Production output directories
# ============================================================================

figdir <- file.path(
  repo_root,
  "outputs",
  "SMR_HEIDI",
  "Figure10"
)

srcdir <- file.path(
  repo_root,
  "processed_results",
  "06_SMR_HEIDI",
  "figure10"
)

dir.create(
  figdir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  srcdir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================================
# Load and construct cohort-specific SMR resources
# ============================================================================

modules <- fread(module_file)
deg     <- fread(deg_file)
smr_all <- fread(smr_file)

required_module_cols <- c(
  "module",
  "gene_name"
)

required_deg_cols <- c(
  "gene",
  "celltype",
  "p_val_adj",
  "avg_log2FC"
)

required_smr_cols <- c(
  "study",
  "tissue",
  "probeID",
  "Gene",
  "topSNP",
  "b_SMR",
  "se_SMR",
  "p_SMR",
  "p_HEIDI",
  "nsnp_HEIDI"
)

check_cols <- function(
    x,
    req,
    label
) {

  miss <- setdiff(
    req,
    names(x)
  )

  if (length(miss)) {
    stop(
      label,
      " is missing column(s): ",
      paste(miss, collapse = ", ")
    )
  }
}

check_cols(
  modules,
  required_module_cols,
  "Module table"
)

check_cols(
  deg,
  required_deg_cols,
  "Canonical revised DEG table"
)

check_cols(
  smr_all,
  required_smr_cols,
  "Canonical combined SMR table"
)

# Normalize tissue labels to the historical underscore convention used
# throughout Figure 10 source data and plotting.
normalize_tissue <- function(x) {
  gsub(
    " ",
    "_",
    x,
    fixed = TRUE
  )
}

prepare_cohort <- function(
    x,
    cohort_name
) {

  z <- copy(
    x[
      study == cohort_name
    ]
  )

  if (nrow(z) == 0L) {
    stop(
      "No SMR rows found for cohort: ",
      cohort_name
    )
  }

  z[
    ,
    tissue :=
      normalize_tissue(tissue)
  ]

  z[
    ,
    valid_p :=
      is.finite(p_SMR) &
      p_SMR >= 0 &
      p_SMR <= 1
  ]

  z[
    ,
    FDR_cohortwide :=
      NA_real_
  ]

  z[
    valid_p == TRUE,
    FDR_cohortwide :=
      p.adjust(
        p_SMR,
        method = "BH"
      )
  ]

  z
}

ukb <- prepare_cohort(
  smr_all,
  "UKB"
)

fin <- prepare_cohort(
  smr_all,
  "FinnGen"
)

# ============================================================================
# Figure 10A: revised discovery-stage integration
# ============================================================================

yellow_genes <- sort(unique(
  na.omit(modules[tolower(module) == "yellow", gene_name])
))

brown_genes <- sort(unique(
  na.omit(modules[tolower(module) == "brown", gene_name])
))

blue_genes <- sort(unique(
  na.omit(modules[tolower(module) == "blue", gene_name])
))

deg_genes <- sort(unique(na.omit(deg$gene)))

# Revision validation
if (length(deg_genes) != 194L) {
  warning(
    "Expected 194 unique revised FAP genes, but found ",
    length(deg_genes)
  )
}

# Corrected UKB discovery rows:
# cohort-wide BH across UKB gene-tissue tests + HEIDI criterion
ukb_discovery_rows <- ukb[
  !is.na(FDR_cohortwide) &
  FDR_cohortwide < 0.05 &
  !is.na(p_HEIDI) &
  p_HEIDI > 0.01
]

smr_ukb_bh_genes <- sort(unique(
  na.omit(ukb_discovery_rows$Gene)
))

module_sets <- list(
  Yellow = yellow_genes,
  Brown = brown_genes,
  Blue = blue_genes
)

candidate_list <- lapply(
  names(module_sets),
  function(mod_name) {

    genes <- Reduce(
      intersect,
      list(
        module_sets[[mod_name]],
        deg_genes,
        smr_ukb_bh_genes
      )
    )

    data.table(
      Gene = sort(genes),
      Module = mod_name
    )
  }
)

discovery_candidates <- rbindlist(
  candidate_list,
  use.names = TRUE
)

# Add exact corrected UKB evidence for each candidate.
# If >1 passing SMR row exists for a candidate, retain all rows.
discovery_candidate_details <- merge(
  discovery_candidates,
  ukb_discovery_rows[
    ,
    .(
      Gene,
      tissue,
      probeID,
      topSNP,
      b_SMR,
      se_SMR,
      p_SMR,
      FDR_cohortwide,
      p_HEIDI,
      nsnp_HEIDI
    )
  ],
  by = "Gene",
  all.x = TRUE,
  allow.cartesian = TRUE
)

setorder(
  discovery_candidate_details,
  Module,
  Gene,
  FDR_cohortwide,
  p_SMR
)

cat("\n============================================================\n")
cat("FIGURE 10A — UKB DISCOVERY\n")
cat("============================================================\n\n")

cat("Yellow module :", length(yellow_genes), "\n")
cat("Brown module  :", length(brown_genes), "\n")
cat("Blue module   :", length(blue_genes), "\n")
cat("FAP DEGs      :", length(deg_genes), "unique genes\n")
cat("UKB SMR BH    :", length(smr_ukb_bh_genes), "unique genes\n")

cat("\nModule-specific convergent candidates:\n")
print(discovery_candidates)

expected_candidates <- c("PLEKHA6", "SMAD3")
observed_candidates <- sort(unique(discovery_candidates$Gene))

if (!setequal(observed_candidates, expected_candidates)) {
  warning(
    "Candidate set differs from the validated audit. Observed: ",
    paste(observed_candidates, collapse = ", ")
  )
}

# --------------------------------------------------------------------------
# UpSet
# --------------------------------------------------------------------------

sets <- list(
  "Yellow module" = yellow_genes,
  "Brown module" = brown_genes,
  "Blue module" = blue_genes,
  "FAP1 to FAP3 DEGs" = deg_genes,
  "UKB SMR BH and HEIDI" = smr_ukb_bh_genes
)

all_genes <- sort(unique(unlist(sets)))

fig10a_source <- data.table(
  gene = all_genes,
  Yellow_Module = as.integer(all_genes %in% yellow_genes),
  Brown_Module = as.integer(all_genes %in% brown_genes),
  Blue_Module = as.integer(all_genes %in% blue_genes),
  FAP_DEGs_BH = as.integer(all_genes %in% deg_genes),
  UKB_SMR_BH = as.integer(all_genes %in% smr_ukb_bh_genes)
)

# Number of actually observed non-empty membership patterns.
# Use this to avoid plotting unnecessary zero-intersection columns.
pattern_counts <- fig10a_source[
  ,
  .N,
  by = .(
    Yellow_Module,
    Brown_Module,
    Blue_Module,
    FAP_DEGs_BH,
    UKB_SMR_BH
  )
][N > 0]

n_intersections_to_plot <- max(
  1L,
  min(26L, nrow(pattern_counts))
)

cat(
  "\nObserved non-empty UpSet membership patterns:",
  nrow(pattern_counts),
  "\n"
)

venn <- Venn(sets)

p10a <- plot_upset(
  venn,
  nintersects = n_intersections_to_plot,
  order.intersect.by = "size",
  relative_height = 2,
  relative_width = 0.3
)

ggsave(
  file.path(figdir, "Figure10A.png"),
  p10a,
  width = 8,
  height = 5,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(figdir, "Figure10A.pdf"),
  p10a,
  width = 8,
  height = 5,
  units = "in",
  bg = "white"
)

fwrite(
  fig10a_source,
  file.path(srcdir, "Figure10A_Upset_source_data.csv")
)

fwrite(
  pattern_counts,
  file.path(srcdir, "Figure10A_nonzero_intersection_patterns.csv")
)

fwrite(
  discovery_candidates,
  file.path(srcdir, "Figure10A_discovery_candidates.csv")
)

fwrite(
  discovery_candidate_details,
  file.path(srcdir, "Figure10A_discovery_candidate_details.csv")
)

fwrite(
  ukb_discovery_rows,
  file.path(srcdir, "Figure10A_UKB_corrected_discovery_rows.tsv"),
  sep = "\t"
)

# ============================================================================
# Canonical targeted FinnGen replication
# ============================================================================
#
# Carry forward the exact UKB discovery hypotheses represented by the two
# convergent candidates and match them in FinnGen by:
#   Gene + tissue + probeID.
#
# BH is applied across these two targeted replication hypotheses.
# ============================================================================

if (
  nrow(discovery_candidate_details) != 2L ||
  !setequal(
    discovery_candidate_details$Gene,
    c(
      "SMAD3",
      "PLEKHA6"
    )
  )
) {
  stop(
    "Expected exactly two convergent UKB discovery candidates: ",
    "SMAD3 and PLEKHA6."
  )
}

rep_primary <- merge(
  discovery_candidate_details[
    ,
    .(
      Gene,
      Module,
      tissue,
      probeID,
      topSNP_UKB = topSNP,
      p_SMR_UKB = p_SMR,
      FDR_UKB = FDR_cohortwide,
      p_HEIDI_UKB = p_HEIDI
    )
  ],
  fin[
    ,
    .(
      Gene,
      tissue,
      probeID,
      topSNP_FinnGen = topSNP,
      b_SMR_FinnGen = b_SMR,
      se_SMR_FinnGen = se_SMR,
      p_SMR_FinnGen = p_SMR,
      p_HEIDI_FinnGen = p_HEIDI,
      nsnp_HEIDI_FinnGen = nsnp_HEIDI
    )
  ],
  by = c(
    "Gene",
    "tissue",
    "probeID"
  ),
  all.x = TRUE,
  allow.cartesian = FALSE
)

if (nrow(rep_primary) != 2L) {
  stop(
    "Expected exactly two exact FinnGen replication rows; found ",
    nrow(rep_primary)
  )
}

if (
  anyNA(rep_primary$p_SMR_FinnGen) ||
  anyNA(rep_primary$p_HEIDI_FinnGen)
) {
  stop(
    "Missing FinnGen evidence for a carried-forward hypothesis."
  )
}

rep_primary[
  ,
  replication_BH :=
    p.adjust(
      p_SMR_FinnGen,
      method = "BH"
    )
]

rep_primary[
  ,
  replication_pass :=
    replication_BH < 0.05 &
    p_HEIDI_FinnGen > 0.01
]

rep_primary[
  ,
  framework :=
    "Targeted FinnGen replication of corrected UKB discovery"
]

rep_primary[
  ,
  revised_analysis_role :=
    paste0(
      "Primary: UKB cohort-wide BH-FDR < 0.05 + HEIDI; ",
      "FinnGen targeted BH replication + HEIDI"
    )
]

setorder(
  rep_primary,
  Gene
)

# Exact safeguards established in the independent C7 validation.
smad3_rep_check <- rep_primary[
  Gene == "SMAD3"
]

plekha6_rep_check <- rep_primary[
  Gene == "PLEKHA6"
]

if (
  nrow(smad3_rep_check) != 1L ||
  abs(
    smad3_rep_check$p_SMR_FinnGen -
      0.001071558
  ) > 1e-12 ||
  abs(
    smad3_rep_check$replication_BH -
      0.002143116
  ) > 1e-12 ||
  !isTRUE(
    smad3_rep_check$replication_pass
  )
) {
  stop(
    "Canonical SMAD3 FinnGen replication does not match the validated result."
  )
}

if (
  nrow(plekha6_rep_check) != 1L ||
  abs(
    plekha6_rep_check$replication_BH -
      0.373608
  ) > 1e-12 ||
  isTRUE(
    plekha6_rep_check$replication_pass
  )
) {
  stop(
    "Canonical PLEKHA6 FinnGen replication does not match the validated result."
  )
}

cat("\n============================================================\n")
cat("TARGETED FINNGEN REPLICATION HYPOTHESES\n")
cat("============================================================\n\n")

print(
  rep_primary[
    ,
    .(
      Gene,
      Module,
      tissue,
      probeID,
      p_SMR_UKB,
      FDR_UKB,
      p_HEIDI_UKB,
      p_SMR_FinnGen,
      replication_BH,
      p_HEIDI_FinnGen,
      replication_pass
    )
  ]
)

fwrite(
  rep_primary,
  file.path(
    srcdir,
    "Figure10_targeted_FinnGen_replication_primary.csv"
  )
)

# Also retain a general canonical replication resource outside the
# Figure-10-specific directory.
fwrite(
  rep_primary,
  file.path(
    repo_root,
    "processed_results",
    "06_SMR_HEIDI",
    "targeted_FinnGen_replication.csv"
  )
)

# ============================================================================
# Figure 10B-C helper
# ============================================================================

pretty_tissue_name <- function(x) {
  out <- x
  out[x == "Adipose_Subcutaneous"] <- "Adipose subcutaneous"
  out[x == "Cells_Cultured_fibroblasts"] <- "Cultured fibroblasts"
  out[x == "Muscle_Skeletal"] <- "Skeletal muscle"
  out[x == "Adipose_Visceral_Omentum"] <- "Adipose visceral omentum"
  out
}

prepare_context_rows <- function(
  smr,
  gene,
  target_tissue,
  target_probe
) {

  z <- copy(
    smr[
      Gene == gene
    ]
  )

  if (nrow(z) == 0L) {
    stop("No SMR rows found for ", gene)
  }

  z[
    ,
    target_hypothesis :=
      tissue == target_tissue &
      probeID == target_probe
  ]

  # For the exact carried hypothesis, force the exact probe to be shown.
  # For other tissues retain the smallest raw SMR P as contextual display.
  setorderv(
    z,
    c("tissue", "target_hypothesis", "p_SMR"),
    c(1L, -1L, 1L),
    na.last = TRUE
  )

  z <- z[
    ,
    .SD[1],
    by = tissue
  ]

  z
}

make_smr_heidi_plot <- function(
  dat,
  study_label,
  subtitle_text,
  highlight_fill
) {

  d <- copy(dat)

  d[
    ,
    log_pSMR := -log10(
      pmax(p_SMR, .Machine$double.xmin)
    )
  ]

  d[
    ,
    log_pHEIDI := -log10(
      pmax(p_HEIDI, .Machine$double.xmin)
    )
  ]

  setorder(d, -log_pSMR)

  d[
    ,
    tissue_display := pretty_tissue_name(tissue)
  ]

  d[
    ,
    tissue_snp_chr := paste(
      tissue_display,
      topSNP,
      sep = "\n"
    )
  ]

  d[
    ,
    tissue_snp := factor(
      tissue_snp_chr,
      levels = tissue_snp_chr
    )
  ]

  d[
    ,
    fill_group := ifelse(
      plot_highlight,
      "highlight",
      "other"
    )
  ]

  fill_values <- c(
    highlight = highlight_fill,
    other = "grey80"
  )

  ymax_smr <- max(
    d$log_pSMR,
    na.rm = TRUE
  ) + 1.2

  thr_heidi <- -log10(0.01)

  ymax_heidi <- max(
    c(
      d$log_pHEIDI,
      thr_heidi
    ),
    na.rm = TRUE
  ) + 0.5

  # ------------------------------------------------------------------------
  # SMR
  # ------------------------------------------------------------------------

  p_smr <- ggplot(
    d,
    aes(
      x = tissue_snp,
      y = log_pSMR,
      fill = fill_group
    )
  ) +
    geom_col(width = 0.8) +
    geom_label(
      aes(
        label = round(log_pSMR, 2)
      ),
      fill = "white",
      color = "black",
      linewidth = 0.25,
      size = 3,
      vjust = -0.2
    ) +
    geom_label(
      data = d[
        !is.na(plot_annotation) &
        plot_annotation != ""
      ],
      aes(
        x = tissue_snp,
        y = log_pSMR,
        label = plot_annotation
      ),
      inherit.aes = FALSE,
      vjust = -1.5,
      size = 2.8,
      linewidth = 0.25,
      fill = "white"
    ) +
    scale_fill_manual(
      values = fill_values,
      guide = "none"
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.05))
    ) +
    coord_cartesian(
      ylim = c(0, ymax_smr),
      clip = "off"
    ) +
    labs(
      title = paste0(
        "SMR results (",
        study_label,
        ")"
      ),
      subtitle = subtitle_text,
      x = NULL,
      y = expression(-log[10](P[SMR]))
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(
        face = "italic",
        hjust = 0
      ),
      plot.subtitle = element_text(
        size = 8
      ),
      axis.text.x.bottom = element_text(
        angle = 30,
        hjust = 1,
        vjust = 1
      ),
      axis.text.y = element_text(size = 9),
      plot.margin = margin(
        t = 5,
        r = 20,
        b = 0,
        l = 5
      )
    )

  # ------------------------------------------------------------------------
  # HEIDI
  # ------------------------------------------------------------------------

  thr_lab <- data.table(
    tissue_snp = tail(
      levels(d$tissue_snp),
      1
    ),
    y = thr_heidi
  )

  p_heidi <- ggplot(
    d,
    aes(
      x = tissue_snp,
      y = log_pHEIDI,
      fill = fill_group
    )
  ) +
    geom_col(width = 0.8) +
    geom_label(
      aes(
        label = round(log_pHEIDI, 2)
      ),
      fill = "white",
      color = "black",
      linewidth = 0.25,
      size = 3,
      vjust = 1.2
    ) +
    geom_hline(
      yintercept = thr_heidi,
      linetype = "dashed",
      linewidth = 0.6
    ) +
    geom_label(
      data = thr_lab,
      aes(
        x = tissue_snp,
        y = y,
        label = "HEIDI exclusion threshold: p = 0.01"
      ),
      inherit.aes = FALSE,
      hjust = 1.05,
      vjust = 1.3,
      size = 2.7,
      linewidth = 0.25,
      fill = "white"
    ) +
    scale_fill_manual(
      values = fill_values,
      guide = "none"
    ) +
    scale_x_discrete(
      position = "top"
    ) +
    scale_y_reverse(
      expand = expansion(mult = c(0, 0))
    ) +
    coord_cartesian(
      ylim = c(ymax_heidi, 0),
      clip = "off"
    ) +
    labs(
      x = NULL,
      y = expression(-log[10](P[HEIDI]))
    ) +
    theme_classic() +
    theme(
      axis.text.x.top = element_blank(),
      axis.text.y = element_text(size = 9),
      plot.margin = margin(
        t = 0,
        r = 20,
        b = 5,
        l = 5
      )
    )

  p_smr / p_heidi +
    plot_layout(
      heights = c(2, 1)
    )
}

# ============================================================================
# Exact SMAD3 carried hypothesis
# ============================================================================

smad3_rep <- rep_primary[
  Gene == "SMAD3"
]

if (nrow(smad3_rep) != 1L) {
  stop(
    "Expected exactly one SMAD3 carried-forward hypothesis; found ",
    nrow(smad3_rep)
  )
}

smad3_tissue <- smad3_rep$tissue[1]
smad3_probe  <- smad3_rep$probeID[1]

# ============================================================================
# Figure 10B — UKB discovery
# ============================================================================

ukb_smad3 <- prepare_context_rows(
  ukb,
  "SMAD3",
  smad3_tissue,
  smad3_probe
)

if (
  !any(
    ukb_smad3$tissue == smad3_tissue &
    ukb_smad3$probeID == smad3_probe
  )
) {
  stop(
    "Exact UKB SMAD3 discovery hypothesis not present in Figure 10B rows."
  )
}

ukb_smad3[
  ,
  corrected_discovery :=
    !is.na(FDR_cohortwide) &
    FDR_cohortwide < 0.05 &
    !is.na(p_HEIDI) &
    p_HEIDI > 0.01
]

ukb_smad3[
  ,
  plot_highlight := corrected_discovery
]

ukb_smad3[
  ,
  plot_annotation := ifelse(
    corrected_discovery,
    paste0(
      "BH FDR = ",
      formatC(
        FDR_cohortwide,
        format = "f",
        digits = 4
      )
    ),
    ""
  )
]

# Validate exact audit values
ukb_target <- ukb_smad3[
  tissue == smad3_tissue &
  probeID == smad3_probe
]

if (
  nrow(ukb_target) != 1L ||
  abs(
    ukb_target$p_SMR -
    smad3_rep$p_SMR_UKB[1]
  ) > 1e-12
) {
  warning(
    "UKB plotted SMAD3 target does not exactly match the replication audit row."
  )
}

p10b <- make_smr_heidi_plot(
  ukb_smad3,
  "UKB",
  "Highlighted: cohort-wide BH FDR < 0.05 and HEIDI p > 0.01",
  "#6EC7D4"
)

ggsave(
  file.path(figdir, "Figure10B_UKB_SMAD3.png"),
  p10b,
  width = 6,
  height = 5,
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(figdir, "Figure10B_UKB_SMAD3.pdf"),
  p10b,
  width = 6,
  height = 5,
  bg = "white"
)

fwrite(
  ukb_smad3,
  file.path(srcdir, "Figure10B_UKB_SMAD3_source_data.csv")
)

# ============================================================================
# Figure 10C — FinnGen targeted replication
# ============================================================================

fin_smad3 <- prepare_context_rows(
  fin,
  "SMAD3",
  smad3_tissue,
  smad3_probe
)

if (
  !any(
    fin_smad3$tissue == smad3_tissue &
    fin_smad3$probeID == smad3_probe
  )
) {
  stop(
    "Exact FinnGen SMAD3 replication hypothesis not present in Figure 10C rows."
  )
}

fin_smad3[
  ,
  targeted_replication :=
    tissue == smad3_tissue &
    probeID == smad3_probe
]

fin_smad3[
  ,
  replication_BH_plot := NA_real_
]

fin_smad3[
  targeted_replication == TRUE,
  replication_BH_plot :=
    smad3_rep$replication_BH[1]
]

fin_smad3[
  ,
  replication_pass_plot := FALSE
]

fin_smad3[
  targeted_replication == TRUE,
  replication_pass_plot :=
    as.logical(
      smad3_rep$replication_pass[1]
    )
]

fin_smad3[
  ,
  plot_highlight := replication_pass_plot
]

fin_smad3[
  ,
  plot_annotation := ifelse(
    targeted_replication,
    paste0(
      "Targeted replication\nBH = ",
      formatC(
        replication_BH_plot,
        format = "f",
        digits = 5
      )
    ),
    ""
  )
]

# Validate exact FinnGen P against replication audit
fin_target <- fin_smad3[
  targeted_replication == TRUE
]

if (
  nrow(fin_target) != 1L ||
  abs(
    fin_target$p_SMR -
    smad3_rep$p_SMR_FinnGen[1]
  ) > 1e-12
) {
  warning(
    "FinnGen plotted SMAD3 target does not exactly match the replication audit row."
  )
}

p10c <- make_smr_heidi_plot(
  fin_smad3,
  "FinnGen",
  "Highlighted: UKB-carried hypothesis replicated after BH correction and HEIDI",
  "#F5D4A4"
)

ggsave(
  file.path(figdir, "Figure10C_FinnGen_SMAD3.png"),
  p10c,
  width = 6,
  height = 5,
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(figdir, "Figure10C_FinnGen_SMAD3.pdf"),
  p10c,
  width = 6,
  height = 5,
  bg = "white"
)

fwrite(
  fin_smad3,
  file.path(srcdir, "Figure10C_FinnGen_SMAD3_source_data.csv")
)

# ============================================================================
# Combined B-C preview
# ============================================================================

p10bc <- p10b | p10c

ggsave(
  file.path(figdir, "Figure10BC_preview.png"),
  p10bc,
  width = 12,
  height = 5,
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(figdir, "Figure10BC_preview.pdf"),
  p10bc,
  width = 12,
  height = 5,
  bg = "white"
)

# ============================================================================
# Console validation summary
# ============================================================================

cat("\n============================================================\n")
cat("FIGURE 10B — UKB SMAD3\n")
cat("============================================================\n\n")

print(
  ukb_smad3[
    ,
    .(
      tissue,
      probeID,
      topSNP,
      p_SMR,
      FDR_cohortwide,
      p_HEIDI,
      corrected_discovery
    )
  ]
)

cat("\n============================================================\n")
cat("FIGURE 10C — FINNGEN SMAD3\n")
cat("============================================================\n\n")

print(
  fin_smad3[
    ,
    .(
      tissue,
      probeID,
      topSNP,
      p_SMR,
      p_HEIDI,
      targeted_replication,
      replication_BH_plot,
      replication_pass_plot
    )
  ]
)

cat("\n============================================================\n")
cat("FILES WRITTEN\n")
cat("============================================================\n\n")

cat("Figures:\n", figdir, "\n")
cat("Source data:\n", srcdir, "\n")

cat("\nProduction Figure 10A-C panels and processed source resources were generated.\n")
