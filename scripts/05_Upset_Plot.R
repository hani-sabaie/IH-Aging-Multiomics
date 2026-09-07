# ===== Clean environment =====
rm(list = ls(all.names = TRUE))
gc()

# ===== Loading relevant libraries =====
library(ggVennDiagram)
library(ggplot2)
library(data.table)
library(dplyr)

# Locate repository root from this script
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (length(file_arg) == 1) {
  script_path <- sub("^--file=", "", file_arg)
  script_dir <- dirname(normalizePath(script_path))
  repo_root <- normalizePath(file.path(script_dir, ".."))
} else {
  repo_root <- normalizePath(".")
}

# ===== Helpers =====
outdir <- file.path(repo_root, "outputs", "Venn")
figdir <- outdir

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

save_gg <- function(p, filename, w=7, h=5, dpi=300){
  if (inherits(p, "ggplot") || inherits(p, "patchwork")){
    ggsave(file.path(figdir, filename), p, width=w, height=h, units="in", dpi=dpi, bg="white")
  } else {
    png(file.path(figdir, filename), width=w, height=h, units="in", res=dpi)
    print(p); dev.off()
  }
}

wcsv <- function(x, path) write.csv(x, path, row.names = TRUE)

# ============================================================================ #
# ===== UpsetPlot =====
modules <- read.csv(
  file.path(
    repo_root,
    "processed_results",
    "05_hdWGCNA",
    "module_assignment_table.csv"
  ),
  stringsAsFactors = FALSE
)
mods <- sort(unique(modules$module))
modules_genes <- setNames(vector("list", length(mods)), mods)

for (m in mods) {
  genes <- modules$gene_name[modules$module == m]
  genes <- unique(na.omit(trimws(genes)))
  
  write.csv(
    data.frame(gene_name = genes),
    file = file.path(outdir, paste0("module_", m, "_genes.csv")),
    row.names = FALSE
  )
  
  safe_m <- make.names(m)
  if (grepl("^[0-9]", safe_m)) safe_m <- paste0("X", safe_m)
  var_name <- paste0(safe_m, "_module_genes")
  
  assign(var_name, genes, envir = .GlobalEnv)
  
  modules_genes[[m]] <- genes
}


# ============================================================================
# Revised FAP-DEG / hdWGCNA / UKB-SMR integration
#
# Discovery-stage SMR significance:
#   - UKB only
#   - BH correction across all valid UKB gene-tissue SMR tests
#   - FDR < 0.05
#   - HEIDI P > 0.01
#
# FinnGen is not included in the discovery UpSet. It is used separately for
# targeted replication of discovery candidates.
# ============================================================================

# ---------------------------------------------------------------------------
# FAP DEGs
# ---------------------------------------------------------------------------
bulk_de_sig_faps <- data.table::fread(
  file.path(
    repo_root,
    "processed_results",
    "02_differential_expression",
    "bulk_de_sig_faps.csv"
  )
)

required_deg_cols <- c(
  "gene",
  "celltype",
  "p_val_adj",
  "avg_log2FC"
)

missing_deg_cols <- setdiff(
  required_deg_cols,
  names(bulk_de_sig_faps)
)

if (length(missing_deg_cols) > 0L) {
  stop(
    "FAP DEG table missing required columns: ",
    paste(missing_deg_cols, collapse = ", ")
  )
}

DEGs_FAPs_genes <- sort(
  unique(
    na.omit(
      bulk_de_sig_faps$gene
    )
  )
)

if (length(DEGs_FAPs_genes) != 194L) {
  stop(
    "Expected 194 unique revised FAP DEG genes; found ",
    length(DEGs_FAPs_genes)
  )
}

# ---------------------------------------------------------------------------
# Canonical combined SMR results
# ---------------------------------------------------------------------------
smr_all <- data.table::fread(
  file.path(
    repo_root,
    "processed_results",
    "06_SMR_HEIDI",
    "SMR_all_studies_all_tissues.csv"
  )
)

required_smr_cols <- c(
  "study",
  "tissue",
  "probeID",
  "Gene",
  "topSNP",
  "p_SMR",
  "p_HEIDI"
)

missing_smr_cols <- setdiff(
  required_smr_cols,
  names(smr_all)
)

if (length(missing_smr_cols) > 0L) {
  stop(
    "Combined SMR table missing required columns: ",
    paste(missing_smr_cols, collapse = ", ")
  )
}

ukb <- data.table::copy(
  smr_all[
    study == "UKB"
  ]
)

ukb[
  ,
  valid_p :=
    is.finite(p_SMR) &
    p_SMR >= 0 &
    p_SMR <= 1
]

ukb[
  ,
  FDR_cohortwide :=
    NA_real_
]

ukb[
  valid_p == TRUE,
  FDR_cohortwide :=
    p.adjust(
      p_SMR,
      method = "BH"
    )
]

ukb[
  ,
  pass_cohortwide :=
    !is.na(FDR_cohortwide) &
    FDR_cohortwide < 0.05 &
    !is.na(p_HEIDI) &
    p_HEIDI > 0.01
]

ukb_discovery_rows <- ukb[
  pass_cohortwide == TRUE
]

smr_ukb_bh_genes <- sort(
  unique(
    na.omit(
      ukb_discovery_rows$Gene
    )
  )
)

if (nrow(ukb_discovery_rows) != 138L) {
  stop(
    "Expected 138 corrected UKB discovery gene-tissue rows; found ",
    nrow(ukb_discovery_rows)
  )
}

if (length(smr_ukb_bh_genes) != 84L) {
  stop(
    "Expected 84 unique corrected UKB discovery genes; found ",
    length(smr_ukb_bh_genes)
  )
}

# ---------------------------------------------------------------------------
# Canonical corrected UKB discovery resources
# ---------------------------------------------------------------------------
smr_processed_dir <- file.path(
  repo_root,
  "processed_results",
  "06_SMR_HEIDI"
)

data.table::fwrite(
  ukb_discovery_rows[
    order(
      FDR_cohortwide,
      p_SMR,
      Gene,
      tissue
    )
  ],
  file.path(
    smr_processed_dir,
    "UKB_discovery_cohortwide_BH.csv"
  )
)

data.table::fwrite(
  data.table::data.table(
    Gene = smr_ukb_bh_genes
  ),
  file.path(
    smr_processed_dir,
    "UKB_discovery_cohortwide_BH_genes.csv"
  )
)

# ---------------------------------------------------------------------------
# hdWGCNA module sets
# ---------------------------------------------------------------------------
get_module_genes <- function(module_name) {

  sort(
    unique(
      na.omit(
        modules$gene_name[
          tolower(modules$module) ==
            tolower(module_name)
        ]
      )
    )
  )
}

yellow_module_genes <- get_module_genes("Yellow")
brown_module_genes  <- get_module_genes("Brown")
blue_module_genes   <- get_module_genes("Blue")

# ---------------------------------------------------------------------------
# Revised discovery-stage integration
# ---------------------------------------------------------------------------
module_sets <- list(
  Yellow = yellow_module_genes,
  Brown = brown_module_genes,
  Blue = blue_module_genes
)

integration_candidates <- data.table::rbindlist(
  lapply(
    names(module_sets),
    function(mod_name) {

      hits <- Reduce(
        intersect,
        list(
          module_sets[[mod_name]],
          DEGs_FAPs_genes,
          smr_ukb_bh_genes
        )
      )

      if (length(hits) == 0L) {
        return(
          data.table::data.table(
            Gene = character(),
            Module = character()
          )
        )
      }

      data.table::data.table(
        Gene = sort(hits),
        Module = mod_name
      )
    }
  ),
  fill = TRUE
)

# Expected study-specific convergence pattern.
yellow_hits <- integration_candidates[
  Module == "Yellow",
  Gene
]

brown_hits <- integration_candidates[
  Module == "Brown",
  Gene
]

blue_hits <- integration_candidates[
  Module == "Blue",
  Gene
]

if (length(yellow_hits) != 0L) {
  stop(
    "Unexpected Yellow-module convergent candidate(s): ",
    paste(yellow_hits, collapse = ", ")
  )
}

if (!identical(sort(brown_hits), "SMAD3")) {
  stop(
    "Expected Brown-module candidate SMAD3; found: ",
    paste(brown_hits, collapse = ", ")
  )
}

if (!identical(sort(blue_hits), "PLEKHA6")) {
  stop(
    "Expected Blue-module candidate PLEKHA6; found: ",
    paste(blue_hits, collapse = ", ")
  )
}

# Attach corrected UKB SMR evidence to convergent candidates.
candidate_details <- merge(
  integration_candidates,
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

candidate_details <- candidate_details[
  order(
    Module,
    Gene,
    FDR_cohortwide
  )
]

data.table::fwrite(
  integration_candidates,
  file.path(
    smr_processed_dir,
    "integration_FAP_hdWGCNA_UKB_candidates.csv"
  )
)

data.table::fwrite(
  candidate_details,
  file.path(
    smr_processed_dir,
    "integration_FAP_hdWGCNA_UKB_candidate_details.csv"
  )
)

# ---------------------------------------------------------------------------
# Revised UpSet
# ---------------------------------------------------------------------------
x <- list(
  Yellow_Module = yellow_module_genes,
  Brown_Module = brown_module_genes,
  Blue_Module = blue_module_genes,
  DEGs_FAPs = DEGs_FAPs_genes,
  SMR_UKB = smr_ukb_bh_genes
)

venn <- Venn(x)

p <- plot_upset(
  venn,
  nintersects = 26,
  order.intersect.by = "size",
  relative_height = 2,
  relative_width = 0.3
)

save_gg(
  p,
  "plot_upset.png",
  w = 8,
  h = 5
)

# ---------------------------------------------------------------------------
# Machine-readable UpSet source resources
# ---------------------------------------------------------------------------
set_sizes <- data.table::data.table(
  set = names(x),
  n_genes = vapply(
    x,
    function(z) length(unique(z)),
    integer(1)
  )
)

data.table::fwrite(
  set_sizes,
  file.path(
    outdir,
    "plot_upset_set_sizes.csv"
  )
)

all_upset_genes <- sort(
  unique(
    unlist(
      x,
      use.names = FALSE
    )
  )
)

membership <- data.table::data.table(
  Gene = all_upset_genes
)

for (nm in names(x)) {
  membership[
    ,
    (nm) :=
      Gene %in% x[[nm]]
  ]
}

data.table::fwrite(
  membership,
  file.path(
    outdir,
    "plot_upset_membership.csv"
  )
)

data.table::fwrite(
  integration_candidates,
  file.path(
    outdir,
    "plot_upset_convergent_candidates.csv"
  )
)

# ---------------------------------------------------------------------------
# Console summary
# ---------------------------------------------------------------------------
cat("\n============================================================\n")
cat("REVISED DISCOVERY-STAGE INTEGRATION\n")
cat("============================================================\n\n")

cat("Yellow-module genes:", length(yellow_module_genes), "\n")
cat("Brown-module genes :", length(brown_module_genes), "\n")
cat("Blue-module genes  :", length(blue_module_genes), "\n")
cat("FAP DEG genes      :", length(DEGs_FAPs_genes), "\n")
cat("UKB SMR genes      :", length(smr_ukb_bh_genes), "\n")
cat("UKB discovery rows :", nrow(ukb_discovery_rows), "\n\n")

cat("Yellow ∩ FAP ∩ UKB:\n")
print(yellow_hits)

cat("\nBrown ∩ FAP ∩ UKB:\n")
print(brown_hits)

cat("\nBlue ∩ FAP ∩ UKB:\n")
print(blue_hits)

cat("\nCandidate details:\n")
print(candidate_details)

cat("\nFinnGen is reserved for targeted replication and is not included in this discovery UpSet.\n")
