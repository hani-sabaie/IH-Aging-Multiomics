rm(list = ls(all.names = TRUE))
gc()

suppressPackageStartupMessages(library(data.table))

de <- fread(
  "processed_results/02_differential_expression/full_precision_audit/all_celltypes_full_precision.tsv"
)

modules <- fread(
  "processed_results/05_hdWGCNA/module_assignment_table.csv"
)

ukb <- fread(
  "processed_results/06_SMR_HEIDI/multiple_testing_audit/discovery_replication/UKB_all_tissues_with_BH.tsv"
)

fin <- fread(
  "processed_results/06_SMR_HEIDI/multiple_testing_audit/discovery_replication/FinnGen_all_tissues_with_BH.tsv"
)

yellow <- unique(modules[tolower(module) == "yellow", gene_name])
brown  <- unique(modules[tolower(module) == "brown",  gene_name])
blue   <- unique(modules[tolower(module) == "blue",   gene_name])

mods <- list(
  Yellow = yellow,
  Brown = brown,
  Blue = blue
)

ukb_sig <- ukb[
  FDR_cohortwide < 0.05 &
  p_HEIDI > 0.01
]

run_framework <- function(label, fap_types, use_fc_cutoff) {

  x <- de[celltype %in% fap_types]

  if (use_fc_cutoff) {
    x <- x[
      p_val_adj < 0.05 &
      abs(avg_log2FC) >= 1
    ]
  } else {
    x <- x[
      p_val_adj < 0.05
    ]
  }

  deg_genes <- unique(x$gene)

  cat("\n============================================================\n")
  cat(label, "\n")
  cat("============================================================\n")
  cat("Significant rows:", nrow(x), "\n")
  cat("Unique genes:", length(deg_genes), "\n\n")

  candidates <- character()

  for (m in names(mods)) {

    z <- sort(
      intersect(
        intersect(deg_genes, mods[[m]]),
        unique(ukb_sig$Gene)
      )
    )

    cat(
      m,
      "∩ FAP-DEG ∩ UKB-SMR-BH:",
      length(z),
      "\n"
    )

    if (length(z)) print(z)

    candidates <- union(candidates, z)
  }

  candidates <- sort(unique(candidates))

  cat("\nAll discovery candidate genes:\n")
  print(candidates)

  if (!length(candidates)) return(NULL)

  disc <- ukb_sig[
    Gene %in% candidates,
    .(
      Gene,
      tissue,
      p_SMR_UKB = p_SMR,
      FDR_UKB = FDR_cohortwide,
      p_HEIDI_UKB = p_HEIDI
    )
  ]

  fg <- fin[
    ,
    .(
      Gene,
      tissue,
      p_SMR_FinnGen = p_SMR,
      p_HEIDI_FinnGen = p_HEIDI
    )
  ]

  rep <- merge(
    disc,
    fg,
    by = c("Gene", "tissue"),
    all.x = TRUE
  )

  rep[, replication_BH :=
    p.adjust(p_SMR_FinnGen, method = "BH")
  ]

  rep[, replication_pass :=
    replication_BH < 0.05 &
    p_HEIDI_FinnGen > 0.01
  ]

  cat("\nFinnGen targeted replication:\n")
  print(rep)

  invisible(rep)
}

run_framework(
  "A: FAP1-FAP4, within-celltype BH + |log2FC| >= 1",
  c("FAP1","FAP2","FAP3","FAP4"),
  TRUE
)

run_framework(
  "B: FAP1-FAP4, within-celltype BH, FDR-only",
  c("FAP1","FAP2","FAP3","FAP4"),
  FALSE
)

run_framework(
  "C: FAP1-FAP3, within-celltype BH + |log2FC| >= 1",
  c("FAP1","FAP2","FAP3"),
  TRUE
)

run_framework(
  "D: FAP1-FAP3, within-celltype BH, FDR-only",
  c("FAP1","FAP2","FAP3"),
  FALSE
)
