rm(list = ls(all.names = TRUE))
gc()

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Provide annotated human RDS path.")
}

rds_file <- args[1]
obj <- readRDS(rds_file)

md <- as.data.table(obj@meta.data)

cell_n <- md[
  ,
  .N,
  by = .(
    sample,
    condition,
    celltype = skeletal_muscle
  )
]

thresholds <- c(1, 5, 10, 20)

out <- rbindlist(lapply(thresholds, function(thr) {

  rbindlist(lapply(sort(unique(cell_n$celltype)), function(ct) {

    x <- cell_n[
      celltype == ct &
      N >= thr
    ]

    data.table(
      min_nuclei = thr,
      celltype = ct,
      Young_donors = uniqueN(x[condition == "Young", sample]),
      Aged_donors  = uniqueN(x[condition == "Aged", sample]),
      analyzable_ge2_each =
        uniqueN(x[condition == "Young", sample]) >= 2 &
        uniqueN(x[condition == "Aged", sample]) >= 2,
      analyzable_ge3_each =
        uniqueN(x[condition == "Young", sample]) >= 3 &
        uniqueN(x[condition == "Aged", sample]) >= 3
    )
  }))
}))

cat("\n============================================================\n")
cat("DONOR ADEQUACY BY MINIMUM NUCLEI PER PSEUDOBULK\n")
cat("============================================================\n\n")

print(out)

cat("\n============================================================\n")
cat("CELL TYPES FAILING >=3 DONORS PER CONDITION\n")
cat("============================================================\n\n")

print(
  out[
    analyzable_ge3_each == FALSE
  ]
)

cat("\n============================================================\n")
cat("RAW DONOR x CELL-TYPE NUCLEUS COUNTS\n")
cat("============================================================\n\n")

wide <- dcast(
  cell_n,
  celltype ~ sample,
  value.var = "N",
  fill = 0
)

print(wide)

outdir <- file.path(
  "processed_results",
  "02_differential_expression",
  "full_precision_audit"
)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

fwrite(
  out,
  file.path(outdir, "pseudobulk_sample_adequacy_thresholds.csv")
)

fwrite(
  cell_n,
  file.path(outdir, "donor_celltype_nucleus_counts.csv")
)
