rm(list = ls(all.names = TRUE))
gc()

suppressPackageStartupMessages({
  library(CellChat)
})

obj_file <- file.path(
  "processed_results",
  "12_CellChat",
  "cellchat_merged.rds"
)

x <- readRDS(obj_file)

cat("\n===== CLASS / SIZE =====\n")
print(class(x))
print(object.size(x), units = "MB")

cat("\n===== DATA SLOTS =====\n")
cat("data dim          : ")
print(dim(x@data))

cat("data.signaling dim: ")
print(dim(x@data.signaling))

cat("\n===== META =====\n")
cat("meta dim: ")
print(dim(x@meta))

cat("meta columns:\n")
print(colnames(x@meta))

cat("\nFirst rows of meta:\n")
print(head(x@meta))

cat("\n===== DATASET COLUMN CANDIDATES =====\n")
for (nm in c(
  "datasets",
  "dataset",
  "condition",
  "group_cellchat",
  "cc_group",
  "labels"
)) {
  if (nm %in% colnames(x@meta)) {
    cat("\n", nm, ":\n", sep = "")
    print(table(x@meta[[nm]], useNA = "ifany"))
  }
}

cat("\n===== IDENTITIES =====\n")
cat("class(x@idents):\n")
print(class(x@idents))

if (is.list(x@idents)) {
  cat("names(x@idents):\n")
  print(names(x@idents))

  for (nm in names(x@idents)) {
    cat("\nidents$", nm, ":\n", sep = "")
    print(table(x@idents[[nm]], useNA = "ifany"))
  }
} else {
  print(table(x@idents, useNA = "ifany"))
}

cat("\n===== CELL NAME CONSISTENCY =====\n")
cat("ncol(data.signaling): ", ncol(x@data.signaling), "\n", sep = "")
cat("nrow(meta)          : ", nrow(x@meta), "\n", sep = "")

if (ncol(x@data.signaling) > 0 && nrow(x@meta) > 0) {
  cat(
    "colnames(data.signaling) identical to rownames(meta): ",
    identical(
      colnames(x@data.signaling),
      rownames(x@meta)
    ),
    "\n",
    sep = ""
  )
}

cat("\n===== SIGNALING GENES =====\n")
cat(
  "Number of signaling genes: ",
  nrow(x@data.signaling),
  "\n",
  sep = ""
)

cat("First 30:\n")
print(head(rownames(x@data.signaling), 30))

cat("\n===== STORED NETWORKS =====\n")
print(names(x@net))

for (ds in names(x@net)) {
  cat("\n", ds, " prob dim: ", sep = "")
  print(dim(x@net[[ds]]$prob))

  cat(ds, " pval dim: ")
  print(dim(x@net[[ds]]$pval))
}

cat("\nDone.\n")
