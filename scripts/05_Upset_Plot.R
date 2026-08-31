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


bulk_de_sig_faps <- read.csv(
  file.path(
    repo_root,
    "processed_results",
    "02_differential_expression",
    "bulk_de_sig_faps.csv"
  )
)
DEGs_FAPs_genes <- unique(bulk_de_sig_faps$gene)

smr_res_sig <- read.csv(
  file = file.path(
  repo_root,
  "processed_results",
  "06_SMR_HEIDI",
  "SMR_all_studies_all_tissues_sig.csv"
),
  header = TRUE
)
smr_res_sig_UKBB <- smr_res_sig %>% subset(study == 'UKB')
smr_res_sig_Finn <- smr_res_sig %>% subset(study == "FinnGen")
smr_res_sig_UKBB_genes <- unique(smr_res_sig_UKBB$Gene)
smr_res_sig_Finn_genes <- unique(smr_res_sig_Finn$Gene)



x <- list(Yellow_Module = yellow_module_genes, 
          Brown_Module = brown_module_genes,
          Blue_Module = blue_module_genes, 
          DEGs_FAPs = DEGs_FAPs_genes,
          SMR_UKB = smr_res_sig_UKBB_genes,
          SMR_FinnGen = smr_res_sig_Finn_genes
          )

venn = Venn(x)
p <- plot_upset(venn, 
           nintersects = 26,
           order.intersect.by = "size",
           relative_height = 2, 
           relative_width = 0.3)

save_gg(p,"plot_upset.png", w = 8, h = 5)

int1 <- intersect(DEGs_FAPs_genes,brown_module_genes)
int2 <- intersect(int1,smr_res_sig_UKBB_genes)
int3 <- intersect(int2,smr_res_sig_Finn_genes)
int3
