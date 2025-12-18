# ========================
# Setting up environment
# ========================
# ===== Clean environment =====
rm(list = ls(all.names = TRUE))
gc()

# ===== Loading relevant libraries =====
library(Seurat)
library(ggplot2)
library(scProportionTest)

# ===== set seed =====
set.seed(1234)

# ===== Helpers =====
figdir <- "../outputs/sc/figs"
save_gg <- function(p, filename, w=7, h=5, dpi=300){
  if (inherits(p, "ggplot") || inherits(p, "patchwork")){
    ggsave(file.path(figdir, filename), p, width=w, height=h, units="in", dpi=dpi, bg="white")
  } else {
    png(file.path(figdir, filename), width=w, height=h, units="in", res=dpi)
    print(p); dev.off()
  }
}

# ============================================================================ #

# ========================
# Load the object
# ========================
obj <- readRDS("../outputs/decont_merged_filt_nodoub_cc_sct_reduc_clust_integ_annot_obj.rds")

# ========================
# Single Cell Proportion Test 
# ========================
prop_test <- sc_utils(obj)

prop_test <- permutation_test(
  prop_test, cluster_identity = "skeletal_muscle",
  sample_1 = "Young", sample_2 = "Aged",
  sample_identity = "condition"
)

p1 <- permutation_plot(prop_test)

save_gg(p1, "Cell_Proportion.png", w= 8, h=5)
