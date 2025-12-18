# Integrative Multi-Omics Analysis of Aging-Associated Inguinal Hernia

This repository contains the full analysis code supporting a manuscript
investigating the role of aging-associated fibro-adipogenic progenitors (FAPs),
fibrosis, and genetic regulation in inguinal hernia.

The analyses integrate single-cell multi-omics profiling, network biology,
and statistical genetics to identify key regulatory mechanisms linked to aging.

---

## 🔬 Overview of Analyses

The repository includes scripts for:

- Single-cell RNA-seq and ATAC-seq preprocessing and QC
- Cell-type composition and trajectory inference (Monocle3)
- High-dimensional WGCNA (hdWGCNA)
- GWAS data preparation
- Summary-data-based Mendelian randomization (SMR/HEIDI)
- Conditional analysis (GCTA-COJO)
- Colocalization and fine-mapping (COLOC, SuSiE)
- Transcription factor network and motif activity (chromVAR)
- Cell–cell communication analysis (CellChat)
- Mouse validation analyses
- Spatial transcriptomics and K-cross spatial statistics

---

## 📂 Repository Structure

Scripts are organized by analytical step under the `scripts/` directory and
are numbered to reflect execution order.

No individual-level or raw sequencing data are included in this repository.

---

## 🧪 Data Sources

All datasets used in this study are publicly available:

- Human skeletal muscle single-cell multi-omics: GEO 
- GWAS summary statistics: UK Biobank, FinnGen
- Mouse datasets: GEO
- Spatial transcriptomics: GEO

Detailed accession numbers are provided in the manuscript.

---

## ▶️ Reproducibility

Analyses were conducted using R (v4.5.2) alongside command-line tools.
All scripts are provided in this repository, and relevant parameters are explicitly
specified within each script. Random seeds were set where applicable.
