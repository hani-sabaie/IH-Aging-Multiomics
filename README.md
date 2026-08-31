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
- Cell-cell communication analysis (CellChat)
- Mouse validation analyses
- Spatial transcriptomics and K-cross spatial statistics

---

## 📂 Repository Structure

Scripts are organized by analytical step under the `scripts/` directory and
are numbered to reflect execution order.

No individual-level or raw sequencing data are included in this repository.

---

## 🧪 Data Sources

The analyses in this study use data obtained from the original repositories
and data providers described below.

### Human skeletal muscle single-nucleus multi-omics

Paired single-nucleus RNA-seq and ATAC-seq data from young and aged human
skeletal muscle were obtained from GEO accession **GSE268953**.

The ATAC fragment files associated with GSE268953 are not publicly available
and were provided directly by the investigators of the original study,
Dr. Huating Wang and Dr. Yang Li. These third-party files are therefore not
redistributed through this repository.

### Genome-wide association data

Inguinal hernia GWAS summary statistics were obtained from:

- UK Biobank-based inguinal hernia GWAS (discovery dataset)
- FinnGen (independent replication dataset)

### eQTL data

GTEx v8 cis-eQTL summary statistics were used for the following tissues:

- Adipose Subcutaneous
- Adipose Visceral Omentum
- Muscle Skeletal
- Cells Cultured fibroblasts

### Linkage disequilibrium reference

European-ancestry samples from the 1000 Genomes Project Phase 3 were used
as the LD reference panel.

### Mouse inguinal hernia data

Mouse single-nucleus multi-omics and spatial transcriptomic datasets were
obtained from GEO accessions **GSE288662** and **GSE288663**.

Detailed information on data acquisition and required input files will be
provided in `DATA_SOURCES.md`. Third-party datasets are not redistributed
when their original access or redistribution conditions do not permit it.

---

## ▶️ Reproducibility

Analyses were conducted using R (v4.5.2) together with command-line tools.
The repository contains the analysis scripts used for the computational
workflows described in the manuscript. As part of the revised reproducibility
package, the repository is being organized to provide explicit documentation
of data acquisition, software dependencies, analysis execution order, and
author-generated processed results underlying the reported analyses.

Third-party datasets and resources subject to access or redistribution
restrictions are not redistributed through this repository.


## License

The analysis code in this repository is released under the MIT License.
See [LICENSE](LICENSE) for details.

Third-party datasets and resources used by these workflows remain subject
to their original licenses, access conditions, and redistribution terms.
See [DATA_SOURCES.md](DATA_SOURCES.md) for data-source details.
