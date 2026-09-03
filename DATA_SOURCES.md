# Data Sources

This document describes the external datasets and reference resources required
to reproduce the analyses reported in the IH-Aging-Multiomics study.

Raw or individual-level data are not redistributed through this repository.
Users should obtain third-party datasets directly from their original
repositories or data providers and follow the applicable access and
redistribution conditions.

---

## 1. Human skeletal muscle single-nucleus multi-omics

### Dataset

- GEO accession: **GSE268953**
- Biological material: human skeletal muscle
- Conditions: young and aged
- Number of biological replicates: 5 young and 5 aged
- Platform: 10x Genomics Chromium Single Cell Multiome
- Modalities:
  - single-nucleus RNA sequencing (snRNA-seq)
  - single-nucleus ATAC sequencing (snATAC-seq)

### Use in this study

This dataset was used for:

- RNA and ATAC quality control
- ambient RNA decontamination
- doublet detection
- cell-cycle assessment
- multimodal integration and WNN clustering
- cell-type annotation
- pseudobulk differential expression
- cell-type composition analysis
- FAP trajectory inference
- hdWGCNA
- transcription factor regulatory network analysis
- chromVAR motif activity analysis
- CellChat analysis
- SMAD3 pseudotime analysis

### Data availability

The gene-expression and peak-level resources associated with GSE268953 should
be obtained from the original GEO record.

The ATAC fragment files used in this study are not publicly available. These
files were provided directly by Dr. Huating Wang and Dr. Yang Li, investigators
of the original study. Because these are third-party resources, they are not
redistributed through this repository.

### Genome build

The human single-nucleus multi-omics analyses were harmonized to **hg38**.

---

## 2. Inguinal hernia GWAS: discovery dataset

### Dataset

- Source: UK Biobank-based inguinal hernia GWAS
- Role in this study: discovery GWAS
- Reference publication: Fadista et al., Nature Communications, 2022
- Total sample size used in the processed summary statistics: **371,810**

### Use in this study

The UK Biobank-based GWAS summary statistics were used for:

- SMR/HEIDI
- regional association analysis at the SMAD3 locus
- GCTA-COJO
- standard colocalization
- SuSiE fine-mapping
- coloc.susie
- locus and eQTL visualization

### Processing performed in this study

The publicly available GWAS file contained genomic-control-adjusted P-values.
Unadjusted test statistics were reconstructed as described in the Supplemental
Methods using a genomic inflation factor of:

- lambda GC = **1.084**

For each variant, the analysis retained or derived:

- chromosome
- genomic position
- effect allele
- non-effect allele
- effect allele frequency
- regression coefficient
- standard error
- P-value
- sample size

The processed summary statistics used by the analysis scripts are
author-generated derived resources and are provided separately from the
original third-party GWAS data where redistribution is appropriate.

---

## 3. Inguinal hernia GWAS: replication dataset

### Dataset

- Source: **FinnGen**
- Role in this study: independent replication GWAS
- Reference publication: Kurki et al., Nature, 2023
- Total sample size used in the processed summary statistics: **207,653**

### Use in this study

FinnGen summary statistics were used for:

- SMR/HEIDI replication
- regional association analysis at the SMAD3 locus
- GCTA-COJO
- standard colocalization
- SuSiE fine-mapping
- coloc.susie
- locus and eQTL visualization

### Processing performed in this study

Variant-level effect sizes, standard errors, -log10 P-values, and allele
frequencies were parsed from the source summary-statistic file as described in
the Supplemental Methods.

When the primary allele-frequency field was missing, the AF value from the
INFO field was used. P-values were reconstructed from the reported -log10
P-values.

The alternate allele was used as the effect allele and the reference allele as
the non-effect allele in the harmonized analysis files.

---

## 4. GTEx v8 cis-eQTL resources

### Source

- Resource: **GTEx v8**
- Molecular trait: cis-eQTL summary statistics

### Tissues analyzed

The following four tissues were included:

1. Adipose Subcutaneous
2. Adipose Visceral Omentum
3. Muscle Skeletal
4. Cells Cultured Fibroblasts

### Use in this study

GTEx v8 cis-eQTL data were used for:

- SMR/HEIDI
- standard colocalization at the SMAD3 locus
- SuSiE-based colocalization

For SMR analyses, cis-eQTLs were restricted to variants within 1 Mb of the
transcription start site, as described in the manuscript.

Third-party GTEx resources are not redistributed through this repository and
should be obtained from the original provider.

---

## 5. Linkage disequilibrium reference panel

### Source

- Resource: **1000 Genomes Project Phase 3**
- Population: European ancestry

### Use in this study

The European-ancestry reference panel was used for:

- SMR/HEIDI linkage disequilibrium calculations
- GCTA-COJO
- SuSiE fine-mapping
- regional LD analyses
- locus visualization

For the SMAD3 regional analyses, variants were restricted to those overlapping
the relevant GWAS locus and the LD reference dataset.

Third-party 1000 Genomes genotype files are not redistributed through this
repository.

---

## 6. Mouse inguinal hernia single-nucleus multi-omics dataset

### Dataset

- GEO accession: **GSE288662**
- Species: mouse
- Disease context: experimental inguinal hernia model

### Use in this study

The processed single-cell resource was used for:

- WNN-based visualization and cell annotation
- pseudobulk SMAD3 expression analysis
- chromVAR analysis of:
  - ETS1
  - ETS2
  - FOS
  - SMAD3

Comparisons were performed between:

- EP: estradiol + progesterone
- Veh: vehicle

The original dataset should be obtained from GEO and is not redistributed
through this repository.

---

## 7. Mouse inguinal hernia spatial transcriptomics dataset

### Dataset

- GEO accession: **GSE288663**
- Species: mouse

### Samples analyzed

- 1450Veh
- 1603Veh
- 1460EP
- 1461EP

### Use in this study

Spatial transcriptomic data were used to assess:

- spatial localization of SMAD3 transcripts
- spatial localization of predicted positive and negative SMAD3 targets
- Kcross spatial association between SMAD3 and its positive targets
- spatial density of SMAD3-positive target transcripts

Transcript-level data were filtered using:

- quality value (qv) > 20

The original spatial dataset should be obtained from GEO and is not
redistributed through this repository.

---

## 8. Genome annotation and regulatory resources

The analyses additionally used the following external reference resources:

### Human genome annotation

- UCSC hg38
- EnsDb.Hsapiens.v86

### Transcription factor motifs

- JASPAR 2024

### ENCODE blacklist

ENCODE blacklist regions were used during ATAC-seq quality control.

These third-party reference resources are not redistributed unless explicitly
permitted by their respective providers.

---

## 9. Data redistribution policy

This repository distinguishes between:

### Third-party input data

Original datasets obtained from GEO, GTEx, FinnGen, the UK Biobank-based GWAS,
the 1000 Genomes Project, and other external providers are subject to the
access and redistribution conditions of their original sources.

Where redistribution is not appropriate, this repository provides
documentation describing the required resource rather than redistributing the
original files.

### Author-generated processed resources

Author-generated summary-level results underlying the analyses, figures, and
tables are provided in the `processed_results/` and `source_data/`
directories where redistribution is appropriate.

These resources are intended to facilitate reproduction of the reported
results without duplicating restricted or externally maintained datasets.

---

## 10. Expected local data organization

Analysis scripts assume that third-party input datasets are stored locally
outside version control.

A recommended organization is:

```text
data/
├── GSE268953/
│   ├── RNA/
│   ├── ATAC/
│   └── Fragments/
├── GWAS/
│   ├── UKB/
│   └── FinnGen/
├── GTEx_v8/
├── 1000Genomes/
├── GSE288662/
└── GSE288663/
```

The `data/` directory should not be committed to the repository.

Exact input filenames required by individual analyses are documented in the
analysis scripts and are also summarized in `RUN_ORDER.md`.
