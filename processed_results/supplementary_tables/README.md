# Supplementary Tables

This directory contains canonical CSV versions of Supplementary Tables
S2-S6 assembled from the reproducible processed-result files.

## Table S2

File:
`Table_S2_significant_DEGs_all_cell_types.csv`

Source:
`processed_results/02_differential_expression/bulk_de_sig.csv`

Rows: 1154 significant gene/cell-type records across cell types meeting the formal pseudobulk representation criterion.

## Table S3

File:
`Table_S3_significant_DEGs_FAP_subclusters.csv`

Source:
`processed_results/02_differential_expression/bulk_de_sig_faps.csv`

Rows: 210 significant FAP gene/cell-type records, representing 194 unique genes.

Table S3 is exactly the FAP1-FAP3 subset of Table S2:
- FAP1: 93
- FAP2: 3
- FAP3: 114

FAP4 is excluded from formal pseudobulk differential-expression testing because it does not meet the donor-representation criterion used for the revised analysis.

## Table S4

File:
`Table_S4_hdWGCNA_module_assignments.csv`

Source:
`processed_results/05_hdWGCNA/module_assignment_table.csv`

Rows: 3386 gene/module records.

## Table S5

File:
`Table_S5_SMR_HEIDI_all_studies_all_tissues.csv`

Source:
`processed_results/06_SMR_HEIDI/SMR_all_studies_all_tissues.csv`

Rows: 58215 combined SMR/HEIDI records across the eight
GWAS-by-tissue analyses.

## Table S6

File:
`Table_S6_mouse_Smad3_pseudobulk_results.csv`

Source:
`processed_results/13_mouse_validation/bulk_de_all_contrasts_mouse.csv`

Rows: 30 Smad3 pseudobulk results:
- Veh_vs_EP: 10
- EPR_vs_EP: 10
- Veh_vs_EPR: 10

Mouse pseudobulk differential-expression P values were adjusted using the
Benjamini-Hochberg method across all genes tested within each cell type x
contrast family. The primary significance criterion is BH FDR < 0.05 without
a hard fold-change threshold. The historical BH FDR < 0.05 and
|log2FC| >= 1 criterion is retained separately as an effect-size sensitivity
analysis.

For Smad3 in type IIa myofibers, the primary three-condition model is
non-significant after BH correction (EP versus Veh: log2FC = 1.352,
FDR = 0.060). A separate sensitivity analysis used the same mouse-level
inclusion rule and the same four EP and three Veh mice contributing to the
primary analysis, but restricted the fitted model to the direct EP-versus-Veh
comparison. The effect estimate was nearly identical and reached significance
(log2FC = 1.359, FDR = 0.033). No additional minimum-nuclei eligibility
threshold was applied in this sensitivity analysis. This sensitivity result
does not replace the primary three-condition analysis.

The primary-versus-sensitivity comparison is available in:
`processed_results/13_mouse_validation/Smad3_TypeIIa_primary_vs_EP_Veh_sensitivity.csv`.

Table S1 is the STROBE-MR checklist and is handled separately from
these analysis-result tables.
