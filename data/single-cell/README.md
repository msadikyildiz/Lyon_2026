# Single-cell data and reconstruction

Samuel Cooke's analysis files, supplied by Adam Rosenthal on 20 September 2026, contain the count matrices, R/Seurat code, GEO metadata templates, full expression tables and Excel volcano plots. `input_manifest.json` identifies each selected file in the original archive, its checksum, and the duplicate copies verified before import. Original inputs remain unchanged. The repository contains 31 selected files (195,521,222 bytes); repeated plots, project copies, Git history and editor/session files are omitted.

## Samples and calculations

The final analysis uses `counts/Filtered_Features/aggr_5removed_filtered_feature_bc_matrix.h5`. Barcode suffixes 1–6 map to 3a, 3b, 4a, 4b, 7a and 7b. Samples 3a/3b are the parent; 4a/4b and 7a/7b are evolved cultures 4 and 7. Each pair is two technical replicates from the same sampling, separately probed and processed through microfluidics. `plot_tables/sample_metadata.csv` retains these identities. Original GEO template files are preserved as supplied; the analysis crosswalk uses the source Rmd sample mapping.

Minimum UMI counts are 27, 27, 60, 70, 65 and 35, respectively. Retained cells are 7,827, 10,774, 4,764, 2,764, 9,619 and 13,135, totaling 48,883. The supplied code retains the highest-total-count probe per gene, log-normalizes, selects variable features, scales, and uses principal components 1–7 for neighbor finding and UMAP. Louvain clustering uses resolution 0.7. The original cluster IDs 0–8 map to published clusters 1–9; their sizes are 9,335, 7,034, 6,008, 5,989, 5,092, 4,967, 4,961, 3,052 and 2,445.

`scripts/reproduce_single_cell.R` executes the relevant original Rmd chunks and exports the computed object under ignored `runs/single-cell/`, plus portable CSV results under `generated/`. The reconstruction uses R 4.4.2, Seurat 5.2.1 and SeuratObject 5.0.2. The original report identifies Seurat 5 but does not record its exact minor version. The complete reconstruction environment is in `single-cell-renv.lock`; `generated/R-session.txt` records the executed session.

## Validation and plots

`scripts/reproduce_single_cell.py` verifies every selected input hash and compares the recomputed full cluster and sample DGE tables with the supplied originals. It requires exact row identities and cluster sizes, fold-change/detection-fraction agreement within 1e-12, identical zero-P patterns, and log10-P agreement within 1e-7. The largest observed difference is 7.61e-9 log10 units. Results are in `plot_tables/validation.json`.

Figure 6 volcano points use absolute log2 fold change > 0.5 and Bonferroni-adjusted P < 0.05. The original shared-gene highlights use a stricter absolute log2 fold change > 0.6 with adjusted P < 0.05 and concordant direction in both default culture-versus-parent comparisons. This rule exactly matches all 91 workbook annotations (82 higher in both evolved cultures, nine higher in the parent); targeted checks do not alter the highlights. Supplementary Figures 12/13 use > 0.25 and the same significance threshold. Their 674, 427, 734 and 1,175 points are matched individually to the recomputed cluster tables. Highlight colors are recovered from the supplied Excel chart point annotations, including manual gene-display aliases. Values below 10^-304 are capped only for plotting. UMAP scatter points are rasterized within otherwise vector PDFs; no cells are subsampled.

Targeted checks of the two genes discussed individually, hipA and rplJ, remove the 1% detection-frequency filter while retaining the supplied Bonferroni correction over 21,701 assay features. The original code selects 4,200 probes, one per gene, and zeros the unused rows rather than removing them. hipA is excluded from the default culture-7 DGE table because it is detected in 0.6% of culture-7 cells and 0.2% of parent cells. Its targeted comparison gives log2 fold change 0.721 and adjusted P = 0.00377. rplJ has higher mean expression in the parent than either evolved culture. Its cluster-1/2 comparisons have nominal P = 0.00626/0.00526 and adjusted P = 1/1.

Both clusters 4 and 9 have a larger fraction of hipA-positive cells than the remaining clusters. Mean normalized expression is higher in cluster 9 and slightly lower in cluster 4 (log2 fold changes 0.597 and −0.136). Detection frequency and mean expression describe different aspects of these sparse distributions; the original qualitative distribution statement is consistent with the higher positive-cell fractions. The cluster comparisons are in Supplementary Table 3; culture comparisons and targeted checks are in Supplementary Table 4. Source Data include both.

## Supplementary tables

| Table | Contents |
|---|---|
| 3 | Sample statistics, full cluster DGE and enrichment tables for clusters 1, 2, 4 and 9 |
| 4 | Full culture-4/7-versus-parent DGE, shared-gene enrichment tables and separately identified targeted marker checks |

The exporter preserves all cells in the relevant supplied worksheets, within floating-point serialization precision, adds the targeted checks as a separate worksheet, and extends the parent-up enrichment table from 19 to 133 rows using its matching original EcoCyc export. `plot_tables/table_mapping.json` makes the numbering explicit.

The [EcoCyc enrichment record](enrichment/README.md) documents five original exports, their worksheet mappings, exact result checks, and the analyst-reported usual Fisher Exact/Benjamini–Hochberg settings with P < 0.1. The historical database snapshot and analysis-specific reference sets are not saved in those exports. EcoCyc 29.0 is consistent with the May 2025 dates but is not a verified version assignment. GSE314756 is scheduled for public release on 21 December 2026; no reviewer credential is stored here.

The cluster and sample marker exports also retain nominal P < 0.01 (`FindAllMarkers(return.thresh=0.01)`), matching the supplied tables. The culture-versus-parent `FindMarkers` tables have no equivalent return-P filter. The enrichment sheet named “up in 1 and gsea w benj hoch” includes rplJ, but its title does not establish the gene-selection or gene-level adjustment method. Those enrichment inputs remain preserved; EcoCyc method provenance is documented separately.
