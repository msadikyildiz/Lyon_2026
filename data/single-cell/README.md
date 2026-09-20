# Single-cell data and reconstruction

Samuel Cooke's analysis files, supplied by Adam Rosenthal on 20 September 2026, contain the count matrices, R/Seurat code, GEO metadata templates, full expression tables and Excel volcano plots. `input_manifest.json` identifies each selected file in the original archive, its checksum, and the duplicate copies verified before import. Original inputs remain unchanged. The repository contains 31 selected files (195,521,222 bytes); repeated plots, project copies, Git history and editor/session files are omitted.

## Samples and calculations

The final analysis uses `counts/Filtered_Features/aggr_5removed_filtered_feature_bc_matrix.h5`. Barcode suffixes 1–6 map to 3a, 3b, 4a, 4b, 7a and 7b. Samples 3a/3b are the parent; 4a/4b and 7a/7b are evolved cultures 4 and 7. Each pair is two technical replicates from the same sampling, separately probed and processed through microfluidics. `plot_tables/sample_metadata.csv` retains these identities. Original GEO template files are preserved as supplied; the analysis crosswalk uses the source Rmd sample mapping.

Minimum UMI counts are 27, 27, 60, 70, 65 and 35, respectively. Retained cells are 7,827, 10,774, 4,764, 2,764, 9,619 and 13,135, totaling 48,883. The supplied code retains the highest-total-count probe per gene, log-normalizes, selects variable features, scales, and uses principal components 1–7 for neighbor finding and UMAP. Louvain clustering uses resolution 0.7. The original cluster IDs 0–8 map to published clusters 1–9; their sizes are 9,335, 7,034, 6,008, 5,989, 5,092, 4,967, 4,961, 3,052 and 2,445.

`scripts/reproduce_single_cell.R` executes the relevant original Rmd chunks and exports the computed object under ignored `runs/single-cell/`, plus portable CSV results under `generated/`. The reconstruction uses R 4.4.2, Seurat 5.2.1 and SeuratObject 5.0.2. The original report identifies Seurat 5 but does not record its exact minor version. The complete reconstruction environment is in `single-cell-renv.lock`; `generated/R-session.txt` records the executed session.

## Validation and plots

`scripts/reproduce_single_cell.py` verifies every selected input hash and compares the recomputed full cluster and sample DGE tables with the supplied originals. It requires exact row identities and cluster sizes, fold-change/detection-fraction agreement within 1e-12, identical zero-P patterns, and log10-P agreement within 1e-7. The largest observed difference is 7.61e-9 log10 units. Results are in `plot_tables/validation.json`.

Figure 6 volcano points use absolute log2 fold change > 0.5 and Bonferroni-adjusted P < 0.05. Supplementary Figures 12/13 use > 0.25 and the same significance threshold. Their 674, 427, 734 and 1,175 points are matched individually to the recomputed cluster tables. Highlight colors are recovered from the supplied Excel chart point annotations, including manual gene-display aliases. Values below 10^-304 are capped only for plotting. UMAP scatter points are rasterized within otherwise vector PDFs; no cells are subsampled.

Targeted checks of the two genes discussed individually, hipA and rplJ, remove the 1% detection-frequency filter while retaining Bonferroni correction over all assay features. hipA is excluded from the default culture-7 DGE table because it is detected in 0.6% of culture-7 cells and 0.2% of parent cells. Its targeted comparison gives log2 fold change 0.721 and adjusted P = 0.00377. rplJ has higher expression in the parent than either evolved culture, but its cluster-1/2 comparisons have adjusted P = 1. hipA is higher in cluster 9 and has a small negative fold change in cluster 4. The exact comparisons are included in Source Data and Supplementary Table 4.

## Supplementary tables

| Table | Contents |
|---|---|
| 3 | Sample statistics, full cluster DGE and enrichment tables for clusters 1, 2, 4 and 9 |
| 4 | Full culture-4/7-versus-parent DGE, shared-gene enrichment tables and separately identified targeted marker checks |

The exporter preserves all cells in the relevant supplied worksheets, within floating-point serialization precision, and adds the targeted checks as a separate worksheet. `plot_tables/table_mapping.json` makes the numbering explicit.

The enrichment worksheets include results and matched gene lists. Neither supplied Rmd file contains their generating method, database version or tested background gene list. These are required to regenerate enrichment, rather than only reproduce the supplied results. GSE314756's public endpoint reported a 21 December 2026 release date when checked on 20 September 2026; no reviewer credential is stored here.
