# Verification

Validation uses Python 3.11.14 on Apple Silicon macOS with `environment-osx-arm64.lock` plus `requirements-pip.txt`. The environment was installed independently of the analysis workstation environment. Times New Roman is installed separately.

The current cached workflow passes 13 failure/invariance tests, the 49-value historical check, 36-contrast primary regression check and 369 analysis checks. It regenerates all 19 figures, Source Data, Supplementary Tables 3 and 4, and their numerical plotting tables. The figures' internal checks validate panel identity, aspect ratios, placement, source hashes and content. The completed reconstruction measurements appear below under 20 September 2026.

A separate refit creates new caches from all six sets of raw OD/plate-layout workbooks and starts with empty generated-output directories, then runs the same analysis and figure workflow. It compares the resulting CSV tables, PNGs and workbook cells with the reference checkout. The cache comparator checks every column, including all nested fit/bootstrap arrays, at relative tolerance 1e-8 and absolute tolerance 1e-10, with matching NaN locations and exact nonnumeric values. Refit outputs are isolated from the supplied caches.

The pip-only environment was also tested. Different native FreeType and numerical-library builds changed glyph rendering and some optimizer results, despite identical Python package versions. The documented default therefore uses the native-library lock. Cross-platform builds from `environment.yml` have not been validated.

Detailed completed-run measurements are recorded in `validation.json`. Runtime logs are kept under `runs/` and are not part of the distribution. Figure coverage and unresolved experimental records are stated in the README and `FIGURE_COVERAGE.md`.

## Initial validation, 17 September 2026

The cached build and isolated raw refit at `f6026fc` passed the automated comparison of 24 CSV tables, 125 PNGs and 414,279 workbook cells, including 138,959 numeric cells. A clean-checkout build independently confirmed regeneration from empty output directories.

The publication-language update preserved all 24 tables, all numeric workbook values and all ten figure assemblies. Twenty-six workbook text cells use publication terminology. Of the 125 PNGs, 124 are byte-identical to the baseline; the standalone Supplementary Figure 11 plot has a wording change confined to its plating-volume note. The plotted observations and geometry are unchanged.

## Source corrections, 18 September 2026

The confirmed S3/S6 pairing changes six primary contrasts; all other 30 primary contrasts are unchanged. The six revised results were independently checked from paired log10 ratios. The recovered Figure 1f observation is verified in both supplied workbooks and included in the canonical inputs, per-culture export, daily summary and plotted mean.

A complete cached rebuild reproduced all 26 CSV tables, all 126 PNGs and 417,798 workbook cells, including 140,622 numeric cells. This run includes the new `survival.py` step. Raw dose-response data and fitted caches are unchanged; the earlier isolated-refit validation remains applicable to those fits. No new raw refit was needed for these corrections.

The subsequent Supplementary Figure 6 day-label correction changes only the retained heatmap label from 81 to 82; the assembly checks passed after rebuilding with that input. The numerical heatmap is unchanged. [Current validation measurements](validation_2026-09-18.json) and [source corrections](DATA_CORRECTIONS.md) record the checks.

## Genomic input integration, 20 September 2026

The supplied archive contains every requested file: 149 mutation-call TSVs and 17 processed CSVs, with no missing or extra genomic inputs. All input hashes remain unchanged after regeneration. Eight tables regenerate from the original notebook data-processing cells. The PbEc table contains all 83 rows; the other five endpoint tables contain 43, 64, 80, 56 and 38 rows. PLAC contains 601 rows after the notebook's allele corrections, and the combined-lineage table contains 963 rows.

All 17 processed exports agree with the rerun calculations before later notebook corrections, and all 331 previously saved endpoint rows agree in mutation identity, frequencies and untruncated annotations. Display-only ellipses, empty-field representations and numerical serialization are normalized for this comparison. The figure-specific table mappings are checked against duplicate notebook data-processing cells. These checks are executable in `scripts/reproduce_genomics.py`; detailed results are in `data/genomics/generated_tables/validation.json`.

The complete cached workflow passed in 156 seconds. All 126 existing PNG files and all non-genomic numerical CSV outputs remain unchanged. The Source Data exporter verified 147,296 numeric cells, 80,279 text cells and 16,220 blank cells across 25 sheets and 89 blocks. Checks of both workbooks confirm that all pre-existing non-genomic data blocks are unchanged. A second genomic run with a different Python hash seed produced 18 identical output files. Six representative genomic table views were inspected with wrapping enabled for full annotations. [Validation measurements](validation_2026-09-20.json) record these results.

## Complete numerical figures and single-cell reconstruction, 20 September 2026

The completed workflow passes all 13 failure/invariance tests, 369 analysis checks, the historical notebook check and the 36-contrast primary regression check. The final cached run completed in 208 seconds and generated all 19 figures across 26 PDF pages, including 50 numerical placements in the ten established assemblies. Numerical and schematic sources are distinguished in the figure map. Genomic figure text boundaries, matrix identities, original chart annotations, assembly geometry and source hashes are checked automatically.

The final statistics, group summaries, sensitivity results and Figure 5 descriptive tables are byte-identical to the preceding committed analysis. Figure 3 zero-count handling and S11 provenance are updated as documented in DATA_CORRECTIONS.md. Source Data verifies 499,906 numeric cells, 273,910 text cells and 30,072 blank cells across 25 sheets and 111 blocks.

The independent R count-matrix execution reproduces 48,883 cells and every cluster count, followed by complete cluster/sample differential-expression calculations. All 31 selected source-file hashes are verified. The 6,543 cluster DGE rows agree within 1.1e-14 in log2 fold changes; the largest culture-comparison discrepancy is 7.61e-9 log10-P units. Targeted hipA/rplJ checks were also executed from the saved Seurat object. The full package lock and R session record are included. EcoCyc method provenance was subsequently supplied and is recorded below.

[Execution measurements](validation_execution_2026-09-20.json) record the completed checks. The earlier raw dose-response refit validation still applies to the unchanged source OD/plate-layout workbooks and fitted caches.

## Clean output build and additional checks, 20 September 2026

A separate checkout with empty generated-output directories passed in 214.33 seconds after correcting the build order: genomic inputs are now regenerated before the tests that consume them. The build reproduced all 28 numerical CSV tables, 186 PNGs, 52 generated data files and 1,322,761 Source Data worksheet cells, including 499,906 numeric cells. The worksheet-cell total includes empty spacers outside the explicitly exported blocks.

An independent calculation directly from the H5 matrix verified all 48,883 retained cell identities, the 4,200 selected probes and eight available marker comparisons using SciPy rank-sum tests rather than Seurat marker functions. It retained the original Bonferroni denominator of 21,701 assay rows.

The final Figure 6 legend distinguishes the 0.5 plotting cutoff from the original 0.6 shared-gene cutoff. An executable check exactly recovers all 91 Shared annotations from the two culture comparisons. Both Figure 6 PNG changes are confined to the legend; the other 184 PNGs, all 28 numerical CSV tables and Source Data are unchanged. All 82,275 cells in Supplementary Tables 3/4 are unchanged. The added worksheet-dimension and gene-name checks pass, as do all 13 failure/invariance tests. [Measurements](validation_clean_build_2026-09-20.json) record the results.

## Reference and EcoCyc records, 21 September 2026

The reference/run-record verifier checks all 600 archive members across 453 distinct stored files and all 149 study-run identities. Archived GenBank sequences match every saved FASTA: 141 U00096.3 runs and eight ATEC runs. Both saved summaries agree for every run. The exact executable version is unrecorded.

All five original EcoCyc export hashes pass. The 746 rows overlapping the source workbook agree in term, P-value and matched genes. The complete parent-up export adds 114 rows, yielding 133 results in Table 4 and Source Data. Existing worksheet cells remain unchanged. No enrichment tests were recalculated.

The full cached workflow passed in 210.01 seconds, including all 13 tests, the primary regression and historical checks. All 201 tracked PNG files and 28 numerical CSV outputs are byte-identical to the preceding commit. Source Data validates 500,020 numeric, 274,137 text and 31,099 blank exported cells across 25 sheets and 111 blocks. All 97 top-level input checksums pass. Unchanged figure containers are retained to avoid changes solely from export timestamps and generated object identifiers.
