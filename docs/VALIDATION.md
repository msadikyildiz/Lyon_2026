# Verification

## Environment and workflow

The reference validation uses Python 3.11.14 on Apple Silicon macOS with `environment-osx-arm64.lock` and `requirements-pip.txt`, installed independently of the analysis workstation environment. Times New Roman is installed separately. The portable setup in `environment.yml` pins OpenBLAS as well as the Python packages; [GitHub Actions](https://github.com/msadikyildiz/Lyon_2026/actions/workflows/reproduce.yml) records full Python workflow checks on Windows, Linux and macOS. R count-matrix reconstruction has been validated on Apple Silicon macOS with R 4.4.2 and `single-cell-renv.lock`.

Numerical validation is separate from rendering validation. PNG file sets and assembled-figure/print-panel dimensions must match; byte differences from fonts or native rendering libraries are listed explicitly. `--strict-images` additionally requires identical PNG bytes. See [Reproduction](REPRODUCING.md) for environment and comparison commands.

The cached workflow runs 24 failure/invariance tests, 368 analysis checks, a 49-value notebook reference check and a 36-contrast primary regression check. It generates all 19 figures across 26 PDF pages, numerical plotting tables, Source Data and Supplementary Tables 3 and 4. Figure checks cover panel identity, aspect ratios, placement, source hashes and content. Clean-output builds verify regeneration without pre-existing outputs.

## Dose-response and statistical checks

An isolated raw refit reconstructs all six datasets from OD and plate-layout workbooks, compares every cache column including nested fit/bootstrap arrays, and runs the full analysis and figure workflow. Numerical comparisons use relative tolerance 1e-8 and absolute tolerance 1e-10, matching NaN locations and exact nonnumeric values. The source caches remain unchanged.

Primary tests verify matched culture IDs, adjustment families and all 36 IC50 contrasts. The notebook reference check independently covers 49 values, including six expected NaNs. Sensitivity outputs identify their resampling units and comparison families. Source Data verification checks values and cell types against the exported blocks.

## Genomic checks

All 149 mutation-call TSVs and 17 processed exports are checked against input hashes. Eight tables are regenerated from notebook data-processing cells. Endpoint row counts are 43, 64, 80, 56, 83 and 38; PLAC contains 601 rows and the combined-lineage table 963. All 17 processed exports agree before subsequent notebook annotation/allele transformations. All 331 available notebook-display rows agree in variant identity, frequencies and full annotations after normalizing display truncation and serialization. A different Python hash seed produces identical genomic outputs. Supplementary Figure 9 uses an explicit display specification so equal gene ranks retain the same order across CPU sorting implementations; shuffled-input tests verify both the matrix and source crosswalk.

The run-record verifier checks 600 archive entries stored as 453 distinct files and all 149 run identities. Archived GenBank sequences match every saved FASTA: 141 U00096.3 runs and eight ATEC runs. Both saved summaries agree for each run. Executable versions are unrecorded; date-based assumptions are documented separately.

## Single-cell and enrichment checks

Independent execution from the count matrix reproduces 48,883 cell identities and every cluster count. All 31 source hashes are checked. All 6,543 cluster DGE rows agree within 1.1e-14 in log2 fold change; the largest culture-comparison discrepancy is 7.61e-9 log10-P units. Exact row identities, zero-P patterns and sample mappings are required.

A separate H5 calculation verifies the 4,200 selected probes and eight available marker comparisons using SciPy rank-sum tests, retaining the Bonferroni denominator of 21,701 assay rows. An executable check matches all 91 Figure 6 shared-gene annotations. Supplementary-table checks compare source cells and worksheet dimensions.

All five EcoCyc export hashes are verified. The 746 rows overlapping the workbook agree in term, P-value and matched genes. Table 4 and Source Data include the complete 133-row parent-up export. Enrichment tests are not rerun.

## Records

Machine-readable validation measurements are retained in `validation*.json` and the genomic and single-cell output directories. They identify the scope and date of each check. Runtime logs under `runs/` are excluded from the distribution. [Figure coverage](FIGURE_COVERAGE.md) and [data-processing notes](DATA_CORRECTIONS.md) document scientific scope and assumptions; the README TODO table lists remaining repository inputs.
