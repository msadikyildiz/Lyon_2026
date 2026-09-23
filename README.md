# Lyon et al. (2026): data, analysis and figures

Data and analysis code for *Recurrent Extinction of Resistance Mutations Leads to Convergent Multidrug Resistance in Sequential Antibiotic Treatment*. Maintained by Adam Lyon and Muhammed Sadik Yildiz.

## TODO

| Task | Input or work remaining | Affected outputs |
|---|---|---|
| Add the sequencing accession | Add the WGS accession. All 149 runs include commands, saved options and sequence-verified references. Date-based software-version assumptions are documented with the run records. | Access to sequencing reads |
| Update the schematic | Add the final editable Figure 4d BioRender artwork and export. | Figure 4d |

The repository includes the experimental workbooks, mutation calls, single-cell count matrices and original analysis code. One command regenerates the numerical analyses, Source Data, Supplementary Tables 3 and 4, and all 19 figures. Schematic artwork is retained explicitly; every numerical panel is generated from included data. See the [figure map](docs/FIGURE_COVERAGE.md).

## Run

On Apple Silicon macOS, install [micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html), then run from the repository root:

```sh
micromamba create -n lyon-2026 -f environment-osx-arm64.lock
micromamba run -n lyon-2026 python -m pip install -r requirements-pip.txt
micromamba run -n lyon-2026 python reproduce.py
```

The default run uses the included dose-response fits and computed single-cell embeddings/DGE tables, checks them against the source tables, and rebuilds the statistics, data exports and figures. Genomic tables are recalculated from mutation calls on every run. Logs and an explicit passed/failed record are written under `runs/`.

The native-library lock controls numerical and font-rendering dependencies. `environment.yml` specifies versions for other platforms, which have not been validated. The pinned plategig source is imported directly; its separate upstream installation is unnecessary. Install Times New Roman separately to reproduce the reference typography. Fonts are not redistributed.

### Refit the dose-response data

```sh
micromamba run -n lyon-2026 python reproduce.py --refit --jobs 2
```

This creates an isolated copy under `runs/`, rebuilds six datasets from raw OD and plate-layout workbooks, compares every fitted/cache field, executes the complete workflow, and compares the resulting tables, workbook cells and PNGs with the reference checkout. Supplied caches remain unchanged. `--output /path/to/new-directory` selects a different new location.

### Recompute single-cell results from count matrices

Use R 4.4.2 with the packages recorded in `single-cell-renv.lock`. With `renv` installed, restore an isolated library and run:

```sh
Rscript -e 'renv::restore(lockfile="single-cell-renv.lock", library="runs/R-library", prompt=FALSE)'
R_LIBS_USER="$PWD/runs/R-library" micromamba run -n lyon-2026 python reproduce.py --single-cell
```

This command replaces the generated single-cell tables in the checkout. Use a separate checkout to compare a different R environment. The R step executes the supplied filtered-analysis Rmd calculations with portable input paths and exports cell identities, embeddings, selected probes, normalized expression and differential-expression results. It reproduces 48,883 cells and the exact nine cluster sizes. All 6,543 cluster DGE rows and the 1,222/863 culture-versus-parent rows agree with the reference expression tables within the recorded numerical tolerances. The [single-cell record](data/single-cell/README.md) explains the sample mapping, table numbering, targeted marker checks and enrichment inputs and assumptions.

## Outputs and inputs

- [All figures](working/figures-assembled/All_figures.pdf), with individual PDF and PNG files alongside it. Larger genomic figures use multiple pages at readable type sizes.
- [Source Data](working/source-data/Source%20Data.xlsx), [Supplementary Table 3](working/source-data/Supplementary%20Table%203.xlsx) and [Supplementary Table 4](working/source-data/Supplementary%20Table%204.xlsx).
- [Final statistics](working/analysis/stats-rework/out/final_statistics.csv), [analysis code](working/analysis/stats-rework), [original workbooks/notebooks](working/figures), [genomic inputs](data/genomics) and [single-cell inputs](data/single-cell).
- [Data processing and source records](docs/DATA_CORRECTIONS.md), [figure coverage](docs/FIGURE_COVERAGE.md) and [verification](docs/VALIDATION.md).

The single-cell inputs occupy approximately 196 MB and include count matrices, analysis code, expression tables, alternative aggregates and excluded sample-5 matrices. Archive paths and per-file checksums identify each source file.

## Analysis conventions

Primary IC50 comparisons use log10 values, culture-ID pairing for Figures 2/3 and Supplementary Figures 3/6, Welch tests for Supplementary Figure 4, and Holm adjustment within each panel and antibiotic. MIC and biological-n=1 mutant measurements are descriptive. Fit-perturbation ranges are sensitivity summaries. Fixed seeds and complete paired IDs are checked. A separate reference check verifies 49 notebook results, including six expected NaNs, separately from the primary analysis.

Zero-colony plates are distinguished from missing counts. Figure 3 parent medians and MADs are invariant over the two censored intervals; omitting the potentially mishandled 3-hour record raises the median by 4%. Supplementary Figure 11 uses a single 10 µL plating at every time point, so equal volume cancels in normalized fractions and detection limits. The recorded 40 µL exceptions and ATEC repeat/exclusion reasons are documented in the data-processing record.

Genomic plots use coordinate/allele identity rather than averaging different variants under the same short label. Each plotted entry has a source crosswalk. Single-cell tests use cells as observations; the six libraries comprise two technical replicates from each of three cultures. Four of the six enrichment worksheets have matching EcoCyc text exports; the cluster-1/2 results remain supplied workbook inputs. The [enrichment record](data/single-cell/enrichment/README.md) describes the full parent-up export and database-version and default-background assumptions. GSE314756 is scheduled for public release on 21 December 2026.

## Compare outputs

```sh
micromamba run -n lyon-2026 python scripts/compare_outputs.py /path/to/reference /path/to/rebuilt --report comparison.json
```

The comparator checks CSV values, PNG bytes, generated genomic/single-cell plotting tables, and Source Data cell values/types. Statistical/workbook tolerances are 1e-8 relative and 1e-10 absolute. PDF, SVG and XLSX containers may differ in timestamps. R reconstruction is validated separately against the original expression tables, including exact row identities and cluster counts.

## License

Project-authored analysis and reproduction code and its documentation are available under the [MIT license](LICENSE).

Study-generated experimental data, derived result tables and original figure elements are licensed under [Creative Commons Attribution 4.0 International (CC BY 4.0)](LICENSE-DATA). Reuse requires appropriate credit, a link to the license and an indication of modifications, where applicable. Academic and commercial reuse are permitted. Credit Lyon et al. (2026) and link to this repository.

Third-party software, reference sequences, database annotations and BioRender assets retain their existing terms and are excluded from these grants. For composite figures, CC BY 4.0 covers only the original study content. The vendored plategig source retains its [MIT license](working/analysis/stats-rework/plategig-88839a2/LICENSE).

## Citation

Please cite the associated study and this repository when using the data or code, and include the Git commit identifier for the version used. [CITATION.cff](CITATION.cff) provides the repository citation metadata. This citation guidance does not add conditions to the MIT license.
