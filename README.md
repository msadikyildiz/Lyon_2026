# Lyon et al. (2026): data, analysis and figures

Data and analysis code for *Recurrent Extinction of Resistance Mutations Leads to Convergent Multidrug Resistance in Sequential Antibiotic Treatment*. Maintained by Adam Lyon and Muhammed Sadik Yildiz.

## TODO

| Task | Input or work remaining | Affected outputs |
|---|---|---|
| Complete genomic inputs | Add per-sample mutation calls or full-precision frequency tables, including PLAC traced-allele data and the 33 missing PbEc rows; retain reference and calling provenance. | Figure 4; Supplementary Figures 4, 8 and 9; upstream genomic analyses |
| Add single-cell inputs and code | Add source objects, sample metadata, embeddings, analysis scripts, and complete differential-expression/enrichment tables; resolve table mapping and sample replication. | Figure 6; Supplementary Figures 12 and 13; Supplementary Table 4 |
| Resolve survival source records | Confirm Supplementary Figure 11 plated volumes, the two blank Figure 3 parent counts, and explicit 40 µL MDK entries. | Survivor fractions and MDK summaries |
| Confirm culture provenance and exclusions | Document parent-culture matching in Supplementary Figures 3/6 and the ATEC-C3 exclusion reason. | Comparison design and exclusion records |
| Complete figure generation | Connect the genomic and single-cell inputs to executable plotting code, replace retained image panels where numerical inputs are available, and automate the remaining figure compositions. | Full figure coverage; see [figure map](docs/FIGURE_COVERAGE.md) |
| Complete data access information | Add the whole-genome sequencing accession and confirm GSE314756 access and sample mapping. | Data availability and sample provenance |
| Specify reuse terms | Add licenses for the study code and data. | Repository license files |

The included workbooks support dose-response, growth and survival analyses. One command regenerates the statistical tables, Source Data workbook, 40 size-specific panels, six trajectory plots and ten assembled figures. The [figure map](docs/FIGURE_COVERAGE.md) documents executable coverage and retained image panels.

## Run

On Apple Silicon macOS, install [micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html), then run from the repository root:

```sh
micromamba create -n lyon-2026 -f environment-osx-arm64.lock
micromamba run -n lyon-2026 python -m pip install -r requirements-pip.txt
micromamba run -n lyon-2026 python reproduce.py
```

The lock fixes native numerical and font-rendering libraries as well as Python packages. `environment.yml` provides a version-pinned specification for other platforms; those platforms have not been validated. `requirements-lock.txt` records Python package versions, but pip wheels alone use different native libraries and do not exactly reproduce the reference fits or typography. The pinned plategig source is imported directly; do not install its upstream project separately.

The default run uses the included fitted curves and redoes the statistics, data exports and rendering. To also refit all six datasets from raw OD and plate-layout workbooks:

```sh
micromamba run -n lyon-2026 python reproduce.py --refit --jobs 2
```

Refits run in a new directory under `runs/`, start with empty generated-output directories, compare all regenerated fit/cache tables with the supplied tables, execute the full pipeline there, and compare the regenerated outputs with the reference checkout. The original caches remain unchanged. An existing output directory is rejected. Use `--output /path/to/new-directory` to place a refit elsewhere. Cached runs take several minutes; refits take longer. Logs and a machine-readable result are written under the run's `runs/` directory.

## Compare outputs

Raw-refit runs automatically compare their rebuilt outputs with this checkout. To repeat the comparison separately:

```sh
micromamba run -n lyon-2026 python scripts/compare_outputs.py /path/to/reference /path/to/rebuilt --report comparison.json
```

The comparator checks all 24 CSV tables (including the panel manifest), 125 PNG files, and Source Data cell values/types. Numerical tolerance is 1e-8 relative and 1e-10 absolute. PNG bytes must match. PDF, SVG and XLSX containers include timestamps, so regeneration can change their file hashes without changing rendered content or values. The driver also checks the 36 primary contrasts against the committed numerical reference; each run records `running`, `failed` or `passed` explicitly.

## Results and source files

- [Figures](working/figures-assembled/All_figures.pdf), also individual PDF/PNG files in that directory.
- [Source Data workbook](working/source-data/Source%20Data.xlsx), with source paths, experimental units and missing-data records.
- [Final statistics](working/analysis/stats-rework/out/final_statistics.csv) and adjacent numerical CSVs.
- [Analysis code](working/analysis/stats-rework), including the pinned [plategig source and MIT license](working/analysis/stats-rework/plategig-88839a2).
- [Original workbooks and reference notebooks](working/figures), [fitting workbooks](working/analysis/stats-rework/biohpc-pull/data), and [recovered genomic tables](data/genomics).
- [Figure coverage](docs/FIGURE_COVERAGE.md), [genomic input requirements](data/genomics/README.md), and [verification record](docs/VALIDATION.md).

Paths retain their existing figure names so that notebook, workbook and manuscript references remain traceable. The repository is self-contained for the supported workflow and does not need a manuscript DOCX or access to the original workstation. Exact typography uses Times New Roman from the validation Mac’s system fonts; see the font note below.

## Analysis conventions

Primary IC50 comparisons use log10 values, culture-ID pairing for Figures 2/3, Welch tests elsewhere, and Holm adjustment within each panel and antibiotic. MIC and biological-n=1 mutant measurements are descriptive. Fit-perturbation ranges are sensitivity summaries. Fixed seeds and complete paired IDs are checked by the analysis. The historical p-value gate verifies 49 saved notebook results, including six expected NaNs; it is separate from the primary statistical analysis.

Supplementary Figure 11 absolute fractions provisionally assume equal plated volumes across time. A 10 µL baseline and 20 µL later schedule would halve those fractions. Figure 3 contains two blank parent counts; several MDK formulas use 40 µL. The ATEC-C3 exclusion reason and parent matching in Supplementary Figures 3/6 remain unconfirmed.

Genomic notebooks and recovered tables are included. Static panels retained in the ten assemblies have explicit source hashes. The TODO table and [genomic input specification](data/genomics/README.md) identify the inputs required for complete numerical regeneration.

## Fonts and reuse

Install Times New Roman to reproduce the reference typography. Matplotlib may substitute another font when it is absent, changing line wrapping and panel geometry. Fonts installed with macOS are not redistributed here. The trajectory notebooks use the bundled DejaVu Sans font explicitly where their original Nimbus Roman request fell back on the validation Mac. Numerical results do not depend on fonts.

The vendored plategig source is distributed under its included MIT license.
