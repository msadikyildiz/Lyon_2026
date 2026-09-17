# Lyon et al. (2026): data, analysis and figures

Reproduction repository for the revised Lyon manuscript. The available dose-response, growth and survival analyses rebuild from the included workbooks. One command regenerates the statistical tables, Source Data workbook, 40 size-specific panels, six original trajectory plots and ten assembled figures. Genomic and single-cell coverage is incomplete; the [figure map](docs/FIGURE_COVERAGE.md) identifies the remaining inputs and retained image panels.

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

Refits run in a new directory under `runs/`, compare all regenerated fit/cache tables with the supplied tables, and execute the full pipeline there. The original caches remain unchanged. An existing output directory is rejected. Use `--output /path/to/new-directory` to place a refit elsewhere. Cached runs take several minutes; refits take longer. Logs and a machine-readable result are written under the run's `runs/` directory.

## Results and source files

- [All revised figures](working/figures-assembled/All_revised_figures.pdf), also individual PDF/PNG files in that directory.
- [Source Data workbook](working/source-data/Source%20Data.xlsx), with source paths, experimental units and missing-data records.
- [Final statistics](working/analysis/stats-rework/out/final_statistics.csv) and adjacent numerical CSVs.
- [Analysis code](working/analysis/stats-rework), including the pinned [plategig source and MIT license](working/analysis/stats-rework/plategig-88839a2).
- [Original workbooks and reference notebooks](working/figures), [fitting workbooks](working/analysis/stats-rework/biohpc-pull/data), and [recovered genomic tables](data/genomics).
- [Figure coverage](docs/FIGURE_COVERAGE.md), [genomic input requirements](data/genomics/README.md), and [verification record](docs/VALIDATION.md).

Paths retain their existing figure names so that notebook, workbook and manuscript references remain traceable. The repository is self-contained for the supported workflow and does not need a manuscript DOCX or access to the original workstation. Exact typography uses Times New Roman, supplied by Microsoft Office on the validation Mac; see the font note below.

## Analysis conventions

Primary IC50 comparisons use log10 values, culture-ID pairing for Figures 2/3, Welch tests elsewhere, and Holm adjustment within each panel and antibiotic. MIC and biological-n=1 mutant measurements are descriptive. Fit-perturbation ranges are sensitivity summaries. Fixed seeds and complete paired IDs are checked by the analysis. The historical p-value gate verifies 49 saved notebook results, including six expected NaNs; it is separate from the revised statistical analysis.

Supplementary Figure 11 absolute fractions provisionally assume equal plated volumes across time. A 10 µL baseline and 20 µL later schedule would halve those fractions. Figure 3 contains two blank parent counts; several MDK formulas use 40 µL. The ATEC-C3 exclusion reason and parent matching in Supplementary Figures 3/6 remain unconfirmed. These records and their consequences are retained for author review.

Genomic notebooks and recovered tables are included, but per-sample mutation calls and some full frequency tables are still missing. Single-cell matrices and analysis code are also needed. These gaps prevent a claim that every manuscript figure can be regenerated from raw data. Static panels retained in the ten assemblies have explicit source hashes; figure assembly alone does not regenerate those panels.

## Fonts and reuse

Install Times New Roman to reproduce the reviewed typography. Matplotlib may substitute another font when it is absent, changing line wrapping and panel geometry. Microsoft fonts are not redistributed here. Numerical results do not depend on fonts.

This repository is private for coauthor review. The authors' code/data license and public release remain to be decided. The vendored plategig MIT license applies to that component. No archival deposition is associated with this version.
