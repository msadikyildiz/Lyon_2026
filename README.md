# Lyon et al. (2026): data, analysis and figures

Data and analysis code for *Recurrent Extinction of Resistance Mutations Leads to Convergent Multidrug Resistance in Sequential Antibiotic Treatment*. Maintained by Adam Lyon and Muhammed Sadik Yildiz.

## TODO

| Task | Input or work remaining | Affected outputs |
|---|---|---|
| Complete sequencing provenance | Add the WGS accession and exact breseq executable version. All 149 runs now include command lines, saved options and sequence-verified references; downstream mutation tables and processed exports are included. | Reproduction from sequencing reads |
| Complete enrichment reproduction | Recover the historical EcoCyc database snapshot and analysis-specific reference sets if available. The analyst identified EcoCyc and supplied usual Fisher/BH settings and original exports; exact reruns remain distinct from reproducing the supplied tables. | Enrichment tables in Supplementary Tables 3 and 4 |
| Resolve the handling record | Establish whether pellet loss was documented for the day-2 parent sample at 3 hours. Both zero-colony records are confirmed and retained; omission sensitivity is supplied. | Figure 3e source annotation |
| Update the schematic | Replace the supplied Figure 4d BioRender panel with its final editable artwork/export. Numerical panels and phenotype labels are regenerated. | Figure 4d |
| Specify reuse terms | Add licenses for study code and data. | Repository license files |

The repository includes the experimental workbooks, mutation calls, single-cell count matrices and original analysis code. One command regenerates the numerical analyses, Source Data, Supplementary Tables 3 and 4, and all 19 figures. Schematic artwork is retained explicitly; every numerical panel is generated from included data. See the [figure map](docs/FIGURE_COVERAGE.md).

## Run

On Apple Silicon macOS, install [micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html), then run from the repository root:

```sh
micromamba create -n lyon-2026 -f environment-osx-arm64.lock
micromamba run -n lyon-2026 python -m pip install -r requirements-pip.txt
micromamba run -n lyon-2026 python reproduce.py
```

The default run uses the included dose-response fits and computed single-cell embeddings/DGE tables, checks them against supplied source tables, and rebuilds the statistics, data exports and figures. Genomic tables are recalculated from mutation calls on every run. Logs and an explicit passed/failed record are written under `runs/`.

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

This command replaces the generated single-cell tables in the checkout. Use a separate checkout to compare a different R environment. The R step executes the supplied filtered-analysis Rmd calculations with portable input paths and exports cell identities, embeddings, selected probes, normalized expression and differential-expression results. It reproduces 48,883 cells and the exact nine cluster sizes. All 6,543 cluster DGE rows and the 1,222/863 culture-versus-parent rows agree with the final supplied tables within the recorded numerical tolerances. The [single-cell record](data/single-cell/README.md) explains the sample mapping, table numbering, targeted marker checks and remaining enrichment dependency.

## Outputs and inputs

- [All figures](working/figures-assembled/All_figures.pdf), with individual PDF and PNG files alongside it. Larger genomic figures use multiple pages at readable type sizes.
- [Source Data](working/source-data/Source%20Data.xlsx), [Supplementary Table 3](working/source-data/Supplementary%20Table%203.xlsx) and [Supplementary Table 4](working/source-data/Supplementary%20Table%204.xlsx).
- [Final statistics](working/analysis/stats-rework/out/final_statistics.csv), [analysis code](working/analysis/stats-rework), [original workbooks/notebooks](working/figures), [genomic inputs](data/genomics) and [single-cell inputs](data/single-cell).
- [Source corrections](docs/DATA_CORRECTIONS.md), [figure coverage](docs/FIGURE_COVERAGE.md) and [verification](docs/VALIDATION.md).

The 6.2 GB single-cell delivery contained many duplicate project copies and figure exports. The repository retains one verified copy of the supplied count matrices, original code and final tables, approximately 196 MB. This includes alternative aggregates and excluded sample-5 matrices for provenance, with archive and per-file checksums. It does not require the original workstation or a manuscript DOCX.

## Analysis conventions

Primary IC50 comparisons use log10 values, culture-ID pairing for Figures 2/3 and Supplementary Figures 3/6, Welch tests for Supplementary Figure 4, and Holm adjustment within each panel and antibiotic. MIC and biological-n=1 mutant measurements are descriptive. Fit-perturbation ranges are sensitivity summaries. Fixed seeds and complete paired IDs are checked. The historical p-value gate verifies 49 notebook results, including six expected NaNs, separately from the primary analysis.

Confirmed zero-colony plates are distinguished from missing counts. Figure 3 parent medians and MADs are invariant over the two censored intervals; omitting the potentially mishandled 3-hour record raises the median by 4%. Supplementary Figure 11 uses the confirmed single 10 µL schedule, so equal volume cancels in normalized fractions and detection limits. The recorded 40 µL exceptions and ATEC repeat/exclusion reasons are documented in the source-correction record.

Genomic plots use coordinate/allele identity rather than averaging different variants under the same short label. Each plotted entry has a source crosswalk. Single-cell tests use cells as observations; the six libraries comprise two technical replicates from each of three cultures. EcoCyc enrichment exports and matched gene lists are retained and verified against the supplied workbook. The complete parent-up export extends its 19-row worksheet to 133 rows without changing existing results. The historical database and reference-set limitations are described in the [enrichment record](data/single-cell/enrichment/README.md). GSE314756 is scheduled for public release on 21 December 2026.

## Compare outputs

```sh
micromamba run -n lyon-2026 python scripts/compare_outputs.py /path/to/reference /path/to/rebuilt --report comparison.json
```

The comparator checks CSV values, PNG bytes, generated genomic/single-cell plotting tables, and Source Data cell values/types. Statistical/workbook tolerances are 1e-8 relative and 1e-10 absolute. PDF, SVG and XLSX containers may differ in timestamps. R reconstruction is validated separately against the original expression tables, including exact row identities and cluster counts.

The vendored plategig source retains its [MIT license](working/analysis/stats-rework/plategig-88839a2/LICENSE).
