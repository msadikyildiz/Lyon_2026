# Lyon et al. (2026): data, analysis and figures

Code and data for *Recurrent Extinction of Resistance Mutations Leads to Convergent Multidrug Resistance in Sequential Antibiotic Treatment*.

## TODO

| Item | Status |
|---|---|
| Public sequencing reads | BioProject [PRJNA1534522](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1534522) is assigned. SRA reads are processing; immediate release was requested. |

## Reproduce the study

Download or clone this repository, install [micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html), and open a terminal in the repository folder. The commands are the same on **Windows (PowerShell), Linux and macOS**:

```sh
micromamba create -y -n lyon-2026 -f environment.yml
micromamba run -n lyon-2026 python reproduce.py
```

This rebuilds all **19 figures**, statistical results, Source Data and Supplementary Tables 3–4. It uses the included dose-response fits and single-cell results, while recalculating genomic tables from mutation calls. Every run checks input hashes and numerical results; logs and pass/fail status are saved in `runs/`.

| Analysis | Included inputs → outputs |
|---|---|
| Resistance and survival | Plate-reader and colony-count workbooks → IC50/MIC, survival, MDK and statistical comparisons |
| Genomics | Mutation calls, reference sequences and run records → mutation tables, trajectories and heatmaps |
| Single-cell expression | Count matrices and R/Seurat code → clusters, differential expression and expression plots |

To refit all six dose-response datasets from raw measurements in an isolated copy:

```sh
micromamba run -n lyon-2026 python reproduce.py --refit
```

Exact raw-refit validation uses the recorded Apple Silicon environment. Windows/Linux refits can fail strict numerical comparisons; the default workflow above passes on all three systems ([details](docs/REPRODUCING.md#raw-refits)).

To also recompute single-cell results from counts, install **R 4.4.2** and its build prerequisites ([setup guide](docs/REPRODUCING.md#single-cell-analysis)), then run:

```sh
Rscript scripts/setup_single_cell.R
micromamba run -n lyon-2026 python reproduce.py --refit --single-cell
```

Supplied EcoCyc enrichment results are validated and exported; enrichment tests and raw-read alignment are not rerun. Schematics use included artwork. Figures use Times New Roman when installed, otherwise the bundled STIX serif font. Numerical checks remain strict across platforms; font rendering can differ.

## Find the results

- [All figures (PDF)](working/figures-assembled/All_figures.pdf), with individual PDF/PNG files in the same folder.
- [Source Data](working/source-data/Source%20Data.xlsx), [Supplementary Table 3](working/source-data/Supplementary%20Table%203.xlsx), [Supplementary Table 4](working/source-data/Supplementary%20Table%204.xlsx) and [statistics](working/analysis/stats-rework/out/final_statistics.csv).
- [Figure-by-figure map](docs/FIGURE_COVERAGE.md), [analysis conventions](docs/DATA_CORRECTIONS.md), [validation](docs/VALIDATION.md) and [setup/comparison options](docs/REPRODUCING.md).

Sequencing deposits: WGS [PRJNA1534522](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1534522); single-cell [GSE314756](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE314756), scheduled for release on 21 December 2026. The count matrices and processed genomic inputs needed here are included in the repository.

## Citation and reuse

Cite the study and [repository](CITATION.cff), including the Git commit used. Code and documentation: [MIT](LICENSE). Study-generated data, derived tables and original figure elements: [CC BY 4.0](LICENSE-DATA). Third-party software, reference sequences, database annotations and BioRender assets retain their own terms.
