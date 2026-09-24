# Reproduction guide

Run commands from the repository root. The [README](../README.md) gives the standard workflow; this page covers optional setup and verification.

## Python environment

`environment.yml` installs Python 3.11.14 and the same pinned Python packages on Windows, Linux and macOS. Micromamba runs the environment directly, so shell activation is unnecessary. Windows on ARM can use x64 Python under emulation; native ARM Windows builds are not tested.

Alternatively, use a Python 3.11 virtual environment and install `requirements-lock.txt`:

| Step | Windows PowerShell | Linux / macOS |
|---|---|---|
| Create | `py -3.11 -m venv .venv` | `python3.11 -m venv .venv` |
| Install | `.venv/Scripts/python -m pip install -r requirements-lock.txt` | `.venv/bin/python -m pip install -r requirements-lock.txt` |
| Run | `.venv/Scripts/python reproduce.py` | `.venv/bin/python reproduce.py` |

The pinned plategig source is included and imported directly. No Office application or BioRender installation is needed. Source workbooks and notebooks are preserved byte-for-byte on checkout. Keep enough disk space for an additional repository copy when using `--refit`.

## Single-cell analysis

The default Python run checks and plots the included single-cell results. Recomputing from the H5 count matrix additionally needs [R 4.4.2](https://cran.r-project.org/) on `PATH` and the packages in `single-cell-renv.lock`. Check `Rscript --version` before installing packages.

Older pinned R packages may compile from source. Install the build prerequisites for your system:

| System | Build prerequisites |
|---|---|
| Windows | [Rtools44](https://cran.r-project.org/bin/windows/Rtools/rtools44/rtools.html) for R 4.4 |
| macOS | Xcode Command Line Tools, a Fortran compiler compatible with your R installation, and HDF5; see [CRAN's macOS tools](https://mac.r-project.org/tools/) |
| Ubuntu / Debian | C/C++ and Fortran compilers, plus HDF5, curl, OpenSSL, XML, PNG and font development libraries (command below) |

```sh
sudo apt-get install build-essential gfortran libhdf5-dev libcurl4-openssl-dev libssl-dev libxml2-dev libpng-dev libfontconfig1-dev libfreetype6-dev
```

Restore the recorded packages into `runs/R-library`, then execute the analysis:

```sh
Rscript scripts/setup_single_cell.R
micromamba run -n lyon-2026 python reproduce.py --single-cell
```

The driver supplies the R library path on all platforms. `--r-library PATH` selects another library. Package restoration needs internet access; analysis uses local inputs. The R run processes 48,883 cells and computes differential expression, so allow more time and memory than for plotting the included results.

`--single-cell` replaces generated single-cell tables in the checkout. Combine it with `--refit` to run both reconstructions in an isolated copy. The R library remains in the original checkout and is reused by that copy. The [single-cell record](../data/single-cell/README.md) explains filtering, sample identities, enrichment inputs and validation tolerances.

## Outputs and checks

`--refit` writes to a new directory under `runs/`, preserves supplied caches, refits all six datasets and compares them with the included fits. `--output PATH` chooses another new directory; `--jobs 1` reduces concurrent fitting. Logs identify any failed stage, and `runs/reproduction.json` records the result. Figures are written to `working/figures-assembled/`; tables to `working/source-data/` and `working/analysis/stats-rework/out/`.

To compare two builds:

```sh
micromamba run -n lyon-2026 python scripts/compare_outputs.py PATH_TO_REFERENCE PATH_TO_BUILD --report comparison.json
```

The comparison requires matching table identities and values, workbook cells/types, PNG file sets, and the dimensions of assembled figures and print panels. Statistical values and cache fits use relative tolerance 1e-8 and absolute tolerance 1e-10. Genomic plotting tables remain exact apart from line endings. Input checksums and primary statistical tests are unchanged by the rendering environment. Full scientific checks are described in [Validation](VALIDATION.md).

PNG byte differences and font-dependent crop sizes of intermediate plots are listed without failing numerical reproduction. Add `--strict-images` to the comparator or `reproduce.py --refit` to require identical PNG bytes. Times New Roman must be installed separately to match reference typography; otherwise plots use Matplotlib's bundled STIXGeneral. `LYON_PLOT_FONT` can select an installed font. Inspect regenerated figures when changing fonts or rendering libraries.

For the recorded Apple Silicon renderer and numerical libraries, the optional `environment-osx-arm64.lock` and `requirements-pip.txt` retain the original environment:

```sh
micromamba create -y -n lyon-reference -f environment-osx-arm64.lock
micromamba run -n lyon-reference python -m pip install -r requirements-pip.txt
micromamba run -n lyon-reference python reproduce.py
```

[GitHub Actions](https://github.com/msadikyildiz/Lyon_2026/actions/workflows/reproduce.yml) rebuilds the default workflow on Windows, Linux and macOS, compares numerical outputs, and saves logs. Its manual **Run workflow** option can also refit dose-response data. The R reconstruction is validated separately; see the recorded environment and results in [Validation](VALIDATION.md).
