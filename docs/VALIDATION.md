# Verification

Validation uses Python 3.11.14 on Apple Silicon macOS with `environment-osx-arm64.lock` plus `requirements-pip.txt`. The environment was installed independently of the analysis workstation environment. Times New Roman is installed separately.

The cached workflow passes six failure/invariance tests, the 49-value historical check and 357 revised-analysis checks. It regenerates 23 numerical tables (excluding the path/provenance panel manifest), Source Data, 59 analytical PNGs, 40 size-specific print panels and ten figure assemblies. The figures' internal checks validate panel identity, aspect ratios, placement, source hashes and content.

A separate refit creates new caches from all six sets of raw OD/plate-layout workbooks and then runs the same analysis and figure workflow. The cache comparator checks every column, including all nested fit/bootstrap arrays, at relative tolerance 1e-8 and absolute tolerance 1e-10, with matching NaN locations and exact nonnumeric values. Refit outputs are isolated from the supplied caches.

The pip-only environment was also tested. Different native FreeType and numerical-library builds changed glyph rendering and some optimizer results, despite identical Python package versions. The documented default therefore uses the native-library lock. Cross-platform builds from `environment.yml` have not been validated.

Detailed completed-run measurements are recorded in `validation.json`. Runtime logs are kept under `runs/` and are not part of the distribution. Figure coverage and unresolved experimental records are stated in the README and `FIGURE_COVERAGE.md`.
