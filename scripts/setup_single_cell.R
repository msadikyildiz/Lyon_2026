# Restore the recorded R packages to a project-local library on every platform.
if (getRversion() != '4.4.2') {
  stop('Use R 4.4.2 for the recorded analysis environment; see docs/REPRODUCING.md.')
}
script <- sub('^--file=', '', commandArgs()[grepl('^--file=', commandArgs())][1])
root <- dirname(dirname(normalizePath(script, winslash='/')))
lib <- file.path(root, 'runs/R-library')
dir.create(lib, recursive=TRUE, showWarnings=FALSE)
.libPaths(c(lib, .libPaths()))
if (!requireNamespace('renv', quietly=TRUE)) {
  install.packages('renv', repos='https://cloud.r-project.org', lib=lib)
}
renv::restore(lockfile=file.path(root, 'single-cell-renv.lock'), library=lib, prompt=FALSE)
cat('R library ready. Run: python reproduce.py --single-cell\n')
