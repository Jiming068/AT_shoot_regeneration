# Simplified manuscript analysis scripts

This directory contains a sanitized, compact copy of the analysis scripts used
for the Arabidopsis shoot-regeneration manuscript.

## Scope of the cleanup

- Preserved the figure-specific analysis logic, model parameters, thresholds,
  random seeds, plotting settings, and output calls.
- Removed local usernames, absolute macOS/Windows/HPC paths, private library
  paths, environment-history comments, and interactive diagnostic statements.
- Replaced machine-specific paths with portable environment-variable settings.
- Removed notebook outputs, embedded figures, execution timestamps, execution
  counters, and unused inspection-only cells.
- Normalized text scripts to UTF-8 with Unix line endings.

No expression matrices, RDS objects, credentials, unpublished numerical output,
or other source data are included in this directory.

## Portable path configuration

Scripts that previously contained machine-specific paths now use the following
environment variables:

- `SHOOT_REGEN_PROJECT_DIR`: project root; defaults to the current directory.
- `SHOOT_REGEN_OUTPUT_DIR`: output root; defaults to `02_results` below the
  project root.
- `SHOOT_REGEN_EXTERNAL_DATA_DIR`: location for inputs originally stored outside
  the project.
- `SHOOT_REGEN_WORK_DIR`: working directory for legacy scripts that relied on
  `setwd()`; defaults to the current directory.

Example:

```bash
export SHOOT_REGEN_PROJECT_DIR=/path/to/01_AT_shoot_regeneration
export SHOOT_REGEN_OUTPUT_DIR=/path/to/results
export SHOOT_REGEN_EXTERNAL_DATA_DIR=/path/to/external_data
export SHOOT_REGEN_WORK_DIR=/path/to/current_figure_workdir
Rscript Figure_2E_2F_monocle3.R
```

Scripts that already provided command-line arguments retain those interfaces.

## Main software dependencies

R packages used across the scripts include Seurat, Matrix, monocle3, ggplot2,
tidyverse, clusterProfiler, org.At.tair.db, ComplexHeatmap, pheatmap, harmony,
igraph, ggraph, patchwork, data.table, and related plotting utilities.

Python packages used across the scripts and notebooks include numpy, pandas,
scanpy, scvelo, scipy, scikit-learn, matplotlib, OpenCV, and stereo.

Exact packages used by an individual analysis remain visible in that script.

## Validation

- All standalone R scripts were checked with `parse()`.
- All standalone Python scripts were checked with `py_compile`.
- All notebooks were validated as JSON and contain zero stored outputs and zero
  execution counters.
- `SANITIZATION_AUDIT.json` records the file-level size reduction and notebook
  cleanup statistics.

## Notebook conversions

For command-line execution and code deposition, each notebook also has a plain
script counterpart. Notebook cell order is preserved and Markdown cells are
retained as comments:

- `Figure_3B_cluster.ipynb` -> `Figure_3B_cluster.py`
- `Figure_S7B_in_situ_expression.ipynb` -> `Figure_S7B_in_situ_expression.py`
- `Figure_S9B_cluster.ipynb` -> `Figure_S9B_cluster.R`
- `Figure_S9D_rds.ipynb` -> `Figure_S9D_rds.R`

The original sanitized notebooks remain in this directory for provenance.
`NOTEBOOK_CONVERSION_MANIFEST.tsv` records language, cell counts, line counts,
and MD5 checksums for each conversion.

The sanitized scripts should still be reviewed together with the associated
Methods and source-data manifest before public deposition. Large input datasets
must be obtained from the manuscript data repository and placed at the paths
specified through the environment variables above.
