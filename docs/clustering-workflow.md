# Clustering workflow

`scripts/clustering/run_clustering.R` selects a cell type, removes sparsely
represented experiments, normalizes by size factor, estimates an
experiment-adjusted stress score, and residualizes candidate features against
experiment, stress, and library-size nuisances. It selects high-variance
residual features, runs PCA, screens Gaussian mixtures, and performs stability
and nuisance diagnostics only for structurally supported candidates.

Accepted resolutions receive geometric-driver, marker, cluster-summary, and
UMAP outputs. UMAP is visualization-only. Experiment, detection, feature count,
stress, and size-factor diagnostics can veto a resolution according to the
preserved thresholds; mitochondrial percentage remains diagnostic-only.

Use `--screen-only` for the cheap phase. Existing terminal results are reused
unless `--overwrite` is supplied. The exporter accepts a complete run directory
and writes only accepted results and supported-marker expression.
