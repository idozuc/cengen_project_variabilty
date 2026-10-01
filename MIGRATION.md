# Migration inventory

The source workspace `/Users/idozuckerman/Documents/Scrutiny CeNGen` was read
but not modified. Files were reorganized as follows.

| Original | New location |
|---|---|
| `Clustering/Clustering_pipeline.R` | `src/R/clustering.R` |
| `Clustering/Clustering_visualization.R` | `src/R/clustering_visualization.R` |
| `R/stress_gene_sets.R` | `src/R/stress_gene_sets.R` |
| `Clustering/run_clustering_pipeline.R` | `scripts/clustering/run_clustering.R` |
| `Clustering/export_streamlit_clustering.R` | `scripts/clustering/export_explorer_assets.R` |
| `HVG/basics/prepare.R` | `scripts/basics/01_prepare.R` |
| `HVG/basics/basics.R` | `scripts/basics/02_fit_cell_type.R` |
| `HVG/basics/diagnose.R` | `scripts/basics/03_diagnose.R` |
| `HVG/basics/detect_lvg.R` | `scripts/basics/04_detect_lvg.R` |
| `HVG/basics/assess_experiment_support.R` | `scripts/basics/05_assess_experiment_support.R` |
| `HVG/Analysis/01_build_master.R` | `scripts/basics/06_build_master.R` |
| `HVG/Analysis/02_prepare_mappings.R` | `scripts/analysis/01_prepare_mappings.R` |
| `HVG/Analysis/03_run_enrichment.R` | `scripts/analysis/02_run_enrichment.R` |
| `HVG/Analysis/04_compare_detection_enrichment.R` | `scripts/analysis/03_compare_detection_enrichment.R` |
| `HVG/Analysis/05_add_detection_columns.R` | `scripts/analysis/04_add_detection_columns.R` |
| `HVG/basics/run_basics.sh` | `slurm/run_basics_array.sh` |
| `HVG/basics/run_lvg.sh` | `slurm/run_lvg_array.sh` |
| `tmp/entropy_streamlit/app.py` | `explorer/app.py` |

## Deliberately excluded

- `Clustering/Unimportant`, experimental GLM-PCA comparisons, and historical reports.
- The alternative scran workflow under `HVG/R` and `HVG/run_scran_hvg.R`.
- Figure and manuscript generation.
- `Rlib*`, Python caches, histories, temporary files, and local environments.
- Raw SCE files, BASiCS chains, completed runs, generated figures, and exports.
- The separate `entropy` project; the explorer now has an in-repository asset builder.

Scientific algorithms and thresholds were retained. Changes are limited to
layout, configuration, portability, output contracts, validation, and provenance.
