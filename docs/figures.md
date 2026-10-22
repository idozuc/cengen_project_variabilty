# Manuscript figures

`scripts/figures/` regenerates every data figure of the manuscript draft
(Figures 1–15) from the outputs of the analysis steps. It only draws figures:
BASiCS, enrichment and clustering run in the RUNBOOK steps listed under
[Inputs](#inputs). The multi-panel figures (3, 4, 5, 6, 7) were assembled from
these panels by hand; file names carry the figure number and panel letter.

```bash
cp config/figures.example.yaml config/figures.yaml   # edit the paths
scripts/figures/run_all.sh config/figures.yaml       # about a minute
```

`run_all.sh` runs the four `0*_export_*.R` scripts (they write
`OUTPUT_DIR/data/`), then every `fig*.py` script (they write
`OUTPUT_DIR/panels/`), then compiles the two TikZ workflow diagrams. Run single
scripts from the repository root with `--config config/figures.yaml`, exports
first.

| Config key | Meaning |
|---|---|
| `sce` | CeNGEN L4 neuron `SingleCellExperiment` RDS |
| `basics_dir` | BASiCS run directory (`outputs/basics/RUN_ID`) |
| `clustering_dir` | Clustering run directory (`outputs/clustering/RUN_ID`) |
| `positive_controls` | RDS from `scripts/clustering/run_positive_controls.R` |
| `output_dir` | Figure build directory (`outputs/figures/RUN_ID`) |

## Figures

| Fig. | Content | Script | Panel files |
|---|---|---|---|
| 1 | Workflow: HVGs and enrichment | `workflows/methods_hvg_workflow.tex` | `fig01_methods_hvg_workflow.pdf` |
| 2 | Workflow: within-cell-type clustering | `workflows/methods_clustering_workflow.tex` | `fig02_methods_clustering_workflow.pdf` |
| 3 A–F | Histograms of ASE/AWC markers (*gcy-5, gcy-7, gcy-22, gcy-6, srsx-3, str-2*) | `fig03_05_expression_histograms.py` | `fig03{A–F}_*.svg` |
| 4 A–B | Expected-marker evidence (A: AWC, B: ASE) | `fig04_marker_evidence_bars.py` | `fig04{A,B}_*.svg` |
| 5 A–B | HVG/LVG examples in URX; high-mean, low-zero HVGs in AVK and RIG | `fig03_05_expression_histograms.py` | `fig05{A,B}_*.svg` |
| 6 A–C | HVG calls by zero-expression bin; Mu vs Epsilon | `fig06_hvg_variability.py` | `fig06{A,B,C}_*.svg` |
| 7 A–B | Neuropeptide enrichment in RMD_DV; cell types per WormCat category | `fig07_hvg_enrichment.py` | `fig07{A,B}_*.svg` |
| 8 | Recurrent WormCat enrichment by functional group | `fig08_enrichment_heatmap.py` | `fig08_*.svg` |
| 9 | PCA/UMAP of the ASE and AWC positive controls | `fig09_positive_controls_pca_umap.py` | `fig09_*.svg` |
| 10–11 | Marker expression on the ASE and AWC UMAPs | `fig10_11_marker_expression_umaps.py` | `fig10_*.svg`, `fig11_*.svg` |
| 12 A–C | Marker detection per cluster in HSN, SMD, NSM | `fig12_subgroup_markers.py` | `fig12_subgroup_markers.pdf` |
| 13–15 | HSN/*tph-1*, SMD/*flp-5*, NSM/*del-3* on the clustering UMAP | `fig13_15_gene_umaps.py` | `fig13_*.svg` … `fig15_*.svg` |

## Inputs

| Input | Made by (RUNBOOK step) | Used by |
|---|---|---|
| `sce` | CeNGEN (step 1) | 01, 03 exports |
| `BASICS_DIR/basics_master_with_detection.csv.gz` | `scripts/analysis/04_add_detection_columns.R` (13) | 01 export, fig06 |
| `BASICS_DIR/analysis/neuron_functional_groups.csv` | `scripts/analysis/01_prepare_mappings.R` (10) | fig08 |
| `BASICS_DIR/analysis/wormcat_hvg_cengen_t2.csv` | `scripts/analysis/05_wormcat_cengen_filtered.py` (17) | fig07, fig08 |
| `positive_controls` | `scripts/clustering/run_positive_controls.R` (16) | 01, 02 exports |
| `CLUSTERING_DIR/run_summary.csv`, `{HSN,SMD,NSM}/assignments_by_resolution.csv` | `scripts/clustering/run_clustering.R` (2) | 03 export |
| `CLUSTERING_DIR/{HSN,SMD,NSM}/visualization_bundle.rds` | same | 04 export |

The clustering and positive-control runs must use the default transgene
exclusion, as the manuscript does.

## Exports

| Script | Writes to `OUTPUT_DIR/data/` |
|---|---|
| `01_export_expression_data.R` | Raw and log1p counts of 11 genes per cell, with the BASiCS Mu, Epsilon and % zeros per panel. URX, AVK and RIG cells come from experiments with ≥ 10 cells (the BASiCS inclusion rule); ASE and AWC cells are those retained by the positive controls |
| `02_export_positive_control_pcs.R` | Selected PC scores, GMM clusters and known identities of the two positive controls, and their statistics |
| `03_export_subgroup_markers.R` | Fraction of cells with ≥ 1 UMI per gene and cluster in HSN, SMD, NSM |
| `04_export_gene_umap_data.R` | Stored UMAP, cluster and counts of *tph-1*, *flp-5*, *del-3* |

Two choices are made in `02_export_positive_control_pcs.R`, not in the
positive-control run:

- For ASE it uses the G = 3 solution. The pipeline rejected it only as
  qc-sensitive within the borderline margin, and the borderline rule accepts it
  (Methods).
- It fixes the sign of each PC (long tail towards negative values). irlba
  starts from a random vector, so PC signs are otherwise arbitrary between
  runs; clusters and UMAPs do not depend on them.

## Definitions fixed in the plotting scripts

- **Analysis HVG** (Figure 6): BASiCS HVG call, `analysis_eligible`, and
  detected in more than half of the cell type's experiments.
- **Marker evidence** (Figure 4): a marker counts as expressed at ≥ 1 UMI.
- **Enrichment** (Figures 7B, 8): WormCat category 3, FDR ≤ 0.05, fold > 1;
  shown when enriched in ≥ 3 cell types.
- **Cluster outlines** (Figures 10–11): convex hull of the 95% of each
  cluster's cells closest to its median.
- **UMAP** (Figures 9–11): computed from the selected PCs, 15 neighbours,
  `min_dist` 0.3, seed 20260809, single-threaded. Figures 13–15 use the UMAP
  stored by the clustering run.

## Software and reproducibility

The panels were built with R 4.6.1 (SingleCellExperiment 1.34.0, mclust 6.1.3,
irlba 2.3.7, Matrix 1.7.5) and Python 3.11.7 with the `figures` extra of
`pyproject.toml`. Built from the original analysis outputs, every panel
reproduces the manuscript's text and embedded images exactly, except the
Figure 8 colour bar, which is one pixel shorter (the original was drawn with
matplotlib 3.8.4). Other umap-learn or numba versions can rotate or mirror the
UMAP layouts.

Built from a fresh clustering run, the clusters and every number in every
figure are the same, but the stored UMAP layout of Figures 13–15 can differ.
The residual PCA uses irlba with an unseeded random start, so PC scores vary
between runs at about 1e-7, and UMAP can turn that into a different layout. In
the check run this happened for NSM (Figure 15), not for HSN or SMD.

Text is Arial, kept editable in the SVGs. Without Arial, matplotlib substitutes
another font and warns.

The title-page images (university logo, NeuroPAL image from the Hobert Lab)
are not produced by these scripts.
