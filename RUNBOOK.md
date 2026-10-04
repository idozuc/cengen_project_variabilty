# CeNGEN project runbook

This is the single start-to-finish operating guide for the supported
clustering, BASiCS/HVG, analysis, and Streamlit workflows. Run every command
from the repository root:

```bash
cd /path/to/cengen_project
```

On the current workstation that path is
`/Users/idozuckerman/Documents/cengen_project`. Commands below remain portable
because the scripts resolve project files from the repository root.

Use a unique run ID for every production analysis. Examples below use
`RUN_ID=example`; replace it with a date or release name.

## 0. Install the environments

### R

**Input:** `renv.lock` and an R 4.4 installation.

```bash
R -q -e 'install.packages("renv"); renv::restore()'
```

**Output:** a project-local `renv/library/` containing BASiCS,
SingleCellExperiment, Matrix, arrow, mclust, irlba, coda, uwot, Shiny, Plotly,
readxl, yaml, and their dependencies. The library is ignored by Git.

### Python

**Input:** `pyproject.toml` and Python 3.11 or newer.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e '.[test]'
```

**Output:** an ignored `.venv/` with pandas, PyArrow, Streamlit, Plotly,
PyYAML, pytest, and schema-validation packages.

## 1. Supply the shared input data

Both main pipelines start from the same RDS containing a
`SingleCellExperiment`.

**Required assay**

- `counts`: genes × cells raw UMI counts.
- Count row names must be unique stable WormBase IDs.
- Count column names must be unique cell IDs.

**Required `colData`**

| Field | Clustering | BASiCS | Meaning |
|---|---:|---:|---|
| `Cell.type` | yes | yes | CeNGEN cell-type label |
| `Experiment` | yes | yes | experiment/batch label |
| `Size_Factor` | yes | no | positive normalization factor |
| `Detection` | yes | no | per-cell detection metric |
| `pct_counts_Mito` | yes | no | mitochondrial percentage |
| `total_features_by_counts` | yes | no | detected feature count |

`rowData(sce)$gene_short_name` is strongly recommended and is required for
readable BASiCS and explorer labels.

Place the RDS under the ignored `data/raw/` directory or reference a shared
absolute path from a local YAML file. Do not commit the RDS.

## 2. Run nuisance-adjusted clustering

Entry point: `scripts/clustering/run_clustering.R`.

### Input

- The SCE described above.
- The embedded 199-gene stress set from `src/R/stress_gene_sets.R`.
- Optional `config/clustering.yaml`, copied from the example.

```bash
cp config/clustering.example.yaml config/clustering.yaml
# Edit sce_path and run_id before continuing.
Rscript scripts/clustering/run_clustering.R --config config/clustering.yaml
```

For a cheap one-type check first:

```bash
Rscript scripts/clustering/run_clustering.R \
  --sce data/raw/L4_neuron_sce.rds \
  --cell-type SMD \
  --run-id smoke \
  --screen-only
```

### Processing

The runner filters rare experiments, normalizes counts, estimates stress,
residualizes experiment/stress/library-size nuisances, selects variable
features and PCs, screens Gaussian mixtures, and then runs stability and
nuisance diagnostics for multi-cluster candidates. Accepted candidates receive
marker, driver, feature-summary, and UMAP outputs. UMAP never determines the
clusters.

### Output

By default, `outputs/clustering/RUN_ID/` contains:

```text
resolved_config.rds
resolved_config.yaml
session_info.txt
run_summary.csv
CELL_TYPE/
├── status.rds
├── screen_result.rds
├── gmm_screen.csv
├── result.rds                       # candidates analyzed in depth
├── resolution_diagnostics.csv       # candidates analyzed in depth
├── assignments_by_resolution.csv    # candidates analyzed in depth
├── primary_drivers.csv              # accepted multi-cluster only
├── primary_markers.csv              # accepted multi-cluster only
├── primary_feature_summaries.csv    # accepted multi-cluster only
└── visualization_bundle.rds         # accepted multi-cluster only
```

The runner reuses terminal results. Supply `--overwrite` only when deliberately
recomputing them.

## 3. Export clustering assets for Streamlit

Entry point: `scripts/clustering/export_explorer_assets.R`.

### Input

- A completed clustering run directory.
- Its `run_summary.csv` and accepted cell-type result files.

```bash
Rscript scripts/clustering/export_explorer_assets.R \
  --input-dir outputs/clustering/RUN_ID \
  --output-dir outputs/explorer/clustering
```

### Output

```text
outputs/explorer/clustering/
├── manifest.parquet
└── CELL_TYPE/
    ├── cells.parquet
    ├── expression.parquet
    ├── summary.parquet
    ├── qc_by_resolution.parquet
    └── markers.parquet
```

Only accepted multi-cluster cell types are exported. `expression.parquet`
contains supported markers only, never the complete count matrix.

## 4. Prepare a BASiCS run

Entry point: `scripts/basics/01_prepare.R`.

### Input

- The SCE described in step 1.
- A new or empty run directory.
- BASiCS thresholds and MCMC settings, normally supplied through YAML.

```bash
cp config/basics.example.yaml config/basics.yaml
# Edit sce and run_dir. The production defaults are already present.
Rscript scripts/basics/01_prepare.R --config config/basics.yaml
```

### Output

```text
outputs/basics/RUN_ID/
├── config.rds
├── resolved_config.yaml
├── manifest.txt
├── session_info.txt
├── cell_types.tsv
├── experiment_counts.tsv
├── logs/
└── cell_types/
```

`cell_types.tsv` freezes task IDs, cell types, seeds, evidence tiers, and cell
counts. `config.rds` records the input MD5; later stages refuse a changed SCE.

## 5. Fit BASiCS for every cell type

Entry point: `scripts/basics/02_fit_cell_type.R`.

### Input

- The prepared BASiCS run.
- One `task_id` from `cell_types.tsv`.
- The unchanged SCE registered during preparation.

For a local one-task smoke run:

```bash
Rscript scripts/basics/02_fit_cell_type.R \
  --run-dir outputs/basics/RUN_ID \
  --task-id 1
```

For all tasks on SLURM:

```bash
RUN=outputs/basics/RUN_ID
N=$(( $(wc -l < "$RUN/cell_types.tsv") - 1 ))
sbatch --array="1-$N" \
  --output="$RUN/logs/%A_%a.out" \
  --error="$RUN/logs/%A_%a.err" \
  --export="ALL,BASICS_RUN=$PWD/$RUN" \
  slurm/run_basics_array.sh
```

### Output per task

```text
cell_types/NNN_CELL_TYPE/
├── chain.rds
├── genes.csv
├── experiments.csv
├── fit_info.tsv
└── _COMPLETE
```

`genes.csv` contains stable ID, symbol, Mu, Delta, Epsilon, posterior
probability, HVG call, detection counts, seed, cells, experiments, and evidence
tier. `_COMPLETE` is written last; reruns skip completed tasks.

## 6. Diagnose MCMC convergence

Entry point: `scripts/basics/03_diagnose.R`.

### Input

- The prepared run and completed `chain.rds` files.

```bash
Rscript scripts/basics/03_diagnose.R --run-dir outputs/basics/RUN_ID
```

### Output

- `cell_types/NNN_CELL_TYPE/diagnostics.csv`: epsilon ESS and Geweke values for
  every modeled gene.
- `diagnostic_summary.csv`: completion state, median ESS, and fraction of genes
  below ESS 100 for every task.

Review this summary before continuing. Extend a run only for poor convergence,
using a new run directory and consistently longer MCMC settings.

## 7. Call lowly variable genes

Entry point: `scripts/basics/04_detect_lvg.R`.

### Input

- `chain.rds`, `genes.csv`, and `diagnostics.csv` for one completed task.
- The prepared run configuration.

Local example:

```bash
Rscript scripts/basics/04_detect_lvg.R \
  --run-dir outputs/basics/RUN_ID \
  --task-id 1
```

SLURM array:

```bash
RUN=outputs/basics/RUN_ID
N=$(( $(wc -l < "$RUN/cell_types.tsv") - 1 ))
sbatch --array="1-$N" \
  --output="$RUN/logs/lvg_%A_%a.out" \
  --error="$RUN/logs/lvg_%A_%a.err" \
  --export="ALL,BASICS_RUN=$PWD/$RUN" \
  slurm/run_lvg_array.sh
```

### Output per task

- `lvg.csv`: LVG call, probability, gene metadata, and convergence diagnostics.
- `_LVG_COMPLETE`: written only after `lvg.csv` is finalized.

## 8. Assess cross-experiment HVG support

Entry point: `scripts/basics/05_assess_experiment_support.R`.

### Input

- The exact SCE registered in the BASiCS run.
- Completed `genes.csv` and `experiments.csv` files.

```bash
Rscript scripts/basics/05_assess_experiment_support.R \
  --sce data/raw/L4_neuron_sce.rds \
  --run-dir outputs/basics/RUN_ID
```

### Output

- `hvg_experiment_detection.csv`: one row per HVG and eligible experiment.
- `hvg_experiment_support.csv`: number of supporting experiments and a
  `robust_detection_support` flag per cell-type/gene pair.

Support means detection in at least `max(3 cells, 5% of experiment cells)`;
robust support requires at least two experiments. Original HVG calls are not
changed.

## 9. Build the canonical BASiCS master table

Entry point: `scripts/basics/06_build_master.R`.

### Input

- `cell_types.tsv`.
- Every task’s `genes.csv`, `lvg.csv`, and `diagnostics.csv`.
- The embedded stable-ID stress list.

```bash
Rscript scripts/basics/06_build_master.R \
  --run-dir outputs/basics/RUN_ID \
  --output outputs/basics/RUN_ID/basics_master.csv.gz
```

### Output

`basics_master.csv.gz`, with one unique row per cell type and stable gene ID.
It joins HVG/LVG estimates, detection data, convergence metrics, stress flags,
and `analysis_eligible`. A row is eligible only when epsilon is valid, ESS is at
least 100, and the gene is not in the 199-gene stress list.

## 10. Prepare gene and neuron mappings

Entry point: `scripts/analysis/01_prepare_mappings.R`.

### Input

- The canonical master table.
- WormCat whole-genome CSV.
- CeNGEN Supplement 13 Excel workbook.

```bash
Rscript scripts/analysis/01_prepare_mappings.R \
  --master outputs/basics/RUN_ID/basics_master.csv.gz \
  --wormcat /path/to/wormcat_whole_genome.csv \
  --cengen-supplement /path/to/cengen_supplement_13.xlsx \
  --output-dir outputs/basics/RUN_ID/analysis
```

### Output

- `gene_classes.csv`: stable IDs mapped to the three WormCat category levels.
- `neuron_classes.csv`: analyzed cell types mapped to neuron families.
- `neuron_functional_groups.csv`: complete published neuron mapping.

## 11. Run standard HVG/LVG enrichment

Entry point: `scripts/analysis/02_run_enrichment.R`.

### Input

- Master table, `gene_classes.csv`, and `neuron_classes.csv`.
- `--target HVG` or `--target LVG`.

```bash
for TARGET in HVG LVG; do
  Rscript scripts/analysis/02_run_enrichment.R \
    --master outputs/basics/RUN_ID/basics_master.csv.gz \
    --gene-classes outputs/basics/RUN_ID/analysis/gene_classes.csv \
    --neuron-classes outputs/basics/RUN_ID/analysis/neuron_classes.csv \
    --output-dir outputs/basics/RUN_ID/analysis \
    --target "$TARGET"
done
```

### Output per target

- `wormcat_TARGET_enrichment.csv`: category enrichment within cell types.
- `neuron_family_TARGET_enrichment.csv`: gene calls concentrated in neuron families.
- `neuron_gene_family_TARGET_enrichment.csv`: category enrichment across each family.

Names use lowercase `hvg` or `lvg` in actual filenames.

## 12. Compare full versus experiment-supported enrichment

Entry point: `scripts/analysis/03_compare_detection_enrichment.R`.

### Input

- Master table and original SCE.
- Prepared BASiCS run.
- `gene_classes.csv`.

```bash
Rscript scripts/analysis/03_compare_detection_enrichment.R \
  --master outputs/basics/RUN_ID/basics_master.csv.gz \
  --sce data/raw/L4_neuron_sce.rds \
  --run-dir outputs/basics/RUN_ID \
  --gene-classes outputs/basics/RUN_ID/analysis/gene_classes.csv \
  --output-dir outputs/basics/RUN_ID/analysis/support
```

### Output

- `gene_experiment_support.csv.gz`.
- `wormcat_full.csv.gz`, `wormcat_supported.csv.gz`, and
  `wormcat_comparison.csv.gz`.
- `wormcat_celltype_summary.csv` and
  `manifest.txt` with annotation/input checksums and testing rules.

`wormcat_supported.csv.gz` is the WormCat input used by the explorer.

## 13. Add plotting-friendly detection columns (optional)

Entry point: `scripts/analysis/04_add_detection_columns.R`.

### Input

- Canonical master table.

```bash
Rscript scripts/analysis/04_add_detection_columns.R \
  --input outputs/basics/RUN_ID/basics_master.csv.gz \
  --output outputs/basics/RUN_ID/basics_master_with_detection.csv.gz
```

### Output

The same master rows plus `pct_zero` and a readable
`experiment_detection` ratio. This file is optional and is not required by the
explorer.

## 14. Build BASiCS assets for Streamlit

Entry point: `explorer/build_assets.py`.

### Input

- `basics_master.csv.gz` from step 9.
- `wormcat_supported.csv.gz` from step 12.

```bash
python explorer/build_assets.py \
  --master outputs/basics/RUN_ID/basics_master.csv.gz \
  --wormcat outputs/basics/RUN_ID/analysis/support/wormcat_supported.csv.gz \
  --output-dir outputs/explorer
```

### Output

```text
outputs/explorer/
├── basics.parquet
├── basics_summary.parquet
├── wormcat_supported.parquet
├── neuropeptide_family_signatures.csv
├── gpcr_ligand_receptor_pairs.csv
├── schema_manifest.json
└── clustering/                         # from step 3, when available
```

The builder validates required columns, unique cell-type/gene keys, and that
the WormCat input contains only the supported condition.

These assets use schema version 2; rebuild existing version-1 explorer assets.
Clustering assets retain schema version 1 and do not need re-exporting.
The app displays significant categories only, with category levels 1–3.

## 15. Launch the explorer

### Input

- The versioned assets from step 14.
- Optional clustering assets from step 3. Without them, tabs 1–7 work and the
  Clusters tab reports that no clustering assets were found.

```bash
streamlit run explorer/app.py
```

For assets stored elsewhere:

```bash
CENGEN_OUTPUT_DIR=/absolute/path/to/explorer/assets \
  streamlit run explorer/app.py
```

### Output

A local interactive web application with Browse, Matrix, Gene, Cell type,
Family, WormCat enrichment, Ligand–Receptor, and accepted Clusters tabs. Download
buttons generate tables in the user’s browser; the app does not modify analysis
results.

## 16. Validation before handoff or publication

```bash
# Parse R, compile Python, and check shell syntax.
find src scripts tests -name '*.R' -type f -print0 | \
  while IFS= read -r -d '' f; do Rscript -e "parse(file='$f')" >/dev/null; done
find explorer tests/python -name '*.py' -type f -print0 | \
  xargs -0 -n1 python -m py_compile
find slurm -name '*.sh' -type f -print0 | xargs -0 -n1 bash -n

# Unit and explorer-asset tests.
Rscript tests/testthat.R
python -m pytest -q
python tests/python/smoke_build_assets.py

# Confirm active code has no machine-specific source paths.
rg '/Users/|/sci/|MY_R|Zaslab/(entropy|Entropy)' src scripts slurm explorer
```

The final `rg` command should print nothing. Before making the repository
public, replace the restrictive placeholder `LICENSE` with the approved project
license and verify that no raw or generated data are staged.
