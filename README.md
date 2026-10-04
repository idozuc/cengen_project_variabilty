# CeNGEN within-cell-type analysis

This repository contains the supported CeNGEN workflows for nuisance-adjusted
within-cell-type clustering, BASiCS HVG/LVG analysis, and a Streamlit explorer.
It is source-first: raw data and generated results are intentionally excluded.

## Quick start

1. Restore R packages with `renv::restore()` and install Python dependencies
   with `python -m pip install -e '.[test]'`.
2. Copy an example YAML from `config/`, replace its input paths, and keep the
   real path file untracked.
3. Run a one-cell-type clustering screen:

   ```bash
   Rscript scripts/clustering/run_clustering.R \
     --sce /path/to/L4_neuron_sce.rds --cell-type SMD --run-id smoke --screen-only
   ```

4. Prepare a BASiCS run and fit one task locally:

   ```bash
   Rscript scripts/basics/01_prepare.R --config config/basics.yaml
   Rscript scripts/basics/02_fit_cell_type.R --run-dir outputs/basics/example --task-id 1
   ```

5. Build explorer assets and launch the app:

   ```bash
   python explorer/build_assets.py --config config/explorer.yaml
   streamlit run explorer/app.py
   ```

See `docs/architecture.md`, `docs/clustering-workflow.md`,
`docs/basics-workflow.md`, and `docs/explorer.md` before a production run.
For a single end-to-end sequence with every input and output, use
[`RUNBOOK.md`](RUNBOOK.md).

## Supported products

- `scripts/clustering`: statistical screening, accepted resolutions, markers,
  visualization bundles, and portable explorer exports.
- `scripts/basics`: checkpointed BASiCS fits, convergence diagnostics, LVG
  calls, experiment support, and a canonical master table.
- `scripts/analysis`: mappings and WormCat enrichment.
- `explorer`: versioned asset builder and Streamlit UI.

The R/Shiny cluster viewer remains an optional diagnostic; Streamlit is the
supported combined explorer.
