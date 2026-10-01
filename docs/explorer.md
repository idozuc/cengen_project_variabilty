# Streamlit explorer

Build assets from the current canonical master and supported GO table:

```bash
python explorer/build_assets.py \
  --master outputs/basics/RUN/basics_master.csv.gz \
  --go outputs/basics/RUN/analysis/go_bp_supported.csv.gz \
  --output-dir outputs/explorer
```

Export accepted clustering results separately with
`scripts/clustering/export_explorer_assets.R`. Launch with:

```bash
streamlit run explorer/app.py
```

Set `CENGEN_OUTPUT_DIR` only when assets live outside `outputs/explorer`.
The application refuses unsupported schema versions or missing required
columns. It provides browse, matrix, gene, cell-type, family, GO,
ligand–receptor, and accepted-cluster tabs.
