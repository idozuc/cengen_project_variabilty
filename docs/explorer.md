# Streamlit explorer

Build assets from the current canonical master and supported WormCat table:

```bash
python explorer/build_assets.py \
  --master outputs/basics/RUN/basics_master.csv.gz \
  --wormcat outputs/basics/RUN/analysis/support/wormcat_supported.csv.gz \
  --output-dir outputs/explorer
```

Export accepted clustering results separately with
`scripts/clustering/export_explorer_assets.R`. Launch with:

```bash
streamlit run explorer/app.py
```

Set `CENGEN_OUTPUT_DIR` only when assets live outside `outputs/explorer`.
The application refuses unsupported schema versions or missing required
columns. It provides browse, matrix, gene, cell-type, family, WormCat,
ligand–receptor, and accepted-cluster tabs.

Rebuild assets after upgrading: the explorer requires schema version 2 and
`wormcat_supported.parquet`. No GO annotation files are needed. WormCat views
show significant supported categories and offer levels 1–3; an empty result
shows a message without preventing other tabs from initializing.
