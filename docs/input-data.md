# Input data

The primary input is an RDS containing a `SingleCellExperiment` with a `counts`
assay. Row names must be unique stable WormBase gene IDs. The following cell
metadata are required by clustering: `Cell.type`, `Experiment`, `Size_Factor`,
`Detection`, `pct_counts_Mito`, and `total_features_by_counts`. BASiCS requires
`Cell.type` and `Experiment`; `rowData(sce)$gene_short_name` supplies symbols.

Raw inputs are never committed. Put local files under `data/raw/` or reference
an absolute shared path from an untracked configuration file. Production runs
record the input checksum. WormCat, CeNGEN supplement, GAF, and OBO annotations
are likewise external inputs whose source/version/checksum must be recorded.
