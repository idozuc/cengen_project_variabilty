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

The CeNGEN-filtered enrichment also needs the CeNGEN gene × cell-type
expression call table at threshold 2 (the CeNGEN file
`021821_medium_threshold2.csv`; record its source): one row per gene with `Wormbase_ID`,
`gene_name`, and one column per cell type, where a value above zero means the
gene is called expressed.
