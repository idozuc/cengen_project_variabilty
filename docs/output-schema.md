# Output contracts

## BASiCS master

One row per `cell_type` and `stable_id`. Required explorer fields include
`gene_name`, `Mu`, `Delta`, `Epsilon`, `Prob`, `HVG`, `LVG`, `n_cells`,
`n_experiments`, `stress_gene`, and `analysis_eligible`.

## Supported GO

One row per cell type and term, with `condition=supported`, `term_id`, `term`,
`ontology`, gene counts, fold enrichment, p/q values, and `significant`.

## Clustering explorer assets

`manifest.parquet` contains schema version 1 and paths to each accepted cell
type’s `cells.parquet`, `expression.parquet`, `summary.parquet`,
`qc_by_resolution.parquet`, and `markers.parquet`. Cell IDs must match exactly
between the cells and expression tables.
