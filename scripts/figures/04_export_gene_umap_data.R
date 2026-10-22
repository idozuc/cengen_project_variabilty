#!/usr/bin/env Rscript

# Export the stored UMAP coordinates, GMM clusters and raw counts of one gene
# per cell type from the clustering run's visualization bundles (Figures 13-15).
#
# Writes OUTPUT_DIR/data/gene_umap_expression.csv

script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_file[[1L]]))
project_root <- dirname(dirname(dirname(script_path)))
source(file.path(project_root, "src", "R", "project_io.R"))
source(file.path(project_root, "src", "R", "figures.R"))

options <- figure_options(
  "Usage: Rscript scripts/figures/04_export_gene_umap_data.R --config config/figures.yaml"
)

targets <- list(HSN = "tph-1", SMD = "flp-5", NSM = "del-3")

export_gene <- function(cell_type, gene) {
  bundle <- readRDS(input_file(options$clustering_dir, cell_type, "visualization_bundle.rds"))

  cells <- names(bundle$primary_labels)
  if (length(cells) == 0L) cells <- rownames(bundle$embedding)

  embedding <- as.matrix(bundle$embedding)
  if (nrow(embedding) != length(cells) || ncol(embedding) < 2L) {
    stop(cell_type, ": stored UMAP embedding has unexpected dimensions")
  }
  if (!is.null(rownames(embedding))) embedding <- embedding[cells, , drop = FALSE]

  if (!setequal(colnames(bundle$counts), cells)) {
    stop(cell_type, ": count-matrix cells do not match clustered cells")
  }

  matched_rows <- which(
    bundle$gene_labels$feature == gene | bundle$gene_labels$gene_short_name == gene
  )
  if (length(matched_rows) == 0L) stop(cell_type, ": gene ", gene, " not found")
  raw_count <- as.numeric(
    Matrix::colSums(bundle$counts[matched_rows, cells, drop = FALSE])
  )

  data.frame(
    cell_type = cell_type,
    gene = gene,
    cell = cells,
    cluster = as.character(bundle$primary_labels[cells]),
    UMAP1 = embedding[, 1L],
    UMAP2 = embedding[, 2L],
    raw_count = raw_count,
    log_count = log1p(raw_count),
    stringsAsFactors = FALSE
  )
}

out <- do.call(rbind, Map(export_gene, names(targets), targets))
write.csv(out, file.path(options$data_dir, "gene_umap_expression.csv"), row.names = FALSE)
cat("Exported", nrow(out), "cells for", length(targets), "cell types\n")
