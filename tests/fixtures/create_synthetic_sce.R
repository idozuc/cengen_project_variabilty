#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(SingleCellExperiment))
set.seed(1)
genes <- sprintf("WBGene%08d", seq_len(220L))
cells <- sprintf("cell_%03d", seq_len(60L))
counts <- matrix(rpois(length(genes) * length(cells), lambda = 1.5), nrow = length(genes), dimnames = list(genes, cells))
sce <- SingleCellExperiment(
  assays = list(counts = counts),
  rowData = S4Vectors::DataFrame(gene_short_name = sprintf("gene-%03d", seq_along(genes))),
  colData = S4Vectors::DataFrame(
    Cell.type = rep("TEST", length(cells)),
    Experiment = rep(c("exp1", "exp2", "exp3"), each = 20L),
    Size_Factor = rep(1, length(cells)),
    Detection = colSums(counts > 0),
    pct_counts_Mito = rep(0, length(cells)),
    total_features_by_counts = colSums(counts > 0),
    total_counts = colSums(counts)
  )
)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: create_synthetic_sce.R OUTPUT.rds")
saveRDS(sce, args[[1L]])
