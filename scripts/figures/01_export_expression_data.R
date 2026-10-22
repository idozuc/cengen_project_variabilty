#!/usr/bin/env Rscript

# Export cell-level raw counts for the histogram, marker-evidence and
# marker-UMAP figures (Figures 3, 4, 5, 10, 11).
#
# URX, AVK and RIG cells follow the BASiCS inclusion rule (experiments with
# >= 10 cells). The pooled ASE and AWC cells are exactly those retained in the
# saved positive-control analysis, with their known identities.
#
# Writes OUTPUT_DIR/data/expression_histogram_data.csv.gz and
#        OUTPUT_DIR/data/expression_histogram_manifest.csv

script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_file[[1L]]))
project_root <- dirname(dirname(dirname(script_path)))
source(file.path(project_root, "src", "R", "project_io.R"))
source(file.path(project_root, "src", "R", "figures.R"))

options <- figure_options(
  "Usage: Rscript scripts/figures/01_export_expression_data.R --config config/figures.yaml"
)
data_dir <- options$data_dir

sce <- read_sce(options)
master <- utils::read.csv(
  gzfile(input_file(options$basics_dir, "basics_master_with_detection.csv.gz")),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
positive_controls <- readRDS(input_file(options$positive_controls))

counts <- SummarizedExperiment::assay(sce, "counts")
metadata <- as.data.frame(SummarizedExperiment::colData(sce))
rownames(metadata) <- colnames(sce)

targets <- data.frame(
  dataset = c(
    "URX", "URX", "URX", "ASE", "AWC", "AVK", "RIG",
    "AWC", "ASE", "ASE", "ASE"
  ),
  panel_order = seq_len(11L),
  stable_id = c(
    "WBGene00001458", # flp-15: HVG, lower zero fraction
    "WBGene00003756", # nlp-18: HVG, higher zero fraction
    "WBGene00004429", # rpl-17: LVG, mean matched to both HVGs
    "WBGene00001533", # gcy-6: ASEL marker
    "WBGene00006070", # str-2: AWC_ON marker
    "WBGene00020317", # pdf-1: high-mean, low-zero HVG in AVK
    "WBGene00017671", # pgal-1: lowest-zero HVG in RIG
    "WBGene00022408", # srsx-3: AWC_OFF marker
    "WBGene00001534", # gcy-7: ASEL marker
    "WBGene00001547", # gcy-22: ASER marker
    "WBGene00001532"  # gcy-5: ASER marker
  ),
  gene_name = c(
    "flp-15", "nlp-18", "rpl-17", "gcy-6", "str-2", "pdf-1", "pgal-1",
    "srsx-3", "gcy-7", "gcy-22", "gcy-5"
  ),
  expected_call = c(
    "HVG", "HVG", "LVG", "marker", "marker", "HVG", "HVG",
    "marker", "marker", "marker", "marker"
  ),
  description = c(
    "HVG neuropeptide with lower zero fraction",
    "HVG neuropeptide with higher zero fraction",
    "LVG with comparable mean expression",
    "ASEL marker in pooled ASEL + ASER cells",
    "AWC_ON marker in pooled AWC_ON + AWC_OFF cells",
    "High-mean, low-zero neuropeptide HVG",
    "Lowest-zero HVG across cell types",
    "AWC_OFF marker in pooled AWC_ON + AWC_OFF cells",
    "ASEL marker in pooled ASEL + ASER cells",
    "ASER marker in pooled ASEL + ASER cells",
    "ASER marker in pooled ASEL + ASER cells"
  ),
  stringsAsFactors = FALSE
)

if (!all(targets$stable_id %in% rownames(counts))) {
  missing_genes <- targets$stable_id[!targets$stable_id %in% rownames(counts)]
  stop("Target genes missing from the SCE: ", paste(missing_genes, collapse = ", "))
}

get_master_row <- function(cell_type, stable_id) {
  row <- master[
    master$cell_type == cell_type & master$stable_id == stable_id,
    ,
    drop = FALSE
  ]
  if (nrow(row) != 1L) {
    stop("Expected one master row for ", cell_type, " / ", stable_id)
  }
  row
}

make_rows <- function(cells, target, identity) {
  if (!all(cells %in% colnames(counts))) {
    stop(target$dataset, ": retained cells are missing from the SCE")
  }
  values <- as.numeric(counts[target$stable_id, cells])
  data.frame(
    dataset = target$dataset,
    panel_order = target$panel_order,
    stable_id = target$stable_id,
    gene_name = target$gene_name,
    expected_call = target$expected_call,
    description = target$description,
    cell = cells,
    identity = as.character(identity),
    experiment = as.character(metadata[cells, "Experiment"]),
    raw_count = values,
    log_count = log1p(values),
    stringsAsFactors = FALSE
  )
}

# BASiCS inclusion rule: keep experiments contributing >= 10 cells.
basics_cells <- function(cell_type) {
  all_cells <- rownames(metadata)[metadata$Cell.type == cell_type]
  experiment_sizes <- table(metadata[all_cells, "Experiment"])
  included <- names(experiment_sizes)[experiment_sizes >= 10L]
  all_cells[metadata[all_cells, "Experiment"] %in% included]
}

expression_rows <- list()

# URX matched HVG/LVG examples, then the AVK and RIG high-mean, low-zero HVGs.
for (index in c(1L, 2L, 3L, 6L, 7L)) {
  cell_type <- targets$dataset[index]
  cells <- basics_cells(cell_type)
  expression_rows[[length(expression_rows) + 1L]] <- make_rows(
    cells,
    targets[index, , drop = FALSE],
    rep(cell_type, length(cells))
  )
}

# Positive controls: use the exact saved cell sets and known identities.
control_specs <- list(
  ASE_gcy_6 = list(result = "gmm_ASEL_ASER", target_index = 4L),
  AWC_str_2 = list(result = "gmm_AWC_pair", target_index = 5L),
  AWC_srsx_3 = list(result = "gmm_AWC_pair", target_index = 8L),
  ASE_gcy_7 = list(result = "gmm_ASEL_ASER", target_index = 9L),
  ASE_gcy_22 = list(result = "gmm_ASEL_ASER", target_index = 10L),
  ASE_gcy_5 = list(result = "gmm_ASEL_ASER", target_index = 11L)
)

for (dataset in names(control_specs)) {
  specification <- control_specs[[dataset]]
  control <- positive_controls[[specification$result]]
  if (is.null(control)) stop("Missing positive control: ", specification$result)

  # true_identity is named by cell and covers every retained cell, also when
  # the pipeline selected G = 1.
  truth <- control$true_identity
  cells <- names(truth)
  if (is.null(cells)) stop(dataset, ": saved identities have no cell names")

  expression_rows[[length(expression_rows) + 1L]] <- make_rows(
    cells,
    targets[specification$target_index, , drop = FALSE],
    truth
  )
}

expression_data <- do.call(rbind, expression_rows)

# BASiCS statistics for the non-marker examples. Marker statistics are
# calculated by the Python scripts for each known identity.
targets$Mu <- NA_real_
targets$Epsilon <- NA_real_
targets$pct_zero_master <- NA_real_

for (index in c(1L, 2L, 3L, 6L, 7L)) {
  row <- get_master_row(targets$dataset[index], targets$stable_id[index])
  targets$Mu[index] <- row$Mu
  targets$Epsilon[index] <- row$Epsilon
  targets$pct_zero_master[index] <- row$pct_zero
}

targets$n_cells <- vapply(targets$panel_order, function(panel_order) {
  sum(expression_data$panel_order == panel_order)
}, integer(1L))
targets$pct_zero_observed <- vapply(targets$panel_order, function(panel_order) {
  values <- expression_data$raw_count[expression_data$panel_order == panel_order]
  100 * mean(values == 0)
}, numeric(1L))

connection <- gzfile(file.path(data_dir, "expression_histogram_data.csv.gz"), open = "wt")
utils::write.csv(expression_data, connection, row.names = FALSE)
close(connection)
utils::write.csv(
  targets,
  file.path(data_dir, "expression_histogram_manifest.csv"),
  row.names = FALSE,
  na = ""
)

cat(sprintf(
  "Exported %s expression observations for %s panels.\n",
  format(nrow(expression_data), big.mark = ","),
  nrow(targets)
))
