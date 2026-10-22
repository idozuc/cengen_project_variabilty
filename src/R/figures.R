# Shared options and input helpers for scripts/figures/*.R.
# Paths are relative to the repository root, where the scripts are run.

figure_options <- function(usage) {
  options <- parse_named_arguments(
    commandArgs(trailingOnly = TRUE),
    list(
      sce = NA_character_,
      basics_dir = NA_character_,
      clustering_dir = NA_character_,
      positive_controls = NA_character_,
      output_dir = NA_character_
    ),
    usage
  )
  require_path_option(options$output_dir, "--output-dir", must_exist = FALSE)
  options$data_dir <- file.path(options$output_dir, "data")
  dir.create(options$data_dir, recursive = TRUE, showWarnings = FALSE)
  options
}

input_file <- function(...) {
  path <- file.path(...)
  if (!file.exists(path)) stop("Missing input: ", path, call. = FALSE)
  path
}

read_sce <- function(options) {
  suppressPackageStartupMessages(library(SingleCellExperiment))
  readRDS(require_path_option(options$sce, "--sce"))
}
