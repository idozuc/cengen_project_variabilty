#!/usr/bin/env Rscript

# Check whether each BASiCS HVG is detectably expressed in multiple experiments.
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
})

script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_file[[1L]]))
project_root <- dirname(dirname(dirname(script_path)))
source(file.path(project_root, "src", "R", "project_io.R"))
options <- parse_named_arguments(
  commandArgs(trailingOnly = TRUE),
  list(
    sce = Sys.getenv("BASICS_RDS", unset = NA_character_),
    run_dir = Sys.getenv("BASICS_RUN", unset = NA_character_)
  ),
  "Usage: Rscript scripts/basics/05_assess_experiment_support.R --sce PATH --run-dir PATH [--config PATH]"
)
sce_file <- normalizePath(require_path_option(options$sce, "--sce"))
run_dir <- normalizePath(require_path_option(options$run_dir, "--run-dir"))
config <- readRDS(file.path(run_dir, "config.rds"))
if (unname(tools::md5sum(sce_file)) != config$input_md5) {
  stop("SCE checksum does not match the BASiCS run")
}

sce <- readRDS(sce_file)
counts <- assay(sce, "counts")
cell_type <- as.character(colData(sce)$Cell.type)
experiment <- as.character(colData(sce)$Experiment)
tasks <- read.delim(file.path(run_dir, "cell_types.tsv"), stringsAsFactors = FALSE)

# An experiment supports a gene when at least three or 5% of its cells detect it.
measure_type <- function(task) {
  result_dir <- basics_result_dir(run_dir, task$task_id, task$cell_type)
  genes <- read.csv(file.path(result_dir, "genes.csv"), stringsAsFactors = FALSE)
  genes <- genes[genes$HVG, c("stable_id", "gene_name")]
  if (!nrow(genes)) return(NULL)

  experiment_table <- read.csv(
    file.path(result_dir, "experiments.csv"), stringsAsFactors = FALSE
  )
  included <- experiment_table$experiment[experiment_table$included]

  do.call(rbind, lapply(included, function(batch) {
    cells <- which(cell_type == task$cell_type & experiment == batch)
    n_detected <- Matrix::rowSums(counts[genes$stable_id, cells, drop = FALSE] > 0)
    threshold <- max(3L, ceiling(0.05 * length(cells)))
    data.frame(
      cell_type = task$cell_type,
      stable_id = genes$stable_id,
      gene_name = genes$gene_name,
      experiment = batch,
      n_cells = length(cells),
      n_detected = as.integer(n_detected),
      detection_fraction = as.numeric(n_detected) / length(cells),
      support_threshold = threshold,
      detection_supported = as.numeric(n_detected) >= threshold
    )
  }))
}

detail <- do.call(rbind, lapply(seq_len(nrow(tasks)), function(i) {
  measure_type(tasks[i, ])
}))

# Summarize experiment support without changing the original HVG calls.
groups <- split(detail, interaction(detail$cell_type, detail$stable_id, drop = TRUE))
summary <- do.call(rbind, lapply(groups, function(x) {
  total_detected <- sum(x$n_detected)
  data.frame(
    cell_type = x$cell_type[1],
    stable_id = x$stable_id[1],
    gene_name = x$gene_name[1],
    n_eligible_experiments = nrow(x),
    n_detection_supported_experiments = sum(x$detection_supported),
    max_experiment_detection_share = if (total_detected) {
      max(x$n_detected) / total_detected
    } else {
      NA_real_
    },
    robust_detection_support = sum(x$detection_supported) >= 2L
  )
}))

detail <- detail[order(detail$cell_type, detail$stable_id, detail$experiment), ]
summary <- summary[order(summary$cell_type, summary$stable_id), ]
write.csv(detail, file.path(run_dir, "hvg_experiment_detection.csv"), row.names = FALSE)
write.csv(summary, file.path(run_dir, "hvg_experiment_support.csv"), row.names = FALSE)

cat(sprintf(
  "%d/%d HVGs have detection support in at least two experiments\n",
  sum(summary$robust_detection_support), nrow(summary)
))
