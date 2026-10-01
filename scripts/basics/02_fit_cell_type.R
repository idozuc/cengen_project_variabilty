#!/usr/bin/env Rscript

# Fit BASiCS for the cell type assigned to this SLURM array task.
suppressPackageStartupMessages({
  library(BASiCS)
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
    run_dir = Sys.getenv("BASICS_RUN", unset = NA_character_),
    task_id = Sys.getenv("SLURM_ARRAY_TASK_ID", unset = NA_character_)
  ),
  "Usage: Rscript scripts/basics/02_fit_cell_type.R --run-dir PATH --task-id INTEGER [--config PATH]"
)
run_dir <- require_path_option(options$run_dir, "--run-dir")
task_id <- as.integer(options$task_id)
if (is.na(task_id)) stop("--task-id must be an integer")

config <- readRDS(file.path(run_dir, "config.rds"))
tasks <- read.delim(file.path(run_dir, "cell_types.tsv"), stringsAsFactors = FALSE)
if (task_id < 1L || task_id > nrow(tasks)) stop("Invalid array task ID")
task <- tasks[task_id, ]
cell_type <- task$cell_type
out_dir <- basics_result_dir(run_dir, task_id, cell_type)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (file.exists(file.path(out_dir, "_COMPLETE"))) quit(save = "no")

# Refuse to run if the registered input file changed after preparation.
if (unname(tools::md5sum(config$input)) != config$input_md5) {
  stop("Input checksum differs from the prepared run")
}

sce <- readRDS(config$input)
counts_all <- assay(sce, "counts")
cell_types <- as.character(colData(sce)$Cell.type)
experiments <- as.character(colData(sce)$Experiment)
gene_names <- as.character(rowData(sce)$gene_short_name)

# Retain adequately represented experiments within this cell type.
type_cells <- which(cell_types == cell_type)
experiment_sizes <- table(experiments[type_cells])
included_experiments <- names(experiment_sizes)[
  experiment_sizes >= config$min_cells_per_experiment
]
selected <- type_cells[experiments[type_cells] %in% included_experiments]
batch <- droplevels(factor(experiments[selected]))
if (nlevels(batch) < config$min_experiments) stop("Too few eligible experiments")

# Retain genes observed across cells and in more than one experiment.
detected <- counts_all[, selected, drop = FALSE] > 0
detected_cells <- Matrix::rowSums(detected)
detected_experiments <- Reduce(`+`, lapply(levels(batch), function(x) {
  Matrix::rowSums(detected[, batch == x, drop = FALSE]) > 0
}))
keep <- detected_cells >= config$min_detected_cells &
  detected_experiments >= config$min_detected_experiments

counts <- as.matrix(counts_all[keep, selected, drop = FALSE])
gene_id <- rownames(counts_all)[keep]
gene_name <- gene_names[keep]
rownames(counts) <- gene_id
colnames(counts) <- sprintf("cell_%06d", seq_len(ncol(counts)))

fit_data <- SingleCellExperiment(
  assays = list(counts = counts),
  rowData = S4Vectors::DataFrame(gene_id = gene_id, gene_name = gene_name),
  colData = S4Vectors::DataFrame(BatchInfo = batch)
)

# Fit residual overdispersion after accounting for mean and experiment batch.
set.seed(task$seed)
prior <- BASiCS_PriorParam(fit_data, PriorMu = "EmpiricalBayes")
started <- Sys.time()
chain <- BASiCS_MCMC(
  Data = fit_data,
  N = config$n_iterations,
  Burn = config$burn,
  Thin = config$thin,
  Regression = TRUE,
  WithSpikes = FALSE,
  PriorParam = prior,
  PrintProgress = FALSE,
  Threads = 1
)

# Save the chain before post-processing so a later failure does not lose it.
saveRDS(chain, file.path(out_dir, "chain.rds"))

hvg <- as.data.frame(BASiCS_DetectHVG(
  chain,
  PercentileThreshold = config$percentile_threshold,
  EFDR = config$efdr,
  MinESS = config$min_ess,
  Plot = FALSE
)@Table)
if (!"GeneName" %in% names(hvg)) stop("BASiCS output lacks GeneName")

# GeneName contains the stable ID because those IDs were used as row names.
hvg$stable_id <- as.character(hvg$GeneName)
hvg$gene_name <- gene_name[match(hvg$stable_id, gene_id)]
hvg$detected_cells <- detected_cells[match(hvg$stable_id, rownames(counts_all))]
hvg$detected_experiments <- detected_experiments[match(hvg$stable_id, rownames(counts_all))]
hvg$cell_type <- cell_type
hvg$evidence_tier <- task$evidence_tier
hvg$seed <- task$seed
hvg$n_cells <- length(selected)
hvg$n_experiments <- nlevels(batch)
write.csv(hvg, file.path(out_dir, "genes.csv"), row.names = FALSE)

write.csv(
  data.frame(experiment = names(experiment_sizes),
             n_cells = as.integer(experiment_sizes),
             included = names(experiment_sizes) %in% included_experiments),
  file.path(out_dir, "experiments.csv"), row.names = FALSE
)
writeLines(c(
  paste("cell_type", cell_type, sep = "\t"),
  paste("seed", task$seed, sep = "\t"),
  paste("elapsed_minutes", round(as.numeric(difftime(Sys.time(), started, units = "mins")), 2), sep = "\t")
), file.path(out_dir, "fit_info.tsv"))
writeLines(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), file.path(out_dir, "_COMPLETE"))

cat(sprintf("Completed %s: %d cells, %d genes, %d experiments\n",
            cell_type, ncol(counts), nrow(counts), nlevels(batch)))
