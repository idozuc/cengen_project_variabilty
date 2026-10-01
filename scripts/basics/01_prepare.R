#!/usr/bin/env Rscript

# Freeze the cell-type list and settings before submitting the array.
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
})

script_path <- script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_file[[1L]]))
project_root <- dirname(dirname(dirname(script_path)))
source(file.path(project_root, "src", "R", "project_io.R"))

options <- parse_named_arguments(
  commandArgs(trailingOnly = TRUE),
  list(
    sce = Sys.getenv("BASICS_RDS", unset = NA_character_),
    run_dir = Sys.getenv("BASICS_RUN", unset = NA_character_),
    min_cells_per_experiment = 10L,
    min_experiments = 3L,
    strong_cells_per_experiment = 20L,
    min_detected_cells = 10L,
    min_detected_experiments = 2L,
    n_iterations = 20000L,
    burn = 10000L,
    thin = 10L,
    percentile_threshold = 0.90,
    efdr = 0.10,
    min_ess = 100L
  ),
  paste(
    "Usage: Rscript scripts/basics/01_prepare.R --sce PATH --run-dir PATH",
    "[--config config/basics.yaml] [--n-iterations N --burn N --thin N]"
  )
)
input <- require_path_option(options$sce, "--sce")
run_dir <- require_path_option(options$run_dir, "--run-dir", must_exist = FALSE)
if (dir.exists(run_dir) && length(list.files(run_dir, all.files = TRUE, no.. = TRUE))) {
  stop("BASICS_RUN must be new or empty: ", run_dir)
}

# These settings define one fixed analysis run.
config <- list(
  input = normalizePath(input),
  input_md5 = unname(tools::md5sum(input)),
  min_cells_per_experiment = as.integer(options$min_cells_per_experiment),
  min_experiments = as.integer(options$min_experiments),
  strong_cells_per_experiment = as.integer(options$strong_cells_per_experiment),
  min_detected_cells = as.integer(options$min_detected_cells),
  min_detected_experiments = as.integer(options$min_detected_experiments),
  n_iterations = as.integer(options$n_iterations),
  burn = as.integer(options$burn),
  thin = as.integer(options$thin),
  percentile_threshold = as.numeric(options$percentile_threshold),
  efdr = as.numeric(options$efdr),
  min_ess = as.integer(options$min_ess),
  created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
)
if (anyNA(unlist(config[c("n_iterations", "burn", "thin")]))) {
  stop("BASICS_N, BASICS_BURN, and BASICS_THIN must be integers")
}
if (config$burn >= config$n_iterations ||
    config$n_iterations %% config$thin != 0L || config$burn %% config$thin != 0L) {
  stop("Require burn < iterations and both divisible by thin")
}

sce <- readRDS(input)
if (!inherits(sce, "SingleCellExperiment")) stop("Input is not a SingleCellExperiment")
if (!all(c("Cell.type", "Experiment") %in% colnames(colData(sce)))) {
  stop("Input needs Cell.type and Experiment in colData")
}
counts <- assay(sce, "counts")
gene_id <- rownames(counts)
if (is.null(gene_id) || anyNA(gene_id) || any(!nzchar(gene_id)) || anyDuplicated(gene_id)) {
  stop("Count row names must be stable, unique gene IDs")
}

cell_type <- as.character(colData(sce)$Cell.type)
experiment <- as.character(colData(sce)$Experiment)
tab <- as.data.frame(table(cell_type, experiment), stringsAsFactors = FALSE)
names(tab)[3] <- "n_cells"
tab <- tab[tab$n_cells > 0, ]

# Keep types represented by at least three sufficiently large experiments.
rows <- lapply(sort(unique(cell_type)), function(type) {
  x <- tab[tab$cell_type == type, ]
  eligible <- x$n_cells >= config$min_cells_per_experiment
  strong <- x$n_cells >= config$strong_cells_per_experiment
  data.frame(
    cell_type = type,
    n_cells = sum(x$n_cells),
    n_eligible_cells = sum(x$n_cells[eligible]),
    n_eligible_experiments = sum(eligible),
    evidence_tier = if (sum(strong) >= config$min_experiments) {
      "strong_support"
    } else {
      "moderate_support"
    }
  )
})
tasks <- do.call(rbind, rows)
tasks <- tasks[tasks$n_eligible_experiments >= config$min_experiments, ]
tasks <- tasks[order(tasks$cell_type), ]
tasks$task_id <- seq_len(nrow(tasks))
tasks$seed <- 100000L + tasks$task_id
tasks <- tasks[, c("task_id", "cell_type", "evidence_tier", "seed",
                   "n_cells", "n_eligible_cells", "n_eligible_experiments")]

dir.create(run_dir, recursive = TRUE)
dir.create(file.path(run_dir, "cell_types"))
dir.create(file.path(run_dir, "logs"))
write.table(tasks, file.path(run_dir, "cell_types.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
write.table(tab, file.path(run_dir, "experiment_counts.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
saveRDS(config, file.path(run_dir, "config.rds"))
write_resolved_yaml(config, file.path(run_dir, "resolved_config.yaml"))
writeLines(capture.output(str(config)), file.path(run_dir, "manifest.txt"))
writeLines(capture.output(sessionInfo()), file.path(run_dir, "session_info.txt"))

cat(sprintf("Prepared %d tasks: %d strong, %d moderate\n",
            nrow(tasks), sum(tasks$evidence_tier == "strong_support"),
            sum(tasks$evidence_tier == "moderate_support")))
