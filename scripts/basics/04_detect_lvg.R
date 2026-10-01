#!/usr/bin/env Rscript

# Call lowly variable genes from one saved BASiCS chain.
suppressPackageStartupMessages(library(BASiCS))

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
  "Usage: Rscript scripts/basics/04_detect_lvg.R --run-dir PATH --task-id INTEGER [--config PATH]"
)
run_dir <- require_path_option(options$run_dir, "--run-dir")
task_id <- as.integer(options$task_id)
if (is.na(task_id)) stop("--task-id must be an integer")

config <- readRDS(file.path(run_dir, "config.rds"))
tasks <- read.delim(file.path(run_dir, "cell_types.tsv"), stringsAsFactors = FALSE)
if (task_id < 1L || task_id > nrow(tasks)) stop("Invalid array task ID")

task <- tasks[task_id, ]
out_dir <- basics_result_dir(run_dir, task$task_id, task$cell_type)
if (!file.exists(file.path(out_dir, "_COMPLETE"))) stop("HVG fit is incomplete")
if (file.exists(file.path(out_dir, "_LVG_COMPLETE"))) quit(save = "no")

chain <- readRDS(file.path(out_dir, "chain.rds"))

# Mirror the HVG analysis by testing the lower 10% of residual variability.
lvg <- as.data.frame(BASiCS_DetectLVG(
  chain,
  PercentileThreshold = 0.10,
  EFDR = config$efdr,
  MinESS = config$min_ess,
  Plot = FALSE
)@Table)
if (!"GeneName" %in% names(lvg)) stop("BASiCS output lacks GeneName")

# Restore gene annotations and attach the existing convergence diagnostics.
hvg <- read.csv(file.path(out_dir, "genes.csv"), stringsAsFactors = FALSE)
diagnostics <- read.csv(
  file.path(out_dir, "diagnostics.csv"),
  stringsAsFactors = FALSE
)
lvg$stable_id <- as.character(lvg$GeneName)
if (!setequal(lvg$stable_id, hvg$stable_id)) stop("LVG and HVG genes differ")

metadata <- c(
  "gene_name", "detected_cells", "detected_experiments", "cell_type",
  "evidence_tier", "seed", "n_cells", "n_experiments"
)
index <- match(lvg$stable_id, hvg$stable_id)
for (column in metadata) lvg[[column]] <- hvg[[column]][index]

diagnostic_index <- match(lvg$stable_id, diagnostics$stable_id)
if (anyNA(diagnostic_index)) stop("Missing diagnostics for LVG genes")
lvg$epsilon_ess <- diagnostics$epsilon_ess[diagnostic_index]
lvg$epsilon_geweke_z <- diagnostics$epsilon_geweke_z[diagnostic_index]

# Rename the temporary file atomically; the marker is written last.
output <- file.path(out_dir, "lvg.csv")
temporary <- paste0(output, ".tmp.", Sys.getpid())
on.exit(unlink(temporary), add = TRUE)
write.csv(lvg, temporary, row.names = FALSE)
if (!file.rename(temporary, output)) stop("Could not finalize lvg.csv")
writeLines(
  format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  file.path(out_dir, "_LVG_COMPLETE")
)

cat(sprintf(
  "Completed %s: %d genes, %d LVGs\n",
  task$cell_type,
  nrow(lvg),
  sum(lvg$LVG, na.rm = TRUE)
))
