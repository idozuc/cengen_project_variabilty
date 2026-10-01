#!/usr/bin/env Rscript

# Calculate convergence diagnostics after the expensive fits finish.
suppressPackageStartupMessages({
  library(BASiCS)
  library(coda)
})

script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_file[[1L]]))
project_root <- dirname(dirname(dirname(script_path)))
source(file.path(project_root, "src", "R", "project_io.R"))
options <- parse_named_arguments(
  commandArgs(trailingOnly = TRUE),
  list(run_dir = Sys.getenv("BASICS_RUN", unset = NA_character_)),
  "Usage: Rscript scripts/basics/03_diagnose.R --run-dir PATH [--config PATH]"
)
run_dir <- require_path_option(options$run_dir, "--run-dir")
tasks <- read.delim(file.path(run_dir, "cell_types.tsv"), stringsAsFactors = FALSE)

summaries <- lapply(seq_len(nrow(tasks)), function(i) {
  task <- tasks[i, ]
  out_dir <- basics_result_dir(run_dir, task$task_id, task$cell_type)
  complete <- file.exists(file.path(out_dir, "_COMPLETE"))
  if (!complete) {
    return(data.frame(task_id = task$task_id, cell_type = task$cell_type,
                      complete = FALSE, median_epsilon_ess = NA,
                      fraction_epsilon_ess_below_100 = NA))
  }

  chain <- readRDS(file.path(out_dir, "chain.rds"))
  epsilon <- as.matrix(displayChainBASiCS(chain, "epsilon"))
  ess <- as.numeric(effectiveSize(mcmc(epsilon)))
  geweke <- as.numeric(geweke.diag(mcmc(epsilon))$z)
  diagnostics <- data.frame(
    stable_id = colnames(epsilon),
    epsilon_ess = ess,
    epsilon_geweke_z = geweke,
    epsilon_ess_below_100 = ess < 100,
    epsilon_abs_geweke_above_2 = abs(geweke) > 2
  )
  write.csv(diagnostics, file.path(out_dir, "diagnostics.csv"), row.names = FALSE)

  data.frame(
    task_id = task$task_id,
    cell_type = task$cell_type,
    complete = TRUE,
    median_epsilon_ess = median(ess, na.rm = TRUE),
    fraction_epsilon_ess_below_100 = mean(ess < 100, na.rm = TRUE)
  )
})

summary <- do.call(rbind, summaries)
write.csv(summary, file.path(run_dir, "diagnostic_summary.csv"), row.names = FALSE)
cat(sprintf("Diagnosed %d/%d completed cell types\n",
            sum(summary$complete), nrow(summary)))
