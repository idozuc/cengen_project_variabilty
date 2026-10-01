#!/usr/bin/env Rscript

# Add zero percentage and a readable experiment-detection ratio to a BASiCS master table.
script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_file[[1L]]))
project_root <- dirname(dirname(dirname(script_path)))
source(file.path(project_root, "src", "R", "project_io.R"))
options <- parse_named_arguments(
  commandArgs(trailingOnly = TRUE),
  list(input = NA_character_, output = NA_character_),
  "Usage: Rscript scripts/analysis/04_add_detection_columns.R --input FILE.csv.gz --output FILE.csv.gz [--config PATH]"
)
input_file <- normalizePath(require_path_option(options$input, "--input"), mustWork = TRUE)
output_file <- require_path_option(options$output, "--output", must_exist = FALSE)
master <- read.csv(gzfile(input_file), stringsAsFactors = FALSE)

required <- c(
  "detected_cells", "n_cells", "detected_experiments", "n_experiments"
)
missing <- setdiff(required, names(master))
if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "))

# Calculate the fraction of cells with zero counts for each cell-type/gene pair.
master$pct_zero <- 100 * (1 - master$detected_cells / master$n_cells)

# Keep the existing numerator and denominator, and add a compact label such as 3/6.
master$experiment_detection <- paste0(
  master$detected_experiments, "/", master$n_experiments
)

if (any(!is.finite(master$pct_zero)) || any(master$pct_zero < 0 | master$pct_zero > 100)) {
  stop("Invalid pct_zero values")
}

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
connection <- gzfile(output_file, open = "wt", encoding = "UTF-8")
tryCatch(
  write.csv(master, connection, row.names = FALSE, na = "NA"),
  finally = close(connection)
)

# Confirm that the new file has the same rows and both added columns.
written <- read.csv(gzfile(output_file), stringsAsFactors = FALSE)
if (nrow(written) != nrow(master) || !identical(names(written), names(master))) {
  stop("Saved table failed validation")
}

cat("Saved", format(nrow(master), big.mark = ","), "rows to", output_file, "\n")
