#!/usr/bin/env Rscript

# Build one stress-filtered BASiCS table from a completed run.
script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_file[[1L]]))
project_root <- dirname(dirname(dirname(script_path)))
source(file.path(project_root, "src", "R", "project_io.R"))
source(file.path(project_root, "src", "R", "stress_gene_sets.R"))
options <- parse_named_arguments(
  commandArgs(trailingOnly = TRUE),
  list(
    run_dir = Sys.getenv("BASICS_RUN", unset = NA_character_),
    output = NA_character_
  ),
  "Usage: Rscript scripts/basics/06_build_master.R --run-dir PATH --output FILE.csv.gz [--config PATH]"
)
run_dir <- normalizePath(require_path_option(options$run_dir, "--run-dir"), mustWork = TRUE)
output_file <- require_path_option(options$output, "--output", must_exist = FALSE)

# Stop when a table is missing columns needed by the analysis.
check_columns <- function(x, required, path) {
  missing <- setdiff(required, names(x))
  if (length(missing)) {
    stop("Missing columns in ", path, ": ", paste(missing, collapse = ", "))
  }
}

tasks <- read.delim(
  file.path(run_dir, "cell_types.tsv"),
  stringsAsFactors = FALSE
)

# Read and combine the BASiCS gene tables.
gene_files <- sort(list.files(
  file.path(run_dir, "cell_types"),
  pattern = "^genes\\.csv$",
  recursive = TRUE,
  full.names = TRUE
))
if (length(gene_files) != nrow(tasks)) {
  stop("Expected ", nrow(tasks), " genes.csv files, found ", length(gene_files))
}

gene_columns <- c(
  "stable_id", "gene_name", "cell_type", "Mu", "Delta", "Epsilon",
  "Prob", "HVG", "detected_cells", "detected_experiments",
  "evidence_tier", "seed", "n_cells", "n_experiments"
)
read_genes <- function(path) {
  x <- read.csv(path, stringsAsFactors = FALSE)
  check_columns(x, gene_columns, path)
  x[gene_columns]
}
master <- do.call(rbind, lapply(gene_files, read_genes))
rownames(master) <- NULL

key <- paste(master$cell_type, master$stable_id, sep = "::")
if (anyDuplicated(key)) stop("Duplicate cell-type/gene combinations found")
if (!setequal(unique(master$cell_type), tasks$cell_type)) {
  stop("Combined cell types do not match cell_types.tsv")
}

# Add the low-variability calls produced from the same BASiCS chains.
lvg_files <- sort(list.files(
  file.path(run_dir, "cell_types"),
  pattern = "^lvg\\.csv$",
  recursive = TRUE,
  full.names = TRUE
))
if (length(lvg_files) != nrow(tasks)) {
  stop("Expected ", nrow(tasks), " lvg.csv files, found ", length(lvg_files))
}

read_lvg <- function(path) {
  x <- read.csv(path, stringsAsFactors = FALSE)
  check_columns(x, c("cell_type", "stable_id", "LVG", "Prob"), path)
  names(x)[names(x) == "Prob"] <- "LVG_Prob"
  x[c("cell_type", "stable_id", "LVG", "LVG_Prob")]
}
lvg <- do.call(rbind, lapply(lvg_files, read_lvg))
rownames(lvg) <- NULL

lvg_key <- paste(lvg$cell_type, lvg$stable_id, sep = "::")
if (anyDuplicated(lvg_key)) stop("Duplicate cell-type/gene LVG rows found")
lvg_index <- match(key, lvg_key)
if (anyNA(lvg_index) || length(lvg_key) != length(key)) {
  stop("HVG and LVG tables do not contain the same genes")
}
master[c("LVG", "LVG_Prob")] <- lvg[lvg_index, c("LVG", "LVG_Prob")]
if (any(master$HVG %in% TRUE & master$LVG %in% TRUE)) {
  stop("A cell-type/gene row is called both HVG and LVG")
}

# Read and combine the per-gene MCMC diagnostics.
diagnostic_files <- sort(list.files(
  file.path(run_dir, "cell_types"),
  pattern = "^diagnostics\\.csv$",
  recursive = TRUE,
  full.names = TRUE
))
if (length(diagnostic_files) != nrow(tasks)) {
  stop(
    "Expected ", nrow(tasks), " diagnostics.csv files, found ",
    length(diagnostic_files)
  )
}

diagnostic_columns <- c(
  "epsilon_ess", "epsilon_geweke_z", "epsilon_ess_below_100",
  "epsilon_abs_geweke_above_2"
)
read_diagnostics <- function(path) {
  x <- read.csv(path, stringsAsFactors = FALSE)
  check_columns(x, c("stable_id", diagnostic_columns), path)
  x$cell_type <- sub("^[0-9]+_", "", basename(dirname(path)))
  x[c("cell_type", "stable_id", diagnostic_columns)]
}
diagnostics <- do.call(rbind, lapply(diagnostic_files, read_diagnostics))
rownames(diagnostics) <- NULL

diagnostic_key <- paste(diagnostics$cell_type, diagnostics$stable_id, sep = "::")
if (anyDuplicated(diagnostic_key)) stop("Duplicate diagnostic keys found")
diagnostic_index <- match(key, diagnostic_key)
if (anyNA(diagnostic_index) || length(diagnostic_key) != length(key)) {
  stop("Gene results and diagnostics do not match")
}
master[diagnostic_columns] <- diagnostics[diagnostic_index, diagnostic_columns]
master$good_mcmc <- is.finite(master$Epsilon) &
  is.finite(master$epsilon_ess) & master$epsilon_ess >= 100

# Use the same immutable stress-gene identifiers as the clustering workflow.
stress <- data.frame(stable_id = stress_genes_full_199, stringsAsFactors = FALSE)
if (nrow(stress) != 199L || anyDuplicated(stress$stable_id)) {
  stop("Embedded stress-gene list failed validation")
}
if (any(!grepl("^WBGene[0-9]+$", stress$stable_id))) {
  stop("Invalid WormBase ID in the stress-gene list")
}

# Flag stress genes and define the rows used by the primary analysis.
stress_index <- match(master$stable_id, stress$stable_id)
master$stress_gene <- !is.na(stress_index)
master$stress_list_name <- NA_character_
master$analysis_eligible <- master$good_mcmc & !master$stress_gene

stress_overlap <- intersect(stress$stable_id, unique(master$stable_id))
if (length(stress_overlap) != 109L) {
  stop("Expected 109 modeled stress genes, found ", length(stress_overlap))
}

# Put identifiers, estimates, diagnostics, and analysis flags in a fixed order.
column_order <- c(
  "cell_type", "stable_id", "gene_name", "Mu", "Delta", "Epsilon",
  "Prob", "HVG", "LVG_Prob", "LVG", "detected_cells",
  "detected_experiments", "n_cells", "n_experiments", "evidence_tier",
  "seed", "epsilon_ess",
  "epsilon_geweke_z", "epsilon_ess_below_100",
  "epsilon_abs_geweke_above_2", "stress_gene", "stress_list_name",
  "good_mcmc", "analysis_eligible"
)
master <- master[column_order]

# Write the canonical compressed table and verify that it can be read back.
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
connection <- gzfile(output_file, open = "wt", encoding = "UTF-8")
tryCatch(
  write.csv(master, connection, row.names = FALSE, na = "NA"),
  finally = close(connection)
)

written <- read.csv(gzfile(output_file), stringsAsFactors = FALSE)
if (nrow(written) != nrow(master) || !identical(names(written), names(master))) {
  stop("Saved master table failed validation")
}

cat(
  "Saved", format(nrow(master), big.mark = ","), "rows from",
  length(unique(master$cell_type)), "cell types to", output_file, "\n"
)
cat(
  "Stress list:", nrow(stress), "genes;", length(stress_overlap),
  "appear in BASiCS\n"
)
cat("LVG calls:", sum(master$LVG %in% TRUE), "\n")
cat("Analysis-eligible rows:", sum(master$analysis_eligible), "\n")
