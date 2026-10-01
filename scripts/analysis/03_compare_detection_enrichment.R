#!/usr/bin/env Rscript

# Compare WormCat HVG enrichment before and after multi-experiment detection.
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
  list(master = NA_character_, sce = NA_character_, run_dir = NA_character_, gene_classes = NA_character_, output_dir = NA_character_),
  paste(
    "Usage: Rscript scripts/analysis/03_compare_detection_enrichment.R --master FILE.csv.gz",
    "--sce FILE.rds --run-dir PATH --gene-classes FILE.csv",
    "--output-dir PATH [--config PATH]"
  )
)
master_file <- normalizePath(require_path_option(options$master, "--master"))
sce_file <- normalizePath(require_path_option(options$sce, "--sce"))
run_dir <- normalizePath(require_path_option(options$run_dir, "--run-dir"))
gene_class_file <- normalizePath(require_path_option(options$gene_classes, "--gene-classes"))
output_dir <- require_path_option(options$output_dir, "--output-dir", must_exist = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

master <- read.csv(gzfile(master_file), stringsAsFactors = FALSE)
tasks <- read.delim(file.path(run_dir, "cell_types.tsv"), stringsAsFactors = FALSE)
config <- readRDS(file.path(run_dir, "config.rds"))
if (unname(tools::md5sum(sce_file)) != config$input_md5) {
  stop("SCE checksum does not match the BASiCS run")
}

# Count experiment-level detections for every modeled gene.
sce <- readRDS(sce_file)
counts <- assay(sce, "counts")
cell_type <- as.character(colData(sce)$Cell.type)
experiment <- as.character(colData(sce)$Experiment)

measure_support <- function(task) {
  rows <- which(master$cell_type == task$cell_type)
  ids <- master$stable_id[rows]
  result_dir <- basics_result_dir(run_dir, task$task_id, task$cell_type)
  experiment_table <- read.csv(
    file.path(result_dir, "experiments.csv"), stringsAsFactors = FALSE
  )
  included <- experiment_table$experiment[experiment_table$included]
  detected <- vapply(included, function(batch) {
    cells <- which(cell_type == task$cell_type & experiment == batch)
    Matrix::rowSums(counts[ids, cells, drop = FALSE] > 0)
  }, numeric(length(ids)))
  if (is.null(dim(detected))) detected <- matrix(detected, ncol = 1L)
  experiment_cells <- vapply(
    included,
    function(batch) sum(cell_type == task$cell_type & experiment == batch),
    integer(1)
  )
  thresholds <- pmax(3L, ceiling(0.05 * experiment_cells))
  supported <- sweep(detected, 2L, thresholds, `>=`)
  total_detected <- rowSums(detected)

  data.frame(
    cell_type = task$cell_type,
    stable_id = ids,
    n_eligible_experiments = length(included),
    n_detection_supported_experiments = rowSums(supported),
    max_experiment_detection_share = apply(detected, 1L, max) / total_detected,
    robust_detection_support = rowSums(supported) >= 2L,
    stringsAsFactors = FALSE
  )
}

support <- do.call(rbind, lapply(seq_len(nrow(tasks)), function(i) {
  measure_support(tasks[i, ])
}))
support_key <- paste(support$cell_type, support$stable_id, sep = "::")
master_key <- paste(master$cell_type, master$stable_id, sep = "::")
index <- match(master_key, support_key)
if (anyNA(index) || anyDuplicated(support_key)) stop("Support rows do not match master")
master$robust_detection_support <- support$robust_detection_support[index]

# Convert the existing WormCat hierarchy into the same gene-term format.
gene_classes <- read.csv(gene_class_file, stringsAsFactors = FALSE)
wormcat_pairs <- do.call(rbind, lapply(
  c("category_1", "category_2", "category_3"),
  function(level) {
    keep <- !is.na(gene_classes[[level]]) & nzchar(gene_classes[[level]])
    data.frame(
      stable_id = gene_classes$stable_id[keep],
      term_id = paste(level, gene_classes[[level]][keep], sep = "::"),
      stringsAsFactors = FALSE
    )
  }
))
wormcat_pairs <- unique(wormcat_pairs)
pieces <- strsplit(unique(wormcat_pairs$term_id), "::", fixed = TRUE)
wormcat_terms <- data.frame(
  term_id = vapply(pieces, paste, character(1), collapse = "::"),
  term = vapply(pieces, `[`, character(1), 2L),
  ontology = vapply(pieces, `[`, character(1), 1L),
  stringsAsFactors = FALSE
)

# Run one-sided over-representation tests with a condition-matched universe.
run_enrichment <- function(pairs, term_table, condition, min_size = 5L,
                           max_size = 500L) {
  annotated <- unique(pairs$stable_id)
  output <- vector("list", nrow(tasks))

  for (i in seq_len(nrow(tasks))) {
    type <- tasks$cell_type[i]
    rows <- master$cell_type == type & master$analysis_eligible
    if (condition == "supported") rows <- rows & master$robust_detection_support
    background <- intersect(master$stable_id[rows], annotated)
    selected_total <- master$stable_id[rows & master$HVG]
    selected <- intersect(selected_total, background)
    local <- pairs[pairs$stable_id %in% background, ]
    K <- table(local$term_id)
    k <- table(local$term_id[local$stable_id %in% selected])
    ids <- names(K)[K >= min_size & K <= max_size]
    K <- as.integer(K[ids])
    k <- as.integer(k[ids]); k[is.na(k)] <- 0L
    N <- length(background)
    n <- length(selected)
    p <- phyper(k - 1L, K, N - K, n, lower.tail = FALSE)
    fold <- (k / n) / (K / N)
    annotation <- term_table[match(ids, term_table$term_id), ]
    output[[i]] <- data.frame(
      condition = condition,
      cell_type = type,
      term_id = ids,
      term = annotation$term,
      ontology = annotation$ontology,
      eligible_genes = sum(rows),
      annotated_background_genes = N,
      selected_genes_total = length(selected_total),
      annotated_selected_genes = n,
      term_genes = K,
      term_selected_genes = k,
      fold_enrichment = fold,
      p_value = p,
      stringsAsFactors = FALSE
    )
  }

  result <- do.call(rbind, output)
  group <- interaction(result$cell_type, result$ontology, drop = TRUE)
  result$q_value <- ave(result$p_value, group, FUN = p.adjust, method = "BH")
  result$significant <- result$q_value <= 0.05 &
    result$fold_enrichment > 1 & result$term_selected_genes >= 3L
  result[order(result$cell_type, result$q_value, -result$fold_enrichment), ]
}

# Compare term calls and effect sizes between the two universes.
compare_results <- function(full, supported) {
  columns <- c(
    "cell_type", "term_id", "term", "ontology", "fold_enrichment",
    "q_value", "significant", "term_selected_genes", "term_genes"
  )
  x <- merge(
    full[columns], supported[columns],
    by = c("cell_type", "term_id", "term", "ontology"),
    suffixes = c("_full", "_supported"), all = TRUE
  )
  x$significant_full[is.na(x$significant_full)] <- FALSE
  x$significant_supported[is.na(x$significant_supported)] <- FALSE
  x$status <- ifelse(
    x$significant_full & x$significant_supported, "retained",
    ifelse(x$significant_full, "lost", ifelse(x$significant_supported, "gained", "neither"))
  )
  x
}

# Summarize comparability separately for each cell type.
summarize_comparison <- function(comparison, full, supported) {
  do.call(rbind, lapply(sort(unique(master$cell_type)), function(type) {
    x <- comparison[comparison$cell_type == type, ]
    common <- is.finite(x$fold_enrichment_full) & is.finite(x$fold_enrichment_supported)
    full_sig <- x$significant_full
    supported_sig <- x$significant_supported
    union <- sum(full_sig | supported_sig)
    full_rows <- full[full$cell_type == type, ][1, ]
    supported_rows <- supported[supported$cell_type == type, ][1, ]
    data.frame(
      cell_type = type,
      full_hvgs = full_rows$selected_genes_total,
      supported_hvgs = supported_rows$selected_genes_total,
      hvg_retention = supported_rows$selected_genes_total / full_rows$selected_genes_total,
      significant_full = sum(full_sig),
      significant_supported = sum(supported_sig),
      significant_retained = sum(full_sig & supported_sig),
      significant_lost = sum(full_sig & !supported_sig),
      significant_gained = sum(!full_sig & supported_sig),
      significant_jaccard = if (union) sum(full_sig & supported_sig) / union else NA_real_,
      fold_spearman = if (sum(common) >= 3L) {
        cor(
          x$fold_enrichment_full[common],
          x$fold_enrichment_supported[common],
          method = "spearman"
        )
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  }))
}

wormcat_full <- run_enrichment(wormcat_pairs, wormcat_terms, "full")
wormcat_supported <- run_enrichment(wormcat_pairs, wormcat_terms, "supported")
wormcat_comparison <- compare_results(wormcat_full, wormcat_supported)
wormcat_summary <- summarize_comparison(
  wormcat_comparison, wormcat_full, wormcat_supported
)

# Write compressed detailed tables and small uncompressed summaries.
write_gz <- function(x, name) {
  connection <- gzfile(file.path(output_dir, name), "wt")
  on.exit(close(connection))
  write.csv(x, connection, row.names = FALSE, na = "NA")
}
support <- support[order(support$cell_type, support$stable_id), ]
write_gz(support, "gene_experiment_support.csv.gz")
write_gz(wormcat_full, "wormcat_full.csv.gz")
write_gz(wormcat_supported, "wormcat_supported.csv.gz")
write_gz(wormcat_comparison, "wormcat_comparison.csv.gz")
write.csv(
  wormcat_summary,
  file.path(output_dir, "wormcat_celltype_summary.csv"),
  row.names = FALSE
)

manifest <- c(
  paste("created_at", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), sep = "\t"),
  paste("master", master_file, sep = "\t"),
  paste("gene_classes_md5", unname(tools::md5sum(gene_class_file)), sep = "\t"),
  paste("sce_md5", unname(tools::md5sum(sce_file)), sep = "\t"),
  paste("support_rule", "detected in max(3 cells, 5% of cells) in >=2 experiments", sep = "\t"),
  paste("stress_filter", "CeNGen 199-gene list", sep = "\t"),
  paste("term_size", "5 to 500 annotated background genes", sep = "\t"),
  paste("significance", "BH q<=0.05, fold>1, >=3 selected genes", sep = "\t")
)
writeLines(manifest, file.path(output_dir, "manifest.txt"))

cat(sprintf(
  "WormCat significant terms: %d full, %d supported\n",
  sum(wormcat_full$significant), sum(wormcat_supported$significant)
))
