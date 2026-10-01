#!/usr/bin/env Rscript

# Test WormCat categories and neuron families for HVG or LVG enrichment.
script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_file[[1L]]))
project_root <- dirname(dirname(dirname(script_path)))
source(file.path(project_root, "src", "R", "project_io.R"))
options <- parse_named_arguments(
  commandArgs(trailingOnly = TRUE),
  list(master = NA_character_, gene_classes = NA_character_, neuron_classes = NA_character_, output_dir = NA_character_, target = "HVG"),
  paste(
    "Usage: Rscript scripts/analysis/02_run_enrichment.R --master FILE.csv.gz",
    "--gene-classes FILE.csv --neuron-classes FILE.csv --output-dir PATH [--target HVG|LVG] [--config PATH]"
  )
)
master_file <- normalizePath(require_path_option(options$master, "--master"), mustWork = TRUE)
gene_file <- normalizePath(require_path_option(options$gene_classes, "--gene-classes"), mustWork = TRUE)
neuron_file <- normalizePath(require_path_option(options$neuron_classes, "--neuron-classes"), mustWork = TRUE)
output_dir <- require_path_option(options$output_dir, "--output-dir", must_exist = FALSE)
target <- toupper(options$target)
if (!target %in% c("HVG", "LVG")) stop("Target must be HVG or LVG")

# Stop when a table lacks columns required by the tests.
check_columns <- function(x, required, path) {
  missing <- setdiff(required, names(x))
  if (length(missing)) {
    stop("Missing columns in ", path, ": ", paste(missing, collapse = ", "))
  }
}

# Calculate a one-sided hypergeometric enrichment test.
test_enrichment <- function(k, K, n, N) {
  p_value <- phyper(k - 1L, K, N - K, n, lower.tail = FALSE)
  fold <- if (n > 0L && K > 0L) (k / n) / (K / N) else NA_real_
  c(fold_enrichment = fold, p_value = p_value)
}

master <- read.csv(gzfile(master_file), stringsAsFactors = FALSE)
genes <- read.csv(gene_file, stringsAsFactors = FALSE)
neurons <- read.csv(neuron_file, stringsAsFactors = FALSE)

check_columns(
  master,
  c("cell_type", "stable_id", "gene_name", target, "analysis_eligible"),
  master_file
)
check_columns(
  genes,
  c("stable_id", "category_1", "category_2", "category_3"),
  gene_file
)
check_columns(neurons, c("cell_type", "neuron_family"), neuron_file)

if (anyDuplicated(paste(master$cell_type, master$stable_id))) {
  stop("Master table contains duplicate cell-type/gene rows")
}
if (anyDuplicated(genes$stable_id)) stop("Gene mapping contains duplicate IDs")
if (anyDuplicated(neurons$cell_type)) stop("Neuron mapping contains duplicate types")

eligible <- master$analysis_eligible %in% TRUE
if (anyNA(master[[target]][eligible])) {
  stop("Eligible rows contain missing ", target, " calls")
}

# Test each WormCat level separately within every cell type.
gene_index <- match(master$stable_id, genes$stable_id)
levels <- c("category_1", "category_2", "category_3")
wormcat_results <- list()
row_id <- 1L

for (cell_type in sort(unique(master$cell_type))) {
  use <- eligible & master$cell_type == cell_type
  selected <- master[[target]][use] %in% TRUE
  N <- sum(use)
  n <- sum(selected)

  for (level in levels) {
    annotation <- genes[[level]][gene_index[use]]
    categories <- sort(unique(annotation[!is.na(annotation) & nzchar(annotation)]))

    for (category in categories) {
      member <- annotation == category & !is.na(annotation)
      K <- sum(member)
      k <- sum(member & selected)
      test <- test_enrichment(k, K, n, N)

      wormcat_results[[row_id]] <- data.frame(
        cell_type = cell_type,
        category_level = level,
        category = category,
        background_genes = N,
        category_genes = K,
        selected_genes = n,
        category_selected_genes = k,
        fold_enrichment = test[["fold_enrichment"]],
        p_value = test[["p_value"]],
        stringsAsFactors = FALSE
      )
      row_id <- row_id + 1L
    }
  }
}

wormcat_results <- do.call(rbind, wormcat_results)
group <- interaction(
  wormcat_results$cell_type,
  wormcat_results$category_level,
  drop = TRUE
)
wormcat_results$q_value <- ave(
  wormcat_results$p_value,
  group,
  FUN = function(x) p.adjust(x, method = "BH")
)
wormcat_results$significant <-
  wormcat_results$q_value <= 0.05 & wormcat_results$fold_enrichment > 1
wormcat_results <- wormcat_results[order(
  wormcat_results$q_value,
  -wormcat_results$fold_enrichment,
  wormcat_results$cell_type
), ]
rownames(wormcat_results) <- NULL

# Test whether each gene's selected calls concentrate in a neuron family.
family_index <- match(master$cell_type, neurons$cell_type)
family <- neurons$neuron_family[family_index]
use <- eligible & !is.na(family)
family_sizes <- table(neurons$neuron_family)
gene_rows <- split(which(use), master$stable_id[use])
neuron_results <- list()
row_id <- 1L

for (stable_id in names(gene_rows)) {
  rows <- gene_rows[[stable_id]]
  selected <- master[[target]][rows] %in% TRUE
  N <- length(rows)
  n <- sum(selected)

  for (neuron_family in names(family_sizes)) {
    member <- family[rows] == neuron_family
    K <- sum(member)
    if (K == 0L || K == N) next

    k <- sum(member & selected)
    test <- test_enrichment(k, K, n, N)
    neuron_results[[row_id]] <- data.frame(
      neuron_family = neuron_family,
      stable_id = stable_id,
      gene_name = master$gene_name[rows[1]],
      mapped_family_cell_types = as.integer(family_sizes[[neuron_family]]),
      eligible_cell_types = N,
      eligible_family_cell_types = K,
      selected_cell_types = n,
      selected_family_cell_types = k,
      fold_enrichment = test[["fold_enrichment"]],
      p_value = test[["p_value"]],
      stringsAsFactors = FALSE
    )
    row_id <- row_id + 1L
  }
}

neuron_results <- do.call(rbind, neuron_results)
neuron_results$q_value <- ave(
  neuron_results$p_value,
  neuron_results$neuron_family,
  FUN = function(x) p.adjust(x, method = "BH")
)
neuron_results$significant <-
  neuron_results$q_value <= 0.05 & neuron_results$fold_enrichment > 1
neuron_results <- neuron_results[order(
  neuron_results$q_value,
  -neuron_results$fold_enrichment,
  neuron_results$neuron_family,
  neuron_results$gene_name
), ]
rownames(neuron_results) <- NULL

# Test WormCat categories among distinct selected genes in each neuron family.
family_category_results <- list()
row_id <- 1L

for (neuron_family in names(family_sizes)) {
  rows <- which(use & family == neuron_family)
  selected_by_gene <- tapply(
    master[[target]][rows] %in% TRUE,
    master$stable_id[rows],
    any
  )
  stable_ids <- names(selected_by_gene)
  mapped_genes <- match(stable_ids, genes$stable_id)
  N <- length(stable_ids)
  n <- sum(selected_by_gene)

  for (level in levels) {
    annotation <- genes[[level]][mapped_genes]
    categories <- sort(unique(annotation[!is.na(annotation) & nzchar(annotation)]))

    for (category in categories) {
      member <- annotation == category & !is.na(annotation)
      K <- sum(member)
      k <- sum(member & selected_by_gene)
      test <- test_enrichment(k, K, n, N)

      family_category_results[[row_id]] <- data.frame(
        neuron_family = neuron_family,
        category_level = level,
        category = category,
        mapped_cell_types = as.integer(family_sizes[[neuron_family]]),
        background_genes = N,
        category_genes = K,
        selected_genes = n,
        category_selected_genes = k,
        fold_enrichment = test[["fold_enrichment"]],
        p_value = test[["p_value"]],
        stringsAsFactors = FALSE
      )
      row_id <- row_id + 1L
    }
  }
}

family_category_results <- do.call(rbind, family_category_results)
group <- interaction(
  family_category_results$neuron_family,
  family_category_results$category_level,
  drop = TRUE
)
family_category_results$q_value <- ave(
  family_category_results$p_value,
  group,
  FUN = function(x) p.adjust(x, method = "BH")
)
family_category_results$significant <-
  family_category_results$q_value <= 0.05 &
  family_category_results$fold_enrichment > 1
family_category_results <- family_category_results[order(
  family_category_results$q_value,
  -family_category_results$fold_enrichment,
  family_category_results$neuron_family
), ]
rownames(family_category_results) <- NULL

# Name count columns and output files for the selected target.
rename_target <- function(x) {
  names(x) <- sub("selected", tolower(target), names(x), fixed = TRUE)
  x
}
wormcat_output <- rename_target(wormcat_results)
neuron_output <- rename_target(neuron_results)
family_category_output <- rename_target(family_category_results)
target_name <- tolower(target)

# Save one complete table for each enrichment question.
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(
  wormcat_output,
  file.path(output_dir, paste0("wormcat_", target_name, "_enrichment.csv")),
  row.names = FALSE,
  na = "NA"
)
write.csv(
  neuron_output,
  file.path(output_dir, paste0("neuron_family_", target_name, "_enrichment.csv")),
  row.names = FALSE,
  na = "NA"
)
write.csv(
  family_category_output,
  file.path(
    output_dir,
    paste0("neuron_gene_family_", target_name, "_enrichment.csv")
  ),
  row.names = FALSE,
  na = "NA"
)

cat(
  target, "WormCat:", nrow(wormcat_results), "tests;",
  sum(wormcat_results$significant), "significant\n"
)
cat(
  target, "neuron families:", nrow(neuron_results), "tests;",
  sum(neuron_results$significant), "significant\n"
)
cat(
  target, "neuron x gene families:", nrow(family_category_results), "tests;",
  sum(family_category_results$significant), "significant\n"
)
