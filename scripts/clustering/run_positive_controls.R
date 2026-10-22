#!/usr/bin/env Rscript

# Positive controls for the clustering pipeline. ASEL/ASER and AWC_ON/AWC_OFF
# are known, distinct identities within one anatomical neuron pair. Each pair is
# pooled into one pseudo-cell-type and run through the same steps and
# thresholds as scripts/clustering/run_clustering.R. The known identities are
# used only at the end, to score the discovered clusters (ARI).
#
# Output: one RDS with, per pair, the selected G, its labels, the known
# identities, the ARI, the resolution ladder, and the selected PC scores.

script_file_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
project_root <- dirname(dirname(dirname(normalizePath(sub("^--file=", "", script_file_argument[[1L]])))))
source(file.path(project_root, "src", "R", "project_io.R"))

options <- parse_named_arguments(
  commandArgs(trailingOnly = TRUE),
  list(
    sce = NA_character_,
    output = file.path(project_root, "outputs", "clustering", "positive_controls.rds"),
    exclude_transgenes = "true",
    seed = "20260809"
  ),
  paste(
    "Usage: Rscript scripts/clustering/run_positive_controls.R --sce PATH",
    "[--output FILE.rds] [--exclude-transgenes true|false] [--seed INTEGER] [--config PATH]"
  )
)
sce_path <- normalizePath(require_path_option(options$sce, "--sce"), mustWork = TRUE)
exclude_transgenes <- tolower(as.character(options$exclude_transgenes)) %in% c("true", "1", "yes")
seed <- as.integer(options$seed)

suppressPackageStartupMessages(library(SingleCellExperiment))
source(file.path(project_root, "src", "R", "clustering.R"))
source(file.path(project_root, "src", "R", "stress_gene_sets.R"))
source(file.path(project_root, "src", "R", "transgene_genes.R"))
source(file.path(project_root, "src", "R", "empty_cluster_guard.R"))

pairs <- list(
  ASEL_ASER = c("ASEL", "ASER"),
  AWC_pair = c("AWC_ON", "AWC_OFF")
)

sce_full <- readRDS(sce_path)

# A merged pseudo-cell-type column lets select_cell_type() run unchanged.
true_identity_full <- as.character(colData(sce_full)$Cell.type)
merged_label <- rep(NA_character_, ncol(sce_full))
for (pair_name in names(pairs)) {
  merged_label[true_identity_full %in% pairs[[pair_name]]] <- pair_name
}
colData(sce_full)$Cell.type.merged <- merged_label

# Same steps and thresholds as prepare_cell_type_for_gmm() and the stability
# ladder in run_clustering.R.
run_control <- function(pair_name) {
  cat("\n---- Positive control:", pair_name, "----\n")
  selected <- select_cell_type(sce_full, pair_name, cell_type_col = "Cell.type.merged")
  selected <- filter_rare_experiments(selected, min_cells = 5L)$sce
  true_identity <- setNames(as.character(colData(selected)$Cell.type), colnames(selected))

  counts <- assay(selected, "counts")
  cell_names <- colnames(counts)
  known_group <- setNames(as.character(colData(selected)$Experiment), cell_names)
  size_factor <- setNames(as.numeric(colData(selected)$Size_Factor), cell_names)
  normalized <- normalize_counts(counts, size_factor)

  stress_selection <- select_stress_features(counts, stress_genes_full_199, min_total_count = 10L)
  stress <- estimate_stress_score(normalized, stress_selection$usable, make_group_design(known_group))
  stress_mixture <- fit_stress_mixture(stress$score)
  nuisance_design <- make_nuisance_design(
    known_group = known_group,
    stress_score = stress$score,
    size_factor = size_factor,
    stress_posterior = stress_mixture$posterior,
    observation_names = cell_names
  )

  excluded_features <- stress_selection$matched
  if (exclude_transgenes) excluded_features <- union(excluded_features, transgene_ids)
  candidate_features <- select_candidate_features(
    counts,
    excluded_features = excluded_features,
    min_total_count = 10L
  )
  candidate_residuals <- residualize_features(
    normalized[candidate_features$eligible, , drop = FALSE],
    nuisance_design
  )
  high_variance <- select_high_variance_features(candidate_residuals, max_features = 3000L)
  pca <- compute_residual_pca(high_variance$residuals, max_pcs = 20L)
  pc_selection <- select_clustering_pcs(pca)

  gmm_candidates <- fit_gmm_candidates(
    pc_selection$scores,
    max_clusters = 10L,
    min_cluster_n = 10L,
    min_cluster_prop = 0.05
  )
  ladder <- build_resolution_ladder(
    scores = pc_selection$scores,
    gmm_candidates = gmm_candidates,
    known_group = known_group,
    stress_score = stress$score,
    size_factor = size_factor,
    gating_categorical = list(Detection = colData(selected)$Detection),
    gating_continuous = list(total_features_by_counts = colData(selected)$total_features_by_counts),
    diagnostic_continuous = list(pct_counts_Mito = colData(selected)$pct_counts_Mito),
    seed = seed
  )

  primary_labels <- ladder$primary_labels
  # A one-cluster selection has no partition to compare (ARI 0 by definition).
  ari <- if (ladder$primary_g > 1L) {
    mclust::adjustedRandIndex(true_identity[names(primary_labels)], primary_labels)
  } else {
    0
  }
  print(ladder$diagnostics[, intersect(
    c("g", "bic_delta_from_one", "stability_median", "status"),
    names(ladder$diagnostics)
  )])
  cat("Cells:", length(cell_names), " selected G:", ladder$primary_g, " ARI:", round(ari, 3), "\n")

  list(
    pair = pair_name,
    n = length(cell_names),
    primary_g = ladder$primary_g,
    ari = ari,
    primary_labels = primary_labels,
    true_identity = true_identity,
    ladder = ladder,
    # PC signs are arbitrary (irlba starts from a random vector).
    scores = pc_selection$scores,
    pca_variance_fraction = pca$variance_fraction_total
  )
}

results <- list()
for (pair_name in names(pairs)) {
  results[[paste0("gmm_", pair_name)]] <- run_control(pair_name)
}

dir.create(dirname(options$output), recursive = TRUE, showWarnings = FALSE)
saveRDS(results, options$output)
cat("\nSaved:", options$output, "\n")
