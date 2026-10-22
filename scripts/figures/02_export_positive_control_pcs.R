#!/usr/bin/env Rscript

# Export the selected PC scores, GMM clusters and known identities of the two
# pooled positive controls (Figures 9, 10, 11), from the RDS written by
# scripts/clustering/run_positive_controls.R.
#
# Writes OUTPUT_DIR/data/positive_ASE.csv, positive_AWC.csv and
#        positive_control_manifest.csv

script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_file[[1L]]))
project_root <- dirname(dirname(dirname(script_path)))
source(file.path(project_root, "src", "R", "project_io.R"))
source(file.path(project_root, "src", "R", "figures.R"))

options <- figure_options(
  "Usage: Rscript scripts/figures/02_export_positive_control_pcs.R --config config/figures.yaml"
)
positive_results <- readRDS(input_file(options$positive_controls))

# ASE G=3 was rejected only as qc-sensitive within the borderline margin and is
# accepted under the borderline rule (Methods); use that solution.
ase <- positive_results$gmm_ASEL_ASER
ase$primary_g <- 3L
ase$primary_labels <- setNames(
  ase$ladder$assignments$cluster_g3,
  ase$ladder$assignments$observation
)
ase$ari <- mclust::adjustedRandIndex(
  ase$primary_labels,
  ase$true_identity[names(ase$primary_labels)]
)
positive_results$gmm_ASEL_ASER <- ase

first_or_na <- function(values, index) {
  if (length(values) >= index) as.numeric(values[index]) else NA_real_
}

control_names <- c(ASE = "gmm_ASEL_ASER", AWC = "gmm_AWC_pair")
display_names <- c(ASE = "ASEL + ASER", AWC = "AWC_ON + AWC_OFF")
manifest_rows <- list()

for (pair_name in names(control_names)) {
  control <- positive_results[[control_names[[pair_name]]]]
  if (is.null(control$scores)) {
    stop("No PC scores in ", options$positive_controls,
         "; rerun scripts/clustering/run_positive_controls.R", call. = FALSE)
  }

  cells <- names(control$primary_labels)
  if (!setequal(cells, rownames(control$scores))) {
    stop(pair_name, ": cluster labels and PC scores cover different cells")
  }
  scores <- as.matrix(control$scores[cells, , drop = FALSE])
  colnames(scores) <- paste0("PC", seq_len(ncol(scores)))
  # irlba starts from a random vector, so each PC's sign is arbitrary between
  # runs. Fix it (long tail towards negative values) so the PCA panels match
  # the manuscript. UMAP and clustering are unaffected by the sign.
  scores <- sweep(scores, 2L, -sign(colSums(scores^3)), `*`)

  primary_g <- as.integer(control$primary_g)
  diagnostics <- control$ladder$diagnostics
  primary <- diagnostics[diagnostics$g == primary_g, , drop = FALSE]

  dataset_id <- paste0("positive_", pair_name)
  data_file <- paste0(dataset_id, ".csv")
  utils::write.csv(
    data.frame(
      dataset_id = dataset_id,
      display_name = display_names[[pair_name]],
      cell = cells,
      cluster = as.character(control$primary_labels[cells]),
      truth = as.character(control$true_identity[cells]),
      scores,
      check.names = FALSE,
      stringsAsFactors = FALSE
    ),
    file.path(options$data_dir, data_file),
    row.names = FALSE,
    na = ""
  )

  manifest_rows[[pair_name]] <- data.frame(
    dataset_id = dataset_id,
    display_name = display_names[[pair_name]],
    data_file = data_file,
    n_cells = nrow(scores),
    primary_g = primary_g,
    external_ari = as.numeric(control$ari),
    bic_delta_from_one = as.numeric(primary$bic_delta_from_one),
    stability_median = as.numeric(primary$stability_median),
    cluster_sizes = as.character(primary$cluster_sizes),
    selected_pc_count = ncol(scores),
    pc1_variance_fraction = first_or_na(control$pca_variance_fraction, 1L),
    pc2_variance_fraction = first_or_na(control$pca_variance_fraction, 2L),
    stringsAsFactors = FALSE
  )
}

manifest <- do.call(rbind, manifest_rows)
utils::write.csv(
  manifest,
  file.path(options$data_dir, "positive_control_manifest.csv"),
  row.names = FALSE,
  na = ""
)
print(manifest[, c("dataset_id", "n_cells", "primary_g", "external_ari", "stability_median")])
