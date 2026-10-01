#!/usr/bin/env Rscript

# Run the hand-written clustering pipeline on one or every cell type.
#
# Examples:
#   Rscript scripts/clustering/run_clustering.R --sce data/raw/L4_neuron_sce.rds --run-id full
#   Rscript scripts/clustering/run_clustering.R --sce data/raw/L4_neuron_sce.rds --cell-type SMD --run-id smoke
#
# The runner has two stages:
#   1. A relatively cheap pass through PCA and the GMM fits.
#   2. Stability, nuisance diagnostics, attribution, and UMAP only when the
#      cheap pass finds a structurally supported multi-cluster resolution.

script_argument <- grep(
  "^--file=",
  commandArgs(trailingOnly = FALSE),
  value = TRUE
)

script_path <- if (length(script_argument)) {
  normalizePath(sub("^--file=", "", script_argument[1L]))
} else {
  normalizePath(file.path("scripts", "clustering", "run_clustering.R"))
}

project_root <- dirname(dirname(dirname(script_path)))

source(file.path(project_root, "src", "R", "clustering.R"))
source(file.path(project_root, "src", "R", "clustering_visualization.R"))
source(file.path(project_root, "src", "R", "stress_gene_sets.R"))

required_packages <- c(
  "SingleCellExperiment",
  "Matrix",
  "irlba",
  "mclust",
  "yaml"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1L), quietly = TRUE)
]

if (length(missing_packages)) {
  stop(
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", ")
  )
}

# Clustering_pipeline.R intentionally uses the familiar unqualified SCE
# accessors (colData, assay, and rowData), so attach their package here.
suppressPackageStartupMessages(library(SingleCellExperiment))

runner_defaults <- list(
  sce_path = NA_character_,
  output_root = file.path(project_root, "outputs", "clustering"),
  run_id = "default",
  output_directory = NULL,
  cell_type = "all",
  screen_only = FALSE,
  overwrite = FALSE,
  seed = 20260809L,
  min_experiment_cells = 5L,
  max_features = 3000L,
  max_pcs = 20L,
  max_clusters = 10L,
  bic_threshold = 10,
  min_cluster_n = 10L,
  min_cluster_prop = 0.05,
  n_stability = 50L,
  subsample_fraction = 0.80,
  stability_threshold = 0.75,
  cramer_threshold = 0.20,
  eta_threshold = 0.05,
  size_factor_eta_threshold = 0.05,
  min_total_count = 10L
)

print_usage <- function() {
  cat(
    paste(
      "Usage:",
      "  Rscript scripts/clustering/run_clustering.R --sce PATH [options]",
      "",
      "Options:",
      "  --cell-type VALUE   all, one cell type, or comma-separated types",
      "  --sce PATH          SingleCellExperiment RDS path",
      "  --config PATH       optional YAML configuration",
      "  --output PATH       explicit run output directory",
      "  --output-root PATH  root used with --run-id",
      "  --run-id VALUE      run directory name (default: default)",
      "  --screen-only       stop after the cheap PCA/GMM pass",
      "  --overwrite         recompute terminal results",
      "  --seed INTEGER      random seed",
      "  --stability INTEGER stability repetitions (default 50)",
      "  --help              print this help",
      sep = "\n"
    ),
    "\n"
  )
}

parse_runner_arguments <- function(arguments, defaults = runner_defaults) {
  config <- defaults
  config_positions <- which(arguments == "--config")
  if (length(config_positions)) {
    if (length(config_positions) != 1L || config_positions == length(arguments)) {
      stop("--config requires exactly one path")
    }
    configured <- yaml::read_yaml(arguments[config_positions + 1L])
    unknown <- setdiff(names(configured), names(config))
    if (length(unknown)) stop("Unknown configuration key(s): ", paste(unknown, collapse = ", "))
    config[names(configured)] <- configured
  }
  index <- 1L

  require_value <- function(option) {
    if (index == length(arguments)) {
      stop(option, " requires a value")
    }
    arguments[index + 1L]
  }

  while (index <= length(arguments)) {
    argument <- arguments[index]

    if (argument == "--help") {
      print_usage()
      quit(save = "no", status = 0L)
    } else if (argument == "--config") {
      index <- index + 2L
    } else if (argument == "--cell-type") {
      config$cell_type <- require_value(argument)
      index <- index + 2L
    } else if (argument == "--sce") {
      config$sce_path <- require_value(argument)
      index <- index + 2L
    } else if (argument == "--output") {
      config$output_directory <- require_value(argument)
      index <- index + 2L
    } else if (argument == "--output-root") {
      config$output_root <- require_value(argument)
      index <- index + 2L
    } else if (argument == "--run-id") {
      config$run_id <- require_value(argument)
      index <- index + 2L
    } else if (argument == "--seed") {
      config$seed <- as.integer(require_value(argument))
      index <- index + 2L
    } else if (argument == "--stability") {
      config$n_stability <- as.integer(require_value(argument))
      index <- index + 2L
    } else if (argument == "--screen-only") {
      config$screen_only <- TRUE
      index <- index + 1L
    } else if (argument == "--overwrite") {
      config$overwrite <- TRUE
      index <- index + 1L
    } else {
      stop("Unknown argument: ", argument)
    }
  }

  if (is.na(config$seed)) {
    stop("--seed must be an integer")
  }
  if (is.na(config$n_stability) || config$n_stability < 1L) {
    stop("--stability must be a positive integer")
  }
  if (is.null(config$output_directory) || !nzchar(config$output_directory)) {
    if (!grepl("^[A-Za-z0-9_.-]+$", config$run_id)) {
      stop("--run-id may contain only letters, numbers, dot, underscore, and hyphen")
    }
    config$output_directory <- file.path(config$output_root, config$run_id)
  }

  config
}

safe_file_label <- function(value) {
  label <- gsub("[^A-Za-z0-9_.-]", "_", value)
  if (!nzchar(label)) stop("Cell type cannot be converted to a file label")
  label
}

validate_input_sce <- function(sce) {
  if (!inherits(sce, "SingleCellExperiment")) {
    stop("The RDS must contain a SingleCellExperiment")
  }

  required_metadata <- c(
    "Cell.type",
    "Experiment",
    "Size_Factor",
    "Detection",
    "pct_counts_Mito",
    "total_features_by_counts"
  )

  missing_metadata <- setdiff(
    required_metadata,
    colnames(SummarizedExperiment::colData(sce))
  )

  if (length(missing_metadata)) {
    stop(
      "Missing required cell metadata: ",
      paste(missing_metadata, collapse = ", ")
    )
  }

  if (!"counts" %in% SummarizedExperiment::assayNames(sce)) {
    stop("The SingleCellExperiment does not contain a 'counts' assay")
  }

  counts <- SummarizedExperiment::assay(sce, "counts")
  if (is.null(rownames(counts)) || is.null(colnames(counts))) {
    stop("The count matrix must have feature and cell names")
  }
  if (anyDuplicated(rownames(counts)) || anyDuplicated(colnames(counts))) {
    stop("The count matrix must have unique feature and cell names")
  }

  invisible(sce)
}

choose_cell_types <- function(sce, requested) {
  available <- sort(unique(as.character(
    SummarizedExperiment::colData(sce)$Cell.type
  )))
  available <- available[!is.na(available) & nzchar(available)]

  if (identical(tolower(requested), "all")) {
    return(available)
  }

  selected <- unique(trimws(strsplit(requested, ",", fixed = TRUE)[[1L]]))
  selected <- selected[nzchar(selected)]
  missing <- setdiff(selected, available)

  if (length(missing)) {
    stop("Unknown cell type(s): ", paste(missing, collapse = ", "))
  }

  selected
}

prepare_cell_type_for_gmm <- function(
    sce,
    cell_type,
    stress_genes,
    config
) {
  selected <- select_cell_type(sce, cell_type)
  original_cell_count <- ncol(selected)

  experiment_filter <- filter_rare_experiments(
    selected,
    experiment_col = "Experiment",
    min_cells = config$min_experiment_cells
  )
  selected <- experiment_filter$sce

  if (ncol(selected) < 3L) {
    stop("Fewer than three cells remain after experiment filtering")
  }

  counts <- SummarizedExperiment::assay(selected, "counts")
  cell_names <- colnames(counts)

  known_group <- setNames(
    as.character(SummarizedExperiment::colData(selected)$Experiment),
    cell_names
  )
  size_factor <- setNames(
    as.numeric(SummarizedExperiment::colData(selected)$Size_Factor),
    cell_names
  )

  normalized <- normalize_counts(counts, size_factor)

  stress_selection <- select_stress_features(
    counts = counts,
    stress_genes = stress_genes,
    min_total_count = config$min_total_count
  )

  experiment_design <- make_group_design(known_group)
  stress <- estimate_stress_score(
    normalized = normalized,
    stress_features = stress_selection$usable,
    experiment_design = experiment_design
  )
  stress_mixture <- fit_stress_mixture(stress$score)

  nuisance_design <- make_nuisance_design(
    known_group = known_group,
    stress_score = stress$score,
    size_factor = size_factor,
    stress_posterior = stress_mixture$posterior,
    observation_names = cell_names
  )

  # Exclude every matched stress gene, including stress genes that were too
  # rare to enter the score, so clustering cannot rediscover the stress list.
  candidate_features <- select_candidate_features(
    counts = counts,
    excluded_features = stress_selection$matched,
    min_total_count = config$min_total_count
  )

  candidate_residuals <- residualize_features(
    expression = normalized[candidate_features$eligible, , drop = FALSE],
    design = nuisance_design
  )

  high_variance <- select_high_variance_features(
    residuals = candidate_residuals,
    max_features = config$max_features
  )

  pca <- compute_residual_pca(
    residuals = high_variance$residuals,
    max_pcs = config$max_pcs
  )
  pc_selection <- select_clustering_pcs(pca)

  gmm_candidates <- fit_gmm_candidates(
    scores = pc_selection$scores,
    max_clusters = config$max_clusters,
    min_cluster_n = config$min_cluster_n,
    min_cluster_prop = config$min_cluster_prop
  )

  gmm_structure <- evaluate_gmm_structure(
    gmm_candidates = gmm_candidates,
    observation_names = cell_names,
    bic_threshold = config$bic_threshold
  )

  supported_rows <- which(
    gmm_structure$diagnostics$g > 1L &
      gmm_structure$diagnostics$structurally_supported
  )

  screen_g <- if (!length(supported_rows)) {
    1L
  } else {
    eligible_table <- gmm_structure$diagnostics[supported_rows, , drop = FALSE]
    as.integer(eligible_table$g[which.max(eligible_table$bic)])
  }

  metadata <- as.data.frame(
    SummarizedExperiment::colData(selected),
    stringsAsFactors = FALSE
  )
  rownames(metadata) <- cell_names

  row_data <- SummarizedExperiment::rowData(selected)
  gene_annotations <- data.frame(
    feature = rownames(selected),
    gene_short_name = if ("gene_short_name" %in% colnames(row_data)) {
      as.character(row_data$gene_short_name)
    } else {
      rownames(selected)
    },
    stringsAsFactors = FALSE
  )

  list(
    cell_type = cell_type,
    counts = counts,
    normalized = normalized,
    known_group = known_group,
    size_factor = size_factor,
    metadata = metadata,
    gene_annotations = gene_annotations,
    original_cell_count = original_cell_count,
    experiment_filter = experiment_filter,
    stress_selection = stress_selection,
    stress = stress,
    stress_mixture = stress_mixture,
    nuisance_design = nuisance_design,
    candidate_features = candidate_features,
    candidate_residuals = candidate_residuals,
    high_variance = high_variance,
    pca = pca,
    pc_selection = pc_selection,
    gmm_candidates = gmm_candidates,
    gmm_structure = gmm_structure,
    screen_g = screen_g
  )
}

compact_screen_result <- function(state) {
  list(
    cell_type = state$cell_type,
    original_cells = state$original_cell_count,
    retained_cells = ncol(state$counts),
    excluded_cells = state$experiment_filter$excluded_cells,
    experiment_counts = state$experiment_filter$experiment_counts,
    retained_experiments = state$experiment_filter$retained_experiments,
    stress_features_matched = state$stress_selection$matched,
    stress_features_used = state$stress$features_used,
    stress_mixture_accepted = state$stress_mixture$accepted,
    candidate_feature_count = length(state$candidate_features$eligible),
    clustering_features = state$high_variance$features,
    selected_pc_count = state$pc_selection$number_selected,
    pc_diagnostics = state$pc_selection$diagnostics,
    gmm_diagnostics = state$gmm_structure$diagnostics,
    screen_g = state$screen_g
  )
}

add_gene_annotations <- function(table, gene_annotations) {
  if (is.null(table) || !nrow(table)) return(table)

  symbol <- as.character(gene_annotations$gene_short_name)
  names(symbol) <- gene_annotations$feature

  table$gene_short_name <- unname(symbol[table$feature])
  table
}

finish_cell_type <- function(state, config) {
  ladder <- build_resolution_ladder(
    scores = state$pc_selection$scores,
    gmm_candidates = state$gmm_candidates,
    known_group = state$known_group,
    stress_score = state$stress$score,
    size_factor = state$size_factor,
    gating_categorical = list(
      Detection = state$metadata$Detection
    ),
    gating_continuous = list(
      total_features_by_counts = state$metadata$total_features_by_counts
    ),
    # Mitochondrial percentage is reported but cannot veto a resolution.
    diagnostic_continuous = list(
      pct_counts_Mito = state$metadata$pct_counts_Mito
    ),
    bic_threshold = config$bic_threshold,
    stability_threshold = config$stability_threshold,
    n_stability = config$n_stability,
    subsample_fraction = config$subsample_fraction,
    cramer_threshold = config$cramer_threshold,
    eta_threshold = config$eta_threshold,
    size_factor_eta_threshold = config$size_factor_eta_threshold,
    seed = config$seed
  )

  result <- list(
    cell_type = state$cell_type,
    screen = compact_screen_result(state),
    resolution_ladder = ladder,
    stress = state$stress,
    stress_mixture = state$stress_mixture,
    candidate_features = state$candidate_features,
    high_variance_features = state$high_variance$features,
    pca = state$pca,
    pc_selection = state$pc_selection,
    nuisance_design = state$nuisance_design,
    assignments = ladder$assignments,
    primary_g = ladder$primary_g,
    primary_labels = ladder$primary_labels,
    attribution = NULL,
    visualization = NULL
  )

  if (ladder$primary_g <= 1L) {
    return(result)
  }

  labels <- ladder$primary_labels

  drivers <- rank_geometric_drivers(
    residuals = state$high_variance$residuals,
    cluster_labels = labels
  )

  markers <- rank_cluster_markers(
    counts = state$counts,
    normalized = state$normalized,
    features = state$candidate_features$eligible,
    cluster_labels = labels,
    nuisance_design = state$nuisance_design,
    clustering_features = state$high_variance$features,
    verbose = TRUE
  )

  summaries <- summarize_features_by_cluster(
    counts = state$counts,
    normalized = state$normalized,
    residuals = state$candidate_residuals,
    cluster_labels = labels
  )

  result$attribution <- list(
    drivers = add_gene_annotations(drivers, state$gene_annotations),
    markers = add_gene_annotations(markers, state$gene_annotations),
    cluster_summaries = add_gene_annotations(
      summaries,
      state$gene_annotations
    )
  )

  result$visualization <- prepare_cluster_visualization(
    counts = state$counts,
    size_factor = state$size_factor,
    pc_scores = state$pc_selection$scores,
    resolution_ladder = ladder,
    cell_metadata = state$metadata,
    gene_annotations = state$gene_annotations,
    gene_ranking = markers,
    seed = config$seed
  )

  result
}

write_cell_type_outputs <- function(
    directory,
    screen,
    result = NULL,
    status
) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)

  saveRDS(screen, file.path(directory, "screen_result.rds"))
  utils::write.csv(
    screen$gmm_diagnostics,
    file.path(directory, "gmm_screen.csv"),
    row.names = FALSE
  )

  if (!is.null(result)) {
    # Keep the sparse count matrix in the standalone viewer bundle instead of
    # duplicating it inside result.rds.
    visualization <- result$visualization
    result_to_save <- result
    if (!is.null(visualization)) {
      result_to_save$visualization <- list(
        bundle_file = "visualization_bundle.rds",
        settings = visualization$settings
      )
    }

    saveRDS(result_to_save, file.path(directory, "result.rds"))
    utils::write.csv(
      result$resolution_ladder$diagnostics,
      file.path(directory, "resolution_diagnostics.csv"),
      row.names = FALSE
    )
    utils::write.csv(
      result$assignments,
      file.path(directory, "assignments_by_resolution.csv"),
      row.names = FALSE
    )

    if (!is.null(result$attribution)) {
      utils::write.csv(
        result$attribution$drivers,
        file.path(directory, "primary_drivers.csv"),
        row.names = FALSE
      )
      utils::write.csv(
        result$attribution$markers,
        file.path(directory, "primary_markers.csv"),
        row.names = FALSE
      )
      utils::write.csv(
        result$attribution$cluster_summaries,
        file.path(directory, "primary_feature_summaries.csv"),
        row.names = FALSE
      )
    }

    if (!is.null(visualization)) {
      save_cluster_visualization_bundle(
        visualization,
        file.path(directory, "visualization_bundle.rds")
      )
    }
  }

  saveRDS(status, file.path(directory, "status.rds"))
  invisible(directory)
}

screen_one_cell_type <- function(sce, cell_type, config) {
  output_directory <- file.path(
    config$output_directory,
    safe_file_label(cell_type)
  )
  status_path <- file.path(output_directory, "status.rds")
  state_path <- file.path(output_directory, "screen_state.rds")
  screen_path <- file.path(output_directory, "screen_result.rds")

  if (file.exists(status_path) && !isTRUE(config$overwrite)) {
    previous <- readRDS(status_path)
    completed_statuses <- c(
      "screen_homogeneous",
      "accepted_multicluster",
      "rejected_after_diagnostics"
    )

    if (previous$status %in% completed_statuses) {
      message("[", cell_type, "] cheap pass already complete")
      return(previous)
    }

    resumable <- file.exists(state_path) && file.exists(screen_path) &&
      previous$status %in% c(
        "screen_candidate",
        "in_depth_pending",
        "in_depth_error"
      )

    if (resumable) {
      screen <- readRDS(screen_path)
      previous$status <- "screen_candidate"
      previous$screen_g <- screen$screen_g
      previous$error <- NA_character_
      saveRDS(previous, status_path)
      message("[", cell_type, "] reusing saved GMM candidate state")
      return(previous)
    }
  }

  message("[", cell_type, "] cheap PCA/GMM pass")
  started <- Sys.time()

  state <- prepare_cell_type_for_gmm(
    sce = sce,
    cell_type = cell_type,
    stress_genes = stress_genes_full_199,
    config = config
  )
  screen <- compact_screen_result(state)

  status <- list(
    cell_type = cell_type,
    status = if (state$screen_g > 1L) {
      "screen_candidate"
    } else {
      "screen_homogeneous"
    },
    retained_cells = ncol(state$counts),
    screen_g = state$screen_g,
    primary_g = NA_integer_,
    elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
    error = NA_character_
  )

  write_cell_type_outputs(output_directory, screen, status = status)

  if (state$screen_g > 1L) {
    # Candidate states are kept only between the two phases. This avoids
    # holding dense residual matrices for several cell types in memory.
    saveRDS(state, state_path, compress = FALSE)
  }

  status
}

finish_screened_cell_type <- function(cell_type, config) {
  output_directory <- file.path(
    config$output_directory,
    safe_file_label(cell_type)
  )
  status_path <- file.path(output_directory, "status.rds")
  state_path <- file.path(output_directory, "screen_state.rds")
  screen_path <- file.path(output_directory, "screen_result.rds")

  if (!file.exists(status_path)) {
    stop("No cheap-pass status exists for ", cell_type)
  }

  previous <- readRDS(status_path)
  if (previous$status %in% c(
    "screen_homogeneous",
    "accepted_multicluster",
    "rejected_after_diagnostics"
  )) {
    return(previous)
  }

  if (!file.exists(state_path) || !file.exists(screen_path)) {
    stop("Saved GMM candidate state is missing for ", cell_type)
  }

  state <- readRDS(state_path)
  screen <- readRDS(screen_path)
  started <- Sys.time()

  pending_status <- list(
    cell_type = cell_type,
    status = "in_depth_pending",
    retained_cells = ncol(state$counts),
    screen_g = state$screen_g,
    primary_g = NA_integer_,
    elapsed_seconds = previous$elapsed_seconds,
    error = NA_character_
  )
  write_cell_type_outputs(
    output_directory,
    screen = screen,
    status = pending_status
  )

  message("[", cell_type, "] in-depth stability and nuisance diagnostics")
  result <- finish_cell_type(state, config)

  total_elapsed <- as.numeric(previous$elapsed_seconds) +
    as.numeric(difftime(Sys.time(), started, units = "secs"))

  status <- list(
    cell_type = cell_type,
    status = if (result$primary_g > 1L) {
      "accepted_multicluster"
    } else {
      "rejected_after_diagnostics"
    },
    retained_cells = ncol(state$counts),
    screen_g = state$screen_g,
    primary_g = result$primary_g,
    elapsed_seconds = total_elapsed,
    error = NA_character_
  )

  write_cell_type_outputs(
    output_directory,
    screen = screen,
    result = result,
    status = status
  )

  # This is an intermediate checkpoint created by this runner. The durable
  # compact screen, final result, tables, and viewer bundle remain.
  unlink(state_path)
  status
}

status_row <- function(status) {
  data.frame(
    cell_type = as.character(status$cell_type),
    status = as.character(status$status),
    retained_cells = as.integer(status$retained_cells),
    screen_g = as.integer(status$screen_g),
    primary_g = as.integer(status$primary_g),
    elapsed_seconds = as.numeric(status$elapsed_seconds),
    error = as.character(status$error),
    stringsAsFactors = FALSE
  )
}

write_run_summary <- function(statuses, output_directory) {
  summary_table <- do.call(rbind, lapply(statuses, status_row))
  utils::write.csv(
    summary_table,
    file.path(output_directory, "run_summary.csv"),
    row.names = FALSE
  )
  summary_table
}

error_status <- function(cell_type, error, stage, config) {
  message("[", cell_type, "] ERROR: ", conditionMessage(error))
  directory <- file.path(
    config$output_directory,
    safe_file_label(cell_type)
  )
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)

  previous_path <- file.path(directory, "status.rds")
  previous <- if (file.exists(previous_path)) {
    readRDS(previous_path)
  } else {
    list(
      retained_cells = NA_integer_,
      screen_g = NA_integer_,
      primary_g = NA_integer_,
      elapsed_seconds = NA_real_
    )
  }

  status <- list(
    cell_type = cell_type,
    status = if (stage == "in_depth") "in_depth_error" else "screen_error",
    retained_cells = previous$retained_cells,
    screen_g = previous$screen_g,
    primary_g = previous$primary_g,
    elapsed_seconds = previous$elapsed_seconds,
    error = conditionMessage(error)
  )
  saveRDS(status, previous_path)
  status
}

run_clustering_pipeline <- function(config = runner_defaults) {
  if (is.na(config$sce_path) || !nzchar(config$sce_path) || !file.exists(config$sce_path)) {
    stop("SCE file not found: ", config$sce_path)
  }

  dir.create(config$output_directory, recursive = TRUE, showWarnings = FALSE)
  config$sce_path <- normalizePath(config$sce_path, mustWork = TRUE)
  config$input_md5 <- unname(tools::md5sum(config$sce_path))
  config$created_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  saveRDS(config, file.path(config$output_directory, "resolved_config.rds"))
  yaml::write_yaml(config, file.path(config$output_directory, "resolved_config.yaml"))
  writeLines(capture.output(sessionInfo()), file.path(config$output_directory, "session_info.txt"))

  message("Loading ", config$sce_path)
  sce <- readRDS(config$sce_path)
  validate_input_sce(sce)
  cell_types <- choose_cell_types(sce, config$cell_type)

  message("Selected ", length(cell_types), " cell type(s)")
  message("Phase 1: cheap PCA/GMM screen for every selected cell type")
  statuses <- vector("list", length(cell_types))

  for (index in seq_along(cell_types)) {
    cell_type <- cell_types[index]
    message("Cell type ", index, " of ", length(cell_types), ": ", cell_type)

    statuses[[index]] <- tryCatch(
      screen_one_cell_type(sce, cell_type, config),
      error = function(error) {
        error_status(cell_type, error, "screen", config)
      }
    )

    write_run_summary(
      statuses[seq_len(index)],
      config$output_directory
    )
  }

  if (isTRUE(config$screen_only)) {
    return(write_run_summary(statuses, config$output_directory))
  }

  candidate_indices <- which(
    vapply(
      statuses,
      function(status) identical(status$status, "screen_candidate"),
      logical(1L)
    )
  )

  message(
    "Phase 2: in-depth analysis for ",
    length(candidate_indices),
    " GMM candidate cell type(s)"
  )

  for (candidate_index in candidate_indices) {
    cell_type <- cell_types[candidate_index]
    statuses[[candidate_index]] <- tryCatch(
      finish_screened_cell_type(cell_type, config),
      error = function(error) {
        error_status(cell_type, error, "in_depth", config)
      }
    )
    write_run_summary(statuses, config$output_directory)
  }

  write_run_summary(statuses, config$output_directory)
}

if (sys.nframe() == 0L) {
  config <- parse_runner_arguments(commandArgs(trailingOnly = TRUE))
  summary <- run_clustering_pipeline(config)
  print(summary)
  cat("Outputs:", normalizePath(config$output_directory), "\n")
}
