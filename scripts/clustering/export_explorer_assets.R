#!/usr/bin/env Rscript

# Export accepted clustering results to small, Python-friendly Parquet assets.
#
# The exporter deliberately keeps only expression for supported marker genes.
# It never writes the complete count matrix. Each cell type receives a compact
# set of tables that the BASiCS Streamlit application can read independently.

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_argument)) {
  normalizePath(sub("^--file=", "", script_argument[1L]))
} else {
  normalizePath(file.path("scripts", "clustering", "export_explorer_assets.R"))
}
project_root <- dirname(dirname(dirname(script_path)))

require_export_package <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Package '", package, "' is required", call. = FALSE)
  }
}

parse_export_arguments <- function(arguments = commandArgs(trailingOnly = TRUE)) {
  defaults <- list(
    input_dir = NA_character_,
    output_dir = file.path(project_root, "outputs", "explorer", "clustering")
  )

  config_positions <- which(arguments == "--config")
  if (length(config_positions)) {
    if (length(config_positions) != 1L || config_positions == length(arguments)) {
      stop("--config requires exactly one path", call. = FALSE)
    }
    require_export_package("yaml")
    configured <- yaml::read_yaml(arguments[config_positions + 1L])
    configured <- configured[intersect(names(configured), names(defaults))]
    defaults[names(configured)] <- configured
  }

  index <- 1L
  while (index <= length(arguments)) {
    argument <- arguments[[index]]
    if (argument %in% c("--help", "-h")) {
      cat(
        "Usage:\n",
        "  Rscript scripts/clustering/export_explorer_assets.R \\\n",
        "    --input-dir outputs/clustering/RUN_ID \\\n",
        "    --output-dir outputs/explorer/clustering\n",
        sep = ""
      )
      quit(save = "no", status = 0L)
    }
    if (argument == "--config") {
      index <- index + 2L
      next
    }
    if (!startsWith(argument, "--") || index == length(arguments)) {
      stop("Arguments must be supplied as --name value pairs", call. = FALSE)
    }

    key <- gsub("-", "_", substring(argument, 3L), fixed = TRUE)
    if (!key %in% names(defaults)) {
      stop("Unknown argument: ", argument, call. = FALSE)
    }
    defaults[[key]] <- arguments[[index + 1L]]
    index <- index + 2L
  }

  if (is.na(defaults$input_dir) || !nzchar(defaults$input_dir)) {
    stop("--input-dir is required", call. = FALSE)
  }
  defaults
}

read_required_csv <- function(path) {
  if (!file.exists(path)) {
    stop("Required file does not exist: ", path, call. = FALSE)
  }
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

read_required_rds <- function(path) {
  if (!file.exists(path)) {
    stop("Required file does not exist: ", path, call. = FALSE)
  }
  readRDS(path)
}

clean_arrow_frame <- function(frame) {
  frame <- as.data.frame(frame, stringsAsFactors = FALSE, check.names = FALSE)
  for (column in names(frame)) {
    if (is.factor(frame[[column]])) {
      frame[[column]] <- as.character(frame[[column]])
    }
  }
  rownames(frame) <- NULL
  frame
}

write_parquet_frame <- function(frame, path) {
  arrow::write_parquet(
    clean_arrow_frame(frame),
    sink = path,
    compression = "zstd"
  )
}

align_size_factors <- function(size_factor, cell_ids) {
  size_factor <- as.numeric(size_factor)
  if (length(size_factor) != length(cell_ids) ||
      anyNA(size_factor) || any(!is.finite(size_factor)) ||
      any(size_factor <= 0)) {
    stop("Visualization bundle contains invalid size factors", call. = FALSE)
  }
  size_factor
}

prepare_cells <- function(bundle, cell_type) {
  embedding <- as.matrix(bundle$embedding)
  if (ncol(embedding) != 2L || is.null(rownames(embedding)) ||
      anyDuplicated(rownames(embedding))) {
    stop(cell_type, ": invalid clustering UMAP embedding", call. = FALSE)
  }

  cell_ids <- rownames(embedding)
  if (!identical(colnames(bundle$counts), cell_ids)) {
    stop(cell_type, ": count columns and clustering UMAP rows differ", call. = FALSE)
  }

  labels <- bundle$primary_labels
  if (!is.null(names(labels))) {
    if (!setequal(names(labels), cell_ids)) {
      stop(cell_type, ": named cluster labels do not match cells", call. = FALSE)
    }
    labels <- labels[cell_ids]
  }
  if (length(labels) != length(cell_ids) || anyNA(labels)) {
    stop(cell_type, ": invalid primary cluster labels", call. = FALSE)
  }

  metadata <- as.data.frame(bundle$cell_metadata, stringsAsFactors = FALSE)
  if (!is.null(rownames(metadata))) {
    if (!setequal(rownames(metadata), cell_ids)) {
      stop(cell_type, ": cell metadata does not match cells", call. = FALSE)
    }
    metadata <- metadata[cell_ids, , drop = FALSE]
  }

  wanted_metadata <- intersect(
    c(
      "Experiment", "Detection", "total_features_by_counts",
      "total_counts", "pct_counts_Mito"
    ),
    names(metadata)
  )

  cells <- data.frame(
    cell_id = cell_ids,
    UMAP_1 = as.numeric(embedding[, 1L]),
    UMAP_2 = as.numeric(embedding[, 2L]),
    cluster = as.character(labels),
    size_factor = align_size_factors(bundle$size_factor, cell_ids),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (length(wanted_metadata)) {
    cells <- cbind(
      cells,
      metadata[, wanted_metadata, drop = FALSE]
    )
  }

  if (anyDuplicated(cells$cell_id) || anyNA(cells$UMAP_1) ||
      anyNA(cells$UMAP_2)) {
    stop(cell_type, ": invalid exported cell table", call. = FALSE)
  }
  cells
}

prepare_markers <- function(marker_path, driver_path, cell_type) {
  markers <- read_required_csv(marker_path)
  required <- c(
    "feature", "gene_short_name", "type", "detection_difference",
    "magnitude_difference", "detection_q", "magnitude_q",
    "highest_detection_cluster", "highest_magnitude_cluster",
    "entered_clustering", "score"
  )
  missing <- setdiff(required, names(markers))
  if (length(missing)) {
    stop(
      cell_type, ": marker table lacks: ", paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  markers <- markers[
    !is.na(markers$type) & markers$type %in% c("on/off", "magnitude", "both"),
    ,
    drop = FALSE
  ]
  markers <- markers[order(markers$score, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
  marker_order <- markers$feature
  if (!nrow(markers) || anyDuplicated(markers$feature)) {
    stop(cell_type, ": supported marker features are empty or duplicated", call. = FALSE)
  }

  gene_name <- as.character(markers$gene_short_name)
  gene_name[is.na(gene_name) | !nzchar(gene_name)] <- markers$feature[
    is.na(gene_name) | !nzchar(gene_name)
  ]
  markers$display_name <- paste0(gene_name, " — ", markers$feature)
  markers$target_cluster <- ifelse(
    markers$type == "on/off",
    markers$highest_detection_cluster,
    markers$highest_magnitude_cluster
  )
  missing_target <- is.na(markers$target_cluster) | !nzchar(markers$target_cluster)
  markers$target_cluster[missing_target] <- markers$highest_detection_cluster[missing_target]

  drivers <- read_required_csv(driver_path)
  if (nrow(drivers)) {
    drivers$driver_rank <- seq_len(nrow(drivers))
    driver_columns <- intersect(
      c(
        "feature", "driver_rank", "between_cluster_variation",
        "driver_eta_squared", "centroid_contribution",
        "cumulative_contribution"
      ),
      names(drivers)
    )
    markers <- merge(
      markers,
      drivers[, driver_columns, drop = FALSE],
      by = "feature",
      all.x = TRUE,
      sort = FALSE
    )
    markers <- markers[match(marker_order, markers$feature), , drop = FALSE]
  }

  rownames(markers) <- NULL
  markers
}

prepare_expression <- function(bundle, cells, markers, cell_type) {
  features <- markers$feature
  missing <- setdiff(features, rownames(bundle$counts))
  if (length(missing)) {
    stop(
      cell_type, ": marker features absent from counts: ",
      paste(utils::head(missing, 5L), collapse = ", "),
      call. = FALSE
    )
  }

  raw <- as.matrix(bundle$counts[features, cells$cell_id, drop = FALSE])
  normalized <- log1p(sweep(raw, 2L, cells$size_factor, FUN = "/"))
  expression <- as.data.frame(
    t(normalized),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  expression <- cbind(
    data.frame(cell_id = cells$cell_id, stringsAsFactors = FALSE),
    expression
  )

  if (!identical(expression$cell_id, cells$cell_id) ||
      !identical(names(expression)[-1L], features) ||
      anyNA(expression[-1L]) || any(as.matrix(expression[-1L]) < 0)) {
    stop(cell_type, ": normalized marker-expression validation failed", call. = FALSE)
  }
  expression
}

prepare_summary <- function(
    run_row,
    result,
    diagnostics,
    markers,
    cells
) {
  primary_g <- as.integer(result$primary_g)
  primary <- diagnostics[diagnostics$g == primary_g, , drop = FALSE]
  if (nrow(primary) != 1L) {
    stop(result$cell_type, ": primary resolution is missing from diagnostics", call. = FALSE)
  }

  data.frame(
    cell_type = result$cell_type,
    retained_cells = nrow(cells),
    selected_pc_count = as.integer(result$screen$selected_pc_count),
    screen_g = as.integer(result$screen$screen_g),
    primary_g = primary_g,
    cluster_sizes = primary$cluster_sizes,
    bic_delta_from_one = primary$bic_delta_from_one,
    stability_median = primary$stability_median,
    status = as.character(run_row$status),
    marker_count = nrow(markers),
    elapsed_seconds = run_row$elapsed_seconds,
    stringsAsFactors = FALSE
  )
}

export_cell_type <- function(cell_type, run_row, input_dir, output_dir) {
  source_dir <- file.path(input_dir, cell_type)
  target_dir <- file.path(output_dir, cell_type)
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)

  bundle <- read_required_rds(file.path(source_dir, "visualization_bundle.rds"))
  result <- read_required_rds(file.path(source_dir, "result.rds"))
  diagnostics <- read_required_csv(file.path(source_dir, "resolution_diagnostics.csv"))
  markers <- prepare_markers(
    file.path(source_dir, "primary_markers.csv"),
    file.path(source_dir, "primary_drivers.csv"),
    cell_type
  )
  cells <- prepare_cells(bundle, cell_type)

  if (isTRUE(bundle$skipped) || as.integer(bundle$primary_g) <= 1L ||
      as.integer(result$primary_g) <= 1L) {
    stop(cell_type, ": attempted to export a non-accepted clustering", call. = FALSE)
  }
  if (as.integer(bundle$primary_g) != as.integer(result$primary_g)) {
    stop(cell_type, ": bundle and result disagree about primary G", call. = FALSE)
  }

  expression <- prepare_expression(bundle, cells, markers, cell_type)
  diagnostics$is_primary <- diagnostics$g == as.integer(result$primary_g)
  diagnostics$mitochondrial_policy <- "diagnostic_only"
  summary <- prepare_summary(run_row, result, diagnostics, markers, cells)

  write_parquet_frame(cells, file.path(target_dir, "cells.parquet"))
  write_parquet_frame(expression, file.path(target_dir, "expression.parquet"))
  write_parquet_frame(summary, file.path(target_dir, "summary.parquet"))
  write_parquet_frame(diagnostics, file.path(target_dir, "qc_by_resolution.parquet"))
  write_parquet_frame(markers, file.path(target_dir, "markers.parquet"))

  data.frame(
    schema_version = 1L,
    cell_type = cell_type,
    retained_cells = summary$retained_cells,
    screen_g = summary$screen_g,
    primary_g = summary$primary_g,
    marker_count = summary$marker_count,
    status = summary$status,
    asset_directory = cell_type,
    cells_file = file.path(cell_type, "cells.parquet"),
    expression_file = file.path(cell_type, "expression.parquet"),
    summary_file = file.path(cell_type, "summary.parquet"),
    qc_file = file.path(cell_type, "qc_by_resolution.parquet"),
    markers_file = file.path(cell_type, "markers.parquet"),
    stringsAsFactors = FALSE
  )
}

export_streamlit_clustering <- function(input_dir, output_dir) {
  require_export_package("arrow")
  require_export_package("Matrix")

  input_dir <- normalizePath(input_dir, mustWork = TRUE)
  output_dir <- normalizePath(
    output_dir,
    mustWork = FALSE
  )
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  run_summary <- read_required_csv(file.path(input_dir, "run_summary.csv"))
  accepted <- run_summary[
    run_summary$status == "accepted_multicluster" &
      !is.na(run_summary$primary_g) & run_summary$primary_g > 1L,
    ,
    drop = FALSE
  ]
  accepted <- accepted[order(accepted$cell_type), , drop = FALSE]
  if (!nrow(accepted)) {
    stop("No accepted multi-cluster cell types were found", call. = FALSE)
  }

  manifest_rows <- vector("list", nrow(accepted))
  for (index in seq_len(nrow(accepted))) {
    cell_type <- accepted$cell_type[[index]]
    message("Exporting ", cell_type, " (", index, "/", nrow(accepted), ")")
    manifest_rows[[index]] <- export_cell_type(
      cell_type,
      accepted[index, , drop = FALSE],
      input_dir,
      output_dir
    )
  }

  manifest <- do.call(rbind, manifest_rows)
  write_parquet_frame(manifest, file.path(output_dir, "manifest.parquet"))
  message(
    "Exported ", nrow(manifest), " accepted cell types and ",
    sum(manifest$marker_count), " marker/cell-type records to ", output_dir
  )
  invisible(manifest)
}

if (sys.nframe() == 0L) {
  arguments <- parse_export_arguments()
  export_streamlit_clustering(arguments$input_dir, arguments$output_dir)
}
