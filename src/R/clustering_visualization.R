# Interactive visualization for an already-fitted clustering result.
#
# This file is deliberately downstream-only: it never changes PC selection,
# mixture fitting, cluster labels, or any clustering diagnostic.

.require_visualization_package <- function(package, purpose) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Package '", package, "' is required for ", purpose)
  }
}

.with_visualization_seed <- function(seed, expression) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(as.integer(seed))
  force(expression)
}

.validate_visualization_size_factors <- function(size_factor, cell_names) {
  if (length(size_factor) != length(cell_names) ||
      anyNA(size_factor) ||
      any(!is.finite(size_factor)) ||
      any(size_factor <= 0)) {
    stop("size_factor must contain one finite positive value per cell")
  }

  if (!is.null(names(size_factor))) {
    if (anyDuplicated(names(size_factor)) ||
        !setequal(names(size_factor), cell_names)) {
      stop("Named size_factor values must match the visualization cells")
    }
    size_factor <- size_factor[cell_names]
  }

  as.numeric(size_factor)
}

.prepare_gene_labels <- function(features, gene_annotations) {
  labels <- data.frame(
    feature = features,
    gene_short_name = features,
    stringsAsFactors = FALSE
  )

  if (is.null(gene_annotations)) {
    return(labels)
  }

  if (is.data.frame(gene_annotations)) {
    if (!"feature" %in% colnames(gene_annotations)) {
      stop("gene_annotations must contain a 'feature' column")
    }
    if (anyDuplicated(gene_annotations$feature)) {
      stop("gene_annotations$feature must be unique")
    }
    if ("gene_short_name" %in% colnames(gene_annotations)) {
      annotation_values <- as.character(gene_annotations$gene_short_name)
      names(annotation_values) <- as.character(gene_annotations$feature)
    } else {
      annotation_values <- as.character(gene_annotations$feature)
      names(annotation_values) <- annotation_values
    }
  } else if (is.atomic(gene_annotations) && !is.null(names(gene_annotations))) {
    annotation_values <- as.character(gene_annotations)
  } else {
    stop(
      "gene_annotations must be NULL, a named vector, or a data frame ",
      "with feature and optional gene_short_name columns"
    )
  }

  matched <- annotation_values[features]
  usable <- !is.na(matched) & nzchar(matched)
  labels$gene_short_name[usable] <- matched[usable]
  labels
}

.add_gene_magnitude_ranking <- function(gene_labels, gene_ranking) {
  gene_labels$magnitude_difference <- NA_real_

  if (is.null(gene_ranking)) {
    return(gene_labels)
  }

  gene_ranking <- as.data.frame(gene_ranking, stringsAsFactors = FALSE)
  required <- c("feature", "magnitude_difference")
  missing <- setdiff(required, colnames(gene_ranking))
  if (length(missing)) {
    stop(
      "gene_ranking is missing column(s): ",
      paste(missing, collapse = ", ")
    )
  }
  if (anyDuplicated(gene_ranking$feature)) {
    stop("gene_ranking$feature must be unique")
  }

  magnitude <- suppressWarnings(as.numeric(gene_ranking$magnitude_difference))
  magnitude[!is.finite(magnitude)] <- NA_real_
  names(magnitude) <- as.character(gene_ranking$feature)
  gene_labels$magnitude_difference <- unname(
    magnitude[gene_labels$feature]
  )
  gene_labels
}

.rank_gene_labels <- function(gene_labels) {
  magnitude <- gene_labels$magnitude_difference
  order(
    is.na(magnitude),
    -ifelse(is.na(magnitude), -Inf, magnitude),
    gene_labels$gene_short_name,
    gene_labels$feature
  )
}

.prepare_visualization_metadata <- function(cell_names, cell_metadata) {
  if (is.null(cell_metadata)) {
    return(data.frame(row.names = cell_names))
  }

  cell_metadata <- as.data.frame(cell_metadata, stringsAsFactors = FALSE)
  if (nrow(cell_metadata) != length(cell_names)) {
    stop("cell_metadata must have one row per visualization cell")
  }

  if (!is.null(rownames(cell_metadata))) {
    if (anyDuplicated(rownames(cell_metadata)) ||
        !setequal(rownames(cell_metadata), cell_names)) {
      stop("cell_metadata row names must match the visualization cells")
    }
    cell_metadata <- cell_metadata[cell_names, , drop = FALSE]
  }

  rownames(cell_metadata) <- cell_names
  cell_metadata
}

#' Compute a UMAP strictly for visualizing an accepted clustering.
#'
#' @param counts Sparse raw count matrix, features by cells.
#' @param size_factor Positive normalization factors, one per cell.
#' @param pc_scores The exact selected PC score matrix supplied to the GMM,
#'   with cells in rows (normally pc_selection$scores).
#' @param resolution_ladder Result from build_resolution_ladder().
#' @param cell_metadata Optional cell-level data frame, indexed by cell name.
#' @param gene_annotations Optional named gene-symbol vector or data frame.
#' @param gene_ranking Optional marker table containing feature and
#'   magnitude_difference columns.
#' @return A visualization bundle, or a skipped bundle when primary_g is one.
prepare_cluster_visualization <- function(
    counts,
    size_factor,
    pc_scores,
    resolution_ladder,
    cell_metadata = NULL,
    gene_annotations = NULL,
    gene_ranking = NULL,
    seed = 20260809L,
    n_neighbors = NULL,
    min_dist = 0.30,
    metric = "euclidean"
) {
  if (is.null(rownames(counts)) || is.null(colnames(counts))) {
    stop("counts must have unique feature row names and cell column names")
  }
  if (anyDuplicated(rownames(counts)) || anyDuplicated(colnames(counts))) {
    stop("counts must have unique feature and cell names")
  }
  if (!is.list(resolution_ladder) || is.null(resolution_ladder$primary_g) ||
      is.null(resolution_ladder$primary_labels)) {
    stop("resolution_ladder must contain primary_g and primary_labels")
  }

  primary_g <- as.integer(resolution_ladder$primary_g)
  if (length(primary_g) != 1L || is.na(primary_g) || primary_g < 1L) {
    stop("resolution_ladder$primary_g must be one positive integer")
  }

  if (primary_g <= 1L) {
    return(structure(
      list(
        skipped = TRUE,
        reason = "No accepted multi-cluster solution: primary_g is 1.",
        primary_g = primary_g
      ),
      class = "cluster_visualization_skipped"
    ))
  }

  .require_visualization_package("uwot", "UMAP visualization")

  pc_scores <- as.matrix(pc_scores)
  if (nrow(pc_scores) < 3L || ncol(pc_scores) < 1L ||
      is.null(rownames(pc_scores)) || anyDuplicated(rownames(pc_scores))) {
    stop("pc_scores must have at least three uniquely named cells and one PC")
  }
  if (anyNA(pc_scores) || any(!is.finite(pc_scores))) {
    stop("pc_scores must contain only finite values")
  }

  cell_names <- rownames(pc_scores)
  if (!all(cell_names %in% colnames(counts))) {
    stop("Every pc_scores row must occur among counts columns")
  }
  if (ncol(counts) != length(cell_names) ||
      !setequal(colnames(counts), cell_names)) {
    stop("counts and pc_scores must contain exactly the same cells")
  }
  counts <- counts[, cell_names, drop = FALSE]
  size_factor <- .validate_visualization_size_factors(size_factor, cell_names)
  cell_metadata <- .prepare_visualization_metadata(cell_names, cell_metadata)

  primary_labels <- factor(resolution_ladder$primary_labels)
  if (length(primary_labels) != length(cell_names) || anyNA(primary_labels)) {
    stop("primary_labels must contain one non-missing label per visualization cell")
  }
  if (!is.null(names(resolution_ladder$primary_labels))) {
    label_names <- names(resolution_ladder$primary_labels)
    if (anyDuplicated(label_names) || !setequal(label_names, cell_names)) {
      stop("Named primary_labels must match pc_scores row names")
    }
    primary_labels <- factor(resolution_ladder$primary_labels[cell_names])
  }
  if (nlevels(primary_labels) != primary_g) {
    stop("primary_labels do not contain primary_g clusters")
  }

  if (is.null(n_neighbors)) {
    n_neighbors <- min(15L, nrow(pc_scores) - 1L)
  }
  if (length(n_neighbors) != 1L || !is.finite(n_neighbors) ||
      n_neighbors < 2L || n_neighbors >= nrow(pc_scores) ||
      n_neighbors != as.integer(n_neighbors)) {
    stop("n_neighbors must be an integer from 2 through nrow(pc_scores) - 1")
  }
  if (length(min_dist) != 1L || !is.finite(min_dist) || min_dist < 0) {
    stop("min_dist must be one finite nonnegative number")
  }
  if (length(metric) != 1L || !is.character(metric) || !nzchar(metric)) {
    stop("metric must be one non-empty character value")
  }

  embedding <- .with_visualization_seed(seed, uwot::umap(
    X = pc_scores,
    n_neighbors = as.integer(n_neighbors),
    min_dist = min_dist,
    metric = metric,
    n_components = 2L,
    n_threads = 1L,
    verbose = FALSE,
    ret_model = FALSE
  ))
  rownames(embedding) <- cell_names
  colnames(embedding) <- c("UMAP_1", "UMAP_2")

  structure(list(
    skipped = FALSE,
    primary_g = primary_g,
    primary_labels = primary_labels,
    embedding = embedding,
    counts = counts,
    size_factor = size_factor,
    cell_metadata = cell_metadata,
    gene_labels = .add_gene_magnitude_ranking(
      .prepare_gene_labels(rownames(counts), gene_annotations),
      gene_ranking
    ),
    settings = list(
      seed = as.integer(seed),
      n_neighbors = as.integer(n_neighbors),
      min_dist = min_dist,
      metric = metric,
      input_space = "selected PC scores used for GMM"
    )
  ), class = "cluster_visualization_bundle")
}

print.cluster_visualization_skipped <- function(x, ...) {
  cat(x$reason, "\n")
  invisible(x)
}

.viewer_expression_data <- function(bundle, feature) {
  # Sparse Matrix subsetting is an S4 method. Ensure that the namespace is
  # loaded before the first `[` call in a fresh R/Shiny session.
  .require_visualization_package("Matrix", "sparse expression access")

  raw_count <- as.numeric(bundle$counts[feature, , drop = TRUE])
  normalized_expression <- log1p(raw_count / bundle$size_factor)
  data.frame(
    cell = colnames(bundle$counts),
    UMAP_1 = bundle$embedding[, "UMAP_1"],
    UMAP_2 = bundle$embedding[, "UMAP_2"],
    cluster = bundle$primary_labels,
    raw_count = raw_count,
    normalized_expression = normalized_expression,
    detected = raw_count > 0,
    stringsAsFactors = FALSE
  )
}

.viewer_hover_text <- function(data, gene_label) {
  paste0(
    "Cell: ", htmltools::htmlEscape(data$cell),
    "<br>Cluster: ", htmltools::htmlEscape(as.character(data$cluster)),
    "<br>", htmltools::htmlEscape(gene_label),
    " raw UMI: ", data$raw_count,
    "<br>Normalized expression: ",
    formatC(data$normalized_expression, digits = 3L, format = "f")
  )
}

#' Build a local or deployable Shiny viewer for a visualization bundle.
cluster_viewer_app <- function(bundle) {
  .require_visualization_package("Matrix", "sparse expression access")
  .require_visualization_package("shiny", "the cluster viewer")
  .require_visualization_package("plotly", "the cluster viewer")
  .require_visualization_package("viridisLite", "expression colouring")

  if (inherits(bundle, "cluster_visualization_skipped") ||
      isTRUE(bundle$skipped)) {
    stop("The viewer is available only when primary_g is greater than one")
  }
  if (!inherits(bundle, "cluster_visualization_bundle")) {
    stop("bundle must be created by prepare_cluster_visualization()")
  }

  if (!"magnitude_difference" %in% colnames(bundle$gene_labels)) {
    bundle$gene_labels$magnitude_difference <- NA_real_
  }

  gene_order <- .rank_gene_labels(bundle$gene_labels)
  ordered_genes <- bundle$gene_labels[gene_order, , drop = FALSE]
  magnitude_label <- ifelse(
    is.na(ordered_genes$magnitude_difference),
    "",
    paste0(
      " | magnitude ",
      formatC(
        ordered_genes$magnitude_difference,
        digits = 3L,
        format = "f"
      )
    )
  )

  gene_choices <- setNames(
    ordered_genes$feature,
    paste0(
      ordered_genes$gene_short_name,
      " — ",
      ordered_genes$feature,
      magnitude_label
    )
  )

  ui <- shiny::fluidPage(
    shiny::titlePanel(
      paste0("Accepted clustering: G = ", bundle$primary_g)
    ),
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        shiny::radioButtons(
          "colour_mode",
          "Colour cells by",
          choices = c("Cluster" = "cluster", "Gene expression" = "expression"),
          selected = "cluster"
        ),
        shiny::selectizeInput(
          "gene",
          "Gene for expression view (ranked by cluster magnitude)",
          choices = NULL,
          options = list(
            placeholder = "Type a gene name or feature ID",
            maxOptions = 500L
          )
        ),
        shiny::checkboxInput(
          "clip_expression",
          "Clip colour scale at the 99th percentile of positive cells",
          value = TRUE
        ),
        shiny::helpText(
          paste(
            "Genes are ranked by the largest between-cluster difference in",
            "mean normalized expression among positive cells.",
            "Grey points are zero UMI. Expression is log1p(raw UMI / size factor)."
          )
        )
      ),
      shiny::mainPanel(
        plotly::plotlyOutput("umap_plot", height = "720px"),
        shiny::verbatimTextOutput("viewer_note")
      )
    )
  )

  server <- function(input, output, session) {
      shiny::updateSelectizeInput(
        session,
        "gene",
        choices = gene_choices,
        selected = ordered_genes$feature[1L],
        server = TRUE
      )

      output$viewer_note <- shiny::renderText({
        paste0(
          "UMAP is a readability-only display of the selected PC space used ",
          "for the GMM. It does not alter clusters or their diagnostics."
        )
      })

      output$umap_plot <- plotly::renderPlotly({
        shiny::req(input$gene)

        gene_row <- bundle$gene_labels[
          match(input$gene, bundle$gene_labels$feature),
          ,
          drop = FALSE
        ]
        if (!nrow(gene_row)) {
          shiny::validate(shiny::need(FALSE, "Choose a valid gene."))
        }

        data <- .viewer_expression_data(bundle, input$gene)
        hover_text <- .viewer_hover_text(
          data,
          paste0(gene_row$gene_short_name, " (", gene_row$feature, ")")
        )

        if (identical(input$colour_mode, "cluster")) {
          return(plotly::plot_ly(
            data = data,
            x = ~UMAP_1,
            y = ~UMAP_2,
            type = "scattergl",
            mode = "markers",
            color = ~cluster,
            text = hover_text,
            hoverinfo = "text",
            marker = list(size = 7, opacity = 0.82)
          ) |> plotly::layout(
            xaxis = list(title = "UMAP 1"),
            yaxis = list(title = "UMAP 2"),
            legend = list(title = list(text = "Primary cluster"))
          ))
        }

        zero_data <- data[!data$detected, , drop = FALSE]
        positive_data <- data[data$detected, , drop = FALSE]
        zero_hover <- hover_text[!data$detected]
        positive_hover <- hover_text[data$detected]

        plot <- plotly::plot_ly(
          data = zero_data,
          x = ~UMAP_1,
          y = ~UMAP_2,
          type = "scattergl",
          mode = "markers",
          name = "Zero UMI",
          text = zero_hover,
          hoverinfo = "text",
          marker = list(size = 7, color = "#BDBDBD", opacity = 0.65)
        )

        if (nrow(positive_data)) {
          colour_max <- max(positive_data$normalized_expression)
          if (isTRUE(input$clip_expression) && nrow(positive_data) >= 2L) {
            colour_max <- as.numeric(stats::quantile(
              positive_data$normalized_expression,
              probs = 0.99,
              names = FALSE
            ))
          }
          colour_max <- max(colour_max, .Machine$double.eps)

          plot <- plotly::add_trace(
            plot,
            data = positive_data,
            x = ~UMAP_1,
            y = ~UMAP_2,
            type = "scattergl",
            mode = "markers",
            name = "Positive UMI",
            text = positive_hover,
            hoverinfo = "text",
            marker = list(
              size = 7,
              opacity = 0.86,
              color = positive_data$normalized_expression,
              cmin = 0,
              cmax = colour_max,
              colorscale = viridisLite::viridis(256L),
              showscale = TRUE,
              colorbar = list(title = "log1p(UMI / size factor)")
            )
          )
        }

        plotly::layout(
          plot,
          xaxis = list(title = "UMAP 1"),
          yaxis = list(title = "UMAP 2")
        )
      })
  }

  shiny::shinyApp(ui = ui, server = server)
}

#' Launch the local viewer from a bundle object or an RDS bundle path.
run_cluster_viewer <- function(bundle_or_path, launch.browser = interactive()) {
  bundle_path <- if (
    is.character(bundle_or_path) && length(bundle_or_path) == 1L
  ) {
    bundle_or_path
  } else {
    NULL
  }

  bundle <- if (!is.null(bundle_path)) {
    readRDS(bundle_path)
  } else {
    bundle_or_path
  }

  # Existing bundles predate embedded marker rankings. When launched from a
  # path, automatically reuse the marker table saved beside the bundle.
  if (
    !is.null(bundle_path) &&
      (
        !"magnitude_difference" %in% colnames(bundle$gene_labels) ||
          all(is.na(bundle$gene_labels$magnitude_difference))
      )
  ) {
    marker_path <- file.path(dirname(bundle_path), "primary_markers.csv")
    if (file.exists(marker_path)) {
      marker_table <- utils::read.csv(marker_path, stringsAsFactors = FALSE)
      bundle$gene_labels <- .add_gene_magnitude_ranking(
        bundle$gene_labels,
        marker_table
      )
    }
  }

  shiny::runApp(cluster_viewer_app(bundle), launch.browser = launch.browser)
}

save_cluster_visualization_bundle <- function(bundle, path) {
  if (!inherits(bundle, "cluster_visualization_bundle") &&
      !inherits(bundle, "cluster_visualization_skipped")) {
    stop("bundle must be created by prepare_cluster_visualization()")
  }
  saveRDS(bundle, path)
  invisible(path)
}
