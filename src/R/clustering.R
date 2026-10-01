

select_cell_type <- function(
    sce,
    cell_type,
    cell_type_col = "Cell.type"
) {
  if (!cell_type_col %in% colnames(colData(sce))) {
    stop("Cell-type column not found: ", cell_type_col)
  }
  
  cell_type_labels <- colData(sce)[[cell_type_col]]
  
  keep <- !is.na(cell_type_labels) &
    as.character(cell_type_labels) == cell_type
  
  if (!any(keep)) {
    stop("No cells found for cell type: ", cell_type)
  }
  
  sce[, keep, drop = FALSE]
}

filter_rare_experiments <- function(
    sce,
    experiment_col = "Experiment",
    min_cells = 5
) {
  if (!experiment_col %in% colnames(colData(sce))) {
    stop("Experiment column not found: ", experiment_col)
  }
  
  if (
    length(min_cells) != 1L ||
    !is.finite(min_cells) ||
    min_cells < 1 ||
    min_cells != as.integer(min_cells)
  ) {
    stop("min_cells must be one positive integer")
  }
  
  experiment <- as.character(colData(sce)[[experiment_col]])
  
  if (anyNA(experiment) || any(trimws(experiment) == "")) {
    stop("Every selected cell must have a known experiment")
  }
  
  experiment_counts <- table(experiment)
  
  retained_experiments <- names(experiment_counts)[
    experiment_counts >= min_cells
  ]
  
  keep <- experiment %in% retained_experiments
  
  if (!any(keep)) {
    stop("No cells remain after rare-experiment filtering")
  }
  
  list(
    sce = sce[, keep, drop = FALSE],
    experiment_counts = experiment_counts,
    retained_experiments = retained_experiments,
    excluded_experiments = setdiff(
      names(experiment_counts),
      retained_experiments
    ),
    excluded_cells = colnames(sce)[!keep]
  )
}





normalize_counts <- function(counts, size_factor) {
  if (
    length(size_factor) != ncol(counts) ||
    anyNA(size_factor) ||
    any(!is.finite(size_factor)) ||
    any(size_factor <= 0)
  ) {
    stop("size_factor must contain one finite positive value per cell")
  }

  normalized <- counts %*% Matrix::Diagonal(x = 1 / size_factor)
  normalized <- log1p(normalized)
  dimnames(normalized) <- dimnames(counts)

  normalized
}



select_stress_features <- function(
    counts,
    stress_genes,
    min_detected = max(3L, ceiling(0.01 * ncol(counts))),
    min_total_count = 10L
) {
  if (is.null(rownames(counts))) {
    stop("counts must have gene row names")
  }
  
  stress_genes <- unique(as.character(stress_genes))
  stress_genes <- stress_genes[
    !is.na(stress_genes) & nzchar(stress_genes)
  ]
  
  if (!length(stress_genes)) {
    stop("stress_genes is empty")
  }
  
  matched <- intersect(stress_genes, rownames(counts))
  missing <- setdiff(stress_genes, rownames(counts))
  
  if (!length(matched)) {
    stop("None of the supplied stress genes occur in counts")
  }
  
  stress_counts <- counts[matched, , drop = FALSE]
  
  detected_cells <- Matrix::rowSums(stress_counts > 0)
  total_count <- Matrix::rowSums(stress_counts)
  
  usable <- detected_cells >= min_detected &
    total_count >= min_total_count
  
  if (sum(usable) < 2L) {
    stop("Fewer than two stress genes pass the coverage filters")
  }
  
  audit <- data.frame(
    feature = matched,
    detected_cells = as.integer(detected_cells),
    detection_fraction = as.numeric(detected_cells) / ncol(counts),
    total_count = as.numeric(total_count),
    usable = usable,
    row.names = NULL
  )
  
  list(
    supplied = stress_genes,
    matched = matched,
    missing = missing,
    usable = matched[usable],
    audit = audit
  )
}

make_group_design <- function(known_group) {
  if (!length(known_group)) {
    stop("known_group is empty")
  }
  
  labels <- as.character(known_group)
  
  if (anyNA(labels) || any(trimws(labels) == "")) {
    stop("Every cell must have a known-group label")
  }
  
  group <- droplevels(factor(labels))
  
  if (nlevels(group) == 1L) {
    design <- matrix(
      1,
      nrow = length(group),
      ncol = 1L
    )
    
    colnames(design) <- "(Intercept)"
  } else {
    design <- model.matrix(~ group)
  }
  
  if (!is.null(names(known_group))) {
    rownames(design) <- names(known_group)
  }
  
  if (qr(design)$rank != ncol(design)) {
    stop("Known-group design matrix is not full rank")
  }
  
  attr(design, "group_levels") <- levels(group)
  
  design
}

residualize_features <- function(expression, design) {
  if (length(dim(expression)) != 2L) {
    stop("expression must be a two-dimensional matrix")
  }
  
  if (length(dim(design)) != 2L) {
    stop("design must be a two-dimensional matrix")
  }
  
  if (nrow(expression) == 0L || ncol(expression) == 0L) {
    stop("expression cannot be empty")
  }
  
  if (ncol(expression) != nrow(design)) {
    stop(
      "The number of expression columns must equal ",
      "the number of design rows"
    )
  }
  
  if (
    !is.null(colnames(expression)) &&
    !is.null(rownames(design)) &&
    !identical(colnames(expression), rownames(design))
  ) {
    stop("Expression columns and design rows are not in the same order")
  }
  
  expression_dense <- as.matrix(expression)
  design_dense <- as.matrix(design)
  
  storage.mode(expression_dense) <- "double"
  storage.mode(design_dense) <- "double"
  
  if (any(!is.finite(expression_dense))) {
    stop("expression contains non-finite values")
  }
  
  if (any(!is.finite(design_dense))) {
    stop("design contains non-finite values")
  }
  
  design_qr <- qr(design_dense)
  
  if (design_qr$rank != ncol(design_dense)) {
    stop("design matrix is not full rank")
  }
  
  residuals <- t(
    qr.resid(
      design_qr,
      t(expression_dense)
    )
  )
  
  dimnames(residuals) <- dimnames(expression)
  
  residuals
}


estimate_stress_score <- function(
    normalized,
    stress_features,
    experiment_design
) {
  if (is.null(rownames(normalized))) {
    stop("normalized must have feature row names")
  }
  
  if (ncol(normalized) < 2L) {
    stop("At least two cells are required")
  }
  
  stress_features <- unique(as.character(stress_features))
  missing_features <- setdiff(stress_features, rownames(normalized))
  
  if (length(missing_features)) {
    stop(
      "Stress features missing from normalized matrix: ",
      paste(missing_features, collapse = ", ")
    )
  }
  
  stress_values <- as.matrix(
    normalized[stress_features, , drop = FALSE]
  )
  normalized_stress_features <- rownames(stress_values)
  
  row_variance <- function(x) {
    number_of_cells <- ncol(x)
    means <- rowMeans(x)
    
    pmax(
      0,
      (
        rowSums(x * x) -
          number_of_cells * means * means
      ) / (number_of_cells - 1L)
    )
  }
  
  variance_before <- row_variance(stress_values)
  variable_before <- variance_before > 0
  excluded_after_normalization <- normalized_stress_features[
    !variable_before
  ]
  
  if (sum(variable_before) < 2L) {
    stop("Fewer than two stress features vary after normalization")
  }
  
  stress_values <- stress_values[
    variable_before,
    ,
    drop = FALSE
  ]
  
  adjusted_stress <- residualize_features(
    expression = stress_values,
    design = experiment_design
  )
  adjusted_stress_features <- rownames(adjusted_stress)
  
  variance_after <- row_variance(adjusted_stress)
  variable_after <- variance_after > .Machine$double.eps
  excluded_after_experiment_adjustment <- adjusted_stress_features[
    !variable_after
  ]
  
  if (sum(variable_after) < 2L) {
    stop(
      paste(
        "Fewer than two stress features vary",
        "after experiment adjustment"
      )
    )
  }
  
  adjusted_stress <- adjusted_stress[
    variable_after,
    ,
    drop = FALSE
  ]
  
  if (!requireNamespace("irlba", quietly = TRUE)) {
    stop("Package 'irlba' is required for stress-score PCA")
  }
  
  pca <- irlba::prcomp_irlba(
    t(adjusted_stress),
    n = 1L,
    center = TRUE,
    scale. = FALSE
  )
  
  score <- as.numeric(pca$x[, 1L])
  loading <- as.numeric(pca$rotation[, 1L])
  
  # PCA signs are arbitrary. Orient the score so that a higher score
  # generally corresponds to higher average stress-gene expression.
  expression_anchor <- colMeans(stress_values)
  
  anchor_correlation <- suppressWarnings(
    stats::cor(score, expression_anchor)
  )
  
  if (
    is.finite(anchor_correlation) &&
    anchor_correlation < 0
  ) {
    score <- -score
    loading <- -loading
  }
  
  # Put the score on a mean-zero, standard-deviation-one scale.
  score <- as.numeric(scale(score))
  
  names(score) <- colnames(normalized)
  names(loading) <- rownames(adjusted_stress)
  
  list(
    score = score,
    loadings = loading,
    adjusted_expression = adjusted_stress,
    features_used = rownames(adjusted_stress),
    excluded_after_normalization =
      excluded_after_normalization,
    excluded_after_experiment_adjustment =
      excluded_after_experiment_adjustment,
    anchor_correlation = abs(anchor_correlation),
    pca_sdev = pca$sdev[1L]
  )
}


fit_stress_mixture <- function(
    stress_score,
    bic_improvement = 10,
    min_component_prop = 0.05
) {
  if (!requireNamespace("mclust", quietly = TRUE)) {
    stop("Package 'mclust' is required")
  }
  
  observation_names <- names(stress_score)
  stress_score <- as.numeric(stress_score)
  
  if (length(stress_score) < 2L) {
    stop("stress_score must contain at least two cells")
  }
  
  if (anyNA(stress_score) || any(!is.finite(stress_score))) {
    stop("stress_score must contain only finite values")
  }
  
  if (stats::sd(stress_score) <= .Machine$double.eps) {
    stop("stress_score has no variation")
  }
  
  if (
    length(bic_improvement) != 1L ||
    !is.finite(bic_improvement) ||
    bic_improvement < 0
  ) {
    stop("bic_improvement must be one nonnegative number")
  }
  
  if (
    length(min_component_prop) != 1L ||
    !is.finite(min_component_prop) ||
    min_component_prop <= 0 ||
    min_component_prop > 0.5
  ) {
    stop("min_component_prop must be between 0 and 0.5")
  }
  
  # mclust::Mclust looks for mclustBIC in its calling environment.
  mclustBIC <- mclust::mclustBIC
  
  fit_one <- tryCatch(
    mclust::Mclust(
      data = stress_score,
      G = 1L,
      modelNames = "E",
      verbose = FALSE
    ),
    error = function(error) NULL
  )
  
  fit_two <- tryCatch(
    mclust::Mclust(
      data = stress_score,
      G = 2L,
      modelNames = "V",
      verbose = FALSE
    ),
    error = function(error) NULL
  )
  
  bic_one <- if (is.null(fit_one)) {
    -Inf
  } else {
    as.numeric(fit_one$bic)
  }
  
  bic_two <- if (is.null(fit_two)) {
    -Inf
  } else {
    as.numeric(fit_two$bic)
  }
  
  if (is.null(fit_two)) {
    component_proportions <- numeric()
  } else {
    component_proportions <- as.numeric(
      table(
        factor(
          fit_two$classification,
          levels = 1:2
        )
      ) / length(stress_score)
    )
    
    names(component_proportions) <- c(
      "component_1",
      "component_2"
    )
  }
  
  accepted <- (
    !is.null(fit_one) &&
      !is.null(fit_two) &&
      is.finite(bic_one) &&
      is.finite(bic_two) &&
      bic_two - bic_one >= bic_improvement &&
      length(component_proportions) == 2L &&
      all(component_proportions >= min_component_prop)
  )
  
  posterior <- NULL
  high_stress_component <- NA_integer_
  
  if (accepted) {
    component_means <- as.numeric(
      fit_two$parameters$mean
    )
    
    high_stress_component <- which.max(component_means)
    
    posterior <- as.numeric(
      fit_two$z[, high_stress_component]
    )
    
    names(posterior) <- observation_names
  }
  
  list(
    accepted = accepted,
    posterior = posterior,
    high_stress_component = high_stress_component,
    bic = c(
      one_component = bic_one,
      two_component = bic_two
    ),
    bic_improvement = bic_two - bic_one,
    component_proportions = component_proportions,
    fit_one = fit_one,
    fit_two = fit_two
  )
}



make_nuisance_design <- function(
    known_group,
    stress_score,
    size_factor,
    stress_posterior = NULL,
    observation_names = NULL
) {
  number_of_cells <- length(known_group)
  
  if (number_of_cells == 0L) {
    stop("known_group is empty")
  }
  
  if (length(stress_score) != number_of_cells) {
    stop("stress_score must have one value per cell")
  }
  
  if (length(size_factor) != number_of_cells) {
    stop("size_factor must have one value per cell")
  }
  
  if (
    anyNA(stress_score) ||
    any(!is.finite(stress_score))
  ) {
    stop("stress_score must contain only finite values")
  }
  
  if (
    anyNA(size_factor) ||
    any(!is.finite(size_factor)) ||
    any(size_factor <= 0)
  ) {
    stop("size_factor must contain finite positive values")
  }
  
  if (
    !is.null(stress_posterior) &&
    length(stress_posterior) != number_of_cells
  ) {
    stop("stress_posterior must have one value per cell")
  }
  
  if (
    !is.null(stress_posterior) &&
    (
      anyNA(stress_posterior) ||
      any(!is.finite(stress_posterior)) ||
      any(stress_posterior < 0) ||
      any(stress_posterior > 1)
    )
  ) {
    stop("stress_posterior must contain values between zero and one")
  }
  
  if (
    !is.null(observation_names) &&
    length(observation_names) != number_of_cells
  ) {
    stop("observation_names must have one name per cell")
  }
  
  group_design <- make_group_design(known_group)
  
  design <- cbind(
    group_design,
    stress_score = as.numeric(stress_score),
    log_size_factor = log(as.numeric(size_factor))
  )
  
  if (!is.null(stress_posterior)) {
    design <- cbind(
      design,
      stress_posterior = as.numeric(stress_posterior)
    )
  }
  
  if (!is.null(observation_names)) {
    rownames(design) <- observation_names
  }
  
  design_qr <- qr(design)
  
  if (design_qr$rank != ncol(design)) {
    stop(
      paste(
        "Combined nuisance design is not full rank.",
        "One nuisance variable may be redundant with another."
      )
    )
  }
  
  design
}


select_candidate_features <- function(
    counts,
    excluded_features = character(),
    min_detected = max(3L, ceiling(0.01 * ncol(counts))),
    min_total_count = 10L
) {
  if (is.null(rownames(counts))) {
    stop("counts must have feature row names")
  }
  
  if (anyDuplicated(rownames(counts))) {
    stop("counts must have unique feature row names")
  }
  
  if (ncol(counts) == 0L) {
    stop("counts contains no cells")
  }
  
  excluded_features <- unique(as.character(excluded_features))
  excluded_features <- excluded_features[
    !is.na(excluded_features) &
      nzchar(excluded_features)
  ]
  
  excluded_matched <- intersect(
    excluded_features,
    rownames(counts)
  )
  
  excluded_missing <- setdiff(
    excluded_features,
    rownames(counts)
  )
  
  candidate_features <- setdiff(
    rownames(counts),
    excluded_matched
  )
  
  if (!length(candidate_features)) {
    stop("No features remain after exclusions")
  }
  
  candidate_counts <- counts[
    candidate_features,
    ,
    drop = FALSE
  ]
  
  detected_cells <- Matrix::rowSums(
    candidate_counts > 0
  )
  
  total_count <- Matrix::rowSums(
    candidate_counts
  )
  
  eligible <- (
    detected_cells >= min_detected &
      total_count >= min_total_count
  )
  
  if (sum(eligible) < 2L) {
    stop(
      paste(
        "Fewer than two non-excluded features",
        "pass the count filters"
      )
    )
  }
  
  audit <- data.frame(
    feature = candidate_features,
    detected_cells = as.integer(detected_cells),
    detection_fraction =
      as.numeric(detected_cells) / ncol(counts),
    total_count = as.numeric(total_count),
    eligible = eligible,
    row.names = NULL
  )
  
  list(
    eligible = candidate_features[eligible],
    excluded_matched = excluded_matched,
    excluded_missing = excluded_missing,
    audit = audit,
    thresholds = list(
      min_detected = min_detected,
      min_total_count = min_total_count
    )
  )
}


select_high_variance_features <- function(
    residuals,
    max_features = 3000L
) {
  if (length(dim(residuals)) != 2L) {
    stop("residuals must be a two-dimensional matrix")
  }
  
  if (is.null(rownames(residuals))) {
    stop("residuals must have feature row names")
  }
  
  if (anyDuplicated(rownames(residuals))) {
    stop("residuals must have unique feature row names")
  }
  
  if (ncol(residuals) < 2L) {
    stop("At least two cells are required")
  }
  
  if (
    length(max_features) != 1L ||
    !is.finite(max_features) ||
    max_features < 2L ||
    max_features != as.integer(max_features)
  ) {
    stop("max_features must be one integer of at least two")
  }
  
  residuals <- as.matrix(residuals)
  
  if (any(!is.finite(residuals))) {
    stop("residuals contains non-finite values")
  }
  
  number_of_cells <- ncol(residuals)
  feature_means <- rowMeans(residuals)
  
  residual_variance <- pmax(
    0,
    (
      rowSums(residuals * residuals) -
        number_of_cells * feature_means^2
    ) / (number_of_cells - 1L)
  )
  
  variable <- residual_variance > .Machine$double.eps
  
  if (sum(variable) < 2L) {
    stop(
      paste(
        "Fewer than two features vary",
        "after nuisance adjustment"
      )
    )
  }
  
  variable_index <- which(variable)
  
  # Feature names break exact variance ties deterministically.
  ordered_index <- variable_index[
    order(
      -residual_variance[variable_index],
      rownames(residuals)[variable_index]
    )
  ]
  
  number_selected <- min(
    length(ordered_index),
    as.integer(max_features)
  )
  
  selected_index <- ordered_index[
    seq_len(number_selected)
  ]
  
  selected <- rep(FALSE, nrow(residuals))
  selected[selected_index] <- TRUE
  
  audit <- data.frame(
    feature = rownames(residuals),
    residual_variance = residual_variance,
    variable = variable,
    selected = selected,
    variance_rank = NA_integer_,
    row.names = NULL
  )
  
  audit$variance_rank[ordered_index] <- seq_along(
    ordered_index
  )
  
  list(
    features = rownames(residuals)[selected_index],
    residuals = residuals[selected_index, , drop = FALSE],
    audit = audit,
    number_variable = sum(variable),
    number_selected = number_selected
  )
}

compute_residual_pca <- function(
    residuals,
    max_pcs = 20L
) {
  if (!requireNamespace("irlba", quietly = TRUE)) {
    stop("Package 'irlba' is required for PCA")
  }
  
  if (length(dim(residuals)) != 2L) {
    stop("residuals must be a two-dimensional matrix")
  }
  
  if (nrow(residuals) < 2L) {
    stop("At least two features are required for PCA")
  }
  
  if (ncol(residuals) < 2L) {
    stop("At least two cells are required for PCA")
  }
  
  if (
    length(max_pcs) != 1L ||
    !is.finite(max_pcs) ||
    max_pcs < 1L ||
    max_pcs != as.integer(max_pcs)
  ) {
    stop("max_pcs must be one positive integer")
  }
  
  residuals <- as.matrix(residuals)
  
  if (any(!is.finite(residuals))) {
    stop("residuals contains non-finite values")
  }
  
  # PCA expects cells in rows and features in columns.
  cell_by_feature <- t(residuals)
  
  pca_rank <- min(
    as.integer(max_pcs),
    nrow(cell_by_feature) - 1L,
    ncol(cell_by_feature) - 1L
  )
  
  if (pca_rank < 1L) {
    stop("The residual matrix does not support PCA")
  }
  
  pca <- irlba::prcomp_irlba(
    cell_by_feature,
    n = pca_rank,
    center = TRUE,
    scale. = FALSE
  )
  
  eigenvalues <- pca$sdev^2
  
  # Total variance across all input features.
  number_of_cells <- nrow(cell_by_feature)
  feature_means <- colMeans(cell_by_feature)
  
  total_variance <- (
    sum(cell_by_feature * cell_by_feature) -
      number_of_cells * sum(feature_means^2)
  ) / (number_of_cells - 1L)
  
  total_variance <- max(0, total_variance)
  
  variance_fraction_total <- if (total_variance > 0) {
    eigenvalues / total_variance
  } else {
    rep(0, length(eigenvalues))
  }
  
  variance_fraction_returned <- (
    eigenvalues / sum(eigenvalues)
  )
  
  rownames(pca$x) <- colnames(residuals)
  rownames(pca$rotation) <- rownames(residuals)
  
  list(
    scores = pca$x,
    loadings = pca$rotation,
    standard_deviations = pca$sdev,
    eigenvalues = eigenvalues,
    variance_fraction_total = variance_fraction_total,
    variance_fraction_returned = variance_fraction_returned,
    center = pca$center,
    rank = pca_rank
  )
}

select_clustering_pcs <- function(
    pca_result,
    elbow_search = 10L,
    minimum_pcs = 2L
) {
  eigenvalues <- as.numeric(pca_result$eigenvalues)
  scores <- as.matrix(pca_result$scores)
  
  if (!length(eigenvalues)) {
    stop("pca_result contains no eigenvalues")
  }
  
  if (ncol(scores) != length(eigenvalues)) {
    stop("PCA scores and eigenvalues do not match")
  }
  
  if (anyNA(eigenvalues) || any(!is.finite(eigenvalues))) {
    stop("PCA eigenvalues must be finite")
  }
  
  if (any(eigenvalues < 0)) {
    stop("PCA eigenvalues cannot be negative")
  }
  
  if (
    length(elbow_search) != 1L ||
    !is.finite(elbow_search) ||
    elbow_search < 1L ||
    elbow_search != as.integer(elbow_search)
  ) {
    stop("elbow_search must be one positive integer")
  }
  
  if (
    length(minimum_pcs) != 1L ||
    !is.finite(minimum_pcs) ||
    minimum_pcs < 1L ||
    minimum_pcs != as.integer(minimum_pcs)
  ) {
    stop("minimum_pcs must be one positive integer")
  }
  
  number_available <- length(eigenvalues)
  
  elbow_limit <- min(
    number_available - 1L,
    as.integer(elbow_search)
  )
  
  if (elbow_limit >= 1L) {
    current_eigenvalue <- pmax(
      eigenvalues[seq_len(elbow_limit)],
      .Machine$double.eps
    )
    
    next_eigenvalue <- pmax(
      eigenvalues[seq_len(elbow_limit) + 1L],
      .Machine$double.eps
    )
    
    log_drop <- (
      log(current_eigenvalue) -
        log(next_eigenvalue)
    )
    
    elbow_pc <- which.max(log_drop)
    
    number_selected <- min(
      number_available,
      max(as.integer(minimum_pcs), elbow_pc)
    )
  } else {
    log_drop <- numeric()
    elbow_pc <- 1L
    number_selected <- 1L
  }
  
  selected_scores <- scores[
    ,
    seq_len(number_selected),
    drop = FALSE
  ]
  
  diagnostics <- data.frame(
    pc = seq_len(number_available),
    eigenvalue = eigenvalues,
    selected = seq_len(number_available) <= number_selected,
    log_drop_to_next = NA_real_,
    row.names = NULL
  )
  
  if (length(log_drop)) {
    diagnostics$log_drop_to_next[
      seq_along(log_drop)
    ] <- log_drop
  }
  
  list(
    scores = selected_scores,
    number_selected = number_selected,
    elbow_pc = elbow_pc,
    diagnostics = diagnostics
  )
}

fit_gmm_candidates <- function(
    scores,
    max_clusters = 10L,
    min_cluster_n = 10L,
    min_cluster_prop = 0.05
) {
  if (!requireNamespace("mclust", quietly = TRUE)) {
    stop("Package 'mclust' is required")
  }
  
  scores <- as.matrix(scores)
  
  if (nrow(scores) < 2L || ncol(scores) < 1L) {
    stop("scores must contain at least two cells and one PC")
  }
  
  if (anyNA(scores) || any(!is.finite(scores))) {
    stop("scores must contain only finite values")
  }
  
  if (
    length(max_clusters) != 1L ||
    !is.finite(max_clusters) ||
    max_clusters < 1L ||
    max_clusters != as.integer(max_clusters)
  ) {
    stop("max_clusters must be one positive integer")
  }
  
  if (
    length(min_cluster_n) != 1L ||
    !is.finite(min_cluster_n) ||
    min_cluster_n < 1L ||
    min_cluster_n != as.integer(min_cluster_n)
  ) {
    stop("min_cluster_n must be one positive integer")
  }
  
  if (
    length(min_cluster_prop) != 1L ||
    !is.finite(min_cluster_prop) ||
    min_cluster_prop <= 0 ||
    min_cluster_prop > 1
  ) {
    stop("min_cluster_prop must be between zero and one")
  }
  
  number_of_cells <- nrow(scores)
  
  minimum_cluster_size <- max(
    as.integer(min_cluster_n),
    ceiling(min_cluster_prop * number_of_cells)
  )
  
  maximum_feasible_clusters <- max(
    1L,
    floor(number_of_cells / minimum_cluster_size)
  )
  
  maximum_fitted_clusters <- min(
    as.integer(max_clusters),
    maximum_feasible_clusters
  )
  
  cluster_numbers <- seq_len(maximum_fitted_clusters)
  
  fits <- vector(
    mode = "list",
    length = maximum_fitted_clusters
  )
  
  names(fits) <- as.character(cluster_numbers)
  
  bic <- rep(-Inf, maximum_fitted_clusters)
  names(bic) <- as.character(cluster_numbers)
  
  errors <- rep(NA_character_, maximum_fitted_clusters)
  names(errors) <- as.character(cluster_numbers)
  
  # Required by the internal evaluation used by mclust::Mclust.
  mclustBIC <- mclust::mclustBIC
  
  for (g in cluster_numbers) {
    fit <- tryCatch(
      mclust::Mclust(
        data = scores,
        G = g,
        modelNames = "VVI",
        verbose = FALSE
      ),
      error = function(error) {
        errors[g] <<- conditionMessage(error)
        NULL
      }
    )
    
    if (
      !is.null(fit) &&
      length(fit$bic) &&
      is.finite(as.numeric(fit$bic))
    ) {
      fits[[g]] <- fit
      bic[g] <- as.numeric(fit$bic)
    } else {
      fits[[g]] <- NULL
      
      if (is.na(errors[g])) {
        errors[g] <- "No valid model was returned"
      }
    }
  }
  
  if (!any(is.finite(bic))) {
    stop("All Gaussian mixture candidate fits failed")
  }
  
  best_g <- which.max(bic)
  
  list(
    fits = fits,
    bic = bic,
    errors = errors,
    best_g = best_g,
    best_fit = fits[[best_g]],
    minimum_cluster_size = minimum_cluster_size,
    maximum_fitted_clusters = maximum_fitted_clusters
  )
}

evaluate_gmm_structure <- function(
    gmm_candidates,
    observation_names = NULL,
    bic_threshold = 10
) {
  fits <- gmm_candidates$fits
  bic <- as.numeric(gmm_candidates$bic)
  names(bic) <- names(gmm_candidates$bic)
  
  if (!length(fits)) {
    stop("gmm_candidates contains no fitted resolutions")
  }
  
  if (!is.finite(bic[1L]) || is.null(fits[[1L]])) {
    stop("A valid one-cluster baseline model is required")
  }
  
  number_of_cells <- length(
    fits[[1L]]$classification
  )
  
  if (is.null(observation_names)) {
    observation_names <- names(
      fits[[1L]]$classification
    )
  }
  
  if (is.null(observation_names)) {
    observation_names <- paste0(
      "cell_",
      seq_len(number_of_cells)
    )
  }
  
  if (length(observation_names) != number_of_cells) {
    stop("observation_names must have one value per cell")
  }
  
  assignments <- data.frame(
    observation = observation_names,
    stringsAsFactors = FALSE,
    row.names = observation_names
  )
  
  rows <- vector("list", length(fits))
  
  for (g in seq_along(fits)) {
    fit <- fits[[g]]
    fit_successful <- !is.null(fit) && is.finite(bic[g])
    
    if (!fit_successful) {
      assignments[[paste0("cluster_g", g)]] <- NA_integer_
      
      rows[[g]] <- data.frame(
        g = g,
        fit_successful = FALSE,
        bic = bic[g],
        bic_delta_from_one = NA_real_,
        bic_delta_from_previous = NA_real_,
        number_of_clusters_returned = NA_integer_,
        minimum_cluster_size = NA_integer_,
        cluster_sizes = NA_character_,
        structurally_supported = FALSE,
        stringsAsFactors = FALSE
      )
      
      next
    }
    
    labels <- factor(fit$classification)
    sizes <- table(labels)
    
    assignments[[paste0("cluster_g", g)]] <- labels
    
    bic_delta_from_one <- if (g == 1L) {
      0
    } else {
      bic[g] - bic[1L]
    }
    
    bic_delta_from_previous <- if (
      g == 1L ||
      !is.finite(bic[g - 1L])
    ) {
      NA_real_
    } else {
      bic[g] - bic[g - 1L]
    }
    
    structurally_supported <- if (g == 1L) {
      TRUE
    } else {
      bic_delta_from_one >= bic_threshold &&
        all(
          sizes >=
            gmm_candidates$minimum_cluster_size
        )
    }
    
    rows[[g]] <- data.frame(
      g = g,
      fit_successful = TRUE,
      bic = bic[g],
      bic_delta_from_one = bic_delta_from_one,
      bic_delta_from_previous = bic_delta_from_previous,
      number_of_clusters_returned = length(sizes),
      minimum_cluster_size = min(as.integer(sizes)),
      cluster_sizes = paste(
        as.integer(sizes),
        collapse = ","
      ),
      structurally_supported = structurally_supported,
      stringsAsFactors = FALSE
    )
  }
  
  diagnostics <- do.call(rbind, rows)
  rownames(diagnostics) <- NULL
  
  list(
    diagnostics = diagnostics,
    assignments = assignments,
    minimum_allowed_cluster_size =
      gmm_candidates$minimum_cluster_size
  )
}


estimate_resolution_stability <- function(
    scores,
    full_labels,
    g,
    n_repeats = 50L,
    subsample_fraction = 0.80,
    seed = 20260809L
) {
  if (!requireNamespace("mclust", quietly = TRUE)) {
    stop("Package 'mclust' is required")
  }
  
  scores <- as.matrix(scores)
  full_labels <- factor(full_labels)
  
  if (nrow(scores) != length(full_labels)) {
    stop("full_labels must have one value per score row")
  }
  
  if (anyNA(full_labels)) {
    stop("full_labels cannot contain missing values")
  }
  
  if (
    length(g) != 1L ||
    !is.finite(g) ||
    g < 1L ||
    g != as.integer(g)
  ) {
    stop("g must be one positive integer")
  }
  
  if (nlevels(full_labels) != g) {
    stop("The number of full-label levels does not equal g")
  }
  
  if (
    length(n_repeats) != 1L ||
    !is.finite(n_repeats) ||
    n_repeats < 1L ||
    n_repeats != as.integer(n_repeats)
  ) {
    stop("n_repeats must be one positive integer")
  }
  
  if (
    length(subsample_fraction) != 1L ||
    !is.finite(subsample_fraction) ||
    subsample_fraction <= 0 ||
    subsample_fraction > 1
  ) {
    stop("subsample_fraction must be between zero and one")
  }
  
  if (g == 1L) {
    ari <- rep(1, n_repeats)
    
    return(
      list(
        ari = ari,
        median_ari = 1,
        successful_repeats = n_repeats,
        failed_repeats = 0L,
        subsample_size = nrow(scores)
      )
    )
  }
  
  number_of_cells <- nrow(scores)
  
  subsample_size <- max(
    g * 3L,
    floor(subsample_fraction * number_of_cells)
  )
  
  subsample_size <- min(
    number_of_cells,
    subsample_size
  )
  
  set.seed(as.integer(seed) + 101L)
  
  repeat_seeds <- sample.int(
    .Machine$integer.max,
    n_repeats
  )
  
  fit_diagonal_gmm <- function(data, g) {
    mclustBIC <- mclust::mclustBIC
    
    mclust::Mclust(
      data = data,
      G = g,
      modelNames = "VVI",
      verbose = FALSE
    )
  }
  
  ari <- vapply(
    repeat_seeds,
    function(repeat_seed) {
      set.seed(repeat_seed)
      
      retained_index <- sort(
        sample.int(
          number_of_cells,
          size = subsample_size,
          replace = FALSE
        )
      )
      
      subset_fit <- tryCatch(
        fit_diagonal_gmm(
          data = scores[
            retained_index,
            ,
            drop = FALSE
          ],
          g = g
        ),
        error = function(error) NULL
      )
      
      if (is.null(subset_fit)) {
        return(NA_real_)
      }
      
      mclust::adjustedRandIndex(
        full_labels[retained_index],
        subset_fit$classification
      )
    },
    numeric(1L)
  )
  
  successful <- sum(is.finite(ari))
  
  median_ari <- if (successful > 0L) {
    stats::median(ari, na.rm = TRUE)
  } else {
    NA_real_
  }
  
  list(
    ari = ari,
    median_ari = median_ari,
    successful_repeats = successful,
    failed_repeats = n_repeats - successful,
    subsample_size = subsample_size,
    repeat_seeds = repeat_seeds
  )
}

cramers_v <- function(x, y) {
  if (length(x) != length(y)) {
    stop("x and y must have the same length")
  }
  
  complete <- !is.na(x) & !is.na(y)
  x <- droplevels(factor(x[complete]))
  y <- droplevels(factor(y[complete]))
  
  if (!length(x)) {
    stop("No complete observations remain")
  }
  
  contingency_table <- table(x, y)
  
  if (
    nrow(contingency_table) < 2L ||
    ncol(contingency_table) < 2L
  ) {
    return(0)
  }
  
  chi_squared <- suppressWarnings(
    stats::chisq.test(
      contingency_table,
      correct = FALSE
    )$statistic
  )
  
  denominator <- sum(contingency_table) *
    min(
      nrow(contingency_table) - 1L,
      ncol(contingency_table) - 1L
    )
  
  if (
    !is.finite(chi_squared) ||
    denominator <= 0
  ) {
    return(0)
  }
  
  as.numeric(
    sqrt(chi_squared / denominator)
  )
}

eta_squared <- function(values, groups) {
  if (length(values) != length(groups)) {
    stop("values and groups must have the same length")
  }
  
  complete <- (
    !is.na(values) &
      is.finite(values) &
      !is.na(groups)
  )
  
  values <- as.numeric(values[complete])
  groups <- droplevels(factor(groups[complete]))
  
  if (!length(values)) {
    stop("No complete observations remain")
  }
  
  if (nlevels(groups) < 2L) {
    return(0)
  }
  
  overall_mean <- mean(values)
  
  total_sum_squares <- sum(
    (values - overall_mean)^2
  )
  
  if (
    !is.finite(total_sum_squares) ||
    total_sum_squares <= .Machine$double.eps
  ) {
    return(0)
  }
  
  group_means <- tapply(
    values,
    groups,
    mean
  )
  
  group_sizes <- table(groups)
  
  between_sum_squares <- sum(
    as.numeric(
      group_sizes[names(group_means)]
    ) *
      (group_means - overall_mean)^2
  )
  
  as.numeric(
    between_sum_squares /
      total_sum_squares
  )
}

omega_squared <- function(values, groups) {
  if (length(values) != length(groups)) {
    stop("values and groups must have the same length")
  }
  
  complete <- (
    !is.na(values) &
      is.finite(values) &
      !is.na(groups)
  )
  
  values <- as.numeric(values[complete])
  groups <- droplevels(factor(groups[complete]))
  
  number_of_cells <- length(values)
  number_of_groups <- nlevels(groups)
  
  if (!number_of_cells) {
    stop("No complete observations remain")
  }
  
  if (
    number_of_groups < 2L ||
    number_of_cells <= number_of_groups
  ) {
    return(0)
  }
  
  overall_mean <- mean(values)
  
  group_means <- tapply(
    values,
    groups,
    mean
  )
  
  group_sizes <- table(groups)
  
  between_sum_squares <- sum(
    as.numeric(
      group_sizes[names(group_means)]
    ) *
      (group_means - overall_mean)^2
  )
  
  fitted_group_mean <- group_means[
    as.character(groups)
  ]
  
  within_sum_squares <- sum(
    (values - fitted_group_mean)^2
  )
  
  total_sum_squares <- (
    between_sum_squares +
      within_sum_squares
  )
  
  within_degrees_of_freedom <- (
    number_of_cells -
      number_of_groups
  )
  
  within_mean_square <- (
    within_sum_squares /
      within_degrees_of_freedom
  )
  
  numerator <- (
    between_sum_squares -
      (number_of_groups - 1L) *
      within_mean_square
  )
  
  denominator <- (
    total_sum_squares +
      within_mean_square
  )
  
  if (
    !is.finite(numerator) ||
    !is.finite(denominator) ||
    denominator <= 0
  ) {
    return(0)
  }
  
  as.numeric(
    max(0, numerator / denominator)
  )
}

evaluate_resolution_associations <- function(
    cluster_labels,
    known_group,
    stress_score,
    size_factor,
    gating_categorical = list(),
    gating_continuous = list(),
    diagnostic_continuous = list(),
    cramer_threshold = 0.20,
    eta_threshold = 0.05,
    size_factor_eta_threshold = 0.05,
    borderline_margin = 0.005
) {
  number_of_cells <- length(cluster_labels)
  cluster_labels <- factor(cluster_labels)
  
  if (
    length(known_group) != number_of_cells ||
    length(stress_score) != number_of_cells ||
    length(size_factor) != number_of_cells
  ) {
    stop("Every nuisance variable must have one value per cell")
  }
  
  if (
    anyNA(size_factor) ||
    any(!is.finite(size_factor)) ||
    any(size_factor <= 0)
  ) {
    stop("size_factor must contain finite positive values")
  }
  
  validate_covariate_list <- function(x, argument_name) {
    if (!is.list(x)) {
      stop(argument_name, " must be a list")
    }
    
    if (
      length(x) &&
      (
        is.null(names(x)) ||
        any(!nzchar(names(x)))
      )
    ) {
      stop(argument_name, " must be a named list")
    }
    
    wrong_length <- vapply(
      x,
      length,
      integer(1L)
    ) != number_of_cells
    
    if (any(wrong_length)) {
      stop(
        argument_name,
        " contains variables with incorrect lengths"
      )
    }
  }
  
  validate_covariate_list(
    gating_categorical,
    "gating_categorical"
  )
  
  validate_covariate_list(
    gating_continuous,
    "gating_continuous"
  )
  
  validate_covariate_list(
    diagnostic_continuous,
    "diagnostic_continuous"
  )
  
  experiment_v <- cramers_v(
    cluster_labels,
    known_group
  )
  
  stress_eta <- eta_squared(
    stress_score,
    cluster_labels
  )
  
  stress_omega <- omega_squared(
    stress_score,
    cluster_labels
  )
  
  log_size_factor <- log(size_factor)
  
  size_eta <- eta_squared(
    log_size_factor,
    cluster_labels
  )
  
  size_omega <- omega_squared(
    log_size_factor,
    cluster_labels
  )
  
  metrics <- c(
    cramers_v_experiment = experiment_v,
    eta_squared_stress = stress_eta,
    omega_squared_stress = stress_omega,
    eta_squared_log_size_factor = size_eta,
    omega_squared_log_size_factor = size_omega
  )
  
  categorical_passes <- logical()
  
  for (variable_name in names(gating_categorical)) {
    safe_name <- make.names(variable_name)
    
    value <- cramers_v(
      cluster_labels,
      gating_categorical[[variable_name]]
    )
    
    metrics[
      paste0("cramers_v_", safe_name)
    ] <- value
    
    categorical_passes <- c(
      categorical_passes,
      value <= cramer_threshold
    )
  }
  
  continuous_passes <- logical()
  continuous_gating_values <- c(
    stress_eta,
    size_eta
  )
  
  continuous_gating_thresholds <- c(
    eta_threshold,
    size_factor_eta_threshold
  )
  
  for (variable_name in names(gating_continuous)) {
    safe_name <- make.names(variable_name)
    
    eta <- eta_squared(
      gating_continuous[[variable_name]],
      cluster_labels
    )
    
    omega <- omega_squared(
      gating_continuous[[variable_name]],
      cluster_labels
    )
    
    metrics[
      paste0("eta_squared_", safe_name)
    ] <- eta
    
    metrics[
      paste0("omega_squared_", safe_name)
    ] <- omega
    
    continuous_passes <- c(
      continuous_passes,
      eta <= eta_threshold
    )
    
    continuous_gating_values <- c(
      continuous_gating_values,
      eta
    )
    
    continuous_gating_thresholds <- c(
      continuous_gating_thresholds,
      eta_threshold
    )
  }
  
  for (variable_name in names(diagnostic_continuous)) {
    safe_name <- make.names(variable_name)
    
    eta <- eta_squared(
      diagnostic_continuous[[variable_name]],
      cluster_labels
    )
    
    omega <- omega_squared(
      diagnostic_continuous[[variable_name]],
      cluster_labels
    )
    
    metrics[
      paste0("eta_squared_", safe_name)
    ] <- eta
    
    metrics[
      paste0("omega_squared_", safe_name)
    ] <- omega
  }
  
  core_nuisance_clear <- (
    experiment_v <= cramer_threshold &&
      stress_eta <= eta_threshold &&
      size_eta <= size_factor_eta_threshold
  )
  
  gating_qc_clear <- (
    all(categorical_passes) &&
      all(continuous_passes)
  )
  
  borderline <- any(
    abs(
      continuous_gating_values -
        continuous_gating_thresholds
    ) <= borderline_margin
  )
  
  list(
    metrics = metrics,
    core_nuisance_clear = core_nuisance_clear,
    gating_qc_clear = gating_qc_clear,
    borderline = borderline,
    association_clear =
      core_nuisance_clear && gating_qc_clear
  )
}


evaluate_one_resolution <- function(
    g,
    scores,
    gmm_candidates,
    structure_row,
    association_data,
    settings
) {
  fit <- gmm_candidates$fits[[g]]
  
  if (is.null(fit)) {
    stop("No successful GMM fit exists for G=", g)
  }
  
  labels <- factor(fit$classification)
  
  if (
    g == 1L ||
    isTRUE(structure_row$structurally_supported)
  ) {
    stability <- estimate_resolution_stability(
      scores = scores,
      full_labels = labels,
      g = g,
      n_repeats = settings$n_stability,
      subsample_fraction =
        settings$subsample_fraction,
      seed = settings$seed
    )
  } else {
    stability <- list(
      ari = rep(
        NA_real_,
        settings$n_stability
      ),
      median_ari = NA_real_,
      successful_repeats = 0L,
      failed_repeats =
        settings$n_stability
    )
  }
  
  stable <- if (g == 1L) {
    TRUE
  } else {
    is.finite(stability$median_ari) &&
      stability$median_ari >=
      settings$stability_threshold
  }
  
  association <-
    evaluate_resolution_associations(
      cluster_labels = labels,
      known_group =
        association_data$known_group,
      stress_score =
        association_data$stress_score,
      size_factor =
        association_data$size_factor,
      gating_categorical =
        association_data$gating_categorical,
      gating_continuous =
        association_data$gating_continuous,
      diagnostic_continuous =
        association_data$diagnostic_continuous,
      cramer_threshold =
        settings$cramer_threshold,
      eta_threshold =
        settings$eta_threshold,
      size_factor_eta_threshold =
        settings$size_factor_eta_threshold,
      borderline_margin =
        settings$borderline_margin
    )
  
  primary_eligible <- (
    g > 1L &&
      isTRUE(
        structure_row$structurally_supported
      ) &&
      stable &&
      association$association_clear
  )
  
  if (g == 1L) {
    status <- "homogeneous-baseline"
  } else if (primary_eligible) {
    status <- "primary-eligible"
  } else {
    failures <- character()
    
    if (
      !isTRUE(
        structure_row$structurally_supported
      )
    ) {
      failures <- c(
        failures,
        "structurally-unsupported"
      )
    }
    
    if (!stable) {
      failures <- c(
        failures,
        "unstable"
      )
    }
    
    if (!association$core_nuisance_clear) {
      failures <- c(
        failures,
        "nuisance-associated"
      )
    }
    
    if (!association$gating_qc_clear) {
      failures <- c(
        failures,
        "qc-sensitive"
      )
    }
    
    status <- paste(
      failures,
      collapse = "+"
    )
  }
  
  diagnostic_row <- cbind(
    structure_row,
    data.frame(
      stability_median =
        stability$median_ari,
      stability_successful =
        stability$successful_repeats,
      stability_failed =
        stability$failed_repeats,
      stable = stable,
      core_nuisance_clear =
        association$core_nuisance_clear,
      gating_qc_clear =
        association$gating_qc_clear,
      association_clear =
        association$association_clear,
      borderline =
        association$borderline,
      primary_eligible =
        primary_eligible,
      status = status,
      stringsAsFactors = FALSE
    ),
    as.data.frame(
      as.list(association$metrics),
      check.names = FALSE
    )
  )
  
  list(
    diagnostic_row = diagnostic_row,
    labels = labels,
    stability_ari = stability$ari
  )
}

build_resolution_ladder <- function(
    scores,
    gmm_candidates,
    known_group,
    stress_score,
    size_factor,
    gating_categorical = list(),
    gating_continuous = list(),
    diagnostic_continuous = list(),
    bic_threshold = 10,
    stability_threshold = 0.75,
    n_stability = 50L,
    subsample_fraction = 0.80,
    cramer_threshold = 0.20,
    eta_threshold = 0.05,
    size_factor_eta_threshold = 0.05,
    borderline_margin = 0.005,
    seed = 20260809L
) {
  structure <- evaluate_gmm_structure(
    gmm_candidates = gmm_candidates,
    observation_names = rownames(scores),
    bic_threshold = bic_threshold
  )
  
  association_data <- list(
    known_group = known_group,
    stress_score = stress_score,
    size_factor = size_factor,
    gating_categorical =
      gating_categorical,
    gating_continuous =
      gating_continuous,
    diagnostic_continuous =
      diagnostic_continuous
  )
  
  settings <- list(
    stability_threshold =
      stability_threshold,
    n_stability = n_stability,
    subsample_fraction =
      subsample_fraction,
    cramer_threshold =
      cramer_threshold,
    eta_threshold = eta_threshold,
    size_factor_eta_threshold =
      size_factor_eta_threshold,
    borderline_margin =
      borderline_margin,
    seed = seed
  )
  
  successful <- structure$diagnostics[
    structure$diagnostics$fit_successful,
    ,
    drop = FALSE
  ]
  
  evaluated <- lapply(
    seq_len(nrow(successful)),
    function(index) {
      structure_row <- successful[
        index,
        ,
        drop = FALSE
      ]
      
      evaluate_one_resolution(
        g = structure_row$g,
        scores = scores,
        gmm_candidates =
          gmm_candidates,
        structure_row = structure_row,
        association_data =
          association_data,
        settings = settings
      )
    }
  )
  
  diagnostics <- do.call(
    rbind,
    lapply(
      evaluated,
      function(result) {
        result$diagnostic_row
      }
    )
  )
  
  rownames(diagnostics) <- NULL
  
  eligible <- which(
    diagnostics$primary_eligible
  )
  
  if (!length(eligible)) {
    primary_g <- 1L
    most_stable_g <- 1L
  } else {
    primary_g <- diagnostics$g[
      eligible[
        which.max(
          diagnostics$bic[eligible]
        )
      ]
    ]
    
    most_stable_g <- diagnostics$g[
      eligible[
        which.max(
          diagnostics$stability_median[
            eligible
          ]
        )
      ]
    ]
  }
  
  labels_by_g <- setNames(
    lapply(
      evaluated,
      function(result) result$labels
    ),
    diagnostics$g
  )
  
  stability_by_g <- setNames(
    lapply(
      evaluated,
      function(result) {
        result$stability_ari
      }
    ),
    diagnostics$g
  )
  
  list(
    diagnostics = diagnostics,
    assignments =
      structure$assignments,
    labels_by_g = labels_by_g,
    stability_by_g = stability_by_g,
    primary_g = as.integer(primary_g),
    most_stable_g =
      as.integer(most_stable_g),
    primary_labels =
      labels_by_g[[as.character(primary_g)]],
    failed_fits =
      structure$diagnostics[
        !structure$diagnostics$fit_successful,
        ,
        drop = FALSE
      ]
  )
}


rank_geometric_drivers <- function(
    residuals,
    cluster_labels
) {
  residuals <- as.matrix(residuals)
  cluster_labels <- droplevels(
    factor(cluster_labels)
  )
  
  if (ncol(residuals) != length(cluster_labels)) {
    stop("cluster_labels must have one value per cell")
  }
  
  if (nlevels(cluster_labels) < 2L) {
    stop("At least two clusters are required")
  }
  
  if (is.null(rownames(residuals))) {
    stop("residuals must have feature row names")
  }
  
  if (any(!is.finite(residuals))) {
    stop("residuals must contain only finite values")
  }
  
  number_of_cells <- ncol(residuals)
  overall_mean <- rowMeans(residuals)
  
  between_cluster_variation <- rep(
    0,
    nrow(residuals)
  )
  
  for (cluster in levels(cluster_labels)) {
    in_cluster <- cluster_labels == cluster
    cluster_size <- sum(in_cluster)
    
    cluster_mean <- rowMeans(
      residuals[
        ,
        in_cluster,
        drop = FALSE
      ]
    )
    
    between_cluster_variation <-
      between_cluster_variation +
      cluster_size *
      (cluster_mean - overall_mean)^2
  }
  
  total_variation <- pmax(
    0,
    rowSums(residuals * residuals) -
      number_of_cells * overall_mean^2
  )
  
  driver_eta_squared <- ifelse(
    total_variation > 0,
    between_cluster_variation /
      total_variation,
    0
  )
  
  total_centroid_separation <- sum(
    between_cluster_variation
  )
  
  centroid_contribution <- if (
    total_centroid_separation > 0
  ) {
    between_cluster_variation /
      total_centroid_separation
  } else {
    rep(0, nrow(residuals))
  }
  
  drivers <- data.frame(
    feature = rownames(residuals),
    between_cluster_variation =
      between_cluster_variation,
    driver_eta_squared =
      driver_eta_squared,
    centroid_contribution =
      centroid_contribution,
    stringsAsFactors = FALSE
  )
  
  drivers <- drivers[
    order(
      -drivers$centroid_contribution,
      -drivers$driver_eta_squared,
      drivers$feature
    ),
    ,
    drop = FALSE
  ]
  
  drivers$cumulative_contribution <- cumsum(
    drivers$centroid_contribution
  )
  
  rownames(drivers) <- NULL
  
  drivers
}


summarize_features_by_cluster <- function(
    counts,
    normalized,
    residuals,
    cluster_labels
) {
  cluster_labels <- droplevels(
    factor(cluster_labels)
  )
  
  if (
    ncol(counts) != length(cluster_labels) ||
    ncol(normalized) != length(cluster_labels) ||
    ncol(residuals) != length(cluster_labels)
  ) {
    stop("All matrices must contain one column per cluster label")
  }
  
  features <- rownames(residuals)
  
  if (
    !all(features %in% rownames(counts)) ||
    !all(features %in% rownames(normalized))
  ) {
    stop("Residual features must occur in counts and normalized")
  }
  
  counts <- counts[
    features,
    ,
    drop = FALSE
  ]
  
  normalized <- normalized[
    features,
    ,
    drop = FALSE
  ]
  
  rows <- lapply(
    levels(cluster_labels),
    function(cluster) {
      in_cluster <- cluster_labels == cluster
      cluster_size <- sum(in_cluster)
      
      cluster_counts <- counts[
        ,
        in_cluster,
        drop = FALSE
      ]
      
      cluster_normalized <- normalized[
        ,
        in_cluster,
        drop = FALSE
      ]
      
      detected_cells <- Matrix::rowSums(
        cluster_counts > 0
      )
      
      raw_count_sum <- Matrix::rowSums(
        cluster_counts
      )
      
      normalized_sum <- Matrix::rowSums(
        cluster_normalized
      )
      
      data.frame(
        feature = features,
        cluster = cluster,
        cells = cluster_size,
        
        detected_cells =
          as.integer(detected_cells),
        
        detection_fraction =
          as.numeric(detected_cells) /
          cluster_size,
        
        positive_umi_mean = ifelse(
          detected_cells > 0,
          raw_count_sum / detected_cells,
          NA_real_
        ),
        
        positive_normalized_mean = ifelse(
          detected_cells > 0,
          normalized_sum / detected_cells,
          NA_real_
        ),
        
        overall_normalized_mean =
          normalized_sum / cluster_size,
        
        residual_mean = rowMeans(
          residuals[
            ,
            in_cluster,
            drop = FALSE
          ]
        ),
        
        stringsAsFactors = FALSE
      )
    }
  )
  
  summaries <- do.call(rbind, rows)
  rownames(summaries) <- NULL
  
  summaries
}

nested_logistic_test <- function(
    detected,
    nuisance_design,
    cluster_design
) {
  detected <- as.integer(detected)
  nuisance_design <- as.matrix(
    nuisance_design
  )
  cluster_design <- as.matrix(
    cluster_design
  )
  
  number_of_cells <- length(detected)
  
  if (
    nrow(nuisance_design) != number_of_cells ||
    nrow(cluster_design) != number_of_cells
  ) {
    stop("Design matrices must have one row per cell")
  }
  
  if (
    anyNA(detected) ||
    any(!detected %in% c(0L, 1L))
  ) {
    stop("detected must contain only zero and one")
  }
  
  if (length(unique(detected)) < 2L) {
    return(NA_real_)
  }
  
  full_design <- cbind(
    nuisance_design,
    cluster_design
  )
  
  reduced_fit <- suppressWarnings(
    stats::glm.fit(
      x = nuisance_design,
      y = detected,
      family = stats::binomial()
    )
  )
  
  full_fit <- suppressWarnings(
    stats::glm.fit(
      x = full_design,
      y = detected,
      family = stats::binomial()
    )
  )
  
  deviance_difference <- (
    reduced_fit$deviance -
      full_fit$deviance
  )
  
  degrees_of_freedom <- (
    full_fit$rank -
      reduced_fit$rank
  )
  
  if (
    !is.finite(deviance_difference) ||
    degrees_of_freedom <= 0L
  ) {
    return(NA_real_)
  }
  
  stats::pchisq(
    max(0, deviance_difference),
    df = degrees_of_freedom,
    lower.tail = FALSE
  )
}



nested_positive_test <- function(
    values,
    detected,
    nuisance_design,
    cluster_design
) {
  values <- as.numeric(values)
  detected <- as.logical(detected)
  nuisance_design <- as.matrix(
    nuisance_design
  )
  cluster_design <- as.matrix(
    cluster_design
  )
  
  number_of_cells <- length(values)
  
  if (length(detected) != number_of_cells) {
    stop("detected and values must have the same length")
  }
  
  if (
    nrow(nuisance_design) != number_of_cells ||
    nrow(cluster_design) != number_of_cells
  ) {
    stop("Design matrices must have one row per cell")
  }
  
  if (
    anyNA(values) ||
    any(!is.finite(values)) ||
    anyNA(detected)
  ) {
    stop("values and detected cannot contain missing values")
  }
  
  positive_index <- which(detected)
  
  full_design <- cbind(
    nuisance_design,
    cluster_design
  )
  
  minimum_positive_cells <- (
    ncol(full_design) + 2L
  )
  
  if (
    length(positive_index) <
    minimum_positive_cells
  ) {
    return(NA_real_)
  }
  
  reduced_fit <- stats::lm.fit(
    x = nuisance_design[
      positive_index,
      ,
      drop = FALSE
    ],
    y = values[positive_index]
  )
  
  full_fit <- stats::lm.fit(
    x = full_design[
      positive_index,
      ,
      drop = FALSE
    ],
    y = values[positive_index]
  )
  
  numerator_degrees_of_freedom <- (
    full_fit$rank -
      reduced_fit$rank
  )
  
  denominator_degrees_of_freedom <- (
    length(positive_index) -
      full_fit$rank
  )
  
  reduced_rss <- sum(
    reduced_fit$residuals^2
  )
  
  full_rss <- sum(
    full_fit$residuals^2
  )
  
  if (
    numerator_degrees_of_freedom <= 0L ||
    denominator_degrees_of_freedom <= 0L ||
    !is.finite(full_rss) ||
    full_rss <= 0
  ) {
    return(NA_real_)
  }
  
  f_statistic <- (
    (
      reduced_rss -
        full_rss
    ) /
      numerator_degrees_of_freedom
  ) /
    (
      full_rss /
        denominator_degrees_of_freedom
    )
  
  stats::pf(
    max(0, f_statistic),
    df1 = numerator_degrees_of_freedom,
    df2 = denominator_degrees_of_freedom,
    lower.tail = FALSE
  )
}


rank_cluster_markers <- function(
    counts,
    normalized,
    features,
    cluster_labels,
    nuisance_design,
    clustering_features = character(),
    fdr_threshold = 0.05,
    min_detection_difference = 0.15,
    min_magnitude_difference = 0.25,
    verbose = TRUE
) {
  features <- unique(as.character(features))
  cluster_labels <- droplevels(
    factor(cluster_labels)
  )
  
  if (nlevels(cluster_labels) < 2L) {
    stop("At least two clusters are required")
  }
  
  if (
    ncol(counts) != length(cluster_labels) ||
    ncol(normalized) != length(cluster_labels) ||
    nrow(nuisance_design) != length(cluster_labels)
  ) {
    stop("All inputs must contain the same cells")
  }
  
  missing_features <- setdiff(
    features,
    rownames(counts)
  )
  
  if (length(missing_features)) {
    stop(
      "Features missing from counts: ",
      paste(head(missing_features), collapse = ", ")
    )
  }
  
  if (
    length(
      setdiff(features, rownames(normalized))
    )
  ) {
    stop("Features are missing from normalized")
  }
  
  cluster_design <- model.matrix(
    ~ cluster_labels
  )[
    ,
    -1L,
    drop = FALSE
  ]
  
  output <- vector(
    mode = "list",
    length = length(features)
  )
  
  for (feature_index in seq_along(features)) {
    feature <- features[feature_index]
    
    if (
      isTRUE(verbose) &&
      (
        feature_index == 1L ||
        feature_index %% 500L == 0L
      )
    ) {
      message(
        "Testing marker ",
        feature_index,
        " of ",
        length(features)
      )
    }
    
    raw_values <- as.numeric(
      counts[feature, ]
    )
    
    normalized_values <- as.numeric(
      normalized[feature, ]
    )
    
    detected <- raw_values > 0
    
    detection_rate <- tapply(
      detected,
      cluster_labels,
      mean
    )
    
    detection_difference <- diff(
      range(detection_rate)
    )
    
    positive_means <- tapply(
      normalized_values[detected],
      cluster_labels[detected],
      mean
    )
    
    positive_means <- positive_means[
      is.finite(positive_means)
    ]
    
    magnitude_difference <- if (
      length(positive_means) >= 2L
    ) {
      diff(range(positive_means))
    } else {
      0
    }
    
    output[[feature_index]] <- data.frame(
      feature = feature,
      
      detection_p = nested_logistic_test(
        detected = detected,
        nuisance_design = nuisance_design,
        cluster_design = cluster_design
      ),
      
      magnitude_p = nested_positive_test(
        values = normalized_values,
        detected = detected,
        nuisance_design = nuisance_design,
        cluster_design = cluster_design
      ),
      
      detection_difference =
        as.numeric(detection_difference),
      
      magnitude_difference =
        as.numeric(magnitude_difference),
      
      highest_detection_cluster =
        names(which.max(detection_rate)),
      
      highest_magnitude_cluster =
        if (length(positive_means)) {
          names(which.max(positive_means))
        } else {
          NA_character_
        },
      
      entered_clustering =
        feature %in% clustering_features,
      
      stringsAsFactors = FALSE
    )
  }
  
  markers <- do.call(rbind, output)
  
  markers$detection_q <- p.adjust(
    markers$detection_p,
    method = "BH"
  )
  
  markers$magnitude_q <- p.adjust(
    markers$magnitude_p,
    method = "BH"
  )
  
  detection_supported <- (
    !is.na(markers$detection_q) &
      markers$detection_q <=
      fdr_threshold &
      markers$detection_difference >=
      min_detection_difference
  )
  
  magnitude_supported <- (
    !is.na(markers$magnitude_q) &
      markers$magnitude_q <=
      fdr_threshold &
      markers$magnitude_difference >=
      min_magnitude_difference
  )
  
  markers$type <- ifelse(
    detection_supported &
      magnitude_supported,
    "both",
    ifelse(
      detection_supported,
      "on/off",
      ifelse(
        magnitude_supported,
        "magnitude",
        "unsupported"
      )
    )
  )
  
  safe_detection_q <- ifelse(
    is.na(markers$detection_q),
    1,
    markers$detection_q
  )
  
  safe_magnitude_q <- ifelse(
    is.na(markers$magnitude_q),
    1,
    markers$magnitude_q
  )
  
  markers$score <- -log10(
    pmax(
      pmin(
        safe_detection_q,
        safe_magnitude_q
      ),
      1e-300
    )
  )
  
  markers <- markers[
    order(
      markers$type == "unsupported",
      -markers$score,
      -markers$detection_difference,
      -markers$magnitude_difference,
      markers$feature
    ),
    ,
    drop = FALSE
  ]
  
  rownames(markers) <- NULL
  
  markers
}
