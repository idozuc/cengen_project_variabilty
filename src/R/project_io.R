# Shared command-line and path helpers for CeNGEN workflows.

script_file <- function(fallback) {
  argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(argument)) normalizePath(sub("^--file=", "", argument[[1L]])) else normalizePath(fallback)
}

project_root_from_script <- function(path) {
  # All executable scripts live exactly three levels below the project root.
  dirname(dirname(dirname(normalizePath(path))))
}

parse_named_arguments <- function(arguments, defaults, usage) {
  if (any(arguments %in% c("--help", "-h"))) {
    cat(usage, "\n")
    quit(save = "no", status = 0L)
  }

  config_positions <- which(arguments == "--config")
  if (length(config_positions)) {
    if (length(config_positions) != 1L || config_positions == length(arguments)) {
      stop("--config requires exactly one YAML path", call. = FALSE)
    }
    if (!requireNamespace("yaml", quietly = TRUE)) {
      stop("Package 'yaml' is required for --config", call. = FALSE)
    }
    configured <- yaml::read_yaml(arguments[[config_positions + 1L]])
    unknown <- setdiff(names(configured), names(defaults))
    if (length(unknown)) {
      stop("Unknown configuration key(s): ", paste(unknown, collapse = ", "), call. = FALSE)
    }
    defaults[names(configured)] <- configured
  }

  index <- 1L
  while (index <= length(arguments)) {
    key <- arguments[[index]]
    if (key == "--config") {
      index <- index + 2L
      next
    }
    if (!startsWith(key, "--") || index == length(arguments)) {
      stop("Arguments must be supplied as --name value pairs", call. = FALSE)
    }
    normalized <- gsub("-", "_", substring(key, 3L), fixed = TRUE)
    if (!normalized %in% names(defaults)) stop("Unknown argument: ", key, call. = FALSE)
    defaults[[normalized]] <- arguments[[index + 1L]]
    index <- index + 2L
  }
  defaults
}

require_path_option <- function(value, label, must_exist = TRUE) {
  if (is.null(value) || length(value) != 1L || is.na(value) || !nzchar(value)) {
    stop(label, " is required", call. = FALSE)
  }
  if (must_exist && !file.exists(value) && !dir.exists(value)) {
    stop(label, " does not exist: ", value, call. = FALSE)
  }
  value
}

basics_result_dir <- function(run_dir, task_id, cell_type) {
  file.path(run_dir, "cell_types", sprintf("%03d_%s", as.integer(task_id), cell_type))
}

write_resolved_yaml <- function(value, path) {
  if (!requireNamespace("yaml", quietly = TRUE)) stop("Package 'yaml' is required", call. = FALSE)
  yaml::write_yaml(value, path)
}
