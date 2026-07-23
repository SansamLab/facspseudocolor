# Input manifests, reading, and validation -----------------------------------

#' Build the sample manifest for an analysis
#'
#' Converts the configured samples or replicates into one row per sample while
#' preserving the configured order.
#'
#' @param config A validated `facs_config` object or a configuration list.
#'
#' @return A data frame describing sample, condition, replicate, prefix, and
#'   reference status.
#' @export
build_sample_manifest <- function(config) {
  if (!inherits(config, "facs_config")) config <- validate_facs_config(config)
  make_sample_manifest(samples = config$samples, replicates = config$replicates)
}

resolve_facs_directory <- function(config, data_dir = NULL) {
  directory <- if (is.null(data_dir)) config$data_dir else data_dir
  if (!config_scalar_string(directory)) {
    stop("`data_dir` must be one explicit, nonempty directory path.", call. = FALSE)
  }
  directory <- path.expand(directory)
  if (!grepl("^(/|[A-Za-z]:[/\\\\])", directory)) {
    config_dir <- attr(config, "config_dir")
    if (is.null(data_dir) && config_scalar_string(config_dir)) {
      directory <- file.path(config_dir, directory)
    } else if (is.null(data_dir)) {
      stop(
        "Relative `data_dir` requires a configuration read from a file; otherwise supply an explicit `data_dir`.",
        call. = FALSE
      )
    }
  }
  normalizePath(directory, mustWork = FALSE)
}

required_population_keys <- function(config, manifest) {
  if (identical(config$plot_type, "poi")) return("complete")
  if (identical(config$plot_type, "ph3")) {
    return(c("complete", "g1", "ph3_positive"))
  }

  require_all_positive <- isTRUE(config$show_edu_apex_line) ||
    isTRUE(config$quantify_phase_median) ||
    isTRUE(config$quantify_whole_median)
  list(
    common = c("complete", "g1"),
    positive_rows = if (require_all_positive) {
      seq_len(nrow(manifest))
    } else {
      which(manifest$is_reference)
    }
  )
}

facs_input_files <- function(config, data_dir = NULL) {
  manifest <- build_sample_manifest(config)
  directory <- resolve_facs_directory(config, data_dir)
  rows <- list()
  add_row <- function(i, key, required = TRUE) {
    suffix <- config$suffixes[[key]]
    data.frame(
      replicate = manifest$replicate[[i]],
      replicate_index = manifest$replicate_index[[i]],
      condition = manifest$condition[[i]],
      prefix = manifest$prefix[[i]],
      population = key,
      required = required,
      path = file.path(directory, paste0(manifest$prefix[[i]], suffix)),
      stringsAsFactors = FALSE
    )
  }

  if (identical(config$plot_type, "poi")) {
    rows <- lapply(seq_len(nrow(manifest)), add_row, key = "complete")
  } else if (identical(config$plot_type, "ph3")) {
    for (i in seq_len(nrow(manifest))) {
      rows[[length(rows) + 1L]] <- add_row(i, "complete")
      rows[[length(rows) + 1L]] <- add_row(i, "g1")
      rows[[length(rows) + 1L]] <- add_row(i, "ph3_positive")
    }
  } else {
    for (i in seq_len(nrow(manifest))) {
      rows[[length(rows) + 1L]] <- add_row(i, "complete")
      rows[[length(rows) + 1L]] <- add_row(i, "g1")
      positive_required <- i %in% required_population_keys(config, manifest)$positive_rows
      rows[[length(rows) + 1L]] <- add_row(i, "edu_positive", positive_required)
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Read and validate one event-level FACS table
#'
#' Reads a CSV path or validates an in-memory data frame. The function does not
#' filter, repair, transform, or substitute event values.
#'
#' @param input An explicit CSV path or an in-memory data frame.
#' @param dna_channel Exact name of the DNA channel column.
#' @param target_channel Exact name of the target channel column.
#' @param sample_id Label used in validation errors and provenance.
#' @param minimum_finite_events Minimum number of rows in which both required
#'   channels are finite.
#'
#' @return The original event table with validation and provenance attributes.
#' @export
read_facs_sample <- function(
    input, dna_channel, target_channel, sample_id = "sample",
    minimum_finite_events = 2L
) {
  if (!config_scalar_string(dna_channel) || !config_scalar_string(target_channel)) {
    stop("`dna_channel` and `target_channel` must be explicit nonempty strings.",
         call. = FALSE)
  }
  if (identical(dna_channel, target_channel)) {
    stop("DNA and target channels must be different columns.", call. = FALSE)
  }
  if (!config_scalar_number(as.numeric(minimum_finite_events)) ||
      minimum_finite_events < 0 || minimum_finite_events %% 1 != 0) {
    stop("`minimum_finite_events` must be a nonnegative integer.", call. = FALSE)
  }

  source_path <- NULL
  if (is.data.frame(input)) {
    events <- input
    input_type <- "data.frame"
  } else if (config_scalar_string(input)) {
    if (!file.exists(input)) {
      stop("Input CSV not found for ", sample_id, ": ", input, call. = FALSE)
    }
    source_path <- normalizePath(input, mustWork = TRUE)
    events <- utils::read.csv(source_path, check.names = FALSE)
    input_type <- "csv"
  } else {
    stop("`input` must be an explicit CSV path or a data frame.", call. = FALSE)
  }

  missing_columns <- setdiff(c(dna_channel, target_channel), names(events))
  if (length(missing_columns)) {
    stop(
      "Missing required channel(s) for ", sample_id, ": ",
      paste(missing_columns, collapse = ", "), ". Available columns: ",
      paste(names(events), collapse = ", "), call. = FALSE
    )
  }
  valid_required_type <- function(value) {
    is.numeric(value) || (nrow(events) == 0L && is.logical(value))
  }
  nonnumeric <- c(dna_channel, target_channel)[
    !vapply(events[c(dna_channel, target_channel)],
            valid_required_type, logical(1))
  ]
  if (length(nonnumeric)) {
    stop("Required channel(s) must be numeric for ", sample_id, ": ",
         paste(nonnumeric, collapse = ", "), call. = FALSE)
  }

  finite <- is.finite(events[[dna_channel]]) & is.finite(events[[target_channel]])
  finite_n <- sum(finite)
  if (finite_n < minimum_finite_events) {
    stop(
      "Too few finite events for ", sample_id, ": found ", finite_n,
      ", require at least ", minimum_finite_events, ".", call. = FALSE
    )
  }

  attr(events, "facs_input") <- list(
    sample_id = sample_id,
    input_type = input_type,
    source_path = source_path,
    event_n = nrow(events),
    finite_event_n = finite_n,
    nonfinite_event_n = nrow(events) - finite_n,
    dna_channel = dna_channel,
    target_channel = target_channel
  )
  events
}

#' Validate all required input CSVs for an analysis
#'
#' Resolves the configured sample files, checks that every required file exists,
#' and validates required channel columns and numeric types. Event values are
#' never filtered or repaired.
#'
#' @param config A validated `facs_config` object or configuration list.
#' @param data_dir Optional explicit data directory overriding `config$data_dir`.
#'
#' @return A data frame with one row per expected population file, including
#'   existence and event-count diagnostics.
#' @export
validate_facs_inputs <- function(config, data_dir = NULL) {
  if (!inherits(config, "facs_config")) config <- validate_facs_config(config)
  files <- facs_input_files(config, data_dir)
  files$exists <- file.exists(files$path)

  missing <- files$required & !files$exists
  if (any(missing)) {
    details <- paste0(
      files$replicate[missing], " / ", files$condition[missing], " / ",
      files$population[missing], ": ", files$path[missing]
    )
    stop(
      paste0("Missing required FACS input file(s):\n- ",
             paste(details, collapse = "\n- ")), call. = FALSE
    )
  }

  files$event_n <- NA_integer_
  files$finite_event_n <- NA_integer_
  files$nonfinite_event_n <- NA_integer_
  files$subset_membership_validated <- NA
  for (i in which(files$exists)) {
    minimum <- if (files$population[[i]] == "complete") {
      10L
    } else if (files$population[[i]] == "ph3_positive") {
      0L
    } else {
      2L
    }
    events <- read_facs_sample(
      files$path[[i]], config$dna_channel, config$target_channel,
      sample_id = paste(files$replicate[[i]], files$condition[[i]],
                        files$population[[i]], sep = " / "),
      minimum_finite_events = minimum
    )
    info <- attr(events, "facs_input")
    files$event_n[[i]] <- info$event_n
    files$finite_event_n[[i]] <- info$finite_event_n
    files$nonfinite_event_n[[i]] <- info$nonfinite_event_n
  }
  if (identical(config$plot_type, "ph3")) {
    grouped <- split(files, files$prefix)
    count_errors <- vapply(grouped, function(x) {
      complete_n <- x$event_n[x$population == "complete"]
      positive_n <- x$event_n[x$population == "ph3_positive"]
      length(complete_n) == 1L && length(positive_n) == 1L &&
        is.finite(complete_n) && is.finite(positive_n) &&
        positive_n > complete_n
    }, logical(1))
    if (any(count_errors)) {
      stop(
        "pH3-positive event count exceeds the Single Cell count for: ",
        paste(names(grouped)[count_errors], collapse = ", "),
        ". Verify that the exported pH3 gate is nested within Single Cells.",
        call. = FALSE
      )
    }
    for (prefix in names(grouped)) {
      group <- grouped[[prefix]]
      complete_path <- group$path[group$population == "complete"]
      positive_path <- group$path[group$population == "ph3_positive"]
      complete <- utils::read.csv(complete_path, check.names = FALSE)
      positive <- utils::read.csv(positive_path, check.names = FALSE)
      validated <- FALSE
      if ("event_index" %in% names(complete) &&
          "event_index" %in% names(positive)) {
        if (anyNA(complete$event_index) || anyNA(positive$event_index) ||
            anyDuplicated(complete$event_index) ||
            anyDuplicated(positive$event_index)) {
          stop(
            "Invalid or duplicate event_index values for PH3 sample ",
            prefix, ".", call. = FALSE
          )
        }
        outside <- setdiff(positive$event_index, complete$event_index)
        if (length(outside)) {
          stop(
            "The pH3-positive population contains event_index values outside ",
            "the Single Cell population for ", prefix,
            ". Verify the FlowJo gate hierarchy.", call. = FALSE
          )
        }
        validated <- TRUE
      }
      row <- files$prefix == prefix & files$population == "ph3_positive"
      files$subset_membership_validated[row] <- validated
    }
  }
  files
}
