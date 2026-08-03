# Optional repository-level FlowJo orchestration.
#
# This file is intentionally outside R/. It may launch the unchanged Python
# exporter, but it is not installed with or required by the facspseudocolor R
# package. Users who already have event-level CSV files never source this file.

flowjo_or <- function(value, default) if (is.null(value)) default else value

flowjo_safe_name <- function(value) {
  value <- gsub("[^A-Za-z0-9._-]+", "_", flowjo_or(value, "replicate"))
  value <- sub("^_+", "", sub("_+$", "", value))
  if (nzchar(value)) value else "replicate"
}

flowjo_export_column <- function(column_names, channel, source = "raw") {
  exact <- paste0(source, "__", channel)
  if (exact %in% column_names) return(exact)
  hits <- column_names[startsWith(column_names, paste0(exact, " "))]
  if (length(hits) == 1L) return(hits[[1]])
  hits <- column_names[startsWith(column_names, exact)]
  if (length(hits) == 1L) return(hits[[1]])
  stop(
    "Could not find a unique '", source, "__", channel,
    "' column in the FlowJo export. Available: ",
    paste(column_names, collapse = ", "), call. = FALSE
  )
}

flowjo_sample_id <- function(sample_ids, fcs) {
  target <- tolower(sub("\\.fcs$", "", fcs, ignore.case = TRUE))
  normalized <- tolower(sub("\\.fcs$", "", sample_ids, ignore.case = TRUE))
  hits <- which(normalized == target | tolower(sample_ids) == tolower(fcs))
  if (!length(hits)) {
    stop(
      "No FlowJo sample matched FCS file '", fcs,
      "'. Exported sample ids: ", paste(sample_ids, collapse = ", "),
      call. = FALSE
    )
  }
  hits[[1]]
}

prepare_flowjo_csvs_external <- function(
    config,
    exporter = "python/export_flowjo_populations.py",
    verbose = TRUE
) {
  if (!inherits(config, "facs_config")) {
    config <- facspseudocolor::validate_facs_config(config)
  }
  replicates <- config$replicates
  if (is.null(replicates)) {
    stop("FlowJo orchestration currently requires a `replicates` configuration.",
         call. = FALSE)
  }
  data_dir <- config$data_dir
  if (!grepl("^/", data_dir)) {
    config_dir <- attr(config, "config_dir")
    if (is.null(config_dir)) {
      stop("Relative data paths require a file-backed configuration.",
           call. = FALSE)
    }
    data_dir <- file.path(config_dir, data_dir)
  }
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

  for (replicate in replicates) {
    flowjo <- utils::modifyList(
      flowjo_or(config$flowjo, list()), flowjo_or(replicate$flowjo, list())
    )
    if (!length(flowjo)) {
      stop("A top-level or per-replicate `flowjo` block is required.",
           call. = FALSE)
    }
    required <- c("source_dir", "workspace", "dna_source_channel",
                  "target_source_channel")
    missing <- required[vapply(flowjo[required], function(x) {
      is.null(x) || !is.character(x) || length(x) != 1L || !nzchar(x)
    }, logical(1))]
    if (length(missing)) {
      stop("Missing FlowJo setting(s) for ", replicate$label, ": ",
           paste(missing, collapse = ", "), call. = FALSE)
    }

    source_dir <- flowjo$source_dir
    workspace <- file.path(source_dir, flowjo$workspace)
    python <- flowjo_or(flowjo$python, Sys.which("python3"))
    rebuild <- isTRUE(flowjo_or(flowjo$rebuild, TRUE))
    population_map <- flowjo$populations
    if (is.null(population_map)) {
      population_map <- list(
        complete = flowjo_or(flowjo$population, "Single Cells")
      )
    }
    population_map <- population_map[
      intersect(names(population_map), names(config$suffixes))
    ]
    if (!length(population_map)) {
      stop("No configured FlowJo populations match the required file suffixes.",
           call. = FALSE)
    }

    prefixes <- vapply(replicate$samples, `[[`, character(1), "prefix")
    expected <- unlist(lapply(names(population_map), function(key) {
      file.path(data_dir, paste0(prefixes, config$suffixes[[key]]))
    }))
    if (!rebuild && all(file.exists(expected))) {
      if (verbose) message("FlowJo CSVs already exist for ", replicate$label)
      next
    }

    # Python, its packages, the workspace, and the original FCS directory are
    # required only when an export will actually be rebuilt.
    if (!nzchar(python) || !file.exists(python)) {
      stop("Python interpreter not found: ", python, call. = FALSE)
    }
    if (!file.exists(exporter)) {
      stop("Python exporter not found: ", exporter, call. = FALSE)
    }
    if (!file.exists(workspace)) {
      stop("FlowJo workspace not found: ", workspace, call. = FALSE)
    }
    dependency_status <- suppressWarnings(system2(
      python, c("-c", shQuote("import flowkit, pandas, lxml")),
      stdout = FALSE, stderr = FALSE
    ))
    if (dependency_status != 0) {
      stop(
        "The selected Python environment is missing flowkit, pandas, or lxml: ",
        python, call. = FALSE
      )
    }

    export_dir <- file.path(
      data_dir, ".flowjo_export", flowjo_safe_name(replicate$label)
    )
    dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)
    populations <- unique(unlist(population_map, use.names = FALSE))
    if (verbose) message("Running unchanged FlowJo exporter for ", replicate$label)
    status <- system2(python, c(
      shQuote(exporter), shQuote(workspace),
      "--fcs-dir", shQuote(source_dir),
      "--output-dir", shQuote(export_dir),
      "--populations", vapply(populations, shQuote, character(1))
    ))
    if (status != 0) {
      stop("FlowJo export failed for ", replicate$label, call. = FALSE)
    }
    counts_file <- file.path(export_dir, "population_counts.csv")
    if (!file.exists(counts_file)) {
      stop("FlowJo population count report was not created: ", counts_file,
           call. = FALSE)
    }
    population_counts <- utils::read.csv(counts_file, check.names = FALSE)

    for (key in names(population_map)) {
      population <- population_map[[key]]
      export_csv <- file.path(
        export_dir,
        paste0(gsub("[^A-Za-z0-9._-]+", "_", population), ".csv")
      )
      if (!file.exists(export_csv)) {
        stop("Expected FlowJo export was not created: ", export_csv,
             call. = FALSE)
      }
      combined <- utils::read.csv(export_csv, check.names = FALSE)
      if (!"sample_id" %in% names(combined)) {
        stop("FlowJo export lacks sample_id: ", export_csv,
             call. = FALSE)
      }
      count_rows <- population_counts[
        population_counts$gate_name == population, , drop = FALSE
      ]
      sample_ids <- unique(c(as.character(count_rows$sample_id),
                             as.character(combined$sample_id)))
      sample_ids <- sample_ids[nzchar(sample_ids)]
      if (!length(sample_ids)) {
        stop(
          "FlowJo did not report sample identities for population '",
          population, "'.", call. = FALSE
        )
      }
      dna_column <- target_column <- NULL
      if (nrow(combined)) {
        dna_column <- flowjo_export_column(
          names(combined), flowjo$dna_source_channel
        )
        target_column <- flowjo_export_column(
          names(combined), flowjo$target_source_channel
        )
      }
      for (sample in replicate$samples) {
        index <- flowjo_sample_id(sample_ids, sample$fcs)
        selected <- combined$sample_id == sample_ids[[index]]
        if (nrow(combined)) {
          output <- data.frame(
            event_index = if ("event_index" %in% names(combined)) {
              combined$event_index[selected]
            } else {
              seq_len(sum(selected)) - 1L
            },
            combined[[dna_column]][selected],
            combined[[target_column]][selected],
            check.names = FALSE
          )
          names(output) <- c(
            "event_index", config$dna_channel, config$target_channel
          )
        } else {
          output <- data.frame(
            event_index = integer(),
            dna = numeric(),
            target = numeric()
          )
          names(output) <- c(
            "event_index", config$dna_channel, config$target_channel
          )
        }
        output_file <- file.path(
          data_dir, paste0(sample$prefix, config$suffixes[[key]])
        )
        utils::write.csv(output, output_file, row.names = FALSE)
        if (verbose) message("Wrote ", output_file)
      }
    }
  }
  invisible(TRUE)
}

# Export exact FlowJo gate vertices to one experiment-level sidecar. This calls
# a separate geometry extractor and does not alter or replace the population
# exporter used above.
prepare_flowjo_gate_geometry_external <- function(
    config,
    output_file,
    population_key = "ph3_positive",
    extractor = "python/export_flowjo_gate_geometry.py",
    overwrite = FALSE,
    verbose = TRUE
) {
  if (!inherits(config, "facs_config")) {
    config <- facspseudocolor::validate_facs_config(config)
  }
  if (!is.character(output_file) || length(output_file) != 1L ||
      !nzchar(output_file)) {
    stop("`output_file` must be one explicit path.", call. = FALSE)
  }
  output_file <- path.expand(output_file)
  if (!grepl("^(/|[A-Za-z]:[/\\\\])", output_file)) {
    config_dir <- attr(config, "config_dir")
    if (is.null(config_dir)) {
      stop("A relative output path requires a file-backed configuration.",
           call. = FALSE)
    }
    output_file <- file.path(config_dir, output_file)
  }
  if (file.exists(output_file) && !isTRUE(overwrite)) {
    stop("Refusing to overwrite existing geometry file: ", output_file,
         call. = FALSE)
  }
  if (!file.exists(extractor)) {
    stop("FlowJo geometry extractor not found: ", extractor, call. = FALSE)
  }
  if (is.null(config$replicates)) {
    stop("FlowJo geometry extraction requires configured replicates.",
         call. = FALSE)
  }

  all_geometry <- list()
  for (replicate in config$replicates) {
    flowjo <- utils::modifyList(
      flowjo_or(config$flowjo, list()), flowjo_or(replicate$flowjo, list())
    )
    population_map <- flowjo$populations
    population <- population_map[[population_key]]
    if (!is.character(population) || length(population) != 1L ||
        !nzchar(population)) {
      stop("FlowJo population key '", population_key, "' is missing for ",
           replicate$label, ".", call. = FALSE)
    }
    required <- c("source_dir", "workspace", "python")
    missing <- required[vapply(flowjo[required], function(value) {
      is.null(value) || !is.character(value) || length(value) != 1L ||
        !nzchar(value)
    }, logical(1))]
    if (length(missing)) {
      stop("Missing FlowJo geometry setting(s) for ", replicate$label, ": ",
           paste(missing, collapse = ", "), call. = FALSE)
    }
    workspace <- file.path(flowjo$source_dir, flowjo$workspace)
    if (!file.exists(workspace)) {
      stop("FlowJo workspace not found: ", workspace, call. = FALSE)
    }
    if (!file.exists(flowjo$python)) {
      stop("Python interpreter not found: ", flowjo$python, call. = FALSE)
    }
    temporary <- tempfile(fileext = ".csv")
    on.exit(unlink(temporary), add = TRUE)
    if (verbose) {
      message("Extracting FlowJo gate geometry for ", replicate$label)
    }
    status <- system2(flowjo$python, c(
      shQuote(extractor), shQuote(workspace),
      "--fcs-dir", shQuote(flowjo$source_dir),
      "--population", shQuote(population),
      "--output", shQuote(temporary)
    ))
    if (status != 0 || !file.exists(temporary)) {
      stop("FlowJo gate-geometry extraction failed for ", replicate$label,
           call. = FALSE)
    }
    geometry <- utils::read.csv(temporary, check.names = FALSE)
    expected_channels <- c(
      flowjo$dna_source_channel, flowjo$target_source_channel
    )
    observed_channels <- unique(c(geometry$x_channel, geometry$y_channel))
    if (!setequal(observed_channels, expected_channels)) {
      stop(
        "Extracted gate dimensions do not match the configured DNA and target ",
        "channels for ", replicate$label, ". Expected ",
        paste(expected_channels, collapse = " and "), "; found ",
        paste(observed_channels, collapse = " and "), ".", call. = FALSE
      )
    }
    reversed <- geometry$x_channel == flowjo$target_source_channel &
      geometry$y_channel == flowjo$dna_source_channel
    if (any(reversed)) {
      for (suffix in c("transformed", "raw")) {
        x_name <- paste0("x_", suffix)
        y_name <- paste0("y_", suffix)
        temporary_value <- geometry[[x_name]][reversed]
        geometry[[x_name]][reversed] <- geometry[[y_name]][reversed]
        geometry[[y_name]][reversed] <- temporary_value
      }
      temporary_channel <- geometry$x_channel[reversed]
      geometry$x_channel[reversed] <- geometry$y_channel[reversed]
      geometry$y_channel[reversed] <- temporary_channel
    }
    sample_ids <- unique(geometry$sample_id)
    geometry$prefix <- NA_character_
    for (sample in replicate$samples) {
      index <- flowjo_sample_id(sample_ids, sample$fcs)
      geometry$prefix[
        geometry$sample_id == sample_ids[[index]]
      ] <- sample$prefix
    }
    geometry <- geometry[!is.na(geometry$prefix), , drop = FALSE]
    expected_prefixes <- vapply(
      replicate$samples, `[[`, character(1), "prefix"
    )
    missing_prefixes <- setdiff(expected_prefixes, unique(geometry$prefix))
    if (length(missing_prefixes)) {
      stop("Geometry is missing configured samples for ", replicate$label,
           ": ", paste(missing_prefixes, collapse = ", "), ".",
           call. = FALSE)
    }
    geometry$replicate <- replicate$label
    all_geometry[[length(all_geometry) + 1L]] <- geometry
  }
  result <- do.call(rbind, all_geometry)
  rownames(result) <- NULL
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(result, output_file, row.names = FALSE)
  invisible(normalizePath(output_file, mustWork = TRUE))
}
