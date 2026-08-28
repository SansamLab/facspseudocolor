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
  operation_artifacts <- character()

  for (replicate in replicates) {
    flowjo <- utils::modifyList(
      flowjo_or(config$flowjo, list()), flowjo_or(replicate$flowjo, list())
    )
    if (!length(flowjo)) {
      stop("A top-level or per-replicate `flowjo` block is required.",
           call. = FALSE)
    }
    required <- c("source_dir", "workspace", "dna_source_channel",
                  "target_source_channel", "contract_metadata",
                  "export_operation_id")
    missing <- required[vapply(flowjo[required], function(x) {
      is.null(x) || !is.character(x) || length(x) != 1L || !nzchar(x)
    }, logical(1))]
    if (length(missing)) {
      stop("Missing FlowJo setting(s) for ", replicate$label, ": ",
           paste(missing, collapse = ", "), call. = FALSE)
    }

    source_dir <- flowjo$source_dir
    workspace <- file.path(source_dir, flowjo$workspace)
    contract_metadata <- file.path(source_dir, flowjo$contract_metadata)
    if (!isTRUE(flowjo$direct_index_semantics_verified)) {
      stop("Production FlowJo orchestration requires the pinned SYNTHETIC direct-index verification.",
           call. = FALSE)
    }
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
      stop(
        "Production FlowJo orchestration cannot bypass manifest and artifact ",
        "verification with `rebuild: false`; validate the completed operation ",
        "explicitly or request a new operation ID/output directory.", call. = FALSE
      )
    }

    # Python, its packages, the workspace, and the original FCS directory are
    # required only when an export will actually be rebuilt.
    if (!nzchar(python) || !file.exists(python)) {
      stop("Python interpreter not found: ", python, call. = FALSE)
    }
    if (!file.exists(exporter)) {
      stop("Python exporter not found: ", exporter, call. = FALSE)
    }
    contract_verifier <- file.path(dirname(exporter), "export_contract.py")
    if (!file.exists(contract_verifier)) {
      stop("Export-contract verifier not found: ", contract_verifier, call. = FALSE)
    }
    if (!file.exists(workspace)) {
      stop("FlowJo workspace not found: ", workspace, call. = FALSE)
    }
    if (!file.exists(contract_metadata)) {
      stop("FlowJo contract metadata not found: ", contract_metadata, call. = FALSE)
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
    populations <- unlist(population_map, use.names = FALSE)
    if (verbose) message("Running contract-aware FlowJo exporter for ", replicate$label)
    status <- system2(python, c(
      shQuote(exporter), shQuote(workspace),
      "--fcs-dir", shQuote(source_dir),
      "--output-dir", shQuote(export_dir),
      "--populations", vapply(populations, shQuote, character(1)),
      "--population-keys", vapply(names(population_map), shQuote, character(1)),
      "--output-suffixes", vapply(
        unname(unlist(config$suffixes[names(population_map)])),
        shQuote, character(1)
      ),
      "--profile", "production_direct_identity_v1",
      "--contract-metadata", shQuote(contract_metadata),
      "--export-operation-id", shQuote(flowjo$export_operation_id),
      "--direct-index-semantics-verified"
    ))
    if (status != 0) {
      stop("FlowJo export failed for ", replicate$label, call. = FALSE)
    }
    counts_file <- file.path(export_dir, "population_counts.csv")
    if (!file.exists(counts_file)) {
      stop("FlowJo population count report was not created: ", counts_file,
           call. = FALSE)
    }
    replicate_artifacts <- character()
    for (key in names(population_map)) {
      identity_columns <- c(
        "acquisition_id", "event_index", "event_identity", "identity_source",
        "identity_method_id", "identity_method_version", "duplicate_occurrence",
        "export_profile", "export_operation_id", "export_manifest_digest",
        "export_manifest_reference"
      )
      for (sample in replicate$samples) {
        artifact <- file.path(
          export_dir, paste0(sample$prefix, config$suffixes[[key]])
        )
        if (!file.exists(artifact)) {
          stop("Expected immutable per-acquisition artifact is missing: ", artifact,
               call. = FALSE)
        }
        exported <- utils::read.csv(
          artifact, check.names = FALSE,
          colClasses = stats::setNames(rep("character", length(identity_columns)),
                                       identity_columns)
        )
        missing_identity <- setdiff(identity_columns, names(exported))
        if (length(missing_identity)) {
          stop("Sequential identity fallback is prohibited; required identity fields are missing.",
               call. = FALSE)
        }
        if (any(exported$export_profile != "production_direct_identity_v1")) {
          stop("Legacy or ambiguous FlowJo exports cannot be consumed as production.",
               call. = FALSE)
        }
        operation_artifacts <- c(operation_artifacts, artifact)
        replicate_artifacts <- c(replicate_artifacts, artifact)
      }
    }
    verification_status <- system2(python, c(
      shQuote(contract_verifier),
      "--verify-operation", shQuote(export_dir),
      "--artifacts", vapply(replicate_artifacts, shQuote, character(1))
    ))
    if (verification_status != 0) {
      stop("Finalized manifest or consumed population artifact verification failed for ",
           replicate$label, ".", call. = FALSE)
    }
  }
  invisible(normalizePath(operation_artifacts, mustWork = TRUE))
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

  if (length(config$replicates) != 1L) {
    stop(
      "Verified geometry orchestration is one export operation at a time; ",
      "call it with a single-replicate configuration for each operation.",
      call. = FALSE
    )
  }
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
    required <- c("source_dir", "workspace", "python", "export_operation_id")
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
    data_dir <- config$data_dir
    if (!grepl("^/", data_dir)) {
      config_dir <- attr(config, "config_dir")
      if (is.null(config_dir)) {
        stop("Relative data paths require a file-backed configuration.", call. = FALSE)
      }
      data_dir <- file.path(config_dir, data_dir)
    }
    operation_dir <- file.path(
      data_dir, ".flowjo_export", flowjo_safe_name(replicate$label)
    )
    expected_output <- normalizePath(
      file.path(operation_dir, "gate_geometry.csv"), mustWork = FALSE
    )
    if (!identical(normalizePath(output_file, mustWork = FALSE), expected_output)) {
      stop("Verified geometry output must be the operation artifact: ",
           expected_output, call. = FALSE)
    }
    if (verbose) {
      message("Extracting FlowJo gate geometry for ", replicate$label)
    }
    status <- system2(flowjo$python, c(
      shQuote(extractor), shQuote(workspace),
      "--fcs-dir", shQuote(flowjo$source_dir),
      "--population", shQuote(population),
      "--population-key", shQuote(population_key),
      "--operation-dir", shQuote(operation_dir),
      "--export-operation-id", shQuote(flowjo$export_operation_id),
      "--output", shQuote(output_file)
    ))
    if (status != 0 || !file.exists(output_file) ||
        !file.exists(file.path(operation_dir, "geometry-manifest.json"))) {
      stop("FlowJo gate-geometry extraction failed for ", replicate$label,
           call. = FALSE)
    }
    geometry <- utils::read.csv(output_file, check.names = FALSE)
    expected_channels <- c(
      flowjo$dna_source_channel, flowjo$target_source_channel
    )
    observed_pairs <- unique(geometry[c("x_channel", "y_channel")])
    if (nrow(observed_pairs) != 1L ||
        !identical(as.character(unlist(observed_pairs[1, ], use.names = FALSE)),
                   expected_channels)) {
      stop(
        "Extracted gate dimensions do not match the configured DNA and target ",
        "channels for ", replicate$label, ". Expected ",
        paste(expected_channels, collapse = " and "), "; found ",
        paste(unlist(observed_pairs[1, ], use.names = FALSE), collapse = " and "),
        ".", call. = FALSE
      )
    }
    expected_prefixes <- vapply(
      replicate$samples, `[[`, character(1), "prefix"
    )
    missing_prefixes <- setdiff(expected_prefixes, unique(geometry$prefix))
    if (length(missing_prefixes)) {
      stop("Geometry is missing configured samples for ", replicate$label,
           ": ", paste(missing_prefixes, collapse = ", "), ".",
           call. = FALSE)
    }
  }
  invisible(normalizePath(output_file, mustWork = TRUE))
}
