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

flowjo_physical_candidate_path <- function(path) {
  # First collapse lexical `.`/`..`; then resolve the nearest existing ancestor
  # physically so a symlink in that ancestor cannot evade containment checks.
  path <- normalizePath(path.expand(path), mustWork = FALSE)
  tail <- character()
  candidate <- path
  while (!file.exists(candidate)) {
    parent <- dirname(candidate)
    if (identical(parent, candidate)) {
      stop("Cannot resolve an existing ancestor for path: ", path, call. = FALSE)
    }
    tail <- c(basename(candidate), tail)
    candidate <- parent
  }
  # Destinations may be either directories (for exported event tables) or
  # files (for the configured panel PDF/PNG).  An existing output file is a
  # valid leaf: resolve it physically for the containment check rather than
  # treating it as an invalid ancestor.
  if (!dir.exists(candidate)) {
    if (!length(tail)) return(normalizePath(candidate, mustWork = TRUE))
    stop("Path has a non-directory existing ancestor: ", candidate, call. = FALSE)
  }
  do.call(file.path, c(list(normalizePath(candidate, mustWork = TRUE)), as.list(tail)))
}

flowjo_path_is_within <- function(path, parent) {
  path <- flowjo_physical_candidate_path(path)
  parent <- normalizePath(path.expand(parent), mustWork = TRUE)
  identical(path, parent) || startsWith(path, paste0(parent, .Platform$file.sep))
}

flowjo_config_path <- function(path, config, label, allow_command = FALSE) {
  if (!is.character(path) || length(path) != 1L || !nzchar(path)) {
    stop(label, " must be one nonempty path.", call. = FALSE)
  }
  if (allow_command && !grepl("[/\\\\]", path)) {
    resolved <- Sys.which(path)
    if (!nzchar(resolved)) {
      stop("Python interpreter not found on PATH: ", path, call. = FALSE)
    }
    return(normalizePath(resolved, mustWork = TRUE))
  }
  path <- path.expand(path)
  if (grepl("^(/|[A-Za-z]:[/\\\\])", path)) return(path)
  config_dir <- attr(config, "config_dir")
  if (!is.character(config_dir) || length(config_dir) != 1L || !nzchar(config_dir)) {
    stop("Relative ", label, " requires a file-backed configuration.", call. = FALSE)
  }
  file.path(config_dir, path)
}

resolve_poi_flowjo_config_paths <- function(config) {
  if (!inherits(config, "facs_config")) {
    config <- facspseudocolor::validate_facs_config(config)
  }
  config$data_dir <- flowjo_config_path(config$data_dir, config, "`data_dir`")
  config$output_pdf <- flowjo_config_path(config$output_pdf, config, "`output_pdf`")
  config$output_png <- flowjo_config_path(config$output_png, config, "`output_png`")
  resolve_flowjo <- function(flowjo) {
    if (!length(flowjo)) return(flowjo)
    if (!is.null(flowjo$source_dir)) {
      flowjo$source_dir <- flowjo_config_path(flowjo$source_dir, config, "FlowJo `source_dir`")
    }
    if (!is.null(flowjo$python)) {
      flowjo$python <- flowjo_config_path(flowjo$python, config, "FlowJo `python`", allow_command = TRUE)
    }
    flowjo
  }
  config$flowjo <- resolve_flowjo(flowjo_or(config$flowjo, list()))
  if (!is.null(config$replicates)) {
    for (index in seq_along(config$replicates)) {
      config$replicates[[index]]$flowjo <- resolve_flowjo(
        flowjo_or(config$replicates[[index]]$flowjo, list())
      )
    }
  }
  config
}

flowjo_require_external_destination <- function(path, label, source_dir, package_root) {
  if (!is.character(path) || length(path) != 1L || !grepl("^/", path)) {
    stop(label, " must be one explicit absolute path outside the source and package directories.",
         call. = FALSE)
  }
  if (flowjo_path_is_within(path, source_dir) ||
      flowjo_path_is_within(path, package_root)) {
    stop(label, " must be outside the FlowJo source directory and package repository.",
         call. = FALSE)
  }
  invisible(normalizePath(path.expand(path), mustWork = FALSE))
}

preflight_poi_flowjo_destinations_external <- function(config,
                                                        exporter = "python/export_flowjo_populations.py") {
  if (!inherits(config, "facs_config")) {
    config <- facspseudocolor::validate_facs_config(config)
  }
  config <- resolve_poi_flowjo_config_paths(config)
  if (!identical(config$plot_type, "poi") || is.null(config$replicates) ||
      length(config$replicates) != 1L) {
    stop("POI FlowJo destination preflight requires exactly one POI replicate.",
         call. = FALSE)
  }
  if (!file.exists(exporter)) {
    stop("Python exporter not found: ", exporter, call. = FALSE)
  }
  flowjo <- utils::modifyList(
    flowjo_or(config$flowjo, list()), flowjo_or(config$replicates[[1L]]$flowjo, list())
  )
  source_dir <- flowjo$source_dir
  if (!is.character(source_dir) || length(source_dir) != 1L || !dir.exists(source_dir)) {
    stop("FlowJo `source_dir` must be one existing local directory.", call. = FALSE)
  }
  package_root <- normalizePath(file.path(dirname(exporter), ".."), mustWork = TRUE)
  flowjo_require_external_destination(config$data_dir, "`data_dir`", source_dir, package_root)
  flowjo_require_external_destination(config$output_pdf, "`output_pdf`", source_dir, package_root)
  flowjo_require_external_destination(config$output_png, "`output_png`", source_dir, package_root)
  invisible(TRUE)
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

prepare_edu_g1_inputs_external <- function(
    config,
    frozen_model_caller = "python/apply_frozen_gate_models.py",
    verbose = TRUE
) {
  if (!inherits(config, "facs_config")) {
    config <- facspseudocolor::validate_facs_config(config)
  }
  if (!identical(config$plot_type, "edu")) {
    stop("G1-source orchestration is available only for EdU configurations.",
         call. = FALSE)
  }
  route <- if (identical(config$gating$mode, "flowjo")) {
    "flowjo"
  } else if (identical(config$gating$profile, "documented_edu_g1_v1")) {
    "documented_edu_g1_v1"
  } else if (identical(config$gating$profile, "single_g1_edu_frozen_v1")) {
    "single_g1_edu_frozen_v1"
  } else {
    stop("Validated EdU gating profile is unsupported by orchestration.", call. = FALSE)
  }
  if (identical(route, "flowjo")) {
    return(prepare_flowjo_csvs_external(config, verbose = verbose))
  }
  flowjo <- config$flowjo
  if (identical(route, "documented_edu_g1_v1")) {
    manifest <- path.expand(config$gating$manifest)
    if (!grepl("^/", manifest)) {
      config_path <- attr(config, "config_path")
      if (!is.character(config_path) || length(config_path) != 1L) {
        stop("Relative documented-profile manifest requires a file-backed config.",
             call. = FALSE)
      }
      manifest <- file.path(dirname(config_path), manifest)
    }
    if (!file.exists(manifest)) {
      stop("Documented model profile requires its existing exported manifest; run the pinned operation separately before analysis.",
           call. = FALSE)
    }
    return(invisible(normalizePath(manifest, mustWork = TRUE)))
  }
  if (identical(route, "single_g1_edu_frozen_v1")) {
    python <- flowjo_or(flowjo$python, Sys.which("python3"))
    if (!nzchar(python) || !file.exists(python)) {
      stop("Python interpreter not found: ", python, call. = FALSE)
    }
    config_path <- attr(config, "config_path")
    if (!is.character(config_path) || length(config_path) != 1L ||
        !file.exists(config_path) || !file.exists(frozen_model_caller)) {
      stop("Frozen model gating requires its file-backed config and caller.",
           call. = FALSE)
    }
    if (verbose) message("Running immutable frozen-profile EdU population caller")
    status <- system2(python, c(shQuote(frozen_model_caller), shQuote(config_path)))
    if (status != 0L) stop("Frozen-profile EdU population calling failed.", call. = FALSE)
    output_dir <- config$gating$output_dir
    if (!grepl("^/", output_dir)) {
      output_dir <- file.path(dirname(config_path), output_dir)
    }
    manifest <- file.path(output_dir, "model-gating-manifest.json")
    if (!file.exists(manifest)) {
      stop("Frozen-profile model-gating manifest was not created.", call. = FALSE)
    }
    return(invisible(normalizePath(manifest, mustWork = TRUE)))
  }
}

prepare_flowjo_csvs_external <- function(
    config,
    exporter = "python/export_flowjo_populations.py",
    verbose = TRUE
) {
  if (!inherits(config, "facs_config")) {
    config <- facspseudocolor::validate_facs_config(config)
  }
  config <- resolve_poi_flowjo_config_paths(config)
  replicates <- config$replicates
  if (is.null(replicates)) {
    stop("FlowJo orchestration currently requires a `replicates` configuration.",
         call. = FALSE)
  }
  # The full pH3-grade production-identity contract (contract_metadata /
  # export_operation_id / an explicit direct-index-semantics attestation) was
  # built for the pH3 pipeline's provenance needs and should apply only
  # there. Every other plot_type uses the exporter's existing, lighter-weight
  # legacy profile -- the same shape this orchestration supported before the
  # pH3 contract was added (see docs/CONFIGURATION.md's "Optional FlowJo
  # block").
  strict <- identical(config$plot_type, "ph3")
  export_profile <- if (strict) "production_direct_identity_v1" else "legacy_count_only_unverified_v1"
  data_dir <- config$data_dir
  operation_artifacts <- character()

  for (replicate in replicates) {
    flowjo <- utils::modifyList(
      flowjo_or(config$flowjo, list()), flowjo_or(replicate$flowjo, list())
    )
    if (!length(flowjo)) {
      stop("A top-level or per-replicate `flowjo` block is required.",
           call. = FALSE)
    }
    minimal <- !strict && is.character(flowjo$contract_metadata) &&
      length(flowjo$contract_metadata) == 1L && nzchar(flowjo$contract_metadata)
    required <- c("source_dir", "workspace", "dna_source_channel",
                  "target_source_channel")
    if (strict) {
      required <- c(required, "contract_metadata", "export_operation_id")
    } else if (minimal) {
      required <- c(required, "contract_metadata")
    }
    missing <- required[vapply(flowjo[required], function(x) {
      is.null(x) || !is.character(x) || length(x) != 1L || !nzchar(x)
    }, logical(1))]
    if (length(missing)) {
      stop("Missing FlowJo setting(s) for ", replicate$label, ": ",
           paste(missing, collapse = ", "), call. = FALSE)
    }

    if (minimal && (!identical(flowjo$dna_source_channel, config$dna_channel) ||
        !identical(flowjo$target_source_channel, config$target_channel))) {
      stop("Minimal FlowJo export raw detector channels must exactly match `dna_channel` and `target_channel`.", call. = FALSE)
    }
    source_dir <- flowjo$source_dir
    workspace <- file.path(source_dir, flowjo$workspace)
    contract_metadata <- if (strict) {
      file.path(source_dir, flowjo$contract_metadata)
    } else if (minimal) {
      flowjo_config_path(flowjo$contract_metadata, config, "FlowJo `contract_metadata`")
    } else {
      NULL
    }
    if (strict && !isTRUE(flowjo$direct_index_semantics_verified)) {
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
    if (minimal) {
      fcs_files <- vapply(replicate$samples, `[[`, character(1), "fcs")
      if (any(!grepl("^[^/\\\\]+\\.fcs$", fcs_files, ignore.case = TRUE)) || anyDuplicated(fcs_files)) {
        stop("Minimal FlowJo export requires unique basename-only `.fcs` sample names.", call. = FALSE)
      }
      sample_population_plan <- stats::setNames(rep(list(names(population_map)), length(fcs_files)), fcs_files)
      plan_file <- tempfile("flowjo_sample_population_plan_", fileext = ".json")
      on.exit(unlink(plan_file), add = TRUE)
      jsonlite::write_json(sample_population_plan, plan_file, auto_unbox = FALSE, pretty = TRUE)
    }
    expected <- unlist(lapply(names(population_map), function(key) {
      file.path(data_dir, paste0(prefixes, config$suffixes[[key]]))
    }))
    if (minimal && !rebuild) {
      stop("Minimal FlowJo orchestration always creates a new immutable operation; `rebuild: false` is unsupported.", call. = FALSE)
    }
    if (!minimal && !rebuild && all(file.exists(expected))) {
      stop(
        "FlowJo orchestration cannot bypass manifest and artifact ",
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
    if ((strict || minimal) && !file.exists(contract_metadata)) {
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

    package_root <- normalizePath(file.path(dirname(exporter), ".."), mustWork = TRUE)
    if (minimal) flowjo_require_external_destination(data_dir, "`data_dir`", source_dir, package_root)
    export_dir <- if (minimal) file.path(data_dir, ".flowjo_export") else file.path(data_dir, ".flowjo_export", flowjo_safe_name(replicate$label))
    if (!minimal) dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)
    populations <- unlist(population_map, use.names = FALSE)
    if (verbose) message("Running contract-aware FlowJo exporter for ", replicate$label)
    exporter_args <- c(
      shQuote(exporter), shQuote(workspace),
      "--fcs-dir", shQuote(source_dir),
      "--output-dir", shQuote(export_dir),
      "--populations", vapply(populations, shQuote, character(1)),
      "--population-keys", vapply(names(population_map), shQuote, character(1)),
      "--output-suffixes", vapply(
        unname(unlist(config$suffixes[names(population_map)])),
        shQuote, character(1)
      ),
      "--profile", if (minimal) "minimal_provenance_v1" else export_profile
    )
    if (minimal) {
      exporter_args <- c(exporter_args, "--sample-population-plan", shQuote(plan_file),
                         "--raw-compatible-columns", "--contract-metadata", shQuote(contract_metadata),
                         "--analysis-channels", shQuote(flowjo$dna_source_channel), shQuote(flowjo$target_source_channel))
    } else if (strict) {
      exporter_args <- c(
        exporter_args,
        "--contract-metadata", shQuote(contract_metadata),
        "--export-operation-id", shQuote(flowjo$export_operation_id),
        "--direct-index-semantics-verified"
      )
    }
    exporter_output <- if (minimal) system2(python, exporter_args, stdout = TRUE, stderr = TRUE) else system2(python, exporter_args)
    status <- if (minimal) flowjo_or(attr(exporter_output, "status"), 0L) else exporter_output
    if (status != 0) {
      stop("FlowJo export failed for ", replicate$label, call. = FALSE)
    }
    if (minimal) {
      operation_lines <- grep("^FLOWJO_OPERATION_DIR=", exporter_output, value = TRUE)
      if (length(operation_lines) != 1L) stop("FlowJo exporter did not return exactly one finalized operation directory.", call. = FALSE)
      export_dir <- sub("^FLOWJO_OPERATION_DIR=", "", operation_lines[[1L]])
      if (!grepl("^/", export_dir) || !flowjo_path_is_within(export_dir, file.path(data_dir, ".flowjo_export"))) {
        stop("FlowJo exporter returned an operation directory outside the approved export root.", call. = FALSE)
      }
      export_dir <- normalizePath(export_dir, mustWork = TRUE)
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
        expected_profile <- if (minimal) "minimal_provenance_v1" else export_profile
        if (any(exported$export_profile != expected_profile)) {
          if (strict) {
            stop("Legacy or ambiguous FlowJo exports cannot be consumed as production.",
                 call. = FALSE)
          }
          stop("Unexpected FlowJo export profile; expected the legacy export profile.",
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

# Point a one-replicate POI analysis at the immutable, verified population
# artifacts just created by `prepare_flowjo_csvs_external()`.  The config is
# copied only in memory; it neither moves/copies artifacts nor changes the
# source config, FCS files, or FlowJo workspace.
prepare_poi_flowjo_analysis_config_external <- function(config, artifacts) {
  if (!inherits(config, "facs_config")) {
    config <- facspseudocolor::validate_facs_config(config)
  }
  if (!identical(config$plot_type, "poi")) {
    stop("The FlowJo POI report requires `plot_type: poi`.", call. = FALSE)
  }
  if (is.null(config$replicates) || length(config$replicates) != 1L) {
    stop(
      "The FlowJo POI report supports exactly one configured replicate/export ",
      "operation. Split multiple biological replicates into separate, explicit ",
      "configurations and reports.", call. = FALSE
    )
  }
  artifacts <- normalizePath(artifacts, mustWork = TRUE)
  expected <- file.path(
    dirname(artifacts[[1L]]),
    paste0(
      vapply(config$replicates[[1L]]$samples, `[[`, character(1), "prefix"),
      config$suffixes$complete
    )
  )
  expected <- normalizePath(expected, mustWork = FALSE)
  if (length(artifacts) != length(expected) ||
      !setequal(artifacts, expected)) {
    stop(
      "FlowJo export artifacts do not exactly match the configured POI Single ",
      "Cells inputs; refusing to select a partial or substituted operation.",
      call. = FALSE
    )
  }
  operation_dir <- dirname(artifacts[[1L]])
  if (!all(dirname(artifacts) == operation_dir) ||
      !file.exists(file.path(operation_dir, "export-manifest.json")) ||
      !file.exists(file.path(operation_dir, "export-manifest.sha256"))) {
    stop("Verified FlowJo operation artifacts must share one finalized operation directory.",
         call. = FALSE)
  }
  analysis_config <- config
  analysis_config$data_dir <- operation_dir
  analysis_config
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
