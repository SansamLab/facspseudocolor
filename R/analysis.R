# Experiment orchestration and result objects --------------------------------

analysis_settings <- function(config) {
  list(
    dna_channel = config$dna_channel,
    target_channel = config$target_channel,
    normalization_method = if (config$plot_type == "edu") {
      "reference_negative_regression"
    } else if (config$plot_type == "ph3") {
      "g1_dna_only"
    } else {
      "background_reference_regression"
    },
    normalize_target = config$normalize_target,
    y_log10 = config$y_log10,
    dna_2n_value = config$dna_2n_value,
    g1_anchor = config$g1_anchor,
    x_limits = as.numeric(unlist(config$x_limits)),
    y_limits = if (is.null(config$y_limits)) NULL else
      as.numeric(unlist(config$y_limits)),
    point_size = config$point_size,
    density_bandwidth = config$density_bandwidth,
    density_lower_clip = config$density_lower_clip,
    density_upper_clip = config$density_upper_clip,
    density_gamma = config$density_gamma,
    baseline_fit_x_range = as.numeric(unlist(config$baseline_fit_x_range)),
    baseline_boundary_bins = config$baseline_boundary_bins,
    baseline_minimum_events_per_bin = config$baseline_minimum_events_per_bin,
    baseline_minimum_negative_events = config$baseline_minimum_negative_events,
    background_quantile = config$background_quantile,
    show_cutoff_line = config$show_cutoff_line,
    poi_dna_align = config$poi_dna_align,
    poi_peak_failure = config$poi_peak_failure,
    y_limit_lower_quantile = config$y_limit_lower_quantile,
    y_limit_upper_quantile = config$y_limit_upper_quantile,
    show_edu_apex_line = config$show_edu_apex_line,
    edu_apex_x_range = as.numeric(unlist(config$edu_apex_x_range)),
    edu_apex_density_adjust = config$edu_apex_density_adjust
  )
}

ph3_normalization_inputs <- function(input_report, prefix, profile) {
  if (identical(profile, "production_direct_identity_v1")) {
    verified <- attr(input_report, "ph3_verified_events")[[prefix]]
    if (is.null(verified) ||
        !all(c("complete", "g1", "ph3_positive") %in% names(verified))) {
      stop("Verified PH3 event tables are missing for ", prefix, ".",
           call. = FALSE)
    }
    return(verified[c("complete", "g1", "ph3_positive")])
  }
  input_path <- function(population) {
    paths <- input_report$path[input_report$prefix == prefix &
                                 input_report$population == population]
    if (length(paths) != 1L) {
      stop("Validated PH3 input path is ambiguous for ", prefix, " / ",
           population, ".", call. = FALSE)
    }
    paths[[1L]]
  }
  list(complete = input_path("complete"), g1 = input_path("g1"),
       ph3_positive = if (identical(profile, "legacy_csv_pilot_v1")) {
         NULL
       } else {
         input_path("ph3_positive")
       })
}

new_facs_analysis <- function(
    config, manifest, input_report, normalized_data, models,
    quantitation = list(), warnings = character()
) {
  ph3_mode <- identical(config$plot_type, "ph3")
  containment <- if (ph3_mode) attr(input_report, "ph3_containment") else NULL
  export_manifests <- if (ph3_mode) {
    attr(input_report, "ph3_export_manifests")
  } else NULL
  attr(input_report, "ph3_verified_events") <- NULL
  provenance <- list(
    package = "facspseudocolor",
    package_version = as.character(
      utils::packageVersion("facspseudocolor")
    ),
    r_version = R.version.string,
    config_path = attr(config, "config_path"),
    config_dir = attr(config, "config_dir"),
    input_files = input_report$path[input_report$exists]
  )
  if (ph3_mode) {
    provenance$ph3_containment <- containment
    provenance$ph3_export_manifests <- export_manifests
  }
  structure(
    list(
      artifact_version = 1L,
      config = config,
      sample_manifest = manifest,
      input_report = input_report,
      normalized_data = normalized_data,
      models = models,
      quantitation = quantitation,
      warnings = unique(warnings),
      provenance = provenance
    ),
    class = "facs_analysis"
  )
}

#' Analyze a configured FACS experiment from CSV inputs
#'
#' Validates all required inputs, fits the configured replicate-level reference
#' models, and normalizes every sample. This function performs no plotting and
#' writes no files. For direct in-memory calculations, use [normalize_edu()],
#' [normalize_poi()], or [normalize_ph3()].
#'
#' @param config A validated `facs_config`, an unvalidated configuration list,
#'   or an explicit path to a YAML configuration.
#' @param data_dir Optional explicit data directory overriding `config$data_dir`.
#'
#' @return A structured object of class `facs_analysis` containing configuration,
#'   sample metadata, input diagnostics, normalized event tables, fitted models,
#'   quantitation placeholders, warnings, and provenance.
#' @export
analyze_facs_experiment <- function(config, data_dir = NULL) {
  if (config_scalar_string(config)) config <- read_facs_config(config)
  if (!inherits(config, "facs_config")) config <- validate_facs_config(config)

  manifest <- build_sample_manifest(config)
  directory <- resolve_facs_directory(config, data_dir)
  input_report <- validate_facs_inputs(config, directory)
  settings <- analysis_settings(config)

  model_gate_inputs <- NULL
  model_gate_manifest_path <- NULL
  model_verified_frames <- NULL
  if (identical(config$plot_type, "edu") &&
      identical(config$gating$mode, "model_experimental") &&
      identical(config$gating$profile, "documented_edu_g1_v1")) {
    model_gate_manifest_path <- path.expand(config$gating$manifest)
    if (!grepl("^(/|[A-Za-z]:[/\\\\])", model_gate_manifest_path)) {
      config_dir <- attr(config, "config_dir")
      if (!config_scalar_string(config_dir)) {
        stop("Relative documented-model manifest requires a file-backed configuration.",
             call. = FALSE)
      }
      model_gate_manifest_path <- file.path(config_dir, model_gate_manifest_path)
    }
    provisional <- structure(list(
      config = config, sample_manifest = manifest, input_report = input_report
    ), class = "facs_analysis")
    model_gate_inputs <- validate_edu_model_gate_manifest(
      provisional, model_gate_manifest_path
    )
    model_verified_frames <- attr(model_gate_inputs, "verified_frames")
  }

  if (config$plot_type == "edu") {
    models <- fit_sample_reference_models(
      manifest, directory, config$suffixes, settings,
      verified_frames = model_verified_frames
    )
    normalized_data <- lapply(seq_len(nrow(manifest)), function(i) {
      model <- models[[manifest$prefix[[i]]]]
      sample <- read_and_normalize_sample(
        prefix = manifest$prefix[[i]],
        condition_label = manifest$condition[[i]],
        data_dir = directory,
        file_suffixes = config$suffixes,
        settings = settings,
        baseline_slope = model$slope,
        verified_frames = model_verified_frames
      )
      add_boundary <- function(data) {
        data$edu_boundary_bgsub <- rep(
          model$positive_boundary_bgsub_cutoff, nrow(data)
        )
        data$edu_boundary_raw <- data$baseline + data$edu_boundary_bgsub
        data$edu_boundary_norm <- data$edu_boundary_raw / data$baseline *
          config$dna_2n_value
        data$edu_computed_positive <- is.finite(data$target_raw) &
          is.finite(data$edu_boundary_raw) &
          data$target_raw >= data$edu_boundary_raw
        data
      }
      sample$data <- add_boundary(sample$data)
      if (!is.null(sample$edu_positive)) {
        sample$edu_positive <- add_boundary(sample$edu_positive)
      }
      sample
    })
  } else if (config$plot_type == "poi") {
    models <- fit_replicate_background_models(
      manifest, directory, config$suffixes, settings
    )
    normalized_data <- lapply(seq_len(nrow(manifest)), function(i) {
      model <- models[[manifest$model_group[[i]]]]
      read_and_normalize_sample(
        prefix = manifest$prefix[[i]],
        condition_label = manifest$condition[[i]],
        data_dir = directory,
        file_suffixes = config$suffixes,
        settings = settings,
        background_model = model
      )
    })
  } else {
    ph3_containment <- attr(input_report, "ph3_containment")
    ph3_export_manifests <- attr(input_report, "ph3_export_manifests")
    normalized_data <- lapply(seq_len(nrow(manifest)), function(i) {
      prefix <- manifest$prefix[[i]]
      label <- manifest$condition[[i]]
      inputs <- ph3_normalization_inputs(input_report, prefix,
                                         config$ph3_input_profile)
      normalized <- normalize_ph3(
        events = inputs$complete,
        g1_events = inputs$g1,
        ph3_positive_events = inputs$ph3_positive,
        dna_channel = config$dna_channel,
        target_channel = config$target_channel,
        dna_2n_value = config$dna_2n_value,
        g1_anchor = config$g1_anchor,
        sample_id = label
      )
      if (identical(config$ph3_input_profile,
                    "production_direct_identity_v1")) {
        normalized$ph3_event_classification <- build_ph3_event_classification(
          normalized = normalized, validated_inputs = inputs, config = config,
          manifest_row = manifest[i, , drop = FALSE],
          containment = ph3_containment,
          export_manifests = ph3_export_manifests
        )
      }
      normalized
    })
    models <- lapply(normalized_data, function(x) {
      list(
        normalization_method = x$normalization_method,
        g1_dna_anchor = x$g1_dna_anchor,
        dna_normalization_factor = x$dna_normalization_factor,
        g1_anchor_method = x$g1_anchor_method
      )
    })
    names(models) <- manifest$prefix
  }
  names(normalized_data) <- manifest$prefix

  analysis <- new_facs_analysis(
    config = config,
    manifest = manifest,
    input_report = input_report,
    normalized_data = normalized_data,
    models = models
  )
  if (identical(config$plot_type, "edu")) {
    analysis$provenance$gating <- if (
      identical(config$gating$mode, "model_experimental") &&
      identical(config$gating$profile, "documented_edu_g1_v1")
    ) {
      list(
        source = "model", profile = config$gating$profile,
        manifest = normalizePath(model_gate_manifest_path, mustWork = TRUE),
        verified_inputs = model_gate_inputs,
        verified_artifacts = attr(model_gate_inputs, "model_artifacts"),
        qc_warnings = attr(model_gate_inputs, "qc_warnings"),
        populations = c("Single Cells", "G1", "EdU Positive"),
        flowjo_workspace_gate_used = FALSE,
        interpretation = paste(
          "experimental model-derived research output;",
          "not expert gating or biological ground truth"
        )
      )
    } else {
      list(source = config$gating$mode, profile = config$gating$profile %||% NULL)
    }
    if (!is.null(model_gate_inputs)) {
      analysis$warnings <- unique(c(
        analysis$warnings,
        paste0("MODEL_QC_REVIEW_REQUIRED: ",
               attr(model_gate_inputs, "qc_warnings"))
      ))
    }
  }
  if (config$plot_type == "ph3" &&
      identical(config$ph3_input_profile,
                "production_direct_identity_v1")) {
    if (identical(config$ph3_positivity_method,
                  "ph3_raw_4n_density_cutoff_v1")) {
      analysis <- ph3_apply_raw_4n_density_cutoff(analysis)
    }
    analysis <- quantify_ph3_production_acquisitions(analysis)
    analysis$ph3_output_model <- derive_ph3_condition_report_model(
      derive_ph3_signal_outcomes(
        apply_ph3_background_regression(build_ph3_output_model(analysis))
      )
    )
    analysis
  } else if (config$plot_type == "ph3" &&
             identical(config$ph3_input_profile, "legacy_csv_pilot_v1")) {
    ph3_build_legacy_csv_pilot(analysis)
  } else if (config$plot_type == "ph3") {
    quantify_ph3(analysis)
  } else {
    quantify_cell_cycle(
      analysis,
      include = c("phase_median", "whole_median", "phase_percent"),
      signal = config$quant_signal,
      reference_condition = config$quant_reference_condition
    )
  }
}

#' Validate provenance for experimental model-derived EdU gate inputs
#'
#' @param analysis A completed `facs_analysis` in EdU mode.
#' @param manifest_path Explicit path to `model-gating-manifest.json`.
#' @return An invisible data frame describing the verified input artifacts.
facs_sha256_file <- function(path) {
  paste0(as.character(openssl::sha256(file(path))), collapse = "")
}

facs_read_file_bytes <- function(path, size) {
  readBin(path, what = "raw", n = size)
}

facs_read_hashed_csv <- function(path, maximum_bytes = 512 * 1024^2) {
  size <- file.info(path)$size
  if (!is.finite(size) || size <= 0 || size > maximum_bytes) {
    stop("Model-gating CSV size is invalid.", call. = FALSE)
  }
  bytes <- facs_read_file_bytes(path, size)
  if (length(bytes) != size) {
    stop("Model-gating CSV changed while it was being read.", call. = FALSE)
  }
  list(
    sha256 = paste0(openssl::sha256(bytes), collapse = ""),
    data = utils::read.csv(text = rawToChar(bytes), check.names = FALSE,
                           stringsAsFactors = FALSE)
  )
}

#' @export
validate_edu_frozen_model_gate_manifest <- function(analysis, manifest_path) {
  if (!inherits(analysis, "facs_analysis") ||
      !identical(analysis$config$plot_type, "edu")) {
    stop("`analysis` must be a completed EdU facs_analysis.", call. = FALSE)
  }
  if (!config_scalar_string(manifest_path) || !file.exists(manifest_path)) {
    stop("Model-gating manifest path must name one existing file.", call. = FALSE)
  }
  manifest_path <- normalizePath(manifest_path, mustWork = TRUE)
  manifest <- jsonlite::fromJSON(manifest_path, simplifyVector = FALSE)
  if (!identical(manifest$status_label,
                 "EXPERIMENTAL MODEL-DERIVED NON-PRODUCTION")) {
    stop("Model-gating manifest lacks the required non-production status.",
         call. = FALSE)
  }
  if (!identical(manifest$gating_mode, "model_experimental") ||
      !identical(manifest$profile, "single_g1_edu_frozen_v1") ||
      !identical(analysis$config$gating$mode, "model_experimental") ||
      !identical(analysis$config$gating$profile, manifest$profile)) {
    stop("Model-gating manifest/config mode or frozen profile differs.", call. = FALSE)
  }
  expected_thresholds <- list(single_cells = 0.205, g1 = 0.255,
                              edu_positive = 0.385)
  expected_g1_dna_minimum <- list(
    method = "median_positive_dna_area_of_preliminary_edu_negative_g1",
    fraction = 0.35
  )
  observed_g1_dna_minimum <- manifest$g1_dna_minimum
  g1_dna_policy_ok <- is.list(observed_g1_dna_minimum) &&
    length(observed_g1_dna_minimum) == 2L &&
    !anyDuplicated(names(observed_g1_dna_minimum)) &&
    setequal(names(observed_g1_dna_minimum),
             names(expected_g1_dna_minimum)) &&
    identical(observed_g1_dna_minimum$method,
              expected_g1_dna_minimum$method) &&
    is.numeric(observed_g1_dna_minimum$fraction) &&
    length(observed_g1_dna_minimum$fraction) == 1L &&
    is.finite(observed_g1_dna_minimum$fraction) &&
    abs(observed_g1_dna_minimum$fraction -
        expected_g1_dna_minimum$fraction) <= 1e-12
  if (!g1_dna_policy_ok) {
    stop("Model-gating G1 DNA minimum policy differs.", call. = FALSE)
  }
  observed_thresholds <- manifest$frozen_thresholds
  threshold_names_ok <- is.list(observed_thresholds) &&
    length(observed_thresholds) == length(expected_thresholds) &&
    !anyDuplicated(names(observed_thresholds)) &&
    setequal(names(observed_thresholds), names(expected_thresholds))
  threshold_values_ok <- threshold_names_ok && all(vapply(
    names(expected_thresholds),
    function(target) {
      observed <- observed_thresholds[[target]]
      is.numeric(observed) && length(observed) == 1L && is.finite(observed) &&
        abs(observed - expected_thresholds[[target]]) <= 1e-12
    }, logical(1)
  ))
  if (!threshold_values_ok) {
    stop("Model-gating manifest frozen thresholds differ from the approved profile.",
         call. = FALSE)
  }
  expected_predictor <- list(
    format = "PORTABLE_HGB_V1",
    implementation = "embedded dependency-free reference",
    derived_from_reference_sha256 = "2edf1bbd529159e8752d01731e87f19efbc63b514409b1d246d7a4ffc266a2c3"
  )
  observed_predictor <- manifest$edu_predictor
  predictor_names_ok <- is.list(observed_predictor) &&
    length(observed_predictor) == length(expected_predictor) &&
    !anyDuplicated(names(observed_predictor)) &&
    setequal(names(observed_predictor), names(expected_predictor))
  predictor_values_ok <- predictor_names_ok && all(vapply(
    names(expected_predictor),
    function(field) identical(observed_predictor[[field]],
                               expected_predictor[[field]]),
    logical(1)
  ))
  if (!predictor_values_ok) {
    stop("Model-gating manifest portable EdU predictor provenance differs.",
         call. = FALSE)
  }
  configured_gating <- analysis$config$gating
  configured_output <- normalizePath(configured_gating$output_dir, mustWork = TRUE)
  observed_output <- normalizePath(dirname(manifest_path), mustWork = TRUE)
  if (!identical(observed_output, configured_output)) {
    stop("Model-gating manifest directory differs from configured output_dir.",
         call. = FALSE)
  }
  acquisitions <- manifest$acquisitions
  if (!is.list(acquisitions) || !length(acquisitions)) {
    stop("Model-gating manifest has no acquisitions.", call. = FALSE)
  }
  artifacts <- manifest$artifacts
  if (!is.list(artifacts) ||
      !setequal(names(artifacts), c("single_cells", "g1", "edu_positive"))) {
    stop("Model-gating manifest must name exactly the Single Cells, G1, and EdU-positive artifacts.",
         call. = FALSE)
  }
  artifact_rows <- list()
  for (target in names(artifacts)) {
    artifact <- artifacts[[target]]
    configured_artifact <- configured_gating$artifacts[[target]]
    is_edu <- identical(target, "edu_positive")
    artifact_fields <- if (is_edu) c("path", "byte_sha256") else
      c("path", "byte_sha256", "semantic_sha256")
    if (!identical(artifact[artifact_fields], configured_artifact[artifact_fields])) {
      stop("Model-gating manifest artifact differs from config for ", target, ".",
           call. = FALSE)
    }
    if (!config_scalar_string(artifact$path) ||
        !config_scalar_string(artifact$byte_sha256) ||
        (!is_edu && (!config_scalar_string(artifact$semantic_sha256))) ||
        !file.exists(artifact$path)) {
      stop("Model artifact provenance is incomplete for ", target, ".",
           call. = FALSE)
    }
    observed_artifact_hash <- facs_sha256_file(artifact$path)
    if (!identical(observed_artifact_hash, artifact$byte_sha256)) {
      stop("Model artifact byte SHA-256 differs for ", target, ".", call. = FALSE)
    }
    if (is_edu) {
      if (!identical(observed_artifact_hash,
                     "9e6ed466687efe5f7523dea633e9e9244381c0e613c86f0944f861c06047fdf4"))
        stop("Portable EdU model artifact provenance differs from the frozen profile.", call. = FALSE)
      artifact_rows[[length(artifact_rows) + 1L]] <- data.frame(
        target = target, byte_sha256 = observed_artifact_hash,
        semantic_sha256 = NA_character_, threshold = 0.385,
        feature_schema = "EDU_POSITIVE_MODEL002/PORTABLE_HGB_V1", stringsAsFactors = FALSE)
      next
    }
    rule <- jsonlite::fromJSON(artifact$path, simplifyVector = FALSE)
    if (!config_scalar_string(rule$sha256) ||
        !config_scalar_string(rule$schema_version) ||
        !config_scalar_number(as.numeric(rule$threshold))) {
      stop("Model artifact lacks required semantic hash, schema, or threshold for ",
           target, ".", call. = FALSE)
    }
    expected_target <- c(single_cells = "single_cells_reference",
                         g1 = "g1_reference")[[target]]
    expected_threshold <- expected_thresholds[[target]]
    observed_threshold <- as.numeric(rule$threshold)
    if (!identical(rule$target, expected_target) ||
        length(observed_threshold) != 1L || !is.finite(observed_threshold) ||
        abs(observed_threshold - expected_threshold) > 1e-12) {
      stop("Model artifact target or threshold differs from the frozen profile for ",
           target, ".", call. = FALSE)
    }
    if (!identical(rule$sha256, artifact$semantic_sha256) ||
        !identical(rule$schema_version, manifest$feature_schema_version)) {
      stop("Model artifact semantic provenance differs for ", target, ".",
           call. = FALSE)
    }
    artifact_rows[[length(artifact_rows) + 1L]] <- data.frame(
      target = target, byte_sha256 = observed_artifact_hash,
      semantic_sha256 = rule$sha256, threshold = observed_threshold,
      feature_schema = rule$schema_version, stringsAsFactors = FALSE
    )
  }
  prefixes <- vapply(acquisitions, `[[`, character(1), "prefix")
  configured <- as.character(analysis$sample_manifest$prefix)
  if (anyDuplicated(prefixes) || !setequal(prefixes, configured) ||
      length(prefixes) != length(configured)) {
    stop("Model-gating manifest prefixes do not exactly match the analysis.",
         call. = FALSE)
  }
  input_paths <- normalizePath(analysis$input_report$path, mustWork = TRUE)
  input_keys <- paste(analysis$input_report$prefix,
                      analysis$input_report$population, sep = "\r")
  population_keys <- c(single_cells = "complete", g1 = "g1",
                       edu_positive = "edu_positive")
  verified <- list()
  population_frames <- list()
  all_event_frames <- list()
  for (acquisition in acquisitions) {
    configured_matches <- Filter(
      function(x) identical(x$prefix, acquisition$prefix),
      configured_gating$acquisitions
    )
    if (length(configured_matches) != 1L) {
      stop("Configured acquisition mapping is missing or ambiguous for ",
           acquisition$prefix, ".", call. = FALSE)
    }
    configured_acquisition <- configured_matches[[1L]]
    observed_roles <- acquisition$channel_roles
    configured_roles <- configured_acquisition$channel_roles
    role_names_ok <- is.list(observed_roles) && is.list(configured_roles) &&
      length(observed_roles) == length(configured_roles) &&
      !anyDuplicated(names(observed_roles)) &&
      !anyDuplicated(names(configured_roles)) &&
      setequal(names(observed_roles), names(configured_roles))
    role_values_ok <- role_names_ok && all(vapply(
      names(configured_roles),
      function(role) identical(observed_roles[[role]], configured_roles[[role]]),
      logical(1)
    ))
    if (!identical(acquisition$acquisition_id,
                   configured_acquisition$acquisition_id) ||
        !identical(normalizePath(acquisition$source_fcs, mustWork = TRUE),
                   normalizePath(configured_acquisition$fcs_path, mustWork = TRUE)) ||
        !identical(acquisition$source_fcs_sha256,
                   configured_acquisition$fcs_sha256) ||
        !role_values_ok ||
        !identical(acquisition$workspace_used, FALSE)) {
      stop("Model-gating manifest acquisition mapping differs from config for ",
           acquisition$prefix, ".", call. = FALSE)
    }
    observed_source_hash <- facs_sha256_file(configured_acquisition$fcs_path)
    if (!identical(observed_source_hash, configured_acquisition$fcs_sha256) ||
        !identical(observed_source_hash, acquisition$source_fcs_sha256)) {
      stop("Configured source FCS SHA-256 differs for ", acquisition$prefix, ".",
           call. = FALSE)
    }
    population_frames[[acquisition$prefix]] <- list()
    output_names <- names(acquisition$outputs)
    if (!setequal(output_names,
                  c("all_events", "single_cells", "g1", "edu_positive"))) {
      stop("Model-gating manifest outputs are invalid for ",
           acquisition$prefix, ".", call. = FALSE)
    }
    for (output_name in names(population_keys)) {
      key <- paste(acquisition$prefix, population_keys[[output_name]], sep = "\r")
      positions <- which(input_keys == key)
      if (length(positions) != 1L) {
        stop("Analysis input is missing or ambiguous for ", key, ".",
             call. = FALSE)
      }
      output <- acquisition$outputs[[output_name]]
      if (is.null(output$path) || is.null(output$sha256)) {
        stop("Model-gating manifest output record is incomplete for ", key, ".",
             call. = FALSE)
      }
      if (!config_scalar_string(output$path) ||
          !identical(basename(output$path), output$path) ||
          output$path %in% c(".", "..")) {
        stop("Model-gating manifest output path must be one path-free basename for ",
             key, ".", call. = FALSE)
      }
      expected_path <- normalizePath(
        file.path(dirname(manifest_path), output$path), mustWork = TRUE
      )
      if (!identical(dirname(expected_path), dirname(manifest_path))) {
        stop("Model-gating manifest output path escapes its manifest directory for ",
             key, ".", call. = FALSE)
      }
      if (!identical(input_paths[[positions]], expected_path)) {
        stop("Analysis input path differs from the model-gating manifest for ",
             key, ".", call. = FALSE)
      }
      verified_csv <- facs_read_hashed_csv(expected_path)
      observed_hash <- verified_csv$sha256
      if (!identical(observed_hash, output$sha256)) {
        stop("Analysis input SHA-256 differs from the model-gating manifest for ",
             key, ".", call. = FALSE)
      }
      output_frame <- verified_csv$data
      if (!is.numeric(output$rows) || length(output$rows) != 1L ||
          as.integer(output$rows) != nrow(output_frame)) {
        stop("Model-gating manifest row count differs for ", key, ".", call. = FALSE)
      }
      identity_columns <- c("acquisition_id", "event_identity", "event_index")
      if (!all(identity_columns %in% names(output_frame))) {
        stop("Model-gating CSV lacks direct event identity columns for ", key, ".",
             call. = FALSE)
      }
      acquisition_ids <- unique(output_frame$acquisition_id)
      identities <- output_frame$event_identity
      indices <- suppressWarnings(as.numeric(output_frame$event_index))
      expected_identities <- paste0(acquisition$acquisition_id,
                                    ":event_index:", indices)
      if (length(acquisition_ids) != 1L ||
          !identical(acquisition_ids[[1L]], acquisition$acquisition_id) ||
          anyNA(indices) || any(indices < 0) || any(indices %% 1 != 0) ||
          anyDuplicated(indices) || anyNA(identities) || any(!nzchar(identities)) ||
          anyDuplicated(identities) || !identical(identities, expected_identities)) {
        stop("Model-gating CSV lacks unique direct event identities for ", key, ".",
             call. = FALSE)
      }
      population_frames[[acquisition$prefix]][[output_name]] <- identities
      verified[[length(verified) + 1L]] <- data.frame(
        prefix = acquisition$prefix,
        population = population_keys[[output_name]], path = expected_path,
        sha256 = observed_hash, stringsAsFactors = FALSE
      )
    }
    all_output <- acquisition$outputs$all_events
    if (!is.null(all_output) &&
        (!is.list(all_output) || !config_scalar_string(all_output$path) ||
        !identical(basename(all_output$path), all_output$path) ||
        all_output$path %in% c(".", "..") ||
        !config_scalar_string(all_output$sha256))) {
      stop("Model-gating all-events output record is incomplete for ",
           acquisition$prefix, ".", call. = FALSE)
    }
    if (!is.null(all_output)) {
      all_path <- normalizePath(file.path(dirname(manifest_path), all_output$path),
                                mustWork = TRUE)
      maximum_all_events_bytes <- 512 * 1024^2
      maximum_all_events_rows <- 5000000L
      if (!identical(dirname(all_path), dirname(manifest_path)) ||
          !file.exists(all_path)) {
        stop("Model-gating all-events path or size is invalid for ",
             acquisition$prefix, ".", call. = FALSE)
      }
      all_csv <- tryCatch(
        facs_read_hashed_csv(all_path, maximum_all_events_bytes),
        error = function(e) stop(
          "Model-gating all-events path or size is invalid for ",
          acquisition$prefix, ".", call. = FALSE
        )
      )
      if (!identical(all_csv$sha256, all_output$sha256)) {
        stop("Model-gating all-events path or SHA-256 is invalid for ",
             acquisition$prefix, ".", call. = FALSE)
      }
      all_frame <- all_csv$data
      identity_columns <- c("acquisition_id", "event_identity", "event_index")
      indices <- suppressWarnings(as.numeric(all_frame$event_index))
      expected_identities <- paste0(acquisition$acquisition_id,
                                    ":event_index:", indices)
      if (!all(identity_columns %in% names(all_frame)) ||
          !identical(as.integer(all_output$rows), nrow(all_frame)) ||
          nrow(all_frame) > maximum_all_events_rows ||
          !identical(as.integer(acquisition$audit$all_event_count), nrow(all_frame)) ||
          length(unique(all_frame$acquisition_id)) != 1L ||
          !identical(unique(all_frame$acquisition_id)[[1L]],
                     acquisition$acquisition_id) ||
          anyNA(indices) || any(indices < 0) || any(indices %% 1 != 0) ||
          anyDuplicated(indices) || anyDuplicated(all_frame$event_identity) ||
          !identical(all_frame$event_identity, expected_identities)) {
        stop("Model-gating all-events identity/count audit differs for ",
             acquisition$prefix, ".", call. = FALSE)
      }
      all_event_frames[[acquisition$prefix]] <- all_frame
    }
  }
  all_event_frames <- all_event_frames[configured]
  for (acquisition in acquisitions) {
    audit <- acquisition$audit
    identities <- population_frames[[acquisition$prefix]]
    all_identities <- all_event_frames[[acquisition$prefix]]$event_identity
    if ((length(all_identities) &&
         !all(identities$single_cells %in% all_identities)) ||
        !all(identities$g1 %in% identities$single_cells) ||
        !all(identities$edu_positive %in% identities$single_cells) ||
        length(intersect(identities$g1, identities$edu_positive))) {
      stop("Model-gating CSV child identities are not contained in Single Cells for ",
           acquisition$prefix, ".", call. = FALSE)
    }
    if (!isTRUE(audit$g1_subset_single_cells) ||
        !isTRUE(audit$edu_positive_subset_single_cells) ||
        !isTRUE(audit$g1_excludes_edu_positive) ||
        !identical(as.integer(audit$model_single_cells_count),
                   as.integer(acquisition$outputs$single_cells$rows)) ||
        !identical(as.integer(audit$model_g1_count),
                   as.integer(acquisition$outputs$g1$rows)) ||
        !identical(as.integer(audit$model_edu_positive_count),
                   as.integer(acquisition$outputs$edu_positive$rows))) {
      stop("Model-gating containment/count audit is incomplete or inconsistent.", call. = FALSE)
    }
    all_frame <- all_event_frames[[acquisition$prefix]]
    if (is.data.frame(all_frame) && nrow(all_frame)) {
      required_scores <- c(
        "DNA content", "model_single_cells_probability",
        "model_g1_probability", "model_edu_positive_probability"
      )
      if (!all(required_scores %in% names(all_frame))) {
        stop("Model-gating all-events table lacks G1 DNA minimum inputs.",
             call. = FALSE)
      }
      dna <- suppressWarnings(as.numeric(all_frame[["DNA content"]]))
      preliminary <- (
        all_frame$model_single_cells_probability >= expected_thresholds$single_cells &
        all_frame$model_g1_probability >= expected_thresholds$g1 &
        all_frame$model_edu_positive_probability < expected_thresholds$edu_positive
      )
      center_values <- dna[preliminary & is.finite(dna) & dna > 0]
      center <- stats::median(center_values)
      minimum <- expected_g1_dna_minimum$fraction * center
      expected_g1 <- all_frame$event_identity[
        preliminary & is.finite(dna) & dna >= minimum
      ]
      numeric_equal <- function(observed, expected) {
        is.numeric(observed) && length(observed) == 1L &&
          isTRUE(all.equal(observed, expected, tolerance = 1e-12))
      }
      if (!length(center_values) || !is.finite(center) ||
          !identical(identities$g1, expected_g1) ||
          !identical(as.integer(audit$preliminary_g1_count), sum(preliminary)) ||
          !identical(as.integer(audit$g1_dna_positive_candidate_count),
                     length(center_values)) ||
          !numeric_equal(audit$g1_dna_2n_center_raw, center) ||
          !numeric_equal(audit$g1_dna_minimum_fraction,
                         expected_g1_dna_minimum$fraction) ||
          !numeric_equal(audit$g1_dna_minimum_raw, minimum) ||
          !identical(as.integer(audit$g1_dna_minimum_excluded_count),
                     sum(preliminary) - length(expected_g1))) {
        stop("Model-gating adaptive G1 DNA minimum audit differs for ",
             acquisition$prefix, ".", call. = FALSE)
      }
    }
  }
  result <- do.call(rbind, verified)
  attr(result, "model_artifacts") <- do.call(rbind, artifact_rows)
  attr(result, "all_events") <- all_event_frames
  invisible(result)
}

#' @export
print.facs_analysis <- function(x, ...) {
  cat("<facs_analysis>\n")
  cat("  mode:       ", x$config$plot_type, "\n", sep = "")
  cat("  samples:    ", nrow(x$sample_manifest), "\n", sep = "")
  cat("  replicates: ", length(unique(x$sample_manifest$replicate_index)),
      "\n", sep = "")
  cat("  models:     ", length(x$models), "\n", sep = "")
  cat("  warnings:   ", length(x$warnings), "\n", sep = "")
  invisible(x)
}

validate_documented_model_parent_membership <- function(parent, g1, edu_positive) {
  parent_ids <- as.character(parent$event_identity)
  for (entry in list(g1 = g1, edu_positive = edu_positive)) {
    child_ids <- as.character(entry$event_identity)
    if (any(!child_ids %in% parent_ids) ||
        !identical(child_ids, parent_ids[parent_ids %in% child_ids])) {
      stop("Verified model child identities are not an order-preserving subset of predicted Single Cells.",
           call. = FALSE)
    }
  }
  invisible(TRUE)
}

format_documented_model_qc_warnings <- function(prefix, warnings) {
  if (!length(warnings)) return(character())
  paste(prefix, warnings, sep = ":")
}

validate_documented_model_manifest_header <- function(manifest, expected_mapping_hash) {
  expected_scopes <- list(
    single_cells = "complete_source_frame",
    g1 = "predicted_single_cells_parent",
    edu_positive = "complete_source_frame_then_predicted_single_cells_parent"
  )
  if (!identical(manifest$schema_version, 3L) ||
      !identical(manifest$status_label,
                 "EXPERIMENTAL MODEL-DERIVED NON-PRODUCTION") ||
      !identical(manifest$predicted_parent_authorization,
        "OWNER_AUTHORIZED_PREDICTED_PARENT_EXPERIMENTAL_BRANCH_ONLY") ||
      !setequal(names(manifest$feature_scopes), names(expected_scopes)) ||
      !all(vapply(names(expected_scopes), function(name) {
        identical(manifest$feature_scopes[[name]], expected_scopes[[name]])
      }, logical(1))) ||
      !identical(manifest$approved_mapping$sha256, expected_mapping_hash) ||
      !identical(manifest$single_cells_source, "model") ||
      !identical(manifest$g1_source, "model") ||
      !identical(manifest$edu_positive_source, "model") ||
      !identical(manifest$flowjo_workspace_gate_used, FALSE)) {
    stop("Model-only schema, authorization, sources, feature scopes, or approved mapping differs.",
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Validate provenance for experimental model-derived EdU gate inputs
#'
#' @param analysis A completed `facs_analysis` in EdU mode.
#' @param manifest_path Explicit path to `model-gating-manifest.json`.
#' @return An invisible data frame describing the verified input artifacts.
#' @export
validate_edu_documented_model_gate_manifest <- function(analysis, manifest_path) {
  read_once <- function(path, maximum_bytes = 512 * 1024^2) {
    info <- file.info(path)
    if (is.na(info$size) || info$isdir || info$size <= 0 ||
        info$size > maximum_bytes) {
      stop("Expected one bounded, nonempty regular provenance file.", call. = FALSE)
    }
    value <- facs_read_file_bytes(path, info$size)
    if (length(value) != info$size) {
      stop("Provenance file changed while it was being read.", call. = FALSE)
    }
    value
  }
  if (!inherits(analysis, "facs_analysis") ||
      !identical(analysis$config$plot_type, "edu")) {
    stop("`analysis` must be a completed EdU facs_analysis.", call. = FALSE)
  }
  if (!config_scalar_string(manifest_path) || !file.exists(manifest_path)) {
    stop("Model-gating manifest path must name one existing file.", call. = FALSE)
  }
  manifest_path <- normalizePath(manifest_path, mustWork = TRUE)
  manifest_bytes <- read_once(manifest_path, 16 * 1024^2)
  manifest <- jsonlite::fromJSON(rawToChar(manifest_bytes), simplifyVector = FALSE)
  expected_scopes <- list(
    single_cells = "complete_source_frame",
    g1 = "predicted_single_cells_parent",
    edu_positive = "complete_source_frame_then_predicted_single_cells_parent"
  )
  mapping_path <- system.file(
    "config", "model_only_edu_figure6_mapping.json",
    package = "facspseudocolor", mustWork = TRUE
  )
  expected_mapping_hash <-
    "7aa668a606a02ea77a0d9a0a69d90ffdac91a37bda5642351283173b2ad6ab7d"
  mapping_bytes <- read_once(mapping_path)
  observed_mapping_hash <- paste0(
    as.character(openssl::sha256(mapping_bytes)), collapse = ""
  )
  if (!identical(observed_mapping_hash, expected_mapping_hash))
    stop("Installed approved mapping digest differs.", call. = FALSE)
  validate_documented_model_manifest_header(manifest, expected_mapping_hash)
  approved_mapping <- jsonlite::fromJSON(rawToChar(mapping_bytes), simplifyVector = FALSE)
  acquisitions <- manifest$acquisitions
  if (!is.list(acquisitions) || !length(acquisitions)) {
    stop("Model-gating manifest has no acquisitions.", call. = FALSE)
  }
  artifacts <- manifest$artifacts
  if (!is.list(artifacts) ||
      !setequal(names(artifacts), c("single_cells", "g1", "edu_positive"))) {
    stop("Model-gating manifest must name exactly three population-model artifacts.",
         call. = FALSE)
  }
  artifact_rows <- list()
  for (target in names(artifacts)) {
    artifact <- artifacts[[target]]
    artifact_path <- if (identical(target, "g1")) artifact$portable_path else artifact$path
    if (!config_scalar_string(artifact_path) ||
        !config_scalar_string(artifact$byte_sha256) ||
        !file.exists(artifact_path)) {
      stop("Model artifact provenance is incomplete for ", target, ".",
           call. = FALSE)
    }
    artifact_bytes <- read_once(artifact_path)
    observed_artifact_hash <- paste0(
      as.character(openssl::sha256(artifact_bytes)), collapse = ""
    )
    if (!identical(target, "g1") &&
        !identical(observed_artifact_hash, artifact$byte_sha256)) {
      stop("Model artifact byte SHA-256 differs for ", target, ".", call. = FALSE)
    }
    expected <- list(
      single_cells = list(id = "SINGLE_CELLS_FROZEN_DEVELOPMENT_RULE",
        threshold = 0.205, byte = "71b459d8fb00930f31a6d289a21f587226fd2d6be4b31ebc782b7bae7839a376"),
      g1 = list(id = "EDU_DOCUMENTED_EDU_G1_RUN001_EDUPANEL001_GENERIC_FL2_FL4",
        threshold = 0.29,
        byte = "1af119ba5b64c845ec99430667b202e50cf7d2fe66b53648eae25a086b833182"),
      edu_positive = list(id = "EDU_POSITIVE_MODEL002", threshold = 0.385,
        byte = "6ef5603555a660b8503379cbbcab131b616336a64330ffd42f45dfaa182042cc")
    )[[target]]
    if (!identical(artifact$model_id, expected$id) ||
        !identical(artifact$byte_sha256, expected$byte) ||
        !isTRUE(all.equal(as.numeric(artifact$threshold), expected$threshold,
                          tolerance = 1e-15))) {
      stop("Model identity or threshold differs from the approved contract for ",
           target, ".", call. = FALSE)
    }
    expected_features <- if (identical(target, "g1")) c(
      "rank_fsc_a", "rank_fsc_h", "rank_ssc_a", "rank_ssc_h",
      "rank_fl2_a", "rank_fl2_h", "dna_signed_log_area_minus_pulse",
      "dna_rank_area_minus_pulse", "log_density_fsc_ssc",
      "log_density_dna_geometry", "rank_edu_area", "log_density_dna_edu"
    ) else if (identical(target, "edu_positive")) c(
      "dna_area_q", "dna_pulse_q", "dna_tls_position", "dna_tls_distance",
      "dna_area_pulse_log_density", "edu_area_q", "dna_edu_log_density"
    ) else NULL
    if (!is.null(expected_features) &&
        !identical(unlist(artifact$feature_schema, use.names = FALSE),
                   expected_features)) {
      stop("Frozen feature schema/order differs for ", target, ".", call. = FALSE)
    }
    if (identical(target, "single_cells")) {
      rule <- jsonlite::fromJSON(rawToChar(artifact_bytes), simplifyVector = FALSE)
      if (!config_scalar_string(artifact$semantic_sha256) ||
          !config_scalar_string(artifact$feature_schema) ||
          !identical(rule$sha256, artifact$semantic_sha256) ||
          !identical(rule$schema_version, artifact$feature_schema) ||
          !isTRUE(all.equal(as.numeric(rule$threshold), 0.205,
                            tolerance = 1e-15))) {
        stop("Single Cells semantic/schema/threshold provenance differs.", call. = FALSE)
      }
    } else {
      metadata <- artifact$metadata_sha256
      expected_metadata <- if (identical(target, "g1")) list() else list(
        `SUMMARY.json` = "dcd16784a80b68070eeb42126d6411db148646af6acc69ded7fb6a529b39d22f",
        `PREFLIGHT.json` = "9bee03f4b039d7645140df43ec03cd586ab7f6d620c0f359646caa65f49914db"
      )
      metadata_matches <- identical(target, "g1") || (is.list(metadata) &&
        setequal(names(metadata), names(expected_metadata)) &&
        all(vapply(names(expected_metadata), function(name) {
          identical(metadata[[name]], expected_metadata[[name]])
        }, logical(1))))
      if (!metadata_matches) {
        stop("Frozen model metadata provenance is missing for ", target, ".",
             call. = FALSE)
      }
      for (name in names(expected_metadata)) {
        metadata_path <- file.path(dirname(artifact_path), name)
        if (!file.exists(metadata_path)) {
          stop("Frozen metadata file is absent for ", target, ".", call. = FALSE)
        }
        observed <- paste0(as.character(openssl::sha256(
          read_once(metadata_path))), collapse = "")
        if (!identical(observed, expected_metadata[[name]])) {
          stop("Frozen metadata SHA-256 differs for ", target, ".", call. = FALSE)
        }
      }
      if (identical(target, "g1") &&
          (!identical(artifact$panel_id, "EDUPANEL001_generic_fl2_fl4") ||
           !identical(artifact$config_sha256,
             "f1a0e6e5327008a1da06ac912a2c25b95b060b579ea3bdf47540b7ca0dafbd88"))) {
        stop("Documented-EdU G1 panel/configuration provenance differs.", call. = FALSE)
      }
      if (identical(target, "edu_positive")) {
        if (!config_scalar_string(artifact$config_path) ||
            !file.exists(artifact$config_path) ||
            !identical(artifact$config_sha256,
              "624a839a38a464976b214afcf8964a5f54009ef7694d16d71a91ba20677d44f3")) {
          stop("Frozen EdU configuration provenance is incomplete.", call. = FALSE)
        }
        config_hash <- paste0(as.character(openssl::sha256(
          read_once(artifact$config_path))), collapse = "")
        if (!identical(config_hash, artifact$config_sha256)) {
          stop("Frozen EdU configuration SHA-256 differs.", call. = FALSE)
        }
      }
      expected_portable <- if (identical(target, "g1"))
        "ea1ccd90328674fb8cc1b6802a66e194f49b3f9f9c1e03c6577ad01fe880c96e" else
        "9e6ed466687efe5f7523dea633e9e9244381c0e613c86f0944f861c06047fdf4"
      if (!config_scalar_string(artifact$portable_path) ||
          !file.exists(artifact$portable_path) ||
          !identical(artifact$portable_sha256, expected_portable)) {
        stop("Portable frozen model provenance is incomplete for ", target,
             ".", call. = FALSE)
      }
      portable_hash <- paste0(as.character(openssl::sha256(
        read_once(artifact$portable_path))), collapse = "")
      if (!identical(portable_hash, expected_portable)) {
        stop("Portable frozen model SHA-256 differs for ", target, ".",
             call. = FALSE)
      }
      portable_root <- dirname(artifact$portable_path)
      package_records <- if (identical(target, "g1")) list(
        `MANIFEST.json` = "1964c6d976b4f25e9b6cfe9def0fb56d53b168fc5f8767ad8a29e7a4c430fe63",
        `COMPATIBILITY_REPORT.json` = "336a7ee3ec3516b6bdd2af1efd11e6300669a68bb335d7b2460096045f4db64f",
        `EXPORT_PROVENANCE.json` = "31d0e75b8d43da8a432a2e10580e7e2b620fcc9124288435822689bf019240a5",
        `NO_PROTECTED_DATA_PROOF.json` = "d81148613b5b31225537e16dd5fb74d70c4225df09203e111e8249fa382d8edb",
        `SYNTHETIC_EQUIVALENCE_FIXTURES.json` = "596125818978609b43c20fefe7a95acd739a01d7487b461074fde1b883a78793",
        `portable_hgb_predictor.py` = "2edf1bbd529159e8752d01731e87f19efbc63b514409b1d246d7a4ffc266a2c3",
        `run_portable_documented_edu_g1_selftest.py` = "8d4c0ff2a93add4302e33dcd9be2d4144bc3849dc12a70a785c333653104c37d"
      ) else list(
        `MANIFEST.json` = "c3885b0c49497672ad83657e2cd5a08beb219201039f5b26cffa002a4dd4fd7f",
        `COMPATIBILITY_REPORT.json` = "60695bb0a5d3d728437aa89aa570b9ee9d2c7be7cacfbbb564ab43b4309f3b38",
        `EXPORT_PROVENANCE.json` = "e0b80bfc44faf26886f04f9093eac1d8ed6474f89a43a90a5a8184819bbd40f8",
        `NO_PROTECTED_DATA_PROOF.json` = "f0b0b3aa0088c70e0f3e29cd75c5bc7e9fc8c809986b86f2ae73dbcb6301111b",
        `portable_hgb_predictor.py` = "2edf1bbd529159e8752d01731e87f19efbc63b514409b1d246d7a4ffc266a2c3"
      )
      for (name in names(package_records)) {
        package_path <- file.path(portable_root, name)
        if (!file.exists(package_path)) {
          stop("Portable conversion provenance file is absent: ", name, ".",
               call. = FALSE)
        }
        package_hash <- paste0(as.character(openssl::sha256(
          read_once(package_path))), collapse = "")
        if (!identical(package_hash, package_records[[name]])) {
          stop("Portable conversion provenance digest differs: ", name, ".",
               call. = FALSE)
        }
      }
    }
    artifact_rows[[length(artifact_rows) + 1L]] <- data.frame(
      target = target, byte_sha256 = artifact$byte_sha256,
      model_id = artifact$model_id, threshold = as.numeric(artifact$threshold),
      stringsAsFactors = FALSE
    )
  }
  prefixes <- vapply(acquisitions, `[[`, character(1), "prefix")
  approved_prefixes <- vapply(approved_mapping$acquisitions, `[[`,
                              character(1), "prefix")
  configured <- as.character(analysis$sample_manifest$prefix)
  if (length(prefixes) != 8L || !identical(prefixes, approved_prefixes) ||
      anyDuplicated(prefixes) || !setequal(prefixes, configured) ||
      length(prefixes) != length(configured)) {
    stop("Model-gating manifest prefixes do not exactly match the analysis.",
         call. = FALSE)
  }
  input_paths <- normalizePath(analysis$input_report$path, mustWork = TRUE)
  input_keys <- paste(analysis$input_report$prefix,
                      analysis$input_report$population, sep = "\r")
  population_keys <- c(single_cells = "complete", g1 = "g1",
                       edu_positive = "edu_positive")
  verified <- list()
  verified_frames <- list()
  qc_warnings <- character()
  for (acquisition in acquisitions) {
    mapped <- Filter(function(x) identical(x$prefix, acquisition$prefix),
                     approved_mapping$acquisitions)
    required_roles <- approved_mapping$channel_roles
    if (length(mapped) != 1L ||
        !identical(acquisition$acquisition_id, mapped[[1]]$acquisition_id) ||
        !identical(acquisition$source_fcs_basename, mapped[[1]]$fcs_basename) ||
        !identical(acquisition$source_fcs_sha256, mapped[[1]]$fcs_sha256) ||
        !identical(acquisition$consumed_fcs_snapshot_sha256,
                   mapped[[1]]$fcs_sha256) ||
        !identical(acquisition$functional_panel,
                   approved_mapping$functional_panel) ||
        !identical(
          unlist(acquisition$channel_roles[sort(names(acquisition$channel_roles))],
                 use.names = TRUE),
          unlist(required_roles[sort(names(required_roles))], use.names = TRUE)
        )) {
      stop("Acquisition identity, FCS digest, panel, or channel mapping differs from the approved mapping.",
           call. = FALSE)
    }
    audit <- acquisition$audit
    required_audit <- c("all_event_count", "model_single_cells_count",
      "model_g1_count", "model_edu_positive_count",
      "g1_subset_predicted_single_cells",
      "edu_positive_subset_predicted_single_cells",
      "parent_event_identity_exact_alignment", "all_event_identity_sha256",
      "predicted_single_cells_identity_sha256", "model_g1_identity_sha256",
      "model_edu_positive_identity_sha256", "g1_feature_scope",
      "edu_positive_feature_scope", "qc_warnings")
    if (!setequal(names(audit), required_audit) ||
        !identical(audit$g1_subset_predicted_single_cells, TRUE) ||
        !identical(audit$edu_positive_subset_predicted_single_cells, TRUE) ||
        !identical(audit$parent_event_identity_exact_alignment, TRUE) ||
        !identical(audit$g1_feature_scope, expected_scopes$g1) ||
        !identical(audit$edu_positive_feature_scope,
                   expected_scopes$edu_positive)) {
      stop("Model-only acquisition audit is incomplete or inconsistent.", call. = FALSE)
    }
    allowed_qc <- paste0(rep(c("single_cells", "g1", "edu_positive"), each = 3L),
      ":", rep(c("LOW_EVENT_SUPPORT", "EXTREME_PREDICTED_FRACTION",
                  "ZERO_PREDICTED_POSITIVES"), times = 3L))
    acquisition_qc <- as.character(unlist(audit$qc_warnings, use.names = FALSE))
    if (any(!acquisition_qc %in% allowed_qc) ||
        anyDuplicated(acquisition_qc)) {
      stop("Model QC warning record contains an unknown or duplicate code.", call. = FALSE)
    }
    qc_warnings <- c(qc_warnings, format_documented_model_qc_warnings(
      acquisition$prefix, acquisition_qc
    ))
    for (output_name in names(population_keys)) {
      key <- paste(acquisition$prefix, population_keys[[output_name]], sep = "\r")
      positions <- which(input_keys == key)
      if (length(positions) != 1L) {
        stop("Analysis input is missing or ambiguous for ", key, ".",
             call. = FALSE)
      }
      output <- acquisition$outputs[[output_name]]
      if (is.null(output$path) || is.null(output$sha256)) {
        stop("Model-gating manifest output record is incomplete for ", key, ".",
             call. = FALSE)
      }
      if (!config_scalar_string(output$path) ||
          !identical(basename(output$path), output$path) ||
          output$path %in% c(".", "..")) {
        stop("Model-gating manifest output path must be one path-free basename for ",
             key, ".", call. = FALSE)
      }
      expected_path <- normalizePath(
        file.path(dirname(manifest_path), output$path), mustWork = TRUE
      )
      if (!identical(dirname(expected_path), dirname(manifest_path))) {
        stop("Model-gating manifest output path escapes its manifest directory for ",
             key, ".", call. = FALSE)
      }
      if (!identical(input_paths[[positions]], expected_path)) {
        stop("Analysis input path differs from the model-gating manifest for ",
             key, ".", call. = FALSE)
      }
      input_bytes <- read_once(expected_path)
      observed_hash <- paste0(
        as.character(openssl::sha256(input_bytes)), collapse = ""
      )
      if (!identical(observed_hash, output$sha256)) {
        stop("Analysis input SHA-256 differs from the model-gating manifest for ",
             key, ".", call. = FALSE)
      }
      frame <- utils::read.csv(text = rawToChar(input_bytes),
                               stringsAsFactors = FALSE, check.names = FALSE)
      expected_rows <- switch(output_name,
        single_cells = audit$model_single_cells_count,
        g1 = audit$model_g1_count,
        edu_positive = audit$model_edu_positive_count)
      expected_identity_hash <- switch(output_name,
        single_cells = audit$predicted_single_cells_identity_sha256,
        g1 = audit$model_g1_identity_sha256,
        edu_positive = audit$model_edu_positive_identity_sha256)
      identities <- paste0(frame$event_identity, collapse = "\n")
      if (nrow(frame)) identities <- paste0(identities, "\n")
      observed_identity_hash <- paste0(as.character(openssl::sha256(
        charToRaw(identities))), collapse = "")
      event_index <- suppressWarnings(as.numeric(frame$event_index))
      expected_identity <- paste0(acquisition$acquisition_id,
                                  ":event_index:", event_index)
      if (!identical(nrow(frame), as.integer(output$rows)) ||
          !identical(nrow(frame), as.integer(expected_rows)) ||
          anyDuplicated(frame$event_identity) ||
          any(!is.finite(event_index)) || any(event_index < 0) ||
          any(event_index %% 1 != 0) || is.unsorted(event_index,
            strictly = TRUE) ||
          !identical(as.character(frame$event_identity), expected_identity) ||
          !identical(observed_identity_hash, output$event_identity_sha256) ||
          !identical(observed_identity_hash, expected_identity_hash)) {
        stop("Model output row count or event identity provenance differs for ",
             key, ".", call. = FALSE)
      }
      verified_frames[[key]] <- frame
      verified[[length(verified) + 1L]] <- data.frame(
        prefix = acquisition$prefix,
        population = population_keys[[output_name]], path = expected_path,
        sha256 = observed_hash, stringsAsFactors = FALSE
      )
    }
    parent <- verified_frames[[paste(acquisition$prefix, "complete", sep = "\r")]]
    validate_documented_model_parent_membership(
      parent,
      verified_frames[[paste(acquisition$prefix, "g1", sep = "\r")]],
      verified_frames[[paste(acquisition$prefix, "edu_positive", sep = "\r")]]
    )
  }
  result <- do.call(rbind, verified)
  attr(result, "model_artifacts") <- do.call(rbind, artifact_rows)
  attr(result, "verified_frames") <- verified_frames
  attr(result, "qc_warnings") <- unique(qc_warnings)
  invisible(result)
}

#' Validate provenance for an explicitly selected EdU model-gating profile
#'
#' @param analysis A completed or provisional EdU analysis.
#' @param manifest_path Explicit model-gating manifest path.
#' @return Verified model-gating input records.
#' @export
validate_edu_model_gate_manifest <- function(analysis, manifest_path) {
  profile <- analysis$config$gating$profile
  if (identical(profile, "single_g1_edu_frozen_v1")) {
    return(validate_edu_frozen_model_gate_manifest(analysis, manifest_path))
  }
  if (identical(profile, "documented_edu_g1_v1")) {
    return(validate_edu_documented_model_gate_manifest(analysis, manifest_path))
  }
  stop("The EdU model-gating profile is absent or unsupported.", call. = FALSE)
}
