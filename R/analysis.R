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

  if (config$plot_type == "edu") {
    models <- fit_sample_reference_models(
      manifest, directory, config$suffixes, settings
    )
    normalized_data <- lapply(seq_len(nrow(manifest)), function(i) {
      model <- models[[manifest$prefix[[i]]]]
      sample <- read_and_normalize_sample(
        prefix = manifest$prefix[[i]],
        condition_label = manifest$condition[[i]],
        data_dir = directory,
        file_suffixes = config$suffixes,
        settings = settings,
        baseline_slope = model$slope
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
validate_edu_model_gate_manifest <- function(analysis, manifest_path) {
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
  if (!identical(dirname(manifest_path), configured_output)) {
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
