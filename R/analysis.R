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

  model_gate_inputs <- NULL
  model_gate_manifest_path <- NULL
  model_verified_frames <- NULL
  if (identical(config$plot_type, "edu") &&
      identical(config$g1_source, "model")) {
    model_gate_manifest_path <- path.expand(config$g1_model_manifest)
    if (!grepl("^(/|[A-Za-z]:[/\\\\])", model_gate_manifest_path)) {
      config_dir <- attr(config, "config_dir")
      if (!config_scalar_string(config_dir)) {
        stop("Relative `g1_model_manifest` requires a file-backed configuration.",
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
    analysis$provenance$g1_calling <- if (identical(config$g1_source, "model")) {
      list(source = "model", manifest = normalizePath(
        model_gate_manifest_path, mustWork = TRUE
      ), verified_inputs = model_gate_inputs,
      verified_artifacts = attr(model_gate_inputs, "model_artifacts"),
      qc_warnings = attr(model_gate_inputs, "qc_warnings"),
      populations = c("Single Cells", "G1", "EdU Positive"),
      flowjo_workspace_gate_used = FALSE,
      interpretation = "experimental model-derived research output; not expert gating or biological ground truth")
    } else {
      list(source = "flowjo", manifest = NULL)
    }
    if (identical(config$g1_source, "model")) {
      analysis$warnings <- unique(c(
        analysis$warnings,
        paste0("MODEL_QC_REVIEW_REQUIRED: ", attr(model_gate_inputs, "qc_warnings"))
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

validate_model_parent_membership <- function(parent, g1, edu_positive) {
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

format_model_qc_warnings <- function(prefix, warnings) {
  if (!length(warnings)) return(character())
  paste(prefix, warnings, sep = ":")
}

validate_model_manifest_header <- function(manifest, expected_mapping_hash) {
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
validate_edu_model_gate_manifest <- function(analysis, manifest_path) {
  read_once <- function(path) {
    info <- file.info(path)
    if (is.na(info$size) || info$isdir) stop("Expected one regular provenance file.", call. = FALSE)
    connection <- file(path, open = "rb")
    value <- readBin(connection, what = "raw", n = info$size)
    close(connection)
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
  manifest_bytes <- read_once(manifest_path)
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
    "a874815c83816d6d4291054fceeb08fa0df4b0c864499176dafe3a6a400597fb"
  mapping_bytes <- read_once(mapping_path)
  observed_mapping_hash <- paste0(
    as.character(openssl::sha256(mapping_bytes)), collapse = ""
  )
  if (!identical(observed_mapping_hash, expected_mapping_hash))
    stop("Installed approved mapping digest differs.", call. = FALSE)
  validate_model_manifest_header(manifest, expected_mapping_hash)
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
    if (!config_scalar_string(artifact$path) ||
        !config_scalar_string(artifact$byte_sha256) ||
        !file.exists(artifact$path)) {
      stop("Model artifact provenance is incomplete for ", target, ".",
           call. = FALSE)
    }
    artifact_bytes <- read_once(artifact$path)
    observed_artifact_hash <- paste0(
      as.character(openssl::sha256(artifact_bytes)), collapse = ""
    )
    if (!identical(observed_artifact_hash, artifact$byte_sha256)) {
      stop("Model artifact byte SHA-256 differs for ", target, ".", call. = FALSE)
    }
    expected <- list(
      single_cells = list(id = "SINGLE_CELLS_FROZEN_DEVELOPMENT_RULE",
        threshold = 0.205, byte = "71b459d8fb00930f31a6d289a21f587226fd2d6be4b31ebc782b7bae7839a376"),
      g1 = list(id = "DNA_PROJECTED_G1_CLUSTER001", threshold = 0.42,
        byte = "9ce7e448455e2d42090a86c80beb10d925a7c0e64e6261e678335385d9714af6"),
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
      "dna_area_q", "dna_pulse_q", "tls_position", "tls_signed_distance",
      "tls_abs_distance", "dna_2d_log_density", "pulse_area_log_ratio"
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
      expected_metadata <- if (identical(target, "g1")) list(
        `SUMMARY.json` = "7a7446d61ca175e56964fc5ab6edbee081f2eb3f16e7925542157bc3719ddb10",
        `PREFLIGHT.json` = "bac957a8837774afee0faba2ff2e832d003f1e37586a68e6af4037083448ec80",
        `per_unit_metrics.csv` = "650076ecf722f165549731a0c84f1cbdde13d7cd36b53653da4a5b7b34b03792",
        `model_comparison.csv` = "54f538ddb00b28a932181d5cd57dda0b3bfe37786a7f4d3357e0a71e9c5c67d2",
        `derived_label_summary.csv` = "e60a3c24bd2ef5b80913d1cbb5387b74df0eef2ac4d67b1dfc5945162943380f",
        `SCOPE.csv` = "264b0ab4ef806933fc03eb66f644ed4405fe77d196bee516a7b32a98acc6c608"
      ) else list(
        `SUMMARY.json` = "dcd16784a80b68070eeb42126d6411db148646af6acc69ded7fb6a529b39d22f",
        `PREFLIGHT.json` = "9bee03f4b039d7645140df43ec03cd586ab7f6d620c0f359646caa65f49914db"
      )
      metadata_matches <- is.list(metadata) &&
        setequal(names(metadata), names(expected_metadata)) &&
        all(vapply(names(expected_metadata), function(name) {
          identical(metadata[[name]], expected_metadata[[name]])
        }, logical(1)))
      if (!metadata_matches) {
        stop("Frozen model metadata provenance is missing for ", target, ".",
             call. = FALSE)
      }
      for (name in names(expected_metadata)) {
        metadata_path <- file.path(dirname(artifact$path), name)
        if (!file.exists(metadata_path)) {
          stop("Frozen metadata file is absent for ", target, ".", call. = FALSE)
        }
        observed <- paste0(as.character(openssl::sha256(
          read_once(metadata_path))), collapse = "")
        if (!identical(observed, expected_metadata[[name]])) {
          stop("Frozen metadata SHA-256 differs for ", target, ".", call. = FALSE)
        }
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
        "7fe256e892adfef6415e3957231bf8f68b7ad401dc909ffd2a801fbde8f4b75d" else
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
      package_records <- list(
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
      target = target, byte_sha256 = observed_artifact_hash,
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
    qc_warnings <- c(qc_warnings, format_model_qc_warnings(
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
    validate_model_parent_membership(
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
