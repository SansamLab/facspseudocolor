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

new_facs_analysis <- function(
    config, manifest, input_report, normalized_data, models,
    quantitation = list(), warnings = character()
) {
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
      provenance = list(
        package = "facspseudocolor",
        package_version = as.character(
          utils::packageVersion("facspseudocolor")
        ),
        r_version = R.version.string,
        config_path = attr(config, "config_path"),
        config_dir = attr(config, "config_dir"),
        input_files = input_report$path[input_report$exists]
      )
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
    normalized_data <- lapply(seq_len(nrow(manifest)), function(i) {
      prefix <- manifest$prefix[[i]]
      label <- manifest$condition[[i]]
      complete_file <- file.path(
        directory, paste0(prefix, config$suffixes$complete)
      )
      g1_file <- file.path(directory, paste0(prefix, config$suffixes$g1))
      positive_file <- file.path(
        directory, paste0(prefix, config$suffixes$ph3_positive)
      )
      normalize_ph3(
        events = complete_file,
        g1_events = g1_file,
        ph3_positive_events = positive_file,
        dna_channel = config$dna_channel,
        target_channel = config$target_channel,
        dna_2n_value = config$dna_2n_value,
        g1_anchor = config$g1_anchor,
        sample_id = label
      )
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
  if (config$plot_type == "ph3") {
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
