# Configuration parsing and validation ---------------------------------------

`%||%` <- function(value, default) if (is.null(value)) default else value

facs_config_keys <- function() {
  c(
    "plot_type", "data_dir", "dna_channel", "target_channel", "target_name",
    "target_display_mode",
    "suffixes", "dna_2n_value", "normalize_target", "g1_anchor",
    "g1_source", "edu_positive_source", "g1_model_manifest",
    "baseline_fit_x_range", "baseline_boundary_bins",
    "baseline_minimum_events_per_bin", "baseline_minimum_negative_events",
    "show_edu_apex_line", "edu_apex_x_range", "edu_apex_density_adjust",
    "background_quantile", "show_cutoff_line", "poi_dna_align", "poi_peak_failure",
    "ph3_boundary_sensitivity_fraction",
    "show_reference_panel",
    "y_limit_lower_quantile", "y_limit_upper_quantile", "y_limits",
    "palette", "y_log10", "pseudocolor_signal",
    "background_subtracted_offset", "x_limits", "point_size", "density_bandwidth",
    "density_lower_clip", "density_upper_clip", "density_gamma",
    "show_phase_gates", "show_gate_assignment", "gate_assignment_max_points",
    "gate_show_labels", "gate_color", "gate_linetype", "gate_linewidth",
    "gate_label_size", "s_phase_bins", "g1_x_range", "g2m_x_range",
    "negative_y_range", "s_phase_y_range", "quantify_phase_median",
    "quantify_whole_median", "quantify_phase_percent", "quant_error_bar",
    "quant_signal", "quantify_reference_normalized",
    "quant_reference_condition", "quant_show_points", "quant_phase_lineplot",
    "bar_colors", "layout", "layout_options", "pdf_width",
    "pdf_height_per_row", "output_pdf", "output_png", "samples",
    "replicates", "flowjo", "gating", "report", "ph3_input_profile", "ph3_export_operation_dirs",
    "ph3_positivity_method",
    "ph3_output_contract", "ph3_pilot", "synchronized_dna_strategy",
    "synchronized_reference_prefix", "synchronized_minimum_g1_events"
  )
}

validate_edu_report_config <- function(report) {
  if (is.null(report)) return(NULL)
  if (!is.list(report)) stop("`report` must be a mapping.")
  allowed <- c(
    "all_events_dir", "dna_height_channel", "show_apex_comparison_toggle",
    "show_phase_gate_toggle", "embed_pseudocolor_pdf_downloads"
  )
  unknown <- setdiff(names(report), allowed)
  if (length(unknown)) {
    stop("Unknown EdU `report` setting(s): ", paste(unknown, collapse = ", "), ".")
  }
  for (key in c("all_events_dir", "dna_height_channel")) {
    if (!config_scalar_string(report[[key]])) {
      stop("EdU `report.", key, "` must be one non-empty string.")
    }
  }
  for (key in c("show_apex_comparison_toggle", "show_phase_gate_toggle",
                "embed_pseudocolor_pdf_downloads")) {
    if (!is.logical(report[[key]]) || length(report[[key]]) != 1L ||
        is.na(report[[key]])) {
      stop("EdU `report.", key, "` must be true or false.")
    }
  }
  report
}

facs_config_defaults <- function(plot_type) {
  list(
    suffixes = if (identical(plot_type, "edu")) {
      list(complete = "_single_cells.csv", g1 = "_g1.csv",
           edu_positive = "_edu_positive.csv")
    } else if (identical(plot_type, "synchronized")) {
      list(complete = "_single_cells.csv", g1 = "_g1.csv")
    } else if (identical(plot_type, "ph3")) {
      list(complete = "_single_cells.csv", g1 = "_g1.csv",
           ph3_positive = "_ph3_positive.csv")
    } else {
      list(complete = "_single_cells.csv")
    },
    dna_2n_value = 1000,
    normalize_target = !identical(plot_type, "synchronized"),
    g1_anchor = "median",
    synchronized_dna_strategy = "per_sample_g1",
    synchronized_reference_prefix = NULL,
    synchronized_minimum_g1_events = 100,
    baseline_fit_x_range = c(1000, 2000),
    baseline_boundary_bins = 20,
    baseline_minimum_events_per_bin = 20,
    baseline_minimum_negative_events = 100,
    show_edu_apex_line = FALSE,
    edu_apex_x_range = c(1400, 1600),
    edu_apex_density_adjust = 1,
    background_quantile = 0.95,
    show_cutoff_line = TRUE,
    poi_dna_align = "per_sample",
    poi_peak_failure = "error",
    ph3_boundary_sensitivity_fraction = 0.05,
    ph3_positivity_method = "flowjo_legacy_v1",
    y_limit_lower_quantile = 0.001,
    y_limit_upper_quantile = 0.999,
    palette = "refined",
    y_log10 = TRUE,
    # PH3-legacy mode (ph3_positivity_method: flowjo_legacy_v1) fits no
    # background model -- target_bgsub is always identical to target_raw --
    # so "background_subtracted" display has no real background to
    # subtract and only adds a display offset that the exact FlowJo gate
    # geometry (plotted from unshifted raw/y_raw coordinates; see
    # read_ph3_gate_geometry()) does not share, making the gate box appear
    # displaced from the data it is meant to outline. "normalized" (offset
    # always 0) is the correct default for ph3; edu/poi keep the prior
    # default unchanged.
    pseudocolor_signal = if (identical(plot_type, "ph3")) "normalized" else "background_subtracted",
    background_subtracted_offset = "auto",
    x_limits = c(700, 2250),
    point_size = 0.3,
    density_bandwidth = 0.5,
    density_lower_clip = 0.05,
    density_upper_clip = 0.90,
    density_gamma = 0.35,
    show_phase_gates = FALSE,
    show_gate_assignment = FALSE,
    gate_assignment_max_points = 4000,
    gate_show_labels = FALSE,
    gate_color = "black",
    gate_linetype = "dashed",
    gate_linewidth = 0.5,
    gate_label_size = 2.2,
    quantify_phase_median = FALSE,
    quantify_whole_median = FALSE,
    quantify_phase_percent = FALSE,
    quant_error_bar = "sd",
    quant_signal = "background_subtracted",
    quantify_reference_normalized = FALSE,
    quant_show_points = TRUE,
    quant_phase_lineplot = TRUE,
    layout = "cowplot",
    layout_options = list(),
    pdf_width = 11.5,
    pdf_height_per_row = 2.5
  )
}

validate_edu_gating_config <- function(gating, replicates, data_dir, config_path) {
  if (is.null(gating)) gating <- list(mode = "flowjo")
  if (!is.list(gating) || !config_scalar_string(gating$mode) ||
      !gating$mode %in% c("flowjo", "model_experimental")) {
    stop("EdU `gating.mode` must be `flowjo` or `model_experimental`.")
  }
  if (identical(gating$mode, "flowjo")) {
    if (!identical(names(gating), "mode"))
      stop("FlowJo gating accepts only `gating.mode`; model settings cannot be mixed in.")
    return(gating)
  }
  if (identical(gating$profile, "documented_edu_g1_v1")) {
    required <- c("mode", "profile", "status_label", "manifest")
    if (!setequal(names(gating), required)) {
      stop("The documented EdU G1 profile requires exactly: ",
           paste(required, collapse = ", "), ".")
    }
    if (!identical(gating$status_label,
                   "EXPERIMENTAL MODEL-DERIVED NON-PRODUCTION") ||
        !config_scalar_string(gating$manifest)) {
      stop("The documented EdU G1 profile requires its exact non-production status and manifest path.")
    }
    return(gating)
  }
  required <- c("mode", "profile", "status_label", "output_dir", "artifacts", "acquisitions")
  if (!setequal(names(gating), required))
    stop("Experimental model gating requires exactly: ", paste(required, collapse = ", "), ".")
  if (!identical(gating$profile, "single_g1_edu_frozen_v1") ||
      !identical(gating$status_label, "EXPERIMENTAL MODEL-DERIVED NON-PRODUCTION"))
    stop("Experimental model gating must select the approved frozen profile and exact status label.")
  if (!config_scalar_string(gating$output_dir) || !config_scalar_string(data_dir))
    stop("Experimental model gating requires explicit output and data directories.")
  base <- if (config_scalar_string(config_path)) dirname(normalizePath(config_path, mustWork = FALSE)) else getwd()
  resolve <- function(path) normalizePath(if (grepl("^/", path)) path else file.path(base, path), mustWork = FALSE)
  if (!identical(resolve(gating$output_dir), resolve(data_dir)))
    stop("In model mode, `data_dir` must equal `gating.output_dir`; no fallback input is allowed.")
  if (!is.list(gating$artifacts) ||
      !setequal(names(gating$artifacts), c("single_cells", "g1", "edu_positive")))
    stop("The frozen model profile requires exactly Single Cells, G1, and EdU-positive artifacts.")
  frozen_hashes <- c(
    single_cells = "71b459d8fb00930f31a6d289a21f587226fd2d6be4b31ebc782b7bae7839a376",
    g1 = "d653d388d15af3bd22185cb9c8e1addea132e0b5d1d18e60e37c65a7548163b7",
    edu_positive = "9e6ed466687efe5f7523dea633e9e9244381c0e613c86f0944f861c06047fdf4"
  )
  frozen_semantic_hashes <- c(
    single_cells = "f9a53e9fd78d2f39cf5980b8f470c3c44e7a187fa942cbcc0324c563fa68a31a",
    g1 = "32723ccc768ed517affa040f2ce62a125dad465399d869f471bc331549f47f6a"
  )
  for (target in names(frozen_hashes)) {
    artifact <- gating$artifacts[[target]]
    required_artifact <- if (identical(target, "edu_positive"))
      c("path", "byte_sha256") else
      c("path", "byte_sha256", "semantic_sha256")
    if (!is.list(artifact) || !identical(names(artifact), required_artifact) ||
        !all(vapply(artifact, config_scalar_string, logical(1))) ||
        !identical(artifact$byte_sha256, unname(frozen_hashes[[target]])) ||
        (!identical(target, "edu_positive") &&
         !identical(artifact$semantic_sha256,
                    unname(frozen_semantic_hashes[[target]]))))
      stop("Model artifact configuration differs from frozen profile for ", target, ".")
  }
  manifest <- make_sample_manifest(replicates = replicates)
  acquisitions <- gating$acquisitions
  if (!is.list(acquisitions) || length(acquisitions) != nrow(manifest))
    stop("Model gating requires exactly one acquisition mapping per configured sample.")
  prefixes <- vapply(acquisitions, function(x) x$prefix %||% "", character(1))
  if (anyDuplicated(prefixes) || !setequal(prefixes, manifest$prefix))
    stop("Model-gating acquisition prefixes must exactly match configured samples.")
  required_roles <- c("fsc_area", "ssc_area", "dna_area", "dna_height", "edu")
  for (item in acquisitions) {
    if (!setequal(names(item), c("prefix", "acquisition_id", "fcs_path", "fcs_sha256", "channel_roles")) ||
        !all(vapply(item[c("prefix", "acquisition_id", "fcs_path", "fcs_sha256")], config_scalar_string, logical(1))) ||
        !grepl("^[0-9a-f]{64}$", item$fcs_sha256) || !is.list(item$channel_roles) ||
        !identical(names(item$channel_roles), required_roles) ||
        !all(vapply(item$channel_roles, config_scalar_string, logical(1))) ||
        anyDuplicated(unlist(item$channel_roles)))
      stop("Each model-gating acquisition requires exact identity, FCS SHA-256, and five distinct ordered channel roles.")
  }
  gating
}

config_scalar_string <- function(value) {
  is.character(value) && length(value) == 1L && !is.na(value) && nzchar(value)
}

config_scalar_number <- function(value) {
  is.numeric(value) && length(value) == 1L && is.finite(value)
}

config_numeric_range <- function(value, increasing = TRUE) {
  value <- as.numeric(unlist(value))
  length(value) == 2L && all(is.finite(value)) &&
    (!increasing || value[[1]] < value[[2]])
}

config_add_error <- function(errors, message) c(errors, message)

#' Validate a flow-cytometry analysis configuration
#'
#' Validates configuration structure and mode-specific settings without reading
#' experimental data. Defaults are applied centrally and returned in the
#' validated configuration. Unknown keys are rejected to catch misspellings.
#'
#' @param config A configuration list, usually produced by
#'   [read_facs_config()].
#' @param config_path Optional path to the source YAML file. Relative data and
#'   output paths are interpreted relative to this file's directory.
#'
#' @return A validated object of class `facs_config`.
#' @export
validate_facs_config <- function(config, config_path = attr(config, "config_path")) {
  if (!is.list(config) || is.data.frame(config)) {
    stop("`config` must be a named list.", call. = FALSE)
  }
  if (is.null(names(config)) || any(!nzchar(names(config)))) {
    stop("Every top-level configuration setting must have a name.", call. = FALSE)
  }

  errors <- character()
  legacy_edu_fields <- c("g1_source", "edu_positive_source", "g1_model_manifest")
  supplied_legacy_edu <- legacy_edu_fields[!vapply(
    config[legacy_edu_fields], is.null, logical(1)
  )]
  unknown <- setdiff(names(config), facs_config_keys())
  if (length(unknown)) {
    errors <- config_add_error(
      errors, paste0("Unknown configuration setting(s): ",
                     paste(unknown, collapse = ", "), ".")
    )
  }

  required_strings <- c("plot_type", "data_dir", "dna_channel",
                        "target_channel", "target_name", "output_pdf",
                        "output_png")
  for (key in required_strings) {
    if (!config_scalar_string(config[[key]])) {
      errors <- config_add_error(errors,
        paste0("`", key, "` must be a nonempty string."))
    }
  }

  plot_type <- config$plot_type
  if (config_scalar_string(plot_type) &&
      !plot_type %in% c("edu", "poi", "ph3", "synchronized")) {
    errors <- config_add_error(
      errors, "`plot_type` must be 'edu', 'poi', 'ph3', or 'synchronized'."
    )
  }
  if (!config_scalar_string(plot_type) ||
      !plot_type %in% c("edu", "poi", "ph3", "synchronized")) {
    plot_type <- "edu"
  }

  if (identical(plot_type, "edu") && length(supplied_legacy_edu)) {
    if (!is.null(config$gating)) {
      errors <- config_add_error(
        errors,
        "Do not mix legacy EdU source fields with the explicit `gating` block."
      )
    } else if (!config_scalar_string(config$g1_source) ||
               !config_scalar_string(config$edu_positive_source) ||
               !identical(config$g1_source, config$edu_positive_source) ||
               !config$g1_source %in% c("model", "flowjo")) {
      errors <- config_add_error(
        errors,
        "Legacy EdU source fields must explicitly select the same `model` or `flowjo` source."
      )
    } else if (identical(config$g1_source, "model")) {
      if (!config_scalar_string(config$g1_model_manifest)) {
        errors <- config_add_error(
          errors,
          "Legacy model-derived EdU configuration requires `g1_model_manifest`."
        )
      } else {
        config$gating <- list(
          mode = "model_experimental",
          profile = "documented_edu_g1_v1",
          status_label = "EXPERIMENTAL MODEL-DERIVED NON-PRODUCTION",
          manifest = config$g1_model_manifest
        )
      }
    } else if (!is.null(config$g1_model_manifest)) {
      errors <- config_add_error(
        errors, "Legacy FlowJo EdU configuration must not declare `g1_model_manifest`."
      )
    } else {
      config$gating <- list(mode = "flowjo")
    }
    config$g1_source <- NULL
    config$edu_positive_source <- NULL
    config$g1_model_manifest <- NULL
  }

  config <- utils::modifyList(facs_config_defaults(plot_type), config)

  if (identical(plot_type, "edu")) {
    gating_error <- tryCatch({
      config$gating <- validate_edu_gating_config(
        config$gating, config$replicates, config$data_dir, config_path
      ); NULL
    }, error = function(e) conditionMessage(e))
    if (!is.null(gating_error)) errors <- config_add_error(errors, gating_error)
    report_error <- tryCatch({
      config$report <- validate_edu_report_config(config$report)
      NULL
    }, error = function(e) conditionMessage(e))
    if (!is.null(report_error)) errors <- config_add_error(errors, report_error)
    if (!is.null(config$target_display_mode) &&
        !config$target_display_mode %in% c("edu", "poi_cdc45")) {
      errors <- config_add_error(
        errors, "`target_display_mode` must be 'edu' or 'poi_cdc45'."
      )
    }
  } else if (!is.null(config$report)) {
    errors <- config_add_error(errors, "`report` settings are supported only for EdU configurations.")
  } else if (!is.null(config$target_display_mode)) {
    errors <- config_add_error(
      errors, "`target_display_mode` is supported only for EdU configurations."
    )
  }

  if (identical(plot_type, "synchronized")) {
    if (!identical(config$normalize_target, FALSE)) {
      errors <- config_add_error(
        errors, "Synchronized mode requires `normalize_target: false`; target values are not transformed."
      )
    }
    if (!config_scalar_string(config$synchronized_dna_strategy) ||
        !config$synchronized_dna_strategy %in%
          c("per_sample_g1", "shared_asynchronous_g1")) {
      errors <- config_add_error(
        errors, "`synchronized_dna_strategy` must be 'per_sample_g1' or 'shared_asynchronous_g1'."
      )
    }
    if (!config_scalar_number(config$synchronized_minimum_g1_events) ||
        config$synchronized_minimum_g1_events < 2 ||
        config$synchronized_minimum_g1_events %% 1 != 0) {
      errors <- config_add_error(
        errors, "`synchronized_minimum_g1_events` must be an integer of at least 2."
      )
    }
    if (identical(config$synchronized_dna_strategy, "shared_asynchronous_g1")) {
      if (!config_scalar_string(config$synchronized_reference_prefix)) {
        errors <- config_add_error(
          errors, "Shared-asynchronous normalization requires an explicit `synchronized_reference_prefix`."
        )
      } else if (!is.null(config$replicates)) {
        prefixes <- unlist(lapply(config$replicates, function(rep) {
          vapply(rep$samples, `[[`, character(1), "prefix")
        }), use.names = FALSE)
        if (sum(prefixes == config$synchronized_reference_prefix) != 1L) {
          errors <- config_add_error(
            errors, "`synchronized_reference_prefix` must identify exactly one configured sample."
          )
        }
      }
    }
  }

  if (!xor(is.null(config$samples), is.null(config$replicates))) {
    errors <- config_add_error(
      errors, "Specify exactly one of `samples` or `replicates`."
    )
  }
  if (!is.null(config$samples)) {
    errors <- config_add_error(
      errors,
      "A flat `samples` list cannot express acquisition structure; use a `replicates` list. POI mode also requires an explicit background-control `reference`."
    )
  }

  if (!is.null(config$replicates)) {
    manifest_error <- tryCatch({
      manifest <- make_sample_manifest(replicates = config$replicates)
      duplicate_prefixes <- unique(manifest$prefix[duplicated(manifest$prefix)])
      if (length(duplicate_prefixes)) {
        stop("Duplicate sample prefix(es): ", paste(duplicate_prefixes, collapse = ", "))
      }
      reference_counts <- vapply(
        split(manifest, manifest$model_group),
        function(x) sum(x$is_reference), integer(1)
      )
      reference_declared <- vapply(
        split(manifest, manifest$model_group),
        function(x) any(!is.na(x$reference_condition)), logical(1)
      )
      if (plot_type == "edu" && any(reference_declared) &&
          any(reference_counts != 1L)) {
        stop(
          "EdU replicate `reference` samples must be absent in every biological/technical replicate pair, or each declared reference must name exactly one matching sample in every pair."
        )
      }
      if (plot_type == "ph3" && any(reference_counts != 0L)) {
        stop("PH3 mode does not use replicate `reference` samples.")
      }
      if (plot_type == "poi" && any(reference_counts != 1L)) {
        stop("Each biological/technical replicate pair must contain exactly one matching `reference` condition.")
      }
      NULL
    }, error = function(e) conditionMessage(e))
    if (!is.null(manifest_error)) errors <- config_add_error(errors, manifest_error)
  }

  positive_numbers <- c(
    "dna_2n_value", "baseline_boundary_bins",
    "baseline_minimum_events_per_bin", "baseline_minimum_negative_events",
    "edu_apex_density_adjust", "point_size", "density_bandwidth",
    "density_gamma", "pdf_width", "pdf_height_per_row"
  )
  for (key in positive_numbers) {
    if (!config_scalar_number(config[[key]]) || config[[key]] <= 0) {
      errors <- config_add_error(errors,
        paste0("`", key, "` must be a positive finite number."))
    }
  }

  for (key in c("x_limits", "baseline_fit_x_range", "edu_apex_x_range")) {
    if (!config_numeric_range(config[[key]])) {
      errors <- config_add_error(errors,
        paste0("`", key, "` must contain two increasing finite numbers."))
    }
  }
  if (!is.null(config$y_limits) && !config_numeric_range(config$y_limits)) {
    errors <- config_add_error(
      errors, "`y_limits` must contain two increasing finite numbers."
    )
  }

  if (!config$g1_anchor %in% c("median", "mode")) {
    errors <- config_add_error(errors, "`g1_anchor` must be 'median' or 'mode'.")
  }
  if (!config$poi_dna_align %in% c("per_sample", "shared_background")) {
    errors <- config_add_error(
      errors, "`poi_dna_align` must be 'per_sample' or 'shared_background'."
    )
  }
  if (!config$poi_peak_failure %in% c("error", "use_background")) {
    errors <- config_add_error(
      errors, "`poi_peak_failure` must be 'error' or 'use_background'."
    )
  }
  if (!config$layout %in% c("plotgardener", "cowplot")) {
    errors <- config_add_error(errors, "`layout` must be 'plotgardener' or 'cowplot'.")
  }
  if (!config$quant_error_bar %in% c("sd", "sem", "none")) {
    errors <- config_add_error(errors, "`quant_error_bar` must be 'sd', 'sem', or 'none'.")
  }
  if (!config$quant_signal %in% c("background_subtracted", "normalized")) {
    errors <- config_add_error(
      errors, "`quant_signal` must be 'background_subtracted' or 'normalized'."
    )
  }
  if (!config$pseudocolor_signal %in% c("background_subtracted", "normalized")) {
    errors <- config_add_error(
      errors,
      "`pseudocolor_signal` must be 'background_subtracted' or 'normalized'."
    )
  }
  offset_is_auto <- config_scalar_string(config$background_subtracted_offset) &&
    identical(config$background_subtracted_offset, "auto")
  offset_is_number <- config_scalar_number(config$background_subtracted_offset) &&
    config$background_subtracted_offset >= 0
  if (!offset_is_auto && !offset_is_number) {
    errors <- config_add_error(
      errors,
      "`background_subtracted_offset` must be 'auto' or a nonnegative finite number."
    )
  }

  probabilities <- c("background_quantile", "density_lower_clip",
                     "density_upper_clip", "y_limit_lower_quantile",
                     "y_limit_upper_quantile",
                     "ph3_boundary_sensitivity_fraction")
  for (key in probabilities) {
    value <- config[[key]]
    if (!config_scalar_number(value) || value < 0 || value > 1) {
      errors <- config_add_error(errors,
        paste0("`", key, "` must be a finite number between 0 and 1."))
    }
  }
  if (config_scalar_number(config$density_lower_clip) &&
      config_scalar_number(config$density_upper_clip) &&
      config$density_lower_clip >= config$density_upper_clip) {
    errors <- config_add_error(errors,
      "`density_lower_clip` must be less than `density_upper_clip`.")
  }
  if (config_scalar_number(config$y_limit_lower_quantile) &&
      config_scalar_number(config$y_limit_upper_quantile) &&
      config$y_limit_lower_quantile >= config$y_limit_upper_quantile) {
    errors <- config_add_error(errors,
      "`y_limit_lower_quantile` must be less than `y_limit_upper_quantile`.")
  }

  if (plot_type == "ph3") {
    if (!config_scalar_string(config$ph3_positivity_method) ||
        !config$ph3_positivity_method %in% c(
          "flowjo_legacy_v1", "ph3_raw_4n_density_cutoff_v1"
        )) {
      errors <- config_add_error(
        errors,
        "`ph3_positivity_method` must be `flowjo_legacy_v1` or `ph3_raw_4n_density_cutoff_v1`."
      )
    }
    if (!config_scalar_string(config$ph3_input_profile) ||
        !config$ph3_input_profile %in% c(
        "production_direct_identity_v1", "legacy_count_only_unverified_v1",
        "legacy_csv_pilot_v1"
        )) {
      errors <- config_add_error(
        errors,
        "PH3 mode requires explicit `ph3_input_profile`: 'production_direct_identity_v1', 'legacy_count_only_unverified_v1', or 'legacy_csv_pilot_v1'."
      )
    } else if (identical(config$ph3_input_profile,
                         "production_direct_identity_v1")) {
      raw_directories <- config$ph3_export_operation_dirs
      valid_directories <-
        (is.character(raw_directories) && length(raw_directories) > 0L &&
         all(vapply(as.list(raw_directories), config_scalar_string, logical(1)))) ||
        (is.list(raw_directories) && length(raw_directories) > 0L &&
         is.null(names(raw_directories)) &&
         all(vapply(raw_directories, config_scalar_string, logical(1))))
      if (!valid_directories) {
        errors <- config_add_error(
          errors,
          "Production PH3 mode requires one or more explicit nonempty `ph3_export_operation_dirs`."
        )
      }
      contract_error <- tryCatch({
        validated_contract <- validate_ph3_output_contract_config(
          config$ph3_output_contract, config$replicates
        )
        attr(validated_contract, "condition_table") <- NULL
        attr(validated_contract, "comparison_table") <- NULL
        attr(validated_contract, "replicate_set_ids") <- NULL
        config$ph3_output_contract <- validated_contract
        NULL
      }, error = function(e) conditionMessage(e))
      if (!is.null(contract_error)) {
        errors <- config_add_error(errors, contract_error)
      }
    } else if (identical(config$ph3_input_profile, "legacy_csv_pilot_v1")) {
      if (!identical(config$ph3_positivity_method,
                     "ph3_raw_4n_density_cutoff_v1")) {
        errors <- config_add_error(errors,
          "Legacy CSV pilot requires `ph3_positivity_method: ph3_raw_4n_density_cutoff_v1`.")
      }
      pilot <- config$ph3_pilot
      if (!is.list(pilot) || !identical(names(pilot), c("control_label", "provenance_label")) ||
          !identical(pilot$control_label, "Untreated") ||
          !identical(pilot$provenance_label, "PILOT / LIMITED-PROVENANCE \u2014 NOT PUBLICATION-GRADE")) {
        errors <- config_add_error(errors,
          "`ph3_pilot` must exactly declare control_label: 'Untreated' and provenance_label: 'PILOT / LIMITED-PROVENANCE \u2014 NOT PUBLICATION-GRADE'.")
      }
      if (!is.list(config$flowjo) || !identical(config$flowjo$rebuild, FALSE)) {
        errors <- config_add_error(errors,
          "Legacy CSV pilot requires `flowjo.rebuild: false`; it never regenerates historical exports.")
      }
    }
    ph3_ranges <- list(
      g1 = config$g1_x_range,
      early = if (is.list(config$s_phase_bins)) config$s_phase_bins$early else NULL,
      mid = if (is.list(config$s_phase_bins)) config$s_phase_bins$mid else NULL,
      late = if (is.list(config$s_phase_bins)) config$s_phase_bins$late else NULL,
      g2m = config$g2m_x_range
    )
    invalid_ranges <- names(ph3_ranges)[
      !vapply(ph3_ranges, config_numeric_range, logical(1))
    ]
    if (length(invalid_ranges)) {
      errors <- config_add_error(
        errors,
        paste0(
          "PH3 mode requires explicit increasing DNA ranges for: ",
          paste(invalid_ranges, collapse = ", "), "."
        )
      )
    } else {
      boundaries <- c(
        as.numeric(unlist(ph3_ranges$g1)),
        as.numeric(unlist(ph3_ranges$early)),
        as.numeric(unlist(ph3_ranges$mid)),
        as.numeric(unlist(ph3_ranges$late)),
        as.numeric(unlist(ph3_ranges$g2m))
      )
      joins <- c(boundaries[2] - boundaries[3],
                 boundaries[4] - boundaries[5],
                 boundaries[6] - boundaries[7],
                 boundaries[8] - boundaries[9])
      if (any(abs(joins) > sqrt(.Machine$double.eps))) {
        errors <- config_add_error(
          errors,
          "PH3 DNA gates must be contiguous: G1, Early S, Mid S, Late S, and G2/M boundaries must meet exactly."
        )
      }
      sensitivity_delta <- config$ph3_boundary_sensitivity_fraction *
        config$dna_2n_value
      g2m_width <- diff(as.numeric(unlist(ph3_ranges$g2m)))
      if (!is.finite(sensitivity_delta) || sensitivity_delta <= 0 ||
          sensitivity_delta >= g2m_width) {
        errors <- config_add_error(
          errors,
          "`ph3_boundary_sensitivity_fraction * dna_2n_value` must be positive and smaller than the configured G2/M gate width."
        )
      }
    }
  }

  required_suffixes <- if (plot_type == "edu") {
    c("complete", "g1", "edu_positive")
  } else if (plot_type == "synchronized") {
    c("complete", "g1")
  } else if (plot_type == "ph3") {
    if (identical(config$ph3_input_profile, "legacy_csv_pilot_v1")) {
      c("complete", "g1")
    } else {
      c("complete", "g1", "ph3_positive")
    }
  } else {
    "complete"
  }
  if (!is.list(config$suffixes) ||
      !all(required_suffixes %in% names(config$suffixes)) ||
      any(!vapply(config$suffixes[required_suffixes], config_scalar_string,
                  logical(1)))) {
    errors <- config_add_error(
      errors,
      paste0("`suffixes` must define nonempty ",
             paste(required_suffixes, collapse = ", "), " suffixes for ",
             plot_type, " mode.")
    )
  }

  if (length(errors)) {
    stop(
      paste0("Invalid FACS configuration:\n- ", paste(unique(errors), collapse = "\n- ")),
      call. = FALSE
    )
  }

  if (!is.null(config_path)) {
    config_path <- normalizePath(config_path, mustWork = FALSE)
    attr(config, "config_path") <- config_path
    attr(config, "config_dir") <- dirname(config_path)
  }
  class(config) <- c("facs_config", "list")
  config
}

#' Read and validate a FACS analysis configuration
#'
#' Reads a YAML configuration from an explicit path. The package does not
#' search the working directory or substitute an example configuration.
#'
#' @param path Path to a YAML configuration file.
#' @param validate If `TRUE`, apply defaults and validate the complete
#'   configuration with [validate_facs_config()].
#'
#' @return A configuration list. When `validate = TRUE`, the object also has
#'   class `facs_config` and records the normalized source path.
#' @export
read_facs_config <- function(path, validate = TRUE) {
  if (!config_scalar_string(path)) {
    stop("`path` must be one explicit, nonempty file path.", call. = FALSE)
  }
  if (!file.exists(path)) {
    stop("Configuration file not found: ", path, call. = FALSE)
  }
  path <- normalizePath(path, mustWork = TRUE)
  config <- yaml::read_yaml(path)
  if (!is.list(config) || is.null(config)) {
    stop("Configuration YAML must contain a top-level mapping.", call. = FALSE)
  }
  attr(config, "config_path") <- path
  attr(config, "config_dir") <- dirname(path)
  if (isTRUE(validate)) validate_facs_config(config, path) else config
}
