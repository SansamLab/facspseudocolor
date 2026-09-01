# pH3 raw 4N density-cutoff positivity --------------------------------------

ph3_raw_4n_cutoff_method <- function() {
  list(
    method_id = "ph3_raw_4n_density_cutoff_v1",
    kernel = "gaussian", bandwidth = "nrd0", adjust = 1,
    grid_size = 2048L, minimum_control_4n_events = 100L,
    event_predicate = "identity_valid & eligible_2to4n & four_n_member",
    positivity_predicate = "target_raw > frozen_control_cutoff"
  )
}

ph3_raw_4n_cutoff_fail <- function(reason, detail) {
  stop("PH3 raw-4N cutoff failed [", reason, "]: ", detail, ".",
       call. = FALSE)
}

ph3_raw_4n_density_cutoff <- function(raw_signal) {
  method <- ph3_raw_4n_cutoff_method()
  if (!is.numeric(raw_signal) || length(raw_signal) < method$minimum_control_4n_events) {
    ph3_raw_4n_cutoff_fail("insufficient_control_4n_events",
                            paste0("at least ", method$minimum_control_4n_events,
                                   " finite untreated 4N events are required"))
  }
  if (any(!is.finite(raw_signal))) {
    ph3_raw_4n_cutoff_fail("nonfinite_control_raw_signal",
                            "every qualifying untreated 4N raw pH3 value must be finite")
  }
  signal_range <- range(raw_signal)
  if (!all(is.finite(signal_range)) || diff(signal_range) <= 0) {
    ph3_raw_4n_cutoff_fail("zero_or_invalid_control_signal_range",
                            "untreated 4N raw pH3 values must span a positive finite range")
  }
  density_error <- NULL
  fit <- tryCatch(stats::density(raw_signal, kernel = method$kernel,
                                 bw = method$bandwidth, adjust = method$adjust,
                                 n = method$grid_size, from = signal_range[[1L]],
                                 to = signal_range[[2L]]), error = function(e) {
                                   density_error <<- conditionMessage(e)
                                   NULL
                                 })
  if (is.null(fit) || !is.numeric(fit$x) || !is.numeric(fit$y) ||
      length(fit$x) != method$grid_size || length(fit$y) != method$grid_size ||
      any(!is.finite(fit$x)) || any(!is.finite(fit$y)) || any(diff(fit$x) <= 0)) {
    ph3_raw_4n_cutoff_fail("invalid_density", paste0(
      "the fixed raw-signal density grid is invalid",
      if (is.null(density_error)) "" else paste0(" (", density_error, ")")
    ))
  }
  maxima <- which(fit$y[2:(length(fit$y) - 1L)] > fit$y[1:(length(fit$y) - 2L)] &
                    fit$y[2:(length(fit$y) - 1L)] >= fit$y[3:length(fit$y)]) + 1L
  if (!length(maxima)) {
    ph3_raw_4n_cutoff_fail("no_interior_local_maximum",
                            "the untreated 4N density has no qualifying dominant peak")
  }
  greatest <- max(fit$y[maxima])
  peak <- maxima[fit$y[maxima] == greatest]
  if (length(peak) != 1L) {
    ph3_raw_4n_cutoff_fail("ambiguous_dominant_peak",
                            "the untreated 4N density has more than one equal dominant local maximum")
  }
  minima <- which(fit$y[2:(length(fit$y) - 1L)] < fit$y[1:(length(fit$y) - 2L)] &
                    fit$y[2:(length(fit$y) - 1L)] <= fit$y[3:length(fit$y)]) + 1L
  right_minima <- minima[minima > peak[[1L]]]
  if (!length(right_minima)) {
    ph3_raw_4n_cutoff_fail("no_right_side_local_minimum",
                            "the dominant untreated 4N peak has no qualifying right-side local minimum")
  }
  cutoff_index <- right_minima[[1L]]
  if (is.na(cutoff_index) ||
      !is.finite(fit$x[[cutoff_index]]) || fit$x[[cutoff_index]] < signal_range[[1L]] ||
      fit$x[[cutoff_index]] > signal_range[[2L]]) {
    ph3_raw_4n_cutoff_fail("no_right_side_local_minimum",
                            "the dominant untreated 4N peak has no qualifying right-side local minimum")
  }
  list(method = method, density_x = fit$x, density_y = fit$y,
       control_event_count = as.integer(length(raw_signal)), raw_min = signal_range[[1L]],
       raw_max = signal_range[[2L]], peak_index = as.integer(peak[[1L]]),
       peak_raw_signal = fit$x[[peak[[1L]]]], peak_density = fit$y[[peak[[1L]]]],
       cutoff_index = as.integer(cutoff_index), cutoff = fit$x[[cutoff_index]],
       cutoff_density = fit$y[[cutoff_index]])
}

ph3_apply_raw_4n_density_cutoff <- function(analysis) {
  if (!inherits(analysis, "facs_analysis") ||
      !identical(analysis$config$plot_type, "ph3") ||
      !identical(analysis$config$ph3_input_profile, "production_direct_identity_v1")) {
    ph3_raw_4n_cutoff_fail("invalid_analysis_profile",
                            "computed pH3 cutoff requires production direct-identity pH3 analysis")
  }
  contract <- validate_ph3_output_contract_config(
    analysis$config$ph3_output_contract, analysis$config$replicates
  )
  conditions <- attr(contract, "condition_table")
  replicate_ids <- attr(contract, "replicate_set_ids")
  manifest <- analysis$sample_manifest
  classifications <- lapply(analysis$normalized_data, `[[`, "ph3_event_classification")
  if (length(classifications) != nrow(manifest) ||
      any(!vapply(classifications, is.data.frame, logical(1))) ||
      anyDuplicated(manifest$prefix)) {
    ph3_raw_4n_cutoff_fail("missing_classification",
                            "exactly one retained classification is required for every configured acquisition")
  }
  records <- vector("list", length(replicate_ids))
  application <- vector("list", nrow(manifest))
  density_curves <- vector("list", length(replicate_ids))
  for (i in seq_along(replicate_ids)) {
    set_id <- replicate_ids[[i]]
    rows <- which(manifest$replicate_index == i)
    condition_id <- conditions$condition_id[match(manifest$condition[rows], conditions$condition_label)]
    control_rows <- rows[conditions$role[match(manifest$condition[rows], conditions$condition_label)] == "control"]
    if (length(control_rows) != 1L || length(rows) < 2L || anyNA(condition_id)) {
      ph3_raw_4n_cutoff_fail("ambiguous_control_mapping",
                              paste0("replicate set `", set_id,
                                     "` must contain exactly one explicit control and mapped conditions"))
    }
    control <- classifications[[control_rows[[1L]]]]
    required <- c("identity_valid", "eligible_2to4n", "four_n_member",
                  "target_raw", "event_identity", "acquisition_id", "sample_id")
    if (!all(required %in% names(control)) || nrow(control) == 0L ||
        any(!control$identity_valid)) {
      ph3_raw_4n_cutoff_fail("invalid_control_classification",
                              paste0("control classification for `", set_id,
                                     "` lacks valid retained direct identities"))
    }
    selected <- control$identity_valid & control$eligible_2to4n & control$four_n_member
    density <- tryCatch(ph3_raw_4n_density_cutoff(control$target_raw[selected]),
                        error = function(e) e)
    if (inherits(density, "error")) {
      reason <- sub("^PH3 raw-4N cutoff failed \\[([^]]+)\\].*$", "\\1",
                    conditionMessage(density))
      expected_failures <- c(
        "insufficient_control_4n_events", "nonfinite_control_raw_signal",
        "zero_or_invalid_control_signal_range", "invalid_density",
        "no_interior_local_maximum", "ambiguous_dominant_peak",
        "no_right_side_local_minimum"
      )
      if (!reason %in% expected_failures) stop(density)
      records[[i]] <- data.frame(
        experiment_id = contract$experiment_id, replicate_set_id = set_id,
        control_condition_id = conditions$condition_id[match(manifest$condition[[control_rows[[1L]]]], conditions$condition_label)],
        control_sample_id = control$sample_id[[1L]], control_acquisition_id = control$acquisition_id[[1L]],
        analysis_id = control$analysis_id[[1L]], config_digest = control$config_digest[[1L]],
        control_export_operation_id = control$export_operation_id[[1L]], control_input_manifest_key = control$input_manifest_key[[1L]],
        raw_channel_id = analysis$config$target_channel,
        event_selection_predicate = ph3_raw_4n_cutoff_method()$event_predicate,
        target_sample_ids = ph3_canonical_string_array(vapply(rows, function(row) classifications[[row]]$sample_id[[1L]], character(1))),
        target_acquisition_ids = ph3_canonical_string_array(vapply(rows, function(row) classifications[[row]]$acquisition_id[[1L]], character(1))),
        cutoff_method_id = ph3_raw_4n_cutoff_method()$method_id, kernel = ph3_raw_4n_cutoff_method()$kernel,
        bandwidth = ph3_raw_4n_cutoff_method()$bandwidth, adjust = ph3_raw_4n_cutoff_method()$adjust,
        grid_size = ph3_raw_4n_cutoff_method()$grid_size, qualifying_control_4n_event_count = as.integer(sum(selected)),
        control_raw_min = NA_real_, control_raw_max = NA_real_, dominant_peak_index = NA_integer_,
        dominant_peak_raw_signal = NA_real_, dominant_peak_density = NA_real_, cutoff_index = NA_integer_,
        cutoff_raw_signal = NA_real_, cutoff_density = NA_real_,
        cutoff_status = "unavailable_cutoff_failure", cutoff_reason_code = reason,
        cutoff_reason_detail = conditionMessage(density), stringsAsFactors = FALSE
      )
      density_curves[[i]] <- data.frame(replicate_set_id = set_id, raw_pH3 = numeric(), density = numeric(), stringsAsFactors = FALSE)
      for (row in rows) {
        classification <- classifications[[row]]
        classification$flowjo_ph3_positive_member <- classification$ph3_positive_member
        classification$flowjo_positivity_method_id <- classification$positivity_method_id
        classification$ph3_positive_member <- rep(FALSE, nrow(classification))
        classification$positivity_call_status <- rep("unavailable_cutoff_failure", nrow(classification))
        classification$positivity_call_reason_code <- rep(reason, nrow(classification))
        classification$positivity_method_id <- ph3_raw_4n_cutoff_method()$method_id
        classification$raw_4n_cutoff <- rep(NA_real_, nrow(classification))
        classifications[[row]] <- classification
        application[[row]] <- data.frame(experiment_id = contract$experiment_id, replicate_set_id = set_id,
          sample_id = classification$sample_id[[1L]], acquisition_id = classification$acquisition_id[[1L]],
          cutoff_method_id = ph3_raw_4n_cutoff_method()$method_id, cutoff_raw_signal = NA_real_,
          retained_event_count = as.integer(nrow(classification)), finite_raw_event_count = as.integer(sum(is.finite(classification$target_raw))),
          nonfinite_raw_event_count = as.integer(sum(!is.finite(classification$target_raw))),
          called_positive_count = 0L, called_negative_count = 0L,
          cutoff_status = "unavailable_cutoff_failure", cutoff_reason_code = reason,
          stringsAsFactors = FALSE)
      }
      next
    }
    target_samples <- vapply(rows, function(row) {
      classifications[[row]]$sample_id[[1L]]
    }, character(1))
    target_acquisitions <- vapply(rows, function(row) {
      classifications[[row]]$acquisition_id[[1L]]
    }, character(1))
    records[[i]] <- data.frame(
      experiment_id = contract$experiment_id, replicate_set_id = set_id,
      control_condition_id = conditions$condition_id[match(manifest$condition[[control_rows[[1L]]]], conditions$condition_label)],
      control_sample_id = control$sample_id[[1L]], control_acquisition_id = control$acquisition_id[[1L]],
      analysis_id = control$analysis_id[[1L]], config_digest = control$config_digest[[1L]],
      control_export_operation_id = control$export_operation_id[[1L]],
      control_input_manifest_key = control$input_manifest_key[[1L]],
      raw_channel_id = analysis$config$target_channel,
      event_selection_predicate = density$method$event_predicate,
      target_sample_ids = ph3_canonical_string_array(target_samples),
      target_acquisition_ids = ph3_canonical_string_array(target_acquisitions),
      cutoff_method_id = density$method$method_id, kernel = density$method$kernel,
      bandwidth = density$method$bandwidth, adjust = density$method$adjust,
      grid_size = density$method$grid_size, qualifying_control_4n_event_count = density$control_event_count,
      control_raw_min = density$raw_min, control_raw_max = density$raw_max,
      dominant_peak_index = density$peak_index, dominant_peak_raw_signal = density$peak_raw_signal,
      dominant_peak_density = density$peak_density, cutoff_index = density$cutoff_index,
      cutoff_raw_signal = density$cutoff, cutoff_density = density$cutoff_density,
      cutoff_status = "available", cutoff_reason_code = NA_character_, cutoff_reason_detail = NA_character_, stringsAsFactors = FALSE
    )
    density_curves[[i]] <- data.frame(
      replicate_set_id = set_id, raw_pH3 = density$density_x,
      density = density$density_y, stringsAsFactors = FALSE
    )
    for (row in rows) {
      classification <- classifications[[row]]
      required_target <- c("identity_valid", "target_raw", "event_identity", "acquisition_id", "sample_id",
                           "eligible_2to4n", "sub_4n_member", "four_n_member")
      if (!all(required_target %in% names(classification)) || any(!classification$identity_valid) ||
          anyDuplicated(classification$event_identity) ||
          anyNA(classification$eligible_2to4n) || anyNA(classification$sub_4n_member) ||
          anyNA(classification$four_n_member) ||
          any(classification$sub_4n_member & classification$four_n_member) ||
          any(classification$eligible_2to4n !=
                (classification$sub_4n_member | classification$four_n_member))) {
        ph3_raw_4n_cutoff_fail("invalid_target_classification",
                                paste0("target classification for `", set_id,
                                       "` lacks valid retained direct identities"))
      }
      finite <- is.finite(classification$target_raw)
      classification$flowjo_ph3_positive_member <- classification$ph3_positive_member
      classification$flowjo_positivity_method_id <- classification$positivity_method_id
      classification$ph3_positive_member <- finite & classification$target_raw > density$cutoff
      classification$positivity_call_status <- ifelse(finite, "called", "unavailable_nonfinite_raw_signal")
      classification$positivity_call_reason_code <- ifelse(finite, NA_character_, "nonfinite_raw_signal")
      classification$positivity_method_id <- density$method$method_id
      classification$raw_4n_cutoff <- rep(density$cutoff, nrow(classification))
      classifications[[row]] <- classification
      application[[row]] <- data.frame(
        experiment_id = contract$experiment_id, replicate_set_id = set_id,
        sample_id = classification$sample_id[[1L]], acquisition_id = classification$acquisition_id[[1L]],
        cutoff_method_id = density$method$method_id, cutoff_raw_signal = density$cutoff,
        retained_event_count = as.integer(nrow(classification)), finite_raw_event_count = as.integer(sum(finite)),
        nonfinite_raw_event_count = as.integer(sum(!finite)), called_positive_count = as.integer(sum(classification$ph3_positive_member)),
        called_negative_count = as.integer(sum(finite & !classification$ph3_positive_member)),
        cutoff_status = "available", cutoff_reason_code = NA_character_, stringsAsFactors = FALSE
      )
    }
  }
  for (i in seq_along(classifications)) analysis$normalized_data[[i]]$ph3_event_classification <- classifications[[i]]
  analysis$provenance$ph3_raw_4n_density_cutoff <- list(
    schema_version = "ph3-raw-4n-density-cutoff-1.0.0",
    method = ph3_raw_4n_cutoff_method(), records = do.call(rbind, records),
    application_qc = do.call(rbind, application),
    density_curves = do.call(rbind, density_curves)
  )
  analysis
}
