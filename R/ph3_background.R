# pH3 event-specific background regression --------------------------------

ph3_background_fail <- function(reason, detail) {
  stop(
    "PH3 background-regression validation failed [", reason, "]: ",
    detail, ".", call. = FALSE
  )
}

ph3_invalid_background_fit <- function(reason, negative_count,
                                       positive_count, negative_range,
                                       positive_range, coverage_status,
                                       design_rank = NA_integer_) {
  list(
    fit_status = "invalid", validity_reason_code = reason,
    negative_event_count = as.integer(negative_count),
    positive_event_count = as.integer(positive_count),
    design_rank = as.integer(design_rank), intercept = NA_real_,
    slope = NA_real_, negative_dna_min = negative_range[[1L]],
    negative_dna_max = negative_range[[2L]],
    positive_dna_min = positive_range[[1L]],
    positive_dna_max = positive_range[[2L]],
    coverage_status = coverage_status, predictions = NULL
  )
}

ph3_validate_background_solution <- function(coefficients, prediction_dna) {
  if (!is.numeric(coefficients) || length(coefficients) != 2L ||
      any(!is.finite(coefficients))) {
    return(list(
      status = "invalid", reason_code = "nonfinite_coefficients",
      predictions = NULL
    ))
  }
  if (!is.numeric(prediction_dna) || any(!is.finite(prediction_dna))) {
    ph3_background_fail(
      "invalid_prediction_vector",
      "prediction DNA supplied after fit validation must be finite numeric data"
    )
  }
  predictions <- coefficients[[1L]] + coefficients[[2L]] * prediction_dna
  if (any(!is.finite(predictions))) {
    return(list(
      status = "invalid", reason_code = "nonfinite_predictions",
      predictions = predictions
    ))
  }
  list(status = "valid", reason_code = "valid", predictions = predictions)
}

ph3_fit_background_model <- function(dna, signal, negative_member,
                                     positive_member) {
  if (!is.numeric(dna) || !is.numeric(signal) ||
      !is.logical(negative_member) || !is.logical(positive_member) ||
      length(dna) != length(signal) ||
      length(dna) != length(negative_member) ||
      length(dna) != length(positive_member) ||
      anyNA(negative_member) || anyNA(positive_member) ||
      any(negative_member & positive_member)) {
    ph3_background_fail(
      "invalid_fit_vectors",
      "DNA, signal, and disjoint negative/positive membership must align exactly"
    )
  }
  negative_count <- sum(negative_member)
  positive_count <- sum(positive_member)
  empty_range <- c(NA_real_, NA_real_)
  negative_range <- empty_range
  positive_range <- empty_range
  if (negative_count > 0L && all(is.finite(dna[negative_member]))) {
    negative_range <- range(dna[negative_member])
  }
  if (positive_count > 0L && all(is.finite(dna[positive_member]))) {
    positive_range <- range(dna[positive_member])
  }
  if (negative_count < 100L) {
    return(ph3_invalid_background_fit(
      "insufficient_negative_events", negative_count, positive_count,
      negative_range, positive_range, "not_evaluated"
    ))
  }
  if (any(!is.finite(dna[negative_member])) ||
      any(!is.finite(signal[negative_member]))) {
    return(ph3_invalid_background_fit(
      "nonfinite_fit_input", negative_count, positive_count,
      negative_range, positive_range, "not_evaluated"
    ))
  }
  if (positive_count > 0L && any(!is.finite(dna[positive_member]))) {
    return(ph3_invalid_background_fit(
      "nonfinite_prediction_input", negative_count, positive_count,
      negative_range, positive_range, "not_evaluated"
    ))
  }
  dna_variation <- diff(negative_range)
  if (!is.finite(dna_variation)) {
    return(ph3_invalid_background_fit(
      "nonfinite_dna_variation", negative_count, positive_count,
      negative_range, positive_range, "not_evaluated"
    ))
  }
  if (dna_variation <= 0) {
    return(ph3_invalid_background_fit(
      "zero_dna_variation", negative_count, positive_count,
      negative_range, positive_range, "not_evaluated"
    ))
  }
  design <- cbind(intercept = 1, dna = dna[negative_member])
  design_rank <- qr(design)$rank
  if (!identical(as.integer(design_rank), 2L)) {
    return(ph3_invalid_background_fit(
      "rank_deficient_design", negative_count, positive_count,
      negative_range, positive_range, "not_evaluated", design_rank
    ))
  }
  fit <- tryCatch(
    stats::lm.fit(x = design, y = signal[negative_member]),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(ph3_invalid_background_fit(
      "fit_error", negative_count, positive_count, negative_range,
      positive_range, "not_evaluated", design_rank
    ))
  }
  coefficients <- unname(fit$coefficients)
  coefficient_check <- ph3_validate_background_solution(
    coefficients, numeric()
  )
  if (!identical(coefficient_check$status, "valid")) {
    return(ph3_invalid_background_fit(
      coefficient_check$reason_code, negative_count, positive_count,
      negative_range, positive_range, "not_evaluated", design_rank
    ))
  }
  coverage_status <- if (positive_count == 0L) {
    "not_applicable_no_positive_events"
  } else if (positive_range[[1L]] >= negative_range[[1L]] &&
             positive_range[[2L]] <= negative_range[[2L]]) {
    "spans_positive_range"
  } else {
    "does_not_span_positive_range"
  }
  if (identical(coverage_status, "does_not_span_positive_range")) {
    return(ph3_invalid_background_fit(
      "positive_dna_extrapolation_required", negative_count, positive_count,
      negative_range, positive_range, coverage_status, design_rank
    ))
  }
  required_member <- negative_member | positive_member
  predictions <- rep(NA_real_, length(dna))
  prediction_check <- ph3_validate_background_solution(
    coefficients, dna[required_member]
  )
  predictions[required_member] <- prediction_check$predictions
  if (!identical(prediction_check$status, "valid")) {
    return(ph3_invalid_background_fit(
      prediction_check$reason_code, negative_count, positive_count,
      negative_range, positive_range, coverage_status, design_rank
    ))
  }
  list(
    fit_status = "valid", validity_reason_code = "valid",
    negative_event_count = as.integer(negative_count),
    positive_event_count = as.integer(positive_count),
    design_rank = as.integer(design_rank), intercept = coefficients[[1L]],
    slope = coefficients[[2L]], negative_dna_min = negative_range[[1L]],
    negative_dna_max = negative_range[[2L]],
    positive_dna_min = positive_range[[1L]],
    positive_dna_max = positive_range[[2L]],
    coverage_status = coverage_status, predictions = predictions
  )
}

ph3_background_event_table <- function(model) {
  classifications <- model$source$event_classifications
  mapping <- unique(
    model$source$quantitation$ph3_metrics_acquisition[c(
      "acquisition_id", "sample_id", "prefix"
    )]
  )
  required <- c(
    "analysis_id", "acquisition_id", "sample_id", "prefix",
    "event_identity", "identity_valid", "dna_finite",
    "ph3_positive_member", "eligible_2to4n",
    "sub_4n_member", "four_n_member", "dna_raw", "dna_norm", "target_raw"
  )
  if (!is.list(classifications) || !length(classifications) ||
      any(!vapply(classifications, is.data.frame, logical(1)))) {
    ph3_background_fail(
      "missing_event_classification",
      "the validated model must retain classification tables"
    )
  }
  classification_acquisition_ids <- vapply(classifications, function(data) {
    if (!is.data.frame(data) || !"acquisition_id" %in% names(data) ||
        !nrow(data) || length(unique(data$acquisition_id)) != 1L) {
      return(NA_character_)
    }
    as.character(data$acquisition_id[[1L]])
  }, character(1))
  if (length(classifications) != nrow(mapping) ||
      anyNA(classification_acquisition_ids) ||
      anyDuplicated(classification_acquisition_ids) ||
      !setequal(classification_acquisition_ids, mapping$acquisition_id)) {
    ph3_background_fail(
      "classification_acquisition_mismatch",
      "exactly one classification table is required for every authoritative acquisition"
    )
  }
  rows <- lapply(seq_along(classifications), function(i) {
    data <- classifications[[i]]
    if (!all(required %in% names(data)) || !nrow(data)) {
      ph3_background_fail(
        "missing_event_field",
        paste0("classification ", i, " lacks a required Slice 2 event field")
      )
    }
    scalar_identity <- c("analysis_id", "acquisition_id", "sample_id", "prefix")
    if (any(vapply(data[scalar_identity], function(x) {
      length(unique(x)) != 1L || is.na(x[[1L]]) || !nzchar(x[[1L]])
    }, logical(1))) || anyNA(data$event_identity) ||
        anyDuplicated(data$event_identity) ||
        !is.numeric(data$dna_raw) || !is.numeric(data$dna_norm) ||
        !is.numeric(data$target_raw) ||
        any(vapply(data[c(
          "ph3_positive_member", "eligible_2to4n", "sub_4n_member",
          "four_n_member", "identity_valid", "dna_finite"
        )], function(x) !is.logical(x) || anyNA(x), logical(1)))) {
      ph3_background_fail(
        "invalid_event_classification",
        paste0("classification ", i, " has ambiguous identity or event types")
      )
    }
    event_row <- as.integer(seq_len(nrow(data)))
    data.frame(
      analysis_id = data$analysis_id,
      acquisition_id = data$acquisition_id,
      sample_id = data$sample_id, prefix = data$prefix,
      event_identity = data$event_identity,
      event_row = event_row,
      identity_valid = data$identity_valid, dna_finite = data$dna_finite,
      ph3_positive_member = data$ph3_positive_member,
      eligible_2to4n = data$eligible_2to4n,
      sub_4n_member = data$sub_4n_member,
      four_n_member = data$four_n_member,
      dna_raw = data$dna_raw, dna_norm = data$dna_norm,
      target_raw = data$target_raw,
      stringsAsFactors = FALSE
    )
  })
  events <- do.call(rbind, rows)
  rownames(events) <- NULL
  if (any(events$eligible_2to4n &
          (!events$identity_valid | !events$dna_finite)) ||
      any(events$sub_4n_member & events$four_n_member) ||
      !identical(events$eligible_2to4n,
                 events$sub_4n_member | events$four_n_member)) {
    ph3_background_fail(
      "altered_event_membership",
      "eligible membership must retain validated identity and the exact DNA partition"
    )
  }
  event_mapping <- unique(events[c("acquisition_id", "sample_id", "prefix")])
  expected_mapping <- mapping[
    match(event_mapping$acquisition_id, mapping$acquisition_id), , drop = FALSE
  ]
  rownames(expected_mapping) <- NULL
  rownames(event_mapping) <- NULL
  if (anyNA(expected_mapping$acquisition_id) ||
      !identical(event_mapping, expected_mapping) ||
      any(events$analysis_id != model$experiment$analysis_id[[1L]]) ||
      anyDuplicated(events[c("acquisition_id", "event_identity")])) {
    ph3_background_fail(
      "altered_event_mapping",
      "event, acquisition, sample, prefix, and analysis identities must remain exact"
    )
  }
  ph3_reconcile_background_counts(
    events, model$source$quantitation$ph3_metrics_acquisition
  )
  sample_map <- model$samples[c("sample_id", "replicate_set_id")]
  if (anyDuplicated(sample_map$sample_id) ||
      !setequal(unique(events$sample_id), sample_map$sample_id)) {
    ph3_background_fail(
      "sample_mapping_mismatch",
      "classification sample IDs must map exactly to validated model samples"
    )
  }
  events$replicate_set_id <- sample_map$replicate_set_id[
    match(events$sample_id, sample_map$sample_id)
  ]
  events$experiment_id <- model$experiment$experiment_id[[1L]]
  events
}

ph3_reconcile_background_counts <- function(events, metrics) {
  metric_ids <- c(
    "ph3_2to4n_positivity_percent",
    "ph3_within_4n_positivity_percent",
    "ph3_4n_positive_prevalence_within_2to4n_percent",
    "ph3_within_sub_4n_positivity_percent",
    "ph3_sub_4n_positive_prevalence_within_2to4n_percent"
  )
  required_columns <- c(
    "acquisition_id", "metric_id", "numerator_count", "denominator_count"
  )
  if (!all(required_columns %in% names(metrics)) ||
      !identical(typeof(metrics$numerator_count), "integer") ||
      !identical(typeof(metrics$denominator_count), "integer")) {
    ph3_background_fail(
      "invalid_authoritative_counts",
      "authoritative acquisition metrics require exact integer count fields"
    )
  }
  for (acquisition_id in unique(events$acquisition_id)) {
    selected_rows <- metrics[
      metrics$acquisition_id == acquisition_id &
        metrics$metric_id %in% metric_ids, , drop = FALSE
    ]
    if (nrow(selected_rows) != length(metric_ids) ||
        anyNA(selected_rows$metric_id) ||
        anyDuplicated(selected_rows$metric_id) ||
        !setequal(selected_rows$metric_id, metric_ids)) {
      ph3_background_fail(
        "invalid_authoritative_counts",
        paste0("acquisition `", acquisition_id,
               "` must have exactly the five authoritative pH3 metrics")
      )
    }
    metric_rows <- selected_rows[
      match(metric_ids, selected_rows$metric_id), , drop = FALSE
    ]
    member <- events$acquisition_id == acquisition_id
    eligible <- events$eligible_2to4n[member]
    positive <- events$ph3_positive_member[member]
    sub_four <- events$sub_4n_member[member]
    four_n <- events$four_n_member[member]
    eligible_n <- as.integer(sum(eligible))
    positive_eligible_n <- as.integer(sum(positive & eligible))
    sub_four_n <- as.integer(sum(sub_four))
    four_n_n <- as.integer(sum(four_n))
    positive_sub_four_n <- as.integer(sum(positive & sub_four))
    positive_four_n_n <- as.integer(sum(positive & four_n))
    expected_numerators <- c(
      positive_eligible_n, positive_four_n_n, positive_four_n_n,
      positive_sub_four_n, positive_sub_four_n
    )
    expected_denominators <- c(
      eligible_n, four_n_n, eligible_n, sub_four_n, eligible_n
    )
    if (!identical(metric_rows$numerator_count, expected_numerators) ||
        !identical(metric_rows$denominator_count, expected_denominators)) {
      ph3_background_fail(
        "classification_metric_count_mismatch",
        paste0("classification membership counts changed for acquisition `",
               acquisition_id, "`")
      )
    }
  }
  invisible(events)
}

ph3_background_fit_row <- function(fit, experiment_id, replicate_set_id,
                                   sample_id, selected_for_analysis) {
  data.frame(
    experiment_id = experiment_id, replicate_set_id = replicate_set_id,
    sample_id = sample_id, fit_status = fit$fit_status,
    validity_reason_code = fit$validity_reason_code,
    negative_event_count = fit$negative_event_count,
    positive_event_count = fit$positive_event_count,
    design_rank = fit$design_rank, intercept = fit$intercept,
    slope = fit$slope, negative_dna_min = fit$negative_dna_min,
    negative_dna_max = fit$negative_dna_max,
    positive_dna_min = fit$positive_dna_min,
    positive_dna_max = fit$positive_dna_max,
    coverage_status = fit$coverage_status,
    selected_for_analysis = selected_for_analysis,
    stringsAsFactors = FALSE
  )
}

ph3_not_attempted_pooled_fit <- function() {
  list(
    fit_status = "not_attempted",
    validity_reason_code = "pooled_fit_not_required",
    negative_event_count = NA_integer_, positive_event_count = NA_integer_,
    design_rank = NA_integer_, intercept = NA_real_, slope = NA_real_,
    negative_dna_min = NA_real_, negative_dna_max = NA_real_,
    positive_dna_min = NA_real_, positive_dna_max = NA_real_,
    coverage_status = "not_evaluated", predictions = NULL
  )
}

ph3_background_analysis_id <- function(value) {
  is.character(value) && length(value) == 1L && !is.na(value) &&
    grepl("^ph3-analysis-sha256:[0-9a-f]{64}$", value)
}

ph3_validate_background_model <- function(model) {
  if (!inherits(model, "ph3_output_model") ||
      !is.list(model$schema) ||
      !identical(model$schema$schema_version,
                 "ph3-output-contract-model-1.0.0") ||
      !is.data.frame(model$experiment) || nrow(model$experiment) != 1L ||
      !is.data.frame(model$samples) || !nrow(model$samples) ||
      !is.data.frame(model$replicate_sets) || !nrow(model$replicate_sets) ||
      !is.data.frame(model$correction)) {
    ph3_background_fail(
      "invalid_output_model",
      "Slice 2 requires one validated Slice 1 pH3 output model"
    )
  }
  ph3_validate_model_table(
    model$experiment, model$schema$experiment_columns,
    model$schema$experiment_types, "experiment"
  )
  ph3_validate_model_table(
    model$replicate_sets, model$schema$replicate_set_columns,
    model$schema$replicate_set_types, "replicate_sets"
  )
  ph3_validate_model_table(
    model$samples, model$schema$sample_columns,
    model$schema$sample_types, "samples"
  )
  ph3_validate_model_table(
    model$correction, model$schema$correction_columns,
    model$schema$correction_types, "correction"
  )
  experiment_id <- model$experiment$experiment_id[[1L]]
  analysis_id <- model$experiment$analysis_id[[1L]]
  set_ids <- model$replicate_sets$replicate_set_id
  sample_ids <- model$samples$sample_id
  if (!ph3_output_contract_key(experiment_id) ||
      !ph3_background_analysis_id(analysis_id) ||
      !ph3_output_contract_string_vector(set_ids) ||
      !ph3_output_contract_string_vector(sample_ids) ||
      any(model$replicate_sets$experiment_id != experiment_id) ||
      any(model$samples$experiment_id != experiment_id) ||
      any(model$correction$experiment_id != experiment_id) ||
      !identical(model$correction$replicate_set_id, set_ids) ||
      any(!model$samples$replicate_set_id %in% set_ids) ||
      !all(model$correction$status == "not_computed") ||
      any(!is.na(model$correction$signal_basis))) {
    ph3_background_fail(
      "invalid_model_identity_mapping",
      "experiment, replicate-set, sample, and Slice 1 correction identities must be exact"
    )
  }
  manifest <- model$source$sample_manifest
  metrics <- model$source$quantitation$ph3_metrics_acquisition
  mapping_columns <- c("acquisition_id", "sample_id", "prefix")
  if (!is.data.frame(manifest) ||
      !all(c("prefix", "replicate_set_id") %in% names(manifest)) ||
      anyDuplicated(manifest$prefix) ||
      !is.data.frame(metrics) ||
      !all(mapping_columns %in% names(metrics))) {
    ph3_background_fail(
      "missing_source_mapping",
      "validated manifest and acquisition identity mappings are required"
    )
  }
  acquisition_mapping <- unique(metrics[mapping_columns])
  if (anyDuplicated(acquisition_mapping$acquisition_id) ||
      anyDuplicated(acquisition_mapping$prefix) ||
      !setequal(acquisition_mapping$prefix, manifest$prefix)) {
    ph3_background_fail(
      "invalid_source_mapping",
      "every acquisition and manifest prefix must map exactly once"
    )
  }
  mapped_sets <- manifest$replicate_set_id[
    match(acquisition_mapping$prefix, manifest$prefix)
  ]
  expected_sets <- model$samples$replicate_set_id[
    match(acquisition_mapping$sample_id, model$samples$sample_id)
  ]
  if (anyNA(expected_sets) || !identical(mapped_sets, expected_sets)) {
    ph3_background_fail(
      "invalid_source_mapping",
      "acquisition samples and prefixes must retain the validated replicate-set mapping"
    )
  }
  invisible(model)
}

apply_ph3_background_regression <- function(model) {
  ph3_validate_background_model(model)
  events <- ph3_background_event_table(model)
  experiment_id <- model$experiment$experiment_id[[1L]]
  sample_ids <- model$samples$sample_id
  individual_fits <- lapply(sample_ids, function(sample_id) {
    member <- events$sample_id == sample_id
    eligible <- events$eligible_2to4n[member]
    positive <- eligible & events$ph3_positive_member[member]
    negative <- eligible & !events$ph3_positive_member[member]
    ph3_fit_background_model(
      events$dna_raw[member], events$target_raw[member], negative, positive
    )
  })
  names(individual_fits) <- sample_ids
  replicate_set_ids <- model$replicate_sets$replicate_set_id
  pooled_fits <- setNames(vector("list", length(replicate_set_ids)),
                         replicate_set_ids)
  invalid_samples <- sample_ids[vapply(individual_fits, function(x) {
    !identical(x$fit_status, "valid")
  }, logical(1))]
  individual_trigger_sample_id <- if (length(invalid_samples)) {
    invalid_samples[[1L]]
  } else NA_character_
  individual_trigger_reason_code <- if (length(invalid_samples)) {
    individual_fits[[individual_trigger_sample_id]]$validity_reason_code
  } else NA_character_
  if (!length(invalid_samples)) {
    for (replicate_set_id in replicate_set_ids) {
      pooled_fits[[replicate_set_id]] <- ph3_not_attempted_pooled_fit()
    }
    selected_basis <- "individual_corrected"
    selected_reason <- "all_individual_fits_valid"
    pooled_failure_set_id <- NA_character_
    pooled_failure_reason_code <- NA_character_
  } else {
    for (replicate_set_id in replicate_set_ids) {
      member <- events$replicate_set_id == replicate_set_id
      eligible <- events$eligible_2to4n[member]
      positive <- eligible & events$ph3_positive_member[member]
      negative <- eligible & !events$ph3_positive_member[member]
      pooled_fits[[replicate_set_id]] <- ph3_fit_background_model(
        events$dna_raw[member], events$target_raw[member], negative, positive
      )
    }
    invalid_pooled_sets <- replicate_set_ids[vapply(
      pooled_fits, function(x) !identical(x$fit_status, "valid"), logical(1)
    )]
    if (!length(invalid_pooled_sets)) {
      selected_basis <- "pooled_corrected"
      selected_reason <- "individual_fit_failure_all_pooled_fits_valid"
      pooled_failure_set_id <- NA_character_
      pooled_failure_reason_code <- NA_character_
    } else {
      selected_basis <- "raw"
      selected_reason <- "pooled_fit_failure_experiment_raw_fallback"
      pooled_failure_set_id <- invalid_pooled_sets[[1L]]
      pooled_failure_reason_code <-
        pooled_fits[[pooled_failure_set_id]]$validity_reason_code
    }
  }
  decisions <- do.call(rbind, lapply(replicate_set_ids, function(set_id) {
    data.frame(
      experiment_id = experiment_id, replicate_set_id = set_id,
      signal_basis = selected_basis, reason_code = selected_reason,
      individual_trigger_sample_id = individual_trigger_sample_id,
      individual_trigger_reason_code = individual_trigger_reason_code,
      pooled_failure_replicate_set_id = pooled_failure_set_id,
      pooled_failure_reason_code = pooled_failure_reason_code,
      stringsAsFactors = FALSE
    )
  }))
  rownames(decisions) <- NULL
  individual_rows <- do.call(rbind, lapply(sample_ids, function(sample_id) {
    sample <- model$samples[model$samples$sample_id == sample_id, , drop = FALSE]
    basis <- decisions$signal_basis[
      decisions$replicate_set_id == sample$replicate_set_id
    ]
    ph3_background_fit_row(
      individual_fits[[sample_id]], experiment_id, sample$replicate_set_id,
      sample_id, identical(basis, "individual_corrected")
    )
  }))
  rownames(individual_rows) <- NULL
  pooled_rows <- do.call(rbind, lapply(replicate_set_ids, function(set_id) {
    basis <- decisions$signal_basis[decisions$replicate_set_id == set_id]
    ph3_background_fit_row(
      pooled_fits[[set_id]], experiment_id, set_id, NA_character_,
      identical(basis, "pooled_corrected")
    )
  }))
  rownames(pooled_rows) <- NULL

  events$predicted_background <- NA_real_
  events$analytical_signal <- NA_real_
  events$signal_basis <- decisions$signal_basis[
    match(events$replicate_set_id, decisions$replicate_set_id)
  ]
  events$fit_scope <- ifelse(
    events$signal_basis == "individual_corrected", "individual",
    ifelse(events$signal_basis == "pooled_corrected", "pooled", "not_applicable")
  )
  events$event_signal_status <- ifelse(
    events$eligible_2to4n, "available", "not_applicable"
  )
  events$event_signal_reason_code <- ifelse(
    events$eligible_2to4n, NA_character_, "event_not_eligible_2to4n"
  )
  for (i in seq_along(sample_ids)) {
    sample_id <- sample_ids[[i]]
    member <- events$sample_id == sample_id
    eligible <- member & events$eligible_2to4n
    basis <- unique(events$signal_basis[member])
    if (identical(basis, "raw")) {
      events$analytical_signal[eligible] <- events$target_raw[eligible]
    } else {
      fit <- if (identical(basis, "individual_corrected")) {
        individual_fits[[sample_id]]
      } else pooled_fits[[model$samples$replicate_set_id[
        model$samples$sample_id == sample_id
      ]]]
      predictions <- fit$intercept + fit$slope * events$dna_raw[eligible]
      events$predicted_background[eligible] <- predictions
      events$analytical_signal[eligible] <-
        events$target_raw[eligible] - predictions
    }
    nonfinite <- eligible & !is.finite(events$analytical_signal)
    events$event_signal_status[nonfinite] <- "unavailable"
    events$event_signal_reason_code[nonfinite] <- "nonfinite_raw_or_corrected_signal"
  }
  event_columns <- c(
    "experiment_id", "replicate_set_id", "sample_id", "analysis_id",
    "acquisition_id", "prefix", "event_identity", "event_row",
    "ph3_positive_member", "eligible_2to4n", "sub_4n_member",
    "four_n_member", "dna_raw", "dna_norm", "target_raw",
    "predicted_background",
    "analytical_signal", "signal_basis", "fit_scope",
    "event_signal_status", "event_signal_reason_code"
  )
  event_signals <- events[event_columns]
  correction <- decisions[c(
    "experiment_id", "replicate_set_id", "signal_basis", "reason_code"
  )]
  correction$status <- "selected"
  correction$reason_detail <- if (identical(selected_basis,
                                             "individual_corrected")) {
    rep(NA_character_, nrow(correction))
  } else {
    detail <- paste0(
      "individual_trigger_sample_id=", individual_trigger_sample_id,
      "; individual_trigger_reason_code=", individual_trigger_reason_code
    )
    if (identical(selected_basis, "raw")) {
      detail <- paste0(
        detail, "; pooled_failure_replicate_set_id=", pooled_failure_set_id,
        "; pooled_failure_reason_code=", pooled_failure_reason_code
      )
    }
    rep(detail, nrow(correction))
  }
  correction <- correction[c(
    "experiment_id", "replicate_set_id", "status", "signal_basis",
    "reason_code", "reason_detail"
  )]
  model$correction <- correction
  model$schema$correction_reason_codes <- c(
    "all_individual_fits_valid",
    "individual_fit_failure_all_pooled_fits_valid",
    "pooled_fit_failure_experiment_raw_fallback"
  )
  model$background_regression <- list(
    schema_version = "ph3-background-regression-1.0.0",
    minimum_negative_event_count = 100L,
    set_decisions = decisions, individual_fits = individual_rows,
    pooled_fits = pooled_rows, event_signals = event_signals
  )
  model
}
