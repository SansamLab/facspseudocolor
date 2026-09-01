# pH3 owner-confirmed signal outcomes ----------------------------------------

ph3_signal_fail <- function(reason, detail) {
  stop(
    "PH3 signal-outcome validation failed [", reason, "]: ", detail, ".",
    call. = FALSE
  )
}

ph3_signal_outcome_specification <- function() {
  data.frame(
    outcome_id = c("C", "D"),
    population_id = c("ph3_positive_4n", "ph3_positive_below_4n"),
    stringsAsFactors = FALSE
  )
}

ph3_signal_cutoff_failure_by_acquisition <- function(model, events) {
  cutoff <- model$source$raw_4n_density_cutoff
  classifications <- model$source$event_classifications
  has_classifications <- is.list(classifications) && length(classifications)
  computed_cutoff <- has_classifications && all(vapply(
    classifications,
    function(data) {
      is.data.frame(data) &&
        identical(ph3_background_positivity_method(data),
                  "ph3_raw_4n_density_cutoff_v1")
    },
    logical(1)
  ))

  if (!computed_cutoff && is.null(cutoff)) {
    if (!is.character(events$positivity_call_status) ||
        !identical(events$positivity_call_status,
                   rep("called", nrow(events))) ||
        !is.character(events$positivity_call_reason_code) ||
        !all(is.na(events$positivity_call_reason_code))) {
      ph3_signal_fail(
        "invalid_legacy_positivity_call_state",
        "a legacy non-computed model must retain explicit called membership"
      )
    }
    return(stats::setNames(rep(FALSE, length(unique(events$acquisition_id))),
                           unique(events$acquisition_id)))
  }

  required <- c("schema_version", "method", "records", "application_qc",
                "density_curves")
  required_records <- c("replicate_set_id", "cutoff_status")
  required_application <- c("replicate_set_id", "acquisition_id",
                            "cutoff_status")
  if (!computed_cutoff || !is.list(cutoff) ||
      !identical(names(cutoff), required) ||
      !identical(cutoff$schema_version,
                 "ph3-raw-4n-density-cutoff-1.0.0") ||
      !identical(cutoff$method, ph3_raw_4n_cutoff_method()) ||
      !is.data.frame(cutoff$records) ||
      !is.data.frame(cutoff$application_qc) ||
      !all(required_records %in% names(cutoff$records)) ||
      !all(required_application %in% names(cutoff$application_qc)) ||
      anyNA(cutoff$records$cutoff_status) ||
      anyNA(cutoff$application_qc$cutoff_status) ||
      any(!cutoff$records$cutoff_status %in%
            c("available", "unavailable_cutoff_failure")) ||
      any(!cutoff$application_qc$cutoff_status %in%
            c("available", "unavailable_cutoff_failure")) ||
      anyDuplicated(cutoff$records$replicate_set_id) ||
      anyDuplicated(cutoff$application_qc$acquisition_id) ||
      !setequal(cutoff$application_qc$acquisition_id,
                unique(events$acquisition_id))) {
    ph3_signal_fail(
      "invalid_raw_4n_cutoff_provenance",
      "computed raw-4N positivity requires complete frozen cutoff provenance"
    )
  }
  if (!is.character(events$positivity_call_status) ||
      anyNA(events$positivity_call_status) ||
      !is.character(events$positivity_call_reason_code) ||
      any(!events$positivity_call_status %in% c(
        "called", "unavailable_nonfinite_raw_signal",
        "unavailable_cutoff_failure"
      ))) {
    ph3_signal_fail(
      "invalid_computed_cutoff_call_state",
      "computed raw-4N positivity must retain explicit event call states"
    )
  }

  result <- logical(length(unique(events$acquisition_id)))
  names(result) <- unique(events$acquisition_id)
  for (acquisition_id in names(result)) {
    application <- cutoff$application_qc[
      cutoff$application_qc$acquisition_id == acquisition_id, , drop = FALSE
    ]
    acquisition_events <- events[events$acquisition_id == acquisition_id, ,
                                 drop = FALSE]
    if (nrow(application) != 1L) {
      ph3_signal_fail(
        "raw_4n_cutoff_application_mismatch",
        "each computed acquisition must bind to exactly one cutoff application record"
      )
    }
    record <- cutoff$records[
      cutoff$records$replicate_set_id == application$replicate_set_id[[1L]],
      , drop = FALSE
    ]
    failed <- identical(application$cutoff_status[[1L]],
                         "unavailable_cutoff_failure")
    if (nrow(record) != 1L ||
        !identical(application$replicate_set_id[[1L]],
                   acquisition_events$replicate_set_id[[1L]]) ||
        !identical(application$cutoff_status[[1L]],
                   record$cutoff_status[[1L]]) ||
        (failed && !all(acquisition_events$positivity_call_status ==
                          "unavailable_cutoff_failure")) ||
        (!failed && any(acquisition_events$positivity_call_status ==
                          "unavailable_cutoff_failure"))) {
      ph3_signal_fail(
        "raw_4n_cutoff_application_mismatch",
        "computed cutoff status must exactly reconcile to every acquisition"
      )
    }
    result[[acquisition_id]] <- failed
  }
  result
}

ph3_validate_signal_model <- function(model) {
  required <- c("experiment", "samples", "correction", "reference", "source",
                "background_regression")
  if (!inherits(model, "ph3_output_model") ||
      !all(required %in% names(model)) ||
      !is.data.frame(model$experiment) || nrow(model$experiment) != 1L ||
      !is.data.frame(model$samples) || !is.data.frame(model$correction) ||
      !is.data.frame(model$reference)) {
    ph3_signal_fail(
      "invalid_output_model",
      "Slice 3 requires the validated Slice 1/2 output model"
    )
  }
  required_events <- c(
    "experiment_id", "replicate_set_id", "sample_id", "acquisition_id",
    "event_identity", "eligible_2to4n", "ph3_positive_member",
    "sub_4n_member", "four_n_member", "analytical_signal", "signal_basis",
    "event_signal_status", "event_signal_reason_code", "positivity_call_status",
    "positivity_call_reason_code"
  )
  events <- model$background_regression$event_signals
  if (!is.data.frame(events) || !all(required_events %in% names(events)) ||
      !nrow(events) || anyNA(events$event_identity) ||
      anyDuplicated(events$event_identity) ||
      !is.character(events$experiment_id) || anyNA(events$experiment_id) ||
      !is.character(events$replicate_set_id) || anyNA(events$replicate_set_id) ||
      !is.character(events$sample_id) || anyNA(events$sample_id) ||
      !is.character(events$acquisition_id) || anyNA(events$acquisition_id) ||
      !is.numeric(events$analytical_signal) ||
      !is.logical(events$eligible_2to4n) ||
      !is.logical(events$ph3_positive_member) ||
      !is.logical(events$sub_4n_member) ||
      !is.logical(events$four_n_member) ||
      anyNA(events$eligible_2to4n) || anyNA(events$ph3_positive_member) ||
      anyNA(events$sub_4n_member) || anyNA(events$four_n_member) ||
      any(!events$eligible_2to4n &
          (events$sub_4n_member | events$four_n_member)) ||
      any(events$sub_4n_member & events$four_n_member) ||
      !identical(events$eligible_2to4n,
                 events$sub_4n_member | events$four_n_member)) {
    ph3_signal_fail(
      "invalid_event_signal_provenance",
      "event signals must preserve exact eligible 2N-4N and population membership"
    )
  }
  experiment_id <- model$experiment$experiment_id
  if (!is.character(experiment_id) || length(experiment_id) != 1L ||
      is.na(experiment_id) || !nzchar(experiment_id) ||
      !identical(unique(events$experiment_id), experiment_id) ||
      !all(c("experiment_id", "replicate_set_id", "sample_id", "condition_id") %in%
           names(model$samples)) ||
      any(vapply(model$samples[c(
        "experiment_id", "replicate_set_id", "sample_id", "condition_id"
      )], function(value) !is.character(value) || anyNA(value) || any(!nzchar(value)),
      logical(1)))) {
    ph3_signal_fail(
      "invalid_sample_provenance",
      "experiment and configured sample identity fields must be complete and exact"
    )
  }
  if (!is.character(events$signal_basis) || anyNA(events$signal_basis) ||
      !is.character(model$correction$signal_basis) ||
      anyNA(model$correction$signal_basis) ||
      any(!events$signal_basis %in% c(
        "individual_corrected", "pooled_corrected", "raw"
      )) || any(!model$correction$signal_basis %in% c(
        "individual_corrected", "pooled_corrected", "raw"
      ))) {
    ph3_signal_fail(
      "mixed_or_invalid_signal_basis",
      "every event and replicate set must retain one approved signal basis"
    )
  }
  expected_samples <- unique(events[c("replicate_set_id", "sample_id")])
  actual_samples <- model$samples[c("replicate_set_id", "sample_id")]
  acquisition_mapping <- unique(events[c("acquisition_id", "sample_id")])
  expected_samples <- expected_samples[order(expected_samples$sample_id), ,
                                       drop = FALSE]
  actual_samples <- actual_samples[order(actual_samples$sample_id), ,
                                   drop = FALSE]
  rownames(expected_samples) <- NULL
  rownames(actual_samples) <- NULL
  if (anyDuplicated(actual_samples$sample_id) ||
      anyDuplicated(acquisition_mapping$acquisition_id) ||
      !identical(expected_samples, actual_samples)) {
    ph3_signal_fail(
      "sample_provenance_mismatch",
      "event signals must map exactly to the configured sample table"
    )
  }
  event_sample_mapping <- unique(events[c(
    "experiment_id", "replicate_set_id", "sample_id"
  )])
  configured_sample_mapping <- model$samples[c(
    "experiment_id", "replicate_set_id", "sample_id"
  )]
  event_sample_mapping <- event_sample_mapping[
    order(event_sample_mapping$sample_id), , drop = FALSE
  ]
  configured_sample_mapping <- configured_sample_mapping[
    order(configured_sample_mapping$sample_id), , drop = FALSE
  ]
  rownames(event_sample_mapping) <- NULL
  rownames(configured_sample_mapping) <- NULL
  if (!identical(
    event_sample_mapping,
    configured_sample_mapping
  )) {
    ph3_signal_fail(
      "sample_provenance_mismatch",
      "event experiment and replicate-set identity must match the configured samples"
    )
  }
  source_mapping <- model$source$sample_mapping
  expected_mapping_columns <- c(
    "experiment_id", "replicate_set_id", "sample_id", "condition_id"
  )
  if (!is.data.frame(source_mapping) ||
      !identical(names(source_mapping), expected_mapping_columns) ||
      !identical(
        source_mapping[order(source_mapping$sample_id), , drop = FALSE],
        model$samples[order(model$samples$sample_id), expected_mapping_columns,
                      drop = FALSE]
      )) {
    ph3_signal_fail(
      "condition_provenance_mismatch",
      "sample condition identity must match the retained immutable source mapping"
    )
  }
  expected_sets <- sort(unique(actual_samples$replicate_set_id))
  if (anyDuplicated(model$reference$replicate_set_id) ||
      !identical(sort(model$reference$replicate_set_id), expected_sets) ||
      any(!model$reference$status %in% c("configured", "not_configured"))) {
    ph3_signal_fail(
      "reference_resolution_mismatch",
      "one valid reference record is required for every configured replicate set"
    )
  }
  for (set_id in expected_sets) {
    selected <- model$correction$signal_basis[
      model$correction$replicate_set_id == set_id
    ]
    observed <- unique(events$signal_basis[events$replicate_set_id == set_id])
    if (length(selected) != 1L || length(observed) != 1L ||
        !identical(observed, selected)) {
      ph3_signal_fail(
        "mixed_or_invalid_signal_basis",
        "each replicate set must retain exactly its selected signal basis"
      )
    }
  }
  for (i in seq_len(nrow(model$reference))) {
    reference <- model$reference[i, , drop = FALSE]
    if (identical(reference$status[[1L]], "configured")) {
      matching <- model$samples[
        model$samples$replicate_set_id == reference$replicate_set_id &
          model$samples$sample_id == reference$sample_id &
          model$samples$condition_id == reference$condition_id,
        , drop = FALSE
      ]
      if (nrow(matching) != 1L) {
        ph3_signal_fail(
          "reference_resolution_mismatch",
          "configured reference sample, condition, and replicate set must match exactly"
        )
      }
    } else if (!is.na(reference$sample_id[[1L]]) ||
               !is.na(reference$condition_id[[1L]]) ||
               !identical(reference$reason_code[[1L]], "reference_not_configured")) {
      ph3_signal_fail(
        "reference_resolution_mismatch",
        "an unconfigured reference must retain only the explicit not-configured state"
      )
    }
  }
  if (!is.character(events$event_signal_status) ||
      !is.character(events$event_signal_reason_code) ||
      anyNA(events$event_signal_status) ||
      !all(events$event_signal_status %in% c(
        "available", "unavailable", "not_applicable"
      )) ||
      any(events$eligible_2to4n & events$event_signal_status == "not_applicable") ||
      any(!events$eligible_2to4n & events$event_signal_status != "not_applicable") ||
      any(events$event_signal_status == "available" &
          (!is.finite(events$analytical_signal) |
           !is.na(events$event_signal_reason_code))) ||
      any(events$event_signal_status == "unavailable" &
          (is.finite(events$analytical_signal) |
           is.na(events$event_signal_reason_code) |
           events$event_signal_reason_code !=
             "nonfinite_raw_or_corrected_signal")) ||
      any(events$event_signal_status == "not_applicable" &
          (!is.na(events$analytical_signal) |
           is.na(events$event_signal_reason_code) |
           events$event_signal_reason_code != "event_not_eligible_2to4n"))) {
    ph3_signal_fail(
      "invalid_event_signal_status",
      "event signal status and reason must retain the authoritative availability semantics"
    )
  }
  invisible(events)
}

ph3_signal_acquisition_table <- function(model, events) {
  outcomes <- ph3_signal_outcome_specification()
  samples <- model$samples
  cutoff_failure_by_acquisition <- ph3_signal_cutoff_failure_by_acquisition(
    model, events
  )
  rows <- list()
  row_index <- 0L
  for (i in seq_len(nrow(samples))) {
    sample <- samples[i, , drop = FALSE]
    sample_events <- events[events$sample_id == sample$sample_id, , drop = FALSE]
    acquisition_ids <- unique(sample_events$acquisition_id)
    for (acquisition_id in acquisition_ids) {
      acquisition_events <- sample_events[
        sample_events$acquisition_id == acquisition_id, , drop = FALSE
      ]
      for (j in seq_len(nrow(outcomes))) {
        outcome <- outcomes[j, , drop = FALSE]
        population <- if (identical(outcome$outcome_id, "C")) {
          acquisition_events$four_n_member
        } else acquisition_events$sub_4n_member
        cutoff_failure <- cutoff_failure_by_acquisition[[acquisition_id]]
        member <- acquisition_events$eligible_2to4n &
          acquisition_events$ph3_positive_member & population
        count <- as.integer(sum(member))
        values <- acquisition_events$analytical_signal[member]
        unavailable <- acquisition_events$event_signal_status[member] != "available" |
          !is.finite(values)
        if (cutoff_failure) {
          status <- "unavailable"
          reason <- "unavailable_cutoff_failure"
          value <- NA_real_
        } else if (!count) {
          status <- "unavailable"
          reason <- "no_qualifying_positive_events"
          value <- NA_real_
        } else if (any(unavailable)) {
          status <- "unavailable"
          reason <- "nonfinite_qualifying_positive_signal"
          value <- NA_real_
        } else {
          status <- "available"
          reason <- NA_character_
          value <- stats::median(values)
        }
        row_index <- row_index + 1L
        rows[[row_index]] <- data.frame(
          experiment_id = sample$experiment_id,
          replicate_set_id = sample$replicate_set_id,
          sample_id = sample$sample_id, condition_id = sample$condition_id,
          acquisition_id = acquisition_id, outcome_id = outcome$outcome_id,
          population_id = outcome$population_id,
          signal_basis = unique(acquisition_events$signal_basis),
          qualifying_event_count = count, acquisition_median = value,
          status = status, reason_code = reason,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

ph3_signal_sample_table <- function(model, acquisition) {
  samples <- model$samples
  reference <- model$reference
  rows <- list()
  row_index <- 0L
  for (i in seq_len(nrow(samples))) {
    sample <- samples[i, , drop = FALSE]
    for (outcome_id in c("C", "D")) {
      source <- acquisition[
        acquisition$sample_id == sample$sample_id &
          acquisition$outcome_id == outcome_id, , drop = FALSE
      ]
      valid <- source$status == "available"
      total <- as.integer(nrow(source))
      valid_count <- as.integer(sum(valid))
      cutoff_failure <- isTRUE(any(
        !is.na(source$reason_code) &
          source$reason_code == "unavailable_cutoff_failure"
      ))
      if (cutoff_failure) {
        direct_status <- "unavailable"
        direct_reason <- "unavailable_cutoff_failure"
        direct_value <- NA_real_
      } else if (!valid_count) {
        direct_status <- "unavailable"
        direct_reason <- "no_valid_acquisition_medians"
        direct_value <- NA_real_
      } else {
        direct_status <- if (valid_count == total) "available" else "available_partial_acquisition_coverage"
        direct_reason <- if (valid_count == total) NA_character_ else "one_or_more_acquisition_medians_unavailable"
        direct_value <- mean(source$acquisition_median[valid])
      }
      reference_row <- reference[
        reference$replicate_set_id == sample$replicate_set_id, , drop = FALSE
      ]
      if (nrow(reference_row) != 1L) {
        ph3_signal_fail(
          "reference_resolution_mismatch",
          "exactly one pre-resolved reference record is required per replicate set"
        )
      }
      row_index <- row_index + 1L
      rows[[row_index]] <- data.frame(
        experiment_id = sample$experiment_id,
        replicate_set_id = sample$replicate_set_id,
        sample_id = sample$sample_id, condition_id = sample$condition_id,
        outcome_id = outcome_id,
        signal_basis = unique(source$signal_basis),
        acquisition_count = total, valid_acquisition_count = valid_count,
        direct_value = direct_value, direct_status = direct_status,
        direct_reason_code = direct_reason,
        reference_condition_id = reference_row$condition_id[[1L]],
        reference_sample_id = reference_row$sample_id[[1L]],
        reference_value = NA_real_, reference_ratio = NA_real_,
        reference_status = NA_character_, reference_reason_code = NA_character_,
        stringsAsFactors = FALSE
      )
    }
  }
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  for (i in seq_len(nrow(result))) {
    reference_row <- reference[
      reference$replicate_set_id == result$replicate_set_id[[i]], , drop = FALSE
    ]
    if (identical(reference_row$status[[1L]], "not_configured")) {
      result$reference_status[[i]] <- "not_applicable"
      result$reference_reason_code[[i]] <- "reference_not_configured"
      next
    }
    if (identical(result$direct_reason_code[[i]], "unavailable_cutoff_failure")) {
      result$reference_status[[i]] <- "unavailable"
      result$reference_reason_code[[i]] <- "unavailable_cutoff_failure"
      next
    }
    reference_source <- result[
      result$sample_id == reference_row$sample_id[[1L]] &
        result$outcome_id == result$outcome_id[[i]], , drop = FALSE
    ]
    if (nrow(reference_source) != 1L ||
        !reference_source$direct_status[[1L]] %in% c(
          "available", "available_partial_acquisition_coverage"
        )) {
      result$reference_status[[i]] <- "unavailable"
      result$reference_reason_code[[i]] <- "reference_value_unavailable"
      next
    }
    reference_value <- reference_source$direct_value[[1L]]
    result$reference_value[[i]] <- reference_value
    if (!is.finite(reference_value)) {
      result$reference_status[[i]] <- "unavailable"
      result$reference_reason_code[[i]] <- "reference_value_nonfinite"
    } else if (reference_value == 0) {
      result$reference_status[[i]] <- "unavailable"
      result$reference_reason_code[[i]] <- "reference_value_zero"
    } else if (!result$direct_status[[i]] %in% c(
      "available", "available_partial_acquisition_coverage"
    )) {
      result$reference_status[[i]] <- "unavailable"
      result$reference_reason_code[[i]] <- "sample_value_unavailable"
    } else {
      result$reference_ratio[[i]] <- result$direct_value[[i]] / reference_value
      result$reference_status[[i]] <- if (
        identical(result$direct_status[[i]], "available") &&
          identical(reference_source$direct_status[[1L]], "available")
      ) "available" else "available_partial_acquisition_coverage"
      result$reference_reason_code[[i]] <- if (
        identical(result$reference_status[[i]], "available")
      ) NA_character_ else "sample_or_reference_partial_acquisition_coverage"
    }
  }
  result
}

#' Derive owner-confirmed pH3 signal outcomes C and D
#'
#' This internal Slice 3 helper preserves the validated event memberships and
#' derives only acquisition and biological-sample signal values. It performs no
#' plotting, writing, statistics, or display offset calculation.
derive_ph3_signal_outcomes <- function(model) {
  events <- ph3_validate_signal_model(model)
  acquisition <- ph3_signal_acquisition_table(model, events)
  samples <- ph3_signal_sample_table(model, acquisition)
  model$signal_outcomes <- list(
    schema_version = "ph3-signal-outcomes-1.0.0",
    acquisition = acquisition,
    sample = samples,
    basis_qc = model$correction,
    reference_qc = model$reference
  )
  model
}
