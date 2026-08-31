# pH3 report-ready condition model ------------------------------------------

ph3_report_model_fail <- function(reason, detail) {
  stop(
    "PH3 report-model validation failed [", reason, "]: ", detail, ".",
    call. = FALSE
  )
}

ph3_report_model_columns <- function() {
  list(
    biological_replicate_values = c(
      "experiment_id", "analysis_id", "condition_id", "condition_label",
      "condition_index", "replicate_set_id", "replicate_label",
      "replicate_index", "outcome_id", "outcome_label", "value_kind",
      "value_mode", "signal_basis", "value", "value_status", "reason_code",
      "technical_acquisition_count", "valid_acquisition_count",
      "source_acquisition_ids"
    ),
    condition_summary = c(
      "experiment_id", "analysis_id", "condition_id", "condition_label",
      "condition_index", "outcome_id", "outcome_label", "value_kind",
      "value_mode", "signal_basis", "mean_value", "sd_value",
      "biological_replicate_count", "finite_biological_replicate_count",
      "undefined_biological_replicate_count", "summary_status", "reason_code",
      "source_replicate_set_ids"
    ),
    signal_basis_strata = c(
      "experiment_id", "analysis_id", "condition_id", "condition_label",
      "condition_index", "outcome_id", "outcome_label", "signal_basis",
      "value_mode", "mean_value", "sd_value", "biological_replicate_count",
      "finite_biological_replicate_count", "undefined_biological_replicate_count",
      "summary_status", "reason_code", "source_replicate_set_ids"
    ),
    qc_flags = c(
      "experiment_id", "analysis_id", "condition_id", "outcome_id",
      "signal_basis_status", "summary_status", "reference_mode",
      "reason_code"
    )
  )
}

ph3_report_model_empty <- function(columns) {
  data.frame(
    setNames(replicate(length(columns), character(), simplify = FALSE), columns),
    stringsAsFactors = FALSE
  )
}

ph3_report_model_validate_table <- function(data, columns, context) {
  if (!is.data.frame(data) || !identical(names(data), columns)) {
    ph3_report_model_fail(
      "invalid_report_schema",
      paste0("`", context, "` must retain its exact ordered schema")
    )
  }
  invisible(data)
}

ph3_report_model_prevalence_values <- function(model) {
  source <- model$source$quantitation$ph3_metrics_biological_replicate
  required <- c(
    "analysis_id", "condition", "condition_index", "replicate",
    "replicate_index", "metric_id", "value_percent", "result_status",
    "technical_acquisition_count", "finite_technical_acquisition_count",
    "source_acquisition_ids"
  )
  metric_map <- c(
    A = "ph3_4n_positive_prevalence_within_2to4n_percent",
    B = "ph3_sub_4n_positive_prevalence_within_2to4n_percent"
  )
  if (!is.data.frame(source) || any(!required %in% names(source))) {
    ph3_report_model_fail(
      "invalid_prevalence_source", "validated Slice 5 biological-replicate metrics are required"
    )
  }
  condition <- model$conditions
  sets <- model$replicate_sets
  outcomes <- model$outcomes
  rows <- list()
  index <- 0L
  for (outcome_id in names(metric_map)) {
    selected <- source[source$metric_id == metric_map[[outcome_id]], , drop = FALSE]
    if (nrow(selected) != nrow(condition) * nrow(sets) ||
        anyDuplicated(selected[c("condition_index", "replicate_index")])) {
      ph3_report_model_fail(
        "invalid_prevalence_source", "A/B values must contain exactly one row per condition and replicate set"
      )
    }
    for (i in seq_len(nrow(selected))) {
      value <- selected[i, , drop = FALSE]
      condition_row <- condition[condition$condition_index == value$condition_index,
                                 , drop = FALSE]
      set_row <- sets[sets$manifest_replicate_index == value$replicate_index,
                      , drop = FALSE]
      outcome <- outcomes[outcomes$outcome_id == outcome_id, , drop = FALSE]
      if (nrow(condition_row) != 1L || nrow(set_row) != 1L ||
          nrow(outcome) != 1L ||
          !identical(value$analysis_id[[1L]], model$experiment$analysis_id[[1L]]) ||
          !identical(value$condition[[1L]], condition_row$condition_label[[1L]]) ||
          !identical(value$replicate[[1L]], set_row$replicate_label[[1L]]) ||
          !value$result_status[[1L]] %in% c("ok", "ok_partial_undefined", "undefined_no_finite_values")) {
        ph3_report_model_fail(
          "prevalence_provenance_mismatch", "A/B values must retain their configured condition and replicate-set provenance"
        )
      }
      available <- is.finite(value$value_percent[[1L]])
      if (available != value$result_status[[1L]] %in% c("ok", "ok_partial_undefined")) {
        ph3_report_model_fail(
          "invalid_prevalence_status", "A/B finite values and statuses must agree"
        )
      }
      index <- index + 1L
      rows[[index]] <- data.frame(
        experiment_id = model$experiment$experiment_id[[1L]],
        analysis_id = model$experiment$analysis_id[[1L]],
        condition_id = condition_row$condition_id[[1L]],
        condition_label = condition_row$condition_label[[1L]],
        condition_index = as.integer(condition_row$condition_index[[1L]]),
        replicate_set_id = set_row$replicate_set_id[[1L]],
        replicate_label = set_row$replicate_label[[1L]],
        replicate_index = as.integer(set_row$manifest_replicate_index[[1L]]),
        outcome_id = outcome_id, outcome_label = outcome$label[[1L]],
        value_kind = outcome$value_kind[[1L]], value_mode = "prevalence_percent",
        signal_basis = "not_applicable", value = value$value_percent[[1L]],
        value_status = if (!available) {
          "unavailable"
        } else if (identical(value$result_status[[1L]], "ok_partial_undefined")) {
          "available_partial_acquisition_coverage"
        } else {
          "available"
        },
        reason_code = if (!available) {
          "undefined_zero_denominator"
        } else if (identical(value$result_status[[1L]], "ok_partial_undefined")) {
          "one_or_more_technical_acquisition_values_undefined"
        } else {
          NA_character_
        },
        technical_acquisition_count = as.integer(value$technical_acquisition_count[[1L]]),
        valid_acquisition_count = as.integer(value$finite_technical_acquisition_count[[1L]]),
        source_acquisition_ids = value$source_acquisition_ids[[1L]],
        stringsAsFactors = FALSE
      )
    }
  }
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

ph3_report_model_signal_values <- function(model) {
  outcomes <- model$signal_outcomes
  if (!is.list(outcomes) || !identical(names(outcomes), c(
    "schema_version", "acquisition", "sample", "basis_qc", "reference_qc"
  )) || !identical(outcomes$schema_version, "ph3-signal-outcomes-1.0.0") ||
      !is.data.frame(outcomes$sample) ||
      !identical(outcomes$basis_qc, model$correction) ||
      !identical(outcomes$reference_qc, model$reference)) {
    ph3_report_model_fail(
      "missing_signal_outcomes", "validated Slice 3 signal outcomes are required"
    )
  }
  required <- c(
    "experiment_id", "replicate_set_id", "sample_id", "condition_id",
    "outcome_id", "signal_basis", "acquisition_count",
    "valid_acquisition_count", "direct_value", "direct_status",
    "direct_reason_code", "reference_ratio", "reference_status",
    "reference_reason_code"
  )
  values <- outcomes$sample
  if (any(!required %in% names(values)) || any(!values$outcome_id %in% c("C", "D")) ||
      anyDuplicated(values[c("sample_id", "outcome_id")])) {
    ph3_report_model_fail(
      "invalid_signal_source", "C/D values must retain one complete row per sample and outcome"
    )
  }
  samples <- model$samples
  condition <- model$conditions
  sets <- model$replicate_sets
  outcomes_spec <- model$outcomes
  reference_configured <- all(model$reference$status == "configured")
  reference_unconfigured <- all(model$reference$status == "not_configured")
  if (!reference_configured && !reference_unconfigured) {
    ph3_report_model_fail(
      "mixed_reference_configuration", "reference configuration must be experiment-wide and explicit"
    )
  }
  rows <- lapply(seq_len(nrow(values)), function(i) {
    value <- values[i, , drop = FALSE]
    sample <- samples[samples$sample_id == value$sample_id, , drop = FALSE]
    condition_row <- condition[condition$condition_id == value$condition_id,
                               , drop = FALSE]
    set_row <- sets[sets$replicate_set_id == value$replicate_set_id,
                    , drop = FALSE]
    outcome <- outcomes_spec[outcomes_spec$outcome_id == value$outcome_id,
                             , drop = FALSE]
    if (nrow(sample) != 1L || nrow(condition_row) != 1L || nrow(set_row) != 1L ||
        nrow(outcome) != 1L || !identical(sample$condition_id[[1L]], value$condition_id[[1L]]) ||
        !identical(sample$replicate_set_id[[1L]], value$replicate_set_id[[1L]]) ||
        !identical(value$experiment_id[[1L]], model$experiment$experiment_id[[1L]])) {
      ph3_report_model_fail(
        "signal_provenance_mismatch", "C/D values must map exactly to validated sample provenance"
      )
    }
    correction <- model$correction[
      model$correction$replicate_set_id == value$replicate_set_id, , drop = FALSE
    ]
    if (nrow(correction) != 1L || !identical(correction$status[[1L]], "selected") ||
        !identical(value$signal_basis[[1L]], correction$signal_basis[[1L]])) {
      ph3_report_model_fail(
        "signal_basis_provenance_mismatch", "C/D values must retain the selected correction basis for their replicate set"
      )
    }
    mode <- if (reference_configured) "reference_ratio" else "direct_median"
    analytical_value <- if (reference_configured) value$reference_ratio[[1L]] else value$direct_value[[1L]]
    status <- if (reference_configured) value$reference_status[[1L]] else value$direct_status[[1L]]
    reason <- if (reference_configured) value$reference_reason_code[[1L]] else value$direct_reason_code[[1L]]
    available <- status %in% c("available", "available_partial_acquisition_coverage")
    if (available != is.finite(analytical_value) ||
        !status %in% c("available", "available_partial_acquisition_coverage", "unavailable")) {
      ph3_report_model_fail(
        "invalid_signal_status", "C/D analytical value, status, and selected reference mode must agree"
      )
    }
    data.frame(
      experiment_id = value$experiment_id[[1L]], analysis_id = model$experiment$analysis_id[[1L]],
      condition_id = value$condition_id[[1L]], condition_label = condition_row$condition_label[[1L]],
      condition_index = as.integer(condition_row$condition_index[[1L]]),
      replicate_set_id = value$replicate_set_id[[1L]], replicate_label = set_row$replicate_label[[1L]],
      replicate_index = as.integer(set_row$manifest_replicate_index[[1L]]),
      outcome_id = value$outcome_id[[1L]], outcome_label = outcome$label[[1L]],
      value_kind = outcome$value_kind[[1L]], value_mode = mode,
      signal_basis = value$signal_basis[[1L]], value = analytical_value,
      value_status = if (available) status else "unavailable", reason_code = reason,
      technical_acquisition_count = as.integer(value$acquisition_count[[1L]]),
      valid_acquisition_count = as.integer(value$valid_acquisition_count[[1L]]),
      source_acquisition_ids = sample$source_acquisition_ids[[1L]],
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  if (nrow(result) != nrow(samples) * 2L ||
      anyDuplicated(result[c("replicate_set_id", "condition_id", "outcome_id")])) {
    ph3_report_model_fail(
      "invalid_signal_source", "C/D report values must retain exactly one sample per condition and replicate set"
    )
  }
  result
}

ph3_report_model_summary_row <- function(data, by_basis = FALSE) {
  finite <- is.finite(data$value)
  total <- as.integer(nrow(data))
  finite_n <- as.integer(sum(finite))
  undefined_n <- as.integer(total - finite_n)
  bases <- unique(data$signal_basis)
  signal <- data$outcome_id[[1L]] %in% c("C", "D")
  mixed <- signal && length(bases) != 1L
  status <- if (mixed) "unavailable_mixed_signal_basis" else if (!finite_n) {
    "unavailable_no_finite_values"
  } else if (undefined_n) {
    "available_partial_biological_replicate_coverage"
  } else {
    "available"
  }
  list(
    mean = if (mixed || !finite_n) NA_real_ else mean(data$value[finite]),
    sd = if (mixed || finite_n < 3L) NA_real_ else stats::sd(data$value[finite]),
    total = total, finite_n = finite_n, undefined_n = undefined_n,
    status = status,
    reason = if (mixed) "mixed_signal_basis" else if (!finite_n) {
      "no_finite_biological_replicate_values"
    } else if (undefined_n) "one_or_more_biological_replicate_values_unavailable" else NA_character_,
    basis = if (signal && length(bases) == 1L) bases[[1L]] else if (signal) "mixed" else "not_applicable"
  )
}

#' Build validated report-ready pH3 condition summaries
#'
#' This internal helper derives only the four owner-confirmed panel data model.
#' It does not render plots, write files, calculate display offsets, or run
#' inferential statistics.
derive_ph3_condition_report_model <- function(model) {
  if (!inherits(model, "ph3_output_model") || is.null(model$signal_outcomes) ||
      !is.data.frame(model$outcomes) || !identical(model$outcomes$outcome_id, c("A", "B", "C", "D"))) {
    ph3_report_model_fail(
      "invalid_output_model", "Slice 4 requires the validated Slice 1-3 output model"
    )
  }
  columns <- ph3_report_model_columns()
  values <- rbind(
    ph3_report_model_prevalence_values(model),
    ph3_report_model_signal_values(model)
  )
  values <- values[order(values$condition_index, values$replicate_index,
                         match(values$outcome_id, c("A", "B", "C", "D"))), , drop = FALSE]
  rownames(values) <- NULL
  ph3_report_model_validate_table(values, columns$biological_replicate_values,
                                  "biological_replicate_values")
  groups <- unique(values[c("condition_id", "condition_label", "condition_index",
                            "outcome_id", "outcome_label", "value_kind", "value_mode")])
  groups <- groups[order(groups$condition_index, match(groups$outcome_id, c("A", "B", "C", "D"))), , drop = FALSE]
  summaries <- lapply(seq_len(nrow(groups)), function(i) {
    group <- groups[i, , drop = FALSE]
    selected <- values[values$condition_id == group$condition_id &
                         values$outcome_id == group$outcome_id, , drop = FALSE]
    result <- ph3_report_model_summary_row(selected)
    data.frame(
      experiment_id = model$experiment$experiment_id[[1L]], analysis_id = model$experiment$analysis_id[[1L]],
      condition_id = group$condition_id[[1L]], condition_label = group$condition_label[[1L]],
      condition_index = as.integer(group$condition_index[[1L]]), outcome_id = group$outcome_id[[1L]],
      outcome_label = group$outcome_label[[1L]], value_kind = group$value_kind[[1L]],
      value_mode = group$value_mode[[1L]], signal_basis = result$basis,
      mean_value = result$mean, sd_value = result$sd,
      biological_replicate_count = result$total, finite_biological_replicate_count = result$finite_n,
      undefined_biological_replicate_count = result$undefined_n,
      summary_status = result$status, reason_code = result$reason,
      source_replicate_set_ids = ph3_canonical_string_array(selected$replicate_set_id),
      stringsAsFactors = FALSE
    )
  })
  summaries <- do.call(rbind, summaries)
  rownames(summaries) <- NULL
  ph3_report_model_validate_table(summaries, columns$condition_summary, "condition_summary")
  signal_values <- values[values$outcome_id %in% c("C", "D"), , drop = FALSE]
  strata_groups <- unique(signal_values[c("condition_id", "condition_label", "condition_index",
                                          "outcome_id", "outcome_label", "signal_basis", "value_mode")])
  strata_groups <- strata_groups[order(strata_groups$condition_index,
                                       match(strata_groups$outcome_id, c("C", "D")),
                                       strata_groups$signal_basis), , drop = FALSE]
  strata <- lapply(seq_len(nrow(strata_groups)), function(i) {
    group <- strata_groups[i, , drop = FALSE]
    selected <- signal_values[signal_values$condition_id == group$condition_id &
                                signal_values$outcome_id == group$outcome_id &
                                signal_values$signal_basis == group$signal_basis, , drop = FALSE]
    result <- ph3_report_model_summary_row(selected, by_basis = TRUE)
    data.frame(
      experiment_id = model$experiment$experiment_id[[1L]], analysis_id = model$experiment$analysis_id[[1L]],
      condition_id = group$condition_id[[1L]], condition_label = group$condition_label[[1L]],
      condition_index = as.integer(group$condition_index[[1L]]), outcome_id = group$outcome_id[[1L]],
      outcome_label = group$outcome_label[[1L]], signal_basis = group$signal_basis[[1L]],
      value_mode = group$value_mode[[1L]], mean_value = result$mean, sd_value = result$sd,
      biological_replicate_count = result$total, finite_biological_replicate_count = result$finite_n,
      undefined_biological_replicate_count = result$undefined_n,
      summary_status = result$status, reason_code = result$reason,
      source_replicate_set_ids = ph3_canonical_string_array(selected$replicate_set_id),
      stringsAsFactors = FALSE
    )
  })
  strata <- if (length(strata)) do.call(rbind, strata) else ph3_report_model_empty(columns$signal_basis_strata)
  rownames(strata) <- NULL
  ph3_report_model_validate_table(strata, columns$signal_basis_strata, "signal_basis_strata")
  flags <- summaries[c("experiment_id", "analysis_id", "condition_id", "outcome_id", "summary_status")]
  flags$signal_basis_status <- ifelse(flags$outcome_id %in% c("C", "D"),
                                       ifelse(flags$summary_status == "unavailable_mixed_signal_basis", "mixed", "single_or_unavailable"),
                                       "not_applicable")
  flags$reference_mode <- summaries$value_mode
  flags$reason_code <- summaries$reason_code
  flags <- flags[c("experiment_id", "analysis_id", "condition_id", "outcome_id",
                   "signal_basis_status", "summary_status", "reference_mode", "reason_code")]
  ph3_report_model_validate_table(flags, columns$qc_flags, "qc_flags")
  model$condition_report <- list(
    schema_version = "ph3-condition-report-model-1.0.0",
    biological_replicate_values = values,
    condition_summary = summaries,
    signal_basis_strata = strata,
    qc_flags = flags,
    provenance = list(
      source_output_model_schema_version = model$schema$schema_version,
      source_signal_outcomes_schema_version = model$signal_outcomes$schema_version,
      analysis_id = model$experiment$analysis_id[[1L]],
      experiment_id = model$experiment$experiment_id[[1L]]
    )
  )
  model
}
