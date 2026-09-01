# Every constructed value in this file is explicitly SYNTHETIC and test-only.
# These fixtures exercise approved accounting and unavailable states, not biology.

synthetic_ph3_signal_events <- function(sample_id, acquisition_id, basis,
                                        four_n_value, below_four_n_value,
                                        nonfinite_four_n = FALSE) {
  values <- c(four_n_value, below_four_n_value, 1, 1)
  if (isTRUE(nonfinite_four_n)) values[[1L]] <- Inf
  data.frame(
    experiment_id = "SYNTHETIC-experiment",
    replicate_set_id = "SYNTHETIC-set-1",
    sample_id = sample_id, acquisition_id = acquisition_id,
    event_identity = paste0(acquisition_id, ":event:", seq_along(values)),
    eligible_2to4n = TRUE,
    ph3_positive_member = c(TRUE, TRUE, FALSE, FALSE),
    sub_4n_member = c(FALSE, TRUE, TRUE, FALSE),
    four_n_member = c(TRUE, FALSE, FALSE, TRUE),
    analytical_signal = values, signal_basis = basis,
    event_signal_status = ifelse(is.finite(values), "available", "unavailable"),
    event_signal_reason_code = ifelse(
      is.finite(values), NA_character_, "nonfinite_raw_or_corrected_signal"
    ),
    # The established pre-computed-cutoff fixture is explicit about its
    # already-called positive membership.  It is not a raw-4N cutoff fixture.
    positivity_call_status = "called",
    positivity_call_reason_code = NA_character_,
    stringsAsFactors = FALSE
  )
}

synthetic_ph3_signal_model <- function(reference = TRUE,
                                       basis = "individual_corrected",
                                       nonfinite_treatment = FALSE) {
  events <- do.call(rbind, list(
    synthetic_ph3_signal_events(
      "SYNTHETIC-reference", "SYNTHETIC-reference-a", basis, 10, 5
    ),
    synthetic_ph3_signal_events(
      "SYNTHETIC-reference", "SYNTHETIC-reference-b", basis, 30, 15
    ),
    synthetic_ph3_signal_events(
      "SYNTHETIC-treatment", "SYNTHETIC-treatment-a", basis, 40, 20,
      nonfinite_treatment
    ),
    synthetic_ph3_signal_events(
      "SYNTHETIC-treatment", "SYNTHETIC-treatment-b", basis, 60, 30
    )
  ))
  samples <- data.frame(
    experiment_id = "SYNTHETIC-experiment",
    replicate_set_id = "SYNTHETIC-set-1",
    sample_id = c("SYNTHETIC-reference", "SYNTHETIC-treatment"),
    condition_id = c("SYNTHETIC-control", "SYNTHETIC-treatment"),
    acquisition_count = c(2L, 2L),
    source_acquisition_ids = c(
      "[\"SYNTHETIC-reference-a\",\"SYNTHETIC-reference-b\"]",
      "[\"SYNTHETIC-treatment-a\",\"SYNTHETIC-treatment-b\"]"
    ), stringsAsFactors = FALSE
  )
  reference_table <- data.frame(
    experiment_id = "SYNTHETIC-experiment",
    replicate_set_id = "SYNTHETIC-set-1",
    status = if (reference) "configured" else "not_configured",
    condition_id = if (reference) "SYNTHETIC-control" else NA_character_,
    sample_id = if (reference) "SYNTHETIC-reference" else NA_character_,
    reason_code = if (reference) NA_character_ else "reference_not_configured",
    reason_detail = NA_character_, stringsAsFactors = FALSE
  )
  structure(list(
    experiment = data.frame(
      experiment_id = "SYNTHETIC-experiment", analysis_id = "SYNTHETIC-analysis",
      stringsAsFactors = FALSE
    ),
    samples = samples,
    correction = data.frame(
      experiment_id = "SYNTHETIC-experiment", replicate_set_id = "SYNTHETIC-set-1",
      status = "selected", signal_basis = basis,
      reason_code = "SYNTHETIC", reason_detail = NA_character_,
      stringsAsFactors = FALSE
    ),
    reference = reference_table,
    source = list(sample_mapping = samples[c(
      "experiment_id", "replicate_set_id", "sample_id", "condition_id"
    )]),
    background_regression = list(event_signals = events)
  ), class = "ph3_output_model")
}

test_that("SYNTHETIC C/D outcomes use equal acquisition-median weighting", {
  model <- synthetic_ph3_signal_model()
  events_before <- serialize(model$background_regression$event_signals, NULL,
                             version = 3L)
  result <- derive_ph3_signal_outcomes(model)
  acquisition <- result$signal_outcomes$acquisition
  samples <- result$signal_outcomes$sample
  expect_identical(result$signal_outcomes$schema_version,
                   "ph3-signal-outcomes-1.0.0")
  expect_identical(serialize(model$background_regression$event_signals, NULL,
                             version = 3L), events_before)
  expect_equal(acquisition$acquisition_median[
    acquisition$acquisition_id == "SYNTHETIC-reference-a" &
      acquisition$outcome_id == "C"
  ], 10)
  expect_equal(samples$direct_value[
    samples$sample_id == "SYNTHETIC-reference" & samples$outcome_id == "C"
  ], 20)
  expect_equal(samples$direct_value[
    samples$sample_id == "SYNTHETIC-treatment" & samples$outcome_id == "C"
  ], 50)
  expect_equal(samples$reference_ratio[
    samples$sample_id == "SYNTHETIC-treatment" & samples$outcome_id == "C"
  ], 2.5)
  expect_equal(samples$reference_ratio[
    samples$sample_id == "SYNTHETIC-treatment" & samples$outcome_id == "D"
  ], 2.5)
  expect_true(all(samples$signal_basis == "individual_corrected"))
})

test_that("SYNTHETIC legacy called membership does not mistake NA reasons for cutoff failure", {
  model <- synthetic_ph3_signal_model()
  failure <- ph3_signal_cutoff_failure_by_acquisition(
    model, model$background_regression$event_signals
  )
  expect_identical(
    failure,
    stats::setNames(
      c(FALSE, FALSE, FALSE, FALSE),
      c(
        "SYNTHETIC-reference-a", "SYNTHETIC-reference-b",
        "SYNTHETIC-treatment-a", "SYNTHETIC-treatment-b"
      )
    )
  )

  computed_without_cutoff <- model
  computed_without_cutoff$source$event_classifications <- list(data.frame(
    positivity_method_id = "ph3_raw_4n_density_cutoff_v1",
    stringsAsFactors = FALSE
  ))
  expect_error(
    ph3_signal_cutoff_failure_by_acquisition(
      computed_without_cutoff,
      computed_without_cutoff$background_regression$event_signals
    ),
    "invalid_raw_4n_cutoff_provenance"
  )
})

test_that("SYNTHETIC unavailable acquisition is explicit and produces partial coverage", {
  result <- derive_ph3_signal_outcomes(
    synthetic_ph3_signal_model(nonfinite_treatment = TRUE)
  )
  acquisition <- result$signal_outcomes$acquisition
  samples <- result$signal_outcomes$sample
  unavailable <- acquisition[
    acquisition$acquisition_id == "SYNTHETIC-treatment-a" &
      acquisition$outcome_id == "C", , drop = FALSE
  ]
  expect_identical(unavailable$status, "unavailable")
  expect_identical(unavailable$reason_code,
                   "nonfinite_qualifying_positive_signal")
  partial <- samples[
    samples$sample_id == "SYNTHETIC-treatment" & samples$outcome_id == "C",
    , drop = FALSE
  ]
  expect_identical(partial$direct_status,
                   "available_partial_acquisition_coverage")
  expect_equal(partial$direct_value, 60)
  expect_identical(partial$reference_status,
                   "available_partial_acquisition_coverage")
})

test_that("SYNTHETIC partial reference coverage retains an explicitly partial ratio", {
  model <- synthetic_ph3_signal_model()
  model$background_regression$event_signals$analytical_signal[
    model$background_regression$event_signals$acquisition_id ==
      "SYNTHETIC-reference-a" &
      model$background_regression$event_signals$four_n_member
  ] <- Inf
  model$background_regression$event_signals$event_signal_status[
    model$background_regression$event_signals$acquisition_id ==
      "SYNTHETIC-reference-a" &
      model$background_regression$event_signals$four_n_member
  ] <- "unavailable"
  model$background_regression$event_signals$event_signal_reason_code[
    model$background_regression$event_signals$acquisition_id ==
      "SYNTHETIC-reference-a" &
      model$background_regression$event_signals$four_n_member
  ] <- "nonfinite_raw_or_corrected_signal"
  result <- derive_ph3_signal_outcomes(model)$signal_outcomes$sample
  treatment <- result[
    result$sample_id == "SYNTHETIC-treatment" & result$outcome_id == "C",
    , drop = FALSE
  ]
  expect_identical(treatment$reference_status,
                   "available_partial_acquisition_coverage")
  expect_identical(treatment$reference_reason_code,
                   "sample_or_reference_partial_acquisition_coverage")
  expect_equal(treatment$reference_value, 30)
  expect_equal(treatment$reference_ratio, 50 / 30)
})

test_that("SYNTHETIC reference is optional and invalid denominators are unavailable", {
  no_reference <- derive_ph3_signal_outcomes(
    synthetic_ph3_signal_model(reference = FALSE)
  )$signal_outcomes$sample
  expect_true(all(no_reference$reference_status == "not_applicable"))
  expect_true(all(is.na(no_reference$reference_ratio)))

  zero_reference <- synthetic_ph3_signal_model()
  events <- zero_reference$background_regression$event_signals
  events$analytical_signal[
    events$sample_id == "SYNTHETIC-reference" & events$four_n_member
  ] <- 0
  zero_reference$background_regression$event_signals <- events
  zero_result <- derive_ph3_signal_outcomes(zero_reference)$signal_outcomes$sample
  row <- zero_result[
    zero_result$sample_id == "SYNTHETIC-treatment" &
      zero_result$outcome_id == "C", , drop = FALSE
  ]
  expect_identical(row$reference_status, "unavailable")
  expect_identical(row$reference_reason_code, "reference_value_zero")
  expect_true(is.na(row$reference_ratio))
})

test_that("SYNTHETIC mixed signal basis is rejected before any C/D calculation", {
  model <- synthetic_ph3_signal_model()
  model$background_regression$event_signals$signal_basis[[1L]] <- "raw"
  expect_error(
    derive_ph3_signal_outcomes(model), "mixed_or_invalid_signal_basis"
  )
})

test_that("SYNTHETIC tampered provenance and status are rejected fail closed", {
  wrong_reference <- synthetic_ph3_signal_model()
  wrong_reference$reference$sample_id <- "SYNTHETIC-treatment"
  expect_error(
    derive_ph3_signal_outcomes(wrong_reference), "reference_resolution_mismatch"
  )

  wrong_experiment <- synthetic_ph3_signal_model()
  wrong_experiment$background_regression$event_signals$experiment_id[[1L]] <-
    "SYNTHETIC-other-experiment"
  expect_error(
    derive_ph3_signal_outcomes(wrong_experiment), "invalid_sample_provenance"
  )

  wrong_status <- synthetic_ph3_signal_model()
  wrong_status$background_regression$event_signals$event_signal_status[[1L]] <-
    "SYNTHETIC-unknown"
  expect_error(
    derive_ph3_signal_outcomes(wrong_status), "invalid_event_signal_status"
  )

  tampered_condition <- synthetic_ph3_signal_model()
  tampered_condition$samples$condition_id[[1L]] <- "SYNTHETIC-wrong"
  expect_error(
    derive_ph3_signal_outcomes(tampered_condition),
    "condition_provenance_mismatch"
  )

  finite_unavailable <- synthetic_ph3_signal_model()
  finite_unavailable$background_regression$event_signals$
    event_signal_status[[1L]] <- "unavailable"
  finite_unavailable$background_regression$event_signals$
    event_signal_reason_code[[1L]] <- "nonfinite_raw_or_corrected_signal"
  expect_error(
    derive_ph3_signal_outcomes(finite_unavailable),
    "invalid_event_signal_status"
  )
})
