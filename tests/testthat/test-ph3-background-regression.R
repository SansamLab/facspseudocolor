# Every constructed event below is explicitly SYNTHETIC and test-only.
# These fixtures exercise numerical validity and fallback, not biology.

synthetic_ph3_background_analysis_id <- paste0(
  "ph3-analysis-sha256:", paste(rep("a", 64L), collapse = "")
)

synthetic_ph3_background_events <- function(sample_id, acquisition_id,
                                            negative_n = 100L,
                                            negative_dna = NULL,
                                            positive_dna = c(0, 10),
                                            nonfinite_negative = FALSE,
                                            nonfinite_positive_signal = FALSE) {
  if (is.null(negative_dna)) {
    negative_dna <- seq(0, 10, length.out = negative_n)
  } else {
    negative_n <- length(negative_dna)
  }
  positive_n <- length(positive_dna)
  dna <- c(negative_dna, positive_dna)
  target <- c(2 + 3 * negative_dna, 12 + 3 * positive_dna)
  if (nonfinite_negative) target[[1L]] <- Inf
  if (nonfinite_positive_signal) target[[negative_n + 1L]] <- Inf
  positive <- c(rep(FALSE, negative_n), rep(TRUE, positive_n))
  data.frame(
    analysis_id = synthetic_ph3_background_analysis_id,
    acquisition_id = acquisition_id, sample_id = sample_id,
    prefix = acquisition_id,
    event_identity = paste0(acquisition_id, ":event:", seq_along(dna)),
    identity_valid = TRUE, dna_finite = is.finite(dna),
    ph3_positive_member = positive, eligible_2to4n = TRUE,
    sub_4n_member = dna < 5, four_n_member = dna >= 5,
    dna_raw = dna, dna_norm = dna, target_raw = target,
    stringsAsFactors = FALSE
  )
}

synthetic_ph3_background_model <- function(specifications) {
  samples <- do.call(rbind, lapply(seq_along(specifications), function(i) {
    specification <- specifications[[i]]
    data.frame(
      experiment_id = "SYNTHETIC-experiment",
      replicate_set_id = specification$replicate_set_id,
      sample_id = specification$sample_id,
      condition_id = paste0("SYNTHETIC-condition-", i),
      acquisition_count = 1L,
      source_acquisition_ids = paste0("[\"SYNTHETIC-acquisition-", i, "\"]"),
      stringsAsFactors = FALSE
    )
  }))
  replicate_set_ids <- unique(samples$replicate_set_id)
  classifications <- lapply(seq_along(specifications), function(i) {
    specification <- specifications[[i]]
    do.call(
      synthetic_ph3_background_events,
      c(list(
        sample_id = specification$sample_id,
        acquisition_id = paste0("SYNTHETIC-acquisition-", i)
      ), specification$event_arguments)
    )
  })
  correction <- data.frame(
    experiment_id = "SYNTHETIC-experiment",
    replicate_set_id = replicate_set_ids, status = "not_computed",
    signal_basis = NA_character_,
    reason_code = "background_correction_not_computed_slice_1",
    reason_detail = NA_character_, stringsAsFactors = FALSE
  )
  sample_manifest <- do.call(rbind, lapply(seq_along(specifications), function(i) {
    data.frame(
      prefix = paste0("SYNTHETIC-acquisition-", i),
      replicate_set_id = specifications[[i]]$replicate_set_id,
      stringsAsFactors = FALSE
    )
  }))
  acquisition_metrics <- do.call(rbind, lapply(seq_along(specifications), function(i) {
    classification <- classifications[[i]]
    eligible <- classification$eligible_2to4n
    positive <- classification$ph3_positive_member
    sub_four <- classification$sub_4n_member
    four_n <- classification$four_n_member
    eligible_n <- as.integer(sum(eligible))
    positive_eligible_n <- as.integer(sum(positive & eligible))
    sub_four_n <- as.integer(sum(sub_four))
    four_n_n <- as.integer(sum(four_n))
    positive_sub_four_n <- as.integer(sum(positive & sub_four))
    positive_four_n_n <- as.integer(sum(positive & four_n))
    data.frame(
      acquisition_id = paste0("SYNTHETIC-acquisition-", i),
      sample_id = specifications[[i]]$sample_id,
      prefix = paste0("SYNTHETIC-acquisition-", i),
      metric_id = c(
        "ph3_2to4n_positivity_percent",
        "ph3_within_4n_positivity_percent",
        "ph3_4n_positive_prevalence_within_2to4n_percent",
        "ph3_within_sub_4n_positivity_percent",
        "ph3_sub_4n_positive_prevalence_within_2to4n_percent"
      ),
      numerator_count = c(
        positive_eligible_n, positive_four_n_n, positive_four_n_n,
        positive_sub_four_n, positive_sub_four_n
      ),
      denominator_count = c(
        eligible_n, four_n_n, eligible_n, sub_four_n, eligible_n
      ),
      stringsAsFactors = FALSE
    )
  }))
  structure(list(
    schema = ph3_output_model_schema(),
    experiment = data.frame(
      experiment_id = "SYNTHETIC-experiment",
      analysis_id = synthetic_ph3_background_analysis_id,
      stringsAsFactors = FALSE
    ),
    replicate_sets = data.frame(
      experiment_id = "SYNTHETIC-experiment",
      replicate_set_id = replicate_set_ids,
      replicate_label = paste0("SYNTHETIC set ", seq_along(replicate_set_ids)),
      manifest_replicate_index = as.integer(seq_along(replicate_set_ids)),
      stringsAsFactors = FALSE
    ),
    samples = samples, correction = correction,
    source = list(
      sample_manifest = sample_manifest,
      event_classifications = classifications,
      quantitation = list(
        ph3_metrics_acquisition = acquisition_metrics
      )
    )
  ), class = "ph3_output_model")
}

synthetic_ph3_background_spec <- function(sample_id, replicate_set_id,
                                          event_arguments = list()) {
  list(
    sample_id = sample_id, replicate_set_id = replicate_set_id,
    event_arguments = event_arguments
  )
}

test_that("SYNTHETIC pre-cutoff classifications receive explicit compatibility fields", {
  model <- synthetic_ph3_background_model(list(
    synthetic_ph3_background_spec("SYNTHETIC-sample", "SYNTHETIC-set")
  ))
  result <- apply_ph3_background_regression(model)
  events <- result$background_regression$event_signals
  expect_true(all(events$positivity_call_status == "called"))
  expect_true(all(is.na(events$positivity_call_reason_code)))
  expect_true(all(is.na(events$config_digest)))
  expect_true(all(is.na(events$export_operation_id)))
  expect_true(all(is.na(events$input_manifest_key)))
})

test_that("SYNTHETIC compatibility preserves legacy fit membership and rejects malformed new fields", {
  legacy <- synthetic_ph3_background_model(list(
    synthetic_ph3_background_spec(
      "SYNTHETIC-sample", "SYNTHETIC-set",
      list(nonfinite_negative = TRUE)
    )
  ))
  legacy_result <- apply_ph3_background_regression(legacy)
  expect_identical(legacy_result$correction$signal_basis, "raw")
  expect_identical(
    legacy_result$background_regression$individual_fits$validity_reason_code,
    "nonfinite_fit_input"
  )

  noncomputed <- synthetic_ph3_background_model(list(
    synthetic_ph3_background_spec("SYNTHETIC-sample", "SYNTHETIC-set")
  ))
  noncomputed$source$event_classifications[[1L]]$positivity_method_id <-
    rep("SYNTHETIC_flowjo_method", nrow(noncomputed$source$event_classifications[[1L]]))
  expect_silent(apply_ph3_background_regression(noncomputed))

  malformed <- synthetic_ph3_background_model(list(
    synthetic_ph3_background_spec("SYNTHETIC-sample", "SYNTHETIC-set")
  ))
  malformed$source$event_classifications[[1L]]$positivity_method_id <-
    rep(c("SYNTHETIC_flowjo_method", "ph3_raw_4n_density_cutoff_v1"),
        length.out = nrow(malformed$source$event_classifications[[1L]]))
  expect_error(apply_ph3_background_regression(malformed),
               "invalid_positivity_method")

  computed_missing <- synthetic_ph3_background_model(list(
    synthetic_ph3_background_spec("SYNTHETIC-sample", "SYNTHETIC-set")
  ))
  computed_missing$source$event_classifications[[1L]]$positivity_method_id <-
    rep("ph3_raw_4n_density_cutoff_v1", nrow(computed_missing$source$event_classifications[[1L]]))
  expect_error(apply_ph3_background_regression(computed_missing),
               "missing_computed_cutoff_event_field")

  invalid_computed_status <- synthetic_ph3_background_model(list(
    synthetic_ph3_background_spec("SYNTHETIC-sample", "SYNTHETIC-set")
  ))
  invalid_data <- invalid_computed_status$source$event_classifications[[1L]]
  invalid_data$positivity_method_id <- rep("ph3_raw_4n_density_cutoff_v1", nrow(invalid_data))
  invalid_data$config_digest <- rep("SYNTHETIC-config", nrow(invalid_data))
  invalid_data$export_operation_id <- rep("SYNTHETIC-operation", nrow(invalid_data))
  invalid_data$input_manifest_key <- rep("SYNTHETIC-manifest", nrow(invalid_data))
  invalid_data$positivity_call_status <- rep("not_a_valid_status", nrow(invalid_data))
  invalid_data$positivity_call_reason_code <- rep(NA_character_, nrow(invalid_data))
  invalid_computed_status$source$event_classifications[[1L]] <- invalid_data
  expect_error(apply_ph3_background_regression(invalid_computed_status),
               "invalid_computed_cutoff_call_state")

  malformed_field <- synthetic_ph3_background_model(list(
    synthetic_ph3_background_spec("SYNTHETIC-sample", "SYNTHETIC-set")
  ))
  malformed_field$source$event_classifications[[1L]]$positivity_call_status <-
    factor(rep("called", nrow(malformed_field$source$event_classifications[[1L]])))
  expect_error(apply_ph3_background_regression(malformed_field),
               "invalid_optional_event_field")
})

test_that("SYNTHETIC individual validity enforces the confirmed fit rules", {
  fit_99 <- ph3_fit_background_model(
    dna = c(seq(0, 10, length.out = 99), 5),
    signal = c(2 + 3 * seq(0, 10, length.out = 99), 17),
    negative_member = c(rep(TRUE, 99), FALSE),
    positive_member = c(rep(FALSE, 99), TRUE)
  )
  expect_identical(fit_99$validity_reason_code,
                   "insufficient_negative_events")

  fit_100 <- ph3_fit_background_model(
    dna = c(seq(0, 10, length.out = 100), 0, 10),
    signal = c(2 + 3 * seq(0, 10, length.out = 100), 12, 42),
    negative_member = c(rep(TRUE, 100), FALSE, FALSE),
    positive_member = c(rep(FALSE, 100), TRUE, TRUE)
  )
  expect_identical(fit_100$fit_status, "valid")
  expect_equal(fit_100$intercept, 2, tolerance = 1e-10)
  expect_equal(fit_100$slope, 3, tolerance = 1e-10)
  expect_identical(fit_100$coverage_status, "spans_positive_range")

  nonfinite <- ph3_fit_background_model(
    dna = seq(0, 10, length.out = 100),
    signal = c(Inf, 2 + 3 * seq(0, 10, length.out = 99)),
    negative_member = rep(TRUE, 100), positive_member = rep(FALSE, 100)
  )
  expect_identical(nonfinite$validity_reason_code, "nonfinite_fit_input")

  nonfinite_dna <- ph3_fit_background_model(
    dna = c(Inf, seq(0, 10, length.out = 99)),
    signal = 2 + 3 * c(0, seq(0, 10, length.out = 99)),
    negative_member = rep(TRUE, 100), positive_member = rep(FALSE, 100)
  )
  expect_identical(nonfinite_dna$validity_reason_code, "nonfinite_fit_input")

  zero_variation <- ph3_fit_background_model(
    dna = rep(5, 100), signal = seq_len(100),
    negative_member = rep(TRUE, 100), positive_member = rep(FALSE, 100)
  )
  expect_identical(zero_variation$validity_reason_code,
                   "zero_dna_variation")
  # With an intercept and one finite DNA predictor, exact rank deficiency
  # implies zero DNA variation. The ordered validator reports that more
  # specific cause before retaining the separate numerical-rank guard.

  nonfinite_prediction_input <- ph3_fit_background_model(
    dna = c(seq(0, 10, length.out = 100), Inf),
    signal = c(2 + 3 * seq(0, 10, length.out = 100), 12),
    negative_member = c(rep(TRUE, 100), FALSE),
    positive_member = c(rep(FALSE, 100), TRUE)
  )
  expect_identical(nonfinite_prediction_input$validity_reason_code,
                   "nonfinite_prediction_input")

  no_positive_events <- ph3_fit_background_model(
    dna = seq(0, 10, length.out = 100),
    signal = 2 + 3 * seq(0, 10, length.out = 100),
    negative_member = rep(TRUE, 100), positive_member = rep(FALSE, 100)
  )
  expect_identical(no_positive_events$fit_status, "valid")
  expect_identical(no_positive_events$coverage_status,
                   "not_applicable_no_positive_events")

  nonfinite_coefficients <- ph3_validate_background_solution(
    coefficients = c(Inf, 1), prediction_dna = 1
  )
  expect_identical(nonfinite_coefficients$status, "invalid")
  expect_identical(nonfinite_coefficients$reason_code,
                   "nonfinite_coefficients")

  nonfinite_predictions <- ph3_validate_background_solution(
    coefficients = rep(.Machine$double.xmax, 2L), prediction_dna = 1
  )
  expect_identical(nonfinite_predictions$status, "invalid")
  expect_identical(nonfinite_predictions$reason_code,
                   "nonfinite_predictions")
})

test_that("SYNTHETIC coverage equality is valid and a coverage gap is invalid", {
  equality <- ph3_fit_background_model(
    dna = c(seq(0, 10, length.out = 100), 0, 10),
    signal = c(2 + 3 * seq(0, 10, length.out = 100), 12, 42),
    negative_member = c(rep(TRUE, 100), FALSE, FALSE),
    positive_member = c(rep(FALSE, 100), TRUE, TRUE)
  )
  expect_identical(equality$fit_status, "valid")

  gap <- ph3_fit_background_model(
    dna = c(seq(1, 9, length.out = 100), 0, 10),
    signal = c(2 + 3 * seq(1, 9, length.out = 100), 12, 42),
    negative_member = c(rep(TRUE, 100), FALSE, FALSE),
    positive_member = c(rep(FALSE, 100), TRUE, TRUE)
  )
  expect_identical(gap$validity_reason_code,
                   "positive_dna_extrapolation_required")
})

test_that("SYNTHETIC one failure selects pooled correction only within its replicate set", {
  model <- synthetic_ph3_background_model(list(
    synthetic_ph3_background_spec("SYNTHETIC-s1", "SYNTHETIC-set-1"),
    synthetic_ph3_background_spec(
      "SYNTHETIC-s2", "SYNTHETIC-set-1", list(negative_n = 99L)
    ),
    synthetic_ph3_background_spec("SYNTHETIC-s3", "SYNTHETIC-set-2"),
    synthetic_ph3_background_spec("SYNTHETIC-s4", "SYNTHETIC-set-2")
  ))
  membership_before <- serialize(lapply(
    model$source$event_classifications,
    function(x) x[c(
      "event_identity", "ph3_positive_member", "eligible_2to4n",
      "sub_4n_member", "four_n_member"
    )]
  ), NULL, version = 3L)
  ab_before <- serialize(model$source$quantitation, NULL, version = 3L)
  result <- apply_ph3_background_regression(model)
  decision <- result$background_regression$set_decisions
  expect_identical(decision$signal_basis,
                   c("pooled_corrected", "individual_corrected"))
  set_one <- result$background_regression$event_signals$replicate_set_id ==
    "SYNTHETIC-set-1"
  expect_true(all(result$background_regression$event_signals$fit_scope[set_one] ==
                    "pooled"))
  expect_true(result$background_regression$pooled_fits$selected_for_analysis[[1L]])
  expect_identical(
    result$background_regression$individual_fits$selected_for_analysis,
    c(FALSE, FALSE, TRUE, TRUE)
  )
  expect_identical(
    result$background_regression$pooled_fits$selected_for_analysis,
    c(TRUE, FALSE)
  )
  expect_match(
    result$correction$reason_detail[[1L]],
    paste0(
      "individual_trigger_sample_id=SYNTHETIC-s2; ",
      "individual_trigger_reason_code=insufficient_negative_events"
    ),
    fixed = TRUE
  )
  expect_identical(decision$individual_trigger_sample_id,
                   c("SYNTHETIC-s2", NA_character_))
  expect_identical(decision$individual_trigger_reason_code,
                   c("insufficient_negative_events", NA_character_))
  expect_true(all(is.na(decision$pooled_failure_replicate_set_id)))
  expect_true(all(is.na(decision$pooled_failure_reason_code)))
  set_one_positive <- set_one &
    result$background_regression$event_signals$ph3_positive_member
  expect_equal(
    result$background_regression$event_signals$analytical_signal[
      set_one_positive
    ],
    rep(10, sum(set_one_positive)), tolerance = 1e-10
  )
  expect_identical(
    serialize(lapply(
      result$source$event_classifications,
      function(x) x[c(
        "event_identity", "ph3_positive_member", "eligible_2to4n",
        "sub_4n_member", "four_n_member"
      )]
    ), NULL, version = 3L), membership_before
  )
  expect_identical(serialize(result$source$quantitation, NULL, version = 3L),
                   ab_before)
})

test_that("SYNTHETIC pooled failure selects raw only within its replicate set", {
  model <- synthetic_ph3_background_model(list(
    synthetic_ph3_background_spec(
      "SYNTHETIC-s1", "SYNTHETIC-set-1",
      list(negative_dna = rep(5, 100))
    ),
    synthetic_ph3_background_spec(
      "SYNTHETIC-s2", "SYNTHETIC-set-1",
      list(negative_dna = rep(5, 100))
    ),
    synthetic_ph3_background_spec("SYNTHETIC-s3", "SYNTHETIC-set-2"),
    synthetic_ph3_background_spec("SYNTHETIC-s4", "SYNTHETIC-set-2")
  ))
  result <- apply_ph3_background_regression(model)
  expect_identical(result$correction$signal_basis,
                   c("raw", "individual_corrected"))
  expect_identical(result$correction$reason_code, c(
    "pooled_fit_failure_replicate_set_raw_fallback",
    "all_individual_fits_valid"
  ))
  expect_match(result$correction$reason_detail[[1L]],
               "pooled_failure_replicate_set_id=SYNTHETIC-set-1",
               fixed = TRUE)
  decision <- result$background_regression$set_decisions
  expect_identical(decision$pooled_failure_replicate_set_id,
                   c("SYNTHETIC-set-1", NA_character_))
  expect_identical(decision$pooled_failure_reason_code,
                   c("zero_dna_variation", NA_character_))
  expect_identical(
    result$background_regression$pooled_fits$fit_status,
    c("invalid", "not_attempted")
  )
  expect_identical(
    result$background_regression$pooled_fits$selected_for_analysis,
    c(FALSE, FALSE)
  )
  signals <- result$background_regression$event_signals
  set_one <- signals$replicate_set_id == "SYNTHETIC-set-1"
  expect_identical(signals$analytical_signal[set_one & signals$eligible_2to4n],
                   signals$target_raw[set_one & signals$eligible_2to4n])
  expect_true(all(signals$fit_scope[set_one] == "not_applicable"))
  expect_true(all(signals$fit_scope[!set_one] == "individual"))
})

test_that("SYNTHETIC all-valid experiment leaves pooled counts inapplicable", {
  model <- synthetic_ph3_background_model(list(
    synthetic_ph3_background_spec("SYNTHETIC-s1", "SYNTHETIC-set-1"),
    synthetic_ph3_background_spec("SYNTHETIC-s2", "SYNTHETIC-set-2")
  ))
  result <- apply_ph3_background_regression(model)
  pooled <- result$background_regression$pooled_fits
  expect_true(all(result$correction$signal_basis == "individual_corrected"))
  expect_true(all(pooled$fit_status == "not_attempted"))
  expect_type(pooled$negative_event_count, "integer")
  expect_type(pooled$positive_event_count, "integer")
  expect_true(all(is.na(pooled$negative_event_count)))
  expect_true(all(is.na(pooled$positive_event_count)))
})

test_that("SYNTHETIC pooled cross-sample coverage failure forces raw within that set", {
  model <- synthetic_ph3_background_model(list(
    synthetic_ph3_background_spec(
      "SYNTHETIC-s1", "SYNTHETIC-set-1",
      list(negative_dna = seq(2, 8, length.out = 50), positive_dna = 0)
    ),
    synthetic_ph3_background_spec(
      "SYNTHETIC-s2", "SYNTHETIC-set-1",
      list(negative_dna = seq(2, 8, length.out = 50), positive_dna = 10)
    ),
    synthetic_ph3_background_spec("SYNTHETIC-s3", "SYNTHETIC-set-2")
  ))
  result <- apply_ph3_background_regression(model)
  pooled <- result$background_regression$pooled_fits
  expect_identical(pooled$validity_reason_code[[1L]],
                   "positive_dna_extrapolation_required")
  expect_identical(result$correction$signal_basis,
                   c("raw", "individual_corrected"))
  expect_identical(result$background_regression$set_decisions$
                   pooled_failure_replicate_set_id,
                   c("SYNTHETIC-set-1", NA_character_))
})

test_that("SYNTHETIC empty positive compartments remain explicit", {
  empty_four_n <- synthetic_ph3_background_model(list(
    synthetic_ph3_background_spec(
      "SYNTHETIC-s1", "SYNTHETIC-set-1", list(positive_dna = c(1, 2))
    )
  ))
  four_result <- apply_ph3_background_regression(empty_four_n)
  four_events <- four_result$background_regression$event_signals
  expect_equal(sum(four_events$ph3_positive_member &
                     four_events$four_n_member), 0L)
  expect_equal(nrow(four_events), 102L)
  expect_identical(four_result$correction$signal_basis,
                   "individual_corrected")

  empty_below_four_n <- synthetic_ph3_background_model(list(
    synthetic_ph3_background_spec(
      "SYNTHETIC-s1", "SYNTHETIC-set-1", list(positive_dna = c(8, 9))
    )
  ))
  below_result <- apply_ph3_background_regression(empty_below_four_n)
  below_events <- below_result$background_regression$event_signals
  expect_equal(sum(below_events$ph3_positive_member &
                     below_events$sub_4n_member), 0L)
  expect_equal(nrow(below_events), 102L)
  expect_identical(below_result$correction$signal_basis,
                   "individual_corrected")
})

test_that("SYNTHETIC nonfinite positive signal is explicit and not omitted", {
  model <- synthetic_ph3_background_model(list(
    synthetic_ph3_background_spec(
      "SYNTHETIC-s1", "SYNTHETIC-set-1",
      list(nonfinite_positive_signal = TRUE)
    )
  ))
  result <- apply_ph3_background_regression(model)
  signals <- result$background_regression$event_signals
  expect_equal(nrow(signals), 102L)
  expect_identical(result$correction$signal_basis, "individual_corrected")
  expect_identical(
    signals$event_signal_status[!is.finite(signals$target_raw)], "unavailable"
  )
  expect_identical(
    signals$event_signal_reason_code[!is.finite(signals$target_raw)],
    "nonfinite_raw_or_corrected_signal"
  )
})

test_that("SYNTHETIC altered Slice 1 identity or membership fails closed", {
  model <- synthetic_ph3_background_model(list(
    synthetic_ph3_background_spec("SYNTHETIC-s1", "SYNTHETIC-set-1"),
    synthetic_ph3_background_spec("SYNTHETIC-s2", "SYNTHETIC-set-1")
  ))
  duplicate_sample <- model
  duplicate_sample$samples$sample_id[[2L]] <- "SYNTHETIC-s1"
  expect_error(
    apply_ph3_background_regression(duplicate_sample),
    "invalid_model_identity_mapping"
  )

  malformed_analysis <- model
  malformed_analysis$experiment$analysis_id <- "SYNTHETIC-analysis"
  expect_error(
    apply_ph3_background_regression(malformed_analysis),
    "invalid_model_identity_mapping"
  )

  deleted_acquisition <- model
  deleted_acquisition$source$event_classifications <-
    deleted_acquisition$source$event_classifications[-1L]
  expect_error(
    apply_ph3_background_regression(deleted_acquisition),
    "classification_acquisition_mismatch"
  )

  duplicated_acquisition <- model
  duplicated_acquisition$source$event_classifications <- c(
    duplicated_acquisition$source$event_classifications,
    duplicated_acquisition$source$event_classifications[1L]
  )
  expect_error(
    apply_ph3_background_regression(duplicated_acquisition),
    "classification_acquisition_mismatch"
  )

  remapped <- model
  remapped$source$sample_manifest$replicate_set_id[[1L]] <-
    "SYNTHETIC-remapped-set"
  expect_error(
    apply_ph3_background_regression(remapped), "invalid_source_mapping"
  )

  altered_membership <- model
  altered_membership$source$event_classifications[[1L]]$
    four_n_member[[1L]] <- TRUE
  altered_membership$source$event_classifications[[1L]]$
    sub_4n_member[[1L]] <- TRUE
  expect_error(
    apply_ph3_background_regression(altered_membership),
    "altered_event_membership"
  )

  altered_positivity <- model
  altered_positivity$source$event_classifications[[1L]]$
    ph3_positive_member[[1L]] <- TRUE
  expect_error(
    apply_ph3_background_regression(altered_positivity),
    "classification_metric_count_mismatch"
  )

  altered_eligibility <- model
  altered_eligibility$source$event_classifications[[1L]]$
    eligible_2to4n[[1L]] <- FALSE
  altered_eligibility$source$event_classifications[[1L]]$
    sub_4n_member[[1L]] <- FALSE
  altered_eligibility$source$event_classifications[[1L]]$
    four_n_member[[1L]] <- FALSE
  expect_error(
    apply_ph3_background_regression(altered_eligibility),
    "classification_metric_count_mismatch"
  )
})
