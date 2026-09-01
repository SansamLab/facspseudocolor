# Every constructed value in this file is explicitly SYNTHETIC and test-only.

synthetic_ph3_condition_report_model <- function(reference = TRUE) {
  conditions <- data.frame(
    condition_id = c("SYNTHETIC-control", "SYNTHETIC-treatment"),
    condition_label = c("SYNTHETIC Control", "SYNTHETIC Treatment"),
    role = c("control", "treatment"), condition_index = 1:2,
    stringsAsFactors = FALSE
  )
  sets <- data.frame(
    experiment_id = "SYNTHETIC-experiment",
    replicate_set_id = c("SYNTHETIC-set-1", "SYNTHETIC-set-2"),
    replicate_label = c("SYNTHETIC R1", "SYNTHETIC R2"),
    manifest_replicate_index = 1:2, stringsAsFactors = FALSE
  )
  samples <- do.call(rbind, lapply(seq_len(nrow(sets)), function(i) {
    data.frame(
      experiment_id = "SYNTHETIC-experiment",
      replicate_set_id = sets$replicate_set_id[[i]],
      sample_id = paste0("SYNTHETIC-s", i, "-", c("control", "treatment")),
      condition_id = conditions$condition_id, acquisition_count = 1L,
      source_acquisition_ids = paste0("[\"SYNTHETIC-a", i, "-", c("c", "t"), "\"]"),
      stringsAsFactors = FALSE
    )
  }))
  outcomes <- data.frame(
    outcome_id = c("A", "B", "C", "D"),
    label = c(
      "4N pH3-positive prevalence", "Below-4N pH3-positive prevalence",
      "pH3 signal in 4N pH3-positive cells",
      "pH3 signal in below-4N pH3-positive cells"
    ),
    value_kind = c("prevalence_percent", "prevalence_percent",
                   "population_signal_median", "population_signal_median"),
    stringsAsFactors = FALSE
  )
  prevalence <- do.call(rbind, lapply(c("A", "B"), function(outcome_id) {
    do.call(rbind, lapply(seq_len(nrow(sets)), function(i) {
      data.frame(
        analysis_id = "SYNTHETIC-analysis",
        condition = conditions$condition_label,
        condition_index = conditions$condition_index,
        replicate = sets$replicate_label[[i]], replicate_index = i,
        metric_id = if (outcome_id == "A") {
          "ph3_4n_positive_prevalence_within_2to4n_percent"
        } else "ph3_sub_4n_positive_prevalence_within_2to4n_percent",
        value_percent = c(10, 20) + i,
        result_status = "ok", technical_acquisition_count = 1L,
        finite_technical_acquisition_count = 1L,
        source_acquisition_ids = samples$source_acquisition_ids[(i - 1L) * 2L + 1:2],
        stringsAsFactors = FALSE
      )
    }))
  }))
  correction <- data.frame(
    experiment_id = "SYNTHETIC-experiment", replicate_set_id = sets$replicate_set_id,
    status = "selected", signal_basis = "individual_corrected",
    reason_code = "all_individual_fits_valid", reason_detail = NA_character_,
    stringsAsFactors = FALSE
  )
  reference_table <- data.frame(
    experiment_id = "SYNTHETIC-experiment", replicate_set_id = sets$replicate_set_id,
    status = if (reference) "configured" else "not_configured",
    condition_id = if (reference) "SYNTHETIC-control" else NA_character_,
    sample_id = if (reference) samples$sample_id[c(1L, 3L)] else NA_character_,
    reason_code = if (reference) NA_character_ else "reference_not_configured",
    reason_detail = NA_character_, stringsAsFactors = FALSE
  )
  signal <- do.call(rbind, lapply(seq_len(nrow(samples)), function(i) {
    sample <- samples[i, , drop = FALSE]
    do.call(rbind, lapply(c("C", "D"), function(outcome_id) {
      direct <- if (outcome_id == "C") 100 + i else 50 + i
      data.frame(
        experiment_id = "SYNTHETIC-experiment",
        replicate_set_id = sample$replicate_set_id, sample_id = sample$sample_id,
        condition_id = sample$condition_id, outcome_id = outcome_id,
        signal_basis = "individual_corrected", acquisition_count = 1L,
        valid_acquisition_count = 1L, direct_value = direct,
        direct_status = "available", direct_reason_code = NA_character_,
        reference_ratio = if (reference) ifelse(sample$condition_id == "SYNTHETIC-control", 1, 2) else NA_real_,
        reference_status = if (reference) "available" else "not_applicable",
        reference_reason_code = if (reference) NA_character_ else "reference_not_configured",
        stringsAsFactors = FALSE
      )
    }))
  }))
  structure(list(
    schema = list(schema_version = "ph3-output-contract-model-1.0.0"),
    experiment = data.frame(experiment_id = "SYNTHETIC-experiment",
                            analysis_id = "SYNTHETIC-analysis", stringsAsFactors = FALSE),
    outcomes = outcomes, conditions = conditions, replicate_sets = sets,
    samples = samples, correction = correction, reference = reference_table,
    source = list(quantitation = list(ph3_metrics_biological_replicate = prevalence)),
    signal_outcomes = list(
      schema_version = "ph3-signal-outcomes-1.0.0", acquisition = data.frame(),
      sample = signal, basis_qc = correction, reference_qc = reference_table
    )
  ), class = "ph3_output_model")
}

test_that("SYNTHETIC report model retains the four owner-confirmed outcomes", {
  result <- derive_ph3_condition_report_model(
    synthetic_ph3_condition_report_model()
  )$condition_report
  expect_identical(result$schema_version, "ph3-condition-report-model-1.0.0")
  expect_identical(names(result), c(
    "schema_version", "biological_replicate_values", "condition_summary",
    "signal_basis_strata", "qc_flags", "provenance"
  ))
  expect_identical(result$biological_replicate_values$outcome_id, rep(
    c("A", "B", "C", "D"), times = 4L
  ))
  expect_identical(nrow(result$biological_replicate_values), 16L)
  expect_identical(nrow(result$condition_summary), 8L)
  expect_identical(result$condition_summary$value_mode[
    result$condition_summary$outcome_id %in% c("A", "B")
  ], rep("prevalence_percent", 4L))
  expect_identical(result$condition_summary$value_mode[
    result$condition_summary$outcome_id %in% c("C", "D")
  ], rep("reference_ratio", 4L))
  expect_true(all(result$condition_summary$biological_replicate_count == 2L))
  expect_true(all(is.na(result$condition_summary$sd_value)))
  expect_identical(
    result$provenance$source_signal_outcomes_schema_version,
    "ph3-signal-outcomes-1.0.0"
  )
})

test_that("SYNTHETIC signal report values use direct medians only without reference", {
  result <- derive_ph3_condition_report_model(
    synthetic_ph3_condition_report_model(reference = FALSE)
  )$condition_report
  signal <- result$biological_replicate_values[
    result$biological_replicate_values$outcome_id %in% c("C", "D"), , drop = FALSE
  ]
  expect_true(all(signal$value_mode == "direct_median"))
  expect_true(all(signal$value_status %in% c(
    "available", "available_partial_acquisition_coverage", "unavailable"
  )))
})

test_that("SYNTHETIC partial prevalence and mixed signal bases remain explicit", {
  partial <- synthetic_ph3_condition_report_model()
  partial$source$quantitation$ph3_metrics_biological_replicate$result_status[[1L]] <-
    "ok_partial_undefined"
  partial_result <- derive_ph3_condition_report_model(partial)$condition_report
  expect_identical(
    partial_result$biological_replicate_values$value_status[[1L]],
    "available_partial_acquisition_coverage"
  )
  expect_identical(
    partial_result$biological_replicate_values$reason_code[[1L]],
    "one_or_more_technical_acquisition_values_undefined"
  )

  mixed <- synthetic_ph3_condition_report_model()
  set_two <- "SYNTHETIC-set-2"
  mixed$correction$signal_basis[
    mixed$correction$replicate_set_id == set_two
  ] <- "raw"
  mixed$signal_outcomes$basis_qc <- mixed$correction
  mixed$signal_outcomes$sample$signal_basis[
    mixed$signal_outcomes$sample$replicate_set_id == set_two
  ] <- "raw"
  mixed_result <- derive_ph3_condition_report_model(mixed)$condition_report
  combined <- mixed_result$condition_summary[
    mixed_result$condition_summary$outcome_id %in% c("C", "D"), , drop = FALSE
  ]
  expect_true(all(combined$summary_status ==
                    "unavailable_mixed_signal_basis"))
  expect_true(all(is.na(combined$mean_value)))
  expect_identical(nrow(mixed_result$signal_basis_strata), 8L)
})

test_that("SYNTHETIC report model rejects tampered basis and prevalence provenance", {
  mixed <- synthetic_ph3_condition_report_model()
  mixed$signal_outcomes$sample$signal_basis[[1L]] <- "raw"
  expect_error(
    derive_ph3_condition_report_model(mixed),
    "signal_basis_provenance_mismatch"
  )

  remapped <- synthetic_ph3_condition_report_model()
  remapped$source$quantitation$ph3_metrics_biological_replicate$condition[[1L]] <-
    "SYNTHETIC wrong"
  expect_error(
    derive_ph3_condition_report_model(remapped),
    "prevalence_provenance_mismatch"
  )
})
