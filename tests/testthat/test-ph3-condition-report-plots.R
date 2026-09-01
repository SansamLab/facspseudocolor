# Every constructed value in this file is explicitly SYNTHETIC and test-only.

synthetic_ph3_plot_model <- function(reference = TRUE, mixed_basis = FALSE) {
  columns <- ph3_report_model_columns()
  conditions <- data.frame(
    condition_id = c("SYNTHETIC-control", "SYNTHETIC-treatment"),
    condition_label = c("SYNTHETIC Control", "SYNTHETIC Treatment"),
    condition_index = 1:2, stringsAsFactors = FALSE
  )
  outcomes <- data.frame(
    outcome_id = c("A", "B", "C", "D"),
    label = c(
      "4N pH3-positive prevalence", "Below-4N pH3-positive prevalence",
      "pH3 signal in 4N pH3-positive cells",
      "pH3 signal in below-4N pH3-positive cells"
    ), stringsAsFactors = FALSE
  )
  values <- do.call(rbind, lapply(seq_len(nrow(conditions)), function(condition_i) {
    do.call(rbind, lapply(1:2, function(replicate_i) {
      do.call(rbind, lapply(c("A", "B", "C", "D"), function(outcome_id) {
        signal <- outcome_id %in% c("C", "D")
        basis <- if (!signal) "not_applicable" else if (mixed_basis &&
          replicate_i == 2L) "raw" else "individual_corrected"
        data.frame(
          experiment_id = "SYNTHETIC-experiment", analysis_id = "SYNTHETIC-analysis",
          condition_id = conditions$condition_id[[condition_i]],
          condition_label = conditions$condition_label[[condition_i]],
          condition_index = as.integer(condition_i),
          replicate_set_id = paste0("SYNTHETIC-set-", replicate_i),
          replicate_label = paste("SYNTHETIC R", replicate_i),
          replicate_index = as.integer(replicate_i), outcome_id = outcome_id,
          outcome_label = outcomes$label[outcomes$outcome_id == outcome_id],
          value_kind = if (signal) "population_signal_median" else "prevalence_percent",
          value_mode = if (!signal) "prevalence_percent" else if (reference) {
            "reference_ratio"
          } else "direct_median",
          signal_basis = basis, value = as.numeric(condition_i * 10 + replicate_i),
          value_status = "available", reason_code = NA_character_,
          technical_acquisition_count = 1L, valid_acquisition_count = 1L,
          source_acquisition_ids = paste0("[\"SYNTHETIC-", condition_i, "-", replicate_i, "\"]"),
          stringsAsFactors = FALSE
        )
      }))
    }))
  }))
  values <- values[columns$biological_replicate_values]
  summary <- do.call(rbind, lapply(seq_len(nrow(conditions)), function(condition_i) {
    do.call(rbind, lapply(c("A", "B", "C", "D"), function(outcome_id) {
      signal <- outcome_id %in% c("C", "D")
      mixed <- signal && mixed_basis
      data.frame(
        experiment_id = "SYNTHETIC-experiment", analysis_id = "SYNTHETIC-analysis",
        condition_id = conditions$condition_id[[condition_i]],
        condition_label = conditions$condition_label[[condition_i]],
        condition_index = as.integer(condition_i), outcome_id = outcome_id,
        outcome_label = outcomes$label[outcomes$outcome_id == outcome_id],
        value_kind = if (signal) "population_signal_median" else "prevalence_percent",
        value_mode = if (!signal) "prevalence_percent" else if (reference) {
          "reference_ratio"
        } else "direct_median",
        signal_basis = if (!signal) "not_applicable" else if (mixed) "mixed" else "individual_corrected",
        mean_value = if (mixed) NA_real_ else condition_i * 10 + 1.5,
        sd_value = NA_real_, biological_replicate_count = 2L,
        finite_biological_replicate_count = 2L,
        undefined_biological_replicate_count = 0L,
        summary_status = if (mixed) "unavailable_mixed_signal_basis" else "available",
        reason_code = if (mixed) "mixed_signal_basis" else NA_character_,
        source_replicate_set_ids = "[\"SYNTHETIC-set-1\",\"SYNTHETIC-set-2\"]",
        stringsAsFactors = FALSE
      )
    }))
  }))
  summary <- summary[columns$condition_summary]
  qc <- summary[c("experiment_id", "analysis_id", "condition_id", "outcome_id", "summary_status")]
  qc$signal_basis_status <- ifelse(qc$outcome_id %in% c("C", "D"), "single_or_unavailable", "not_applicable")
  qc$reference_mode <- summary$value_mode
  qc$reason_code <- summary$reason_code
  qc <- qc[columns$qc_flags]
  structure(list(
    experiment = data.frame(experiment_id = "SYNTHETIC-experiment",
                            analysis_id = "SYNTHETIC-analysis", stringsAsFactors = FALSE),
    outcomes = outcomes, conditions = conditions,
    condition_report = list(
      schema_version = "ph3-condition-report-model-1.0.0",
      biological_replicate_values = values, condition_summary = summary,
      signal_basis_strata = ph3_report_model_empty(columns$signal_basis_strata),
      qc_flags = qc,
      provenance = list(analysis_id = "SYNTHETIC-analysis",
                        experiment_id = "SYNTHETIC-experiment")
    )
  ), class = "ph3_output_model")
}

test_that("SYNTHETIC Slice 5 panels retain exactly the owner-confirmed A-D model", {
  model <- synthetic_ph3_plot_model()
  report <- ph3_report_plot_validate_model(model)
  colours <- ph3_report_plot_colours(model$conditions)
  panels <- stats::setNames(lapply(ph3_report_plot_outcome_ids(), function(id) {
    ph3_report_plot_panel(report, model$outcomes, id, colours)
  }), ph3_report_plot_outcome_ids())
  expect_identical(names(panels), c("A", "B", "C", "D"))
  expect_true(all(vapply(panels, inherits, logical(1), "ggplot")))
  expect_match(panels$A$labels$y, "Analysis singlets")
  expect_match(panels$B$labels$y, "Analysis singlets")
  expect_match(panels$C$labels$y, "matched reference")
  expect_match(panels$D$labels$y, "matched reference")
})

test_that("SYNTHETIC Slice 5 signal panels state direct and unavailable bases", {
  direct <- synthetic_ph3_plot_model(reference = FALSE)
  direct_panel <- ph3_report_plot_panel(
    ph3_report_plot_validate_model(direct), direct$outcomes, "C",
    ph3_report_plot_colours(direct$conditions)
  )
  expect_match(direct_panel$labels$y, "background-corrected")

  mixed <- synthetic_ph3_plot_model(mixed_basis = TRUE)
  mixed_panel <- ph3_report_plot_panel(
    ph3_report_plot_validate_model(mixed), mixed$outcomes, "C",
    ph3_report_plot_colours(mixed$conditions)
  )
  expect_match(mixed_panel$labels$caption, "Unavailable condition summary")
  expect_match(mixed_panel$labels$caption, "mixed_signal_basis")
  expect_match(mixed_panel$labels$caption, "RAW—NOT BACKGROUND SUBTRACTED")
  expect_match(mixed_panel$labels$y, "RAW—NOT BACKGROUND SUBTRACTED")
})

test_that("SYNTHETIC Slice 5 rejects a tampered condition report", {
  model <- synthetic_ph3_plot_model()
  model$condition_report$provenance$analysis_id <- "SYNTHETIC-tampered"
  expect_error(ph3_report_plot_validate_model(model), "invalid_condition_report_model")
})
