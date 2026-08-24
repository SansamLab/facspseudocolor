# All constructed events in this file are SYNTHETIC and test-only. They encode
# arithmetic and boundary invariants, not biological observations or claims.

synthetic_edu_contract_analysis <- function(display_offset = 10000) {
  config <- validate_facs_config(minimal_config("edu"))
  config$show_reference_panel <- TRUE
  config$y_limits <- c(100, 30000)
  config$background_subtracted_offset <- display_offset
  manifest <- build_sample_manifest(config)

  dna <- c(
    800, 1000, 1224, 1225, 1300, 1619, 1620, 1674.5, 1675, 2000,
    800, rep(1000, 10), rep(1400, 10), rep(1800, 10), 2000, 2200, 2200
  )
  positive <- c(
    rep(FALSE, 10), TRUE, rep(TRUE, 30), TRUE, FALSE, TRUE
  )
  make_sample <- function(signal_offset) {
    target_bgsub <- seq_along(dna) * 10 + signal_offset
    list(
      data = data.frame(
        dna_norm = dna,
        target_norm = target_bgsub + 10000,
        target_bgsub = target_bgsub,
        edu_computed_positive = positive
      ),
      edu_positive = NULL,
      normalization_method = "reference_negative_regression",
      g1_anchor_target = 6000 + signal_offset
    )
  }
  normalized <- list(reference = make_sample(0), treatment = make_sample(100))
  new_facs_analysis(
    config, manifest,
    data.frame(path = character(), exists = logical()),
    normalized, models = list(list())
  )
}

quantify_synthetic_edu_contract <- function(display_offset = 10000) {
  expect_warning(
    result <- quantify_cell_cycle(
      synthetic_edu_contract_analysis(display_offset)
    ),
    "aliases.*deprecated"
  )
  result
}

test_that("EdU contract exposes only the seven approved canonical tables", {
  result <- quantify_synthetic_edu_contract()
  expected <- c(
    "edu_assigned_phase_composition",
    "edu_single_cells_phase_composition",
    "edu_six_gate_phase_composition",
    "edu_regional_positivity",
    "edu_overall_positivity",
    "edu_positive_cell_regional_intensity",
    "edu_positive_population_intensity"
  )
  expect_true(all(expected %in% names(result$quantitation)))
  expect_true(all(paste0(expected, "_acquisition") %in%
                    names(result$quantitation)))
  expect_false(any(grepl(
    "all_cell_regional_intensity|single_cells_intensity|figure1",
    names(result$quantitation), ignore.case = TRUE
  )))
  expect_identical(result$provenance$output_schema_version, 2L)
  expect_true(all(vapply(
    result$quantitation[expected],
    function(x) all(x$output_schema_version == 2L), logical(1)
  )))
})

test_that("Negative S uses approved half-open boundaries", {
  result <- quantify_synthetic_edu_contract()
  classified <- result$normalized_data$reference$edu_event_classification

  expect_identical(
    classified$historical_five_gate_assignment[c(3, 9)],
    c("g1", "g2m")
  )
  expect_identical(classified$negative_s[c(3, 4, 8, 9)],
                   c(FALSE, TRUE, TRUE, FALSE))
  expect_identical(classified$six_gate_assignment[c(4, 8)],
                   c("negative_s", "negative_s"))
})

test_that("composition denominators and Unassigned QC reconcile exactly", {
  result <- quantify_synthetic_edu_contract()
  q <- result$quantitation
  assigned <- q$edu_assigned_phase_composition_acquisition
  single <- q$edu_single_cells_phase_composition_acquisition
  six <- q$edu_six_gate_phase_composition_acquisition
  first <- assigned$condition == "Reference"

  expect_equal(assigned$numerator_n[first], c(3, 10, 10, 10, 2))
  expect_equal(assigned$denominator_n[first], rep(35, 5))
  expect_equal(assigned$unassigned_n[first], rep(9, 5))
  expect_equal(
    assigned$assigned_phase_composition_pct[first],
    100 * c(3, 10, 10, 10, 2) / 35
  )
  expect_equal(single$denominator_n[single$condition == "Reference"],
               rep(44, 6))
  expect_equal(single$unassigned_n[single$condition == "Reference"],
               rep(4, 6))
  expect_equal(
    unique(single$unassigned_pct_of_eligible_single_cells[
      single$condition == "Reference"
    ]),
    100 * 4 / 44
  )
  expect_equal(sum(single$single_cells_phase_composition_pct[
    single$condition == "Reference"
  ]), 100 * 40 / 44)
  expect_equal(sum(six$six_gate_phase_composition_pct[
    six$condition == "Reference"
  ]), 100)
})

test_that("regional and overall positivity retain every eligible denominator event", {
  result <- quantify_synthetic_edu_contract()
  regional <- result$quantitation$edu_regional_positivity_acquisition
  regional <- regional[regional$condition == "Reference", ]
  overall <- result$quantitation$edu_overall_positivity_acquisition
  overall <- overall[overall$condition == "Reference", ]

  expect_equal(regional$numerator_n, c(10, 10, 10))
  expect_equal(regional$denominator_n, c(13, 12, 13))
  expect_equal(regional$regional_edu_positive_pct,
               100 * c(10 / 13, 10 / 12, 10 / 13))
  expect_equal(overall$numerator_n, 32)
  expect_equal(overall$denominator_n, 42)
  expect_equal(overall$phase_unassigned_n, 2)
  expect_equal(overall$overall_edu_positive_pct, 100 * 32 / 42)
  expect_equal(overall$dna_interval_min, 775)
  expect_equal(overall$dna_interval_max, 2125)
  expect_identical(overall$dna_interval_lower_inclusive, TRUE)
  expect_identical(overall$dna_interval_upper_inclusive, FALSE)
})

test_that("regional and overall intervals are exactly half-open", {
  analysis <- synthetic_edu_contract_analysis()
  config <- analysis$config
  boundary_data <- data.frame(
    dna_norm = c(774.999, 775, 1264.999, 1265, 1619.999, 1620,
                 2124.999, 2125),
    target_norm = rep(10000, 8),
    target_bgsub = seq_len(8),
    edu_computed_positive = rep(TRUE, 8)
  )
  classified <- facspseudocolor:::classify_edu_events(
    boundary_data, config, 10000
  )
  manifest <- analysis$sample_manifest[1, , drop = FALSE]
  regional <- facspseudocolor:::collect_edu_regional_positivity(
    list(classified), manifest, config, 10000
  )
  overall <- facspseudocolor:::collect_edu_overall_positivity(
    list(classified), manifest, config, 10000
  )

  expect_equal(regional$denominator_n, c(1, 2, 1))
  expect_equal(regional$numerator_n, c(1, 2, 1))
  expect_equal(overall$denominator_n, 6)
  expect_equal(overall$numerator_n, 6)
})

test_that("zero regional denominators are structured and do not invent values", {
  analysis <- synthetic_edu_contract_analysis()
  for (i in seq_along(analysis$normalized_data)) {
    keep <- analysis$normalized_data[[i]]$data$dna_norm < 1265
    analysis$normalized_data[[i]]$data <-
      analysis$normalized_data[[i]]$data[keep, , drop = FALSE]
  }
  expect_warning(result <- quantify_cell_cycle(analysis), "deprecated")
  regional <- result$quantitation$edu_regional_positivity_acquisition
  empty <- regional$region %in% c("mid", "late")
  expect_true(all(is.na(regional$regional_edu_positive_pct[empty])))
  expect_true(all(regional$metric_status[empty] == "zero_denominator"))
  intensity <- result$quantitation$
    edu_positive_cell_regional_intensity_acquisition
  insufficient <- intensity$phase %in% c("mid", "late")
  expect_true(all(is.na(
    intensity$positive_cell_regional_edu_bgsub_median[insufficient]
  )))
  expect_true(all(intensity$metric_status[insufficient] ==
                    "insufficient_events"))
})

test_that("canonical intensities preserve positive-only legacy values", {
  result <- quantify_synthetic_edu_contract()
  q <- result$quantitation
  regional <- q$edu_positive_cell_regional_intensity_acquisition
  regional <- regional[regional$condition == "Reference", ]
  whole <- q$edu_positive_population_intensity_acquisition
  whole <- whole[whole$condition == "Reference", ]

  expect_equal(regional$positive_cell_regional_edu_bgsub_median,
               c(165, 265, 365))
  expect_equal(whole$positive_population_edu_bgsub_median, 270)
  expect_equal(
    q$phase_medians_acquisition$median_signal,
    q$edu_positive_cell_regional_intensity_acquisition$
      positive_cell_regional_edu_bgsub_median
  )
  expect_equal(
    q$whole_medians_acquisition$median_signal,
    q$edu_positive_population_intensity_acquisition$
      positive_population_edu_bgsub_median
  )
  expect_equal(
    q$phase_percentages_acquisition$phase_percent,
    q$edu_assigned_phase_composition_acquisition$
      assigned_phase_composition_pct
  )
  expect_equal(
    q$phase_medians$median_signal,
    q$edu_positive_cell_regional_intensity$
      positive_cell_regional_edu_bgsub_median
  )
  expect_equal(
    q$whole_medians$median_signal,
    q$edu_positive_population_intensity$
      positive_population_edu_bgsub_median
  )
  expect_equal(
    q$phase_percentages$phase_percent,
    q$edu_assigned_phase_composition$assigned_phase_composition_pct
  )
})

test_that("display offset is recorded but cannot change quantitative values", {
  first <- quantify_synthetic_edu_contract(10000)
  second <- quantify_synthetic_edu_contract(20000)
  first_q <- first$quantitation$edu_positive_population_intensity_acquisition
  second_q <- second$quantitation$edu_positive_population_intensity_acquisition

  expect_equal(first_q$positive_population_edu_bgsub_median,
               second_q$positive_population_edu_bgsub_median)
  first_events <- first$normalized_data$reference$edu_event_classification
  second_events <- second$normalized_data$reference$edu_event_classification
  expect_equal(second_events$display_signal - first_events$display_signal,
               rep(10000, nrow(first_events)))
  expect_equal(first_events$quantitative_signal,
               second_events$quantitative_signal)
  expect_false(first$provenance$edu_output_contract$
                 quantitative_values_include_display_offset)
})

test_that("legacy divided display records its actual coordinate and zero offset", {
  analysis <- synthetic_edu_contract_analysis(10000)
  analysis$config$pseudocolor_signal <- "normalized"
  for (i in seq_along(analysis$normalized_data)) {
    analysis$normalized_data[[i]]$g1_anchor_target <- NA_real_
  }
  expect_warning(result <- quantify_cell_cycle(analysis), "deprecated")
  classified <- result$normalized_data$reference$edu_event_classification
  canonical <- result$quantitation$edu_positive_population_intensity_acquisition

  expect_equal(classified$display_signal,
               analysis$normalized_data$reference$data$target_norm)
  expect_true(all(classified$display_transform ==
                    "legacy_background_divided"))
  expect_true(all(classified$display_offset == 0))
  expect_false(any(classified$display_offset_applied))
  expect_true(all(canonical$display_transform ==
                    "legacy_background_divided"))
  expect_true(all(canonical$display_offset == 0))
  expect_identical(result$provenance$edu_output_contract$display_offset, 0)
})

test_that("deprecated aliases warn once and plot switches cannot change calculations", {
  analysis <- synthetic_edu_contract_analysis()
  analysis$config$quantify_phase_median <- FALSE
  analysis$config$quantify_whole_median <- FALSE
  analysis$config$quantify_phase_percent <- FALSE
  expect_warning(first <- quantify_cell_cycle(analysis), "deprecated")
  expect_warning(second <- quantify_cell_cycle(first), NA)
  expect_equal(sum(grepl("aliases.*deprecated", second$warnings)), 1)

  enabled <- synthetic_edu_contract_analysis()
  enabled$config$quantify_phase_median <- TRUE
  enabled$config$quantify_whole_median <- TRUE
  enabled$config$quantify_phase_percent <- TRUE
  expect_warning(enabled <- quantify_cell_cycle(enabled), "deprecated")
  canonical <- c(
    "edu_assigned_phase_composition", "edu_single_cells_phase_composition",
    "edu_six_gate_phase_composition", "edu_regional_positivity",
    "edu_overall_positivity", "edu_positive_cell_regional_intensity",
    "edu_positive_population_intensity"
  )
  expect_equal(first$quantitation[canonical], enabled$quantitation[canonical])
})

test_that("technical acquisitions retain unweighted summary aggregation", {
  acquisition <- data.frame(
    replicate = c("SYNTHETIC Bio 1", "SYNTHETIC Bio 1"),
    replicate_index = c(1L, 1L),
    technical_replicate = c("SYNTHETIC T1", "SYNTHETIC T2"),
    condition = c("SYNTHETIC condition", "SYNTHETIC condition"),
    condition_index = c(1L, 1L),
    region = c("early", "early"), region_label = c("Early S", "Early S"),
    region_index = c(1L, 1L), numerator_n = c(1L, 90L),
    denominator_n = c(10L, 100L), regional_edu_positive_pct = c(10, 90),
    source_population = "eligible_single_cells_in_dna_region",
    dna_interval_min = 910, dna_interval_max = 1265,
    dna_interval_lower_inclusive = TRUE,
    dna_interval_upper_inclusive = FALSE,
    signal_transform = "background_subtracted",
    display_transform = "background_subtracted_plus_offset",
    display_offset = 10000, display_offset_applied = TRUE,
    reference_normalization_status = "not_applied",
    output_schema_version = 2L, metric_status = "ok"
  )
  averaged <- facspseudocolor:::average_edu_metric_table(
    acquisition, "regional_edu_positive_pct",
    c("region", "region_label", "region_index")
  )
  expect_equal(averaged$regional_edu_positive_pct, 50)
  expect_equal(averaged$numerator_n, 91)
  expect_equal(averaged$denominator_n, 110)
  expect_equal(averaged$technical_n, 2)
})

test_that("technical aggregation preserves mixed acquisition status", {
  acquisition <- data.frame(
    replicate = c("SYNTHETIC Bio 1", "SYNTHETIC Bio 1"),
    replicate_index = c(1L, 1L),
    technical_replicate = c("SYNTHETIC T1", "SYNTHETIC T2"),
    condition = c("SYNTHETIC condition", "SYNTHETIC condition"),
    condition_index = c(1L, 1L),
    region = c("early", "early"), region_label = c("Early S", "Early S"),
    region_index = c(1L, 1L), numerator_n = c(1L, 0L),
    denominator_n = c(10L, 0L), regional_edu_positive_pct = c(10, NA),
    source_population = "eligible_single_cells_in_dna_region",
    dna_interval_min = 910, dna_interval_max = 1265,
    dna_interval_lower_inclusive = TRUE,
    dna_interval_upper_inclusive = FALSE,
    signal_transform = "background_subtracted",
    display_transform = "background_subtracted_plus_offset",
    display_offset = 10000, display_offset_applied = TRUE,
    reference_normalization_status = "not_applied",
    output_schema_version = 2L,
    metric_status = c("ok", "zero_denominator")
  )
  averaged <- facspseudocolor:::average_edu_metric_table(
    acquisition, "regional_edu_positive_pct",
    c("region", "region_label", "region_index")
  )
  expect_equal(averaged$regional_edu_positive_pct, 10)
  expect_identical(averaged$metric_status, "partial")
  expect_equal(averaged$non_ok_technical_n, 1)
  expect_identical(averaged$technical_metric_statuses,
                   "ok;zero_denominator")
})

test_that("modern and explicitly requested legacy signal populations remain separate", {
  analysis <- synthetic_edu_contract_analysis()
  expect_warning(
    result <- quantify_cell_cycle(analysis, signal = "normalized"),
    "deprecated"
  )
  expect_true("normalized" %in% names(result$quantitation$by_signal))
  expect_true("legacy_background_divided" %in% names(result$quantitation))
  expect_identical(
    result$quantitation$edu_positive_population_intensity$
      source_population,
    rep("whole_computed_positive_eligible_population", 2)
  )
  expect_false(any(grepl("flowjo|figure1", names(result$quantitation),
                         ignore.case = TRUE)))
})

test_that("CSV export uses canonical names and requires explicit legacy selection", {
  result <- quantify_synthetic_edu_contract()
  directory <- withr::local_tempdir()
  written <- save_facs_results(
    result, output_pdf = NULL, output_png = NULL,
    output_csv_dir = directory
  )
  expected <- c(
    facspseudocolor:::edu_canonical_table_names(),
    paste0(facspseudocolor:::edu_canonical_table_names(), "_acquisition")
  )
  expect_true(all(file.exists(file.path(directory, paste0(expected, ".csv")))))
  expect_false(file.exists(file.path(directory, "phase_percentages.csv")))
  expect_true(all(paste0("csv_", expected) %in% names(written)))

  legacy_directory <- withr::local_tempdir()
  save_facs_results(
    result, output_pdf = NULL, output_png = NULL,
    output_csv_dir = legacy_directory, include_deprecated_csv = TRUE
  )
  expect_true(file.exists(file.path(legacy_directory, "phase_percentages.csv")))
  expect_true(file.exists(file.path(legacy_directory, "phase_medians.csv")))
  expect_true(file.exists(file.path(legacy_directory, "whole_medians.csv")))
})
