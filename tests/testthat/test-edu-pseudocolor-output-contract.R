# All events in this file are unmistakably SYNTHETIC and test-only.

synthetic_edu_display_analysis <- function() {
  raw_config <- minimal_config("edu")
  raw_config$replicates[[1L]]$reference <- "Reference"
  config <- validate_facs_config(raw_config)
  config$y_log10 <- FALSE
  config$x_limits <- c(700, 2250)
  manifest <- build_sample_manifest(config)
  make_sample <- function(raw_shift, baseline) {
    dna <- seq(800, 2200, length.out = 24L)
    raw <- seq(100, 2400, length.out = 24L) + raw_shift
    bgsub <- raw - baseline
    negative <- seq_len(12L)
    list(
      sample = list(data = data.frame(
        dna_norm = dna, target_raw = raw, target_bgsub = bgsub,
        target_norm = raw, baseline = rep(baseline, length(raw)),
        edu_computed_positive = seq_along(raw) > 12L
      )),
      model = list(
        negative_event_index = negative,
        negative_event_fingerprint = edu_negative_event_fingerprint(dna[negative], raw[negative])
      )
    )
  }
  reference <- make_sample(100, 60)
  treatment <- make_sample(600, 160)
  samples <- list(reference$sample, treatment$sample)
  names(samples) <- manifest$prefix
  models <- list(reference$model, treatment$model)
  for (i in seq_along(models)) models[[i]]$sample_prefix <- manifest$prefix[[i]]
  names(models) <- manifest$prefix
  analysis <- new_facs_analysis(
    config, manifest, data.frame(path = character(), exists = logical()),
    samples, models,
    quantitation = list(SYNTHETIC_quantitation = data.frame(value = c(1, 2)))
  )
  # The report consumes the existing canonical EdU positivity tables and the
  # already-classified event rows. The warning is the explicitly tested legacy
  # alias deprecation, not an outcome of this SYNTHETIC fixture.
  suppressWarnings(quantify_cell_cycle(analysis))
}

test_that("SYNTHETIC EdU panels use exact per-sample negative restoration offsets", {
  analysis <- synthetic_edu_display_analysis()
  before <- serialize(analysis$quantitation, NULL)
  result <- build_edu_pseudocolor_output_contract(analysis)
  qc <- result$display_offset_qc

  expect_identical(names(result$panels), analysis$sample_manifest$prefix)
  expect_true(all(qc$display_status == "available"))
  expect_equal(qc$display_offset, c(60, 160))
  expect_equal(qc$raw_negative_median - qc$corrected_negative_median,
               qc$display_offset)
  expect_identical(serialize(analysis$quantitation, NULL), before)
  expect_true(all(result$panel_qc$displayed_event_n > 0L))
  expect_true(all(result$panel_qc$y_limit_lower < result$panel_qc$y_limit_upper))
  expect_identical(result$provenance$analytical_values_mutated, FALSE)
})

test_that("SYNTHETIC EdU display offset cannot be reassigned after model reordering", {
  analysis <- synthetic_edu_display_analysis()
  analysis$models <- analysis$models[rev(names(analysis$models))]
  names(analysis$models) <- analysis$sample_manifest$prefix
  result <- build_edu_pseudocolor_output_contract(analysis)
  expect_true(all(result$display_offset_qc$display_status == "suppressed"))
  expect_true(all(result$display_offset_qc$suppression_reason_code == "sample_identity_mismatch"))
})

test_that("SYNTHETIC EdU missing or mutated membership proof suppresses only that panel", {
  analysis <- synthetic_edu_display_analysis()
  analysis$models$reference$negative_event_fingerprint <- "SYNTHETIC-mutated"
  result <- build_edu_pseudocolor_output_contract(analysis)
  qc <- result$display_offset_qc
  expect_identical(qc$display_status, c("suppressed", "available"))
  expect_identical(qc$suppression_reason_code[[1L]], "negative_membership_proof_mismatch")
  expect_identical(result$panel_qc$panel_status, c("suppressed", "available"))
})

test_that("SYNTHETIC EdU display labels retain condition replicate and technical acquisition", {
  analysis <- synthetic_edu_display_analysis()
  result <- build_edu_pseudocolor_output_contract(analysis)
  labels <- vapply(result$panels, function(plot) plot$labels$subtitle, character(1))
  expect_true(all(grepl("Condition:|Biological replicate:|Technical acquisition:", labels)))
})

test_that("SYNTHETIC EdU report retains approved 2N-4N and canonical regional positivity", {
  analysis <- synthetic_edu_display_analysis()
  result <- build_edu_pseudocolor_output_contract(analysis)
  positivity <- result$positivity

  expect_identical(result$schema_version, "edu-pseudocolor-output-contract-1.2.0")
  expect_identical(names(positivity$panels), c("two_to_four_n", "early_mid_late_s"))
  expect_identical(positivity$overall_2to4n_biological_replicate$region_label,
                   c("2N\u20134N", "2N\u20134N"))
  expect_true(all(positivity$overall_2to4n_acquisition$dna_interval_lower_inclusive))
  expect_true(all(positivity$overall_2to4n_acquisition$dna_interval_upper_inclusive))
  # The user-approved report denominator includes the 4N endpoint, unlike the
  # half-open Early/Mid/Late-S regional source tables.
  expected_denominators <- vapply(analysis$normalized_data, function(sample) {
    x <- sample$data$dna_norm
    positive_known <- !is.na(sample$data$edu_computed_positive)
    sum(positive_known & x >= analysis$config$dna_2n_value &
          x <= 2 * analysis$config$dna_2n_value)
  }, integer(1))
  expect_identical(positivity$overall_2to4n_acquisition$denominator_n,
                   unname(expected_denominators))
  expect_setequal(positivity$regional_biological_replicate$region,
                  c("early", "mid", "late"))
  expect_match(positivity$provenance$overall_definition, "inclusive")
})

test_that("SYNTHETIC EdU positivity report rejects nonreconciling canonical regions", {
  analysis <- synthetic_edu_display_analysis()
  analysis$quantitation$edu_regional_positivity$regional_edu_positive_pct[[1L]] <-
    analysis$quantitation$edu_regional_positivity$regional_edu_positive_pct[[1L]] + 1
  expect_error(
    build_edu_pseudocolor_output_contract(analysis),
    "regional_positivity_reconciliation_failed"
  )
})

test_that("SYNTHETIC EdU regional positivity rejects substituted acquisition identity", {
  analysis <- synthetic_edu_display_analysis()
  analysis$quantitation$edu_regional_positivity_acquisition$
    technical_replicate[[1L]] <- "SYNTHETIC substituted acquisition"
  expect_error(
    build_edu_pseudocolor_output_contract(analysis),
    "acquisition_manifest_identity_mismatch"
  )
})

test_that("SYNTHETIC EdU regional positivity rejects omitted acquisition identity", {
  analysis <- synthetic_edu_display_analysis()
  analysis$quantitation$edu_regional_positivity_acquisition <-
    analysis$quantitation$edu_regional_positivity_acquisition[-1L, , drop = FALSE]
  expect_error(
    build_edu_pseudocolor_output_contract(analysis),
    "acquisition_manifest_identity_mismatch"
  )
})

test_that("SYNTHETIC EdU regional positivity rejects duplicate acquisition identity", {
  analysis <- synthetic_edu_display_analysis()
  acquisition <- analysis$quantitation$edu_regional_positivity_acquisition
  analysis$quantitation$edu_regional_positivity_acquisition <-
    rbind(acquisition, acquisition[1L, , drop = FALSE])
  expect_error(
    build_edu_pseudocolor_output_contract(analysis),
    "acquisition_manifest_identity_mismatch"
  )
})

test_that("SYNTHETIC EdU intensity report uses the configured reference after canonical aggregation", {
  analysis <- synthetic_edu_display_analysis()
  result <- build_edu_pseudocolor_output_contract(analysis)
  intensity <- result$intensity
  canonical_all <- analysis$quantitation$edu_positive_population_intensity
  canonical_regional <- analysis$quantitation$edu_positive_cell_regional_intensity

  expect_identical(names(intensity$panels), c("all_computed_positive", "early_mid_late_s"))
  expect_identical(
    intensity$overall_canonical_biological_replicate$positive_population_edu_bgsub_median,
    canonical_all$positive_population_edu_bgsub_median
  )
  expect_identical(
    intensity$regional_canonical_biological_replicate$positive_cell_regional_edu_bgsub_median,
    canonical_regional$positive_cell_regional_edu_bgsub_median
  )
  expect_true(all(intensity$regional_biological_replicate$signal_transform ==
                    "background_subtracted"))
  available_all <- intensity$overall_biological_replicate$
    reference_normalization_status == "available"
  available_regional <- intensity$regional_biological_replicate$
    reference_normalization_status == "available"
  expect_true(all(intensity$overall_biological_replicate$
                    reference_relative_intensity[available_all &
                      intensity$overall_biological_replicate$condition == "Reference"] == 1))
  expect_true(all(intensity$regional_biological_replicate$
                    reference_relative_intensity[available_regional &
                      intensity$regional_biological_replicate$condition == "Reference"] == 1))
  expect_true(all(intensity$overall_biological_replicate$
                    reference_normalization_status %in% c(
                      "available", "unavailable_invalid_reference_intensity",
                      "unavailable_nonfinite_canonical_intensity"
                    )))
  expect_identical(intensity$provenance$signal, "background_subtracted")
  expect_match(intensity$provenance$reference_normalization, "Untreated")
})

test_that("SYNTHETIC EdU intensity report rejects absent references and exposes an invalid denominator", {
  analysis <- synthetic_edu_display_analysis()
  analysis$sample_manifest$is_reference <- FALSE
  expect_error(
    build_edu_pseudocolor_output_contract(analysis),
    "invalid_reference_identity"
  )

  analysis <- synthetic_edu_display_analysis()
  reference_row <- which(analysis$quantitation$edu_positive_population_intensity$condition == "Reference")
  analysis$quantitation$edu_positive_population_intensity$
    positive_population_edu_bgsub_median[[reference_row]] <- 0
  analysis$quantitation$edu_positive_population_intensity_acquisition$
    positive_population_edu_bgsub_median[[reference_row]] <- 0
  result <- build_edu_pseudocolor_output_contract(analysis)
  rows <- result$intensity$overall_biological_replicate
  expect_true(all(rows$reference_normalization_status ==
                    "unavailable_invalid_reference_intensity"))
  expect_true(all(is.na(rows$reference_relative_intensity)))
  expect_true(all(rows$reference_intensity == 0))
})

test_that("SYNTHETIC EdU intensity report exposes nonfinite canonical numerators without fallback", {
  analysis <- synthetic_edu_display_analysis()
  treatment_row <- which(analysis$quantitation$edu_positive_population_intensity$condition == "Treatment")
  analysis$quantitation$edu_positive_population_intensity_acquisition$
    positive_population_edu_bgsub_median[[treatment_row]] <- NA_real_
  analysis$quantitation$edu_positive_population_intensity <-
    average_edu_metric_table(
      analysis$quantitation$edu_positive_population_intensity_acquisition,
      "positive_population_edu_bgsub_median"
    )
  result <- build_edu_pseudocolor_output_contract(analysis)
  rows <- result$intensity$overall_biological_replicate
  treatment <- rows[rows$condition == "Treatment", , drop = FALSE]
  expect_identical(treatment$reference_normalization_status,
                   "unavailable_nonfinite_canonical_intensity")
  expect_true(is.na(treatment$reference_relative_intensity))
  expect_true(is.finite(treatment$reference_intensity))
})

test_that("SYNTHETIC EdU intensity report rejects nonreconciling canonical values", {
  analysis <- synthetic_edu_display_analysis()
  analysis$quantitation$edu_positive_population_intensity$
    positive_population_edu_bgsub_median[[1L]] <- 12345
  expect_error(
    build_edu_pseudocolor_output_contract(analysis),
    "intensity_reconciliation_failed"
  )
})

test_that("SYNTHETIC EdU intensity report rejects substituted acquisition identity", {
  analysis <- synthetic_edu_display_analysis()
  analysis$quantitation$edu_positive_population_intensity_acquisition$
    condition[[1L]] <- "SYNTHETIC substituted condition"
  expect_error(
    build_edu_pseudocolor_output_contract(analysis),
    "acquisition_manifest_identity_mismatch"
  )
})

test_that("SYNTHETIC EdU regional intensity rejects substituted acquisition identity", {
  analysis <- synthetic_edu_display_analysis()
  analysis$quantitation$edu_positive_cell_regional_intensity_acquisition$
    technical_replicate[[1L]] <- "SYNTHETIC substituted acquisition"
  expect_error(
    build_edu_pseudocolor_output_contract(analysis),
    "acquisition_manifest_identity_mismatch"
  )
})

test_that("SYNTHETIC EdU regional intensity rejects omitted acquisition identity", {
  analysis <- synthetic_edu_display_analysis()
  analysis$quantitation$edu_positive_cell_regional_intensity_acquisition <-
    analysis$quantitation$edu_positive_cell_regional_intensity_acquisition[-1L, , drop = FALSE]
  expect_error(
    build_edu_pseudocolor_output_contract(analysis),
    "acquisition_manifest_identity_mismatch"
  )
})

test_that("SYNTHETIC EdU regional intensity rejects duplicate acquisition identity", {
  analysis <- synthetic_edu_display_analysis()
  acquisition <- analysis$quantitation$edu_positive_cell_regional_intensity_acquisition
  analysis$quantitation$edu_positive_cell_regional_intensity_acquisition <-
    rbind(acquisition, acquisition[1L, , drop = FALSE])
  expect_error(
    build_edu_pseudocolor_output_contract(analysis),
    "acquisition_manifest_identity_mismatch"
  )
})
