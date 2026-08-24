test_that("high-level EdU analysis returns a structured non-saving result", {
  config <- read_facs_config(system.file(
    "config", "config_edu.yml", package = "facspseudocolor"
  ))
  expect_warning(
    result <- analyze_facs_experiment(config),
    "aliases.*deprecated"
  )

  expect_s3_class(result, "facs_analysis")
  expect_equal(nrow(result$sample_manifest), 6)
  expect_length(result$models, 6)
  expect_length(result$normalized_data, 6)
  expect_identical(names(result$normalized_data), result$sample_manifest$prefix)
  expect_equal(
    vapply(result$normalized_data, function(x) {
      stats::median(x$data$target_norm, na.rm = TRUE)
    }, numeric(1)),
    c(rep1_NT = 1752.82011390766, rep1_2h = 1504.51043111236,
      rep1_4h = 1368.51065076182, rep2_NT = 2349.03387094802,
      rep2_2h = 1643.1981974707, rep2_4h = 1472.14509202877),
    tolerance = 1e-8
  )
  expect_identical(names(result$quantitation$by_signal),
                   "background_subtracted")
  expect_true(all(c("phase_medians", "whole_medians") %in%
                    names(result$quantitation$by_signal$background_subtracted)))
  expect_true("phase_percentages" %in% names(result$quantitation))
  expect_identical(result$provenance$package, "facspseudocolor")
})

test_that("high-level POI analysis preserves current normalized results", {
  config <- read_facs_config(system.file(
    "config", "config_poi.yml", package = "facspseudocolor"
  ))
  result <- analyze_facs_experiment(config)

  expect_s3_class(result, "facs_analysis")
  expect_length(result$models, 2)
  expect_equal(
    vapply(result$normalized_data, function(x) {
      stats::median(x$data$target_norm, na.rm = TRUE)
    }, numeric(1)),
    c(rep1_bg = 956.742596862670, rep1_NT = 4991.65220845589,
      rep1_2h = 1638.55373071234, rep1_4h = 1826.48029208043,
      rep2_bg = 968.992295650424, rep2_NT = 4098.33453848040,
      rep2_2h = 1648.79195919279, rep2_4h = 1689.28981136966),
    tolerance = 1e-8
  )
})

test_that("analysis print method summarizes without exposing event data", {
  object <- structure(list(
    config = list(plot_type = "poi"),
    sample_manifest = data.frame(replicate_index = c(1, 1)),
    models = list(list()), warnings = character()
  ), class = "facs_analysis")
  expect_output(print(object), "mode:       poi")
  expect_output(print(object), "samples:    2")
})
