test_that("POI cutoff uses the type-7 quantile of corrected background", {
  data_dir <- withr::local_tempdir()

  # SYNTHETIC: the target residuals sum to zero and have zero covariance with
  # DNA, so the fitted background is exactly 10 with slope zero. The corrected
  # values are therefore 100 times the target values.
  SYNTHETIC_background <- data.frame(
    DNA = seq(95, 106),
    Target = c(11, rep(10, 4), 49 / 6, rep(10, 5), 65 / 6)
  )
  utils::write.csv(
    SYNTHETIC_background,
    file.path(data_dir, "SYNTHETIC_background_single.csv"),
    row.names = FALSE
  )

  # Peak estimation is unrelated to the cutoff convention under test. Fixing
  # it at 100 makes the DNA normalization independently hand-checkable.
  testthat::local_mocked_bindings(
    estimate_lower_dna_peak = function(x, adjust = 0.75) 100,
    .package = "facspseudocolor"
  )

  model <- fit_background_model(
    prefix = "SYNTHETIC_background",
    condition_label = "SYNTHETIC background",
    replicate_label = "SYNTHETIC replicate",
    data_dir = data_dir,
    file_suffixes = list(complete = "_single.csv"),
    settings = list(
      dna_channel = "DNA",
      target_channel = "Target",
      dna_2n_value = 1000,
      background_quantile = 0.95
    )
  )

  # For 12 sorted corrected values, type 7 places p = 0.95 at index 11.45.
  # Interpolating 45% from 3250/3 to 1100 gives 6545/6 = 1090.8333...
  expect_equal(model$intercept, 10, tolerance = 1e-12)
  expect_equal(model$slope, 0, tolerance = 1e-12)
  expect_equal(model$cutoff, 6545 / 6, tolerance = 1e-12)
})
