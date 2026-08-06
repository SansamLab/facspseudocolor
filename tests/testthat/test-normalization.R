test_that("G1-median normalization anchors DNA and target at 1000", {
  data_dir <- withr::local_tempdir()
  write_test_csv(file.path(data_dir, "sample_single.csv"), c(100, 150, 200), c(10, 20, 30))
  write_test_csv(file.path(data_dir, "sample_g1.csv"), c(90, 100, 110), c(8, 10, 12))
  result <- read_and_normalize_sample(
    "sample", "Sample", data_dir, list(complete = "_single.csv", g1 = "_g1.csv"),
    list(dna_channel = "DNA", target_channel = "Target", dna_2n_value = 1000, normalize_target = TRUE),
    normalization_method = "g1_median")
  expect_equal(result$g1_median_dna, 100)
  expect_equal(result$g1_anchor_target, 10)
  expect_equal(result$data$dna_norm, c(1000, 1500, 2000))
  expect_equal(result$data$target_norm, c(1000, 2000, 3000))
})

test_that("EdU normalization uses the supplied slope through the sample G1 anchor", {
  data_dir <- withr::local_tempdir()
  write_test_csv(file.path(data_dir, "sample_single.csv"), c(100, 150, 200), c(10, 15, 20))
  write_test_csv(file.path(data_dir, "sample_g1.csv"), c(90, 100, 110), c(9, 10, 11))
  write_test_csv(file.path(data_dir, "sample_positive.csv"), c(125, 175), c(25, 35))
  result <- read_and_normalize_sample(
    "sample", "Sample", data_dir,
    list(complete = "_single.csv", g1 = "_g1.csv", edu_positive = "_positive.csv"),
    list(dna_channel = "DNA", target_channel = "Target", dna_2n_value = 1000, g1_anchor = "median"),
    baseline_slope = 0.01, normalization_method = "reference_negative_regression")
  expect_equal(result$data$dna_norm, c(1000, 1500, 2000))
  expect_equal(result$data$baseline, c(10, 15, 20))
  expect_equal(result$data$target_norm, rep(1000, 3))
  expect_equal(result$data$target_bgsub, c(0, 0, 0))
  expect_equal(result$edu_positive$target_bgsub, c(12.5, 17.5))
  expect_error(read_and_normalize_sample(
    "sample", "Sample", data_dir, list(complete = "_single.csv", g1 = "_g1.csv"),
    list(dna_channel = "DNA", target_channel = "Target"),
    normalization_method = "reference_negative_regression"), "baseline_slope.*finite number")
})

test_that("EdU boundary bins include endpoints and assign internal breaks right", {
  SYNTHETIC_lower_and_midpoint <- data.frame(
    dna_norm = c(1000, 1499, 1500, 2000),
    target = c(1, 5, 2, 4)
  )
  lower_and_midpoint <- fit_positive_minimum_boundary(
    SYNTHETIC_lower_and_midpoint,
    y_column = "target",
    fit_x_range = c(1000, 2000),
    boundary_bins = 2,
    minimum_events_per_bin = 1
  )

  expect_identical(lower_and_midpoint$points$bin, c(1L, 2L))
  expect_equal(lower_and_midpoint$points$boundary_x, c(1000, 1500))
  expect_equal(lower_and_midpoint$points$boundary_y, c(1, 2))

  SYNTHETIC_upper_and_outside <- data.frame(
    dna_norm = c(999, 1000, 1499, 1500, 2000, 2001),
    target = c(-100, 5, 4, 3, 2, -200)
  )
  upper_and_outside <- fit_positive_minimum_boundary(
    SYNTHETIC_upper_and_outside,
    y_column = "target",
    fit_x_range = c(1000, 2000),
    boundary_bins = 2,
    minimum_events_per_bin = 1
  )

  expect_identical(upper_and_outside$points$bin, c(1L, 2L))
  expect_equal(upper_and_outside$points$boundary_x, c(1499, 2000))
  expect_equal(upper_and_outside$points$boundary_y, c(4, 2))
})

test_that("EdU boundary bins select the first row when minima are tied", {
  SYNTHETIC_tied_minima <- data.frame(
    dna_norm = c(1100, 1200, 1300, 1400, 1600, 1700),
    target = c(20, 10, 10, 30, 40, 50)
  )
  boundary <- fit_positive_minimum_boundary(
    SYNTHETIC_tied_minima,
    y_column = "target",
    fit_x_range = c(1000, 2000),
    boundary_bins = 2,
    minimum_events_per_bin = 1
  )

  # DNA 1200 and 1300 tie at target 10; which.min() selects the first row.
  expect_identical(boundary$points$bin, c(1L, 2L))
  expect_equal(boundary$points$boundary_x, c(1200, 1600))
  expect_equal(boundary$points$boundary_y, c(10, 40))
})

test_that("EdU boundary fitting excludes nonfinite rows and requires enough valid events", {
  SYNTHETIC_with_nonfinite <- data.frame(
    dna_norm = c(1100, 1200, Inf, -Inf, 1600, 1700, 1800),
    target = c(10, NA_real_, 5, 6, 20, NaN, -Inf)
  )
  boundary <- fit_positive_minimum_boundary(
    SYNTHETIC_with_nonfinite,
    y_column = "target",
    fit_x_range = c(1000, 2000),
    boundary_bins = 2,
    minimum_events_per_bin = 1
  )

  expect_identical(boundary$points$bin, c(1L, 2L))
  expect_equal(boundary$points$boundary_x, c(1100, 1600))
  expect_equal(boundary$points$boundary_y, c(10, 20))

  SYNTHETIC_too_few_finite <- data.frame(
    dna_norm = c(1100, 1200, Inf, 1600),
    target = c(10, NA_real_, 20, NaN)
  )
  expect_error(
    fit_positive_minimum_boundary(
      SYNTHETIC_too_few_finite,
      y_column = "target",
      fit_x_range = c(1000, 2000),
      boundary_bins = 2,
      minimum_events_per_bin = 1
    ),
    "Too few EdU\\+ events for boundary fitting"
  )
})

test_that("EdU negative regression uses strict boundary and inclusive DNA endpoints", {
  data_dir <- withr::local_tempdir()
  SYNTHETIC_complete <- data.frame(
    DNA = c(100, 100, 100, 90, 200, 210),
    Target = c(99, 100, 101, 50, 98, 50)
  )
  SYNTHETIC_g1 <- data.frame(
    DNA = c(90, 100, 110),
    Target = c(90, 100, 110)
  )
  SYNTHETIC_edu_positive <- data.frame(
    DNA = c(100, 150),
    Target = c(100, 100)
  )
  utils::write.csv(
    SYNTHETIC_complete,
    file.path(data_dir, "SYNTHETIC_sample_single.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    SYNTHETIC_g1,
    file.path(data_dir, "SYNTHETIC_sample_g1.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    SYNTHETIC_edu_positive,
    file.path(data_dir, "SYNTHETIC_sample_positive.csv"),
    row.names = FALSE
  )

  model <- fit_reference_negative_model(
    prefix = "SYNTHETIC_sample",
    condition_label = "SYNTHETIC condition",
    replicate_label = "SYNTHETIC replicate",
    data_dir = data_dir,
    file_suffixes = list(
      complete = "_single.csv",
      g1 = "_g1.csv",
      edu_positive = "_positive.csv"
    ),
    settings = list(
      dna_channel = "DNA",
      target_channel = "Target",
      dna_2n_value = 1000,
      g1_anchor = "median",
      baseline_fit_x_range = c(1000, 2000),
      baseline_boundary_bins = 2,
      baseline_minimum_events_per_bin = 1,
      baseline_minimum_negative_events = 2
    )
  )

  expect_equal(model$positive_boundary_intercept, 100)
  expect_equal(model$positive_boundary_slope, 0)
  expect_identical(model$negative_event_n, 2L)
  expect_equal(model$negative_intercept, 100)
  expect_equal(model$slope, -0.001)
})

test_that("POI normalization applies the supplied background regression", {
  data_dir <- withr::local_tempdir()
  write_test_csv(file.path(data_dir, "sample_single.csv"), c(100, 150, 200), c(20, 25, 30))
  model <- list(dna_peak = 100, intercept = 10, slope = 0.01, floor = 1, cutoff = 1200)
  result <- read_and_normalize_sample(
    "sample", "Sample", data_dir, list(complete = "_single.csv"),
    list(dna_channel = "DNA", target_channel = "Target", dna_2n_value = 1000, poi_dna_align = "shared_background"),
    background_model = model, normalization_method = "background_reference_regression")
  expect_equal(result$data$dna_norm, c(1000, 1500, 2000))
  expect_equal(result$data$baseline, c(20, 25, 30))
  expect_equal(result$data$target_norm, rep(1000, 3))
  expect_equal(result$cutoff, 1200)
})

test_that("normalization reports missing files and channels", {
  data_dir <- withr::local_tempdir()
  settings <- list(dna_channel = "DNA", target_channel = "Target")
  suffixes <- list(complete = "_single.csv", g1 = "_g1.csv")
  expect_error(read_and_normalize_sample("missing", "Missing sample", data_dir, suffixes, settings,
                                         normalization_method = "g1_median"), "Missing files")
  utils::write.csv(data.frame(DNA = 1:3), file.path(data_dir, "bad_single.csv"), row.names = FALSE)
  write_test_csv(file.path(data_dir, "bad_g1.csv"), 1:3, 1:3)
  expect_error(read_and_normalize_sample("bad", "Bad sample", data_dir, suffixes, settings,
                                         normalization_method = "g1_median"), "Missing channels")
})
