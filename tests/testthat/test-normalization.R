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

test_that("EdU regression uses a shared slope through the sample G1 anchor", {
  data_dir <- withr::local_tempdir()
  write_test_csv(file.path(data_dir, "sample_single.csv"), c(100, 150, 200), c(10, 15, 20))
  write_test_csv(file.path(data_dir, "sample_g1.csv"), c(90, 100, 110), c(9, 10, 11))
  write_test_csv(file.path(data_dir, "sample_positive.csv"), c(125, 175), c(25, 35))
  result <- read_and_normalize_sample(
    "sample", "Sample", data_dir,
    list(complete = "_single.csv", g1 = "_g1.csv", edu_positive = "_positive.csv"),
    list(dna_channel = "DNA", target_channel = "Target", dna_2n_value = 1000, g1_anchor = "median"),
    baseline_slope = 0.01, normalization_method = "reference_negative_regression")
  expect_equal(result$data$baseline, c(10, 15, 20))
  expect_equal(result$data$target_norm, rep(1000, 3))
  expect_equal(result$edu_positive$target_bgsub, c(12.5, 17.5))
  expect_error(read_and_normalize_sample(
    "sample", "Sample", data_dir, list(complete = "_single.csv", g1 = "_g1.csv"),
    list(dna_channel = "DNA", target_channel = "Target"),
    normalization_method = "reference_negative_regression"), "baseline_slope.*finite number")
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
