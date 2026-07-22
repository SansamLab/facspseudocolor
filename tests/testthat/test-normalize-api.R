test_that("normalize_edu performs calculation without file I/O", {
  events <- data.frame(DNA = c(100, 150, 200), Target = c(10, 15, 20))
  g1 <- data.frame(DNA = c(90, 100, 110), Target = c(9, 10, 11))
  positive <- data.frame(DNA = c(125, 175), Target = c(25, 35))

  result <- normalize_edu(
    events, g1, "DNA", "Target", baseline_slope = 0.01,
    edu_positive_events = positive, sample_id = "test"
  )

  expect_equal(result$data$dna_norm, c(1000, 1500, 2000))
  expect_equal(result$data$baseline, c(10, 15, 20))
  expect_equal(result$data$target_norm, rep(1000, 3))
  expect_equal(result$edu_positive$target_bgsub, c(12.5, 17.5))
  expect_false("dna_norm" %in% names(events))
  expect_false("dna_norm" %in% names(g1))
})

test_that("normalize_poi performs calculation without file I/O", {
  events <- data.frame(
    DNA = c(100, 120, 140, 160, 180, 200, 220, 240, 260, 280),
    Target = c(20, 22, 24, 26, 28, 30, 32, 34, 36, 38)
  )
  model <- list(dna_peak = 100, intercept = 10, slope = 0.01,
                floor = 1, cutoff = 1200)
  result <- normalize_poi(
    events, model, "DNA", "Target", dna_align = "shared_background"
  )

  expect_equal(result$data$dna_norm, events$DNA * 10)
  expect_equal(result$data$baseline, 10 + 0.01 * result$data$dna_norm)
  expect_equal(result$cutoff, 1200)
  expect_false("dna_norm" %in% names(events))
})

test_that("POI peak failure stops unless legacy fallback is explicit", {
  events <- data.frame(DNA = rep(100, 10), Target = seq_len(10) + 10)
  model <- list(dna_peak = 80, intercept = 10, slope = 0.01, floor = 1)

  expect_error(
    normalize_poi(events, model, "DNA", "Target", sample_id = "constant DNA"),
    "peak detection failed.*constant DNA"
  )
  fallback <- normalize_poi(
    events, model, "DNA", "Target", peak_failure = "use_background"
  )
  expect_equal(fallback$dna_peak, 80)
  expect_identical(fallback$peak_failure, "use_background")
})
