# Every value in this file is explicitly SYNTHETIC and test-only.

test_that("SYNTHETIC raw-4N cutoff uses the fixed density method and first valley", {
  negative <- seq(100, 130, length.out = 160L)
  positive <- seq(260, 290, length.out = 40L)
  result <- ph3_raw_4n_density_cutoff(c(negative, positive))
  expect_identical(result$method, ph3_raw_4n_cutoff_method())
  expect_identical(result$control_event_count, 200L)
  expect_true(result$cutoff_index > result$peak_index)
  expect_true(is.finite(result$cutoff))
  expect_true(result$cutoff >= result$raw_min)
  expect_true(result$cutoff <= result$raw_max)
})

test_that("SYNTHETIC raw-4N cutoff fails closed for insufficient, nonfinite, and flat inputs", {
  expect_error(
    ph3_raw_4n_density_cutoff(seq_len(99L)),
    "insufficient_control_4n_events"
  )
  expect_error(
    ph3_raw_4n_density_cutoff(c(seq_len(99L), Inf)),
    "nonfinite_control_raw_signal"
  )
  expect_error(
    ph3_raw_4n_density_cutoff(rep(100, 100L)),
    "zero_or_invalid_control_signal_range"
  )
})

test_that("SYNTHETIC raw-4N cutoff never substitutes a FlowJo membership label", {
  values <- c(seq(100, 130, length.out = 160L), seq(260, 290, length.out = 40L))
  result <- ph3_raw_4n_density_cutoff(values)
  computed <- values > result$cutoff
  flowjo_label <- rep(c(TRUE, FALSE), length.out = length(values))
  expect_false(identical(computed, flowjo_label))
  expect_identical(computed, is.finite(values) & values > result$cutoff)
})
