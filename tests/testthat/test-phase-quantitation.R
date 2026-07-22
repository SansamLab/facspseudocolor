test_that("default phase windows are ordered and scaled from the 2N value", {
  windows <- make_phase_windows(dna_2n_value = 900)
  expect_identical(windows$phase_label, c("G1", "Early S", "Mid S", "Late S", "G2/M"))
  expect_equal(windows$x_min[2:4], c(900, 1200, 1500))
  expect_equal(windows$x_max[2:4], c(1200, 1500, 1800))
  expect_error(validate_s_phase_bins(list(early = 1:2)), "early, mid, and late")
})

test_that("rectangular EdU gates separate low-signal and S-phase bands", {
  gates <- make_rectangular_phase_gates(dna_2n_value = 1000, y_limits = c(500, 80000))
  expect_equal(nrow(gates), 5)
  expect_equal(gates$ymax[c(1, 5)], c(2500, 2500))
  expect_equal(gates$ymin[2:4], rep(2500, 3))
})

test_that("events use half-open gate boundaries and retain unmatched events", {
  gates <- data.frame(gate = c("A", "B"), xmin = c(0, 1), xmax = c(1, 2), ymin = c(0, 0), ymax = c(2, 2))
  events <- data.frame(dna_norm = c(0, 0.999, 1, 2, NA), target_norm = rep(1, 5))
  assigned <- assign_events_to_gates(events, gates)
  expect_identical(as.character(assigned$gate), c("A", "A", "B", "ungated", "ungated"))
})

test_that("phase medians enforce the minimum event count", {
  windows <- data.frame(phase = c("early", "late"), phase_label = c("Early S", "Late S"),
                        x_min = c(0, 2), x_max = c(2, 4), phase_index = 1:2)
  dat <- data.frame(dna_norm = c(0.5, 1.5, 2.5), target_norm = c(10, 20, 100))
  result <- calculate_phase_signal_medians(dat, windows, minimum_events = 2)
  expect_equal(result$median_signal[[1]], 15)
  expect_true(is.na(result$median_signal[[2]]))
  expect_identical(result$n, c(2L, 1L))
})

test_that("replicate summaries and reference ratios are computed within groups", {
  dat <- data.frame(replicate_index = c(1, 1, 2, 2), condition = rep(c("Control", "Drug"), 2),
                    value = c(10, 20, 5, 15))
  ratio <- add_reference_ratio(dat, "value", "Control")
  expect_equal(ratio$ratio, c(1, 2, 1, 3))
  summary <- summarize_across_replicates(ratio, "ratio", "condition")
  drug <- summary[summary$condition == "Drug", ]
  expect_equal(drug$replicate_n, 2)
  expect_equal(drug$mean, 2.5)
  expect_error(add_reference_ratio(dat, "value", "Missing"), "not among the quantified samples")
})
