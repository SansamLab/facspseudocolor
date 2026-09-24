# All event tables in this file are unmistakably SYNTHETIC and test-only.

test_that("synchronized defaults preserve the target and require G1 inputs", {
  config <- validate_facs_config(minimal_config("synchronized"))
  expect_identical(config$synchronized_dna_strategy, "per_sample_g1")
  expect_equal(config$synchronized_minimum_g1_events, 100)
  expect_false(config$normalize_target)
  expect_identical(config$suffixes,
                   list(complete = "_single_cells.csv", g1 = "_g1.csv"))
})

test_that("synchronized mode validates an explicit shared reference", {
  raw <- minimal_config("synchronized")
  raw$synchronized_dna_strategy <- "shared_asynchronous_g1"
  raw$synchronized_reference_prefix <- "reference"
  config <- validate_facs_config(raw)
  inputs <- facspseudocolor:::facs_input_files(config, withr::local_tempdir())
  expect_equal(sum(inputs$population == "complete"), 2)
  expect_identical(inputs$prefix[inputs$population == "g1"], "reference")
  raw$synchronized_reference_prefix <- "absent"
  expect_error(validate_facs_config(raw), "identify exactly one configured sample")
})

test_that("synchronized normalization preserves zero and negative target values", {
  events <- data.frame(DNA = c(100, 150, 200), Target = c(-5, 0, 10))
  result <- facspseudocolor:::synchronized_normalize_table(
    events, "DNA", "Target", anchor = 100, dna_2n_value = 1000
  )
  expect_equal(result$dna_norm, c(1000, 1500, 2000))
  expect_identical(result$target_norm, events$Target)
  expect_identical(result$target_bgsub, events$Target)
})

test_that("synchronized G1 containment rejects an event outside Single Cells", {
  complete <- data.frame(
    event_index = seq_len(120),
    DNA = rep(c(100, 200), each = 60), Target = seq_len(120)
  )
  g1 <- complete[seq_len(100), , drop = FALSE]
  expect_identical(
    facspseudocolor:::synchronized_validate_g1_containment(
      complete, g1, "SYNTHETIC sample", 100,
      required_columns = "event_index"
    )$status,
    "validated"
  )
  unmatched_child <- g1
  unmatched_child$event_index[[100]] <- 121
  expect_error(
    facspseudocolor:::synchronized_validate_g1_containment(
      complete, unmatched_child, "SYNTHETIC sample", 100,
      required_columns = "event_index"
    ), "not an exact multiset subset"
  )
})

test_that("synchronized G1 containment requires a unique stable event index", {
  complete <- data.frame(event_index = seq_len(100), DNA = seq_len(100), Target = 1)
  g1 <- complete
  g1$event_index[[1]] <- NA
  expect_error(
    facspseudocolor:::synchronized_validate_g1_containment(
      complete, g1, "SYNTHETIC sample", 100
    ), "Invalid, missing, or duplicate event_index"
  )
  expect_error(
    facspseudocolor:::synchronized_validate_g1_containment(
      complete[, c("DNA", "Target")], g1[, c("DNA", "Target")],
      "SYNTHETIC sample", 100
    ), "Required event identity/export columns are missing"
  )
})
