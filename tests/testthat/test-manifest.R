test_that("a flat sample list produces a one-row manifest", {
  samples <- list(list(label = "Control", prefix = "control"), list(label = "Treatment", prefix = "treatment"))
  manifest <- make_sample_manifest(samples = samples)
  expect_equal(nrow(manifest), 2)
  expect_identical(manifest$condition, c("Control", "Treatment"))
  expect_identical(manifest$condition_index, 1:2)
  expect_false(any(manifest$is_reference))
})

test_that("replicate manifests identify exactly the named reference", {
  replicates <- lapply(1:2, function(i) list(
    label = paste("Replicate", i), reference = "Control",
    samples = list(list(label = "Control", prefix = paste0("r", i, "_control")),
                   list(label = "Treatment", prefix = paste0("r", i, "_treatment")))))
  manifest <- make_sample_manifest(replicates = replicates)
  expect_equal(nrow(manifest), 4)
  expect_identical(which(manifest$is_reference), c(1L, 3L))
  expect_identical(manifest$replicate_index, c(1L, 1L, 2L, 2L))
})

test_that("invalid sample layouts fail with actionable errors", {
  expect_error(make_sample_manifest(), "No samples")
  expect_error(make_sample_manifest(samples = list(list(label = "Control", prefix = ""))), "nonempty label and prefix")
  mismatched <- list(
    list(label = "R1", samples = list(list(label = "A", prefix = "r1a"), list(label = "B", prefix = "r1b"))),
    list(label = "R2", samples = list(list(label = "B", prefix = "r2b"), list(label = "A", prefix = "r2a"))))
  expect_error(make_sample_manifest(replicates = mismatched), "same conditions in the same order")
})

test_that("G1 anchors support median and mode with validation", {
  expect_equal(g1_target_anchor(c(1, 2, 100), "median"), 2)
  set.seed(2)
  expect_equal(g1_target_anchor(c(stats::rnorm(200, 5, 0.1), 10), "mode"), 5, tolerance = 0.1)
  expect_error(g1_target_anchor(1, "median"), "Too few")
  expect_error(g1_target_anchor(1:3, "mean"), "Unknown g1_anchor")
})
