# All events in this file are unmistakably SYNTHETIC and test-only.

synthetic_edu_display_analysis <- function() {
  config <- validate_facs_config(minimal_config("edu"))
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
        target_norm = raw, baseline = rep(baseline, length(raw))
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
  new_facs_analysis(
    config, manifest, data.frame(path = character(), exists = logical()),
    samples, models,
    quantitation = list(SYNTHETIC_quantitation = data.frame(value = c(1, 2)))
  )
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
