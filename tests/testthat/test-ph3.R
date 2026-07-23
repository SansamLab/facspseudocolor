ph3_test_config <- function() {
  config <- minimal_config("ph3")
  config$g1_x_range <- c(750, 1250)
  config$s_phase_bins <- list(
    early = c(1250, 1450),
    mid = c(1450, 1650),
    late = c(1650, 1800)
  )
  config$g2m_x_range <- c(1800, 2400)
  validate_facs_config(config)
}

test_that("PH3 configuration requires no reference and explicit contiguous gates", {
  config <- ph3_test_config()
  expect_identical(config$plot_type, "ph3")
  expect_identical(
    names(config$suffixes), c("complete", "g1", "ph3_positive")
  )
  expect_false(any(build_sample_manifest(config)$is_reference))

  invalid <- minimal_config("ph3")
  expect_error(validate_facs_config(invalid), "explicit increasing DNA ranges")

  noncontiguous <- unclass(config)
  noncontiguous$s_phase_bins$early <- c(1300, 1450)
  expect_error(validate_facs_config(noncontiguous), "must be contiguous")

  invalid_sensitivity <- unclass(config)
  invalid_sensitivity$ph3_boundary_sensitivity_fraction <- 0
  expect_error(
    validate_facs_config(invalid_sensitivity),
    "ph3_boundary_sensitivity_fraction"
  )
})

test_that("PH3 normalization uses the G1 DNA gate and does not fit positivity", {
  complete <- data.frame(DNA = rep(c(100, 200), 10), PH3 = seq_len(20))
  g1 <- data.frame(DNA = c(95, 100, 100, 100, 105), PH3 = 1:5)
  positive <- data.frame(
    DNA = c(100, 130, 150, 170, 200, 260, NA_real_),
    PH3 = 101:107
  )
  result <- normalize_ph3(
    complete, g1, positive, "DNA", "PH3",
    dna_2n_value = 1000, g1_anchor = "median", sample_id = "test"
  )
  expect_equal(result$g1_dna_anchor, 100)
  expect_equal(result$dna_normalization_factor, 10)
  expect_equal(result$data$dna_norm, complete$DNA * 10)
  expect_identical(result$data$target_norm, complete$PH3)
  expect_equal(nrow(result$ph3_positive), 7)
  expect_identical(result$normalization_method, "g1_dna_only")
})

test_that("PH3 percentages always use all Single Cell events", {
  config <- ph3_test_config()
  manifest <- build_sample_manifest(config)
  complete <- data.frame(DNA = rep(c(100, 200), 10), PH3 = seq_len(20))
  g1 <- data.frame(DNA = c(95, 100, 100, 100, 105), PH3 = 1:5)
  positive <- data.frame(
    DNA = c(100, 130, 150, 170, 200, 260, NA_real_),
    PH3 = 101:107
  )
  normalized <- normalize_ph3(
    complete, g1, positive, "DNA", "PH3", sample_id = "Reference"
  )
  normalized2 <- normalize_ph3(
    complete, g1, positive[0, ], "DNA", "PH3", sample_id = "Treatment"
  )
  analysis <- new_facs_analysis(
    config, manifest,
    data.frame(path = character(), exists = logical()),
    list(reference = normalized, treatment = normalized2),
    models = list(list(), list())
  )
  result <- quantify_ph3(analysis)
  summary <- result$quantitation$ph3$sample_summary
  phases <- result$quantitation$ph3$phase_percentages

  expect_equal(summary$ph3_positive_percent, c(35, 0))
  expect_equal(summary$unassigned_percent, c(10, 0))
  expect_equal(
    phases$phase_percent[phases$prefix == "reference"],
    c(5, 5, 5, 5, 5, 10)
  )
  expect_equal(
    sum(phases$phase_percent[phases$prefix == "reference"]),
    summary$ph3_positive_percent[[1]]
  )
  expect_identical(result$quantitation$ph3$denominator,
                   "all_single_cell_events")
  sensitivity <- result$quantitation$ph3$boundary_sensitivity
  expect_equal(
    sensitivity$g2m_percent[
      sensitivity$prefix == "reference" &
        sensitivity$variant == "Configured"
    ],
    5
  )
})

test_that("empty PH3-positive CSVs are valid and counts cannot exceed denominator", {
  empty <- data.frame(DNA = numeric(), PH3 = numeric())
  expect_silent(read_facs_sample(empty, "DNA", "PH3", minimum_finite_events = 0))

  config <- ph3_test_config()
  directory <- withr::local_tempdir()
  for (prefix in c("reference", "treatment")) {
    write_test_csv(
      file.path(directory, paste0(prefix, "_single_cells.csv")),
      rep(c(100, 200), 5), seq_len(10)
    )
    write_test_csv(
      file.path(directory, paste0(prefix, "_g1.csv")),
      rep(100, 3), 1:3
    )
    utils::write.csv(
      data.frame(DNA = numeric(), Target = numeric()),
      file.path(directory, paste0(prefix, "_ph3_positive.csv")),
      row.names = FALSE
    )
  }
  expect_silent(validate_facs_inputs(config, directory))
})

test_that("PH3 input validation verifies event membership when indices exist", {
  config <- ph3_test_config()
  directory <- withr::local_tempdir()
  for (prefix in c("reference", "treatment")) {
    complete <- data.frame(
      event_index = 0:9,
      DNA = rep(c(100, 200), 5),
      Target = seq_len(10)
    )
    g1 <- complete[1:3, ]
    positive <- complete[c(2, 8), ]
    utils::write.csv(
      complete, file.path(directory, paste0(prefix, "_single_cells.csv")),
      row.names = FALSE
    )
    utils::write.csv(
      g1, file.path(directory, paste0(prefix, "_g1.csv")),
      row.names = FALSE
    )
    utils::write.csv(
      positive, file.path(directory, paste0(prefix, "_ph3_positive.csv")),
      row.names = FALSE
    )
  }
  report <- validate_facs_inputs(config, directory)
  expect_true(all(
    report$subset_membership_validated[
      report$population == "ph3_positive"
    ]
  ))

  invalid <- utils::read.csv(
    file.path(directory, "treatment_ph3_positive.csv"),
    check.names = FALSE
  )
  invalid$event_index[[1]] <- 999
  utils::write.csv(
    invalid, file.path(directory, "treatment_ph3_positive.csv"),
    row.names = FALSE
  )
  expect_error(validate_facs_inputs(config, directory), "outside")
})

test_that("PH3 plotting APIs return editable plots", {
  config <- ph3_test_config()
  manifest <- build_sample_manifest(config)
  complete <- data.frame(
    DNA = rep(seq(80, 220, length.out = 20), 2),
    PH3 = seq(100, 4000, length.out = 40)
  )
  g1 <- data.frame(DNA = seq(95, 105, length.out = 10), PH3 = 1:10)
  positive <- complete[c(5, 12, 18, 25, 32, 38), ]
  normalized <- lapply(c("Reference", "Treatment"), function(label) {
    normalize_ph3(complete, g1, positive, "DNA", "PH3", sample_id = label)
  })
  analysis <- quantify_ph3(new_facs_analysis(
    config, manifest, data.frame(path = character(), exists = logical()),
    stats::setNames(normalized, manifest$prefix), models = list(list(), list())
  ))
  expect_s3_class(plot_ph3_overall(analysis), "ggplot")
  expect_s3_class(plot_ph3_phase(analysis), "ggplot")
  expect_s3_class(plot_ph3_phase(analysis, "stacked"), "ggplot")
  expect_s3_class(plot_ph3_diagnostic(analysis), "ggplot")
  expect_s3_class(plot_ph3_boundary_sensitivity(analysis), "ggplot")
  plots <- plot_facs_quantitation(analysis)
  expect_true(all(c("ph3_overall", "phase_percent", "ph3_phase_stacked",
                    "ph3_diagnostic", "ph3_boundary_sensitivity") %in%
                  names(plots)))
  bundle <- build_facs_figure_bundle(analysis)
  expect_s3_class(bundle, "facs_figure_bundle")
  expect_true("ph3_overall" %in% names(bundle$quantitation))
})

test_that("high-level PH3 analysis runs from the three explicit CSV populations", {
  config <- ph3_test_config()
  directory <- withr::local_tempdir()
  for (prefix in c("reference", "treatment")) {
    write_test_csv(
      file.path(directory, paste0(prefix, "_single_cells.csv")),
      rep(seq(80, 220, length.out = 20), 2),
      seq(100, 4000, length.out = 40)
    )
    write_test_csv(
      file.path(directory, paste0(prefix, "_g1.csv")),
      seq(95, 105, length.out = 10), seq_len(10)
    )
    write_test_csv(
      file.path(directory, paste0(prefix, "_ph3_positive.csv")),
      c(100, 130, 150, 170, 200, 260),
      seq(1000, 6000, length.out = 6)
    )
  }
  analysis <- analyze_facs_experiment(config, data_dir = directory)
  expect_s3_class(analysis, "facs_analysis")
  expect_identical(analysis$config$plot_type, "ph3")
  expect_length(analysis$normalized_data, 2)
  expect_length(analysis$models, 2)
  expect_true("ph3" %in% names(analysis$quantitation))
  expect_equal(
    analysis$quantitation$ph3$sample_summary$ph3_positive_percent,
    c(15, 15)
  )
})
