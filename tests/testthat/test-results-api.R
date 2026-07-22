small_analysis_object <- function(plot_type = "poi") {
  config <- validate_facs_config(minimal_config(plot_type))
  config$show_reference_panel <- TRUE
  config$y_limits <- c(500, 5000)
  manifest <- build_sample_manifest(config)
  make_sample <- function(offset) list(
    data = data.frame(
      dna_norm = seq(800, 2200, length.out = 60),
      target_norm = seq(700, 4000, length.out = 60) + offset,
      target_bgsub = seq(100, 1000, length.out = 60) + offset
    ),
    edu_positive = if (plot_type == "edu") data.frame(
      dna_norm = seq(1050, 1950, length.out = 30),
      target_norm = seq(2800, 4500, length.out = 30) + offset,
      target_bgsub = seq(500, 1500, length.out = 30) + offset
    ) else NULL,
    normalization_method = if (plot_type == "edu") {
      "reference_negative_regression"
    } else {
      "background_reference_regression"
    },
    cutoff = if (plot_type == "poi") 1200 else NA_real_
  )
  normalized <- list(reference = make_sample(0), treatment = make_sample(100))
  new_facs_analysis(
    config, manifest,
    data.frame(path = character(), exists = logical()),
    normalized, models = list(list())
  )
}

test_that("quantitation adds tables without altering normalized events", {
  analysis <- small_analysis_object("poi")
  original_names <- names(analysis$normalized_data$reference$data)
  result <- quantify_cell_cycle(analysis, reference_condition = "Reference")

  expect_s3_class(result, "facs_analysis")
  expect_true(all(c("phase_medians", "whole_medians", "phase_percentages") %in%
                    names(result$quantitation)))
  expect_true("ratio" %in% names(result$quantitation$phase_medians_reference))
  expect_identical(names(analysis$normalized_data$reference$data), original_names)
  expect_identical(analysis$quantitation, list())
})

test_that("plotting consumes analysis data and returns plot objects", {
  analysis <- small_analysis_object("poi")
  analysis$config$layout <- "cowplot"
  plots <- plot_pseudocolor_panels(analysis)

  expect_s3_class(plots, "facs_plot_set")
  expect_length(plots$decorated_plots, 2)
  expect_true(all(vapply(plots$decorated_plots, inherits, logical(1), "ggplot")))
  expect_s3_class(plots$grid, "ggplot")
})

test_that("quantitation plotting returns reusable ggplots", {
  analysis <- quantify_cell_cycle(
    small_analysis_object("poi"), reference_condition = "Reference"
  )
  plots <- plot_facs_quantitation(analysis)
  expect_true(all(c("phase_median", "whole_median", "phase_percent") %in%
                    names(plots)))
  expect_true(all(vapply(plots, inherits, logical(1), "ggplot")))
})

test_that("saving is explicit and refuses accidental overwrite", {
  analysis <- small_analysis_object("poi")
  analysis$config$layout <- "cowplot"
  plots <- plot_pseudocolor_panels(analysis)
  directory <- withr::local_tempdir()
  rds <- file.path(directory, "analysis.rds")

  written <- save_facs_results(
    analysis, plots, output_pdf = NULL, output_png = NULL,
    output_rds = rds
  )
  expect_true(file.exists(rds))
  expect_identical(unname(written), rds)
  expect_error(
    save_facs_results(analysis, plots, NULL, NULL, rds),
    "Refusing to overwrite"
  )
  expect_silent(save_facs_results(
    analysis, plots, NULL, NULL, rds, overwrite = TRUE
  ))
})
