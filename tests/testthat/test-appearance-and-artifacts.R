artifact_test_analysis <- function() {
  config <- validate_facs_config(minimal_config("poi"))
  config$show_reference_panel <- TRUE
  config$y_limits <- c(500, 5000)
  manifest <- build_sample_manifest(config)
  make_sample <- function(offset) list(
    data = data.frame(
      dna_norm = seq(800, 2200, length.out = 60),
      target_norm = seq(700, 4000, length.out = 60) + offset,
      target_bgsub = seq(100, 1000, length.out = 60) + offset
    ),
    normalization_method = "background_reference_regression",
    cutoff = 1200
  )
  facspseudocolor:::new_facs_analysis(
    config, manifest, data.frame(path = character(), exists = logical()),
    list(reference = make_sample(0), treatment = make_sample(100)),
    models = list(list())
  )
}

test_that("appearance defaults are dynamic and null overrides are ignored", {
  analysis <- artifact_test_analysis()
  style <- resolve_facs_appearance(
    analysis,
    list(target_axis_label = NULL)
  )
  expect_s3_class(style, "facs_appearance")
  expect_match(style$target_axis_label, "background-normalized")
  expect_named(style$condition_colors, c("Reference", "Treatment"))
  expect_identical(unname(style$condition_colors), c("#440154", "#FDE725"))
})

test_that("partial condition colors override only known labels", {
  analysis <- artifact_test_analysis()
  style <- resolve_facs_appearance(
    analysis, list(condition_colors = list(Treatment = "#112233"))
  )
  expect_identical(unname(style$condition_colors[["Treatment"]]), "#112233")
  expect_true(nzchar(style$condition_colors[["Reference"]]))
  expect_error(
    resolve_facs_appearance(
      analysis, list(condition_colors = list(Misspelled = "#112233"))
    ),
    "Unknown condition color label"
  )
})

test_that("analysis and figure bundle artifacts round trip", {
  analysis <- quantify_cell_cycle(artifact_test_analysis())
  bundle <- build_facs_figure_bundle(
    analysis, appearance = list(phase_lineplot = FALSE)
  )
  expect_s3_class(bundle, "facs_figure_bundle")
  expect_s3_class(bundle$quantitation_page, "facs_panel_bundle")
  expect_equal(bundle$quantitation_page$layout$panel, 3.2)
  expect_true(all(vapply(bundle$quantitation, inherits, logical(1), "ggplot")))
  expect_false("normalized_data" %in% names(bundle))
  expect_true(all(vapply(bundle$pseudocolor$panel_results, function(x) {
    is.null(x$sample_data)
  }, logical(1))))
  expect_s3_class(get_facs_panel(bundle, "treatment"), "ggplot")

  directory <- withr::local_tempdir()
  analysis_path <- file.path(directory, "analysis.rds")
  bundle_path <- file.path(directory, "figures.rds")
  save_facs_analysis(analysis, analysis_path)
  save_facs_figure_bundle(bundle, bundle_path)
  expect_s3_class(read_facs_analysis(analysis_path), "facs_analysis")
  expect_s3_class(read_facs_figure_bundle(bundle_path), "facs_figure_bundle")
  expect_error(save_facs_analysis(analysis, analysis_path), "Refusing to overwrite")
})

test_that("focused plotting functions return individual editable ggplots", {
  analysis <- quantify_cell_cycle(
    artifact_test_analysis(), reference_condition = "Reference"
  )
  expect_s3_class(plot_facs_phase_signal(analysis), "ggplot")
  expect_s3_class(plot_facs_phase_signal(analysis, "line"), "ggplot")
  expect_s3_class(plot_facs_whole_signal(analysis), "ggplot")
  expect_s3_class(plot_facs_phase_percent(analysis), "ggplot")
  expect_s3_class(plot_facs_gate_assignments(analysis), "ggplot")
})

test_that("Quarto templates use exactly-one input setup", {
  helper <- system.file("quarto", "_report-setup.R", package = "facspseudocolor")
  expect_true(nzchar(helper))
  environment <- new.env(parent = globalenv())
  sys.source(helper, envir = environment)
  expect_error(environment$facs_report_input(), "exactly one")
  expect_error(environment$facs_report_input("a.yml", "a.rds"), "exactly one")
})
