# All analysis objects used here contain unmistakably SYNTHETIC test-only data.

edu_report_resource <- function(...) {
  relative <- file.path(...)
  installed <- system.file("quarto", relative, package = "facspseudocolor")
  if (nzchar(installed)) return(installed)
  source_path <- testthat::test_path("..", "..", "inst", "quarto", relative)
  if (file.exists(source_path)) return(source_path)
  stop("Required installed or source-tree report resource is unavailable: ",
       relative, call. = FALSE)
}

sys.source(
  edu_report_resource("_report-setup.R"),
  envir = environment()
)

synthetic_edu_report_analysis <- function() {
  raw_config <- minimal_config("edu")
  raw_config$replicates[[1L]]$reference <- "Reference"
  config <- validate_facs_config(raw_config)
  config$y_log10 <- FALSE
  config$x_limits <- c(700, 2250)
  manifest <- build_sample_manifest(config)
  make_sample <- function(raw_shift, baseline) {
    dna <- seq(800, 2200, length.out = 24L)
    raw <- seq(100, 2400, length.out = 24L) + raw_shift
    bgsub <- raw - baseline
    negative <- seq_len(12L)
    sample_data <- data.frame(
        dna_norm = dna, target_raw = raw, target_bgsub = bgsub,
        target_norm = raw, baseline = rep(baseline, length(raw)),
        edu_computed_positive = seq_along(raw) > 12L
    )
    list(
      sample = list(
        data = sample_data,
        edu_positive = sample_data[seq_along(raw) > 12L, , drop = FALSE]
      ),
      model = list(
        negative_event_index = negative,
        negative_event_fingerprint = edu_negative_event_fingerprint(
          dna[negative], raw[negative]
        )
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
  analysis <- new_facs_analysis(
    config, manifest, data.frame(
      replicate = character(), condition = character(), prefix = character(),
      population = character(), required = logical(), exists = logical(),
      event_n = integer(), finite_event_n = integer(),
      nonfinite_event_n = integer(), path = character()
    ), samples, models
  )
  suppressWarnings(quantify_cell_cycle(analysis))
}

test_that("standard EdU report preserves the reader-facing section contract", {
  path <- edu_report_resource("facs_edu_standard_v2.qmd")
  text <- readLines(path, warn = FALSE)
  headings <- trimws(text[grepl("^## ", text)])
  expect_identical(headings, c(
    "## All-condition visual overview", "## Quantitation overview",
    "## Experiment overview", "## Methods", "## Provenance"
  ))
  expect_false(any(grepl("output_csv_dir:", text, fixed = TRUE)))
  expect_false(any(grepl("show_quantitative_tables:", text, fixed = TRUE)))
  expect_true(any(grepl("show_audit_tables: false", text, fixed = TRUE)))
  expect_true(any(grepl("show_apex_comparison_toggle: true", text,
                        fixed = TRUE)))
  expect_true(any(grepl("show_phase_gate_toggle: true", text,
                        fixed = TRUE)))
  expect_true(any(grepl("embed_pseudocolor_pdf_downloads: true", text,
                        fixed = TRUE)))
  expect_true(any(grepl("download_plot, filename, width = 3, height = 3", text,
                        fixed = TRUE)))
  expect_true(any(grepl('id="apex-overview-toggle"', text, fixed = TRUE)))
  expect_true(any(grepl('id="phase-gate-toggle"', text, fixed = TRUE)))
  expect_true(any(grepl('"download-all-pdf-plain"', text, fixed = TRUE)))
  expect_true(any(grepl('"download-all-pdf-apex"', text, fixed = TRUE)))
  expect_true(any(grepl('"download-all-pdf-phase"', text, fixed = TRUE)))
  expect_true(any(grepl('"download-all-pdf-apex-phase"', text, fixed = TRUE)))
  expect_true(any(grepl('"all-pseudocolor-plots-cell-cycle-gates.pdf"',
                        text, fixed = TRUE)))
  expect_true(any(grepl('"all-pseudocolor-plots-apex-and-cell-cycle-gates.pdf"',
                        text, fixed = TRUE)))
  expect_true(any(grepl('label = "Download all PDF"', text, fixed = TRUE)))
  expect_true(any(grepl('apexButton.setAttribute("aria-pressed"', text,
                        fixed = TRUE)))
  expect_true(any(grepl("model_gate_manifest: null", text, fixed = TRUE)))
  expect_true(any(grepl("validate_edu_model_gate_manifest", text, fixed = TRUE)))
  expect_true(any(grepl("EXPERIMENTAL MODEL-DERIVED NON-PRODUCTION", text,
                        fixed = TRUE)))
  expect_true(any(grepl("share_model_gate_inputs$path <- basename", text,
                        fixed = TRUE)))
  expect_true(any(grepl("may be supplied only with model-derived gating", text,
                        fixed = TRUE)))
  validation_line <- grep("model_gate_inputs <- validate_edu_model_gate_manifest",
                          text, fixed = TRUE)
  report_line <- grep("report <- build_edu_pseudocolor_output_contract",
                      text, fixed = TRUE)
  expect_length(validation_line, 1L)
  expect_length(report_line, 1L)
  expect_lt(validation_line, report_line)
  expect_true(any(grepl("max_condition_columns: 4", text, fixed = TRUE)))
  expect_true(any(grepl("min_panel_width: 220", text, fixed = TRUE)))
  expect_true(any(grepl("width = 2.3, height = 2.15", text, fixed = TRUE)))
  expect_true(any(grepl("dpi = 192L", text, fixed = TRUE)))
  expect_true(any(grepl("repeat(auto-fit, minmax", text, fixed = TRUE)))
  expect_true(any(grepl("--panel-min", text, fixed = TRUE)))
  expect_true(any(grepl("facs_report_edu_responsive_groups", text, fixed = TRUE)))
  expect_true(any(grepl("--panel-min: %.0fpx; --panel-max: %.0fpx", text, fixed = TRUE)))
  expect_false(any(grepl("reference-derived EdU-negative", text, fixed = TRUE)))
  expect_true(any(grepl("independent EdU-negative baseline slope", text, fixed = TRUE)))
  expect_true(any(grepl("echo: false", text, fixed = TRUE)))
  expect_true(any(grepl("code-tools: false", text, fixed = TRUE)))
  expect_false(any(grepl("path = unname(exported_tables)", text, fixed = TRUE)))
  expect_false(any(grepl("theme: cosmo", text, fixed = TRUE)))
  expect_false(any(grepl("fonts.googleapis.com", text, fixed = TRUE)))
  expect_false(any(grepl("#| label: every-sample-panels", text,
                         fixed = TRUE)))
  removed_chunks <- c(
    "overview-balance", "overview-inputs", "qc-status", "qc-normalization",
    "qc-visual-review", "qc-software", "comparison-matrices",
    "replicate-variation", "quantitative-tables", "export-status"
  )
  expect_false(any(vapply(removed_chunks, function(label) {
    any(grepl(paste0("#| label: ", label), text, fixed = TRUE))
  }, logical(1))))
  expect_false(any(trimws(text) %in% c(
    "### Design completeness", "### Input populations and event counts",
    "## Quality control", "## Side-by-side comparisons",
    "## Replicate variation", "## Quantitative tables", "### Export"
  )))
  all_heading <- grep("^## All-condition visual overview$", text)
  quant_heading <- grep("^## Quantitation overview$", text)
  expect_identical(text[which(nzchar(trimws(text)) & seq_along(text) > all_heading)[[1L]]],
                   "```{r}")
  expect_true(any(grepl("These cards are the report's quantitative figure set", text,
                        fixed = TRUE)))
  quantitation_start <- grep("#\\| label: quantitation-overview", text)
  expect_length(quantitation_start, 1L)
  quantitation_end <- which(seq_along(text) > quantitation_start & text == "```")[[1L]]
  quantitation_text <- text[quantitation_start:quantitation_end]
  object_paths <- c(
    "report$positivity$panels$two_to_four_n",
    "report$intensity$panels$all_computed_positive",
    "report$positivity$panels$early_mid_late_s",
    "report$intensity$panels$early_mid_late_s"
  )
  observed_paths <- unlist(regmatches(
    quantitation_text,
    gregexpr("report\\$(positivity|intensity)\\$panels\\$[a-z0-9_]+",
             quantitation_text, perl = TRUE)
  ), use.names = FALSE)
  expect_identical(observed_paths[nzchar(observed_paths)], object_paths)
  quantitation_ids <- facs_report_card_ids(
    "quantitation-card", length(object_paths)
  )
  expect_length(quantitation_ids, 4L)
  expect_length(unique(quantitation_ids), 4L)
  expect_identical(quantitation_ids, paste0("quantitation-card-", 1:4))
  expect_true(any(grepl(
    'quantitation_card_ids <- facs_report_card_ids(',
    quantitation_text, fixed = TRUE
  )))
  expect_true(any(grepl(
    "caption_id <- quantitation_card_ids[[card_index]]",
    quantitation_text, fixed = TRUE
  )))
  expect_true(any(grepl(
    "aria-labelledby=\"', caption_id", quantitation_text, fixed = TRUE
  )))
  expect_true(any(grepl(
    "<h4 id=\"', caption_id", quantitation_text, fixed = TRUE
  )))
  expect_equal(sum(grepl("    alt = ", quantitation_text, fixed = TRUE)), 4L)
  expect_true(any(grepl("fig_alt = card$alt", quantitation_text, fixed = TRUE)))
  quant_plot_line <- grep("compact_plot <- card\\$plot", quantitation_text)
  quant_caption_line <- grep("facs_report_heading_label\\(card\\$title", quantitation_text)
  expect_length(quant_plot_line, 1L)
  expect_length(quant_caption_line, 1L)
  expect_lt(quant_plot_line, quant_caption_line)
  expect_true(any(grepl('legend.position = "bottom"', quantitation_text,
                        fixed = TRUE)))
  expect_true(any(grepl("title_size = 8, secondary_size = 7.5, legend_size = 7",
                        quantitation_text, fixed = TRUE)))
  expect_false(any(grepl("element_text(size = 7.5)", quantitation_text,
                         fixed = TRUE)))
  expect_equal(sum(grepl("ncol = quantitation_legend$columns",
                         quantitation_text, fixed = TRUE)), 2L)
  expect_true(any(grepl(
    "height = quantitation_legend$figure_height", quantitation_text,
    fixed = TRUE
  )))
  expect_true(any(grepl("for (row_index in seq_along(quantitation_rows))", quantitation_text,
                        fixed = TRUE)))
  expect_true(any(grepl("facs_report_quantitation_row_layout", quantitation_text,
                        fixed = TRUE)))
  expect_true(any(grepl("facs_report_quantitation_card_width", quantitation_text,
                        fixed = TRUE)))
  expect_true(any(grepl("facs_report_quantitation_plot_width(card_width)",
                        quantitation_text, fixed = TRUE)))
  expect_true(any(grepl("width = plot_width", quantitation_text, fixed = TRUE)))
  expect_true(any(grepl("quantitation_plot_widths", quantitation_text,
                        fixed = TRUE)))
  expect_true(any(grepl('"all-quantitation-plots.pdf"', quantitation_text,
                        fixed = TRUE)))
  expect_true(any(grepl('label = "Download all PDF"', quantitation_text,
                        fixed = TRUE)))
  expect_true(any(grepl('paste0(caption_id, ".pdf")', quantitation_text,
                        fixed = TRUE)))
  expect_true(any(grepl("Gating by sample", quantitation_text,
                        fixed = TRUE)))
  expect_true(any(grepl("facs_report_edu_gating_cards", quantitation_text,
                        fixed = TRUE)))
  expect_true(any(grepl('paste0(gating_card$prefix, "-gating.pdf")',
                        quantitation_text, fixed = TRUE)))
  expect_true(any(grepl("width = 9, height = 3.25", quantitation_text,
                        fixed = TRUE)))
  expect_true(any(grepl("Early-, Mid-, and Late-S populations by sample",
                        quantitation_text, fixed = TRUE)))
  expect_true(any(grepl("facs_report_edu_s_phase_cards", quantitation_text,
                        fixed = TRUE)))
  expect_true(any(grepl('"-s-phase-populations.pdf"', quantitation_text,
                        fixed = TRUE)))
  expect_false(any(grepl("#| label: every-sample-panels", text, fixed = TRUE)))
  expect_false(any(grepl("#| label: s-phase-positivity", text, fixed = TRUE)))
  expect_false(any(grepl("#| label: two-to-four-n-positivity", text, fixed = TRUE)))
  expect_false(any(grepl("#| label: positive-intensity-overall", text, fixed = TRUE)))
  expect_false(any(grepl("#| label: positive-intensity-regional", text, fixed = TRUE)))
  expect_true(any(grepl("quantitation_row_categories <- c(1L, 3L)",
                        quantitation_text, fixed = TRUE)))
  expect_true(any(grepl("analysis$sample_manifest$condition",
                        quantitation_text, fixed = TRUE)))
  expect_true(any(grepl(".overview-grid .quantitation-row", quantitation_text,
                        fixed = TRUE)))
  expect_true(any(grepl("axis.text.x = ggplot2::element_blank()",
                        quantitation_text, fixed = TRUE)))
  expect_true(any(grepl("axis.title.x = ggplot2::element_blank()",
                        quantitation_text, fixed = TRUE)))
  expect_equal(sum(grepl("hide_single_category_x_tick = TRUE",
                         quantitation_text, fixed = TRUE)), 2L)
  expect_equal(sum(grepl("hide_x_title = TRUE", quantitation_text,
                         fixed = TRUE)), 2L)
  expect_true(any(grepl('title = "2N–4N EdU-positive percentage"',
                        quantitation_text, fixed = TRUE)))
  expect_true(any(grepl('title = "All computed-positive EdU intensity"',
                        quantitation_text, fixed = TRUE)))
  expect_equal(sum(grepl('y_label = "% EdU positive"', quantitation_text,
                         fixed = TRUE)), 2L)
  expect_equal(sum(grepl('y_label = "Relative EdU level"', quantitation_text,
                         fixed = TRUE)), 2L)
  expect_false(any(grepl('y_label = "EdU Level"', quantitation_text,
                         fixed = TRUE)))
  expect_true(any(grepl("fill = NULL, colour = NULL", quantitation_text,
                        fixed = TRUE)))
  expect_true(any(grepl("legend.title = ggplot2::element_blank()",
                        quantitation_text, fixed = TRUE)))
  expect_equal(sum(grepl("title = NULL, ncol = quantitation_legend$columns",
                         quantitation_text, fixed = TRUE)), 2L)
  expect_true(any(grepl("caption = NULL, y = card$y_label",
                        quantitation_text, fixed = TRUE)))
  expect_equal(sum(grepl("small_legend_keys = TRUE", quantitation_text,
                         fixed = TRUE)), 2L)
  expect_true(any(grepl("facs_report_compact_legend_key_theme(0.28)",
                        quantitation_text, fixed = TRUE)))
  expect_true(any(grepl("Points are biological-replicate values",
                        quantitation_text, fixed = TRUE)))
  expect_true(any(grepl("Bar heights are condition means",
                        quantitation_text, fixed = TRUE)))
  expect_true(any(grepl("Horizontal marks are condition means",
                        quantitation_text, fixed = TRUE)))
  expect_true(any(grepl("Error bars use the configured",
                        quantitation_text, fixed = TRUE)))
  expect_true(any(grepl("Error bars are disabled by configuration",
                        quantitation_text, fixed = TRUE)))
  expect_true(any(grepl('"quantitation-card"', quantitation_text, fixed = TRUE)))
  pre_qc_line <- grep("#\\| label: pre-results-qc", text)
  all_condition_line <- grep("^## All-condition visual overview$", text)
  quantitation_line <- grep("^## Quantitation overview$", text)
  expect_length(pre_qc_line, 1L)
  expect_lt(pre_qc_line, all_condition_line)
  expect_lt(pre_qc_line, quantitation_line)
  expect_true(any(grepl("pre_results_qc$level != \"CLEAR\"", text, fixed = TRUE)))
  expect_false(any(grepl('href="#quality-control"', text, fixed = TRUE)))
  card_plot_line <- grep("card\\$plot, width = 2.3", text)
  card_caption_line <- grep("facs_report_heading_label\\(card_title", text)
  expect_length(card_plot_line, 1L)
  expect_length(card_caption_line, 1L)
  expect_lt(card_plot_line, card_caption_line)
  expect_true(any(grepl('role="group" aria-labelledby=', text, fixed = TRUE)))
})

test_that("legacy EdU report preserves its established parameter and section contract", {
  path <- edu_report_resource("facs_edu_pseudocolor_output_contract.qmd")
  text <- readLines(path, warn = FALSE)
  headings <- trimws(text[grepl("^## ", text)])
  expect_identical(headings, c(
    "## All-condition visual overview", "## Quantitation overview",
    "## Experiment overview", "## Every-sample DNA-versus-EdU plots",
    "## S-phase views", "## EdU-positive percentages",
    "## Background-subtracted EdU measurements", "## Methods",
    "## Provenance"
  ))
  expect_true(any(grepl("show_audit_tables: false", text, fixed = TRUE)))
  expect_false(any(grepl("show_apex_comparison_toggle", text, fixed = TRUE)))
  expect_false(any(grepl("show_phase_gate_toggle", text, fixed = TRUE)))
  expect_false(any(grepl("embed_pseudocolor_pdf_downloads", text, fixed = TRUE)))
})

test_that("compact quantitation legend space follows explicit condition count", {
  four <- facs_report_quantitation_legend_layout(4L)
  expect_identical(four$condition_count, 4L)
  expect_identical(four$columns, 2L)
  expect_identical(four$rows, 2L)
  expect_equal(four$figure_height, 3)

  seven <- facs_report_quantitation_legend_layout(7L)
  expect_identical(seven$condition_count, 7L)
  expect_identical(seven$columns, 2L)
  expect_identical(seven$rows, 4L)
  expect_equal(seven$figure_height, 3)

  expect_error(facs_report_quantitation_legend_layout(0), "positive integer")
  expect_error(facs_report_quantitation_legend_layout(2.5), "positive integer")

  two_card_row <- facs_report_quantitation_row_layout(2L, 220L)
  expect_identical(two_card_row$columns, 2L)
  expect_identical(two_card_row$min_panel_width, 220L)
  expect_equal(two_card_row$maximum_width, 456)
  expect_identical(facs_report_quantitation_row_layout(3L, 220L)$columns, 2L)
  expect_error(facs_report_quantitation_row_layout(0, 220), "positive integers")
})

test_that("compact quantitation card width follows conditions and categories", {
  expect_identical(facs_report_quantitation_card_width(4L, 1L), 220L)
  expect_identical(facs_report_quantitation_card_width(4L, 3L), 285L)

  expect_identical(facs_report_quantitation_card_width(1L, 1L), 220L)
  expect_identical(facs_report_quantitation_card_width(10L, 3L), 480L)

  expect_error(facs_report_quantitation_card_width(0, 1), "positive integers")
  expect_error(facs_report_quantitation_card_width(2.5, 1), "positive integers")
  expect_error(facs_report_quantitation_card_width(4, 0), "positive integers")
  expect_error(facs_report_quantitation_card_width(4, "3"), "positive integers")
  expect_error(facs_report_quantitation_card_width(Inf, 3), "positive integers")

  expect_equal(facs_report_quantitation_plot_width(220), 220 / 96)
  expect_equal(facs_report_quantitation_plot_width(285), 285 / 96)
  expect_error(facs_report_quantitation_plot_width(0), "finite positive")
  expect_error(facs_report_quantitation_plot_width(Inf), "finite positive")
  expect_error(facs_report_quantitation_plot_width("284"), "finite positive")
})

test_that("compact quantitation typography enforces readable reduced sizes", {
  typography <- facs_report_compact_typography(8, 7.5, 7)
  expect_equal(typography$axis.title$size, 8)
  expect_equal(typography$axis.text$size, 7.5)
  expect_equal(typography$legend.title$size, 8)
  expect_equal(typography$legend.text$size, 7)
  expect_equal(typography$strip.text$size, 7.5)
  expect_error(facs_report_compact_typography(7.9, 7.5), "at least 8 pt")
  expect_error(facs_report_compact_typography(8, 7.4), "at least 7.5 pt")
  expect_error(facs_report_compact_typography(8, NA_real_), "at least 8 pt")
  expect_error(facs_report_compact_typography(8, 7.5, 6.9), "at least 7 pt")

  intensity_plot_body <- paste(deparse(body(edu_intensity_plot)), collapse = "\n")
  expect_match(intensity_plot_body, "theme_classic\\(base_size = 10\\)")
})

test_that("compact legend key theme validates and reduces overview keys", {
  key_theme <- facs_report_compact_legend_key_theme(0.28)
  expect_equal(as.numeric(key_theme$legend.key.size), 0.28)
  expect_equal(as.numeric(key_theme$legend.spacing.x), 0.08)
  expect_error(facs_report_compact_legend_key_theme(0), "positive finite")
  expect_error(facs_report_compact_legend_key_theme(NA_real_), "positive finite")
})

test_that("compact caption removal leaves detailed quantitative plots unchanged", {
  analysis <- synthetic_edu_report_analysis()
  report <- build_edu_pseudocolor_output_contract(analysis)
  sources <- list(
    report$positivity$panels$two_to_four_n,
    report$intensity$panels$all_computed_positive,
    report$positivity$panels$early_mid_late_s,
    report$intensity$panels$early_mid_late_s
  )
  source_captions <- vapply(
    sources, function(plot) plot$labels$caption, character(1)
  )

  compact_copies <- lapply(sources, function(plot) {
    plot + ggplot2::labs(caption = NULL)
  })

  expect_true(all(vapply(
    compact_copies, function(plot) is.null(plot$labels$caption), logical(1)
  )))
  expect_identical(
    vapply(sources, function(plot) plot$labels$caption, character(1)),
    source_captions
  )
  expect_true(all(nzchar(source_captions)))
})

test_that("compact condition guides hide titles but retain keys and source scales", {
  analysis <- synthetic_edu_report_analysis()
  report <- build_edu_pseudocolor_output_contract(analysis)
  cases <- list(
    list(plot = report$positivity$panels$two_to_four_n, aesthetic = "fill"),
    list(plot = report$intensity$panels$all_computed_positive,
         aesthetic = "colour")
  )
  for (case in cases) {
    source_built <- ggplot2::ggplot_build(case$plot)
    source_scale <- source_built$plot$scales$get_scales(case$aesthetic)
    source_labels <- source_scale$get_labels()
    source_limits <- source_scale$get_limits()
    source_values <- source_scale$map(source_limits)
    source_name <- source_scale$name
    guide <- ggplot2::guide_legend(title = NULL, ncol = 2, byrow = TRUE)
    guides <- list(guide)
    names(guides) <- case$aesthetic
    compact <- case$plot +
      ggplot2::theme(legend.title = ggplot2::element_blank()) +
      do.call(ggplot2::guides, guides)
    compact_built <- ggplot2::ggplot_build(compact)
    compact_scale <- compact_built$plot$scales$get_scales(case$aesthetic)

    expect_s3_class(compact_built$plot$theme$legend.title, "element_blank")
    expect_identical(compact_scale$get_labels(), source_labels)
    expect_identical(compact_scale$get_limits(), source_limits)
    expect_identical(compact_scale$map(compact_scale$get_limits()),
                     source_values)
    expect_identical(case$plot$scales$get_scales(case$aesthetic)$name,
                     source_name)
    expect_gt(length(source_labels), 0L)
  }
})

test_that("responsive overview uses explicit model groups and reference-first cards", {
  analysis <- synthetic_edu_report_analysis()
  report <- build_edu_pseudocolor_output_contract(analysis)
  second <- analysis$sample_manifest
  second$model_group <- "1::SYNTHETIC_T2"
  second$technical_replicate <- "SYNTHETIC_T2"
  second$prefix <- paste0(second$prefix, "_SYNTHETIC_T2")
  analysis$sample_manifest <- rbind(analysis$sample_manifest, second)
  report$panels <- c(report$panels, stats::setNames(report$panels, second$prefix))
  report$panel_qc$y_limit_lower <- c(10, 20)
  report$panel_qc$y_limit_upper <- c(100, 200)
  second_qc <- report$panel_qc
  second_qc$prefix <- second$prefix
  second_qc$y_limit_lower <- c(1000, 2000)
  second_qc$y_limit_upper <- c(3000, 4000)
  report$panel_qc <- rbind(report$panel_qc, second_qc)
  original_audit_y <- vapply(report$panels, function(panel) panel$labels$y, character(1))

  overview <- facs_report_edu_responsive_groups(analysis, report, 2L, 300L)
  expect_identical(vapply(overview$groups, `[[`, character(1), "model_group"),
                   c("1::1", "1::SYNTHETIC_T2"))
  expect_true(all(vapply(overview$groups, function(group) {
    isTRUE(group$cards[[1L]]$is_reference)
  }, logical(1))))
  represented <- unlist(lapply(overview$groups, function(group) {
    vapply(group$cards, `[[`, character(1), "prefix")
  }), use.names = FALSE)
  expect_setequal(represented[!is.na(represented)], analysis$sample_manifest$prefix)
  expect_identical(overview$max_condition_columns, 2L)
  expect_identical(overview$min_panel_width, 300L)
  expect_identical(overview$groups[[1L]]$shared_y_limits, c(10, 200))
  expect_identical(overview$groups[[2L]]$shared_y_limits, c(1000, 4000))
  expect_identical(vapply(report$panels, function(panel) panel$labels$y, character(1)),
                   original_audit_y)
  expect_true(all(vapply(overview$groups, function(group) {
    all(vapply(group$cards, function(card) {
      is.null(card$plot$labels$title) && is.null(card$plot$labels$subtitle)
    }, logical(1)))
  }, logical(1))))
  expect_true(all(vapply(report$panels, function(panel) {
    !is.null(panel$labels$title) && !is.null(panel$labels$subtitle)
  }, logical(1))))
})

test_that("responsive overview shares the union of retained y ranges within each group", {
  analysis <- synthetic_edu_report_analysis()
  report <- build_edu_pseudocolor_output_contract(analysis)
  report$panel_qc$y_limit_lower <- c(10, 20)
  report$panel_qc$y_limit_upper <- c(100, 400)
  original_y_limits <- lapply(
    report$panels, function(panel) panel$coordinates$limits$y
  )
  scale_signature <- function(panel) lapply(panel$scales$scales, function(scale) {
    list(aesthetics = scale$aesthetics, transformation = scale$trans$name,
         breaks = scale$breaks)
  })
  original_scale_signatures <- lapply(report$panels, scale_signature)

  overview <- facs_report_edu_responsive_groups(analysis, report)

  expect_identical(overview$groups[[1L]]$shared_y_limits, c(10, 400))
  expect_true(all(vapply(overview$groups[[1L]]$cards, function(card) {
    identical(card$plot$coordinates$limits$y, c(10, 400))
  }, logical(1))))
  expect_identical(
    lapply(report$panels, function(panel) panel$coordinates$limits$y),
    original_y_limits
  )
  expect_identical(lapply(report$panels, scale_signature),
                   original_scale_signatures)
  for (card in overview$groups[[1L]]$cards) {
    source <- report$panels[[card$prefix]]
    expect_identical(scale_signature(card$plot), scale_signature(source))
  }
})

test_that("responsive overview rejects invalid available-panel y-limit pairs", {
  analysis <- synthetic_edu_report_analysis()

  reversed <- build_edu_pseudocolor_output_contract(analysis)
  reversed$panel_qc$y_limit_lower[[1L]] <- 100
  reversed$panel_qc$y_limit_upper[[1L]] <- 10
  expect_error(
    facs_report_edu_responsive_groups(analysis, reversed),
    "strictly below its upper limit"
  )

  equal <- build_edu_pseudocolor_output_contract(analysis)
  equal$panel_qc$y_limit_lower[[1L]] <- 100
  equal$panel_qc$y_limit_upper[[1L]] <- 100
  expect_error(
    facs_report_edu_responsive_groups(analysis, equal),
    "strictly below its upper limit"
  )

  nonnumeric <- build_edu_pseudocolor_output_contract(analysis)
  nonnumeric$panel_qc$y_limit_lower <- as.character(
    nonnumeric$panel_qc$y_limit_lower
  )
  expect_error(
    facs_report_edu_responsive_groups(analysis, nonnumeric),
    "numeric retained lower and upper"
  )
})

test_that("responsive overview validates width and normalization identity", {
  analysis <- synthetic_edu_report_analysis()
  report <- build_edu_pseudocolor_output_contract(analysis)
  expect_error(facs_report_edu_responsive_groups(analysis, report, 3L, 0),
               "min_panel_width")
  expect_error(facs_report_edu_responsive_groups(
    analysis, report, 3L, as.double(.Machine$integer.max) + 1
  ), "min_panel_width")
  analysis$sample_manifest$model_group[[1L]] <- ""
  expect_error(facs_report_edu_responsive_groups(analysis, report),
               "normalization identities")
})

test_that("responsive paired groups require consistent reference designation", {
  analysis <- synthetic_edu_report_analysis()
  report <- build_edu_pseudocolor_output_contract(analysis)
  analysis$sample_manifest$is_reference[[1L]] <- NA
  expect_error(facs_report_edu_responsive_groups(analysis, report),
               "nonmissing is_reference")

  analysis <- synthetic_edu_report_analysis()
  report <- build_edu_pseudocolor_output_contract(analysis)
  analysis$sample_manifest$reference_condition[[2L]] <- "Treatment"
  expect_error(facs_report_edu_responsive_groups(analysis, report),
               "same nonempty reference_condition")

  analysis <- synthetic_edu_report_analysis()
  report <- build_edu_pseudocolor_output_contract(analysis)
  analysis$sample_manifest$is_reference <- rev(analysis$sample_manifest$is_reference)
  expect_error(facs_report_edu_responsive_groups(analysis, report),
               "same row marked is_reference")
})

test_that("responsive group legend fallback retains every card", {
  analysis <- synthetic_edu_report_analysis()
  report <- build_edu_pseudocolor_output_contract(analysis)
  report$panel_qc$panel_status[[1L]] <- "suppressed"
  report$panel_qc$y_limit_lower <- c(-1000, 20)
  report$panel_qc$y_limit_upper <- c(1000, 40)
  overview <- facs_report_edu_responsive_groups(analysis, report)
  expect_identical(overview$max_condition_columns, 4L)
  expect_identical(overview$min_panel_width, 220L)
  expect_identical(overview$groups[[1L]]$legend_source_prefix,
                   analysis$sample_manifest$prefix[[2L]])
  expect_identical(overview$groups[[1L]]$shared_legend_count, 1L)
  expect_length(overview$groups[[1L]]$cards, nrow(analysis$sample_manifest))
  expect_identical(overview$groups[[1L]]$shared_y_limits, c(20, 40))
  suppressed_card <- overview$groups[[1L]]$cards[
    vapply(overview$groups[[1L]]$cards, `[[`, character(1), "prefix") ==
      analysis$sample_manifest$prefix[[1L]]
  ][[1L]]
  expect_identical(suppressed_card$status, "suppressed")
  expect_false(identical(suppressed_card$plot$coordinates$limits$y, c(20, 40)))

  report$panel_qc$panel_status[] <- "suppressed"
  suppressed <- facs_report_edu_responsive_groups(analysis, report)
  expect_identical(suppressed$groups[[1L]]$legend_status,
                   "unavailable_all_group_panels_suppressed")
  expect_identical(suppressed$groups[[1L]]$shared_legend_count, 0L)
  expect_length(suppressed$groups[[1L]]$cards, nrow(analysis$sample_manifest))
})

test_that("QC separates scientific checks, visual review, and software notices", {
  analysis <- synthetic_edu_report_analysis()
  report <- build_edu_pseudocolor_output_contract(analysis)
  analysis$warnings <- c(
    paste0(
      "EdU aliases `phase_percentages`, `phase_medians`, and `whole_medians` ",
      "are deprecated for one major-release compatibility window; use the ",
      "canonical `edu_*` tables. Alias meanings are unchanged."
    ),
    "SYNTHETIC retained data warning",
    "SYNTHETIC software version mismatch requiring scientific review"
  )
  sections <- facs_report_edu_qc_sections(analysis, report)
  expect_true(any(grepl("deprecated", sections$software$finding, ignore.case = TRUE)))
  expect_false(any(grepl("deprecated", sections$data_input$finding, ignore.case = TRUE)))
  expect_identical(sections$visual_review$level, "REQUIRED")
  expect_true(any(grepl("SYNTHETIC retained data warning", sections$data_input$finding, fixed = TRUE)))
  expect_true(any(grepl("software version mismatch", sections$data_input$finding, fixed = TRUE)))
  expect_false(any(grepl("software version mismatch", sections$software$finding, fixed = TRUE)))
})

test_that("all-condition overview retains manifest identities and missing cells", {
  analysis <- synthetic_edu_report_analysis()
  report <- build_edu_pseudocolor_output_contract(analysis)
  overview <- facs_report_edu_all_conditions(analysis, report)
  expect_length(overview$blocks, 1L)
  expect_s3_class(overview$blocks[[1L]]$plot, "ggplot")
  expect_gt(overview$blocks[[1L]]$figure_width, 0)
  expect_gt(overview$blocks[[1L]]$figure_height, 0)
  expect_identical(overview$blocks[[1L]]$shared_legend_count, 1L)
  expect_null(overview$blocks[[1L]]$shared_legend_title)
  expect_false(overview$blocks[[1L]]$panel_legends_visible)
  expect_identical(overview$blocks[[1L]]$y_axis_title,
                   "EdU")
  expect_identical(overview$blocks[[1L]]$row_header_position,
                   "above_full_width")
  expect_true(overview$blocks[[1L]]$replicate_header_height_fixed)
  expect_identical(overview$blocks[[1L]]$replicate_header_rel_height, .12)
  expect_true(all(overview$blocks[[1L]]$replicate_header_rel_heights == .12))
  expect_identical(overview$blocks[[1L]]$replicate_labels,
                   unique(as.character(analysis$sample_manifest$replicate)))
  expect_identical(report$panels[[1L]]$labels$x, "DNA")
  expect_identical(report$panels[[1L]]$labels$y, "EdU")
  expect_false(grepl("Condition: ", report$panels[[1L]]$labels$subtitle,
                     fixed = TRUE))
  density_scale <- report$panels[[1L]]$scales$get_scales("colour")
  expect_null(density_scale$name)
  expect_identical(density_scale$breaks, c(0, 0.5, 1))
  expect_identical(density_scale$labels, c("0", "0.5", "1"))
  guide_theme <- density_scale$guide$params$theme
  expect_equal(as.numeric(guide_theme[["legend.key.height"]]), 0.6)
  expect_equal(as.numeric(guide_theme[["legend.key.width"]]), 0.153)
  expect_identical(sum(overview$balance$`Technical acquisitions`), nrow(analysis$sample_manifest))
  expect_true(all(overview$balance$Status == "Present"))
  reference <- analysis$sample_manifest[analysis$sample_manifest$is_reference, , drop = FALSE]
  reference$replicate <- "SYNTHETIC replicate with missing treatment"
  reference$replicate_index <- max(analysis$sample_manifest$replicate_index) + 1L
  reference$prefix <- "SYNTHETIC_missing_treatment_reference"
  analysis$sample_manifest <- rbind(analysis$sample_manifest, reference)
  report$panels[[reference$prefix]] <- report$panels[[1L]]
  report$panel_qc <- rbind(report$panel_qc, transform(
    report$panel_qc[1L, , drop = FALSE], prefix = reference$prefix,
    panel_status = "available", panel_reason_code = NA_character_
  ))
  overview_missing <- facs_report_edu_all_conditions(analysis, report)
  expect_true(any(overview_missing$balance$Status == "MISSING"))
})

test_that("replicate header allocation is fixed when acquisition stacks differ", {
  analysis <- synthetic_edu_report_analysis()
  report <- build_edu_pseudocolor_output_contract(analysis)
  technical <- analysis$sample_manifest[1L, , drop = FALSE]
  technical$technical_replicate <- "2"
  technical$prefix <- "SYNTHETIC_extra_technical_acquisition"
  analysis$sample_manifest <- rbind(analysis$sample_manifest, technical)
  report$panels[[technical$prefix]] <- report$panels[[1L]]
  report$panel_qc <- rbind(report$panel_qc, transform(
    report$panel_qc[1L, , drop = FALSE], prefix = technical$prefix
  ))
  overview <- facs_report_edu_all_conditions(analysis, report)
  expect_true(all(overview$blocks[[1L]]$replicate_header_rel_heights == .12))
  expect_gt(overview$blocks[[1L]]$replicate_body_rel_heights[[1L]], 1L)
})

test_that("all-condition overview validates wrapping and retains condition order", {
  analysis <- synthetic_edu_report_analysis()
  report <- build_edu_pseudocolor_output_contract(analysis)
  template <- analysis$sample_manifest[1L, , drop = FALSE]
  conditions <- paste("SYNTHETIC condition", 1:5)
  manifest <- do.call(rbind, lapply(seq_along(conditions), function(index) {
    row <- template
    row$condition <- conditions[[index]]
    row$condition_index <- index
    row$prefix <- paste0("SYNTHETIC_condition_", index)
    row$technical_replicate <- "1"
    row
  }))
  analysis$sample_manifest <- manifest
  report$panels <- stats::setNames(
    rep(report$panels[1L], length(conditions)), manifest$prefix
  )
  panel_qc_template <- report$panel_qc[1L, , drop = FALSE]
  report$panel_qc <- do.call(rbind, lapply(
    as.character(analysis$sample_manifest$prefix), function(prefix_value) {
      row <- panel_qc_template
      row$prefix <- as.character(prefix_value)
      row$panel_status <- "available"
      row$panel_reason_code <- NA_character_
      row
    }
  ))
  rownames(report$panel_qc) <- NULL
  expect_identical(
    as.character(report$panel_qc$prefix),
    as.character(analysis$sample_manifest$prefix)
  )
  overview <- facs_report_edu_all_conditions(analysis, report, 2L)
  expect_length(overview$blocks, 3L)
  expect_identical(unlist(lapply(overview$blocks, `[[`, "conditions"), use.names = FALSE),
                   conditions)
  expect_identical(sum(vapply(overview$blocks, function(block) length(block$conditions), integer(1))),
                   length(conditions))
  expect_true(all(vapply(overview$blocks, `[[`, integer(1), "shared_legend_count") == 1L))
  expect_true(all(lengths(lapply(overview$blocks, `[[`, "conditions")) <= 2L))
  expect_error(facs_report_edu_all_conditions(analysis, report, 0),
               "one positive integer")
  expect_error(facs_report_edu_all_conditions(analysis, report, 2.5),
               "one positive integer")
  expect_error(facs_report_edu_all_conditions(analysis, report, c(2, 3)),
               "one positive integer")
  expect_error(facs_report_edu_all_conditions(
    analysis, report, as.double(.Machine$integer.max) + 1
  ), "one positive integer")
})

test_that("all-condition legend uses the first available panel", {
  analysis <- synthetic_edu_report_analysis()
  report <- build_edu_pseudocolor_output_contract(analysis)
  report$panel_qc$panel_status[[1L]] <- "suppressed"
  report$panel_qc$panel_reason_code[[1L]] <- "SYNTHETIC_suppressed"
  overview <- facs_report_edu_all_conditions(analysis, report)
  expect_identical(overview$blocks[[1L]]$legend_source_prefix,
                   analysis$sample_manifest$prefix[[2L]])
  expect_identical(overview$blocks[[1L]]$shared_legend_count, 1L)
  expect_false(overview$all_panels_suppressed)
  expect_identical(sum(overview$balance$`Technical acquisitions`),
                   nrow(analysis$sample_manifest))
})

test_that("all-condition overview reports all panels suppressed without stopping", {
  analysis <- synthetic_edu_report_analysis()
  report <- build_edu_pseudocolor_output_contract(analysis)
  report$panel_qc$panel_status[] <- "suppressed"
  report$panel_qc$panel_reason_code[] <- "SYNTHETIC_suppressed"
  overview <- facs_report_edu_all_conditions(analysis, report)
  expect_true(overview$all_panels_suppressed)
  expect_identical(overview$density_legend_status,
                   "unavailable_all_panels_suppressed")
  expect_true(all(vapply(overview$blocks, `[[`, integer(1),
                         "shared_legend_count") == 0L))
  expect_true(all(is.na(vapply(overview$blocks, `[[`, character(1),
                               "legend_source_prefix"))))
  expect_identical(sum(overview$balance$`Technical acquisitions`),
                   nrow(analysis$sample_manifest))
})

test_that("all-condition overview rejects ambiguous manifest identities", {
  analysis <- synthetic_edu_report_analysis()
  report <- build_edu_pseudocolor_output_contract(analysis)
  analysis$sample_manifest$prefix[[2L]] <- analysis$sample_manifest$prefix[[1L]]
  expect_error(facs_report_edu_all_conditions(analysis, report), "exact manifest-prefix order")

  analysis <- synthetic_edu_report_analysis()
  report <- build_edu_pseudocolor_output_contract(analysis)
  analysis$sample_manifest$condition[[1L]] <- ""
  expect_error(facs_report_edu_all_conditions(analysis, report), "nonmissing condition")

  analysis <- synthetic_edu_report_analysis()
  report <- build_edu_pseudocolor_output_contract(analysis)
  analysis$sample_manifest$replicate_index[[2L]] <- analysis$sample_manifest$replicate_index[[1L]]
  analysis$sample_manifest$replicate[[2L]] <- "SYNTHETIC conflicting label"
  expect_error(facs_report_edu_all_conditions(analysis, report), "exactly one biological-replicate label")
})

test_that("report-generated headings escape markup without changing labels", {
  expect_identical(
    facs_report_heading_label("Condition <A> & prefix '1'"),
    "Condition &lt;A&gt; &amp; prefix &#39;1&#39;"
  )
  expect_identical(
    facs_report_heading_label("# condition\n[link](bad)   *emphasis*"),
    "# condition [link](bad) *emphasis*"
  )
})

test_that("report plots use Knitr child chunks with explicit dimensions", {
  helper <- readLines(
    edu_report_resource("_report-setup.R"),
    warn = FALSE
  )
  qmd <- readLines(
    edu_report_resource("facs_edu_standard_v2.qmd"),
    warn = FALSE
  )
  expect_true(any(grepl("knitr::knit_child", helper, fixed = TRUE)))
  expect_true(any(grepl("#| fig-width:", helper, fixed = TRUE)))
  expect_true(any(grepl("#| fig-height:", helper, fixed = TRUE)))
  expect_true(any(grepl("#| fig-dpi:", helper, fixed = TRUE)))
  expect_true(any(grepl("#| fig-alt:", helper, fixed = TRUE)))
  expect_true(any(grepl("facs_report_knit_plot", qmd, fixed = TRUE)))
  expect_false(any(grepl("knitr::knit_print", qmd, fixed = TRUE)))
  expect_false(any(grepl("opts_current", qmd, fixed = TRUE)))
  analysis <- synthetic_edu_report_analysis()
  report <- build_edu_pseudocolor_output_contract(analysis)
  expect_error(
    facs_report_knit_plot(report$panels[[1L]], 2.3, 2.15, dpi = 0),
    "DPI must be one positive integer"
  )
  expect_error(
    facs_report_knit_plot(report$panels[[1L]], 2.3, 2.15, fig_alt = ""),
    "fig_alt"
  )
  withr::local_dir(withr::local_tempdir(pattern = "SYNTHETIC_grob_render_"))
  rendered_grob <- facs_report_knit_plot(
    grid::textGrob("SYNTHETIC legend"), 2, 0.7
  )
  expect_false(any(grepl("TableGrob", rendered_grob, fixed = TRUE)))
})

test_that("EdU report overview exposes configured identities and input provenance", {
  analysis <- synthetic_edu_report_analysis()
  overview <- facs_report_edu_overview(analysis)

  expect_identical(
    overview$samples$prefix,
    as.character(analysis$sample_manifest$prefix)
  )
  expect_true(all(c("FlowJo Single Cells", "FlowJo G1", "FlowJo EdU Positive") %in%
                    strsplit(overview$summary$value[overview$summary$item == "Population sources"],
                             ", ", fixed = TRUE)[[1L]]))
  expect_true(all(c("prefix", "population", "path") %in% names(overview$inputs)))

  model_analysis <- analysis
  model_analysis$config$gating <- list(mode = "model_experimental")
  model_overview <- facs_report_edu_overview(model_analysis)
  expect_true(all(c("Model-derived Single Cells", "Model-derived G1",
                    "Model-derived EdU Positive") %in%
                  strsplit(model_overview$summary$value[
                    model_overview$summary$item == "Population sources"
                  ], ", ", fixed = TRUE)[[1L]]))
})

test_that("report helpers resolve package internals across the sourced-file boundary", {
  analysis <- synthetic_edu_report_analysis()
  report <- build_edu_pseudocolor_output_contract(analysis)
  helper_environment <- new.env(parent = baseenv())
  sys.source(
    edu_report_resource("_report-setup.R"),
    envir = helper_environment
  )

  overview <- helper_environment$facs_report_edu_overview(analysis)
  qc <- helper_environment$facs_report_edu_qc(analysis, report)

  expect_identical(overview$samples$prefix, analysis$sample_manifest$prefix)
  expect_true(is.data.frame(qc))
  expect_true(all(c("level", "area", "sample", "finding") %in% names(qc)))
})

test_that("EdU report QC makes retained failures visible without changing analysis", {
  analysis <- synthetic_edu_report_analysis()
  report <- build_edu_pseudocolor_output_contract(analysis)
  analysis$warnings <- "SYNTHETIC retained analysis warning"
  report$panel_qc$panel_status[[1L]] <- "suppressed"
  report$panel_qc$panel_reason_code[[1L]] <- "SYNTHETIC_reason"
  before <- serialize(analysis, NULL)

  qc <- facs_report_edu_qc(analysis, report)

  expect_true(any(qc$level == "REVIEW" & qc$area == "Analysis"))
  expect_true(any(qc$level == "BLOCKED" & qc$area == "DNA-versus-EdU panel"))
  expect_identical(serialize(analysis, NULL), before)
})

test_that("apex comparison retains exact references and distinct line roles", {
  analysis <- synthetic_edu_report_analysis()
  analysis$config$edu_apex_x_range <- c(800, 2200)
  report <- build_edu_pseudocolor_output_contract(analysis)
  comparison <- facs_report_edu_apex_comparison(analysis, report)

  expect_identical(comparison$apex$prefix,
                   as.character(analysis$sample_manifest$prefix))
  expect_true(all(comparison$apex$apex_event_n >= 10L))
  expect_identical(unique(comparison$apex$reference_prefix), "reference")
  expect_equal(comparison$apex$current_line_display_signal,
               comparison$apex$apex_background_subtracted_signal +
                 comparison$apex$display_offset)
  expect_equal(comparison$apex$reference_line_display_signal,
               comparison$apex$reference_apex_background_subtracted_signal +
                 comparison$apex$display_offset)
  is_hline <- function(layer) inherits(layer$geom, "GeomHline")
  reference_lines <- Filter(is_hline, comparison$report$panels$reference$layers)
  treatment_lines <- Filter(is_hline, comparison$report$panels$treatment$layers)
  expect_identical(length(reference_lines), 1L)
  expect_identical(length(treatment_lines), 2L)
  expect_identical(reference_lines[[1L]]$aes_params$colour, "#2166AC")
  expect_identical(reference_lines[[1L]]$aes_params$linetype, "dashed")
  expect_identical(treatment_lines[[2L]]$aes_params$colour, "#B2182B")
  expect_identical(treatment_lines[[2L]]$aes_params$linetype, "solid")
  expect_identical(comparison$reference_colour, "#2166AC")
  expect_identical(comparison$sample_colour, "#B2182B")
})

test_that("phase-gate comparison uses all five established EdU gates", {
  analysis <- synthetic_edu_report_analysis()
  analysis$normalized_data <- lapply(analysis$normalized_data, function(sample) {
    sample$data$edu_boundary_bgsub <- rep(2500, nrow(sample$data))
    sample
  })
  report <- build_edu_pseudocolor_output_contract(analysis)
  comparison <- facs_report_edu_phase_gate_comparison(analysis, report)

  expect_identical(names(comparison$report$panels),
                   as.character(analysis$sample_manifest$prefix))
  expect_identical(names(comparison$gates),
                   as.character(analysis$sample_manifest$prefix))
  expect_true(all(vapply(
    comparison$gates,
    function(x) identical(as.character(x$gate),
                           c("G1", "Early S", "Mid S", "Late S", "G2/M")),
    logical(1)
  )))
  expect_true(all(vapply(comparison$report$panels, function(plot) {
    any(vapply(plot$layers, function(layer) inherits(layer$geom, "GeomPolygon"),
               logical(1))) &&
      any(vapply(plot$layers, function(layer) inherits(layer$geom, "GeomText"),
                 logical(1)))
  }, logical(1))))
  expect_identical(comparison$line_type, "dashed")
})

test_that("apex comparison fails closed on malformed reference identities", {
  analysis <- synthetic_edu_report_analysis()
  analysis$config$edu_apex_x_range <- c(800, 2200)
  report <- build_edu_pseudocolor_output_contract(analysis)

  malformed <- analysis
  malformed$sample_manifest$model_group[[1L]] <- ""
  expect_error(facs_report_edu_apex_comparison(malformed, report),
               "nonempty character identities")
  malformed <- analysis
  malformed$sample_manifest$is_reference <- as.integer(
    malformed$sample_manifest$is_reference
  )
  expect_error(facs_report_edu_apex_comparison(malformed, report),
               "strict nonmissing logical")
  malformed <- analysis
  malformed$sample_manifest$is_reference[[1L]] <- NA
  expect_error(facs_report_edu_apex_comparison(malformed, report),
               "strict nonmissing logical")
  malformed <- analysis
  malformed$sample_manifest$reference_condition[[2L]] <- "Different reference"
  expect_error(facs_report_edu_apex_comparison(malformed, report),
               "exactly one correctly identified reference")
})

test_that("embedded plot PDF links are self-contained and filename-safe", {
  link <- facs_report_embedded_plot_pdf_link(
    ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
      ggplot2::geom_point(),
    "../SYNTHETIC unsafe plot.pdf", width = 2, height = 2
  )
  expect_match(link, 'download="SYNTHETIC_unsafe_plot.pdf"', fixed = TRUE)
  expect_match(link, 'href="data:application/pdf;base64,JVBER', fixed = TRUE)
  expect_false(grepl("../", link, fixed = TRUE))
  expect_true(is.numeric(attr(link, "pdf_bytes")))
  escaped <- facs_report_embedded_plot_pdf_link(
    ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
      ggplot2::geom_point(),
    'SYNTHETIC "<unsafe>.pdf', width = 2, height = 2,
    label = "Download <SYNTHETIC> & plot"
  )
  expect_false(grepl('download="[^\"]*<', escaped))
  expect_match(escaped, "Download &lt;SYNTHETIC&gt; &amp; plot", fixed = TRUE)
  expect_error(facs_report_embedded_plot_pdf_link(
    ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
      ggplot2::geom_point(),
    "SYNTHETIC-too-large.pdf", width = 2, height = 2, maximum_bytes = 1
  ), "per-plot safety limit")
})

test_that("all-panel PDF canvas contains every named pseudocolor panel", {
  analysis <- synthetic_edu_report_analysis()
  report <- build_edu_pseudocolor_output_contract(analysis)
  canvas <- facs_report_edu_all_panels_canvas(
    report, max_columns = 2, panel_width = 3, panel_height = 3
  )
  expect_s3_class(canvas$plot, "ggplot")
  expect_identical(canvas$panel_count, length(report$panels))
  expect_identical(canvas$columns, 2L)
  expect_identical(canvas$rows, ceiling(length(report$panels) / 2))
  expect_identical(canvas$width, 6)
  expect_identical(canvas$height, 3 * canvas$rows)
  expect_error(facs_report_edu_all_panels_canvas(list(panels = list())),
               "named ggplot panels")
})

test_that("EdU gating cards preserve every sample and validated population", {
  analysis <- synthetic_edu_report_analysis()
  analysis$config$dna_channel <- "DNA content"
  temporary <- withr::local_tempdir()
  expected_positive_n <- integer(length(analysis$normalized_data))
  for (i in seq_along(analysis$normalized_data)) {
    sample <- analysis$normalized_data[[i]]
    event_n <- nrow(sample$data)
    sample$data$event_identity <- paste0(
      "SYNTHETIC-", i, ":event_index:", seq_len(event_n)
    )
    sample$data[["DNA-H"]] <- seq_len(event_n) * 8
    sample$data[[analysis$config$dna_channel]] <- sample$data$dna_norm
    sample$data[[analysis$config$target_channel]] <- sample$data$target_raw
    sample$edu_positive <- sample$edu_positive[
      seq_len(nrow(sample$edu_positive) - i + 1L), , drop = FALSE
    ]
    positive_n <- nrow(sample$edu_positive)
    expected_positive_n[[i]] <- positive_n
    sample$edu_positive[[analysis$config$dna_channel]] <-
      sample$edu_positive$dna_norm
    sample$edu_positive[[analysis$config$target_channel]] <-
      sample$edu_positive$target_raw
    sample$edu_positive$event_identity <-
      sample$data$event_identity[seq_len(positive_n)]
    g1 <- sample$data[seq_len(8L), c(analysis$config$dna_channel,
                                     analysis$config$target_channel,
                                     "event_identity"),
                                   drop = FALSE]
    sample$g1 <- g1
    sample$g1_file <- file.path(temporary, paste0("SYNTHETIC-g1-", i, ".csv"))
    utils::write.csv(g1, sample$g1_file, row.names = FALSE)
    utils::write.csv(data.frame(invalid = "SYNTHETIC replacement"),
                     sample$g1_file, row.names = FALSE)
    analysis$normalized_data[[i]] <- sample
  }
  analysis$config$gating <- list(
    mode = "model_experimental",
    acquisitions = lapply(seq_len(nrow(analysis$sample_manifest)), function(i) {
      list(prefix = analysis$sample_manifest$prefix[[i]],
           channel_roles = list(dna_area = "DNA-A", dna_height = "DNA-H"))
    })
  )
  all_events <- lapply(analysis$normalized_data, `[[`, "data")
  names(all_events) <- as.character(analysis$sample_manifest$prefix)
  cards <- facs_report_edu_gating_cards(
    analysis, all_events = all_events,
    display_offsets = stats::setNames(
      rep(1000, nrow(analysis$sample_manifest)),
      analysis$sample_manifest$prefix
    ),
    max_points = 10L
  )
  expect_length(cards, nrow(analysis$sample_manifest))
  expect_identical(vapply(cards, `[[`, character(1), "prefix"),
                   as.character(analysis$sample_manifest$prefix))
  expect_true(all(vapply(cards, function(card) inherits(card$plot, "ggplot"),
                         logical(1))))
  expect_identical(cards[[1L]]$single_cells_x_axis, "DNA-A")
  expect_identical(cards[[1L]]$single_cells_y_axis, "DNA-H")
  expect_identical(cards[[1L]]$dna_area_source_channel, "DNA-A")
  expect_identical(cards[[1L]]$dna_height_source_channel, "DNA-H")
  expect_true(all(is.finite(cards[[1L]]$single_cells_x_limits)))
  expect_true(all(is.finite(cards[[1L]]$single_cells_y_limits)))
  expect_true(cards[[1L]]$single_cells_axis_coverage >= 0.99)
  expect_identical(cards[[1L]]$overlay_x_axis, "DNA")
  expect_identical(cards[[1L]]$overlay_y_axis, "EdU")
  expect_identical(cards[[1L]]$background_line_source,
                   "sample_fitted_raw_edu_baseline")
  expect_true(inherits(cards[[1L]]$background_plot, "ggplot"))
  expect_true(is.finite(cards[[1L]]$background_display_offset))
  expect_true(all(vapply(cards, `[[`, integer(1), "g1_n") == 8L))
  expect_identical(vapply(cards, `[[`, integer(1), "edu_positive_n"),
                   expected_positive_n)
  zero_offset_analysis <- analysis
  zero_offset_analysis$normalized_data <- lapply(
    zero_offset_analysis$normalized_data,
    function(sample) {
      sample$data$target_bgsub <- abs(sample$data$target_bgsub) + 1
      sample
    }
  )
  zero_offset_cards <- facs_report_edu_gating_cards(
    zero_offset_analysis, all_events = all_events,
    display_offsets = stats::setNames(
      rep(0, nrow(analysis$sample_manifest)),
      analysis$sample_manifest$prefix
    ),
    max_points = 10L
  )
  expect_length(zero_offset_cards, nrow(analysis$sample_manifest))
  expect_error(facs_report_edu_gating_cards(
                 analysis,
                 display_offsets = stats::setNames(
                   rep(1000, nrow(analysis$sample_manifest)),
                   analysis$sample_manifest$prefix
                 ),
                 max_points = 0),
               "positive point limit")
})

test_that("FlowJo gating cards avoid unconfigured DNA-height assumptions", {
  analysis <- synthetic_edu_report_analysis()
  analysis$config$gating <- list(mode = "flowjo")
  for (i in seq_along(analysis$normalized_data)) {
    sample <- analysis$normalized_data[[i]]
    sample$data[[analysis$config$dna_channel]] <- sample$data$dna_norm
    sample$data[[analysis$config$target_channel]] <- sample$data$target_raw
    sample$edu_positive[[analysis$config$dna_channel]] <-
      sample$edu_positive$dna_norm
    sample$edu_positive[[analysis$config$target_channel]] <-
      sample$edu_positive$target_raw
    sample$g1 <- sample$data[seq_len(8L), c(analysis$config$dna_channel,
                                            analysis$config$target_channel),
                                          drop = FALSE]
    analysis$normalized_data[[i]] <- sample
  }
  cards <- facs_report_edu_gating_cards(
    analysis,
    display_offsets = stats::setNames(
      rep(1000, nrow(analysis$sample_manifest)),
      analysis$sample_manifest$prefix
    ),
    max_points = 10L
  )
  expect_length(cards, nrow(analysis$sample_manifest))
  expect_true(all(vapply(cards, `[[`, character(1), "single_cells_x_axis") ==
                    "DNA"))
  expect_true(all(vapply(cards, `[[`, character(1), "single_cells_y_axis") ==
                    "EdU"))
  expect_true(all(is.na(vapply(cards, `[[`, character(1),
                                   "dna_height_source_channel"))))
})

test_that("FlowJo gating cards use explicit all-events DNA-A and DNA-H", {
  analysis <- synthetic_edu_report_analysis()
  analysis$config$gating <- list(mode = "flowjo")
  all_events <- vector("list", length(analysis$normalized_data))
  names(all_events) <- analysis$sample_manifest$prefix
  for (i in seq_along(analysis$normalized_data)) {
    sample <- analysis$normalized_data[[i]]
    n <- nrow(sample$data)
    sample$data$event_index <- as.character(seq_len(n))
    sample$data[[analysis$config$dna_channel]] <- sample$data$dna_norm
    sample$data[[analysis$config$target_channel]] <- sample$data$target_raw
    sample$data[["DNA-H"]] <- sample$data$dna_norm * 0.9
    sample$edu_positive[[analysis$config$dna_channel]] <-
      sample$edu_positive$dna_norm
    sample$edu_positive[[analysis$config$target_channel]] <-
      sample$edu_positive$target_raw
    sample$edu_positive$event_index <- sample$data$event_index[
      seq_len(nrow(sample$edu_positive))
    ]
    sample$g1 <- sample$data[seq_len(8L), c(
      analysis$config$dna_channel, analysis$config$target_channel,
      "event_index"
    ), drop = FALSE]
    analysis$normalized_data[[i]] <- sample
    extra <- sample$data
    extra$event_index <- as.character(seq_len(n) + n)
    all_events[[i]] <- rbind(sample$data, extra)
  }
  cards <- facs_report_edu_gating_cards(
    analysis, all_events = all_events, dna_height_channel = "DNA-H",
    display_offsets = stats::setNames(
      rep(1000, nrow(analysis$sample_manifest)),
      analysis$sample_manifest$prefix
    ), max_points = 10L
  )
  expect_true(all(vapply(cards, `[[`, character(1), "single_cells_x_axis") ==
                    "DNA-A"))
  expect_true(all(vapply(cards, `[[`, character(1), "single_cells_y_axis") ==
                    "DNA-H"))
  expect_true(all(vapply(cards, `[[`, character(1),
                               "dna_height_source_channel") == "DNA-H"))
  mismatched <- all_events
  mismatched[[1L]]$event_index <- paste0("wrong-", mismatched[[1L]]$event_index)
  expect_error(facs_report_edu_gating_cards(
    analysis, all_events = mismatched, dna_height_channel = "DNA-H",
    display_offsets = stats::setNames(
      rep(1000, nrow(analysis$sample_manifest)),
      analysis$sample_manifest$prefix
    ), max_points = 10L
  ), "not contained")
})

test_that("gating card errors name the missing population and column", {
  analysis <- synthetic_edu_report_analysis()
  analysis$config$gating <- list(mode = "flowjo")
  all_events <- vector("list", length(analysis$normalized_data))
  names(all_events) <- analysis$sample_manifest$prefix
  for (i in seq_along(analysis$normalized_data)) {
    sample <- analysis$normalized_data[[i]]
    n <- nrow(sample$data)
    sample$data$event_index <- as.character(seq_len(n))
    sample$data[[analysis$config$dna_channel]] <- sample$data$dna_norm
    sample$data[[analysis$config$target_channel]] <- sample$data$target_raw
    sample$data[["DNA-H"]] <- sample$data$dna_norm * 0.9
    sample$edu_positive[[analysis$config$dna_channel]] <-
      sample$edu_positive$dna_norm
    sample$edu_positive[[analysis$config$target_channel]] <-
      sample$edu_positive$target_raw
    sample$edu_positive$event_index <- sample$data$event_index[
      seq_len(nrow(sample$edu_positive))
    ]
    sample$g1 <- sample$data[seq_len(8L), c(
      analysis$config$dna_channel, analysis$config$target_channel,
      "event_index"
    ), drop = FALSE]
    analysis$normalized_data[[i]] <- sample
    all_events[[i]] <- sample$data
  }
  display_offsets <- stats::setNames(
    rep(1000, nrow(analysis$sample_manifest)), analysis$sample_manifest$prefix
  )

  missing_single_height <- analysis
  missing_single_height$normalized_data[[1L]]$data[["DNA-H"]] <- NULL
  expect_error(facs_report_edu_gating_cards(
    missing_single_height, all_events = all_events, dna_height_channel = "DNA-H",
    display_offsets = display_offsets, max_points = 10L
  ), "Single Cells missing column\\(s\\): DNA-H")

  missing_all_events_height <- all_events
  missing_all_events_height[[1L]][["DNA-H"]] <- NULL
  expect_error(facs_report_edu_gating_cards(
    analysis, all_events = missing_all_events_height, dna_height_channel = "DNA-H",
    display_offsets = display_offsets, max_points = 10L
  ), "missing column\\(s\\): DNA-H")

  empty_g1 <- analysis
  empty_g1$normalized_data[[1L]]$g1 <- empty_g1$normalized_data[[1L]]$g1[0L, , drop = FALSE]
  expect_error(facs_report_edu_gating_cards(
    empty_g1, all_events = all_events, dna_height_channel = "DNA-H",
    display_offsets = display_offsets, max_points = 10L
  ), "G1 has zero events")
})

test_that("S-phase cards use established regional assignments for every sample", {
  analysis <- synthetic_edu_report_analysis()
  analysis$normalized_data[[1L]]$data$target_bgsub[[13L]] <- Inf
  report <- build_edu_pseudocolor_output_contract(analysis)
  cards <- facs_report_edu_s_phase_cards(analysis, report, max_points = 10L)
  expect_length(cards, nrow(analysis$sample_manifest))
  expect_identical(vapply(cards, `[[`, character(1), "prefix"),
                   as.character(analysis$sample_manifest$prefix))
  expect_true(all(vapply(cards, function(card) inherits(card$plot, "ggplot"),
                         logical(1))))
  expect_true(all(vapply(cards, `[[`, integer(1), "early_s_n") >= 0L))
  expect_true(all(vapply(cards, `[[`, integer(1), "mid_s_n") >= 0L))
  expect_true(all(vapply(cards, `[[`, integer(1), "late_s_n") >= 0L))
  gates <- analysis$quantitation$gates
  early <- gates[gates$gate == "Early S", , drop = FALSE]
  expected_early <- with(analysis$normalized_data[[1L]]$data,
                         sum(edu_computed_positive %in% TRUE &
                               dna_norm >= early$xmin[[1L]] &
                               dna_norm < early$xmax[[1L]]))
  expect_identical(cards[[1L]]$early_s_n, as.integer(expected_early))
  expect_error(facs_report_edu_s_phase_cards(analysis, report, max_points = 0),
               "established regional gates")
})
