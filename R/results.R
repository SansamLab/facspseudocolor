# Result-based quantitation, plotting, and saving ----------------------------

validate_analysis_object <- function(analysis) {
  if (!inherits(analysis, "facs_analysis")) {
    stop("`analysis` must be a facs_analysis object.", call. = FALSE)
  }
  invisible(analysis)
}

analysis_display_rows <- function(analysis) {
  show_reference <- analysis$config$show_reference_panel
  if (is.null(show_reference)) show_reference <- analysis$config$plot_type == "edu"
  if (isTRUE(show_reference)) {
    rep(TRUE, nrow(analysis$sample_manifest))
  } else {
    !analysis$sample_manifest$is_reference
  }
}

analysis_y_limits <- function(analysis, rows = analysis_display_rows(analysis)) {
  if (!is.null(analysis$config$y_limits)) {
    return(as.numeric(unlist(analysis$config$y_limits)))
  }
  values <- unlist(lapply(analysis$normalized_data[rows], function(sample) {
    value <- sample$data$target_norm
    value[is.finite(value) & value > 0]
  }), use.names = FALSE)
  if (length(values) < 2L) {
    stop("Too few positive normalized target values to determine y limits.",
         call. = FALSE)
  }
  stats::quantile(
    values,
    c(analysis$config$y_limit_lower_quantile,
      analysis$config$y_limit_upper_quantile),
    names = FALSE
  )
}

analysis_gate_rectangles <- function(analysis, y_limits) {
  config <- analysis$config
  numeric_or_null <- function(value) {
    if (is.null(value)) NULL else as.numeric(unlist(value))
  }
  if (config$plot_type == "edu") {
    make_rectangular_phase_gates(
      s_phase_bins = config$s_phase_bins,
      g1_x_range = numeric_or_null(config$g1_x_range),
      g2m_x_range = numeric_or_null(config$g2m_x_range),
      negative_y_range = numeric_or_null(config$negative_y_range),
      s_phase_y_range = numeric_or_null(config$s_phase_y_range),
      dna_2n_value = config$dna_2n_value,
      y_limits = y_limits
    )
  } else {
    make_dna_phase_gates(
      s_phase_bins = config$s_phase_bins,
      g1_x_range = numeric_or_null(config$g1_x_range),
      g2m_x_range = numeric_or_null(config$g2m_x_range),
      dna_2n_value = config$dna_2n_value,
      y_limits = y_limits
    )
  }
}

#' Quantify cell-cycle phases from an analysis result
#'
#' Calculates result tables without drawing or saving plots. The returned
#' `facs_analysis` contains the original normalized data plus phase medians,
#' whole-population medians, phase percentages, gate definitions, and optional
#' within-replicate reference ratios.
#'
#' @param analysis A `facs_analysis` object.
#' @param include Any of `"phase_median"`, `"whole_median"`, and
#'   `"phase_percent"`.
#' @param signal Either `"background_subtracted"` or `"normalized"`.
#' @param reference_condition Optional displayed condition used for
#'   within-replicate ratios. If omitted, ratios are added only when requested
#'   by the configuration.
#'
#' @return An updated `facs_analysis` object.
#' @export
quantify_cell_cycle <- function(
    analysis,
    include = c("phase_median", "whole_median", "phase_percent"),
    signal = analysis$config$quant_signal,
    reference_condition = NULL
) {
  validate_analysis_object(analysis)
  allowed <- c("phase_median", "whole_median", "phase_percent")
  unknown <- setdiff(include, allowed)
  if (length(unknown)) {
    stop("Unknown quantitation type(s): ", paste(unknown, collapse = ", "),
         call. = FALSE)
  }
  signal <- match.arg(signal, c("background_subtracted", "normalized"))
  value_column <- if (signal == "normalized") "target_norm" else "target_bgsub"

  rows <- analysis_display_rows(analysis)
  manifest <- analysis$sample_manifest[rows, , drop = FALSE]
  samples <- analysis$normalized_data[rows]
  result_wrappers <- lapply(samples, function(sample) list(sample_data = sample))
  y_limits <- analysis_y_limits(analysis, rows)
  gates <- analysis_gate_rectangles(analysis, y_limits)
  quantitation <- list(
    signal = signal,
    y_limits = y_limits,
    gates = gates
  )

  source_field <- if (analysis$config$plot_type == "edu") {
    "edu_positive"
  } else {
    "data"
  }
  if ("phase_median" %in% include) {
    windows <- make_phase_windows(
      s_phase_bins = analysis$config$s_phase_bins,
      g1_x_range = analysis$config$g1_x_range,
      g2m_x_range = analysis$config$g2m_x_range,
      dna_2n_value = analysis$config$dna_2n_value,
      include_g1_g2m = analysis$config$plot_type == "poi"
    )
    quantitation$phase_medians <- collect_phase_signal_medians(
      result_wrappers, manifest, windows,
      source_field = source_field, y_col = value_column
    )
  }
  if ("whole_median" %in% include) {
    quantitation$whole_medians <- collect_whole_population_medians(
      result_wrappers, manifest, source_field = source_field,
      y_col = value_column
    )
  }
  if ("phase_percent" %in% include) {
    quantitation$phase_percentages <- collect_phase_counts(
      result_wrappers, manifest, gates
    )
  }

  ratio_requested <- !is.null(reference_condition) ||
    isTRUE(analysis$config$quantify_reference_normalized)
  if (ratio_requested) {
    if (is.null(reference_condition)) {
      reference_condition <- analysis$config$quant_reference_condition
    }
    if (is.null(reference_condition) && analysis$config$plot_type == "edu") {
      refs <- unique(manifest$condition[manifest$is_reference])
      if (length(refs) == 1L) reference_condition <- refs[[1]]
    }
    if (!config_scalar_string(reference_condition)) {
      stop(
        "Reference-normalized quantitation requires an explicit displayed `reference_condition`.",
        call. = FALSE
      )
    }
    quantitation$reference_condition <- reference_condition
    if (!is.null(quantitation$phase_medians)) {
      quantitation$phase_medians_reference <- add_reference_ratio(
        quantitation$phase_medians, "median_signal", reference_condition,
        phase_col = "phase"
      )
    }
    if (!is.null(quantitation$whole_medians)) {
      quantitation$whole_medians_reference <- add_reference_ratio(
        quantitation$whole_medians, "median_signal", reference_condition
      )
    }
  }

  analysis$quantitation <- quantitation
  analysis
}

#' Plot configured cell-cycle quantitation results
#'
#' Converts tables produced by [quantify_cell_cycle()] into ggplot objects.
#' This function writes no files.
#'
#' @param analysis A quantified `facs_analysis` object.
#' @param seed Explicit seed used only if gate-assignment points must be
#'   downsampled.
#'
#' @return A named list of ggplot objects.
#' @export
plot_facs_quantitation <- function(analysis, seed = 1L) {
  validate_analysis_object(analysis)
  quantitation <- analysis$quantitation
  if (!length(quantitation) || is.null(quantitation$gates)) {
    stop("Run `quantify_cell_cycle()` before plotting quantitation.",
         call. = FALSE)
  }
  config <- analysis$config
  error_bar <- config$quant_error_bar
  show_points <- config$quant_show_points
  colors <- config$bar_colors
  signal_label <- if (quantitation$signal == "normalized") {
    "normalized"
  } else {
    "background-subtracted"
  }
  plots <- list()

  if (!is.null(quantitation$phase_medians)) {
    y_title <- paste("Median", signal_label, config$target_name)
    plots$phase_median <- make_phase_signal_barplot(
      quantitation$phase_medians, error_bar = error_bar, y_title = y_title,
      show_points = show_points, fill_colors = colors
    )
    if (isTRUE(config$quant_phase_lineplot)) {
      plots$phase_median_line <- make_phase_signal_lineplot(
        quantitation$phase_medians, error_bar = error_bar, y_title = y_title,
        show_points = show_points, fill_colors = colors
      )
    }
  }
  if (!is.null(quantitation$whole_medians)) {
    plots$whole_median <- make_whole_population_barplot(
      quantitation$whole_medians, error_bar = error_bar,
      y_title = paste("Median", signal_label, config$target_name),
      show_points = show_points, fill_colors = colors
    )
  }
  if (!is.null(quantitation$phase_percentages)) {
    plots$phase_percent <- make_phase_percentage_plot(
      quantitation$phase_percentages, error_bar = error_bar,
      show_points = show_points, fill_colors = colors
    )
  }
  if (!is.null(quantitation$phase_medians_reference)) {
    plots$phase_median_reference <- make_phase_signal_barplot(
      quantitation$phase_medians_reference, error_bar = error_bar,
      y_title = "Signal relative to reference", show_points = show_points,
      fill_colors = colors, value_col = "ratio", ref_line = 1
    )
    if (isTRUE(config$quant_phase_lineplot)) {
      plots$phase_median_reference_line <- make_phase_signal_lineplot(
        quantitation$phase_medians_reference, error_bar = error_bar,
        y_title = "Signal relative to reference", show_points = show_points,
        fill_colors = colors, value_col = "ratio", ref_line = 1
      )
    }
  }
  if (!is.null(quantitation$whole_medians_reference)) {
    plots$whole_median_reference <- make_whole_population_barplot(
      quantitation$whole_medians_reference, error_bar = error_bar,
      y_title = "Signal relative to reference", show_points = show_points,
      fill_colors = colors, value_col = "ratio", ref_line = 1
    )
  }

  if (isTRUE(config$show_gate_assignment)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
      stop("`seed` must be one explicit finite number.", call. = FALSE)
    }
    rows <- analysis_display_rows(analysis)
    manifest <- analysis$sample_manifest[rows, , drop = FALSE]
    wrappers <- lapply(
      analysis$normalized_data[rows], function(sample) list(sample_data = sample)
    )
    old_seed_exists <- exists(".Random.seed", envir = .GlobalEnv,
                              inherits = FALSE)
    if (old_seed_exists) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
      if (old_seed_exists) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
    assignments <- collect_gate_assignments(
      wrappers, manifest, quantitation$gates,
      max_points = config$gate_assignment_max_points
    )
    plots$gate_assignment <- make_gate_assignment_plot(
      assignments, quantitation$gates,
      x_limits = as.numeric(unlist(config$x_limits)),
      y_limits = quantitation$y_limits,
      y_log10 = config$y_log10,
      dna_2n_value = config$dna_2n_value
    )
  }
  plots
}

#' Build pseudocolor panels from a normalized analysis result
#'
#' Uses normalized event data already stored in `analysis`; input CSV files are
#' not read again. The function returns plot objects and writes no files.
#'
#' @param analysis A `facs_analysis` object.
#' @param palette A named package palette or a custom color vector. Defaults to
#'   the validated configuration value.
#'
#' @return An object of class `facs_plot_set` containing individual panels,
#'   layout information, and an assembled cowplot grid when applicable.
#' @export
plot_pseudocolor_panels <- function(analysis, palette = analysis$config$palette) {
  validate_analysis_object(analysis)
  rows <- analysis_display_rows(analysis)
  manifest <- analysis$sample_manifest[rows, , drop = FALSE]
  samples <- analysis$normalized_data[rows]
  settings <- analysis_settings(analysis$config)
  settings$y_limits <- analysis_y_limits(analysis, rows)
  settings$y_axis_title <- if (analysis$config$plot_type == "edu") {
    paste0(analysis$config$target_name, "\n(baseline-normalized)")
  } else {
    paste0(analysis$config$target_name, "\n(background-normalized)")
  }
  palette <- resolve_palette(palette)

  panel_results <- lapply(seq_along(samples), function(i) {
    make_pseudocolor_plot_from_data(
      sample_data = samples[[i]],
      condition_label = manifest$condition[[i]],
      settings = settings,
      palette = palette
    )
  })
  gates <- analysis_gate_rectangles(analysis, settings$y_limits)
  if (isTRUE(analysis$config$show_phase_gates)) {
    panel_results <- lapply(panel_results, function(result) {
      result$panel <- add_phase_gates_to_plot(
        result$panel, gates,
        color = analysis$config$gate_color,
        linetype = analysis$config$gate_linetype,
        linewidth = analysis$config$gate_linewidth,
        show_labels = analysis$config$gate_show_labels,
        label_size = analysis$config$gate_label_size,
        y_log10 = analysis$config$y_log10
      )
      result$plot <- add_phase_gates_to_plot(
        result$plot, gates,
        color = analysis$config$gate_color,
        linetype = analysis$config$gate_linetype,
        linewidth = analysis$config$gate_linewidth,
        show_labels = analysis$config$gate_show_labels,
        label_size = analysis$config$gate_label_size,
        y_log10 = analysis$config$y_log10
      )
      result
    })
  }

  n_columns <- length(unique(manifest$condition_index))
  n_rows <- length(unique(manifest$replicate_index))
  replicate_labels <- unique(manifest$replicate)
  condition_labels <- unique(manifest$condition)
  page_layout <- facs_page_layout(
    n_columns, n_rows, analysis$config$layout_options
  )
  grid <- NULL
  if (analysis$config$layout == "cowplot") {
    formatted <- format_panel_plots(
      lapply(panel_results, `[[`, "plot"), n_columns
    )
    grid <- make_panel_grid(
      formatted, n_columns, row_labels = replicate_labels
    )
  }

  structure(list(
    layout = analysis$config$layout,
    panels = lapply(panel_results, `[[`, "panel"),
    decorated_plots = lapply(panel_results, `[[`, "plot"),
    panel_results = panel_results,
    grid = grid,
    page_layout = page_layout,
    condition_labels = condition_labels,
    replicate_labels = replicate_labels,
    n_columns = n_columns,
    n_rows = n_rows,
    x_limits = settings$x_limits,
    y_limits = settings$y_limits,
    y_log10 = settings$y_log10,
    dna_2n_value = settings$dna_2n_value,
    y_axis_title = settings$y_axis_title
  ), class = "facs_plot_set")
}

draw_facs_plot_set <- function(x) {
  if (x$layout == "cowplot") {
    print(x$grid)
  } else {
    assemble_facs_plotgardener(
      panels = x$panels,
      n_cols = x$n_columns,
      n_rows = x$n_rows,
      condition_labels = x$condition_labels,
      replicate_labels = x$replicate_labels,
      x_limits = x$x_limits,
      y_limits = x$y_limits,
      y_log10 = x$y_log10,
      x_breaks = c(x$dna_2n_value, 2 * x$dna_2n_value),
      x_labels = c("2N", "4N"),
      y_axis_title = x$y_axis_title,
      x_axis_title = "DNA content",
      layout = x$page_layout
    )
  }
  invisible(x)
}

#' @export
print.facs_plot_set <- function(x, ...) draw_facs_plot_set(x)

resolve_analysis_output <- function(analysis, path) {
  if (is.null(path)) return(NULL)
  if (!config_scalar_string(path)) {
    stop("Every output path must be a nonempty string.", call. = FALSE)
  }
  path <- path.expand(path)
  if (!grepl("^(/|[A-Za-z]:[/\\\\])", path)) {
    config_dir <- attr(analysis$config, "config_dir")
    if (!config_scalar_string(config_dir)) {
      stop("Relative output paths require a file-backed configuration.",
           call. = FALSE)
    }
    path <- file.path(config_dir, path)
  }
  normalizePath(path, mustWork = FALSE)
}

#' Save explicitly requested analysis outputs
#'
#' This is the package's file-writing layer. Existing files are never
#' overwritten unless `overwrite = TRUE`.
#'
#' @param analysis A `facs_analysis` object.
#' @param plots Optional `facs_plot_set`; generated from `analysis` when a PDF
#'   or PNG is requested and `plots` is `NULL`.
#' @param output_pdf Explicit PDF path, or `NULL` for no PDF. Defaults to the
#'   configured path.
#' @param output_png Explicit PNG path, or `NULL` for no PNG. Defaults to the
#'   configured path.
#' @param output_rds Explicit path for the complete analysis object, or `NULL`
#'   to skip serialization.
#' @param overwrite Whether existing requested outputs may be replaced.
#' @param dpi PNG resolution.
#'
#' @return Invisibly, a named character vector of files written.
#' @export
save_facs_results <- function(
    analysis, plots = NULL,
    output_pdf = analysis$config$output_pdf,
    output_png = analysis$config$output_png,
    output_rds = NULL,
    overwrite = FALSE,
    dpi = 300
) {
  validate_analysis_object(analysis)
  paths <- c(
    pdf = resolve_analysis_output(analysis, output_pdf),
    png = resolve_analysis_output(analysis, output_png),
    rds = resolve_analysis_output(analysis, output_rds)
  )
  paths <- paths[nzchar(paths)]
  if (!length(paths)) {
    stop("At least one explicit output path must be requested.", call. = FALSE)
  }
  duplicates <- unique(paths[duplicated(paths)])
  if (length(duplicates)) {
    stop("Output paths must be distinct: ", paste(duplicates, collapse = ", "),
         call. = FALSE)
  }
  existing <- paths[file.exists(paths)]
  if (length(existing) && !isTRUE(overwrite)) {
    stop(
      "Refusing to overwrite existing output file(s):\n- ",
      paste(existing, collapse = "\n- "),
      "\nSet `overwrite = TRUE` to replace them.", call. = FALSE
    )
  }
  for (directory in unique(dirname(paths))) {
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  }

  figure_requested <- any(names(paths) %in% c("pdf", "png"))
  if (figure_requested) {
    if (is.null(plots)) plots <- plot_pseudocolor_panels(analysis)
    if (!inherits(plots, "facs_plot_set")) {
      stop("`plots` must be a facs_plot_set object.", call. = FALSE)
    }
  }
  if ("pdf" %in% names(paths)) {
    width <- if (plots$layout == "plotgardener") plots$page_layout$width else
      analysis$config$pdf_width
    height <- if (plots$layout == "plotgardener") plots$page_layout$height else
      analysis$config$pdf_height_per_row * plots$n_rows
    grDevices::pdf(paths[["pdf"]], width = width, height = height,
                   bg = "white", useDingbats = FALSE)
    tryCatch(draw_facs_plot_set(plots), finally = grDevices::dev.off())
  }
  if ("png" %in% names(paths)) {
    width <- if (plots$layout == "plotgardener") plots$page_layout$width else
      analysis$config$pdf_width
    height <- if (plots$layout == "plotgardener") plots$page_layout$height else
      analysis$config$pdf_height_per_row * plots$n_rows
    grDevices::png(paths[["png"]], width = width, height = height,
                   units = "in", res = dpi, bg = "white")
    tryCatch(draw_facs_plot_set(plots), finally = grDevices::dev.off())
  }
  if ("rds" %in% names(paths)) saveRDS(analysis, paths[["rds"]])
  invisible(paths)
}
