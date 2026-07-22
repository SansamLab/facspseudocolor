facs_report_input <- function(config = NULL, analysis_rds = NULL) {
  has_config <- is.character(config) && length(config) == 1L && nzchar(config)
  has_rds <- is.character(analysis_rds) && length(analysis_rds) == 1L && nzchar(analysis_rds)
  if (!xor(has_config, has_rds)) {
    stop("Supply exactly one of `config` or `analysis_rds`.", call. = FALSE)
  }
  if (has_config) analyze_facs_experiment(config) else read_facs_analysis(analysis_rds)
}

facs_report_appearance <- function(params) {
  keys <- c(
    "target_axis_label", "dna_axis_label", "palette", "condition_palette",
    "condition_colors", "x_limits", "y_limits", "y_log10", "point_size",
    "density_bandwidth", "density_lower_clip", "density_upper_clip",
    "density_gamma", "show_phase_gates", "gate_color", "gate_linetype",
    "gate_linewidth", "gate_show_labels", "gate_label_size", "base_font_size",
    "error_bar", "show_points", "phase_lineplot", "layout_options"
  )
  values <- params[intersect(keys, names(params))]
  values[!vapply(values, is.null, logical(1))]
}

facs_report_save_artifacts <- function(
    analysis, bundle, analysis_path = NULL, bundle_path = NULL,
    overwrite = FALSE
) {
  if (is.character(analysis_path) && length(analysis_path) == 1L && nzchar(analysis_path)) {
    save_facs_analysis(analysis, analysis_path, overwrite = overwrite)
  }
  if (is.character(bundle_path) && length(bundle_path) == 1L && nzchar(bundle_path)) {
    save_facs_figure_bundle(bundle, bundle_path, overwrite = overwrite)
  }
  invisible(NULL)
}
