# Report-specific appearance settings ---------------------------------------

facs_appearance_keys <- function() {
  c(
    "target_axis_label", "dna_axis_label", "palette", "condition_palette",
    "condition_colors", "x_limits", "y_limits", "y_log10", "point_size",
    "density_bandwidth", "density_lower_clip", "density_upper_clip",
    "density_gamma", "show_phase_gates", "gate_color", "gate_linetype",
    "gate_linewidth", "gate_show_labels", "gate_label_size", "base_font_size",
    "error_bar", "show_points", "phase_lineplot", "layout_options"
  )
}

#' Read optional report appearance overrides
#'
#' @param path Explicit YAML file containing presentation settings only.
#' @return A named list of appearance overrides.
#' @export
read_facs_appearance <- function(path) {
  if (!config_scalar_string(path) || !file.exists(path)) {
    stop("Appearance file not found: ", path, call. = FALSE)
  }
  value <- yaml::read_yaml(normalizePath(path, mustWork = TRUE))
  if (is.null(value)) value <- list()
  if (!is.list(value) || is.data.frame(value) || is.null(names(value))) {
    stop("Appearance YAML must contain a top-level mapping.", call. = FALSE)
  }
  value
}

facs_named_palette <- function(name, n) {
  name <- name %||% "default"
  if (!config_scalar_string(name) || !name %in% c("default", "colorblind", "viridis")) {
    stop("`condition_palette` must be 'default', 'colorblind', or 'viridis'.",
         call. = FALSE)
  }
  if (name == "default") return(condition_fill_palette(n))
  anchors <- if (name == "colorblind") {
    c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
      "#0072B2", "#D55E00", "#CC79A7")
  } else {
    c("#440154", "#414487", "#2A788E", "#22A884", "#7AD151", "#FDE725")
  }
  if (name == "viridis") return(grDevices::colorRampPalette(anchors)(n))
  if (n <= length(anchors)) anchors[seq_len(n)] else
    grDevices::colorRampPalette(anchors)(n)
}

#' Resolve package defaults and optional presentation overrides
#' @param analysis A complete `facs_analysis` object.
#' @param overrides Optional named presentation settings. Null entries are ignored.
#' @param appearance_file Optional appearance YAML file.
#' @return A validated `facs_appearance` object.
#' @export
resolve_facs_appearance <- function(analysis, overrides = NULL, appearance_file = NULL) {
  validate_analysis_object(analysis)
  if (!is.null(appearance_file)) {
    file_values <- read_facs_appearance(appearance_file)
    overrides <- utils::modifyList(file_values, overrides %||% list(), keep.null = TRUE)
  }
  overrides <- overrides %||% list()
  if (!is.list(overrides) || is.data.frame(overrides)) {
    stop("`appearance` must be a named list.", call. = FALSE)
  }
  unknown <- setdiff(names(overrides), facs_appearance_keys())
  if (length(unknown)) {
    stop("Unknown appearance setting(s): ", paste(unknown, collapse = ", "),
         call. = FALSE)
  }
  overrides <- overrides[!vapply(overrides, is.null, logical(1))]
  config <- analysis$config
  defaults <- list(
    target_axis_label = if (config$plot_type == "edu") {
      paste0(config$target_name, "\n(baseline-normalized)")
    } else paste0(config$target_name, "\n(background-normalized)"),
    dna_axis_label = "DNA content", palette = config$palette,
    condition_palette = "viridis", condition_colors = config$bar_colors,
    x_limits = as.numeric(unlist(config$x_limits)),
    y_limits = analysis_y_limits(analysis), y_log10 = config$y_log10,
    point_size = config$point_size, density_bandwidth = config$density_bandwidth,
    density_lower_clip = config$density_lower_clip,
    density_upper_clip = config$density_upper_clip,
    density_gamma = config$density_gamma,
    show_phase_gates = config$show_phase_gates, gate_color = config$gate_color,
    gate_linetype = config$gate_linetype, gate_linewidth = config$gate_linewidth,
    gate_show_labels = config$gate_show_labels,
    gate_label_size = config$gate_label_size, base_font_size = 11,
    error_bar = config$quant_error_bar, show_points = config$quant_show_points,
    phase_lineplot = config$quant_phase_lineplot,
    layout_options = config$layout_options
  )
  out <- utils::modifyList(defaults, overrides, keep.null = TRUE)
  for (key in c("x_limits", "y_limits")) {
    if (!config_numeric_range(out[[key]])) {
      stop("`", key, "` must contain two increasing finite numbers.", call. = FALSE)
    }
    out[[key]] <- as.numeric(unlist(out[[key]]))
  }
  if (isTRUE(out$y_log10) && out$y_limits[[1]] <= 0) {
    stop("Logarithmic y-axis limits must be positive.", call. = FALSE)
  }
  if (!out$error_bar %in% c("sd", "sem", "none")) {
    stop("`error_bar` must be 'sd', 'sem', or 'none'.", call. = FALSE)
  }
  conditions <- unique(analysis$sample_manifest$condition[analysis_display_rows(analysis)])
  colors <- stats::setNames(facs_named_palette(out$condition_palette, length(conditions)),
                            conditions)
  if (!is.null(out$condition_colors)) {
    supplied <- unlist(out$condition_colors)
    if (is.null(names(supplied)) || any(!nzchar(names(supplied)))) {
      stop("`condition_colors` must map condition labels to colors.", call. = FALSE)
    }
    unknown_conditions <- setdiff(names(supplied), conditions)
    if (length(unknown_conditions)) {
      stop("Unknown condition color label(s): ",
           paste(unknown_conditions, collapse = ", "), call. = FALSE)
    }
    tryCatch(grDevices::col2rgb(unname(supplied)), error = function(e) {
      stop("Invalid value in `condition_colors`: ", conditionMessage(e), call. = FALSE)
    })
    colors[names(supplied)] <- unname(supplied)
  }
  out$condition_colors <- colors
  class(out) <- c("facs_appearance", "list")
  out
}
