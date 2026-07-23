# FlowJo pH3 gate geometry and focused report plots -------------------------

#' Read and validate exact pH3 gate geometry exported from FlowJo
#'
#' The sidecar must be produced by the repository-level
#' `python/export_flowjo_gate_geometry.py` script and mapped to configured
#' sample prefixes by `prepare_flowjo_gate_geometry_external()`. Gate boundaries
#' are never inferred from event extrema.
#'
#' @param path Explicit gate-geometry CSV path.
#' @param analysis A completed PH3-mode analysis.
#'
#' @return A validated data frame with raw and normalized plotting coordinates.
#' @export
read_ph3_gate_geometry <- function(path, analysis) {
  validate_analysis_object(analysis)
  if (!identical(analysis$config$plot_type, "ph3")) {
    stop("pH3 gate geometry requires a PH3-mode analysis.", call. = FALSE)
  }
  if (!config_scalar_string(path) || !file.exists(path)) {
    stop("FlowJo gate-geometry file not found: ", path, call. = FALSE)
  }
  path <- normalizePath(path, mustWork = TRUE)
  geometry <- utils::read.csv(path, check.names = FALSE)
  required <- c(
    "prefix", "sample_id", "gate_name", "gate_path", "gate_type",
    "vertex_index", "x_channel", "y_channel", "x_raw", "y_raw", "workspace"
  )
  missing <- setdiff(required, names(geometry))
  if (length(missing)) {
    stop("FlowJo gate geometry is missing column(s): ",
         paste(missing, collapse = ", "), ".", call. = FALSE)
  }
  for (column in c("vertex_index", "x_raw", "y_raw")) {
    if (!is.numeric(geometry[[column]]) || any(!is.finite(geometry[[column]]))) {
      stop("FlowJo gate geometry column `", column,
           "` must contain only finite numeric values.", call. = FALSE)
    }
  }
  if (!nrow(geometry) || anyNA(geometry$prefix) ||
      any(!nzchar(geometry$prefix))) {
    stop("FlowJo gate geometry contains no usable mapped vertices.",
         call. = FALSE)
  }
  expected <- analysis$sample_manifest$prefix
  observed <- unique(geometry$prefix)
  missing_samples <- setdiff(expected, observed)
  unknown_samples <- setdiff(observed, expected)
  if (length(missing_samples) || length(unknown_samples)) {
    details <- c(
      if (length(missing_samples)) paste0(
        "missing configured prefixes: ", paste(missing_samples, collapse = ", ")
      ),
      if (length(unknown_samples)) paste0(
        "unknown prefixes: ", paste(unknown_samples, collapse = ", ")
      )
    )
    stop("FlowJo gate geometry sample mismatch (",
         paste(details, collapse = "; "), ").", call. = FALSE)
  }

  groups <- split(geometry, geometry$prefix)
  for (prefix in names(groups)) {
    group <- groups[[prefix]]
    identity_columns <- c(
      "sample_id", "gate_name", "gate_path", "gate_type",
      "x_channel", "y_channel", "workspace"
    )
    ambiguous <- identity_columns[
      vapply(group[identity_columns], function(x) length(unique(x)) != 1L,
             logical(1))
    ]
    if (length(ambiguous)) {
      stop("Ambiguous FlowJo gate metadata for ", prefix, ": ",
           paste(ambiguous, collapse = ", "), ".", call. = FALSE)
    }
    if (!group$gate_type[[1]] %in% c("PolygonGate", "RectangleGate")) {
      stop("Unsupported FlowJo gate type for ", prefix, ": ",
           group$gate_type[[1]], ".", call. = FALSE)
    }
    ordered <- group[order(group$vertex_index), , drop = FALSE]
    if (nrow(ordered) < 4L ||
        !isTRUE(all.equal(
          unname(as.numeric(ordered[1, c("x_raw", "y_raw")])),
          unname(as.numeric(ordered[nrow(ordered), c("x_raw", "y_raw")])),
          tolerance = 1e-10
        ))) {
      stop("FlowJo gate polygon is not closed for ", prefix, ".",
           call. = FALSE)
    }
  }

  factors <- vapply(expected, function(prefix) {
    model <- analysis$models[[prefix]]
    value <- model$dna_normalization_factor
    if (!config_scalar_number(value) || value <= 0) {
      stop("Missing DNA normalization factor for ", prefix, ".",
           call. = FALSE)
    }
    value
  }, numeric(1))
  geometry$dna_norm <- geometry$x_raw * factors[geometry$prefix]
  geometry$target_norm <- geometry$y_raw
  geometry <- geometry[
    order(match(geometry$prefix, expected), geometry$vertex_index),
    , drop = FALSE
  ]
  rownames(geometry) <- NULL
  attr(geometry, "geometry_path") <- path
  class(geometry) <- c("ph3_gate_geometry", "data.frame")
  geometry
}

ph3_gate_label_position <- function(gate, y_log10) {
  x <- mean(range(gate$dna_norm))
  y_range <- range(gate$target_norm)
  y <- if (isTRUE(y_log10) && all(y_range > 0)) {
    10^(mean(log10(y_range)))
  } else {
    mean(y_range)
  }
  c(x = x, y = y)
}

#' Plot focused pH3 gate pseudocolor panels
#'
#' Overlays the exact FlowJo polygon and annotates the percentage of all Single
#' Cell events contained in the exported pH3-positive population.
#'
#' @param analysis A completed PH3-mode analysis.
#' @param geometry A value returned by [read_ph3_gate_geometry()].
#' @param appearance Optional presentation overrides.
#' @param appearance_file Optional appearance YAML.
#' @param gate_fill Polygon fill color.
#' @param gate_alpha Polygon fill opacity.
#' @param label_color Percentage-label color.
#' @param label_size Percentage-label size.
#' @param label_digits Digits after the decimal point.
#'
#' @return An editable `facs_plot_set`.
#' @export
plot_ph3_4n_gate_panels <- function(
    analysis, geometry, appearance = NULL, appearance_file = NULL,
    gate_fill = "#FFFFFF", gate_alpha = 0.08,
    label_color = "#FFFFFF", label_size = 3.5, label_digits = 1L
) {
  validate_analysis_object(analysis)
  if (!inherits(geometry, "ph3_gate_geometry")) {
    stop("`geometry` must come from `read_ph3_gate_geometry()`.",
         call. = FALSE)
  }
  if (!config_scalar_number(gate_alpha) || gate_alpha < 0 || gate_alpha > 1) {
    stop("`gate_alpha` must be between 0 and 1.", call. = FALSE)
  }
  if (!config_scalar_number(label_size) || label_size <= 0) {
    stop("`label_size` must be positive.", call. = FALSE)
  }
  if (!config_scalar_number(label_digits) || label_digits < 0 ||
      label_digits %% 1 != 0) {
    stop("`label_digits` must be a nonnegative integer.", call. = FALSE)
  }
  style <- if (inherits(appearance, "facs_appearance")) appearance else
    resolve_facs_appearance(analysis, appearance, appearance_file)
  style$show_phase_gates <- FALSE
  plots <- plot_pseudocolor_panels(analysis, appearance = style)
  summary <- analysis$quantitation$ph3$sample_summary

  for (prefix in names(plots$panels)) {
    gate <- geometry[geometry$prefix == prefix, , drop = FALSE]
    value <- summary$ph3_positive_percent[summary$prefix == prefix]
    if (length(value) != 1L || !is.finite(value)) {
      stop("Missing pH3-positive percentage for ", prefix, ".",
           call. = FALSE)
    }
    position <- ph3_gate_label_position(gate, style$y_log10)
    label <- paste0(
      "pH3+: ", formatC(value, format = "f", digits = label_digits), "%"
    )
    overlay <- list(
      ggplot2::geom_polygon(
        data = gate,
        ggplot2::aes(x = dna_norm, y = target_norm),
        inherit.aes = FALSE, fill = gate_fill, alpha = gate_alpha,
        color = style$gate_color, linetype = style$gate_linetype,
        linewidth = style$gate_linewidth
      ),
      ggplot2::annotate(
        "label", x = position[["x"]], y = position[["y"]], label = label,
        color = label_color, fill = scales::alpha("black", 0.55),
        linewidth = 0, size = label_size
      )
    )
    plots$panels[[prefix]] <- plots$panels[[prefix]] + overlay
    plots$decorated_plots[[prefix]] <-
      plots$decorated_plots[[prefix]] + overlay
    index <- match(prefix, analysis$sample_manifest$prefix)
    plots$panel_results[[index]]$panel <- plots$panels[[prefix]]
    plots$panel_results[[index]]$plot <- plots$decorated_plots[[prefix]]
  }
  plots$flowjo_gate_geometry <- geometry
  plots$annotation_denominator <- "all_single_cell_events"
  plots
}

#' Build an editable pH3 gate-focused figure bundle
#' @inheritParams plot_ph3_4n_gate_panels
#' @export
build_ph3_4n_figure_bundle <- function(
    analysis, geometry, appearance = NULL, appearance_file = NULL,
    gate_fill = "#FFFFFF", gate_alpha = 0.08,
    label_color = "#FFFFFF", label_size = 3.5, label_digits = 1L
) {
  resolved <- resolve_facs_appearance(analysis, appearance, appearance_file)
  layout_analysis <- analysis
  layout_analysis$config$layout <- "plotgardener"
  pseudocolor <- plot_ph3_4n_gate_panels(
    layout_analysis, geometry, resolved,
    gate_fill = gate_fill, gate_alpha = gate_alpha,
    label_color = label_color, label_size = label_size,
    label_digits = label_digits
  )
  pseudocolor$panel_results <- lapply(
    pseudocolor$panel_results,
    function(value) {
      value$sample_data <- NULL
      value
    }
  )
  structure(list(
    artifact_version = 1L,
    report_type = "ph3_4n_gate",
    pseudocolor = pseudocolor,
    flowjo_gate_geometry = geometry,
    plot_data = analysis$quantitation$ph3$sample_summary,
    appearance = resolved,
    sample_manifest = analysis$sample_manifest,
    provenance = analysis$provenance,
    analysis_config_path = attr(analysis$config, "config_path")
  ), class = "facs_figure_bundle")
}
