# Durable analysis and figure artifacts -------------------------------------

require_plotgardener <- function(
    available = requireNamespace("plotgardener", quietly = TRUE)
) {
  if (!isTRUE(available)) {
    stop(
      "This figure layout requires the optional plotgardener package. ",
      "Install it with BiocManager::install('plotgardener').",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

artifact_write <- function(object, path, overwrite = FALSE) {
  if (!config_scalar_string(path)) stop("`path` must be explicit.", call. = FALSE)
  path <- normalizePath(path.expand(path), mustWork = FALSE)
  if (file.exists(path) && !isTRUE(overwrite)) {
    stop("Refusing to overwrite existing file: ", path, call. = FALSE)
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(object, path)
  invisible(path)
}

#' Save or read a complete FACS analysis artifact
#' @param analysis A `facs_analysis` object.
#' @param path Explicit RDS path.
#' @param overwrite Whether an existing file may be replaced.
#' @export
save_facs_analysis <- function(analysis, path, overwrite = FALSE) {
  validate_analysis_object(analysis)
  artifact_write(analysis, path, overwrite)
}

#' @rdname save_facs_analysis
#' @export
read_facs_analysis <- function(path) {
  if (!config_scalar_string(path) || !file.exists(path)) {
    stop("Analysis RDS not found: ", path, call. = FALSE)
  }
  value <- readRDS(path)
  validate_analysis_object(value)
  value
}

#' Build a reusable collection of editable FACS plots
#'
#' Event-level data remain in the analysis object. The bundle retains editable
#' ggplots, plot-ready summary tables, dynamic layout metadata, appearance, and
#' provenance needed to rebuild Plotgardener pages.
#' @param analysis A complete `facs_analysis` object.
#' @param appearance Optional presentation-only overrides.
#' @param appearance_file Optional YAML appearance file.
#' @param seed Explicit diagnostic downsampling seed.
#' @export
build_facs_figure_bundle <- function(
    analysis, appearance = NULL, appearance_file = NULL, seed = 1L
) {
  require_plotgardener()
  validate_analysis_object(analysis)
  resolved <- resolve_facs_appearance(analysis, appearance, appearance_file)
  layout_analysis <- analysis
  layout_analysis$config$layout <- "plotgardener"
  pseudocolor <- plot_pseudocolor_panels(layout_analysis, appearance = resolved)
  pseudocolor$panel_results <- lapply(pseudocolor$panel_results, function(value) {
    value$sample_data <- NULL
    value
  })
  quantitation <- plot_facs_quantitation(analysis, seed = seed,
                                         appearance = resolved)
  quantitation_page <- assemble_facs_panels(
    quantitation,
    ncol = min(2L, max(1L, length(quantitation))),
    layout_options = resolved$layout_options
  )
  structure(list(
    artifact_version = 1L,
    pseudocolor = pseudocolor,
    quantitation = quantitation,
    quantitation_page = quantitation_page,
    plot_data = analysis$quantitation,
    appearance = resolved,
    sample_manifest = analysis$sample_manifest,
    provenance = analysis$provenance,
    analysis_config_path = attr(analysis$config, "config_path")
  ), class = "facs_figure_bundle")
}

#' Assemble editable ggplots on a dynamic Plotgardener page
#' @param plots A named list of ggplot objects.
#' @param ncol Number of columns.
#' @param layout_options Optional Plotgardener sizing settings.
#' @export
assemble_facs_panels <- function(plots, ncol = NULL, layout_options = list()) {
  require_plotgardener()
  if (inherits(plots, "ggplot")) plots <- list(plots)
  if (!is.list(plots) || !length(plots) ||
      any(!vapply(plots, inherits, logical(1), "ggplot"))) {
    stop("`plots` must be a nonempty list of ggplot objects.", call. = FALSE)
  }
  if (is.null(ncol)) ncol <- ceiling(sqrt(length(plots)))
  if (!is.numeric(ncol) || length(ncol) != 1L || ncol < 1 || ncol %% 1 != 0) {
    stop("`ncol` must be one positive integer.", call. = FALSE)
  }
  nrow <- ceiling(length(plots) / ncol)
  generic_defaults <- list(
    panel_size = 3.2, panel_gap_x = 0.35, panel_gap_y = 0.35,
    margin_left = 0.25, margin_top = 0.25,
    margin_right = 0.25, margin_bottom = 0.25
  )
  layout_options <- utils::modifyList(generic_defaults, layout_options)
  layout <- facs_page_layout(ncol, nrow, layout_options)
  structure(list(plots = plots, ncol = as.integer(ncol), nrow = nrow,
                 layout = layout), class = "facs_panel_bundle")
}

#' @export
print.facs_panel_bundle <- function(x, ...) {
  require_plotgardener()
  plotgardener::pageCreate(width = x$layout$width, height = x$layout$height,
                           default.units = "inches", showGuides = FALSE,
                           xgrid = 0, ygrid = 0)
  for (i in seq_along(x$plots)) {
    row <- ceiling(i / x$ncol)
    column <- ((i - 1L) %% x$ncol) + 1L
    plotgardener::plotGG(
      plot = x$plots[[i]], x = x$layout$panel_x(column),
      y = x$layout$panel_y(row), width = x$layout$panel,
      height = x$layout$panel, just = c("left", "top"),
      default.units = "inches"
    )
  }
  invisible(x)
}

#' Save or read a FACS figure bundle
#' @param bundle A `facs_figure_bundle` object.
#' @param path Explicit RDS path.
#' @param overwrite Whether an existing file may be replaced.
#' @export
save_facs_figure_bundle <- function(bundle, path, overwrite = FALSE) {
  if (!inherits(bundle, "facs_figure_bundle")) {
    stop("`bundle` must be a facs_figure_bundle object.", call. = FALSE)
  }
  artifact_write(bundle, path, overwrite)
}

#' @rdname save_facs_figure_bundle
#' @export
read_facs_figure_bundle <- function(path) {
  if (!config_scalar_string(path) || !file.exists(path)) {
    stop("Figure-bundle RDS not found: ", path, call. = FALSE)
  }
  value <- readRDS(path)
  if (!inherits(value, "facs_figure_bundle")) {
    stop("RDS does not contain a facs_figure_bundle object.", call. = FALSE)
  }
  value
}

#' Retrieve one editable plot or its plot-ready data
#' @param bundle A `facs_figure_bundle`.
#' @param name Plot name or pseudocolor sample prefix.
#' @export
get_facs_panel <- function(bundle, name) {
  if (!inherits(bundle, "facs_figure_bundle")) stop("Invalid figure bundle.", call. = FALSE)
  candidates <- c(bundle$pseudocolor$decorated_plots, bundle$quantitation)
  if (!name %in% names(candidates)) stop("Plot not found in bundle: ", name, call. = FALSE)
  candidates[[name]]
}

#' @rdname get_facs_panel
#' @export
get_facs_plot_data <- function(bundle, name = NULL) {
  if (!inherits(bundle, "facs_figure_bundle")) stop("Invalid figure bundle.", call. = FALSE)
  if (is.null(name)) return(bundle$plot_data)
  if (!name %in% names(bundle$plot_data)) stop("Plot data not found: ", name, call. = FALSE)
  bundle$plot_data[[name]]
}
