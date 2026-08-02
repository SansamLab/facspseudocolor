# Interactive configurator ---------------------------------------------------

#' Set the data directory used by the interactive configurator
#'
#' Captures an explicit directory—by default the current R working directory—
#' in an environment variable inherited by a subsequently launched Quarto
#' configurator. Call this from the RStudio Console immediately before clicking
#' **Run Document** in `facs_configurator.qmd`.
#'
#' @param path Existing directory containing the exported sample CSV files.
#'
#' @return The normalized directory path, invisibly.
#' @export
set_facs_configurator_directory <- function(path = getwd()) {
  if (!is.character(path) || length(path) != 1L || is.na(path) ||
      !nzchar(path) || !dir.exists(path)) {
    stop("`path` must be one existing directory.", call. = FALSE)
  }
  path <- normalizePath(path, mustWork = TRUE)
  Sys.setenv(FACS_CONFIG_DATA_DIR = path)
  message("FACS configurator data directory: ", path)
  invisible(path)
}
