#' facspseudocolor: Analyze signal versus DNA flow cytometry data
#'
#' The package validates, normalizes, quantifies, and plots event-level flow
#' cytometry data for EdU incorporation or a protein of interest versus DNA
#' content. FlowJo workspace processing remains an optional, external Python
#' preprocessing step.
#'
#' @keywords internal
"_PACKAGE"

utils::globalVariables(c(
  ".data", "background_fitted", "boundary_x", "boundary_y", "condition",
  "density_color", "dna_norm", "error", "gate", "label_x", "label_y",
  "phase_label", "phase_percent", "target_norm", "target_raw", "xmax",
  "xmin", "ymax", "ymin"
))
