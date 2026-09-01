#' facspseudocolor: Analyze signal versus DNA flow cytometry data
#'
#' The package validates, normalizes, quantifies, and plots event-level flow
#' cytometry data for EdU incorporation, a protein of interest, or
#' phospho-histone H3 versus DNA content. FlowJo workspace processing remains
#' an optional, external Python preprocessing step.
#'
#' @keywords internal
"_PACKAGE"

utils::globalVariables(c(
  ".data", "background_fitted", "boundary_x", "boundary_y", "condition",
  "density_color", "dna_norm", "error", "g2m_percent", "gate", "gate_index",
  "label_x",
  "label_y", "ph3_phase", "phase_label", "phase_percent", "target_norm",
  "target_raw", "variant", "x", "xmax", "xmin", "y", "ymax", "ymin",
  "raw_pH3", "density", "condition_label", "value", "signal_basis",
  "mean_value", "sd_value"
))
