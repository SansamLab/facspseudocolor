# Pure normalization functions ----------------------------------------------

#' Normalize EdU signal and DNA content for one sample
#'
#' Applies the existing EdU normalization calculation to in-memory event
#' tables. The reference-derived slope is supplied explicitly; this function
#' does not read files, fit a reference sample, draw plots, or save output.
#'
#' @param events All single-cell events.
#' @param g1_events G1-gated events for the same sample.
#' @param dna_channel Exact DNA channel column name.
#' @param target_channel Exact EdU channel column name.
#' @param baseline_slope Finite slope fitted from the replicate's reference
#'   EdU-negative events.
#' @param dna_2n_value Numeric value to which the G1 DNA anchor is mapped.
#' @param g1_anchor Either `"median"` or `"mode"`.
#' @param edu_positive_events Optional EdU-positive events for the same sample.
#' @param sample_id Label used in validation errors and provenance.
#'
#' @return A list containing normalized all-cell and optional EdU-positive
#'   tables plus the fitted sample anchors.
#' @export
normalize_edu <- function(
    events, g1_events, dna_channel, target_channel, baseline_slope,
    dna_2n_value = 1000, g1_anchor = "median", edu_positive_events = NULL,
    sample_id = "sample"
) {
  events <- read_facs_sample(
    events, dna_channel, target_channel, paste(sample_id, "single cells"), 2L
  )
  g1_events <- read_facs_sample(
    g1_events, dna_channel, target_channel, paste(sample_id, "G1"), 2L
  )
  if (!config_scalar_number(baseline_slope)) {
    stop("`baseline_slope` must be one finite number for ", sample_id, ".",
         call. = FALSE)
  }
  if (!config_scalar_number(dna_2n_value) || dna_2n_value <= 0) {
    stop("`dna_2n_value` must be a positive finite number.", call. = FALSE)
  }
  if (!config_scalar_string(g1_anchor) || !g1_anchor %in% c("median", "mode")) {
    stop("`g1_anchor` must be 'median' or 'mode'.", call. = FALSE)
  }

  g1_median_dna <- stats::median(g1_events[[dna_channel]], na.rm = TRUE)
  g1_anchor_target <- g1_target_anchor(g1_events[[target_channel]], g1_anchor)
  if (!is.finite(g1_median_dna) || g1_median_dna == 0) {
    stop("Invalid G1 DNA median for ", sample_id, ".", call. = FALSE)
  }
  if (!is.finite(g1_anchor_target) || g1_anchor_target == 0) {
    stop("Invalid G1 target anchor for ", sample_id, ".", call. = FALSE)
  }

  normalize_table <- function(table) {
    table$target_raw <- table[[target_channel]]
    table$dna_norm <- table[[dna_channel]] / g1_median_dna * dna_2n_value
    table$baseline <- g1_anchor_target +
      baseline_slope * (table$dna_norm - dna_2n_value)
    table$target_norm <- table$target_raw / table$baseline * dna_2n_value
    table$target_bgsub <- table$target_raw - table$baseline
    table
  }

  normalized_events <- normalize_table(events)
  normalized_positive <- NULL
  if (!is.null(edu_positive_events)) {
    edu_positive_events <- read_facs_sample(
      edu_positive_events, dna_channel, target_channel,
      paste(sample_id, "EdU positive"), 2L
    )
    normalized_positive <- normalize_table(edu_positive_events)
  }

  list(
    data = normalized_events,
    edu_positive = normalized_positive,
    g1_median_dna = g1_median_dna,
    g1_anchor_target = g1_anchor_target,
    baseline_slope = baseline_slope,
    normalization_method = "reference_negative_regression"
  )
}

#' Normalize POI signal and DNA content for one sample
#'
#' Applies the existing protein-of-interest background regression to an
#' in-memory event table. The background model is supplied explicitly; this
#' function does not read files, fit the background, draw plots, or save output.
#'
#' @param events All single-cell events.
#' @param background_model Named list containing `dna_peak`, `intercept`,
#'   `slope`, and `floor`; `cutoff` is optional.
#' @param dna_channel Exact DNA channel column name.
#' @param target_channel Exact POI channel column name.
#' @param dna_2n_value Numeric value to which the G1 DNA peak is mapped.
#' @param dna_align Either `"per_sample"` or `"shared_background"`.
#' @param peak_failure What to do if per-sample DNA peak detection fails.
#'   `"error"` stops with an actionable error. `"use_background"` explicitly
#'   selects the legacy fallback to the background-control peak.
#' @param sample_id Label used in validation errors and provenance.
#'
#' @return A list containing the normalized event table, DNA peak used, cutoff,
#'   and normalization metadata.
#' @export
normalize_poi <- function(
    events, background_model, dna_channel, target_channel,
    dna_2n_value = 1000, dna_align = "per_sample", peak_failure = "error",
    sample_id = "sample"
) {
  events <- read_facs_sample(
    events, dna_channel, target_channel, paste(sample_id, "single cells"), 2L
  )
  required_model <- c("dna_peak", "intercept", "slope", "floor")
  if (!is.list(background_model) ||
      !all(required_model %in% names(background_model)) ||
      any(!vapply(background_model[required_model], config_scalar_number,
                  logical(1)))) {
    stop(
      "`background_model` must contain finite dna_peak, intercept, slope, and floor values for ",
      sample_id, ".", call. = FALSE
    )
  }
  if (background_model$dna_peak <= 0 || background_model$floor <= 0) {
    stop("Background dna_peak and floor must be positive for ", sample_id, ".",
         call. = FALSE)
  }
  if (!config_scalar_number(dna_2n_value) || dna_2n_value <= 0) {
    stop("`dna_2n_value` must be a positive finite number.", call. = FALSE)
  }
  dna_align <- match.arg(dna_align, c("per_sample", "shared_background"))
  peak_failure <- match.arg(peak_failure, c("error", "use_background"))

  if (dna_align == "shared_background") {
    sample_dna_peak <- background_model$dna_peak
  } else {
    peak_result <- tryCatch(
      estimate_lower_dna_peak(events[[dna_channel]]),
      error = identity
    )
    if (inherits(peak_result, "error")) {
      if (peak_failure == "use_background") {
        sample_dna_peak <- background_model$dna_peak
      } else {
        stop(
          "Per-sample DNA peak detection failed for ", sample_id, ": ",
          conditionMessage(peak_result),
          ". Inspect the DNA distribution or explicitly set ",
          "`poi_peak_failure: use_background` to use the legacy fallback.",
          call. = FALSE
        )
      }
    } else {
      sample_dna_peak <- peak_result
    }
  }

  events$target_raw <- events[[target_channel]]
  events$dna_norm <- events[[dna_channel]] / sample_dna_peak * dna_2n_value
  events$baseline <- pmax(
    background_model$intercept + background_model$slope * events$dna_norm,
    background_model$floor
  )
  events$target_norm <- dna_2n_value * events$target_raw / events$baseline
  events$target_bgsub <- events$target_raw - events$baseline

  list(
    data = events,
    dna_peak = sample_dna_peak,
    cutoff = get_setting(background_model, "cutoff", NA_real_),
    normalization_method = "background_reference_regression",
    dna_align = dna_align,
    peak_failure = peak_failure
  )
}
