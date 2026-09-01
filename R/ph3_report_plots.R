# pH3 output-contract condition plots ---------------------------------------

ph3_report_plot_fail <- function(reason, detail) {
  stop(
    "PH3 report-plot validation failed [", reason, "]: ", detail, ".",
    call. = FALSE
  )
}

ph3_report_plot_outcome_ids <- function() c("A", "B", "C", "D")

ph3_report_plot_validate_model <- function(model) {
  report <- model$condition_report
  expected <- c(
    "schema_version", "biological_replicate_values", "condition_summary",
    "signal_basis_strata", "qc_flags", "provenance"
  )
  if (!inherits(model, "ph3_output_model") || !is.list(report) ||
      !identical(names(report), expected) ||
      !identical(report$schema_version, "ph3-condition-report-model-1.0.0") ||
      !identical(report$provenance$analysis_id, model$experiment$analysis_id[[1L]]) ||
      !identical(report$provenance$experiment_id, model$experiment$experiment_id[[1L]])) {
    ph3_report_plot_fail(
      "invalid_condition_report_model",
      "the validated Slice 4 condition report model is required"
    )
  }
  columns <- ph3_report_model_columns()
  ph3_report_model_validate_table(
    report$biological_replicate_values, columns$biological_replicate_values,
    "biological_replicate_values"
  )
  ph3_report_model_validate_table(
    report$condition_summary, columns$condition_summary, "condition_summary"
  )
  ph3_report_model_validate_table(report$qc_flags, columns$qc_flags, "qc_flags")
  outcomes <- model$outcomes
  if (!is.data.frame(outcomes) ||
      !identical(outcomes$outcome_id, ph3_report_plot_outcome_ids()) ||
      nrow(report$condition_summary) != nrow(model$conditions) * 4L ||
      anyDuplicated(report$condition_summary[c("condition_id", "outcome_id")])) {
    ph3_report_plot_fail(
      "invalid_condition_report_model",
      "the report must retain one A-D summary for every configured condition"
    )
  }
  invisible(report)
}

ph3_report_plot_colours <- function(conditions, supplied = NULL) {
  labels <- conditions$condition_label
  colours <- stats::setNames(facs_named_palette("colorblind", length(labels)), labels)
  if (!is.null(supplied)) {
    supplied <- unlist(supplied, use.names = TRUE)
    if (!is.character(supplied) || is.null(names(supplied)) ||
        any(!nzchar(names(supplied))) || any(!names(supplied) %in% labels)) {
      ph3_report_plot_fail(
        "invalid_condition_colours",
        "appearance colours must name only configured pH3 conditions"
      )
    }
    tryCatch(grDevices::col2rgb(unname(supplied)), error = function(error) {
      ph3_report_plot_fail("invalid_condition_colours", conditionMessage(error))
    })
    colours[names(supplied)] <- unname(supplied)
  }
  colours
}

ph3_report_plot_y_label <- function(outcome_id, value_mode, signal_basis) {
  if (outcome_id == "A") {
    return("4N pH3-positive cells (% of Analysis singlets, 2N-4N)")
  }
  if (outcome_id == "B") {
    return("Below-4N pH3-positive cells (% of Analysis singlets, 2N-4N)")
  }
  if (identical(signal_basis, "raw") && identical(value_mode, "reference_ratio")) {
    return("Median pH3 signal relative to matched reference (RAW\u2014NOT BACKGROUND SUBTRACTED)")
  }
  if (identical(value_mode, "reference_ratio")) {
    return("Median pH3 signal relative to matched reference")
  }
  if (identical(signal_basis, "raw")) {
    return("Median raw pH3 signal (RAW\u2014NOT BACKGROUND SUBTRACTED)")
  }
  "Median background-corrected pH3 signal"
}

ph3_report_plot_caption <- function(summary, points, raw_replicate_sets = character()) {
  unavailable <- summary[summary$summary_status %in% c(
    "unavailable_mixed_signal_basis", "unavailable_no_finite_values"
  ), , drop = FALSE]
  partial <- summary[summary$summary_status %in% c(
    "available_partial_biological_replicate_coverage",
    "available_partial_unavailable_cutoff_failure"
  ), , drop = FALSE]
  pieces <- c(
    "Points are biological-replicate values; bars show the condition mean. SD is shown only with at least three finite biological replicates."
  )
  if (nrow(partial)) {
    pieces <- c(pieces, paste0(
      "Partial biological-replicate coverage: ",
      paste(partial$condition_label, collapse = ", "), "."
    ))
  }
  if (nrow(unavailable)) {
    pieces <- c(pieces, paste0(
      "Unavailable condition summary: ",
      paste(paste0(unavailable$condition_label, " (", unavailable$reason_code, ")"),
            collapse = "; "), "."
    ))
  }
  bases <- unique(points$signal_basis[points$signal_basis != "not_applicable"])
  if (length(bases)) {
    pieces <- c(pieces, paste0("Signal basis: ", paste(bases, collapse = ", "), "."))
  }
  raw_sets <- unique(raw_replicate_sets)
  if (length(raw_sets)) {
    pieces <- c(pieces, paste0(
      "RAW\u2014NOT BACKGROUND SUBTRACTED for replicate set(s): ",
      paste(raw_sets, collapse = ", "), "."
    ))
  }
  paste(pieces, collapse = " ")
}

ph3_report_plot_panel <- function(report, outcomes, outcome_id, colours) {
  outcome <- outcomes[outcomes$outcome_id == outcome_id, , drop = FALSE]
  summary <- report$condition_summary[
    report$condition_summary$outcome_id == outcome_id, , drop = FALSE
  ]
  source_points <- report$biological_replicate_values[
    report$biological_replicate_values$outcome_id == outcome_id, , drop = FALSE
  ]
  points <- source_points[is.finite(source_points$value), , drop = FALSE]
  if (nrow(outcome) != 1L || nrow(summary) == 0L) {
    ph3_report_plot_fail("missing_panel_outcome", "every A-D panel requires source rows")
  }
  summary <- summary[order(summary$condition_index), , drop = FALSE]
  points <- points[order(points$condition_index, points$replicate_index), , drop = FALSE]
  display_basis <- unique(source_points$signal_basis[
    source_points$signal_basis != "not_applicable"
  ])
  if (!length(display_basis) && outcome_id %in% c("A", "B")) {
    display_basis <- "not_applicable"
  }
  y_label <- ph3_report_plot_y_label(
    outcome_id, summary$value_mode[[1L]],
    if ("raw" %in% display_basis) "raw" else if (length(display_basis) == 1L) {
      display_basis[[1L]]
    } else "mixed"
  )
  mean_rows <- summary[is.finite(summary$mean_value), , drop = FALSE]
  point_mapping <- if (outcome_id %in% c("C", "D")) {
    ggplot2::aes(x = condition_label, y = value, colour = condition_label,
                 shape = signal_basis)
  } else {
    ggplot2::aes(x = condition_label, y = value, colour = condition_label)
  }
  plot <- ggplot2::ggplot(points, point_mapping) +
    ggplot2::geom_point(
      position = ggplot2::position_jitter(width = 0.08, height = 0, seed = 1L),
      size = 2.4, na.rm = TRUE
    ) +
    ggplot2::geom_errorbar(
      data = mean_rows,
      ggplot2::aes(
        x = condition_label, ymin = mean_value - sd_value,
        ymax = mean_value + sd_value, colour = condition_label
      ),
      width = 0.14, na.rm = TRUE, inherit.aes = FALSE
    ) +
    ggplot2::geom_point(
      data = mean_rows,
      ggplot2::aes(x = condition_label, y = mean_value, colour = condition_label),
      size = 4.0, shape = 95, stroke = 1.2, inherit.aes = FALSE
    ) +
    ggplot2::scale_colour_manual(values = colours, name = "Condition") +
    ggplot2::labs(
      title = outcome$label[[1L]], x = NULL, y = y_label,
      caption = ph3_report_plot_caption(
        summary, points,
        source_points$replicate_set_id[source_points$signal_basis == "raw"]
      )
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      plot.caption = ggplot2::element_text(hjust = 0, size = 8),
      axis.text.x = ggplot2::element_text(angle = 20, hjust = 1)
    )
  if (outcome_id %in% c("C", "D")) {
    plot <- plot + ggplot2::guides(shape = ggplot2::guide_legend(title = "Signal basis"))
  }
  plot
}

ph3_raw_4n_cutoff_qc_plots <- function(cutoff) {
  required <- c("schema_version", "method", "records", "application_qc", "density_curves")
  if (is.null(cutoff)) return(list())
  if (!is.list(cutoff) || !identical(names(cutoff), required) ||
      !identical(cutoff$schema_version, "ph3-raw-4n-density-cutoff-1.0.0") ||
      !is.data.frame(cutoff$records) || !is.data.frame(cutoff$density_curves)) {
    ph3_report_plot_fail("invalid_cutoff_provenance",
                         "computed positivity requires valid raw-4N cutoff provenance")
  }
  stats::setNames(lapply(seq_len(nrow(cutoff$records)), function(i) {
    record <- cutoff$records[i, , drop = FALSE]
    curve <- cutoff$density_curves[
      cutoff$density_curves$replicate_set_id == record$replicate_set_id[[1L]], , drop = FALSE
    ]
    if (identical(record$cutoff_status[[1L]], "unavailable_cutoff_failure")) {
      return(ggplot2::ggplot() + ggplot2::annotate(
        "text", x = 0, y = 0,
        label = paste0("No cutoff: ", record$cutoff_reason_code[[1L]]),
        hjust = 0
      ) + ggplot2::labs(
        title = paste0("Raw 4N control pH3 cutoff: ", record$replicate_set_id[[1L]]),
        subtitle = "Cutoff unavailable; this replicate set is excluded from condition means",
        x = "Raw pH3 fluorescence", y = "Density"
      ) + ggplot2::theme_void())
    }
    if (!nrow(curve) || any(!is.finite(curve$raw_pH3)) || any(!is.finite(curve$density))) {
      ph3_report_plot_fail("invalid_cutoff_density_curve",
                           "each cutoff record requires a finite raw-signal density curve")
    }
    ggplot2::ggplot(curve, ggplot2::aes(x = raw_pH3, y = density)) +
      ggplot2::geom_line() +
      ggplot2::geom_vline(xintercept = record$cutoff_raw_signal[[1L]],
                          linetype = "dashed", colour = "firebrick") +
      ggplot2::geom_point(data = data.frame(raw_pH3 = record$dominant_peak_raw_signal,
                                            density = record$dominant_peak_density),
                          ggplot2::aes(x = raw_pH3, y = density), inherit.aes = FALSE) +
      ggplot2::labs(
        title = paste0("Raw 4N control pH3 cutoff: ", record$replicate_set_id[[1L]]),
        subtitle = "Gaussian KDE (nrd0; adjust 1); dashed line = first right-side local minimum",
        x = "Raw pH3 fluorescence", y = "Density"
      ) + ggplot2::theme_classic(base_size = 11)
  }), cutoff$records$replicate_set_id)
}

#' Build the four owner-confirmed pH3 condition plots
#'
#' This function consumes the validated Slice 4 report model only. It does not
#' recalculate events, gates, thresholds, background correction, reference
#' ratios, or condition summaries, and writes no files.
#'
#' @param analysis A completed direct-identity pH3 `facs_analysis`.
#' @param appearance Optional presentation-only appearance overrides.
#' @param appearance_file Optional YAML appearance file.
#' @return A named list containing four editable ggplots and the validated
#'   report-model tables used to render them.
#' @export
plot_ph3_condition_report <- function(
    analysis, appearance = NULL, appearance_file = NULL
) {
  validate_analysis_object(analysis)
  if (!identical(analysis$config$plot_type, "ph3") ||
      !inherits(analysis$ph3_output_model, "ph3_output_model")) {
    ph3_report_plot_fail(
      "invalid_analysis", "one completed production pH3 output-contract analysis is required"
    )
  }
  model <- analysis$ph3_output_model
  report <- ph3_report_plot_validate_model(model)
  resolved <- resolve_facs_appearance(analysis, appearance, appearance_file)
  colours <- ph3_report_plot_colours(model$conditions, resolved$condition_colors)
  panels <- stats::setNames(lapply(ph3_report_plot_outcome_ids(), function(id) {
    ph3_report_plot_panel(report, model$outcomes, id, colours)
  }), c(
    "ph3_4n_positive_prevalence",
    "ph3_below_4n_positive_prevalence",
    "ph3_4n_positive_signal",
    "ph3_below_4n_positive_signal"
  ))
  cutoff <- model$source$raw_4n_density_cutoff
  cutoff_qc <- ph3_raw_4n_cutoff_qc_plots(cutoff)
  list(
    schema_version = "ph3-condition-report-plots-1.0.0",
    panels = panels,
    biological_replicate_values = report$biological_replicate_values,
    condition_summary = report$condition_summary,
    signal_basis_strata = report$signal_basis_strata,
    qc_flags = report$qc_flags,
    provenance = report$provenance,
    cutoff_qc = cutoff_qc,
    cutoff_records = if (is.null(cutoff)) NULL else cutoff$records,
    cutoff_application_qc = if (is.null(cutoff)) NULL else cutoff$application_qc
  )
}
