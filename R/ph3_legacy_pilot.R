# Legacy CSV pH3 pilot -------------------------------------------------------

# This deliberately separate path permits a real-data smoke test of the pH3
# report using historical CSV exports.  It is never a substitute for the
# direct-identity production profile: no event identity or gate containment is
# claimed, and every output carries the limited-provenance label.

ph3_pilot_fail <- function(reason, detail) {
  stop("PH3 legacy CSV pilot failed [", reason, "]: ", detail, ".", call. = FALSE)
}

ph3_pilot_label <- function(analysis) analysis$config$ph3_pilot$provenance_label

# The legacy pilot has no direct event identity, so the explicitly configured
# manifest prefix is the only permitted sample-to-table join. Do not infer
# this relationship from list position: a reordered list could otherwise
# attach a different sample's events to a clone or condition.
ph3_pilot_validate_normalized_data_alignment <- function(analysis, manifest) {
  expected <- as.character(manifest$prefix)
  actual <- names(analysis$normalized_data)
  if (is.null(actual) || !identical(actual, expected)) {
    ph3_pilot_fail(
      "normalized_data_prefix_alignment",
      "normalized-data names must exactly equal configured manifest prefixes in manifest order"
    )
  }
  invisible(TRUE)
}

ph3_pilot_regions <- function(config, data) {
  lower <- as.numeric(unlist(config$g1_x_range))[[1L]]
  four <- as.numeric(unlist(config$g2m_x_range))
  eligible <- is.finite(data$dna_norm) & data$dna_norm >= lower & data$dna_norm <= four[[2L]]
  list(eligible = eligible, four_n = eligible & data$dna_norm >= four[[1L]],
       below_4n = eligible & data$dna_norm < four[[1L]])
}

ph3_pilot_value <- function(x, denominator) {
  if (!length(denominator) || !any(denominator)) return(NA_real_)
  100 * sum(x, na.rm = TRUE) / sum(denominator)
}

# Display-only metadata for background-subtracted pH3-versus-DNA panels. These rows point back
# to the normalized complete-event tables that the pilot already used; they do
# not recreate or claim any FlowJo gate geometry.
ph3_pilot_pseudocolor_metadata <- function(analysis, cutoff_records,
                                            display_offset_qc, log_display_qc) {
  manifest <- analysis$sample_manifest
  records <- cutoff_records[match(manifest$replicate_index,
                                  cutoff_records$replicate_index), , drop = FALSE]
  if (nrow(records) != nrow(manifest) ||
      any(records$replicate_index != manifest$replicate_index) ||
      anyDuplicated(manifest$prefix)) {
    ph3_pilot_fail("invalid_cutoff_records",
                   "every configured pilot sample must map to one clone cutoff record")
  }
  offsets <- display_offset_qc[match(manifest$prefix, display_offset_qc$prefix), , drop = FALSE]
  if (nrow(offsets) != nrow(manifest) || anyNA(match(manifest$prefix, offsets$prefix)) ||
      !identical(offsets$prefix, manifest$prefix)) {
    ph3_pilot_fail("invalid_display_offset_qc",
                   "every configured pilot sample must map to one display-offset QC row")
  }
  log_display <- log_display_qc[match(manifest$prefix, log_display_qc$prefix), , drop = FALSE]
  if (nrow(log_display) != nrow(manifest) || anyNA(match(manifest$prefix, log_display$prefix)) ||
      !identical(log_display$prefix, manifest$prefix)) {
    ph3_pilot_fail("invalid_log_display_qc",
                   "every configured pilot sample must map to one log-display QC row")
  }
  data.frame(
    prefix = manifest$prefix,
    replicate = manifest$replicate,
    replicate_index = as.integer(manifest$replicate_index),
    condition = manifest$condition,
    condition_index = as.integer(manifest$condition_index),
    cutoff_control_prefix = records$control_prefix,
    cutoff_corrected_signal = records$cutoff_corrected_signal,
    cutoff_status = records$cutoff_status,
    cutoff_reason_code = records$reason_code,
    display_offset_status = offsets$display_offset_status,
    display_offset_reason_code = offsets$display_offset_reason_code,
    shared_raw_negative_median = offsets$shared_raw_negative_median,
    shared_target_status = offsets$shared_target_status,
    shared_target_reason_code = offsets$shared_target_reason_code,
    display_offset = offsets$display_offset,
    displayed_cutoff_signal = offsets$displayed_cutoff_signal,
    log_display_status = log_display$log_display_status,
    log_display_reason_code = log_display$log_display_reason_code,
    visual_window_finite_event_count = log_display$visual_window_finite_event_count,
    positive_domain_display_event_count = log_display$positive_domain_display_event_count,
    nonpositive_display_event_count = log_display$nonpositive_display_event_count,
    signal_basis = "background_subtracted",
    provenance_label = ph3_pilot_label(analysis),
    stringsAsFactors = FALSE
  )
}

# This is a presentation-only audit. It runs after the owner-approved display
# offset and records exactly which finite, visual-window events cannot be shown
# on a log10 axis. It never changes corrected event values, cutoffs, or outcomes.
ph3_pilot_log_display_qc <- function(analysis, display_offset_qc, corrections) {
  manifest <- analysis$sample_manifest
  x_limits <- c(0.8, 2.2) * as.numeric(analysis$config$dna_2n_value)
  rows <- lapply(seq_len(nrow(manifest)), function(i) {
    prefix <- manifest$prefix[[i]]
    offset <- display_offset_qc[display_offset_qc$prefix == prefix, , drop = FALSE]
    base <- data.frame(
      prefix = prefix,
      replicate = manifest$replicate[[i]],
      replicate_index = as.integer(manifest$replicate_index[[i]]),
      condition = manifest$condition[[i]],
      visual_window_finite_event_count = NA_integer_,
      positive_domain_display_event_count = NA_integer_,
      nonpositive_display_event_count = NA_integer_,
      log_display_status = "unavailable",
      log_display_reason_code = NA_character_,
      stringsAsFactors = FALSE
    )
    if (nrow(offset) != 1L || !identical(offset$display_offset_status[[1L]], "available") ||
        !is.finite(offset$display_offset[[1L]])) {
      base$log_display_reason_code <- "display_offset_unavailable"
      return(base)
    }
    sample <- analysis$normalized_data[[prefix]]
    corrected <- corrections[[prefix]]
    if (is.null(sample) || !is.data.frame(sample$data) ||
        !all(c("dna_norm", "target_raw") %in% names(sample$data)) ||
        !is.data.frame(corrected) || !identical(corrected$event_row, seq_len(nrow(sample$data))) ||
        !all(c("corrected_signal", "correction_status") %in% names(corrected)) ||
        !all(corrected$correction_status == "available") ||
        any(!is.finite(corrected$corrected_signal))) {
      base$log_display_reason_code <- "background_correction_unavailable"
      return(base)
    }
    in_window <- is.finite(sample$data$dna_norm) &
      sample$data$dna_norm >= x_limits[[1L]] & sample$data$dna_norm <= x_limits[[2L]]
    displayed <- corrected$corrected_signal[in_window] + offset$display_offset[[1L]]
    # corrected_signal and the offset were already required finite; retain the
    # check as a fail-visible guard against an unexpected arithmetic result.
    if (any(!is.finite(displayed))) {
      base$log_display_reason_code <- "nonfinite_display_signal"
      return(base)
    }
    total <- length(displayed)
    positive <- sum(displayed > 0)
    nonpositive <- sum(displayed <= 0)
    base$visual_window_finite_event_count <- as.integer(total)
    base$positive_domain_display_event_count <- as.integer(positive)
    base$nonpositive_display_event_count <- as.integer(nonpositive)
    if (!positive) {
      base$log_display_reason_code <- "no_positive_domain_display_events"
    } else if (nonpositive) {
      base$log_display_status <- "available_nonpositive_display_values_excluded"
      base$log_display_reason_code <- "nonpositive_display_values_excluded"
    } else {
      base$log_display_status <- "available"
      base$log_display_reason_code <- NA_character_
    }
    base
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

ph3_pilot_display_offset_qc <- function(analysis, cutoff_records, corrections) {
  manifest <- analysis$sample_manifest
  # First retain each sample's independently auditable final-negative medians.
  # A clone's display target is deliberately not calculated from pooled events:
  # every configured sample contributes one median with equal weight.
  local_rows <- lapply(seq_len(nrow(manifest)), function(i) {
    prefix <- manifest$prefix[[i]]
    cutoff <- cutoff_records[cutoff_records$replicate_index == manifest$replicate_index[[i]], , drop = FALSE]
    base <- data.frame(
      prefix = prefix,
      replicate = manifest$replicate[[i]],
      replicate_index = as.integer(manifest$replicate_index[[i]]),
      condition = manifest$condition[[i]],
      cutoff_control_prefix = if (nrow(cutoff) == 1L) cutoff$control_prefix[[1L]] else NA_character_,
      cutoff_corrected_signal = if (nrow(cutoff) == 1L) cutoff$cutoff_corrected_signal[[1L]] else NA_real_,
      eligible_final_negative_event_count = NA_integer_,
      raw_negative_median = NA_real_,
      corrected_negative_median = NA_real_,
      shared_raw_negative_median = NA_real_,
      shared_target_status = "unavailable",
      shared_target_reason_code = NA_character_,
      display_offset = NA_real_,
      displayed_cutoff_signal = NA_real_,
      display_offset_status = "unavailable",
      display_offset_reason_code = NA_character_,
      stringsAsFactors = FALSE
    )
    if (nrow(cutoff) != 1L || !identical(cutoff$cutoff_status[[1L]], "available") ||
        !is.finite(cutoff$cutoff_corrected_signal[[1L]])) {
      base$display_offset_reason_code <- "final_corrected_cutoff_unavailable"
      return(base)
    }
    corrected <- corrections[[prefix]]
    if (is.null(corrected) || !is.data.frame(corrected) ||
        !all(c("dna_norm", "target_raw", "corrected_signal", "correction_status") %in% names(corrected)) ||
        !all(corrected$correction_status == "available") ||
        any(!is.finite(corrected$dna_norm)) || any(!is.finite(corrected$target_raw)) ||
        any(!is.finite(corrected$corrected_signal))) {
      base$display_offset_reason_code <- "background_correction_unavailable"
      return(base)
    }
    regions <- ph3_pilot_regions(analysis$config, corrected)
    final_negative <- regions$eligible &
      corrected$corrected_signal <= cutoff$cutoff_corrected_signal[[1L]]
    count <- sum(final_negative)
    base$eligible_final_negative_event_count <- as.integer(count)
    if (!count) {
      base$display_offset_reason_code <- "no_eligible_final_negative_events"
      return(base)
    }
    raw_median <- stats::median(corrected$target_raw[final_negative])
    corrected_median <- stats::median(corrected$corrected_signal[final_negative])
    if (!is.finite(raw_median) || !is.finite(corrected_median)) {
      base$display_offset_reason_code <- "nonfinite_final_negative_display_offset"
      return(base)
    }
    base$raw_negative_median <- raw_median
    base$corrected_negative_median <- corrected_median
    base
  })
  result <- do.call(rbind, local_rows)
  rownames(result) <- NULL

  for (replicate_index in unique(manifest$replicate_index)) {
    members <- which(result$replicate_index == replicate_index)
    # The required group is the complete explicit manifest group.  There is no
    # partial-group target or fallback to a single usable sample.
    local_valid <- result$display_offset_reason_code[members] == "" |
      is.na(result$display_offset_reason_code[members])
    local_valid <- local_valid & is.finite(result$raw_negative_median[members]) &
      is.finite(result$corrected_negative_median[members])
    if (!all(local_valid)) {
      result$shared_target_status[members] <- "unavailable"
      result$shared_target_reason_code[members] <- "configured_group_member_display_input_unavailable"
      result$display_offset_status[members] <- "unavailable"
      result$display_offset_reason_code[members] <- "clone_group_shared_raw_negative_median_unavailable"
      next
    }
    # stats::median is the owner-approved deterministic equal-sample rule,
    # including its standard arithmetic mean for an even number of samples.
    target <- stats::median(result$raw_negative_median[members])
    offsets <- target - result$corrected_negative_median[members]
    displayed_cutoffs <- result$cutoff_corrected_signal[members] + offsets
    if (!is.finite(target) || any(!is.finite(offsets)) || any(!is.finite(displayed_cutoffs))) {
      result$shared_target_status[members] <- "unavailable"
      result$shared_target_reason_code[members] <- "nonfinite_group_shared_raw_negative_median"
      result$display_offset_status[members] <- "unavailable"
      result$display_offset_reason_code[members] <- "clone_group_shared_raw_negative_median_unavailable"
      next
    }
    result$shared_raw_negative_median[members] <- target
    result$shared_target_status[members] <- "available"
    result$shared_target_reason_code[members] <- NA_character_
    result$display_offset[members] <- offsets
    result$displayed_cutoff_signal[members] <- displayed_cutoffs
    result$display_offset_status[members] <- "available"
    result$display_offset_reason_code[members] <- NA_character_
  }
  rownames(result) <- NULL
  result
}

ph3_pilot_correction_qc_row <- function(analysis, prefix) {
  qc <- analysis$ph3_legacy_pilot$provisional_background_qc
  required <- c("prefix", "correction_status", "correction_reason_code")
  if (!is.data.frame(qc) || !all(required %in% names(qc)) ||
      anyDuplicated(qc$prefix) || !identical(qc$prefix, analysis$sample_manifest$prefix)) {
    ph3_pilot_fail("invalid_background_correction_qc",
                   "one ordered background-correction QC row is required for every configured pilot sample")
  }
  row <- qc[qc$prefix == prefix, , drop = FALSE]
  if (nrow(row) != 1L || !row$correction_status[[1L]] %in% c("available", "unavailable") ||
      !is.character(row$correction_reason_code) || length(row$correction_reason_code) != 1L) {
    ph3_pilot_fail("invalid_background_correction_qc",
                   paste0("background-correction QC is invalid for ", prefix))
  }
  row
}

ph3_pilot_unavailable_correction_panel <- function(analysis, metadata_row, correction_qc) {
  prefix <- metadata_row$prefix[[1L]]
  reason <- correction_qc$correction_reason_code[[1L]]
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0.5, y = 0.62, label = "BACKGROUND CORRECTION UNAVAILABLE",
                      fontface = "bold", colour = "#A00000", size = 3.2) +
    ggplot2::annotate("text", x = 0.5, y = 0.49,
                      label = paste0("Status: ", correction_qc$correction_status[[1L]],
                                     "\nReason: ", reason),
                      size = 2.6) +
    ggplot2::annotate("text", x = 0.5, y = 0.32,
                      label = "No raw events, corrected density, or corrected cutoff are shown.\nPILOT / LIMITED-PROVENANCE \u2014 NOT PUBLICATION-GRADE.",
                      size = 2.3) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE, clip = "off") +
    ggplot2::labs(
      title = paste0(metadata_row$replicate[[1L]], " \u2014 ", metadata_row$condition[[1L]]),
      subtitle = paste0("Sample: ", prefix, ". Background correction unavailable; clone outcomes are unavailable."),
      x = "Normalized DNA content", y = "Background-subtracted pH3 fluorescence",
      caption = paste0(metadata_row$provenance_label[[1L]],
                       ". No corrected event display is fabricated when correction is unavailable.")
    ) +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(
      aspect.ratio = 1,
      plot.caption = ggplot2::element_text(hjust = 0, size = 6),
      plot.title = ggplot2::element_text(size = 8),
      plot.subtitle = ggplot2::element_text(size = 6),
      axis.title = ggplot2::element_text(size = 7),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    )
}

ph3_pilot_unavailable_display_offset_panel <- function(analysis, metadata_row) {
  prefix <- metadata_row$prefix[[1L]]
  reason <- metadata_row$display_offset_reason_code[[1L]]
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0.5, y = 0.62, label = "CORRECTED DISPLAY OFFSET UNAVAILABLE",
                      fontface = "bold", colour = "#A00000", size = 3.2) +
    ggplot2::annotate("text", x = 0.5, y = 0.49,
                      label = paste0("Status: ", metadata_row$display_offset_status[[1L]],
                                     "\nReason: ", reason), size = 2.6) +
    ggplot2::annotate("text", x = 0.5, y = 0.32,
                      label = "No raw events, corrected density, or cutoff are shown.\nPILOT / LIMITED-PROVENANCE \u2014 NOT PUBLICATION-GRADE.",
                      size = 2.3) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE, clip = "off") +
    ggplot2::labs(
      title = paste0(metadata_row$replicate[[1L]], " \u2014 ", metadata_row$condition[[1L]]),
      subtitle = paste0("Sample: ", prefix, ". Corrected display offset unavailable; clone outcomes are unavailable."),
      x = "Normalized DNA content", y = "Background-subtracted pH3 fluorescence (display offset restored)",
      caption = paste0(metadata_row$provenance_label[[1L]],
                       ". No corrected event display is fabricated when display restoration is unavailable.")
    ) +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(
      aspect.ratio = 1,
      plot.caption = ggplot2::element_text(hjust = 0, size = 6),
      plot.title = ggplot2::element_text(size = 8),
      plot.subtitle = ggplot2::element_text(size = 6),
      axis.title = ggplot2::element_text(size = 7),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    )
}

ph3_pilot_unavailable_log_display_panel <- function(analysis, metadata_row) {
  prefix <- metadata_row$prefix[[1L]]
  reason <- metadata_row$log_display_reason_code[[1L]]
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0.5, y = 0.62, label = "LOG10 DISPLAY UNAVAILABLE",
                      fontface = "bold", colour = "#A00000", size = 3.2) +
    ggplot2::annotate("text", x = 0.5, y = 0.49,
                      label = paste0("Status: ", metadata_row$log_display_status[[1L]],
                                     "\nReason: ", reason), size = 2.6) +
    ggplot2::annotate("text", x = 0.5, y = 0.32,
                      label = "No positive-domain display events are available.\nNo raw fallback, density, or cutoff overlay is shown.",
                      size = 2.3) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE, clip = "off") +
    ggplot2::labs(
      title = paste0(metadata_row$replicate[[1L]], " \u2014 ", metadata_row$condition[[1L]]),
      subtitle = paste0("Sample: ", prefix, ". Background-subtracted with display offset; log10 display unavailable."),
      x = "Normalized DNA content", y = "Background-subtracted pH3 fluorescence (display offset restored; log10 display)",
      caption = paste0(metadata_row$provenance_label[[1L]],
                       ". No corrected event display is fabricated when no positive log-domain values exist.")
    ) +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(
      aspect.ratio = 1,
      plot.caption = ggplot2::element_text(hjust = 0, size = 6),
      plot.title = ggplot2::element_text(size = 8),
      plot.subtitle = ggplot2::element_text(size = 6),
      axis.title = ggplot2::element_text(size = 7),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    )
}

ph3_pilot_pseudocolor_panel <- function(analysis, metadata_row) {
  prefix <- metadata_row$prefix[[1L]]
  correction_qc <- ph3_pilot_correction_qc_row(analysis, prefix)
  if (identical(correction_qc$correction_status[[1L]], "unavailable")) {
    return(ph3_pilot_unavailable_correction_panel(analysis, metadata_row, correction_qc))
  }
  if (!identical(metadata_row$cutoff_status[[1L]], "available") ||
      !is.finite(metadata_row$cutoff_corrected_signal[[1L]]) ||
      !identical(metadata_row$display_offset_status[[1L]], "available") ||
      !is.finite(metadata_row$display_offset[[1L]]) ||
      !is.finite(metadata_row$displayed_cutoff_signal[[1L]])) {
    return(ph3_pilot_unavailable_display_offset_panel(analysis, metadata_row))
  }
  if (!metadata_row$log_display_status[[1L]] %in% c("available", "available_nonpositive_display_values_excluded")) {
    return(ph3_pilot_unavailable_log_display_panel(analysis, metadata_row))
  }
  sample <- analysis$normalized_data[[prefix]]
  if (is.null(sample) || !is.data.frame(sample$data) ||
      !all(c("dna_norm", "target_raw") %in% names(sample$data))) {
    ph3_pilot_fail("missing_display_events",
                   paste0("complete-event data are required for ", prefix))
  }
  dat <- sample$data
  corrected <- analysis$ph3_legacy_pilot$event_corrections[[prefix]]
  if (!is.data.frame(corrected) || !identical(corrected$event_row, seq_len(nrow(dat))) ||
      !all(c("corrected_signal", "correction_status") %in% names(corrected)) ||
      !all(corrected$correction_status == "available") ||
      any(!is.finite(corrected$corrected_signal))) {
    ph3_pilot_fail("invalid_display_correction",
                   paste0("available event-level background correction is required for ", prefix))
  }
  # The three values below define the unchanged analytic regions: lower
  # analysis bound, 4N-region lower bound, and upper analysis bound.
  # The visual window deliberately extends beyond those regions to 1.6N--4.4N
  # and never contributes to eligibility, cutoff estimation, or summaries.
  bounds <- c(as.numeric(unlist(analysis$config$g1_x_range))[[1L]],
              as.numeric(unlist(analysis$config$g2m_x_range)))
  x_limits <- c(0.8, 2.2) * as.numeric(analysis$config$dna_2n_value)
  four_n_lower <- bounds[[2L]]
  display <- data.frame(dna_norm = dat$dna_norm,
                        corrected_signal = corrected$corrected_signal,
                        displayed_signal = corrected$corrected_signal + metadata_row$display_offset[[1L]])[
    is.finite(dat$dna_norm) & is.finite(corrected$corrected_signal) &
      dat$dna_norm >= x_limits[[1L]] & dat$dna_norm <= x_limits[[2L]],
    c("dna_norm", "corrected_signal", "displayed_signal"), drop = FALSE
  ]
  if (nrow(display) < 10L) {
    ph3_pilot_fail("insufficient_display_events",
                   paste0("at least 10 finite background-subtracted pH3 events are required for ", prefix))
  }
  positive_domain <- display$displayed_signal > 0
  expected_total <- metadata_row$visual_window_finite_event_count[[1L]]
  expected_positive <- metadata_row$positive_domain_display_event_count[[1L]]
  expected_nonpositive <- metadata_row$nonpositive_display_event_count[[1L]]
  if (!identical(as.integer(nrow(display)), expected_total) ||
      !identical(as.integer(sum(positive_domain)), expected_positive) ||
      !identical(as.integer(sum(!positive_domain)), expected_nonpositive)) {
    ph3_pilot_fail("log_display_qc_mismatch",
                   paste0("log-display QC must exactly reconcile visual events for ", prefix))
  }
  display <- display[positive_domain, , drop = FALSE]
  if (!nrow(display)) {
    return(ph3_pilot_unavailable_log_display_panel(analysis, metadata_row))
  }
  # Estimate pseudocolor density in the same coordinate system used by the
  # log10 panel.  Density in linear fluorescence units would be warped by the
  # later y transformation and can create a display-only horizontal artifact.
  # The plotted events remain on their restored linear signal scale; this does
  # not alter corrections, cutoffs, positivity, or summaries.
  display$density <- compute_point_density(
    display$dna_norm, log10(display$displayed_signal), n = 300L,
    bandwidth_multiplier = 0.5
  )
  display$density_color <- prepare_density_color(display$density)
  display <- display[order(display$density_color), , drop = FALSE]
  y_limits <- range(display$displayed_signal, finite = TRUE)
  if (length(y_limits) != 2L || any(!is.finite(y_limits)) ||
      any(y_limits <= 0) || !isTRUE(y_limits[[1L]] < y_limits[[2L]])) {
    ph3_pilot_fail("invalid_log_display_y_limits",
                   paste0("finite, positive, ordered log10 display limits are required for ", prefix))
  }
  regions <- data.frame(
    xmin = c(bounds[[1L]], four_n_lower), xmax = c(four_n_lower, bounds[[3L]]),
    ymin = rep(y_limits[[1L]], 2L), ymax = rep(y_limits[[2L]], 2L),
    label = c("Below 4N", "4N"), stringsAsFactors = FALSE
  )
  subtitle <- paste0(
    "Sample: ", prefix, ". Background-subtracted pH3 with display offset; log10 display. ",
    "Corrected cutoff from matched Untreated: ", metadata_row$cutoff_control_prefix[[1L]], "."
  )
  if (identical(metadata_row$cutoff_status[[1L]], "available")) {
    subtitle <- paste0(subtitle, " Corrected 4N cutoff = ",
                       format(metadata_row$cutoff_corrected_signal[[1L]], trim = TRUE),
                       "; displayed at ",
                       format(metadata_row$displayed_cutoff_signal[[1L]], trim = TRUE), ".")
  } else {
    subtitle <- paste0(subtitle, " Cutoff unavailable: ",
                       metadata_row$cutoff_reason_code[[1L]], ".")
  }
  if (identical(metadata_row$log_display_status[[1L]], "available_nonpositive_display_values_excluded")) {
    subtitle <- paste0(subtitle, " ", metadata_row$nonpositive_display_event_count[[1L]],
                       " nonpositive display value(s) excluded from this log10 plot only.")
  }
  plot <- ggplot2::ggplot(display, ggplot2::aes(x = dna_norm, y = displayed_signal,
                                                  colour = density_color)) +
    ggplot2::geom_rect(
      data = regions,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                   fill = label),
      inherit.aes = FALSE, alpha = 0.07, colour = NA
    ) +
    ggplot2::geom_point(shape = 16, size = 0.35, stroke = 0, alpha = 0.9) +
    ggplot2::scale_colour_gradientn(
      colours = refined_density_palette(), limits = c(0, 1),
      oob = scales::squish, name = "Relative density"
    ) +
    ggplot2::guides(colour = ggplot2::guide_colourbar(
      barheight = grid::unit(0.32, "in"), barwidth = grid::unit(0.07, "in"),
      title.theme = ggplot2::element_text(size = 6),
      label.theme = ggplot2::element_text(size = 5)
    )) +
    ggplot2::scale_fill_manual(values = c("Below 4N" = "#4C78A8", "4N" = "#F58518"),
                               guide = "none") +
    ggplot2::geom_vline(xintercept = bounds, linetype = c("dashed", "solid", "dashed"),
                        colour = c("grey35", "#F58518", "grey35"), linewidth = 0.45) +
    ggplot2::scale_x_continuous(breaks = c(analysis$config$dna_2n_value,
                                            2 * analysis$config$dna_2n_value),
                                labels = c("2N", "4N"), minor_breaks = NULL) +
    ggplot2::scale_y_log10(labels = scales::label_comma()) +
    ggplot2::coord_cartesian(xlim = x_limits, ylim = y_limits, clip = "off") +
    ggplot2::labs(
      title = paste0(metadata_row$replicate[[1L]], " \u2014 ", metadata_row$condition[[1L]]),
      subtitle = subtitle, x = "Normalized DNA content",
      y = "Background-subtracted pH3 fluorescence (display offset restored; log10 display)",
      caption = paste0(metadata_row$provenance_label[[1L]],
                       ". Display-only analytic regions; no verified FlowJo geometry is claimed.")
    ) +
    ggplot2::theme_classic(base_size = 10) +
    # Keep the data viewport square. The report supplies physical figure
    # dimensions intended to leave an approximately 1.25 x 1.25 inch data
    # viewport; this is display-only and does not alter density, cutoffs, or
    # regions.
    ggplot2::theme(
      aspect.ratio = 1,
      plot.caption = ggplot2::element_text(hjust = 0, size = 6),
      plot.title = ggplot2::element_text(size = 8),
      plot.subtitle = ggplot2::element_text(size = 6),
      axis.title = ggplot2::element_text(size = 7),
      axis.text = ggplot2::element_text(size = 6),
      legend.title = ggplot2::element_text(size = 6),
      legend.text = ggplot2::element_text(size = 5),
      legend.key.height = grid::unit(0.09, "in"),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.box.margin = ggplot2::margin(0, 0, 0, 2)
    )
  if (identical(metadata_row$cutoff_status[[1L]], "available") &&
      metadata_row$displayed_cutoff_signal[[1L]] > 0) {
    plot <- plot + ggplot2::geom_hline(
      yintercept = metadata_row$displayed_cutoff_signal[[1L]], colour = "#C00000",
      linewidth = 0.7
    )
  } else if (identical(metadata_row$cutoff_status[[1L]], "available")) {
    plot <- plot + ggplot2::labs(subtitle = paste0(
      subtitle, " Displayed cutoff is nonpositive; no cutoff overlay is shown on the log10 panel."
    ))
  }
  plot
}

ph3_pilot_expected_cutoff_failure_codes <- function() {
  c(
    "insufficient_control_4n_events", "nonfinite_control_raw_signal",
    "zero_or_invalid_control_signal_range", "invalid_density",
    "no_interior_local_maximum", "ambiguous_dominant_peak",
    "no_right_side_local_minimum"
  )
}

ph3_pilot_cutoff_result <- function(signal, stage) {
  density <- tryCatch(ph3_raw_4n_density_cutoff(signal), error = identity)
  if (!inherits(density, "error")) return(list(available = TRUE, density = density,
                                                reason = NA_character_, detail = NA_character_))
  reason <- sub("^PH3 raw-4N cutoff failed \\[([^]]+)\\].*$", "\\1", conditionMessage(density))
  if (!reason %in% ph3_pilot_expected_cutoff_failure_codes()) stop(density)
  list(available = FALSE, density = NULL, reason = paste0(stage, "_", reason),
       detail = conditionMessage(density))
}

ph3_pilot_correct_sample <- function(data, regions) {
  required <- c("dna_norm", "target_raw")
  if (!is.data.frame(data) || !all(required %in% names(data)) || !nrow(data) ||
      !is.numeric(data$dna_norm) || !is.numeric(data$target_raw) ||
      any(!is.finite(data$dna_norm)) || any(!is.finite(data$target_raw))) {
    return(list(status = "unavailable", reason = "nonfinite_or_missing_sample_data",
                detail = "every complete-event DNA and raw pH3 value must be finite", events = NULL,
                provisional = NULL, fit = NULL))
  }
  provisional <- ph3_pilot_cutoff_result(data$target_raw[regions$four_n], "provisional_raw")
  if (!provisional$available) {
    return(list(status = "unavailable", reason = provisional$reason, detail = provisional$detail,
                events = NULL, provisional = provisional, fit = NULL))
  }
  provisional_positive <- data$target_raw > provisional$density$cutoff
  fit <- ph3_fit_background_model(data$dna_norm, data$target_raw,
                                  !provisional_positive, provisional_positive)
  if (!identical(fit$fit_status, "valid") || length(fit$predictions) != nrow(data) ||
      any(!is.finite(fit$predictions))) {
    return(list(status = "unavailable", reason = paste0("background_fit_", fit$validity_reason_code),
                detail = "the per-sample provisional-negative linear background fit was invalid",
                events = NULL, provisional = provisional, fit = fit))
  }
  list(status = "available", reason = NA_character_, detail = NA_character_,
       provisional = provisional, fit = fit,
       events = data.frame(event_row = seq_len(nrow(data)), dna_norm = data$dna_norm,
                           target_raw = data$target_raw,
                           provisional_negative_member = !provisional_positive,
                           provisional_positive_member = provisional_positive,
                           predicted_background = fit$predictions,
                           corrected_signal = data$target_raw - fit$predictions,
                           correction_status = "available", stringsAsFactors = FALSE))
}

ph3_build_legacy_csv_pilot <- function(analysis) {
  if (!inherits(analysis, "facs_analysis") ||
      !identical(analysis$config$ph3_input_profile, "legacy_csv_pilot_v1")) {
    ph3_pilot_fail("invalid_profile", "legacy CSV pilot profile is required")
  }
  manifest <- analysis$sample_manifest
  if (nrow(manifest) != length(analysis$normalized_data) || anyDuplicated(manifest$prefix)) {
    ph3_pilot_fail("invalid_manifest", "each configured pilot sample must resolve exactly once")
  }
  ph3_pilot_validate_normalized_data_alignment(analysis, manifest)
  control_label <- analysis$config$ph3_pilot$control_label
  if (any(vapply(split(manifest, manifest$replicate_index), function(x) {
    sum(x$condition == control_label) != 1L
  }, logical(1)))) {
    ph3_pilot_fail("ambiguous_control_mapping", "each clone requires exactly one configured Untreated control")
  }
  cutoff_rows <- list(); value_rows <- list(); density_rows <- list()
  correction_rows <- list(); provisional_density_rows <- list(); corrections <- list(); ci <- 0L; vi <- 0L
  for (set_index in sort(unique(manifest$replicate_index))) {
    rows <- which(manifest$replicate_index == set_index)
    control_row <- rows[manifest$condition[rows] == control_label]
    sample_corrections <- lapply(rows, function(row) {
      data <- analysis$normalized_data[[row]]$data
      ph3_pilot_correct_sample(data, ph3_pilot_regions(analysis$config, data))
    })
    names(sample_corrections) <- manifest$prefix[rows]
    for (i in seq_along(rows)) {
      row <- rows[[i]]; corrected <- sample_corrections[[i]]
      correction_rows[[length(correction_rows) + 1L]] <- data.frame(
        prefix = manifest$prefix[[row]], replicate = manifest$replicate[[row]],
        replicate_index = as.integer(set_index), condition = manifest$condition[[row]],
        provisional_cutoff_raw_signal = if (isTRUE(corrected$provisional$available)) corrected$provisional$density$cutoff else NA_real_,
        provisional_cutoff_status = if (isTRUE(corrected$provisional$available)) "available" else "unavailable",
        correction_status = corrected$status, correction_reason_code = corrected$reason,
        negative_event_count = if (is.null(corrected$fit)) NA_integer_ else corrected$fit$negative_event_count,
        positive_event_count = if (is.null(corrected$fit)) NA_integer_ else corrected$fit$positive_event_count,
        intercept = if (is.null(corrected$fit)) NA_real_ else corrected$fit$intercept,
        slope = if (is.null(corrected$fit)) NA_real_ else corrected$fit$slope,
        stringsAsFactors = FALSE)
      if (isTRUE(corrected$provisional$available)) {
        provisional_density_rows[[length(provisional_density_rows) + 1L]] <- data.frame(
          prefix = manifest$prefix[[row]], replicate_index = as.integer(set_index),
          raw_pH3 = corrected$provisional$density$density_x,
          density = corrected$provisional$density$density_y,
          provisional_cutoff_raw_signal = corrected$provisional$density$cutoff,
          stringsAsFactors = FALSE)
      }
      if (identical(corrected$status, "available")) {
        corrections[[manifest$prefix[[row]]]] <- corrected$events
      }
    }
    all_corrected <- all(vapply(sample_corrections, function(x) identical(x$status, "available"), logical(1)))
    control_corrected <- sample_corrections[[which(rows == control_row)]]
    final <- if (all_corrected) {
      control_data <- analysis$normalized_data[[control_row]]$data
      control_regions <- ph3_pilot_regions(analysis$config, control_data)
      ph3_pilot_cutoff_result(control_corrected$events$corrected_signal[control_regions$four_n], "final_corrected")
    } else list(available = FALSE, density = NULL, reason = "unavailable_sample_background_correction",
                detail = "every sample in a clone must have an available provisional cutoff and background fit")
    available <- final$available
    reason <- final$reason
    ci <- ci + 1L
    cutoff_rows[[ci]] <- data.frame(
      replicate = manifest$replicate[[control_row]], replicate_index = as.integer(set_index),
      control_prefix = manifest$prefix[[control_row]], control_condition = control_label,
      cutoff_corrected_signal = if (available) final$density$cutoff else NA_real_,
      cutoff_status = if (available) "available" else "unavailable_cutoff_failure",
      reason_code = reason, cutoff_basis = "background_subtracted_4n_pH3",
      provenance_label = ph3_pilot_label(analysis), stringsAsFactors = FALSE)
    density_rows[[ci]] <- if (available) data.frame(
      replicate_index = as.integer(set_index), corrected_pH3 = final$density$density_x, density = final$density$density_y,
      cutoff_corrected_signal = final$density$cutoff, stringsAsFactors = FALSE
    ) else data.frame(replicate_index = integer(), corrected_pH3 = numeric(), density = numeric(), cutoff_corrected_signal = numeric())
    for (row in rows) {
      data <- analysis$normalized_data[[row]]$data
      regions <- ph3_pilot_regions(analysis$config, data)
      corrected <- sample_corrections[[which(rows == row)]]
      positive <- if (available) corrected$events$corrected_signal > final$density$cutoff else rep(FALSE, nrow(data))
      values <- if (!available) {
        list(A = NA_real_, B = NA_real_, C = NA_real_, D = NA_real_)
      } else {
        list(
          A = ph3_pilot_value(positive & regions$four_n, regions$eligible),
          B = ph3_pilot_value(positive & regions$below_4n, regions$eligible),
          C = if (any(positive & regions$four_n)) stats::median(corrected$events$corrected_signal[positive & regions$four_n]) else NA_real_,
          D = if (any(positive & regions$below_4n)) stats::median(corrected$events$corrected_signal[positive & regions$below_4n]) else NA_real_
        )
      }
      for (outcome in names(values)) {
        vi <- vi + 1L
        value_rows[[vi]] <- data.frame(
          replicate = manifest$replicate[[row]], replicate_index = as.integer(set_index),
          condition = manifest$condition[[row]], condition_index = as.integer(manifest$condition_index[[row]]),
          prefix = manifest$prefix[[row]], outcome_id = outcome,
          value = values[[outcome]], value_status = if (available && is.finite(values[[outcome]])) "available" else "unavailable",
          reason_code = if (!available) "unavailable_cutoff_failure" else if (!is.finite(values[[outcome]])) "no_qualifying_events" else NA_character_,
          signal_basis = "background_subtracted",
          cutoff_corrected_signal = if (available) final$density$cutoff else NA_real_,
          provenance_label = ph3_pilot_label(analysis), stringsAsFactors = FALSE)
      }
    }
  }
  values <- do.call(rbind, value_rows); rownames(values) <- NULL
  summary <- do.call(rbind, lapply(split(values, interaction(values$condition_index, values$outcome_id, drop = TRUE)), function(x) {
    finite <- is.finite(x$value)
    data.frame(condition = x$condition[[1L]], condition_index = x$condition_index[[1L]], outcome_id = x$outcome_id[[1L]],
      mean_value = if (any(finite)) mean(x$value[finite]) else NA_real_,
      biological_replicate_count = nrow(x), finite_biological_replicate_count = sum(finite),
      summary_status = if (all(finite)) "available" else if (any(finite)) "available_partial_coverage" else "unavailable",
      provenance_label = ph3_pilot_label(analysis), stringsAsFactors = FALSE)
  })); rownames(summary) <- NULL
  cutoff_records <- do.call(rbind, cutoff_rows)
  display_offset_qc <- ph3_pilot_display_offset_qc(analysis, cutoff_records, corrections)
  log_display_qc <- ph3_pilot_log_display_qc(analysis, display_offset_qc, corrections)
  analysis$ph3_legacy_pilot <- list(
    schema_version = "ph3-legacy-csv-pilot-1.0.0", provenance_label = ph3_pilot_label(analysis),
    limitations = c("Legacy FlowJo CSV exports: event identity and gate containment are unverified.",
                    "Final cutoff, signal panels, and C/D values use per-sample background-subtracted pH3; provisional raw cutoffs are retained only for correction audit.",
                    "PILOT ONLY \u2014 NOT PUBLICATION-GRADE."),
    cutoff_records = cutoff_records, provisional_background_qc = do.call(rbind, correction_rows),
    provisional_density_curves = if (length(provisional_density_rows)) do.call(rbind, provisional_density_rows) else data.frame(
      prefix = character(), replicate_index = integer(), raw_pH3 = numeric(), density = numeric(), provisional_cutoff_raw_signal = numeric()),
    event_corrections = corrections, density_curves = do.call(rbind, density_rows),
    display_offset_qc = display_offset_qc,
    log_display_qc = log_display_qc,
    pseudocolor_metadata = ph3_pilot_pseudocolor_metadata(
      analysis, cutoff_records, display_offset_qc, log_display_qc
    ),
    biological_replicate_values = values, condition_summary = summary)
  analysis$warnings <- unique(c(analysis$warnings, paste0("PILOT / LIMITED-PROVENANCE: ", ph3_pilot_label(analysis))))
  analysis
}

#' Plot a limited-provenance legacy pH3 pilot report
#'
#' Produces four clearly labelled pilot panels from [analyze_facs_experiment()]
#' run with `ph3_input_profile: legacy_csv_pilot_v1`. Signal panels use
#' background-subtracted pH3, and this output is never publication-grade.
#'
#' @param analysis A completed legacy CSV pH3 pilot analysis.
#' @return A named list of four outcome panels, one background-subtracted pH3-versus-DNA
#'   pseudocolor panel per configured pilot sample, and pilot QC tables.
#' @export
plot_ph3_legacy_pilot_report <- function(analysis) {
  if (!inherits(analysis, "facs_analysis") || is.null(analysis$ph3_legacy_pilot)) {
    ph3_pilot_fail("missing_pilot_model", "a completed legacy CSV pilot analysis is required")
  }
  pilot <- analysis$ph3_legacy_pilot; values <- pilot$biological_replicate_values; summary <- pilot$condition_summary
  labels <- c(A = "4N pH3-positive cells (% of Analysis singlets, 2N-4N)", B = "Below-4N pH3-positive cells (% of Analysis singlets, 2N-4N)", C = "Median background-subtracted pH3 signal in 4N pH3-positive cells", D = "Median background-subtracted pH3 signal below 4N")
  panels <- lapply(names(labels), function(outcome) {
    points <- values[values$outcome_id == outcome & is.finite(values$value), , drop = FALSE]
    means <- summary[summary$outcome_id == outcome & is.finite(summary$mean_value), , drop = FALSE]
    ggplot2::ggplot(points, ggplot2::aes(condition, value, colour = condition)) +
      ggplot2::geom_point(position = ggplot2::position_jitter(width = .08, seed = 1L), size = 2.4) +
      ggplot2::geom_point(data = means, ggplot2::aes(condition, mean_value), inherit.aes = FALSE, shape = 95, size = 4) +
      ggplot2::labs(title = labels[[outcome]], x = NULL,
        y = if (outcome %in% c("A", "B")) "Percent" else "Background-subtracted pH3 fluorescence",
        caption = paste0(pilot$provenance_label, ". Final corrected-density cutoff from matched Untreated control.")) +
      ggplot2::theme_classic(base_size = 11)
  })
  names(panels) <- c("ph3_4n_positive_prevalence", "ph3_below_4n_positive_prevalence", "ph3_4n_positive_signal", "ph3_below_4n_positive_signal")
  metadata <- pilot$pseudocolor_metadata
  expected_metadata <- c("prefix", "replicate", "replicate_index", "condition", "condition_index",
                         "cutoff_control_prefix", "cutoff_corrected_signal", "cutoff_status",
                         "cutoff_reason_code", "display_offset_status", "display_offset_reason_code",
                         "shared_raw_negative_median", "shared_target_status", "shared_target_reason_code",
                         "display_offset", "displayed_cutoff_signal", "log_display_status",
                         "log_display_reason_code", "visual_window_finite_event_count",
                         "positive_domain_display_event_count", "nonpositive_display_event_count",
                         "signal_basis", "provenance_label")
  if (!is.data.frame(metadata) || !identical(names(metadata), expected_metadata) ||
      nrow(metadata) != nrow(analysis$sample_manifest) || anyDuplicated(metadata$prefix) ||
      !identical(metadata$prefix, analysis$sample_manifest$prefix)) {
    ph3_pilot_fail("invalid_pseudocolor_metadata",
                   "one display metadata row is required for every configured pilot sample")
  }
  cutoff_rows <- pilot$cutoff_records[match(metadata$replicate_index,
                                             pilot$cutoff_records$replicate_index), , drop = FALSE]
  if (nrow(cutoff_rows) != nrow(metadata) ||
      !identical(metadata$cutoff_control_prefix, cutoff_rows$control_prefix) ||
      !identical(metadata$cutoff_corrected_signal, cutoff_rows$cutoff_corrected_signal) ||
      !identical(metadata$cutoff_status, cutoff_rows$cutoff_status) ||
      !identical(metadata$cutoff_reason_code, cutoff_rows$reason_code)) {
    ph3_pilot_fail("pseudocolor_cutoff_provenance_mismatch",
                   "display cutoff metadata must exactly match its clone cutoff record")
  }
  offset_qc <- pilot$display_offset_qc
  expected_offset_qc <- c(
    "prefix", "replicate", "replicate_index", "condition", "cutoff_control_prefix",
    "cutoff_corrected_signal", "eligible_final_negative_event_count",
    "raw_negative_median", "corrected_negative_median", "shared_raw_negative_median",
    "shared_target_status", "shared_target_reason_code", "display_offset",
    "displayed_cutoff_signal", "display_offset_status", "display_offset_reason_code"
  )
  if (!is.data.frame(offset_qc) || !identical(names(offset_qc), expected_offset_qc) ||
      nrow(offset_qc) != nrow(metadata) || anyDuplicated(offset_qc$prefix) ||
      !identical(offset_qc$prefix, metadata$prefix) ||
      !identical(offset_qc$display_offset_status, metadata$display_offset_status) ||
      !identical(offset_qc$display_offset_reason_code, metadata$display_offset_reason_code) ||
      !identical(offset_qc$shared_raw_negative_median, metadata$shared_raw_negative_median) ||
      !identical(offset_qc$shared_target_status, metadata$shared_target_status) ||
      !identical(offset_qc$shared_target_reason_code, metadata$shared_target_reason_code) ||
      !identical(offset_qc$display_offset, metadata$display_offset) ||
      !identical(offset_qc$displayed_cutoff_signal, metadata$displayed_cutoff_signal)) {
    ph3_pilot_fail("pseudocolor_display_offset_provenance_mismatch",
                   "display offset metadata must exactly match its per-sample QC record")
  }
  log_display_qc <- pilot$log_display_qc
  expected_log_display_qc <- c(
    "prefix", "replicate", "replicate_index", "condition",
    "visual_window_finite_event_count", "positive_domain_display_event_count",
    "nonpositive_display_event_count", "log_display_status", "log_display_reason_code"
  )
  if (!is.data.frame(log_display_qc) || !identical(names(log_display_qc), expected_log_display_qc) ||
      nrow(log_display_qc) != nrow(metadata) || anyDuplicated(log_display_qc$prefix) ||
      !identical(log_display_qc$prefix, metadata$prefix) ||
      !identical(log_display_qc$log_display_status, metadata$log_display_status) ||
      !identical(log_display_qc$log_display_reason_code, metadata$log_display_reason_code) ||
      !identical(log_display_qc$visual_window_finite_event_count, metadata$visual_window_finite_event_count) ||
      !identical(log_display_qc$positive_domain_display_event_count, metadata$positive_domain_display_event_count) ||
      !identical(log_display_qc$nonpositive_display_event_count, metadata$nonpositive_display_event_count)) {
    ph3_pilot_fail("pseudocolor_log_display_provenance_mismatch",
                   "log-display metadata must exactly match its per-sample QC record")
  }
  pseudocolor_panels <- stats::setNames(lapply(seq_len(nrow(metadata)), function(i) {
    ph3_pilot_pseudocolor_panel(analysis, metadata[i, , drop = FALSE])
  }), metadata$prefix)
  list(panels = panels, pseudocolor_panels = pseudocolor_panels,
       pseudocolor_metadata = metadata, cutoff_records = pilot$cutoff_records,
       provisional_background_qc = pilot$provisional_background_qc,
       display_offset_qc = pilot$display_offset_qc,
       log_display_qc = pilot$log_display_qc,
       provisional_density_curves = pilot$provisional_density_curves, density_curves = pilot$density_curves,
       biological_replicate_values = values, condition_summary = summary, limitations = pilot$limitations)
}
