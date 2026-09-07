# EdU Slice 1 every-sample pseudocolor output contract ----------------------

# This module deliberately constructs a display-only view of an existing EdU
# analysis.  It neither refits a model nor alters event classifications,
# numerical results, or the analysis object supplied by the caller.

edu_negative_event_fingerprint <- function(dna_norm, target_raw) {
  if (length(dna_norm) != length(target_raw) || !length(dna_norm) ||
      any(!is.finite(dna_norm)) || any(!is.finite(target_raw))) {
    return(NA_character_)
  }
  bytes <- serialize(list(as.numeric(dna_norm), as.numeric(target_raw)), NULL,
                     version = 2L)
  paste(sprintf("%02x", as.integer(openssl::sha256(bytes))), collapse = "")
}

edu_display_suppression <- function(prefix, reason, detail) {
  list(prefix = prefix, status = "suppressed", reason_code = reason,
       reason_detail = detail)
}

edu_resolve_sample_display_offset <- function(sample, model, prefix) {
  required_model <- c("sample_prefix", "negative_event_index",
                      "negative_event_fingerprint")
  if (!is.list(model) || !all(required_model %in% names(model))) {
    return(edu_display_suppression(
      prefix, "missing_negative_membership_proof",
      "The established EdU background fit has no retained negative-event membership proof"
    ))
  }
  if (!identical(as.character(model$sample_prefix), as.character(prefix))) {
    return(edu_display_suppression(
      prefix, "sample_identity_mismatch",
      "The established EdU background-fit identity does not match this sample"
    ))
  }
  data <- sample$data
  required_data <- c("dna_norm", "target_raw", "target_bgsub")
  if (!is.data.frame(data) || !all(required_data %in% names(data))) {
    return(edu_display_suppression(
      prefix, "missing_display_signal",
      "The established analysis does not retain the required raw and corrected event signals"
    ))
  }
  index <- model$negative_event_index
  if (!is.numeric(index) || !length(index) || any(!is.finite(index)) ||
      any(index %% 1 != 0) || any(index < 1L) || any(index > nrow(data)) ||
      anyDuplicated(index)) {
    return(edu_display_suppression(
      prefix, "invalid_negative_membership_proof",
      "The retained EdU-negative event membership proof is invalid"
    ))
  }
  index <- as.integer(index)
  fingerprint <- edu_negative_event_fingerprint(
    data$dna_norm[index], data$target_raw[index]
  )
  if (is.na(fingerprint) || !identical(fingerprint, model$negative_event_fingerprint)) {
    return(edu_display_suppression(
      prefix, "negative_membership_proof_mismatch",
      "The retained EdU-negative event membership proof does not match this sample's event table"
    ))
  }
  raw <- data$target_raw[index]
  corrected <- data$target_bgsub[index]
  if (any(!is.finite(raw)) || any(!is.finite(corrected))) {
    return(edu_display_suppression(
      prefix, "invalid_negative_display_values",
      "The established EdU-negative events do not have finite raw and corrected signals"
    ))
  }
  raw_median <- stats::median(raw)
  corrected_median <- stats::median(corrected)
  offset <- raw_median - corrected_median
  if (!is.finite(raw_median) || !is.finite(corrected_median) || !is.finite(offset)) {
    return(edu_display_suppression(
      prefix, "unavailable_display_offset",
      "The required per-sample EdU display offset is not finite"
    ))
  }
  list(
    prefix = prefix, status = "available", reason_code = NA_character_,
    reason_detail = NA_character_, negative_event_n = length(index),
    raw_negative_median = raw_median,
    corrected_negative_median = corrected_median, display_offset = offset,
    membership_fingerprint = fingerprint
  )
}

edu_display_limits <- function(values, config) {
  values <- values[is.finite(values)]
  if (isTRUE(config$y_log10)) values <- values[values > 0]
  if (length(values) < 2L) return(NULL)
  limits <- stats::quantile(
    values,
    c(config$y_limit_lower_quantile, config$y_limit_upper_quantile),
    names = FALSE
  )
  if (length(limits) != 2L || any(!is.finite(limits)) || limits[[1L]] >= limits[[2L]] ||
      (isTRUE(config$y_log10) && limits[[1L]] <= 0)) return(NULL)
  as.numeric(limits)
}

edu_display_panel <- function(sample, manifest_row, offset_record, analysis) {
  prefix <- manifest_row$prefix[[1L]]
  identity <- paste0(
    "Condition: ", manifest_row$condition[[1L]],
    " | Biological replicate: ", manifest_row$replicate[[1L]],
    " | Technical acquisition: ", manifest_row$technical_replicate[[1L]]
  )
  if (!identical(offset_record$status, "available")) {
    return(list(
      plot = ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.55, hjust = 0.5,
                          label = paste("SUPPRESSED:", offset_record$reason_code,
                                        "\n", offset_record$reason_detail)) +
        ggplot2::labs(title = prefix, subtitle = identity,
                      x = "Normalized DNA content", y = "EdU display signal") +
        ggplot2::theme_void() + ggplot2::theme(plot.title = ggplot2::element_text(face = "bold")),
      status = offset_record$status, reason_code = offset_record$reason_code,
      displayed_event_n = 0L, display_limits = c(NA_real_, NA_real_)
    ))
  }
  data <- sample$data
  displayed <- data$target_bgsub + offset_record$display_offset
  limits <- edu_display_limits(displayed, analysis$config)
  if (is.null(limits)) {
    suppressed <- edu_display_suppression(
      prefix, "unavailable_display_domain",
      "The offset display coordinate has too few finite values in the required plotting domain"
    )
    return(edu_display_panel(sample, manifest_row, suppressed, analysis))
  }
  keep <- is.finite(data$dna_norm) & is.finite(displayed) &
    data$dna_norm >= analysis$config$x_limits[[1L]] &
    data$dna_norm <= analysis$config$x_limits[[2L]] &
    displayed >= limits[[1L]] & displayed <= limits[[2L]]
  if (isTRUE(analysis$config$y_log10)) keep <- keep & displayed > 0
  display <- data.frame(dna_norm = data$dna_norm[keep], displayed_signal = displayed[keep])
  if (nrow(display) < 10L) {
    suppressed <- edu_display_suppression(
      prefix, "too_few_display_events",
      "Fewer than ten established events are available in this sample's display coordinate"
    )
    return(edu_display_panel(sample, manifest_row, suppressed, analysis))
  }
  density_y <- if (isTRUE(analysis$config$y_log10)) log10(display$displayed_signal) else display$displayed_signal
  display$density <- compute_point_density(
    display$dna_norm, density_y,
    bandwidth_multiplier = analysis$config$density_bandwidth
  )
  display$density_color <- prepare_density_color(
    display$density, analysis$config$density_lower_clip,
    analysis$config$density_upper_clip, analysis$config$density_gamma
  )
  display <- display[order(display$density_color), , drop = FALSE]
  plot <- ggplot2::ggplot(display, ggplot2::aes(dna_norm, displayed_signal,
                                                colour = density_color)) +
    ggplot2::geom_point(size = analysis$config$point_size, stroke = 0) +
    ggplot2::scale_color_gradientn(colours = resolve_palette(analysis$config$palette),
                                   limits = c(0, 1), oob = scales::squish,
                                   name = "Relative density") +
    ggplot2::scale_x_continuous(breaks = c(analysis$config$dna_2n_value,
                                            2 * analysis$config$dna_2n_value),
                                labels = c("2N", "4N")) +
    ggplot2::coord_cartesian(xlim = as.numeric(analysis$config$x_limits), ylim = limits) +
    ggplot2::labs(title = prefix, subtitle = identity,
                  x = "Normalized DNA content",
                  y = "EdU background-subtracted fluorescence (display offset)") +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(aspect.ratio = 1,
                   plot.title = ggplot2::element_text(face = "bold"))
  if (isTRUE(analysis$config$y_log10)) plot <- plot + ggplot2::scale_y_log10()
  list(plot = plot, status = "available", reason_code = NA_character_,
       displayed_event_n = nrow(display), display_limits = limits)
}

#' Build every-sample EdU pseudocolor panels using per-sample display offsets
#'
#' This display-only report model reconstructs a per-sample EdU display
#' coordinate from the exact negative events used by the already-completed
#' background fit. It never changes the analysis object or numerical results.
#'
#' @param analysis A completed EdU `facs_analysis` object.
#' @return A list of individually identified ggplot panels and per-sample
#'   display-offset/QC records.
#' @export
build_edu_pseudocolor_output_contract <- function(analysis) {
  validate_analysis_object(analysis)
  if (!identical(analysis$config$plot_type, "edu")) {
    stop("EdU pseudocolor output requires an EdU analysis.", call. = FALSE)
  }
  manifest <- analysis$sample_manifest
  expected <- as.character(manifest$prefix)
  if (is.null(names(analysis$normalized_data)) ||
      !identical(names(analysis$normalized_data), expected) ||
      is.null(names(analysis$models)) || !identical(names(analysis$models), expected)) {
    stop("EdU pseudocolor output requires normalized data and background models named in exact manifest-prefix order.", call. = FALSE)
  }
  offsets <- lapply(seq_len(nrow(manifest)), function(i) {
    edu_resolve_sample_display_offset(
      analysis$normalized_data[[i]], analysis$models[[i]], manifest$prefix[[i]]
    )
  })
  panels <- lapply(seq_len(nrow(manifest)), function(i) {
    edu_display_panel(analysis$normalized_data[[i]], manifest[i, , drop = FALSE],
                      offsets[[i]], analysis)
  })
  offset_qc <- do.call(rbind, lapply(offsets, function(x) {
    data.frame(prefix = x$prefix, display_status = x$status,
               suppression_reason_code = x$reason_code,
               suppression_reason_detail = x$reason_detail,
               negative_event_n = x$negative_event_n %||% NA_integer_,
               raw_negative_median = x$raw_negative_median %||% NA_real_,
               corrected_negative_median = x$corrected_negative_median %||% NA_real_,
               display_offset = x$display_offset %||% NA_real_,
               membership_fingerprint = x$membership_fingerprint %||% NA_character_,
               stringsAsFactors = FALSE)
  }))
  panel_qc <- data.frame(
    manifest[, c("prefix", "condition", "replicate", "technical_replicate"), drop = FALSE],
    panel_status = vapply(panels, `[[`, character(1), "status"),
    panel_reason_code = vapply(panels, function(x) x$reason_code %||% NA_character_, character(1)),
    displayed_event_n = vapply(panels, `[[`, integer(1), "displayed_event_n"),
    y_limit_lower = vapply(panels, function(x) x$display_limits[[1L]], numeric(1)),
    y_limit_upper = vapply(panels, function(x) x$display_limits[[2L]], numeric(1)),
    stringsAsFactors = FALSE
  )
  structure(list(
    schema_version = "edu-pseudocolor-output-contract-1.0.0",
    panels = stats::setNames(lapply(panels, `[[`, "plot"), expected),
    display_offset_qc = offset_qc, panel_qc = panel_qc,
    provenance = list(
      display_offset_method = "raw_negative_median_minus_corrected_negative_median_v1",
      analytical_values_mutated = FALSE,
      panel_identity_key = "manifest_prefix"
    )
  ), class = "edu_pseudocolor_output_contract")
}
