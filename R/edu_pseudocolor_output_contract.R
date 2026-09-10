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
    # This report deliberately uses the project-standard refined palette,
    # matching the preferred pH3 report density rendering.  It is a
    # presentation choice only; density values and analytical results are not
    # changed by this palette selection.
    ggplot2::scale_color_gradientn(colours = refined_density_palette(),
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

edu_positivity_report_fail <- function(reason, detail) {
  stop("EdU positivity report validation failed [", reason, "]: ", detail,
       call. = FALSE)
}

edu_positivity_required_columns <- function() {
  c(
    "replicate", "replicate_index", "technical_replicate", "condition",
    "condition_index", "region", "region_label", "region_index",
    "numerator_n", "denominator_n", "regional_edu_positive_pct",
    "source_population", "dna_interval_min", "dna_interval_max",
    "dna_interval_lower_inclusive", "dna_interval_upper_inclusive",
    "metric_status", "aggregation_level", "aggregation_method"
  )
}

edu_validate_acquisition_manifest_identity <- function(
    acquisition, manifest, category_col = NULL, category_values = NULL,
    context
) {
  identity_cols <- c(
    "replicate", "replicate_index", "technical_replicate", "condition",
    "condition_index"
  )
  required <- c(identity_cols, category_col)
  if (!all(required %in% names(acquisition))) {
    edu_positivity_report_fail(
      "acquisition_manifest_identity_mismatch",
      paste0(context, " is missing the manifest identity field(s) required for report provenance")
    )
  }
  expected <- manifest[, identity_cols, drop = FALSE]
  if (!is.null(category_col)) {
    if (!is.character(category_values) || !length(category_values)) {
      edu_positivity_report_fail(
        "invalid_report_categories",
        paste0(context, " requires explicit approved category values")
      )
    }
    expected <- do.call(rbind, lapply(category_values, function(category) {
      out <- expected
      out[[category_col]] <- category
      out
    }))
  }
  key <- function(df, columns) {
    do.call(paste, c(lapply(df[columns], as.character), sep = "\r"))
  }
  actual_key <- key(acquisition, required)
  expected_key <- key(expected, required)
  if (anyDuplicated(actual_key) || length(actual_key) != length(expected_key) ||
      !identical(sort(actual_key), sort(expected_key))) {
    edu_positivity_report_fail(
      "acquisition_manifest_identity_mismatch",
      paste0(context,
             " must contain exactly one current-manifest row for every required acquisition identity")
    )
  }
  invisible(TRUE)
}

edu_validate_regional_positivity_for_report <- function(analysis) {
  q <- analysis$quantitation
  acquisition_name <- "edu_regional_positivity_acquisition"
  aggregate_name <- "edu_regional_positivity"
  if (!is.list(q) || !is.data.frame(q[[acquisition_name]]) ||
      !is.data.frame(q[[aggregate_name]])) {
    edu_positivity_report_fail(
      "missing_canonical_regional_positivity",
      "the completed analysis must retain both canonical regional positivity tables"
    )
  }
  acquisition <- q[[acquisition_name]]
  aggregate <- q[[aggregate_name]]
  required <- edu_positivity_required_columns()
  missing <- setdiff(required, names(acquisition))
  if (length(missing)) {
    edu_positivity_report_fail(
      "invalid_canonical_regional_positivity",
      "the canonical acquisition table is incomplete or has duplicate identity rows"
    )
  }
  expected_regions <- c("early", "mid", "late")
  manifest <- analysis$sample_manifest
  edu_validate_acquisition_manifest_identity(
    acquisition, manifest, category_col = "region",
    category_values = expected_regions,
    context = "canonical regional positivity acquisition table"
  )
  if (nrow(acquisition) != 3L * nrow(manifest) ||
      !identical(as.character(acquisition$region),
                 rep(expected_regions, times = nrow(manifest))) ||
      !identical(as.character(acquisition$source_population), rep(
        "eligible_single_cells_in_dna_region", nrow(acquisition)
      ))) {
    edu_positivity_report_fail(
      "invalid_canonical_regional_positivity",
      "the canonical table must contain Early, Mid, and Late S rows in manifest order"
    )
  }
  expected_aggregate <- average_edu_metric_table(
    acquisition, "regional_edu_positive_pct",
    c("region", "region_label", "region_index")
  )
  if (!identical(aggregate, expected_aggregate)) {
    edu_positivity_report_fail(
      "regional_positivity_reconciliation_failed",
      "the biological-replicate regional table must exactly reconcile to canonical acquisition rows"
    )
  }
  aggregate
}

edu_collect_2to4n_positivity_for_report <- function(analysis) {
  manifest <- analysis$sample_manifest
  samples <- analysis$normalized_data
  if (is.null(names(samples)) || !identical(names(samples), as.character(manifest$prefix))) {
    edu_positivity_report_fail(
      "sample_identity_mismatch",
      "normalized samples must be named in exact manifest-prefix order"
    )
  }
  lower <- analysis$config$dna_2n_value
  upper <- 2 * analysis$config$dna_2n_value
  rows <- lapply(seq_len(nrow(manifest)), function(i) {
    sample <- samples[[i]]
    source <- sample$data
    classified <- sample$edu_event_classification
    required <- c("event_row", "dna_norm", "computed_positive", "positivity_eligible")
    if (!is.data.frame(source) || !is.data.frame(classified) ||
        !all(c("dna_norm", "edu_computed_positive") %in% names(source)) ||
        !all(required %in% names(classified)) || nrow(source) != nrow(classified) ||
        !identical(classified$event_row, seq_len(nrow(source))) ||
        !identical(classified$dna_norm, source$dna_norm) ||
        !identical(classified$computed_positive, source$edu_computed_positive)) {
      edu_positivity_report_fail(
        "classification_identity_mismatch",
        paste0("the retained EdU classification must exactly match source events for ",
               manifest$prefix[[i]])
      )
    }
    positive_known <- !is.na(source$edu_computed_positive) &
      source$edu_computed_positive %in% c(TRUE, FALSE)
    expected_eligible <- is.finite(source$dna_norm) & positive_known
    if (!identical(classified$positivity_eligible, expected_eligible)) {
      edu_positivity_report_fail(
        "classification_eligibility_mismatch",
        paste0("the retained EdU positivity eligibility must reconcile for ",
               manifest$prefix[[i]])
      )
    }
    in_2to4n <- classified$positivity_eligible &
      classified$dna_norm >= lower & classified$dna_norm <= upper
    denominator_n <- sum(in_2to4n)
    numerator_n <- sum(in_2to4n & classified$computed_positive %in% TRUE)
    metadata <- edu_metric_metadata(
      "eligible_single_cells_in_2n_to_4n_inclusive",
      classified$display_offset[[1L]],
      classified$display_offset_applied[[1L]] %in% TRUE,
      lower, upper, display_transform = classified$display_transform[[1L]]
    )
    metadata$dna_interval_upper_inclusive <- TRUE
    cbind(
      edu_sample_metadata(manifest, i),
      data.frame(
        region = "2to4n", region_label = "2N\u20134N", region_index = 1L,
        numerator_n = numerator_n, denominator_n = denominator_n,
        two_to_four_n_edu_positive_pct = if (denominator_n > 0L) {
          100 * numerator_n / denominator_n
        } else NA_real_,
        stringsAsFactors = FALSE
      ),
      metadata,
      data.frame(
        metric_status = if (denominator_n > 0L) "ok" else "zero_denominator",
        stringsAsFactors = FALSE
      )
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

edu_positivity_plot <- function(
    points, value_col, category_col, category_levels, category_label,
    y_label, analysis
) {
  if (!is.data.frame(points) || !all(c(
    value_col, category_col, "condition", "condition_index", "replicate",
    "replicate_index"
  ) %in% names(points))) {
    edu_positivity_report_fail("invalid_plot_input", "validated positivity rows are required")
  }
  categories <- as.character(points[[category_col]])
  if (!identical(sort(unique(categories)), sort(category_levels))) {
    edu_positivity_report_fail("invalid_plot_input", "the approved positivity categories are missing or duplicated")
  }
  summary <- summarize_across_replicates(
    points, value_col, c(category_col, "condition", "condition_index")
  )
  names(summary)[[1L]] <- category_col
  summary[[category_col]] <- factor(summary[[category_col]], levels = category_levels)
  points[[category_col]] <- factor(points[[category_col]], levels = category_levels)
  condition_levels <- unique(summary$condition[order(summary$condition_index)])
  summary$condition <- factor(summary$condition, levels = condition_levels)
  points$condition <- factor(points$condition, levels = condition_levels)
  summary$error <- switch(
    analysis$config$quant_error_bar,
    sd = summary$sd, sem = summary$sem, none = rep(NA_real_, nrow(summary))
  )
  colours <- resolve_condition_colors(condition_levels, analysis$config$bar_colors)
  dodge <- 0.82
  plot <- ggplot2::ggplot(
    summary, ggplot2::aes(
      x = .data[[category_col]], y = .data$mean, fill = .data$condition
    )
  ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(width = dodge), width = 0.76,
      color = "white", linewidth = 0.25, na.rm = TRUE
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = pmax(0, mean - error), ymax = mean + error),
      position = ggplot2::position_dodge(width = dodge), width = 0.14,
      linewidth = 0.45, na.rm = TRUE
    )
  if (isTRUE(analysis$config$quant_show_points)) {
    plot <- plot + ggplot2::geom_point(
      data = points,
      ggplot2::aes(
        x = .data[[category_col]], y = .data[[value_col]], group = .data$condition
      ),
      inherit.aes = FALSE, shape = 21, fill = "white", color = "black",
      size = 1.8, stroke = 0.4,
      position = ggplot2::position_jitterdodge(
        jitter.width = 0.05, dodge.width = dodge, seed = 1L
      ), na.rm = TRUE
    )
  }
  plot +
    ggplot2::scale_fill_manual(values = colours, name = "Condition") +
    ggplot2::scale_y_continuous(
      labels = function(x) paste0(x, "%"),
      expand = ggplot2::expansion(mult = c(0, 0.10))
    ) +
    ggplot2::labs(x = category_label, y = y_label,
                  caption = "Points are biological-replicate values; bars show condition means.") +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))
}

edu_build_positivity_report <- function(analysis) {
  regional <- edu_validate_regional_positivity_for_report(analysis)
  overall_acquisition <- edu_collect_2to4n_positivity_for_report(analysis)
  overall <- average_edu_metric_table(
    overall_acquisition, "two_to_four_n_edu_positive_pct",
    c("region", "region_label", "region_index")
  )
  list(
    schema_version = "edu-positivity-report-1.0.0",
    overall_2to4n_acquisition = overall_acquisition,
    overall_2to4n_biological_replicate = overall,
    regional_biological_replicate = regional,
    panels = list(
      two_to_four_n = edu_positivity_plot(
        overall, "two_to_four_n_edu_positive_pct", "region_label", "2N\u20134N",
        "DNA content", "2N\u20134N EdU-positive percentage", analysis
      ),
      early_mid_late_s = edu_positivity_plot(
        regional, "regional_edu_positive_pct", "region_label",
        c("Early S", "Mid S", "Late S"), "S-phase region",
        "EdU-positive percentage within DNA region", analysis
      )
    ),
    provenance = list(
      overall_definition = paste0(
        "computed EdU-positive events with DNA in inclusive [",
        analysis$config$dna_2n_value, ", ", 2 * analysis$config$dna_2n_value,
        "] divided by all positivity-eligible events in the same interval"
      ),
      regional_source = "canonical_edu_regional_positivity_biological_replicate",
      aggregation = "unweighted technical-acquisition mean within biological replicate"
    )
  )
}

edu_validate_intensity_table_for_report <- function(
    analysis, acquisition_name, aggregate_name, value_col, group_cols,
    expected_source_population, expected_rows_per_acquisition,
    category_col = NULL, category_values = NULL
) {
  q <- analysis$quantitation
  if (!is.list(q) || !is.data.frame(q[[acquisition_name]]) ||
      !is.data.frame(q[[aggregate_name]])) {
    edu_positivity_report_fail(
      "missing_canonical_intensity",
      paste0("the completed analysis must retain `", acquisition_name,
             "` and `", aggregate_name, "`")
    )
  }
  acquisition <- q[[acquisition_name]]
  aggregate <- q[[aggregate_name]]
  required <- c(
    "replicate", "replicate_index", "technical_replicate", "condition",
    "condition_index", "source_population_n", value_col,
    "source_population", "signal_transform", "aggregation_level",
    "aggregation_method", "metric_status"
  )
  if (!all(required %in% names(acquisition))) {
    edu_positivity_report_fail(
      "invalid_canonical_intensity",
      paste0("the canonical `", acquisition_name,
             "` is missing the required intensity fields")
    )
  }
  edu_validate_acquisition_manifest_identity(
    acquisition, analysis$sample_manifest, category_col, category_values,
    paste0("canonical `", acquisition_name, "`")
  )
  if (nrow(acquisition) != expected_rows_per_acquisition * nrow(analysis$sample_manifest) ||
      !identical(as.character(acquisition$source_population),
                 rep(expected_source_population, nrow(acquisition)))) {
    edu_positivity_report_fail(
      "invalid_canonical_intensity",
      paste0("the canonical `", acquisition_name,
             "` rows do not have the approved identity and source population")
    )
  }
  if (!all(acquisition$signal_transform == "background_subtracted")) {
    edu_positivity_report_fail(
      "invalid_canonical_intensity",
      paste0("the canonical `", acquisition_name,
             "` must contain only background-subtracted values")
    )
  }
  expected <- average_edu_metric_table(acquisition, value_col, group_cols)
  if (!identical(aggregate, expected)) {
    edu_positivity_report_fail(
      "intensity_reconciliation_failed",
      paste0("the canonical `", aggregate_name,
             "` must exactly reconcile to its acquisition rows")
    )
  }
  aggregate
}

edu_intensity_plot <- function(
    points, value_col, category_col, category_levels, category_label,
    y_label, analysis
) {
  if (!is.data.frame(points) || !all(c(
    value_col, category_col, "condition", "condition_index", "replicate",
    "replicate_index"
  ) %in% names(points))) {
    edu_positivity_report_fail("invalid_intensity_plot_input", "validated canonical intensity rows are required")
  }
  if (!identical(sort(unique(as.character(points[[category_col]]))),
                 sort(category_levels))) {
    edu_positivity_report_fail("invalid_intensity_plot_input", "the approved intensity categories are missing or duplicated")
  }
  summary <- summarize_across_replicates(
    points, value_col, c(category_col, "condition", "condition_index")
  )
  names(summary)[[1L]] <- category_col
  summary[[category_col]] <- factor(summary[[category_col]], levels = category_levels)
  points[[category_col]] <- factor(points[[category_col]], levels = category_levels)
  condition_levels <- unique(summary$condition[order(summary$condition_index)])
  summary$condition <- factor(summary$condition, levels = condition_levels)
  points$condition <- factor(points$condition, levels = condition_levels)
  summary$error <- switch(
    analysis$config$quant_error_bar,
    sd = summary$sd, sem = summary$sem, none = rep(NA_real_, nrow(summary))
  )
  colours <- resolve_condition_colors(condition_levels, analysis$config$bar_colors)
  dodge <- 0.82
  plot <- ggplot2::ggplot(
    points,
    ggplot2::aes(
      x = .data[[category_col]], y = .data[[value_col]], colour = .data$condition
    )
  )
  if (isTRUE(analysis$config$quant_show_points)) {
    plot <- plot + ggplot2::geom_point(
      shape = 16, size = 2.0,
      position = ggplot2::position_jitterdodge(
        jitter.width = 0.05, dodge.width = dodge, seed = 1L
      ), na.rm = TRUE
    )
  }
  plot +
    ggplot2::geom_errorbar(
      data = summary,
      ggplot2::aes(
        x = .data[[category_col]], ymin = .data$mean - .data$error,
        ymax = .data$mean + .data$error, colour = .data$condition
      ), position = ggplot2::position_dodge(width = dodge), width = 0.14,
      linewidth = 0.45, na.rm = TRUE, inherit.aes = FALSE
    ) +
    ggplot2::geom_point(
      data = summary,
      ggplot2::aes(
        x = .data[[category_col]], y = .data$mean, colour = .data$condition
      ), position = ggplot2::position_dodge(width = dodge), shape = 95,
      size = 5, stroke = 1.1, na.rm = TRUE, inherit.aes = FALSE
    ) +
    ggplot2::scale_colour_manual(values = colours, name = "Condition") +
    ggplot2::labs(
      x = category_label, y = y_label,
      caption = "Points are biological-replicate values; horizontal marks show condition means."
    ) +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))
}

# Derive a report-only fold value from the retained biological-replicate
# intensity table.  This intentionally happens after the established
# acquisition-to-biological-replicate averaging, so each configured matched
# reference sample is the within-biological-replicate denominator. It neither refits a
# background model nor changes any canonical quantitative table.
edu_reference_relative_intensity <- function(
    points, value_col, manifest, category_col = NULL
) {
  required_manifest <- c(
    "replicate", "replicate_index", "technical_replicate", "model_group",
    "condition", "reference_condition", "is_reference"
  )
  if (!is.data.frame(manifest) || !all(required_manifest %in% names(manifest))) {
    edu_positivity_report_fail(
      "missing_reference_provenance",
      "the current sample manifest lacks the reference identity required for reference-relative intensity"
    )
  }
  required_points <- c("replicate", "replicate_index", "condition", value_col)
  if (!is.data.frame(points) || !all(required_points %in% names(points))) {
    edu_positivity_report_fail(
      "invalid_reference_intensity_input",
      "validated biological-replicate intensity rows are required"
    )
  }
  if (!is.null(category_col) && !category_col %in% names(points)) {
    edu_positivity_report_fail(
      "invalid_reference_intensity_input",
      "the required report category is absent from the intensity rows"
    )
  }

  group_rows <- split(seq_len(nrow(manifest)), manifest$model_group)
  reference_by_group <- vapply(group_rows, function(idx) {
    group <- manifest[idx, , drop = FALSE]
    declared <- unique(as.character(group$reference_condition))
    declared <- declared[!is.na(declared) & nzchar(declared)]
    reference <- which(group$is_reference %in% TRUE)
    if (length(declared) != 1L || length(reference) != 1L ||
        !identical(as.character(group$condition[[reference]]), declared[[1L]])) {
      edu_positivity_report_fail(
        "invalid_reference_identity",
        paste0("each biological/technical replicate pair must contain exactly one configured reference; invalid group: ",
               group$model_group[[1L]])
      )
    }
    declared[[1L]]
  }, character(1))
  reference_by_replicate <- lapply(
    split(seq_len(nrow(manifest)), manifest$replicate_index),
    function(idx) {
      values <- unique(unname(reference_by_group[manifest$model_group[idx]]))
      if (length(values) != 1L) {
        edu_positivity_report_fail(
          "mixed_reference_identity",
          paste0("all technical acquisitions must use one configured reference within biological replicate ",
                 manifest$replicate[[idx[[1L]]]])
        )
      }
      values[[1L]]
    }
  )
  reference_by_replicate <- unlist(reference_by_replicate, use.names = TRUE)

  replicate_key <- as.character(points$replicate_index)
  expected_reference <- unname(reference_by_replicate[replicate_key])
  if (anyNA(expected_reference) || any(!nzchar(expected_reference))) {
    edu_positivity_report_fail(
      "reference_coverage_mismatch",
      "every report intensity row must map to a configured biological-replicate reference"
    )
  }
  category_key <- if (is.null(category_col)) rep("", nrow(points)) else
    as.character(points[[category_col]])
  group_key <- paste(replicate_key, category_key, sep = "\r")
  out <- points
  out$reference_condition <- expected_reference
  out$reference_intensity <- NA_real_
  out$reference_relative_intensity <- NA_real_
  for (key in unique(group_key)) {
    idx <- which(group_key == key)
    reference_idx <- idx[out$condition[idx] == out$reference_condition[idx]]
    if (length(reference_idx) != 1L) {
      edu_positivity_report_fail(
        "reference_coverage_mismatch",
        "each biological-replicate intensity/category group must contain exactly one configured reference row"
      )
    }
    reference_value <- out[[value_col]][[reference_idx]]
    values <- out[[value_col]][idx]
    out$reference_intensity[idx] <- reference_value
    if (!is.finite(reference_value) || reference_value <= 0) {
      # A background-subtracted signal can validly be nonpositive.  It cannot
      # define a fold denominator, however, so retain the canonical result but
      # expose no substitute ratio for any condition in this exact group.
      out$reference_normalization_status[idx] <-
        "unavailable_invalid_reference_intensity"
      next
    }
    finite <- is.finite(values)
    out$reference_relative_intensity[idx[finite]] <-
      values[finite] / reference_value
    out$reference_normalization_status[idx[finite]] <- "available"
    out$reference_normalization_status[idx[!finite]] <-
      "unavailable_nonfinite_canonical_intensity"
  }
  if (!"reference_normalization_status" %in% names(out)) {
    out$reference_normalization_status <- NA_character_
  }
  out$reference_normalization_method <-
    "biological_replicate_canonical_median_divided_by_configured_reference_v1"
  out
}

edu_build_intensity_report <- function(analysis) {
  overall <- edu_validate_intensity_table_for_report(
    analysis,
    "edu_positive_population_intensity_acquisition",
    "edu_positive_population_intensity",
    "positive_population_edu_bgsub_median", character(),
    "whole_computed_positive_eligible_population", 1L
  )
  regional <- edu_validate_intensity_table_for_report(
    analysis,
    "edu_positive_cell_regional_intensity_acquisition",
    "edu_positive_cell_regional_intensity",
    "positive_cell_regional_edu_bgsub_median",
    c("phase", "phase_label", "phase_index"),
    "computed_positive_eligible_cells_in_dna_region", 3L,
    category_col = "phase", category_values = c("early", "mid", "late")
  )
  overall$population_label <- "All computed EdU-positive"
  overall_relative <- edu_reference_relative_intensity(
    overall, "positive_population_edu_bgsub_median", analysis$sample_manifest
  )
  regional_relative <- edu_reference_relative_intensity(
    regional, "positive_cell_regional_edu_bgsub_median",
    analysis$sample_manifest, category_col = "phase"
  )
  list(
    schema_version = "edu-intensity-report-1.1.0",
    overall_canonical_biological_replicate = overall,
    regional_canonical_biological_replicate = regional,
    overall_biological_replicate = overall_relative,
    regional_biological_replicate = regional_relative,
    panels = list(
      all_computed_positive = edu_intensity_plot(
        overall_relative, "reference_relative_intensity", "population_label",
        "All computed EdU-positive", "Population",
        "Median background-subtracted EdU fluorescence\n(relative to configured matched reference)", analysis
      ),
      early_mid_late_s = edu_intensity_plot(
        regional_relative, "reference_relative_intensity", "phase_label",
        c("Early S", "Mid S", "Late S"), "S-phase region",
        "Median background-subtracted EdU fluorescence\n(relative to configured matched reference)", analysis
      )
    ),
    provenance = list(
      overall_source = "canonical_edu_positive_population_intensity_biological_replicate",
      regional_source = "canonical_edu_positive_cell_regional_intensity_biological_replicate",
      signal = "background_subtracted",
      aggregation = "unweighted technical-acquisition mean within biological replicate",
      reference_normalization = "canonical biological-replicate median divided by the explicitly configured matched reference"
    )
  )
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
  positivity <- edu_build_positivity_report(analysis)
  intensity <- edu_build_intensity_report(analysis)
  structure(list(
    schema_version = "edu-pseudocolor-output-contract-1.2.0",
    panels = stats::setNames(lapply(panels, `[[`, "plot"), expected),
    display_offset_qc = offset_qc, panel_qc = panel_qc,
    positivity = positivity,
    intensity = intensity,
    provenance = list(
      display_offset_method = "raw_negative_median_minus_corrected_negative_median_v1",
      analytical_values_mutated = FALSE,
      panel_identity_key = "manifest_prefix",
      density_palette = "refined_density_palette_v1"
    )
  ), class = "edu_pseudocolor_output_contract")
}
