# =============================================================================
# phase_quantitation.R
#
# Cell-cycle phase gates and quantitation of the normalized signal, generalized
# from the Figure 1 analysis. Works for both modes (edu / poi) using the
# normalized columns produced by pseudocolor_helpers.R: `dna_norm` (2N at
# dna_2n_value, 4N at 2x) and `target_norm`.
#
# Provides:
#   * phase gates: G1, Early S, Mid S, Late S, G2/M
#   * a gate overlay for the pseudocolor panels
#   * per-phase median of the normalized signal (5 DNA windows)
#   * whole-population median of the normalized signal
#   * phase percentages (fraction of cells within each 2D gate)
#   * across-replicate summaries and Figure 1-style bar plots
#
# Requires ggplot2, scales, grDevices (namespaced).
# =============================================================================


# ---------------------------------------------------------------------------
# S-phase bins and phase windows / gates
# ---------------------------------------------------------------------------

default_s_phase_bins <- function(dna_2n_value = 1000) {
  list(
    early = c(0.910, 1.265) * dna_2n_value,
    mid   = c(1.265, 1.620) * dna_2n_value,
    late  = c(1.620, 1.975) * dna_2n_value
  )
}

# Normalize the s_phase_bins config into an ordered table of early/mid/late S
# DNA windows. Accepts a named list of length-2 ranges.
validate_s_phase_bins <- function(bins, dna_2n_value = 1000) {
  if (is.null(bins)) bins <- default_s_phase_bins(dna_2n_value)
  required <- c("early", "mid", "late")
  if (!is.list(bins) || !all(required %in% names(bins))) {
    stop("s_phase_bins must define named early, mid, and late ranges.")
  }
  bins <- lapply(bins[required], function(r) as.numeric(unlist(r)))
  if (any(vapply(bins, function(r) length(r) != 2 || any(!is.finite(r)) ||
                 r[[1]] >= r[[2]], logical(1)))) {
    stop("Each S-phase bin must contain two finite increasing numbers.")
  }
  data.frame(
    phase = required,
    phase_label = c("Early S", "Mid S", "Late S"),
    x_min = vapply(bins, `[[`, numeric(1), 1),
    x_max = vapply(bins, `[[`, numeric(1), 2),
    stringsAsFactors = FALSE
  )
}

# DNA windows (x-only) for per-phase signal medians. With include_g1_g2m = TRUE
# (poi default) returns G1, Early S, Mid S, Late S, G2/M; with FALSE (edu, where
# only S-phase EdU+ cells are informative) returns only Early/Mid/Late S.
make_phase_windows <- function(
    s_phase_bins = NULL,
    g1_x_range = NULL,
    g2m_x_range = NULL,
    dna_2n_value = 1000,
    include_g1_g2m = TRUE
) {
  s <- validate_s_phase_bins(s_phase_bins, dna_2n_value)
  if (include_g1_g2m) {
    if (is.null(g1_x_range))  g1_x_range  <- c(0.775, 1.225) * dna_2n_value
    if (is.null(g2m_x_range)) g2m_x_range <- c(1.675, 2.125) * dna_2n_value
    g1_x_range <- as.numeric(unlist(g1_x_range))
    g2m_x_range <- as.numeric(unlist(g2m_x_range))
    out <- data.frame(
      phase = c("g1", s$phase, "g2m"),
      phase_label = c("G1", s$phase_label, "G2/M"),
      x_min = c(g1_x_range[[1]], s$x_min, g2m_x_range[[1]]),
      x_max = c(g1_x_range[[2]], s$x_max, g2m_x_range[[2]]),
      stringsAsFactors = FALSE
    )
  } else {
    out <- data.frame(
      phase = s$phase, phase_label = s$phase_label,
      x_min = s$x_min, x_max = s$x_max, stringsAsFactors = FALSE
    )
  }
  out$phase_index <- seq_len(nrow(out))
  out
}

# Five 2D rectangular gates for counting cells (percentages). Non-replicating
# G1/G2M gates use the low-signal band; S gates use the high-signal band.
make_rectangular_phase_gates <- function(
    s_phase_bins = NULL,
    g1_x_range = NULL,
    g2m_x_range = NULL,
    negative_y_range = NULL,
    s_phase_y_range = NULL,
    dna_2n_value = 1000,
    y_limits = c(500, 80000)
) {
  if (is.null(g1_x_range))  g1_x_range  <- c(0.775, 1.225) * dna_2n_value
  if (is.null(g2m_x_range)) g2m_x_range <- c(1.675, 2.125) * dna_2n_value
  threshold <- 2.5 * dna_2n_value
  if (is.null(negative_y_range)) negative_y_range <- c(y_limits[[1]], threshold)
  if (is.null(s_phase_y_range))  s_phase_y_range  <- c(threshold, y_limits[[2]])
  g1_x_range <- as.numeric(unlist(g1_x_range))
  g2m_x_range <- as.numeric(unlist(g2m_x_range))
  negative_y_range <- as.numeric(unlist(negative_y_range))
  s_phase_y_range <- as.numeric(unlist(s_phase_y_range))

  s <- validate_s_phase_bins(s_phase_bins, dna_2n_value)
  n_s <- nrow(s)

  data.frame(
    gate = c("G1", s$phase_label, "G2/M"),
    gate_index = seq_len(n_s + 2),
    xmin = c(g1_x_range[[1]], s$x_min, g2m_x_range[[1]]),
    xmax = c(g1_x_range[[2]], s$x_max, g2m_x_range[[2]]),
    ymin = c(negative_y_range[[1]], rep(s_phase_y_range[[1]], n_s),
             negative_y_range[[1]]),
    ymax = c(negative_y_range[[2]], rep(s_phase_y_range[[2]], n_s),
             negative_y_range[[2]]),
    stringsAsFactors = FALSE
  )
}


# ---------------------------------------------------------------------------
# Gate overlay
# ---------------------------------------------------------------------------

# Add the rectangular gates (and optional labels) to a ggplot panel. Works for
# both the decorated plot and the naked plotgardener panel.
add_phase_gates_to_plot <- function(
    plot, gate_rectangles, color = "black", linetype = "dashed",
    linewidth = 0.5, show_labels = FALSE, label_size = 2.2, y_log10 = TRUE
) {
  gr <- gate_rectangles
  plot <- plot + ggplot2::geom_rect(
    data = gr,
    ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE, fill = NA, color = color,
    linetype = linetype, linewidth = linewidth
  )
  if (show_labels) {
    gr$label_x <- (gr$xmin + gr$xmax) / 2
    gr$label_y <- if (y_log10) sqrt(gr$ymin * gr$ymax) else (gr$ymin + gr$ymax) / 2
    plot <- plot + ggplot2::geom_text(
      data = gr, ggplot2::aes(x = label_x, y = label_y, label = gate),
      inherit.aes = FALSE, color = color, size = label_size, fontface = "bold"
    )
  }
  plot
}

make_computed_edu_gate_polygons <- function(
    sample, gate_rectangles, y_limits, offset = 0
) {
  required <- c("dna_norm", "edu_boundary_bgsub")
  if (!all(required %in% names(sample$data))) {
    stop("Computed EdU boundary columns are unavailable.", call. = FALSE)
  }
  x <- sample$data$dna_norm
  boundary <- sample$data$edu_boundary_bgsub + offset
  fit <- stats::lm(boundary ~ x)
  rows <- lapply(seq_len(nrow(gate_rectangles)), function(i) {
    gate <- gate_rectangles[i, ]
    bx <- c(gate$xmin, gate$xmax)
    by <- as.numeric(stats::predict(fit, newdata = data.frame(x = bx)))
    s_phase <- gate$gate %in% c("Early S", "Mid S", "Late S")
    if (s_phase) {
      px <- c(bx[[1]], bx[[1]], bx[[2]], bx[[2]])
      py <- c(by[[1]], y_limits[[2]], y_limits[[2]], by[[2]])
    } else {
      px <- c(bx[[1]], bx[[1]], bx[[2]], bx[[2]])
      py <- c(y_limits[[1]], by[[1]], by[[2]], y_limits[[1]])
    }
    data.frame(
      gate = gate$gate, gate_index = gate$gate_index,
      vertex = seq_along(px), x = px, y = py
    )
  })
  do.call(rbind, rows)
}

add_phase_gate_polygons_to_plot <- function(
    plot, polygons, color = "black", linetype = "dashed", linewidth = 0.5
) {
  plot + ggplot2::geom_polygon(
    data = polygons,
    ggplot2::aes(x = x, y = y, group = gate_index),
    inherit.aes = FALSE, fill = NA, color = color,
    linetype = linetype, linewidth = linewidth
  )
}


# ---------------------------------------------------------------------------
# Per-phase median of the normalized signal
# ---------------------------------------------------------------------------

calculate_phase_signal_medians <- function(
    data, windows, x_col = "dna_norm", y_col = "target_norm",
    minimum_events = 10, condition_label = "sample"
) {
  x <- data[[x_col]]; y <- data[[y_col]]
  rows <- lapply(seq_len(nrow(windows)), function(i) {
    sel <- is.finite(x) & is.finite(y) &
      x >= windows$x_min[[i]] & x < windows$x_max[[i]]
    vals <- y[sel]
    med <- if (length(vals) >= minimum_events) stats::median(vals) else NA_real_
    data.frame(
      phase = windows$phase[[i]], phase_label = windows$phase_label[[i]],
      phase_index = windows$phase_index[[i]],
      median_signal = med, n = length(vals), stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# Return the per-sample event table used for quantitation. "data" = all cells;
# "edu_positive" = the EdU+ (S gate) subset (required for edu-mode medians).
quant_source_data <- function(result, source_field, prefix) {
  if (source_field == "computed_edu_positive") {
    ev <- result$sample_data$data
    if (!"edu_computed_positive" %in% names(ev)) {
      legacy <- result$sample_data$edu_positive
      if (!is.null(legacy)) return(legacy)
      stop("Computed EdU boundary is unavailable for ", prefix, ".", call. = FALSE)
    }
    ev[ev$edu_computed_positive %in% TRUE, , drop = FALSE]
  } else if (source_field == "edu_positive") {
    ev <- result$sample_data$edu_positive
    if (is.null(ev)) {
      stop(paste0("EdU-positive events are not available for ", prefix,
                  " (needed for edu-mode medians). Ensure the edu_positive ",
                  "file/gate is exported."))
    }
    ev
  } else {
    result$sample_data$data
  }
}

# One row per (sample, phase). `source_field` selects all cells ("data", poi)
# or the EdU+ subset ("edu_positive", edu). `y_col` selects the signal column
# (e.g. "target_bgsub" for background-subtracted, "target_norm" for normalized).
collect_phase_signal_medians <- function(
    plot_results, sample_manifest, windows, source_field = "data",
    y_col = "target_norm") {
  rows <- lapply(seq_along(plot_results), function(i) {
    dat <- quant_source_data(plot_results[[i]], source_field,
                             sample_manifest$prefix[[i]])
    med <- calculate_phase_signal_medians(
      dat, windows, y_col = y_col,
      condition_label = sample_manifest$prefix[[i]])
    data.frame(
      replicate = sample_manifest$replicate[[i]],
      replicate_index = sample_manifest$replicate_index[[i]],
      technical_replicate = sample_manifest$technical_replicate[[i]],
      condition = sample_manifest$condition[[i]],
      condition_index = sample_manifest$condition_index[[i]],
      med, stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows); rownames(out) <- NULL; out
}

# Average independently processed technical acquisitions within each biological
# replicate. This prevents technical reruns from becoming independent points.
average_technical_replicates <- function(df, value_col,
                                         extra_group_cols = character()) {
  group_cols <- c("replicate", "replicate_index", "condition",
                  "condition_index", extra_group_cols)
  key <- do.call(paste, c(df[group_cols], sep = "\r"))
  rows <- lapply(split(seq_len(nrow(df)), key), function(idx) {
    values <- df[[value_col]][idx]
    values <- values[is.finite(values)]
    out <- df[idx[[1]], group_cols, drop = FALSE]
    out[[value_col]] <- if (length(values)) mean(values) else NA_real_
    out$technical_n <- length(values)
    if ("n" %in% names(df)) out$n <- sum(df$n[idx], na.rm = TRUE)
    out
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# Generic across-replicate summary of a per-sample value, grouped by one or two
# keys. `group_cols` are carried through; statistics computed on `value_col`.
summarize_across_replicates <- function(df, value_col, group_cols) {
  key <- do.call(paste, c(df[group_cols], sep = "\r"))
  split_idx <- split(seq_len(nrow(df)), key)
  rows <- lapply(split_idx, function(idx) {
    vals <- df[[value_col]][idx]; vals <- vals[is.finite(vals)]
    n <- length(vals)
    base <- df[idx[[1]], group_cols, drop = FALSE]
    data.frame(
      base,
      replicate_n = n,
      mean = if (n > 0) mean(vals) else NA_real_,
      sd = if (n > 1) stats::sd(vals) else NA_real_,
      sem = if (n > 1) stats::sd(vals) / sqrt(n) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows); rownames(out) <- NULL; out
}

# Normalize a per-sample value to a reference condition WITHIN each biological
# replicate (and, if `phase_col` is given, within each phase). Adds a `ratio`
# column = value / reference value; the reference condition becomes 1.
add_reference_ratio <- function(df, value_col, reference_condition,
                                phase_col = NULL) {
  if (!reference_condition %in% df$condition) {
    stop(paste0("Reference condition '", reference_condition,
                "' is not among the quantified samples. Set ",
                "`quant_reference_condition` to a displayed condition."))
  }
  key_cols <- c("replicate_index", phase_col)
  key <- do.call(paste, c(df[key_cols], sep = "\r"))
  df$ratio <- NA_real_
  for (k in unique(key)) {
    idx <- which(key == k)
    ref_idx <- idx[df$condition[idx] == reference_condition]
    if (length(ref_idx) != 1) next
    ref_val <- df[[value_col]][[ref_idx]]
    if (is.finite(ref_val) && ref_val != 0) {
      df$ratio[idx] <- df[[value_col]][idx] / ref_val
    }
  }
  df
}

# Default palette: evenly spaced viridis-style colors.
condition_fill_palette <- function(n) {
  anchors <- c("#440154", "#414487", "#2A788E", "#22A884", "#7AD151", "#FDE725")
  grDevices::colorRampPalette(anchors)(n)
}

# Resolve one fill color per condition level. `fill_colors` may be a named
# vector/list mapping condition label -> color (any not listed fall back to a
# neutral grey); if NULL, the default palette is used.
resolve_condition_colors <- function(cond_levels, fill_colors = NULL) {
  if (is.null(fill_colors)) {
    return(stats::setNames(condition_fill_palette(length(cond_levels)),
                           cond_levels))
  }
  fc <- unlist(fill_colors)
  out <- unname(fc[cond_levels])
  out[is.na(out)] <- "#B5B5B5"
  stats::setNames(out, cond_levels)
}

# Grouped bar plot: x = phase, fill = condition. `show_points` toggles replicate
# points; `fill_colors` is an optional named condition -> color map.
make_phase_signal_barplot <- function(
    phase_medians, error_bar = c("sd", "sem", "none"),
    y_title = "Median normalized signal", base_font_size = 11,
    show_points = TRUE, fill_colors = NULL,
    value_col = "median_signal", ref_line = NULL
) {
  error_bar <- match.arg(error_bar)
  s <- summarize_across_replicates(
    phase_medians, value_col, c("phase_label", "phase_index",
                                "condition", "condition_index"))
  phase_levels <- unique(s$phase_label[order(s$phase_index)])
  cond_levels <- unique(s$condition[order(s$condition_index)])
  s$phase_label <- factor(s$phase_label, levels = phase_levels)
  s$condition <- factor(s$condition, levels = cond_levels)
  s$error <- switch(error_bar, sd = s$sd, sem = s$sem,
                    none = rep(NA_real_, nrow(s)))
  pts <- phase_medians
  pts$phase_label <- factor(pts$phase_label, levels = phase_levels)
  pts$condition <- factor(pts$condition, levels = cond_levels)
  dodge <- 0.82

  plot <- ggplot2::ggplot(
    s, ggplot2::aes(x = phase_label, y = mean, fill = condition))
  if (!is.null(ref_line)) {
    plot <- plot + ggplot2::geom_hline(
      yintercept = ref_line, linetype = "dashed", color = "grey45",
      linewidth = 0.4)
  }
  plot <- plot +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = dodge),
                      width = 0.76, color = "white", linewidth = 0.25) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = pmax(0, mean - error), ymax = mean + error),
      position = ggplot2::position_dodge(width = dodge),
      width = 0.16, linewidth = 0.45, na.rm = TRUE)
  if (isTRUE(show_points)) {
    plot <- plot + ggplot2::geom_point(
      data = pts, ggplot2::aes(x = phase_label, y = .data[[value_col]],
                               group = condition),
      inherit.aes = FALSE, shape = 3, size = 1.6, stroke = 0.6, color = "black",
      position = ggplot2::position_jitterdodge(
        jitter.width = 0.05, jitter.height = 0, dodge.width = dodge, seed = 1))
  }
  plot +
    ggplot2::scale_fill_manual(
      values = resolve_condition_colors(cond_levels, fill_colors), name = NULL) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.10))) +
    ggplot2::labs(x = "Cell-cycle phase", y = y_title) +
    ggplot2::theme_classic(base_size = base_font_size) +
    ggplot2::theme(legend.position = "right",
                   axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
}

# Line plot: x = phase, one line per condition, error bars across replicates.
# Same summarized data as make_phase_signal_barplot, drawn as connected lines.
make_phase_signal_lineplot <- function(
    phase_medians, error_bar = c("sd", "sem", "none"),
    y_title = "Median normalized signal", base_font_size = 11,
    show_points = TRUE, fill_colors = NULL,
    value_col = "median_signal", ref_line = NULL
) {
  error_bar <- match.arg(error_bar)
  s <- summarize_across_replicates(
    phase_medians, value_col, c("phase_label", "phase_index",
                                "condition", "condition_index"))
  phase_levels <- unique(s$phase_label[order(s$phase_index)])
  cond_levels <- unique(s$condition[order(s$condition_index)])
  s$phase_label <- factor(s$phase_label, levels = phase_levels)
  s$condition <- factor(s$condition, levels = cond_levels)
  s$error <- switch(error_bar, sd = s$sd, sem = s$sem,
                    none = rep(NA_real_, nrow(s)))
  cols <- resolve_condition_colors(cond_levels, fill_colors)

  plot <- ggplot2::ggplot(
    s, ggplot2::aes(x = phase_label, y = mean, color = condition,
                    group = condition))
  if (!is.null(ref_line)) {
    plot <- plot + ggplot2::geom_hline(
      yintercept = ref_line, linetype = "dashed", color = "grey45",
      linewidth = 0.4)
  }
  plot <- plot +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = pmax(0, mean - error), ymax = mean + error),
      width = 0.14, linewidth = 0.5, na.rm = TRUE)
  if (isTRUE(show_points)) {
    plot <- plot + ggplot2::geom_point(size = 2)
  }
  plot +
    ggplot2::scale_color_manual(values = cols, name = NULL) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0.03, 0.08))) +
    ggplot2::labs(x = "Phase", y = y_title) +
    ggplot2::theme_classic(base_size = base_font_size)
}


# ---------------------------------------------------------------------------
# Whole-population median of the normalized signal
# ---------------------------------------------------------------------------

collect_whole_population_medians <- function(
    plot_results, sample_manifest, source_field = "data",
    y_col = "target_norm") {
  rows <- lapply(seq_along(plot_results), function(i) {
    dat <- quant_source_data(plot_results[[i]], source_field,
                             sample_manifest$prefix[[i]])
    v <- dat[[y_col]]
    v <- v[is.finite(v)]
    data.frame(
      replicate = sample_manifest$replicate[[i]],
      replicate_index = sample_manifest$replicate_index[[i]],
      technical_replicate = sample_manifest$technical_replicate[[i]],
      condition = sample_manifest$condition[[i]],
      condition_index = sample_manifest$condition_index[[i]],
      median_signal = if (length(v)) stats::median(v) else NA_real_,
      n = length(v), stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows); rownames(out) <- NULL; out
}

# Bar plot: x = condition, bars colored per condition, optional replicate points.
make_whole_population_barplot <- function(
    whole_medians, error_bar = c("sd", "sem", "none"),
    y_title = "Median normalized signal", base_font_size = 11,
    show_points = TRUE, fill_colors = NULL,
    value_col = "median_signal", ref_line = NULL
) {
  error_bar <- match.arg(error_bar)
  s <- summarize_across_replicates(whole_medians, value_col,
                                   c("condition", "condition_index"))
  cond_levels <- unique(s$condition[order(s$condition_index)])
  s$condition <- factor(s$condition, levels = cond_levels)
  s$error <- switch(error_bar, sd = s$sd, sem = s$sem,
                    none = rep(NA_real_, nrow(s)))
  pts <- whole_medians
  pts$condition <- factor(pts$condition, levels = cond_levels)

  plot <- ggplot2::ggplot(s, ggplot2::aes(x = condition, y = mean,
                                          fill = condition))
  if (!is.null(ref_line)) {
    plot <- plot + ggplot2::geom_hline(
      yintercept = ref_line, linetype = "dashed", color = "grey45",
      linewidth = 0.4)
  }
  plot <- plot +
    ggplot2::geom_col(width = 0.70, color = "black", linewidth = 0.4) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = pmax(0, mean - error), ymax = mean + error),
      width = 0.18, linewidth = 0.5, na.rm = TRUE)
  if (isTRUE(show_points)) {
    plot <- plot + ggplot2::geom_point(
      data = pts, ggplot2::aes(x = condition, y = .data[[value_col]]),
      inherit.aes = FALSE, shape = 21, fill = "white", color = "black",
      size = 2.0, stroke = 0.5,
      position = ggplot2::position_jitter(width = 0.07, height = 0, seed = 1))
  }
  plot +
    ggplot2::scale_fill_manual(
      values = resolve_condition_colors(cond_levels, fill_colors),
      guide = "none") +
    ggplot2::scale_y_continuous(labels = scales::label_comma(),
                                expand = ggplot2::expansion(mult = c(0, 0.10))) +
    ggplot2::labs(x = NULL, y = y_title) +
    ggplot2::theme_classic(base_size = base_font_size) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
}


# ---------------------------------------------------------------------------
# Phase percentages (2D gate counts)
# ---------------------------------------------------------------------------

collect_phase_counts <- function(plot_results, sample_manifest, gate_rectangles) {
  rows <- lapply(seq_along(plot_results), function(i) {
    ev <- plot_results[[i]]$sample_data$data
    x <- ev$dna_norm; y <- ev$target_norm
    computed_edu <- "edu_computed_positive" %in% names(ev)
    counts <- vapply(seq_len(nrow(gate_rectangles)), function(g) {
      gr <- gate_rectangles[g, ]
      selected <- is.finite(x) & is.finite(y) &
        x >= gr$xmin & x < gr$xmax
      if (computed_edu) {
        if (gr$gate %in% c("Early S", "Mid S", "Late S")) {
          selected <- selected & ev$edu_computed_positive %in% TRUE
        } else {
          selected <- selected & ev$edu_computed_positive %in% FALSE
        }
      } else {
        selected <- selected & y >= gr$ymin & y < gr$ymax
      }
      sum(selected)
    }, integer(1))
    denom <- sum(counts)
    if (denom <= 0) {
      stop(paste("No events in the phase gates for",
                 sample_manifest$prefix[[i]]))
    }
    data.frame(
      replicate = sample_manifest$replicate[[i]],
      replicate_index = sample_manifest$replicate_index[[i]],
      technical_replicate = sample_manifest$technical_replicate[[i]],
      condition = sample_manifest$condition[[i]],
      condition_index = sample_manifest$condition_index[[i]],
      gate = gate_rectangles$gate, gate_index = gate_rectangles$gate_index,
      gate_count = counts, phase_denominator = denom,
      phase_percent = 100 * counts / denom, stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows); rownames(out) <- NULL; out
}

# EdU output-contract classification and tables -----------------------------

edu_output_schema_version <- function() 2L

edu_contract_ranges <- function(config) {
  windows <- make_phase_windows(
    s_phase_bins = config$s_phase_bins,
    g1_x_range = config$g1_x_range,
    g2m_x_range = config$g2m_x_range,
    dna_2n_value = config$dna_2n_value,
    include_g1_g2m = TRUE
  )
  list(
    windows = windows,
    g1 = unname(c(windows$x_min[[1]], windows$x_max[[1]])),
    g2m = unname(c(windows$x_min[[5]], windows$x_max[[5]])),
    s = windows[2:4, , drop = FALSE]
  )
}

classify_edu_events <- function(data, config, display_offset) {
  required <- c(
    "dna_norm", "target_norm", "target_bgsub", "edu_computed_positive"
  )
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop(
      "EdU output classification requires event field(s): ",
      paste(missing, collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (!config_scalar_number(display_offset) || display_offset < 0) {
    stop("`display_offset` must be one nonnegative finite number.",
         call. = FALSE)
  }

  ranges <- edu_contract_ranges(config)
  x <- data$dna_norm
  quantitative_signal <- data$target_bgsub
  computed_positive <- data$edu_computed_positive
  positivity_known <- !is.na(computed_positive) &
    computed_positive %in% c(TRUE, FALSE)
  dna_eligible <- is.finite(x) & positivity_known
  composition_eligible <- dna_eligible
  # The historical implementation also required a finite target_norm even
  # though that display/divided coordinate does not define the new metrics.
  historical_composition_eligible <- dna_eligible & is.finite(data$target_norm)
  quantitative_eligible <- dna_eligible & is.finite(quantitative_signal)
  positive_population_intensity_eligible <- positivity_known &
    computed_positive %in% TRUE & is.finite(quantitative_signal)

  region <- rep(NA_character_, nrow(data))
  for (i in seq_len(nrow(ranges$s))) {
    selected <- is.finite(x) & x >= ranges$s$x_min[[i]] &
      x < ranges$s$x_max[[i]]
    region[selected] <- ranges$s$phase[[i]]
  }

  historical_assignment <- rep(NA_character_, nrow(data))
  biological_assignment <- rep(NA_character_, nrow(data))
  negative <- positivity_known & computed_positive %in% FALSE
  positive <- positivity_known & computed_positive %in% TRUE
  historical_assignment[historical_composition_eligible & negative &
    x >= ranges$g1[[1]] & x < ranges$g1[[2]]] <- "g1"
  historical_assignment[historical_composition_eligible & positive &
    region %in% c("early", "mid", "late")] <- region[
      historical_composition_eligible & positive &
        region %in% c("early", "mid", "late")
    ]
  historical_assignment[historical_composition_eligible & negative &
    x >= ranges$g2m[[1]] & x < ranges$g2m[[2]]] <- "g2m"

  biological_assignment[composition_eligible & negative &
    x >= ranges$g1[[1]] & x < ranges$g1[[2]]] <- "g1"
  biological_assignment[composition_eligible & positive &
    region %in% c("early", "mid", "late")] <- region[
      composition_eligible & positive & region %in% c("early", "mid", "late")
    ]
  biological_assignment[composition_eligible & negative &
    x >= ranges$g2m[[1]] & x < ranges$g2m[[2]]] <- "g2m"

  negative_s <- composition_eligible & negative &
    x >= ranges$g1[[2]] & x < ranges$g2m[[1]]
  six_gate_assignment <- biological_assignment
  six_gate_assignment[negative_s] <- "negative_s"
  display_offset_applied <- identical(
    config$pseudocolor_signal, "background_subtracted"
  )
  actual_display_offset <- if (display_offset_applied) display_offset else 0
  display_signal <- if (display_offset_applied) {
    quantitative_signal + actual_display_offset
  } else {
    data$target_norm
  }
  display_transform <- if (display_offset_applied) {
    "background_subtracted_plus_offset"
  } else {
    "legacy_background_divided"
  }

  data.frame(
    event_row = seq_len(nrow(data)),
    composition_eligible = composition_eligible,
    historical_composition_eligible = historical_composition_eligible,
    positivity_eligible = dna_eligible,
    regional_intensity_eligible = quantitative_eligible & positive,
    positive_population_intensity_eligible =
      positive_population_intensity_eligible,
    computed_positive = computed_positive,
    s_region = region,
    historical_five_gate_assignment = historical_assignment,
    negative_s = negative_s,
    six_gate_assignment = six_gate_assignment,
    unassigned = composition_eligible & is.na(six_gate_assignment),
    positivity_unassigned = dna_eligible & is.na(six_gate_assignment),
    quantitative_signal = quantitative_signal,
    display_signal = display_signal,
    display_transform = rep(display_transform, nrow(data)),
    display_offset = rep(actual_display_offset, nrow(data)),
    display_offset_applied = rep(display_offset_applied, nrow(data)),
    dna_norm = x,
    stringsAsFactors = FALSE
  )
}

edu_sample_metadata <- function(manifest, i) {
  data.frame(
    replicate = manifest$replicate[[i]],
    replicate_index = manifest$replicate_index[[i]],
    technical_replicate = manifest$technical_replicate[[i]],
    condition = manifest$condition[[i]],
    condition_index = manifest$condition_index[[i]],
    stringsAsFactors = FALSE
  )
}

edu_metric_metadata <- function(
    source_population, display_offset, display_offset_applied,
    dna_interval_min = NA_real_, dna_interval_max = NA_real_,
    display_transform = if (display_offset_applied) {
      "background_subtracted_plus_offset"
    } else {
      "legacy_background_divided"
    }
) {
  data.frame(
    source_population = source_population,
    dna_interval_min = dna_interval_min,
    dna_interval_max = dna_interval_max,
    dna_interval_lower_inclusive = !is.na(dna_interval_min),
    dna_interval_upper_inclusive = FALSE,
    signal_transform = "background_subtracted",
    display_transform = display_transform,
    display_offset = display_offset,
    display_offset_applied = display_offset_applied,
    reference_normalization_status = "not_applied",
    aggregation_level = "acquisition",
    aggregation_method = "none",
    output_schema_version = edu_output_schema_version(),
    stringsAsFactors = FALSE
  )
}

collect_edu_composition <- function(
    classifications, manifest, assignment_col, categories, value_col,
    source_population, display_offset, denominator = c("assigned", "eligible")
) {
  denominator <- match.arg(denominator)
  labels <- c(
    g1 = "G1", early = "Early S", mid = "Mid S", late = "Late S",
    negative_s = "Negative S", g2m = "G2/M"
  )
  rows <- lapply(seq_along(classifications), function(i) {
    dat <- classifications[[i]]
    counts <- vapply(
      categories,
      function(category) sum(dat[[assignment_col]] %in% category, na.rm = TRUE),
      integer(1)
    )
    denominator_n <- if (denominator == "eligible") {
      sum(dat$composition_eligible)
    } else {
      sum(counts)
    }
    metric_unassigned <- dat$composition_eligible &
      is.na(dat[[assignment_col]])
    unassigned_n <- sum(metric_unassigned)
    value <- if (denominator_n > 0) 100 * counts / denominator_n else
      rep(NA_real_, length(counts))
    metric_status <- if (denominator_n > 0) "ok" else "zero_denominator"
    value_frame <- data.frame(value, stringsAsFactors = FALSE)
    names(value_frame) <- value_col
    cbind(
      edu_sample_metadata(manifest, i),
      data.frame(
        gate = unname(labels[categories]),
        gate_index = seq_along(categories),
        numerator_n = counts,
        denominator_n = rep(denominator_n, length(counts)),
        unassigned_n = rep(unassigned_n, length(counts)),
        unassigned_pct_of_eligible_single_cells = rep(
          if (sum(dat$composition_eligible) > 0) {
            100 * unassigned_n / sum(dat$composition_eligible)
          } else {
            NA_real_
          },
          length(counts)
        ),
        stringsAsFactors = FALSE
      ),
      value_frame,
      edu_metric_metadata(
        source_population, dat$display_offset[[1]],
        dat$display_offset_applied[[1]] %in% TRUE,
        display_transform = dat$display_transform[[1]]
      )[rep(1, length(counts)), , drop = FALSE],
      data.frame(
        metric_status = rep(metric_status, length(counts)),
        stringsAsFactors = FALSE
      )
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

collect_edu_regional_positivity <- function(
    classifications, manifest, config, display_offset
) {
  ranges <- edu_contract_ranges(config)$s
  rows <- lapply(seq_along(classifications), function(i) {
    dat <- classifications[[i]]
    values <- lapply(seq_len(nrow(ranges)), function(j) {
      in_region <- dat$positivity_eligible &
        dat$dna_norm >= ranges$x_min[[j]] & dat$dna_norm < ranges$x_max[[j]]
      denominator_n <- sum(in_region)
      numerator_n <- sum(in_region & dat$computed_positive %in% TRUE)
      data.frame(
        region = ranges$phase[[j]],
        region_label = ranges$phase_label[[j]],
        region_index = ranges$phase_index[[j]] - 1L,
        numerator_n = numerator_n,
        denominator_n = denominator_n,
        regional_edu_positive_pct = if (denominator_n > 0) {
          100 * numerator_n / denominator_n
        } else {
          NA_real_
        },
        edu_metric_metadata(
          "eligible_single_cells_in_dna_region", dat$display_offset[[1]],
          dat$display_offset_applied[[1]] %in% TRUE,
          ranges$x_min[[j]], ranges$x_max[[j]],
          display_transform = dat$display_transform[[1]]
        ),
        metric_status = if (denominator_n > 0) "ok" else "zero_denominator",
        stringsAsFactors = FALSE
      )
    })
    cbind(edu_sample_metadata(manifest, i), do.call(rbind, values))
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

collect_edu_overall_positivity <- function(
    classifications, manifest, config, display_offset
) {
  ranges <- edu_contract_ranges(config)
  lower <- ranges$g1[[1]]
  upper <- ranges$g2m[[2]]
  rows <- lapply(seq_along(classifications), function(i) {
    dat <- classifications[[i]]
    in_span <- dat$positivity_eligible & dat$dna_norm >= lower &
      dat$dna_norm < upper
    denominator_n <- sum(in_span)
    numerator_n <- sum(in_span & dat$computed_positive %in% TRUE)
    unassigned_n <- sum(in_span & dat$positivity_unassigned)
    cbind(
      edu_sample_metadata(manifest, i),
      data.frame(
        region = "g1_through_g2m",
        region_label = "G1 through G2/M",
        numerator_n = numerator_n,
        denominator_n = denominator_n,
        phase_unassigned_n = unassigned_n,
        overall_edu_positive_pct = if (denominator_n > 0) {
          100 * numerator_n / denominator_n
        } else {
          NA_real_
        },
        edu_metric_metadata(
          "eligible_single_cells_in_g1_through_g2m_span",
          dat$display_offset[[1]],
          dat$display_offset_applied[[1]] %in% TRUE, lower, upper,
          display_transform = dat$display_transform[[1]]
        ),
        metric_status = if (denominator_n > 0) "ok" else "zero_denominator",
        stringsAsFactors = FALSE
      )
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

collect_edu_positive_intensity <- function(
    classifications, manifest, config, display_offset,
    regional = TRUE, minimum_events = 10L
) {
  ranges <- edu_contract_ranges(config)$s
  rows <- lapply(seq_along(classifications), function(i) {
    dat <- classifications[[i]]
    display_applied <- dat$display_offset_applied[[1]] %in% TRUE
    if (isTRUE(regional)) {
      values <- lapply(seq_len(nrow(ranges)), function(j) {
        selected <- dat$regional_intensity_eligible &
          dat$dna_norm >= ranges$x_min[[j]] & dat$dna_norm < ranges$x_max[[j]]
        signal <- dat$quantitative_signal[selected]
        n <- length(signal)
        data.frame(
          phase = ranges$phase[[j]],
          phase_label = ranges$phase_label[[j]],
          phase_index = ranges$phase_index[[j]] - 1L,
          source_population_n = n,
          positive_cell_regional_edu_bgsub_median = if (n >= minimum_events) {
            stats::median(signal)
          } else {
            NA_real_
          },
          edu_metric_metadata(
            "computed_positive_eligible_cells_in_dna_region",
            dat$display_offset[[1]], display_applied,
            ranges$x_min[[j]], ranges$x_max[[j]],
            display_transform = dat$display_transform[[1]]
          ),
          metric_status = if (n >= minimum_events) "ok" else
            "insufficient_events",
          stringsAsFactors = FALSE
        )
      })
      value <- do.call(rbind, values)
    } else {
      selected <- dat$positive_population_intensity_eligible
      signal <- dat$quantitative_signal[selected]
      n <- length(signal)
      value <- data.frame(
        source_population_n = n,
        positive_population_edu_bgsub_median = if (n > 0) {
          stats::median(signal)
        } else {
          NA_real_
        },
        edu_metric_metadata(
          "whole_computed_positive_eligible_population",
          dat$display_offset[[1]], display_applied,
          display_transform = dat$display_transform[[1]]
        ),
        metric_status = if (n > 0) "ok" else "insufficient_events",
        stringsAsFactors = FALSE
      )
    }
    cbind(edu_sample_metadata(manifest, i), value)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

average_edu_metric_table <- function(df, value_col, extra_group_cols = character()) {
  group_cols <- c(
    "replicate", "replicate_index", "condition", "condition_index",
    extra_group_cols
  )
  key <- do.call(paste, c(df[group_cols], sep = "\r"))
  count_cols <- names(df)[grepl("(_n|^numerator_n$|^denominator_n$)$", names(df))]
  count_cols <- setdiff(count_cols, c("replicate_index", "condition_index"))
  averaged_metadata_cols <- intersect(
    "unassigned_pct_of_eligible_single_cells", names(df)
  )
  rows <- lapply(split(seq_len(nrow(df)), key), function(idx) {
    values <- df[[value_col]][idx]
    finite <- is.finite(values)
    keep <- setdiff(
      names(df),
      c(
        "technical_replicate", value_col, count_cols, averaged_metadata_cols,
        "metric_status"
      )
    )
    out <- df[idx[[1]], keep, drop = FALSE]
    out[[value_col]] <- if (any(finite)) mean(values[finite]) else NA_real_
    for (name in count_cols) out[[name]] <- sum(df[[name]][idx], na.rm = TRUE)
    for (name in averaged_metadata_cols) {
      metadata_values <- df[[name]][idx]
      metadata_values <- metadata_values[is.finite(metadata_values)]
      out[[name]] <- if (length(metadata_values)) mean(metadata_values) else
        NA_real_
    }
    out$technical_n <- sum(finite)
    statuses <- unique(df$metric_status[idx])
    out$technical_acquisition_n <- length(idx)
    out$non_ok_technical_n <- sum(df$metric_status[idx] != "ok")
    out$technical_metric_statuses <- paste(sort(statuses), collapse = ";")
    out$aggregation_level <- "biological_replicate_condition"
    out$aggregation_method <-
      "unweighted_mean_of_technical_acquisition_values"
    out$metric_status <- if (all(df$metric_status[idx] == "ok")) {
      "ok"
    } else if (any(finite)) {
      "partial"
    } else if (length(statuses) == 1L) {
      statuses[[1]]
    } else {
      "mixed"
    }
    out
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# Grouped bar plot: x = gate, fill = condition. `show_points` toggles replicate
# points; `fill_colors` is an optional named condition -> color map.
make_phase_percentage_plot <- function(
    phase_counts, error_bar = c("sd", "sem", "none"), base_font_size = 11,
    show_points = TRUE, fill_colors = NULL,
    y_title = "Cells within phase gate\n(% of five-gate total)",
    caption = NULL
) {
  error_bar <- match.arg(error_bar)
  s <- summarize_across_replicates(
    phase_counts, "phase_percent", c("gate", "gate_index",
                                    "condition", "condition_index"))
  gate_levels <- unique(s$gate[order(s$gate_index)])
  cond_levels <- unique(s$condition[order(s$condition_index)])
  s$gate <- factor(s$gate, levels = gate_levels)
  s$condition <- factor(s$condition, levels = cond_levels)
  s$error <- switch(error_bar, sd = s$sd, sem = s$sem,
                    none = rep(NA_real_, nrow(s)))
  pts <- phase_counts
  pts$gate <- factor(pts$gate, levels = gate_levels)
  pts$condition <- factor(pts$condition, levels = cond_levels)
  dodge <- 0.86

  plot <- ggplot2::ggplot(s, ggplot2::aes(x = gate, y = mean, fill = condition)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = dodge),
                      width = 0.80, color = "white", linewidth = 0.25) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = pmax(0, mean - error), ymax = mean + error),
      position = ggplot2::position_dodge(width = dodge),
      width = 0.14, linewidth = 0.45, na.rm = TRUE)
  if (isTRUE(show_points)) {
    plot <- plot + ggplot2::geom_point(
      data = pts, ggplot2::aes(x = gate, y = phase_percent, group = condition),
      inherit.aes = FALSE, shape = 3, size = 1.5, stroke = 0.6, color = "black",
      position = ggplot2::position_jitterdodge(
        jitter.width = 0.05, jitter.height = 0, dodge.width = dodge, seed = 1))
  }
  plot +
    ggplot2::scale_fill_manual(
      values = resolve_condition_colors(cond_levels, fill_colors), name = NULL) +
    ggplot2::scale_y_continuous(labels = function(x) paste0(x, "%"),
                                expand = ggplot2::expansion(mult = c(0, 0.10))) +
    ggplot2::labs(x = "Cell-cycle phase", y = y_title, caption = caption) +
    ggplot2::theme_classic(base_size = base_font_size) +
    ggplot2::theme(legend.position = "right",
                   axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
}


# ---------------------------------------------------------------------------
# DNA-content phase gates (poi) and gate-membership visualization
# ---------------------------------------------------------------------------

# For POI there is no S-phase marker on the y-axis, so cell-cycle phase is set by
# DNA content alone. Standard DNA gating: a G1 gate straddling the 2N peak, a
# G2/M gate straddling the 4N peak, and S as the region between them, split into
# equal Early/Mid/Late thirds. Full-height, contiguous, non-overlapping bands.
make_dna_phase_gates <- function(
    s_phase_bins = NULL, g1_x_range = NULL, g2m_x_range = NULL,
    dna_2n_value = 1000, y_limits = c(1, 1e5)
) {
  if (is.null(g1_x_range))  g1_x_range  <- c(0.775, 1.225) * dna_2n_value
  if (is.null(g2m_x_range)) g2m_x_range <- c(1.675, 2.125) * dna_2n_value
  g1_x_range <- as.numeric(unlist(g1_x_range))
  g2m_x_range <- as.numeric(unlist(g2m_x_range))
  s_lo <- g1_x_range[[2]]    # S starts where the G1 gate ends
  s_hi <- g2m_x_range[[1]]   # S ends where the G2/M gate starts
  if (s_lo >= s_hi) {
    stop("The G1 and G2/M DNA gates overlap; widen the gap via g1_x_range / g2m_x_range.")
  }
  b <- seq(s_lo, s_hi, length.out = 4)   # Early / Mid / Late S thirds

  data.frame(
    gate = c("G1", "Early S", "Mid S", "Late S", "G2/M"),
    gate_index = 1:5,
    xmin = c(g1_x_range[[1]], b[[1]], b[[2]], b[[3]], g2m_x_range[[1]]),
    xmax = c(g1_x_range[[2]], b[[2]], b[[3]], b[[4]], g2m_x_range[[2]]),
    ymin = y_limits[[1]], ymax = y_limits[[2]],
    stringsAsFactors = FALSE
  )
}

# Distinct colors for the phase gates (+ grey for ungated events).
phase_gate_colors <- function(gate_levels) {
  base <- c("G1" = "#7F7F7F", "Early S" = "#0B4BFF", "Mid S" = "#12BED0",
            "Late S" = "#A026A3", "G2/M" = "#E08214", "ungated" = "#DADADA")
  out <- base[gate_levels]
  out[is.na(out)] <- "#B0B0B0"
  stats::setNames(unname(out), gate_levels)
}

# Assign each event to the first gate that contains it; unmatched -> "ungated".
assign_events_to_gates <- function(
    data, gate_rectangles, x_col = "dna_norm", y_col = "target_norm"
) {
  x <- data[[x_col]]; y <- data[[y_col]]
  assigned <- rep(NA_character_, length(x))
  for (i in seq_len(nrow(gate_rectangles))) {
    gr <- gate_rectangles[i, ]
    sel <- is.na(assigned) & is.finite(x) & is.finite(y) &
      x >= gr$xmin & x < gr$xmax & y >= gr$ymin & y < gr$ymax
    assigned[sel] <- gr$gate
  }
  assigned[is.na(assigned)] <- "ungated"
  data$gate <- factor(assigned, levels = c(gate_rectangles$gate, "ungated"))
  data
}

# Collect gate-labeled events (downsampled) across displayed samples.
collect_gate_assignments <- function(
    plot_results, sample_manifest, gate_rectangles, source_field = "data",
    max_points = 4000
) {
  rows <- lapply(seq_along(plot_results), function(i) {
    dat <- quant_source_data(plot_results[[i]], source_field,
                             sample_manifest$prefix[[i]])
    if ("edu_computed_positive" %in% names(dat)) {
      assigned <- rep("ungated", nrow(dat))
      for (g in seq_len(nrow(gate_rectangles))) {
        gr <- gate_rectangles[g, ]
        selected <- is.finite(dat$dna_norm) &
          dat$dna_norm >= gr$xmin & dat$dna_norm < gr$xmax
        if (gr$gate %in% c("Early S", "Mid S", "Late S")) {
          selected <- selected & dat$edu_computed_positive %in% TRUE
        } else {
          selected <- selected & dat$edu_computed_positive %in% FALSE
        }
        assigned[selected & assigned == "ungated"] <- gr$gate
      }
      dat$gate <- factor(assigned, levels = c(gate_rectangles$gate, "ungated"))
    } else {
      dat <- assign_events_to_gates(dat, gate_rectangles)
    }
    if (nrow(dat) > max_points) {
      dat <- dat[sample.int(nrow(dat), max_points), , drop = FALSE]
    }
    data.frame(
      replicate = sample_manifest$replicate[[i]],
      replicate_index = sample_manifest$replicate_index[[i]],
      condition = sample_manifest$condition[[i]],
      condition_index = sample_manifest$condition_index[[i]],
      dna_norm = dat$dna_norm, target_norm = dat$target_norm,
      gate = dat$gate, stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# Faceted scatter (replicate x condition) coloring each event by its phase gate.
make_gate_assignment_plot <- function(
    assignments, gate_rectangles, x_limits, y_limits, y_log10 = TRUE,
    dna_2n_value = 1000, point_size = 0.25, base_font_size = 9
) {
  cond_levels <- unique(
    assignments$condition[order(assignments$condition_index)])
  rep_levels <- unique(
    assignments$replicate[order(assignments$replicate_index)])
  assignments$condition <- factor(assignments$condition, levels = cond_levels)
  assignments$replicate <- factor(assignments$replicate, levels = rep_levels)
  gate_levels <- c(gate_rectangles$gate, "ungated")
  assignments$gate <- factor(as.character(assignments$gate), levels = gate_levels)

  p <- ggplot2::ggplot(
    assignments,
    ggplot2::aes(x = dna_norm, y = target_norm, color = gate)) +
    ggplot2::geom_point(size = point_size, stroke = 0, shape = 16) +
    ggplot2::scale_color_manual(values = phase_gate_colors(gate_levels),
                                name = "Phase gate", drop = FALSE) +
    ggplot2::geom_vline(xintercept = dna_2n_value, linetype = "dashed",
                        color = "grey55", linewidth = 0.3) +
    ggplot2::scale_x_continuous(breaks = c(dna_2n_value, 2 * dna_2n_value),
                                labels = c("2N", "4N"), minor_breaks = NULL) +
    ggplot2::coord_cartesian(xlim = x_limits, ylim = y_limits) +
    ggplot2::facet_grid(replicate ~ condition) +
    ggplot2::labs(x = "DNA content", y = "Signal") +
    ggplot2::theme_bw(base_size = base_font_size) +
    ggplot2::theme(
      legend.position = "right",
      panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::guides(color = ggplot2::guide_legend(
      override.aes = list(size = 2.5)))
  if (y_log10) p <- p + ggplot2::scale_y_log10()
  p
}
