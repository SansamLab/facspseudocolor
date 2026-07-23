# PH3 normalization, quantitation, and plotting ------------------------------

#' Normalize DNA for one PH3 sample using its FlowJo G1 gate
#'
#' PH3 positivity is supplied by the user-defined FlowJo population. This
#' function performs no positivity threshold fitting and does not modify target
#' intensities.
#'
#' @param events All Single Cell events.
#' @param g1_events G1-gated events for the same sample.
#' @param ph3_positive_events Events in the user-defined pH3-positive gate.
#' @param dna_channel Exact DNA channel column.
#' @param target_channel Exact pH3 channel column.
#' @param dna_2n_value Value to which the G1 DNA anchor is mapped.
#' @param g1_anchor Either `"median"` or `"mode"`.
#' @param sample_id Label used in errors and provenance.
#'
#' @return Normalized all-cell, G1, and pH3-positive tables plus the DNA anchor.
#' @export
normalize_ph3 <- function(
    events, g1_events, ph3_positive_events, dna_channel, target_channel,
    dna_2n_value = 1000, g1_anchor = "median", sample_id = "sample"
) {
  events <- read_facs_sample(
    events, dna_channel, target_channel, paste(sample_id, "Single Cells"), 2L
  )
  g1_events <- read_facs_sample(
    g1_events, dna_channel, target_channel, paste(sample_id, "G1"), 2L
  )
  ph3_positive_events <- read_facs_sample(
    ph3_positive_events, dna_channel, target_channel,
    paste(sample_id, "pH3 positive"), 0L
  )
  if (!config_scalar_number(dna_2n_value) || dna_2n_value <= 0) {
    stop("`dna_2n_value` must be a positive finite number.", call. = FALSE)
  }
  g1_anchor <- match.arg(g1_anchor, c("median", "mode"))
  anchor <- if (g1_anchor == "median") {
    stats::median(g1_events[[dna_channel]], na.rm = TRUE)
  } else {
    g1_target_anchor(g1_events[[dna_channel]], "mode")
  }
  if (!is.finite(anchor) || anchor <= 0) {
    stop("Invalid G1 DNA anchor for ", sample_id, ".", call. = FALSE)
  }
  normalize_table <- function(table) {
    table$target_raw <- table[[target_channel]]
    table$dna_norm <- table[[dna_channel]] / anchor * dna_2n_value
    table$target_norm <- table$target_raw
    table$target_bgsub <- table$target_raw
    table
  }
  list(
    data = normalize_table(events),
    g1 = normalize_table(g1_events),
    ph3_positive = normalize_table(ph3_positive_events),
    g1_dna_anchor = anchor,
    dna_normalization_factor = dna_2n_value / anchor,
    g1_anchor_method = g1_anchor,
    normalization_method = "g1_dna_only"
  )
}

ph3_gate_table <- function(config) {
  ranges <- list(
    "G1" = config$g1_x_range,
    "Early S" = config$s_phase_bins$early,
    "Mid S" = config$s_phase_bins$mid,
    "Late S" = config$s_phase_bins$late,
    "G2/M" = config$g2m_x_range
  )
  data.frame(
    gate = names(ranges),
    gate_index = seq_along(ranges),
    xmin = vapply(ranges, function(x) as.numeric(unlist(x))[[1]], numeric(1)),
    xmax = vapply(ranges, function(x) as.numeric(unlist(x))[[2]], numeric(1)),
    stringsAsFactors = FALSE
  )
}

assign_ph3_dna_phase <- function(dna, gates) {
  assignment <- rep("Unassigned", length(dna))
  for (i in seq_len(nrow(gates))) {
    upper_test <- if (i == nrow(gates)) dna <= gates$xmax[[i]] else
      dna < gates$xmax[[i]]
    selected <- is.finite(dna) & dna >= gates$xmin[[i]] & upper_test
    assignment[selected] <- gates$gate[[i]]
  }
  factor(assignment, levels = c(gates$gate, "Unassigned"))
}

#' Quantify user-gated pH3-positive events
#'
#' Every percentage uses all Single Cell events from the same sample as its
#' denominator. PH3-positive events outside the explicit DNA gates, including
#' nonfinite DNA values, are retained as `Unassigned`.
#'
#' @param analysis A PH3-mode `facs_analysis`.
#' @return The updated analysis with `quantitation$ph3`.
#' @export
quantify_ph3 <- function(analysis) {
  validate_analysis_object(analysis)
  if (!identical(analysis$config$plot_type, "ph3")) {
    stop("`quantify_ph3()` requires a PH3-mode analysis.", call. = FALSE)
  }
  manifest <- analysis$sample_manifest
  gates <- ph3_gate_table(analysis$config)
  summaries <- vector("list", nrow(manifest))
  phase_rows <- vector("list", nrow(manifest))
  sensitivity_rows <- vector("list", nrow(manifest))
  delta <- analysis$config$ph3_boundary_sensitivity_fraction *
    analysis$config$dna_2n_value
  shift_label <- format(
    100 * analysis$config$ph3_boundary_sensitivity_fraction,
    trim = TRUE, scientific = FALSE
  )

  for (i in seq_len(nrow(manifest))) {
    sample <- analysis$normalized_data[[i]]
    complete_n <- nrow(sample$data)
    positive_n <- nrow(sample$ph3_positive)
    if (complete_n <= 0) {
      stop("Single Cell denominator is zero for ", manifest$prefix[[i]], ".",
           call. = FALSE)
    }
    if (positive_n > complete_n) {
      stop(
        "pH3-positive count exceeds Single Cell count for ",
        manifest$prefix[[i]], ".", call. = FALSE
      )
    }
    assignment <- assign_ph3_dna_phase(sample$ph3_positive$dna_norm, gates)
    analysis$normalized_data[[i]]$ph3_positive$ph3_phase <- assignment
    counts <- table(assignment)
    counts <- as.integer(counts[levels(assignment)])
    summaries[[i]] <- data.frame(
      replicate = manifest$replicate[[i]],
      replicate_index = manifest$replicate_index[[i]],
      condition = manifest$condition[[i]],
      condition_index = manifest$condition_index[[i]],
      prefix = manifest$prefix[[i]],
      single_cell_n = complete_n,
      ph3_positive_n = positive_n,
      assigned_n = sum(counts[seq_len(nrow(gates))]),
      unassigned_n = counts[[nrow(gates) + 1L]],
      ph3_positive_percent = 100 * positive_n / complete_n,
      unassigned_percent = 100 * counts[[nrow(gates) + 1L]] / complete_n,
      stringsAsFactors = FALSE
    )
    phase_rows[[i]] <- data.frame(
      replicate = manifest$replicate[[i]],
      replicate_index = manifest$replicate_index[[i]],
      condition = manifest$condition[[i]],
      condition_index = manifest$condition_index[[i]],
      prefix = manifest$prefix[[i]],
      gate = levels(assignment),
      gate_index = seq_along(levels(assignment)),
      gate_count = counts,
      denominator = complete_n,
      phase_percent = 100 * counts / complete_n,
      stringsAsFactors = FALSE
    )

    g2m <- gates[gates$gate == "G2/M", , drop = FALSE]
    variants <- data.frame(
      variant = c(
        "Configured",
        paste0("Lower -", shift_label, "% of 2N"),
        paste0("Lower +", shift_label, "% of 2N"),
        paste0("Upper -", shift_label, "% of 2N"),
        paste0("Upper +", shift_label, "% of 2N")
      ),
      lower = c(g2m$xmin, g2m$xmin - delta, g2m$xmin + delta,
                g2m$xmin, g2m$xmin),
      upper = c(g2m$xmax, g2m$xmax, g2m$xmax,
                g2m$xmax - delta, g2m$xmax + delta),
      stringsAsFactors = FALSE
    )
    dna <- sample$ph3_positive$dna_norm
    variant_counts <- vapply(seq_len(nrow(variants)), function(j) {
      sum(is.finite(dna) & dna >= variants$lower[[j]] &
            dna <= variants$upper[[j]])
    }, integer(1))
    sensitivity_rows[[i]] <- data.frame(
      replicate = manifest$replicate[[i]],
      replicate_index = manifest$replicate_index[[i]],
      condition = manifest$condition[[i]],
      condition_index = manifest$condition_index[[i]],
      prefix = manifest$prefix[[i]],
      variants,
      g2m_count = variant_counts,
      denominator = complete_n,
      g2m_percent = 100 * variant_counts / complete_n,
      stringsAsFactors = FALSE
    )
  }
  analysis$quantitation$ph3 <- list(
    sample_summary = do.call(rbind, summaries),
    phase_percentages = do.call(rbind, phase_rows),
    gates = gates,
    boundary_sensitivity = do.call(rbind, sensitivity_rows),
    denominator = "all_single_cell_events",
    positivity_source = "user_defined_flowjo_gate"
  )
  rownames(analysis$quantitation$ph3$sample_summary) <- NULL
  rownames(analysis$quantitation$ph3$phase_percentages) <- NULL
  rownames(analysis$quantitation$ph3$boundary_sensitivity) <- NULL
  analysis
}

ph3_style <- function(analysis, appearance = NULL, appearance_file = NULL) {
  if (inherits(appearance, "facs_appearance")) appearance else
    resolve_facs_appearance(analysis, appearance, appearance_file)
}

#' Plot the overall pH3-positive percentage
#' @param analysis A quantified PH3 analysis.
#' @param appearance Optional appearance overrides.
#' @param appearance_file Optional appearance YAML.
#' @export
plot_ph3_overall <- function(analysis, appearance = NULL, appearance_file = NULL) {
  validate_analysis_object(analysis)
  values <- analysis$quantitation$ph3$sample_summary
  if (is.null(values)) stop("PH3 quantitation is unavailable.", call. = FALSE)
  style <- ph3_style(analysis, appearance, appearance_file)
  make_whole_population_barplot(
    values, error_bar = style$error_bar,
    y_title = "pH3-positive cells (% of Single Cells)",
    base_font_size = style$base_font_size, show_points = style$show_points,
    fill_colors = style$condition_colors,
    value_col = "ph3_positive_percent"
  )
}

#' Plot phase-specific pH3-positive percentages
#' @inheritParams plot_ph3_overall
#' @param geometry Either `"grouped"` or `"stacked"`.
#' @export
plot_ph3_phase <- function(
    analysis, geometry = c("grouped", "stacked"),
    appearance = NULL, appearance_file = NULL
) {
  geometry <- match.arg(geometry)
  values <- analysis$quantitation$ph3$phase_percentages
  if (is.null(values)) stop("PH3 quantitation is unavailable.", call. = FALSE)
  style <- ph3_style(analysis, appearance, appearance_file)
  if (geometry == "grouped") {
    return(make_phase_percentage_plot(
      values, error_bar = style$error_bar,
      base_font_size = style$base_font_size,
      show_points = style$show_points,
      fill_colors = style$condition_colors
    ) + ggplot2::labs(
      y = "pH3-positive cells\n(% of all Single Cells)"
    ))
  }
  summary <- summarize_across_replicates(
    values, "phase_percent",
    c("gate", "gate_index", "condition", "condition_index")
  )
  gate_levels <- unique(summary$gate[order(summary$gate_index)])
  condition_levels <- unique(summary$condition[order(summary$condition_index)])
  summary$gate <- factor(summary$gate, levels = gate_levels)
  summary$condition <- factor(summary$condition, levels = condition_levels)
  phase_colors <- stats::setNames(
    c("#4E79A7", "#76B7B2", "#59A14F", "#F28E2B", "#E15759", "#BDBDBD"),
    c("G1", "Early S", "Mid S", "Late S", "G2/M", "Unassigned")
  )
  ggplot2::ggplot(
    summary, ggplot2::aes(x = condition, y = mean, fill = gate)
  ) +
    ggplot2::geom_col(width = 0.72, color = "white", linewidth = 0.25) +
    ggplot2::scale_fill_manual(values = phase_colors, name = NULL) +
    ggplot2::scale_y_continuous(
      labels = function(x) paste0(x, "%"),
      expand = ggplot2::expansion(mult = c(0, 0.08))
    ) +
    ggplot2::labs(
      x = NULL,
      y = "pH3-positive cells\n(% of all Single Cells)"
    ) +
    ggplot2::theme_classic(base_size = style$base_font_size) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
      legend.position = "right"
    )
}

#' Plot PH3 DNA-gate assignments over all Single Cell events
#' @inheritParams plot_ph3_overall
#' @param seed Explicit downsampling seed.
#' @export
plot_ph3_diagnostic <- function(
    analysis, appearance = NULL, appearance_file = NULL, seed = 1L
) {
  validate_analysis_object(analysis)
  if (!identical(analysis$config$plot_type, "ph3")) {
    stop("PH3 diagnostic plotting requires PH3 mode.", call. = FALSE)
  }
  if (!config_scalar_number(as.numeric(seed))) {
    stop("`seed` must be one finite number.", call. = FALSE)
  }
  style <- ph3_style(analysis, appearance, appearance_file)
  manifest <- analysis$sample_manifest
  maximum <- analysis$config$gate_assignment_max_points
  old_seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (old_seed_exists) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (old_seed_exists) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  background <- list()
  positive <- list()
  for (i in seq_len(nrow(manifest))) {
    all_events <- analysis$normalized_data[[i]]$data
    pos_events <- analysis$normalized_data[[i]]$ph3_positive
    if (nrow(all_events) > maximum) {
      all_events <- all_events[sample.int(nrow(all_events), maximum), , drop = FALSE]
    }
    if (nrow(pos_events) > maximum) {
      pos_events <- pos_events[sample.int(nrow(pos_events), maximum), , drop = FALSE]
    }
    labels <- data.frame(
      replicate = manifest$replicate[[i]],
      condition = manifest$condition[[i]],
      stringsAsFactors = FALSE
    )
    background[[i]] <- cbind(labels[rep(1, nrow(all_events)), , drop = FALSE],
                             all_events[c("dna_norm", "target_norm")])
    positive[[i]] <- cbind(labels[rep(1, nrow(pos_events)), , drop = FALSE],
                           pos_events[c("dna_norm", "target_norm", "ph3_phase")])
  }
  background <- do.call(rbind, background)
  positive <- do.call(rbind, positive)
  gates <- analysis$quantitation$ph3$gates
  plot <- ggplot2::ggplot(background, ggplot2::aes(dna_norm, target_norm)) +
    ggplot2::geom_point(color = "grey82", size = 0.35, alpha = 0.45) +
    ggplot2::geom_point(
      data = positive, ggplot2::aes(color = ph3_phase),
      size = 0.65, alpha = 0.8
    ) +
    ggplot2::geom_vline(
      xintercept = unique(c(gates$xmin, gates$xmax)),
      color = style$gate_color, linetype = style$gate_linetype,
      linewidth = style$gate_linewidth
    ) +
    ggplot2::facet_grid(replicate ~ condition) +
    ggplot2::coord_cartesian(xlim = style$x_limits, ylim = style$y_limits) +
    ggplot2::labs(
      x = style$dna_axis_label, y = style$target_axis_label,
      color = "pH3-positive\nDNA assignment"
    ) +
    ggplot2::theme_classic(base_size = style$base_font_size)
  if (isTRUE(style$y_log10)) plot <- plot + ggplot2::scale_y_log10()
  plot
}

#' Plot G2/M boundary-sensitivity diagnostics
#' @inheritParams plot_ph3_overall
#' @export
plot_ph3_boundary_sensitivity <- function(
    analysis, appearance = NULL, appearance_file = NULL
) {
  values <- analysis$quantitation$ph3$boundary_sensitivity
  if (is.null(values)) stop("PH3 boundary sensitivity is unavailable.", call. = FALSE)
  style <- ph3_style(analysis, appearance, appearance_file)
  ggplot2::ggplot(
    values,
    ggplot2::aes(x = variant, y = g2m_percent, color = condition,
                 group = interaction(replicate, condition))
  ) +
    ggplot2::geom_line(alpha = 0.65) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::scale_color_manual(values = style$condition_colors, name = NULL) +
    ggplot2::scale_y_continuous(labels = function(x) paste0(x, "%")) +
    ggplot2::labs(
      x = "G2/M boundary variant",
      y = "pH3-positive G2/M cells\n(% of all Single Cells)"
    ) +
    ggplot2::theme_classic(base_size = style$base_font_size) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
}
