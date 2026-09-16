facs_report_input <- function(config = NULL, analysis_rds = NULL) {
  has_config <- is.character(config) && length(config) == 1L && nzchar(config)
  has_rds <- is.character(analysis_rds) && length(analysis_rds) == 1L && nzchar(analysis_rds)
  if (!xor(has_config, has_rds)) {
    stop("Supply exactly one of `config` or `analysis_rds`.", call. = FALSE)
  }
  if (has_config) analyze_facs_experiment(config) else read_facs_analysis(analysis_rds)
}

facs_report_appearance <- function(params) {
  keys <- c(
    "target_axis_label", "dna_axis_label", "palette", "condition_palette",
    "condition_colors", "x_limits", "y_limits", "y_log10", "point_size",
    "density_bandwidth", "density_lower_clip", "density_upper_clip",
    "density_gamma", "show_phase_gates", "gate_color", "gate_linetype",
    "gate_linewidth", "gate_show_labels", "gate_label_size", "base_font_size",
    "error_bar", "show_points", "phase_lineplot", "layout_options"
  )
  values <- params[intersect(keys, names(params))]
  values[!vapply(values, is.null, logical(1))]
}

facs_report_save_artifacts <- function(
    analysis, bundle, analysis_path = NULL, bundle_path = NULL,
    overwrite = FALSE
) {
  if (is.character(analysis_path) && length(analysis_path) == 1L && nzchar(analysis_path)) {
    save_facs_analysis(analysis, analysis_path, overwrite = overwrite)
  }
  if (is.character(bundle_path) && length(bundle_path) == 1L && nzchar(bundle_path)) {
    save_facs_figure_bundle(bundle, bundle_path, overwrite = overwrite)
  }
  invisible(NULL)
}

facs_report_software_used <- function() {
  package_version <- function(package) {
    if (requireNamespace(package, quietly = TRUE)) {
      as.character(utils::packageVersion(package))
    } else {
      "not available"
    }
  }

  data.frame(
    component = c("facspseudocolor", "R", "Quarto", "ggplot2", "knitr"),
    version = c(
      package_version("facspseudocolor"),
      R.version.string,
      Sys.getenv("QUARTO_VERSION", unset = "not reported to R"),
      package_version("ggplot2"),
      package_version("knitr")
    ),
    stringsAsFactors = FALSE
  )
}

facs_report_compact_table <- function(data, digits = 3L) {
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  numeric <- vapply(data, is.numeric, logical(1))
  data[numeric] <- lapply(data[numeric], function(x) round(x, digits))
  knitr::kable(data, row.names = FALSE, format = "html",
               table.attr = 'class="report-table"')
}

facs_report_html_escape <- function(value) {
  value <- as.character(value)
  value <- gsub("&", "&amp;", value, fixed = TRUE)
  value <- gsub("<", "&lt;", value, fixed = TRUE)
  value <- gsub(">", "&gt;", value, fixed = TRUE)
  value <- gsub('"', "&quot;", value, fixed = TRUE)
  gsub("'", "&#39;", value, fixed = TRUE)
}

facs_report_heading_label <- function(value) {
  value <- paste(as.character(value), collapse = " ")
  value <- gsub("[[:space:]]+", " ", value)
  facs_report_html_escape(trimws(value))
}

facs_report_card_ids <- function(prefix, count) {
  if (!is.character(prefix) || length(prefix) != 1L || is.na(prefix) ||
      !grepl("^[A-Za-z][A-Za-z0-9_-]*$", prefix)) {
    stop("Card ID prefix must start with a letter and contain only letters, numbers, underscores, or hyphens.", call. = FALSE)
  }
  if (!is.numeric(count) || length(count) != 1L || is.na(count) ||
      !is.finite(count) || count < 1 || count > .Machine$integer.max ||
      count != floor(count)) {
    stop("Card ID count must be one positive integer.", call. = FALSE)
  }
  paste0(prefix, "-", seq_len(as.integer(count)))
}

facs_report_quantitation_legend_layout <- function(condition_count) {
  if (!is.numeric(condition_count) || length(condition_count) != 1L ||
      is.na(condition_count) || !is.finite(condition_count) ||
      condition_count < 1 || condition_count > .Machine$integer.max ||
      condition_count != floor(condition_count)) {
    stop("Quantitation legend condition count must be one positive integer.",
         call. = FALSE)
  }
  condition_count <- as.integer(condition_count)
  columns <- min(2L, condition_count)
  rows <- ceiling(condition_count / columns)
  list(
    condition_count = condition_count,
    columns = columns,
    rows = as.integer(rows),
    figure_height = 3
  )
}

facs_report_quantitation_row_layout <- function(card_count, min_panel_width) {
  values <- list(card_count = card_count, min_panel_width = min_panel_width)
  invalid <- vapply(values, function(value) {
    !is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value < 1 || value > .Machine$integer.max ||
      value != floor(value)
  }, logical(1))
  if (any(invalid)) {
    stop("Quantitation row card count and minimum panel width must be positive integers.",
         call. = FALSE)
  }
  columns <- min(2L, as.integer(card_count))
  list(
    columns = columns,
    min_panel_width = as.integer(min_panel_width),
    maximum_width = columns * as.integer(min_panel_width) + 16 * (columns - 1L)
  )
}

facs_report_quantitation_card_width <- function(condition_count, category_count) {
  values <- list(condition_count = condition_count, category_count = category_count)
  invalid <- vapply(values, function(value) {
    !is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value < 1 || value > .Machine$integer.max ||
      value != floor(value)
  }, logical(1))
  if (any(invalid)) {
    stop("Quantitation card condition and displayed category counts must be positive integers.",
         call. = FALSE)
  }

  condition_count <- as.double(condition_count)
  category_count <- as.double(category_count)
  calculated_width <- 140 + 36 * condition_count +
    18 * condition_count * (category_count - 1)
  if (!is.finite(calculated_width)) {
    stop("Calculated quantitation card width must be finite.", call. = FALSE)
  }
  existing_width <- min(720, max(280, calculated_width))
  scaled_width <- min(480, max(220, (2 / 3) * existing_width))
  as.integer(round(scaled_width))
}

facs_report_quantitation_plot_width <- function(card_width_px) {
  if (!is.numeric(card_width_px) || length(card_width_px) != 1L ||
      is.na(card_width_px) || !is.finite(card_width_px) || card_width_px <= 0) {
    stop("Quantitation plot CSS width must be one finite positive number.",
         call. = FALSE)
  }
  plot_width_in <- card_width_px / 96
  if (!is.finite(plot_width_in) || plot_width_in <= 0) {
    stop("Quantitation plot width in inches must be finite and positive.",
         call. = FALSE)
  }
  plot_width_in
}

facs_report_compact_typography <- function(
    title_size = 8, secondary_size = 7.5, legend_size = 7
) {
  sizes <- c(
    title_size = title_size,
    secondary_size = secondary_size,
    legend_size = legend_size
  )
  if (!is.numeric(sizes) || length(sizes) != 3L || anyNA(sizes) ||
      any(!is.finite(sizes)) || title_size < 8 || secondary_size < 7.5 ||
      legend_size < 7) {
    stop("Compact report title size must be at least 8 pt, secondary size at least 7.5 pt, and legend size at least 7 pt.",
         call. = FALSE)
  }
  ggplot2::theme(
    axis.title = ggplot2::element_text(size = title_size),
    axis.text = ggplot2::element_text(size = secondary_size),
    legend.title = ggplot2::element_text(size = title_size),
    legend.text = ggplot2::element_text(size = legend_size),
    strip.text = ggplot2::element_text(size = secondary_size)
  )
}

facs_report_compact_legend_key_theme <- function(key_size_cm = 0.28) {
  if (!is.numeric(key_size_cm) || length(key_size_cm) != 1L ||
      is.na(key_size_cm) || !is.finite(key_size_cm) || key_size_cm <= 0) {
    stop("Compact legend key size must be one positive finite number of centimetres.",
         call. = FALSE)
  }
  ggplot2::theme(
    legend.key.size = grid::unit(key_size_cm, "cm"),
    legend.spacing.x = grid::unit(0.08, "cm")
  )
}

facs_report_knit_plot <- function(
    plot, width, height, envir = parent.frame(), dpi = 96L, fig_alt = NULL
) {
  if (!inherits(plot, c("ggplot", "grob", "gTree", "gtable"))) {
    stop("`plot` must be a printable plot or grid object.", call. = FALSE)
  }
  if (!is.numeric(width) || length(width) != 1L || !is.finite(width) || width <= 0 ||
      !is.numeric(height) || length(height) != 1L || !is.finite(height) || height <= 0) {
    stop("Plot width and height must be finite positive scalars.", call. = FALSE)
  }
  if (!is.numeric(dpi) || length(dpi) != 1L || !is.finite(dpi) ||
      dpi < 1 || dpi > .Machine$integer.max || dpi != floor(dpi)) {
    stop("Plot DPI must be one positive integer.", call. = FALSE)
  }
  if (!is.null(fig_alt) &&
      (!is.character(fig_alt) || length(fig_alt) != 1L || is.na(fig_alt) ||
       !nzchar(trimws(fig_alt)))) {
    stop("`fig_alt` must be null or one nonempty string.", call. = FALSE)
  }
  child_environment <- new.env(parent = envir)
  child_environment$report_plot <- plot
  child <- c(
    "```{r}",
    paste0("#| fig-width: ", format(width, scientific = FALSE, trim = TRUE)),
    paste0("#| fig-height: ", format(height, scientific = FALSE, trim = TRUE)),
    paste0("#| fig-dpi: ", as.integer(dpi)),
    "#| echo: false",
    if (!is.null(fig_alt)) paste0(
      "#| fig-alt: ",
      jsonlite::toJSON(as.character(fig_alt), auto_unbox = TRUE)
    ),
    "if (inherits(report_plot, 'ggplot')) {",
    "  print(report_plot)",
    "} else {",
    "  grid::grid.newpage()",
    "  grid::grid.draw(report_plot)",
    "}",
    "```"
  )
  knitr::knit_child(text = child, envir = child_environment, quiet = TRUE)
}

facs_report_edu_overview <- function(analysis) {
  facspseudocolor:::validate_analysis_object(analysis)
  if (!identical(analysis$config$plot_type, "edu")) {
    stop("The EdU report overview requires an EdU analysis.", call. = FALSE)
  }
  manifest <- analysis$sample_manifest
  input <- analysis$input_report
  reference <- unique(as.character(manifest$reference_condition))
  reference <- reference[!is.na(reference) & nzchar(reference)]
  target_label <- facspseudocolor::facs_report_edu_target_label(analysis)
  population_names <- if (identical(analysis$config$gating$mode,
                                    "model_experimental")) {
    c(complete = "Model-derived Single Cells", g1 = "Model-derived G1",
      edu_positive = paste("Model-derived", target_label, "Positive"))
  } else {
    c(complete = "FlowJo Single Cells", g1 = "FlowJo G1",
      edu_positive = paste("FlowJo", target_label, "Positive"))
  }
  populations <- unique(as.character(input$population[input$exists %in% TRUE]))
  if (!length(populations)) {
    populations <- intersect(names(analysis$config$suffixes), names(population_names))
  }
  populations <- unname(population_names[populations])
  populations <- populations[!is.na(populations)]
  list(
    summary = data.frame(
      item = c(
        "Samples / acquisitions", "Biological replicates", "Conditions",
        "DNA channel", paste(target_label, "channel"), "Population sources",
        "DNA normalization", paste(target_label, "correction"), "Intensity reference",
        "Configuration"
      ),
      value = c(
        nrow(manifest), length(unique(manifest$replicate_index)),
        paste(unique(manifest$condition[order(manifest$condition_index)]), collapse = ", "),
        analysis$config$dna_channel, analysis$config$target_channel,
        paste(populations, collapse = ", "),
        paste0("G1/2N mapped to ", analysis$config$dna_2n_value,
               " using the configured ", analysis$config$g1_anchor, " anchor"),
        paste0(
          "Independent ", target_label, "-negative baseline slope fitted for each acquisition ",
          "and combined with that acquisition's G1-derived anchor using the configured anchor method"
        ),
        if (length(reference)) paste(reference, collapse = ", ") else "Not configured",
        if (is.character(analysis$provenance$config_path) &&
            length(analysis$provenance$config_path) == 1L &&
            !is.na(analysis$provenance$config_path) &&
            nzchar(analysis$provenance$config_path)) {
          basename(analysis$provenance$config_path)
        } else "Completed analysis artifact"
      ), stringsAsFactors = FALSE
    ),
    samples = manifest[, c(
      "replicate", "technical_replicate", "condition", "prefix", "is_reference"
    ), drop = FALSE],
    inputs = input[, intersect(c(
      "replicate", "condition", "prefix", "population", "required", "exists",
      "event_n", "finite_event_n", "nonfinite_event_n", "path"
    ), names(input)), drop = FALSE]
  )
}

facs_report_edu_qc <- function(analysis, report) {
  facspseudocolor:::validate_analysis_object(analysis)
  findings <- list()
  add <- function(level, area, sample, message) {
    findings[[length(findings) + 1L]] <<- data.frame(
      level = level, area = area, sample = sample, finding = message,
      stringsAsFactors = FALSE
    )
  }
  warnings <- as.character(analysis$warnings)
  warnings <- warnings[!is.na(warnings) & nzchar(warnings)]
  for (warning in warnings) add("REVIEW", "Analysis", "Experiment", warning)
  input <- analysis$input_report
  if (is.data.frame(input)) {
    for (i in seq_len(nrow(input))) {
      label <- paste(input$prefix[[i]], input$population[[i]], sep = " / ")
      if (isTRUE(input$required[[i]]) && !isTRUE(input$exists[[i]])) {
        add("BLOCKED", "Input", label, "Required input is missing.")
      }
      if ("nonfinite_event_n" %in% names(input) &&
          is.finite(input$nonfinite_event_n[[i]]) && input$nonfinite_event_n[[i]] > 0) {
        add("REVIEW", "Input", label,
            paste(input$nonfinite_event_n[[i]], "non-finite event rows recorded."))
      }
    }
  }
  panel_qc <- report$panel_qc
  for (i in which(panel_qc$panel_status != "available")) {
    add("BLOCKED", "DNA-versus-EdU panel", panel_qc$prefix[[i]],
        paste("Panel suppressed:", panel_qc$panel_reason_code[[i]]))
  }
  positivity_tables <- list(
    "2N-4N positivity" = report$positivity$overall_2to4n_biological_replicate,
    "regional positivity" = report$positivity$regional_biological_replicate
  )
  for (metric in names(positivity_tables)) {
    values <- positivity_tables[[metric]]
    if ("metric_status" %in% names(values)) {
      for (i in which(values$metric_status != "ok")) {
        add("BLOCKED", "Quantitation", as.character(values$condition[[i]]),
            paste(metric, "status:", values$metric_status[[i]]))
      }
    }
  }
  intensity_tables <- list(
    "all-positive intensity" = report$intensity$overall_biological_replicate,
    "regional intensity" = report$intensity$regional_biological_replicate
  )
  for (metric in names(intensity_tables)) {
    values <- intensity_tables[[metric]]
    if ("reference_normalization_status" %in% names(values)) {
      for (i in which(values$reference_normalization_status != "available")) {
        add("BLOCKED", "Reference normalization",
            as.character(values$condition[[i]]),
            paste(metric, "status:", values$reference_normalization_status[[i]]))
      }
    }
  }
  if (!length(findings)) {
    add("CLEAR", "Report checks", "Experiment",
        "No retained warning, missing input, non-finite input, suppressed panel, zero denominator, or unavailable reference normalization was found.")
  }
  do.call(rbind, findings)
}

facs_report_edu_qc_sections <- function(analysis, report) {
  all <- facs_report_edu_qc(analysis, report)
  known_alias_deprecation <- paste0(
    "EdU aliases `phase_percentages`, `phase_medians`, and `whole_medians` ",
    "are deprecated for one major-release compatibility window; use the ",
    "canonical `edu_*` tables. Alias meanings are unchanged."
  )
  software_notice <- all$finding == known_alias_deprecation
  normalization <- all$area %in% c("Reference normalization", "Quantitation")
  data_input <- !(software_notice | normalization)
  clear_row <- function(area, message) data.frame(
    level = "CLEAR", area = area, sample = "Experiment", finding = message,
    stringsAsFactors = FALSE
  )
  subset_or_clear <- function(keep, area, message) {
    value <- all[keep & all$level != "CLEAR", , drop = FALSE]
    if (nrow(value)) value else clear_row(area, message)
  }
  list(
    data_input = subset_or_clear(
      data_input, "Data and inputs",
      "No retained data/input warning, missing required input, non-finite input, or suppressed panel was found."
    ),
    normalization = subset_or_clear(
      normalization, "Normalization and quantitation",
      "No zero denominator or unavailable configured-reference normalization was found."
    ),
    visual_review = data.frame(
      level = "REQUIRED", area = "Human visual review", sample = "Every acquisition",
      finding = "Inspect DNA-versus-EdU distributions, density shape, clipping, and replicate agreement; no new automated distribution threshold is applied.",
      stringsAsFactors = FALSE
    ),
    software = subset_or_clear(
      software_notice, "Software notices", "No retained software notice was found."
    )
  )
}

facs_report_edu_balance <- function(analysis) {
  manifest <- analysis$sample_manifest
  conditions <- unique(as.character(manifest$condition[order(manifest$condition_index)]))
  replicates <- unique(manifest$replicate_index[order(manifest$replicate_index)])
  rows <- lapply(replicates, function(replicate_index) {
    replicate_rows <- manifest[manifest$replicate_index == replicate_index, , drop = FALSE]
    labels <- unique(as.character(replicate_rows$replicate))
    if (length(labels) != 1L) stop("Each replicate_index must map to one replicate label.", call. = FALSE)
    data.frame(
      `Biological replicate` = labels[[1L]],
      Condition = conditions,
      `Technical acquisitions` = vapply(
        conditions, function(condition) sum(replicate_rows$condition == condition), integer(1)
      ), check.names = FALSE, stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$Status <- ifelse(out$`Technical acquisitions` > 0L, "Present", "MISSING")
  rownames(out) <- NULL
  out
}

facs_report_edu_all_conditions <- function(analysis, report, max_condition_columns = 3L) {
  if (!is.numeric(max_condition_columns) || length(max_condition_columns) != 1L ||
      is.na(max_condition_columns) || !is.finite(max_condition_columns) ||
      max_condition_columns < 1 || max_condition_columns > .Machine$integer.max ||
      max_condition_columns != floor(max_condition_columns)) {
    stop("`max_condition_columns` must be one positive integer.", call. = FALSE)
  }
  max_condition_columns <- as.integer(max_condition_columns)
  manifest <- analysis$sample_manifest
  required <- c("prefix", "condition", "condition_index", "replicate", "replicate_index", "technical_replicate")
  missing_fields <- setdiff(required, names(manifest))
  if (length(missing_fields)) stop(paste("All-condition overview requires:", paste(missing_fields, collapse = ", ")), call. = FALSE)
  prefixes <- as.character(manifest$prefix)
  if (anyNA(prefixes) || any(!nzchar(prefixes)) || anyDuplicated(prefixes) ||
      is.null(names(report$panels)) || !identical(names(report$panels), prefixes)) {
    stop("All-condition panels must be named in exact manifest-prefix order.", call. = FALSE)
  }
  identity_fields <- c("condition", "replicate", "technical_replicate")
  invalid_identity <- vapply(identity_fields, function(field) {
    values <- as.character(manifest[[field]])
    anyNA(values) || any(!nzchar(values))
  }, logical(1))
  if (any(invalid_identity) || anyNA(manifest$condition_index) ||
      anyNA(manifest$replicate_index)) {
    stop("All-condition overview requires nonmissing condition, biological-replicate, technical-acquisition, and ordering identities.", call. = FALSE)
  }
  labels_per_replicate <- vapply(
    split(as.character(manifest$replicate), manifest$replicate_index),
    function(value) length(unique(value)), integer(1)
  )
  if (any(labels_per_replicate != 1L)) {
    stop("Each replicate_index must map to exactly one biological-replicate label.", call. = FALSE)
  }
  conditions <- unique(as.character(manifest$condition[order(manifest$condition_index)]))
  replicates <- unique(manifest$replicate_index[order(manifest$replicate_index)])
  placeholder <- function(label) ggplot2::ggplot() +
    ggplot2::annotate("text", .5, .56, label = "MISSING", fontface = "bold", colour = "#9b2c2c", size = 5) +
    ggplot2::annotate("text", .5, .43, label = label, colour = "#5f3030", size = 3) +
    ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) + ggplot2::theme_void()
  overview_panel <- function(panel) {
    panel + ggplot2::labs(x = "DNA", y = facspseudocolor::facs_report_edu_target_label(analysis)) +
      ggplot2::theme(legend.position = "none")
  }
  stack <- function(rows, label) {
    if (!nrow(rows)) return(placeholder(label))
    panels <- lapply(report$panels[as.character(rows$prefix)], overview_panel)
    cowplot::plot_grid(plotlist = panels, ncol = 1L, align = "v")
  }
  panel_qc <- report$panel_qc
  if (!is.data.frame(panel_qc) ||
      !all(c("prefix", "panel_status") %in% names(panel_qc)) ||
      !identical(as.character(panel_qc$prefix), prefixes)) {
    stop("All-condition overview requires panel QC in exact manifest-prefix order.", call. = FALSE)
  }
  available_index <- which(panel_qc$panel_status == "available")
  legend_source_prefix <- if (length(available_index)) prefixes[[available_index[[1L]]]] else NA_character_
  shared_legend <- if (length(available_index)) {
    cowplot::get_legend(
      report$panels[[available_index[[1L]]]] +
        ggplot2::labs(colour = NULL) +
        ggplot2::theme(legend.position = "bottom")
    )
  } else NULL
  legend_available <- !is.null(shared_legend)
  condition_blocks <- split(conditions, ceiling(seq_along(conditions) / max_condition_columns))
  blocks <- lapply(seq_along(condition_blocks), function(block_index) {
    block_conditions <- unname(condition_blocks[[block_index]])
    row_weights <- vapply(replicates, function(index) {
      rows <- manifest[manifest$replicate_index == index, , drop = FALSE]
      as.integer(max(1L, vapply(block_conditions, function(condition) {
        sum(rows$condition == condition)
      }, integer(1))))
    }, integer(1))
    row_parts <- unlist(lapply(replicates, function(index) {
      rows <- manifest[manifest$replicate_index == index, , drop = FALSE]
      label <- unique(as.character(rows$replicate))[[1L]]
      cells <- lapply(block_conditions, function(condition) stack(
        rows[rows$condition == condition, , drop = FALSE], paste(label, "/", condition)
      ))
      cell_row <- cowplot::plot_grid(plotlist = cells, nrow = 1L, align = "hv")
      row_header <- cowplot::ggdraw() + cowplot::draw_label(
        label, x = 0, hjust = 0, fontface = "bold", size = 9
      )
      list(row_header, cell_row)
    }), recursive = FALSE)
    header <- cowplot::plot_grid(
      plotlist = lapply(block_conditions, function(condition) {
        cowplot::ggdraw() + cowplot::draw_label(condition, fontface = "bold", size = 9)
      }), nrow = 1L
    )
    body <- cowplot::plot_grid(
      plotlist = c(list(header), row_parts), ncol = 1L,
      rel_heights = c(.08, as.vector(rbind(rep(.12, length(replicates)), row_weights)))
    )
    block_plot <- if (legend_available) {
      cowplot::plot_grid(body, shared_legend, ncol = 1L, rel_heights = c(1, .08))
    } else body
    list(
      block_index = block_index,
      conditions = block_conditions,
      plot = block_plot,
      figure_width = max(8, 3.1 * length(block_conditions)),
      figure_height = 1.2 + 3.8 * sum(row_weights),
      shared_legend_count = if (legend_available) 1L else 0L,
      shared_legend_title = NULL,
      legend_source_prefix = legend_source_prefix,
      panel_legends_visible = FALSE,
      y_axis_title = facspseudocolor::facs_report_edu_target_label(analysis),
      row_header_position = "above_full_width",
      replicate_header_rel_height = .12,
      replicate_header_height_fixed = TRUE,
      replicate_header_rel_heights = rep(.12, length(replicates)),
      replicate_body_rel_heights = unname(row_weights),
      replicate_labels = vapply(replicates, function(index) {
        unique(as.character(manifest$replicate[manifest$replicate_index == index]))[[1L]]
      }, character(1))
    )
  })
  list(
    blocks = blocks,
    condition_order = conditions,
    max_condition_columns = max_condition_columns,
    density_legend_status = if (length(available_index)) {
      if (legend_available) "available" else "unavailable_from_available_panel"
    } else "unavailable_all_panels_suppressed",
    all_panels_suppressed = !length(available_index),
    balance = facs_report_edu_balance(analysis)
  )
}

facs_report_edu_apex_comparison <- function(analysis, report) {
  facspseudocolor:::validate_analysis_object(analysis)
  if (!inherits(report, "edu_pseudocolor_output_contract")) {
    stop("A completed EdU pseudocolor output contract is required for apex comparison.",
         call. = FALSE)
  }
  manifest <- analysis$sample_manifest
  required_manifest <- c("prefix", "condition", "model_group", "is_reference",
                         "reference_condition")
  if (!all(required_manifest %in% names(manifest))) {
    stop("Apex comparison requires explicit model-group and reference identities.",
         call. = FALSE)
  }
  character_identity_fields <- c("prefix", "condition", "model_group",
                                 "reference_condition")
  invalid_identity <- vapply(character_identity_fields, function(field) {
    value <- manifest[[field]]
    !is.character(value) || length(value) != nrow(manifest) ||
      anyNA(value) || any(!nzchar(value))
  }, logical(1))
  if (any(invalid_identity) || !is.logical(manifest$is_reference) ||
      length(manifest$is_reference) != nrow(manifest) ||
      anyNA(manifest$is_reference)) {
    stop("Apex comparison requires nonempty character identities and strict nonmissing logical reference flags.",
         call. = FALSE)
  }
  prefixes <- as.character(manifest$prefix)
  if (anyDuplicated(prefixes) || is.null(names(analysis$normalized_data)) ||
      !identical(names(analysis$normalized_data), prefixes) ||
      is.null(names(report$panels)) || !identical(names(report$panels), prefixes)) {
    stop("Apex comparison inputs must follow exact unique manifest-prefix order.",
         call. = FALSE)
  }
  offset_qc <- report$display_offset_qc
  if (!is.data.frame(offset_qc) ||
      !all(c("prefix", "display_status", "display_offset") %in% names(offset_qc)) ||
      !identical(as.character(offset_qc$prefix), prefixes)) {
    stop("Apex comparison requires display offsets in exact manifest-prefix order.",
         call. = FALSE)
  }
  apex_range <- as.numeric(unlist(analysis$config$edu_apex_x_range))
  density_adjust <- analysis$config$edu_apex_density_adjust
  if (length(apex_range) != 2L || any(!is.finite(apex_range)) ||
      apex_range[[1L]] >= apex_range[[2L]] ||
      !is.numeric(density_adjust) || length(density_adjust) != 1L ||
      !is.finite(density_adjust) || density_adjust <= 0) {
    stop("Apex comparison requires a valid configured apex range and density adjustment.",
         call. = FALSE)
  }
  apex_rows <- lapply(seq_along(prefixes), function(i) {
    prefix <- prefixes[[i]]
    if (!identical(offset_qc$display_status[[i]], "available") ||
        !is.finite(offset_qc$display_offset[[i]])) {
      stop("Apex comparison requires an available finite display offset for ",
           prefix, ".", call. = FALSE)
    }
    positive <- analysis$normalized_data[[i]]$edu_positive
    if (!is.data.frame(positive) ||
        !all(c("dna_norm", "target_bgsub") %in% names(positive))) {
      stop("Apex comparison requires retained EdU-positive normalized events for ",
           prefix, ".", call. = FALSE)
    }
    apex <- suppressWarnings(facspseudocolor:::calculate_edu_apex_density(
      positive, y_column = "target_bgsub", mid_x_range = apex_range,
      density_adjust = density_adjust, minimum_events = 10L,
      condition_label = prefix
    ))
    if (!is.finite(apex$y)) {
      stop("Apex comparison requires at least ten finite positive-signal EdU-positive events in the configured apex window for ",
           prefix, "; observed ", apex$n, ".", call. = FALSE)
    }
    data.frame(prefix = prefix, model_group = as.character(manifest$model_group[[i]]),
               condition = as.character(manifest$condition[[i]]),
               is_reference = manifest$is_reference[[i]],
               display_offset = offset_qc$display_offset[[i]],
               apex_background_subtracted_signal = apex$y,
               apex_event_n = apex$n,
               stringsAsFactors = FALSE)
  })
  apex_table <- do.call(rbind, apex_rows)
  reference_by_group <- lapply(unique(apex_table$model_group), function(group_id) {
    group_index <- which(apex_table$model_group == group_id)
    reference_index <- group_index[apex_table$is_reference[group_index]]
    declared <- unique(as.character(manifest$reference_condition[
      as.character(manifest$model_group) == group_id
    ]))
    if (length(reference_index) != 1L || length(declared) != 1L ||
        is.na(declared) || !nzchar(declared) ||
        !identical(apex_table$condition[[reference_index]], declared[[1L]])) {
      stop("Apex comparison requires exactly one correctly identified reference in model group ",
           group_id, ".", call. = FALSE)
    }
    c(group_id = group_id, reference_prefix = apex_table$prefix[[reference_index]],
      reference_apex = apex_table$apex_background_subtracted_signal[[reference_index]])
  })
  reference_table <- do.call(rbind, lapply(reference_by_group, function(value) {
    data.frame(model_group = value[["group_id"]],
               reference_prefix = value[["reference_prefix"]],
               reference_apex_background_subtracted_signal =
                 as.numeric(value[["reference_apex"]]),
               stringsAsFactors = FALSE)
  }))
  apex_table$reference_prefix <- reference_table$reference_prefix[
    match(apex_table$model_group, reference_table$model_group)
  ]
  apex_table$reference_apex_background_subtracted_signal <-
    reference_table$reference_apex_background_subtracted_signal[
      match(apex_table$model_group, reference_table$model_group)
    ]
  if (anyNA(apex_table$reference_prefix) ||
      any(!is.finite(apex_table$reference_apex_background_subtracted_signal))) {
    stop("Apex comparison reference coverage is incomplete.", call. = FALSE)
  }
  apex_table$current_line_display_signal <-
    apex_table$apex_background_subtracted_signal + apex_table$display_offset
  apex_table$reference_line_display_signal <-
    apex_table$reference_apex_background_subtracted_signal +
      apex_table$display_offset
  apex_report <- report
  for (i in seq_along(prefixes)) {
    current <- apex_table[i, , drop = FALSE]
    panel <- apex_report$panels[[prefixes[[i]]]] + ggplot2::geom_hline(
      yintercept = current$reference_line_display_signal[[1L]],
      colour = "#2166AC", linetype = "dashed", linewidth = 0.7
    )
    if (!isTRUE(current$is_reference[[1L]])) {
      panel <- panel + ggplot2::geom_hline(
        yintercept = current$current_line_display_signal[[1L]],
        colour = "#B2182B", linetype = "solid", linewidth = 0.7
      )
    }
    apex_report$panels[[prefixes[[i]]]] <- panel
  }
  group_apex <- split(seq_len(nrow(apex_table)), apex_table$model_group)
  for (indices in group_apex) {
    panel_indices <- match(apex_table$prefix[indices], report$panel_qc$prefix)
    values <- c(apex_table$current_line_display_signal[indices],
                apex_table$reference_line_display_signal[indices])
    lower <- min(c(report$panel_qc$y_limit_lower[panel_indices], values))
    upper <- max(c(report$panel_qc$y_limit_upper[panel_indices], values))
    if (isTRUE(analysis$config$y_log10)) {
      if (lower <= 0) stop("Apex comparison log display limits must be positive.",
                           call. = FALSE)
      padding <- diff(log10(c(lower, upper))) * 0.03
      lower <- 10^(log10(lower) - padding)
      upper <- 10^(log10(upper) + padding)
    } else {
      padding <- (upper - lower) * 0.03
      lower <- lower - padding
      upper <- upper + padding
    }
    apex_report$panel_qc$y_limit_lower[panel_indices] <- lower
    apex_report$panel_qc$y_limit_upper[panel_indices] <- upper
  }
  list(report = apex_report, apex = apex_table,
       reference_colour = "#2166AC", reference_linetype = "dashed",
       sample_colour = "#B2182B", sample_linetype = "solid",
       method = "mode_of_log10_background_subtracted_signal_within_configured_edu_positive_apex_window_then_destination_display_offset")
}

facs_report_edu_phase_gate_comparison <- function(analysis, report) {
  facspseudocolor:::validate_analysis_object(analysis)
  if (!inherits(report, "edu_pseudocolor_output_contract")) {
    stop("A completed EdU pseudocolor output contract is required for phase-gate overlays.",
         call. = FALSE)
  }
  prefixes <- as.character(analysis$sample_manifest$prefix)
  if (!identical(names(report$panels), prefixes) ||
      !identical(as.character(report$panel_qc$prefix), prefixes) ||
      !identical(as.character(report$display_offset_qc$prefix), prefixes)) {
    stop("Phase-gate overlays require exact manifest-prefix report order.",
         call. = FALSE)
  }
  phase_report <- report
  gate_tables <- vector("list", length(prefixes))
  names(gate_tables) <- prefixes
  for (i in seq_along(prefixes)) {
    limits <- c(report$panel_qc$y_limit_lower[[i]],
                report$panel_qc$y_limit_upper[[i]])
    offset <- report$display_offset_qc$display_offset[[i]]
    if (any(!is.finite(limits)) || limits[[1L]] >= limits[[2L]] ||
        !is.finite(offset)) {
      stop("Phase-gate overlays require available display limits and offsets for ",
           prefixes[[i]], ".", call. = FALSE)
    }
    gates <- facspseudocolor:::analysis_display_gate_rectangles(
      analysis, limits, offset
    )
    polygons <- facspseudocolor:::make_computed_edu_gate_polygons(
      analysis$normalized_data[[i]], gates, limits, offset
    )
    labels <- do.call(rbind, lapply(split(polygons, polygons$gate), function(x) {
      positive_y <- x$y[is.finite(x$y) & x$y > 0]
      if (!length(positive_y)) {
        stop("Phase-gate label has no positive display coordinate.",
             call. = FALSE)
      }
      data.frame(
        gate = x$gate[[1L]], x = mean(range(x$x)),
        y = exp(mean(log(range(positive_y)))), stringsAsFactors = FALSE
      )
    }))
    phase_report$panels[[i]] <- facspseudocolor:::add_phase_gate_polygons_to_plot(
      phase_report$panels[[i]], polygons, color = "#1A1A1A",
      linetype = "dashed", linewidth = 0.55
    ) + ggplot2::geom_text(
      data = labels, ggplot2::aes(x = .data$x, y = .data$y, label = .data$gate),
      inherit.aes = FALSE, colour = "#1A1A1A", size = 2.2,
      fontface = "bold"
    )
    gate_tables[[i]] <- gates
  }
  list(report = phase_report, gates = gate_tables,
       line_colour = "#1A1A1A", line_type = "dashed")
}

facs_report_embedded_plot_pdf_link <- function(
    plot, filename, width = 3, height = 3, label = "Download PDF",
    maximum_bytes = 8 * 1024^2
) {
  if (!inherits(plot, "ggplot")) {
    stop("Embedded PDF download requires one ggplot object.", call. = FALSE)
  }
  if (!is.character(filename) || length(filename) != 1L || is.na(filename) ||
      !nzchar(filename) || !is.character(label) || length(label) != 1L ||
      is.na(label) || !nzchar(label) ||
      !is.numeric(width) || length(width) != 1L ||
      !is.finite(width) || width <= 0 || !is.numeric(height) ||
      length(height) != 1L || !is.finite(height) || height <= 0 ||
      !is.numeric(maximum_bytes) || length(maximum_bytes) != 1L ||
      !is.finite(maximum_bytes) || maximum_bytes <= 0) {
    stop("Embedded PDF download requires a filename and positive dimensions.",
         call. = FALSE)
  }
  safe_filename <- gsub("[^A-Za-z0-9._-]", "_", basename(filename))
  if (!grepl("[.]pdf$", safe_filename, ignore.case = TRUE)) {
    safe_filename <- paste0(safe_filename, ".pdf")
  }
  pdf_path <- tempfile("facs-plot-download-", fileext = ".pdf")
  device_open <- FALSE
  on.exit({
    if (device_open) grDevices::dev.off()
    unlink(pdf_path)
  }, add = TRUE)
  grDevices::pdf(pdf_path, width = width, height = height,
                 onefile = TRUE, useDingbats = FALSE)
  device_open <- TRUE
  print(plot)
  grDevices::dev.off()
  device_open <- FALSE
  size <- file.info(pdf_path)$size
  if (!is.finite(size) || size <= 0) {
    stop("Embedded plot PDF generation produced no readable content.",
         call. = FALSE)
  }
  if (size > maximum_bytes) {
    stop("Embedded plot PDF exceeds the configured per-plot safety limit.",
         call. = FALSE)
  }
  bytes <- readBin(pdf_path, what = "raw", n = size)
  encoded <- openssl::base64_encode(bytes)
  safe_label <- gsub("&", "&amp;", label, fixed = TRUE)
  safe_label <- gsub("<", "&lt;", safe_label, fixed = TRUE)
  safe_label <- gsub(">", "&gt;", safe_label, fixed = TRUE)
  result <- paste0(
    '<a class="plot-pdf-download" download="', safe_filename,
    '" href="data:application/pdf;base64,', encoded, '">',
    safe_label, '</a>'
  )
  attr(result, "pdf_bytes") <- as.numeric(size)
  result
}

facs_report_edu_all_panels_canvas <- function(report, max_columns = 4L,
                                              panel_width = 3,
                                              panel_height = 3) {
  if (!is.list(report) || !is.list(report$panels) || !length(report$panels) ||
      is.null(names(report$panels)) || any(!nzchar(names(report$panels))) ||
      any(!vapply(report$panels, inherits, logical(1), what = "ggplot")) ||
      !is.numeric(max_columns) || length(max_columns) != 1L ||
      is.na(max_columns) || !is.finite(max_columns) || max_columns < 1 ||
      max_columns != floor(max_columns) ||
      !is.numeric(panel_width) || length(panel_width) != 1L ||
      !is.finite(panel_width) || panel_width <= 0 ||
      !is.numeric(panel_height) || length(panel_height) != 1L ||
      !is.finite(panel_height) || panel_height <= 0) {
    stop("All-panel PDF canvas requires named ggplot panels and positive dimensions.",
         call. = FALSE)
  }
  columns <- min(as.integer(max_columns), length(report$panels))
  rows <- ceiling(length(report$panels) / columns)
  list(
    plot = cowplot::plot_grid(plotlist = unname(report$panels), ncol = columns,
                              align = "hv"),
    width = columns * panel_width,
    height = rows * panel_height,
    panel_count = length(report$panels),
    columns = columns,
    rows = rows
  )
}

facs_report_edu_gating_cards <- function(analysis, all_events = NULL,
                                         dna_height_channel = NULL,
                                         display_offsets = NULL,
                                         max_points = 3000L) {
  if (!inherits(analysis, "facs_analysis") ||
      !identical(analysis$config$plot_type, "edu") ||
      length(analysis$normalized_data) != nrow(analysis$sample_manifest) ||
      !is.numeric(max_points) || length(max_points) != 1L ||
      is.na(max_points) || !is.finite(max_points) || max_points < 1 ||
      max_points != floor(max_points)) {
    stop("EdU gating cards require one valid EdU analysis and a positive point limit.",
         call. = FALSE)
  }
  max_points <- as.integer(max_points)
  prefixes <- as.character(analysis$sample_manifest$prefix)
  if (!is.numeric(display_offsets) ||
      !identical(names(display_offsets), prefixes) ||
      any(!is.finite(display_offsets)) || any(display_offsets < 0)) {
    stop("Background comparison cards require one nonnegative provenance-bound report display offset in exact sample order.",
         call. = FALSE)
  }
  gating_acquisitions <- analysis$config$gating$acquisitions
  model_derived <- identical(analysis$config$gating$mode, "model_experimental")
  all_event_display <- !is.null(all_events)
  if (model_derived &&
      (!is.list(gating_acquisitions) || !length(gating_acquisitions))) {
    stop("DNA-W gating displays require explicit per-acquisition DNA-height channel roles.",
         call. = FALSE)
  }
  if (all_event_display &&
      (!is.list(all_events) ||
       !identical(names(all_events), as.character(analysis$sample_manifest$prefix)))) {
    stop("DNA-A/DNA-H gating displays require explicit all-events exports in exact sample order.",
         call. = FALSE)
  }
  if (!model_derived && all_event_display &&
      (!is.character(dna_height_channel) || length(dna_height_channel) != 1L ||
       is.na(dna_height_channel) || !nzchar(dna_height_channel))) {
    stop("FlowJo DNA-A/DNA-H gating displays require an explicit DNA-height channel.",
         call. = FALSE)
  }
  retain_points <- function(data) {
    if (nrow(data) <= max_points) return(data)
    data[unique(as.integer(round(seq(1, nrow(data), length.out = max_points)))),
         , drop = FALSE]
  }
  robust_axis_limits <- function(values) {
    values <- values[is.finite(values)]
    if (length(values) < 2L) {
      stop("DNA-A/DNA-H display limits require at least two finite values.",
           call. = FALSE)
    }
    limits <- as.numeric(stats::quantile(
      values, c(0.001, 0.999), names = FALSE, type = 7
    ))
    span <- diff(limits)
    if (!all(is.finite(limits)) || span <= 0) {
      stop("DNA-A/DNA-H display limits are degenerate.", call. = FALSE)
    }
    limits + c(-1, 1) * span * 0.03
  }
  overlay_plot <- function(parent, child, title, child_colour, identity,
                           dna_channel, target_channel) {
    identity_column <- if (model_derived) "event_identity" else "event_index"
    if (identity_column %in% names(parent) && identity_column %in% names(child)) {
      if (anyDuplicated(parent[[identity_column]]) ||
          anyDuplicated(child[[identity_column]])) {
        stop("Exclusive gating overlays require unique event identities.",
             call. = FALSE)
      }
      parent <- parent[
        !parent[[identity_column]] %in% child[[identity_column]], , drop = FALSE
      ]
    }
    parent <- parent[is.finite(parent[[target_channel]]) &
                       parent[[target_channel]] > 0, , drop = FALSE]
    child <- child[is.finite(child[[target_channel]]) &
                     child[[target_channel]] > 0, , drop = FALSE]
    ggplot2::ggplot() +
      ggplot2::geom_point(
        data = retain_points(parent),
        ggplot2::aes(x = .data[[dna_channel]], y = .data[[target_channel]]),
        colour = "#4D4D4D", shape = 1, size = 0.34, stroke = 0.25,
        alpha = 0.65
      ) +
      ggplot2::geom_point(
        data = retain_points(child),
        ggplot2::aes(x = .data[[dna_channel]], y = .data[[target_channel]]),
        colour = child_colour, shape = 16, size = 0.42, stroke = 0,
        alpha = 0.78
      ) +
      ggplot2::labs(title = title, subtitle = identity, x = "DNA", y = facspseudocolor::facs_report_edu_target_label(analysis)) +
      ggplot2::scale_y_log10() +
      ggplot2::theme_classic(base_size = 8) +
      ggplot2::theme(aspect.ratio = 1,
                     plot.title = ggplot2::element_text(face = "bold"))
  }
  manifest <- analysis$sample_manifest
  lapply(seq_len(nrow(manifest)), function(i) {
    sample <- analysis$normalized_data[[i]]
    prefix <- as.character(manifest$prefix[[i]])
    display_offset <- unname(display_offsets[[prefix]])
    condition <- as.character(manifest$condition[[i]])
    identity <- paste0(condition, " | Sample: ", prefix)
    dna_area_source_channel <- NA_character_
    sample_dna_height_channel <- NA_character_
    if (model_derived) {
      acquisition_prefixes <- vapply(
        gating_acquisitions, function(item) as.character(item$prefix %||% ""),
        character(1)
      )
      acquisition_index <- which(acquisition_prefixes == prefix)
      if (length(acquisition_index) != 1L) {
        stop("DNA channel roles are not uniquely configured for ", prefix,
             ".", call. = FALSE)
      }
      roles <- gating_acquisitions[[acquisition_index]]$channel_roles
      dna_area_source_channel <- roles$dna_area
      sample_dna_height_channel <- roles$dna_height
      if (!is.character(dna_area_source_channel) ||
          length(dna_area_source_channel) != 1L ||
          is.na(dna_area_source_channel) || !nzchar(dna_area_source_channel) ||
          !is.character(sample_dna_height_channel) ||
          length(sample_dna_height_channel) != 1L ||
          is.na(sample_dna_height_channel) || !nzchar(sample_dna_height_channel) ||
          !identical(analysis$config$dna_channel, "DNA content")) {
        stop("Model-derived DNA-A/DNA-W display roles or exported DNA-area mapping are unavailable for ",
             prefix, ".", call. = FALSE)
      }
    } else if (all_event_display) {
      dna_area_source_channel <- analysis$config$dna_channel
      sample_dna_height_channel <- dna_height_channel
    }
    if (!is.data.frame(sample$data) || !is.data.frame(sample$g1) ||
        !is.data.frame(sample$edu_positive)) {
      stop("Validated Single Cells, G1, and EdU-positive inputs are required for gating card ",
           prefix, ".", call. = FALSE)
    }
    g1 <- sample$g1
    single <- sample$data
    positive <- sample$edu_positive
    required_single <- c(analysis$config$dna_channel,
                         analysis$config$target_channel, "baseline",
                         "target_bgsub")
    if (all_event_display) required_single <- c(required_single, sample_dna_height_channel)
    required_child <- c(analysis$config$dna_channel,
                        analysis$config$target_channel)
    if (any(!required_single %in% names(single)) ||
        any(!required_child %in% names(g1)) ||
        any(!required_child %in% names(positive)) ||
        !nrow(single) || !nrow(g1) || !nrow(positive)) {
      problems <- character()
      missing_single <- setdiff(required_single, names(single))
      missing_g1 <- setdiff(required_child, names(g1))
      missing_positive <- setdiff(required_child, names(positive))
      if (length(missing_single)) problems <- c(problems, paste0(
        "Single Cells missing column(s): ", paste(missing_single, collapse = ", ")
      ))
      if (length(missing_g1)) problems <- c(problems, paste0(
        "G1 missing column(s): ", paste(missing_g1, collapse = ", ")
      ))
      if (length(missing_positive)) problems <- c(problems, paste0(
        "EdU-positive missing column(s): ", paste(missing_positive, collapse = ", ")
      ))
      if (!nrow(single)) problems <- c(problems, "Single Cells has zero events")
      if (!nrow(g1)) problems <- c(problems, "G1 has zero events")
      if (!nrow(positive)) problems <- c(problems, "EdU-positive has zero events")
      stop("Gating card inputs lack required channels or events for ", prefix,
           ": ", paste(problems, collapse = "; "), ".", call. = FALSE)
    }
    if (all_event_display) {
      all_data <- all_events[[prefix]]
      identity_column <- if (model_derived) "event_identity" else "event_index"
      required_all_events <- c(analysis$config$dna_channel,
                               sample_dna_height_channel, identity_column)
      if (!is.data.frame(all_data) ||
          any(!required_all_events %in% names(all_data)) || !nrow(all_data)) {
        detail <- character()
        if (!is.data.frame(all_data)) {
          detail <- c(detail, "the all-events export is not a data frame")
        } else {
          missing_all_events <- setdiff(required_all_events, names(all_data))
          if (length(missing_all_events)) detail <- c(detail, paste0(
            "missing column(s): ", paste(missing_all_events, collapse = ", ")
          ))
          if (!nrow(all_data)) detail <- c(detail, "zero events")
        }
        stop("Validated all-events DNA-A/DNA-H coordinates are unavailable for ",
             prefix, ": ", paste(detail, collapse = "; "), ".", call. = FALSE)
      }
      single_x_limits <- robust_axis_limits(
        all_data[[analysis$config$dna_channel]]
      )
      single_y_limits <- robust_axis_limits(all_data[[sample_dna_height_channel]])
      single_axis_coverage <- mean(
        is.finite(all_data[[analysis$config$dna_channel]]) &
          is.finite(all_data[[sample_dna_height_channel]]) &
          all_data[[analysis$config$dna_channel]] >= single_x_limits[[1L]] &
          all_data[[analysis$config$dna_channel]] <= single_x_limits[[2L]] &
          all_data[[sample_dna_height_channel]] >= single_y_limits[[1L]] &
          all_data[[sample_dna_height_channel]] <= single_y_limits[[2L]]
      )
      if (!identity_column %in% names(single) ||
          anyDuplicated(all_data[[identity_column]]) ||
          anyDuplicated(single[[identity_column]])) {
        stop("Single Cells display requires unique direct event identities for ",
             prefix, ".", call. = FALSE)
      }
      if (any(!single[[identity_column]] %in% all_data[[identity_column]])) {
        stop("Single Cells event identities are not contained in the explicit all-events acquisition for ",
             prefix, ".", call. = FALSE)
      }
      all_context <- all_data[
        !all_data[[identity_column]] %in% single[[identity_column]], , drop = FALSE
      ]
      single_plot <- ggplot2::ggplot() +
        ggplot2::geom_point(
          data = retain_points(all_context),
          ggplot2::aes(x = .data[[analysis$config$dna_channel]],
                       y = .data[[sample_dna_height_channel]]),
          colour = "#4D4D4D", shape = 1, size = 0.34, stroke = 0.25,
          alpha = 0.65
        ) +
        ggplot2::geom_point(
          data = retain_points(single),
          ggplot2::aes(x = .data[[analysis$config$dna_channel]],
                       y = .data[[sample_dna_height_channel]]),
          colour = "#0072B2", shape = 16, size = 0.42, stroke = 0,
          alpha = 0.78
        ) +
        ggplot2::labs(title = "Single Cells", subtitle = identity,
                      x = "DNA-A", y = "DNA-H") +
        ggplot2::coord_cartesian(
          xlim = single_x_limits, ylim = single_y_limits
        )
      single_x_axis <- "DNA-A"
      single_y_axis <- "DNA-H"
    } else {
      single_plot <- ggplot2::ggplot(
        retain_points(single),
        ggplot2::aes(x = .data[[analysis$config$dna_channel]],
                     y = .data[[analysis$config$target_channel]])) +
        ggplot2::geom_point(colour = "#0072B2", shape = 16,
                            size = 0.42, stroke = 0, alpha = 0.78) +
        ggplot2::labs(title = "Single Cells", subtitle = identity,
                      x = "DNA", y = facspseudocolor::facs_report_edu_target_label(analysis))
      single_x_axis <- "DNA"
      single_y_axis <- facspseudocolor::facs_report_edu_target_label(analysis)
      single_x_limits <- c(NA_real_, NA_real_)
      single_y_limits <- c(NA_real_, NA_real_)
      single_axis_coverage <- 1
    }
    single_plot <- single_plot +
      ggplot2::theme_classic(base_size = 8) +
      ggplot2::theme(aspect.ratio = 1,
                     plot.title = ggplot2::element_text(face = "bold"))
    g1_plot <- overlay_plot(
      single, g1, "G1", "#0072B2", identity,
      analysis$config$dna_channel, analysis$config$target_channel
    )
    positive_plot <- overlay_plot(
      single, positive, paste0(facspseudocolor::facs_report_edu_target_label(analysis), "+"),
      "#D55E00", identity,
      analysis$config$dna_channel, analysis$config$target_channel
    )
    background_curve <- single[
      is.finite(single[[analysis$config$dna_channel]]) &
        is.finite(single$baseline) & single$baseline > 0,
      c(analysis$config$dna_channel, "baseline"), drop = FALSE
    ]
    background_curve <- background_curve[
      order(background_curve[[analysis$config$dna_channel]]), , drop = FALSE
    ]
    if (nrow(background_curve) < 2L) {
      stop("Fitted EdU background coordinates are unavailable for ", prefix,
           ".", call. = FALSE)
    }
    comparison_parent <- single[
      is.finite(single[[analysis$config$dna_channel]]) &
        is.finite(single[[analysis$config$target_channel]]) &
        single[[analysis$config$target_channel]] > 0 &
        is.finite(single$target_bgsub + display_offset) &
        single$target_bgsub + display_offset > 0, , drop = FALSE
    ]
    if (nrow(comparison_parent) < 2L) {
      stop("Too few common log-display-valid background comparison events for ",
           prefix, ".", call. = FALSE)
    }
    raw_background_plot <- ggplot2::ggplot() +
      ggplot2::geom_point(
        data = retain_points(comparison_parent),
        ggplot2::aes(x = .data[[analysis$config$dna_channel]],
                     y = .data[[analysis$config$target_channel]]),
        colour = "#000000", size = 0.25, stroke = 0, alpha = 0.45
      ) +
      ggplot2::geom_line(
        data = background_curve,
        ggplot2::aes(x = .data[[analysis$config$dna_channel]], y = .data$baseline),
        colour = "#542788", linewidth = 0.55, linetype = "dashed"
      ) +
      ggplot2::labs(title = paste("Raw", facspseudocolor::facs_report_edu_target_label(analysis), "+ fitted background"), subtitle = identity,
                    x = "DNA", y = facspseudocolor::facs_report_edu_target_label(analysis)) +
      ggplot2::scale_y_log10() +
      ggplot2::theme_classic(base_size = 8) +
      ggplot2::theme(aspect.ratio = 1,
                     plot.title = ggplot2::element_text(face = "bold"))
    corrected <- comparison_parent
    corrected$displayed_signal <- corrected$target_bgsub + display_offset
    corrected_plot <- ggplot2::ggplot(
      retain_points(corrected),
      ggplot2::aes(x = .data[[analysis$config$dna_channel]],
                   y = .data$displayed_signal)
    ) +
      ggplot2::geom_point(
        colour = "#4D4D4D", size = 0.25, stroke = 0, alpha = 0.55
      ) +
      ggplot2::labs(
        title = paste("Background-subtracted", facspseudocolor::facs_report_edu_target_label(analysis), "+ offset"), subtitle = identity,
        x = "DNA", y = facspseudocolor::facs_report_edu_target_label(analysis)
      ) +
      ggplot2::scale_y_log10() +
      ggplot2::theme_classic(base_size = 8) +
      ggplot2::theme(aspect.ratio = 1,
                     plot.title = ggplot2::element_text(face = "bold"))
    background_comparison_plot <- cowplot::plot_grid(
      raw_background_plot, corrected_plot, nrow = 1L, align = "hv"
    )
    list(
      prefix = prefix, condition = condition,
      plot = cowplot::plot_grid(single_plot, g1_plot, positive_plot,
                                nrow = 1L, align = "hv"),
      background_plot = background_comparison_plot,
      single_cells_n = nrow(single), g1_n = nrow(g1),
      edu_positive_n = nrow(positive), displayed_point_limit = max_points,
      single_cells_x_axis = single_x_axis,
      single_cells_y_axis = single_y_axis,
      dna_area_source_channel = dna_area_source_channel,
      dna_height_source_channel = sample_dna_height_channel,
      single_cells_x_limits = single_x_limits,
      single_cells_y_limits = single_y_limits,
      single_cells_axis_coverage = single_axis_coverage,
      overlay_x_axis = "DNA", overlay_y_axis = facspseudocolor::facs_report_edu_target_label(analysis),
      background_line_source = "sample_fitted_raw_edu_baseline",
      background_display_offset = display_offset
    )
  })
}

facs_report_edu_s_phase_cards <- function(analysis, report,
                                          max_points = 3000L) {
  required_gates <- c("Early S", "Mid S", "Late S")
  gates <- analysis$quantitation$gates
  prefixes <- as.character(analysis$sample_manifest$prefix)
  if (!inherits(analysis, "facs_analysis") ||
      !identical(analysis$config$plot_type, "edu") ||
      !is.list(report) || !is.data.frame(report$panel_qc) ||
      !is.data.frame(report$display_offset_qc) ||
      !is.data.frame(gates) ||
      !all(c("gate", "xmin", "xmax") %in% names(gates)) ||
      !all(required_gates %in% gates$gate) ||
      !identical(as.character(report$panel_qc$prefix), prefixes) ||
      !identical(as.character(report$display_offset_qc$prefix), prefixes) ||
      !is.numeric(max_points) || length(max_points) != 1L ||
      is.na(max_points) || !is.finite(max_points) || max_points < 1 ||
      max_points != floor(max_points)) {
    stop("S-phase cards require validated EdU panels, offsets, and established regional gates.",
         call. = FALSE)
  }
  max_points <- as.integer(max_points)
  retain_points <- function(data) {
    if (nrow(data) <= max_points) return(data)
    data[unique(as.integer(round(seq(1, nrow(data), length.out = max_points)))),
         , drop = FALSE]
  }
  regional_gates <- gates[match(required_gates, gates$gate), , drop = FALSE]
  colours <- c("Early S" = "#0B4BFF", "Mid S" = "#12BED0",
               "Late S" = "#A026A3")
  lapply(seq_along(prefixes), function(i) {
    prefix <- prefixes[[i]]
    data <- analysis$normalized_data[[i]]$data
    offset <- report$display_offset_qc$display_offset[[i]]
    lower <- report$panel_qc$y_limit_lower[[i]]
    upper <- report$panel_qc$y_limit_upper[[i]]
    if (!is.data.frame(data) ||
        !all(c("dna_norm", "target_bgsub", "edu_computed_positive") %in%
             names(data)) ||
        !is.finite(offset) || !is.finite(lower) || !is.finite(upper) ||
        lower >= upper) {
      stop("S-phase card inputs are unavailable for ", prefix, ".",
           call. = FALSE)
    }
    display <- data.frame(
      dna_norm = data$dna_norm,
      displayed_signal = data$target_bgsub + offset,
      edu_computed_positive = data$edu_computed_positive,
      phase = "Other", stringsAsFactors = FALSE
    )
    for (gate_index in seq_len(nrow(regional_gates))) {
      gate <- regional_gates[gate_index, , drop = FALSE]
      selected <- display$phase == "Other" &
        display$edu_computed_positive %in% TRUE &
        is.finite(display$dna_norm) &
        display$dna_norm >= gate$xmin[[1L]] &
        display$dna_norm < gate$xmax[[1L]]
      display$phase[selected] <- as.character(gate$gate[[1L]])
    }
    assignment_counts <- table(factor(display$phase, levels = required_gates))
    keep <- is.finite(display$dna_norm) & is.finite(display$displayed_signal) &
      display$dna_norm >= analysis$config$x_limits[[1L]] &
      display$dna_norm <= analysis$config$x_limits[[2L]] &
      display$displayed_signal >= lower & display$displayed_signal <= upper
    if (isTRUE(analysis$config$y_log10)) {
      keep <- keep & display$displayed_signal > 0
    }
    display <- display[keep, , drop = FALSE]
    background <- retain_points(display[display$phase == "Other", , drop = FALSE])
    regional <- retain_points(display[display$phase %in% required_gates, , drop = FALSE])
    regional$phase <- factor(regional$phase, levels = required_gates)
    identity <- paste0(analysis$sample_manifest$condition[[i]],
                       " | Sample: ", prefix)
    plot <- ggplot2::ggplot() +
      ggplot2::geom_point(
        data = background,
        ggplot2::aes(x = dna_norm, y = displayed_signal),
        colour = "#D0D0D0", size = 0.25, stroke = 0
      ) +
      ggplot2::geom_point(
        data = regional,
        ggplot2::aes(x = dna_norm, y = displayed_signal, colour = phase),
        size = 0.3, stroke = 0
      ) +
      ggplot2::scale_colour_manual(values = colours, breaks = required_gates,
                                   drop = FALSE, name = NULL) +
      ggplot2::scale_x_continuous(
        breaks = c(analysis$config$dna_2n_value,
                   2 * analysis$config$dna_2n_value),
        labels = c("2N", "4N")
      ) +
      ggplot2::coord_cartesian(xlim = as.numeric(analysis$config$x_limits),
                               ylim = c(lower, upper)) +
      ggplot2::labs(title = prefix, subtitle = identity, x = "DNA", y = facspseudocolor::facs_report_edu_target_label(analysis)) +
      ggplot2::theme_classic(base_size = 9) +
      ggplot2::theme(aspect.ratio = 1, legend.position = "bottom",
                     plot.title = ggplot2::element_text(face = "bold"))
    if (isTRUE(analysis$config$y_log10)) {
      plot <- plot + ggplot2::scale_y_log10()
    }
    list(prefix = prefix,
         condition = as.character(analysis$sample_manifest$condition[[i]]),
         plot = plot,
         early_s_n = unname(as.integer(assignment_counts[["Early S"]])),
         mid_s_n = unname(as.integer(assignment_counts[["Mid S"]])),
         late_s_n = unname(as.integer(assignment_counts[["Late S"]])),
         displayed_point_limit = max_points)
  })
}

facs_report_edu_responsive_groups <- function(
    analysis, report, max_condition_columns = 4L, min_panel_width = 220L
) {
  validate_positive_integer <- function(value, name) {
    if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
        !is.finite(value) || value < 1 || value > .Machine$integer.max ||
        value != floor(value)) {
      stop(paste0("`", name, "` must be one positive integer."), call. = FALSE)
    }
    as.integer(value)
  }
  max_condition_columns <- validate_positive_integer(
    max_condition_columns, "max_condition_columns"
  )
  min_panel_width <- validate_positive_integer(min_panel_width, "min_panel_width")
  manifest <- analysis$sample_manifest
  required <- c(
    "prefix", "condition", "condition_index", "replicate", "replicate_index",
    "technical_replicate", "model_group", "is_reference", "reference_condition"
  )
  missing_fields <- setdiff(required, names(manifest))
  if (length(missing_fields)) {
    stop(paste("Responsive overview requires explicit manifest fields:",
               paste(missing_fields, collapse = ", ")), call. = FALSE)
  }
  identity_fields <- c("prefix", "condition", "replicate", "technical_replicate", "model_group")
  invalid <- vapply(identity_fields, function(field) {
    value <- as.character(manifest[[field]])
    anyNA(value) || any(!nzchar(value))
  }, logical(1))
  prefixes <- as.character(manifest$prefix)
  if (any(invalid) || anyNA(manifest$condition_index) ||
      anyNA(manifest$replicate_index) || anyDuplicated(prefixes)) {
    stop("Responsive overview requires unique prefixes and complete explicit normalization identities.", call. = FALSE)
  }
  if (is.null(names(report$panels)) || !identical(names(report$panels), prefixes)) {
    stop("Responsive overview panels must be named in exact manifest-prefix order.", call. = FALSE)
  }
  panel_qc <- report$panel_qc
  if (!is.data.frame(panel_qc) ||
      !all(c("prefix", "panel_status", "y_limit_lower", "y_limit_upper") %in%
           names(panel_qc)) ||
      !identical(as.character(panel_qc$prefix), prefixes) ||
      anyNA(panel_qc$panel_status) ||
      any(!panel_qc$panel_status %in% c("available", "suppressed"))) {
    stop("Responsive overview requires panel QC and display limits in exact manifest-prefix order.", call. = FALSE)
  }
  condition_order <- unique(as.character(
    manifest$condition[order(manifest$condition_index, seq_len(nrow(manifest)))]
  ))
  overview_panel <- function(panel, shared_y_limits = NULL) {
    out <- panel + ggplot2::labs(
      title = NULL, subtitle = NULL,
      x = "DNA", y = facspseudocolor::facs_report_edu_target_label(analysis)
    ) +
    ggplot2::theme(legend.position = "none")
    if (!is.null(shared_y_limits)) {
      if (!is.numeric(shared_y_limits) || length(shared_y_limits) != 2L ||
          any(!is.finite(shared_y_limits)) || shared_y_limits[[1L]] >= shared_y_limits[[2L]]) {
        stop("Shared overview y-axis limits must be two increasing finite numbers.", call. = FALSE)
      }
      # Replace only the Cartesian viewport through ggplot2's supported API.
      # Carry the existing x viewport, expansion, and clipping forward; y
      # scales, transformations, breaks, data, and the source audit plot remain
      # unchanged. `default = TRUE` makes this intentional replacement quiet.
      existing_coordinates <- out$coordinates
      out <- out + ggplot2::coord_cartesian(
        xlim = existing_coordinates$limits$x,
        ylim = as.numeric(shared_y_limits),
        expand = existing_coordinates$expand,
        default = TRUE,
        clip = existing_coordinates$clip
      )
    }
    out
  }
  placeholder <- function(label) ggplot2::ggplot() +
    ggplot2::annotate("text", .5, .56, label = "MISSING", fontface = "bold",
                      colour = "#9b2c2c", size = 5) +
    ggplot2::annotate("text", .5, .43, label = label, colour = "#5f3030", size = 3) +
    ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) + ggplot2::theme_void()
  group_ids <- unique(as.character(manifest$model_group))
  groups <- lapply(seq_along(group_ids), function(group_index) {
    group_id <- group_ids[[group_index]]
    rows <- manifest[as.character(manifest$model_group) == group_id, , drop = FALSE]
    replicate_index <- unique(rows$replicate_index)
    replicate <- unique(as.character(rows$replicate))
    technical <- unique(as.character(rows$technical_replicate))
    if (length(replicate_index) != 1L || length(replicate) != 1L ||
        length(technical) != 1L) {
      stop("Each model_group must identify exactly one biological replicate and technical acquisition.", call. = FALSE)
    }
    if (anyNA(rows$is_reference)) {
      stop("Each paired acquisition group requires nonmissing is_reference values.", call. = FALSE)
    }
    declared_reference <- as.character(rows$reference_condition)
    if (anyNA(declared_reference) || any(!nzchar(declared_reference)) ||
        length(unique(declared_reference)) != 1L) {
      stop("Every row in a paired acquisition group must declare the same nonempty reference_condition.", call. = FALSE)
    }
    reference_match <- as.character(rows$condition) == declared_reference[[1L]]
    if (sum(reference_match) != 1L ||
        !identical(rows$is_reference %in% TRUE, reference_match)) {
      stop("Each paired acquisition group must have exactly one condition matching reference_condition and the same row marked is_reference.", call. = FALSE)
    }
    reference <- rows[reference_match, , drop = FALSE]
    group_panel_index <- match(as.character(rows$prefix), prefixes)
    available_index <- group_panel_index[
      panel_qc$panel_status[group_panel_index] == "available"
    ]
    shared_y_limits <- NULL
    if (length(available_index)) {
      lower <- panel_qc$y_limit_lower[available_index]
      upper <- panel_qc$y_limit_upper[available_index]
      if (!is.numeric(lower) || !is.numeric(upper)) {
        stop("Each available overview panel requires numeric retained lower and upper display y-axis limits.", call. = FALSE)
      }
      if (any(!is.finite(lower)) || any(!is.finite(upper))) {
        stop("Each available overview panel requires finite retained lower and upper display y-axis limits.", call. = FALSE)
      }
      invalid_pair <- lower >= upper
      if (any(invalid_pair)) {
        invalid_prefixes <- panel_qc$prefix[available_index[invalid_pair]]
        stop(
          "Each available overview panel requires a retained lower display y-axis limit strictly below its upper limit; invalid panel(s): ",
          paste(invalid_prefixes, collapse = ", "), ".",
          call. = FALSE
        )
      }
      shared_y_limits <- c(min(lower), max(upper))
    }
    cards <- unlist(lapply(condition_order, function(condition) {
      found <- rows[rows$condition == condition, , drop = FALSE]
      if (!nrow(found)) return(list(list(
        condition = condition, prefix = NA_character_, status = "missing",
        is_reference = identical(condition, as.character(reference$condition[[1L]])),
        plot = placeholder(paste(replicate, technical, condition, sep = " / "))
      )))
      found <- found[order(found$condition_index, match(found$prefix, prefixes)), , drop = FALSE]
      lapply(seq_len(nrow(found)), function(i) {
        prefix <- as.character(found$prefix[[i]])
        panel_index <- match(prefix, prefixes)
        status <- as.character(panel_qc$panel_status[[panel_index]])
        list(
          condition = as.character(found$condition[[i]]),
          prefix = prefix, status = status,
          is_reference = found$is_reference[[i]] %in% TRUE,
          plot = overview_panel(
            report$panels[[prefix]],
            if (identical(status, "available")) shared_y_limits else NULL
          )
        )
      })
    }), recursive = FALSE)
    reference_card <- vapply(cards, `[[`, logical(1), "is_reference")
    cards <- c(cards[reference_card], cards[!reference_card])
    card_prefixes <- vapply(cards, `[[`, character(1), "prefix")
    candidate <- match(card_prefixes[!is.na(card_prefixes)], prefixes)
    candidate <- candidate[panel_qc$panel_status[candidate] == "available"]
    legend_source_prefix <- if (length(candidate)) prefixes[[candidate[[1L]]]] else NA_character_
    legend <- if (length(candidate)) cowplot::get_legend(
      report$panels[[candidate[[1L]]]] +
        ggplot2::labs(colour = NULL) +
        ggplot2::theme(legend.position = "bottom")
    ) else NULL
    list(
      group_index = group_index, model_group = group_id,
      replicate = replicate[[1L]], replicate_index = replicate_index[[1L]],
      technical_replicate = technical[[1L]],
      reference_condition = as.character(reference$condition[[1L]]),
      shared_y_limits = shared_y_limits,
      cards = cards, legend = legend,
      legend_source_prefix = legend_source_prefix,
      legend_status = if (!length(candidate)) "unavailable_all_group_panels_suppressed" else
        if (is.null(legend)) "unavailable_from_available_panel" else "available",
      shared_legend_count = if (is.null(legend)) 0L else 1L
    )
  })
  list(
    groups = groups, condition_order = condition_order,
    max_condition_columns = max_condition_columns,
    min_panel_width = min_panel_width,
    balance = facs_report_edu_balance(analysis)
  )
}

# pH3 Standard v1 report helpers -----------------------------------------
#
# These mirror the EdU Standard v2 report helpers above, adapted for
# plot_type: ph3. Unlike EdU, pH3-legacy mode (ph3_positivity_method:
# flowjo_legacy_v1) has no per-sample fitted background model and no
# reference-relative replicate pairing (PH3 replicates have no `reference`
# field), so these are deliberately simpler than their EdU counterparts:
# no model_group/is_reference/shared-legend logic, and the "normalization"
# shown is DNA-content normalization (raw -> G1-anchored, mapped to the
# configured 2N value), not a background-subtracted target signal.

facs_report_ph3_overview <- function(analysis) {
  facspseudocolor:::validate_analysis_object(analysis)
  if (!identical(analysis$config$plot_type, "ph3")) {
    stop("The pH3 report overview requires a pH3 analysis.", call. = FALSE)
  }
  manifest <- analysis$sample_manifest
  input <- analysis$input_report
  target_name <- analysis$config$target_name
  gates <- analysis$quantitation$ph3$gates
  populations <- c("FlowJo Single Cells", "FlowJo G1", paste0("FlowJo ", target_name, "+"))
  list(
    summary = data.frame(
      item = c(
        "Samples / acquisitions", "Biological replicates", "Conditions",
        "DNA channel", paste(target_name, "channel"), "Population sources",
        "DNA normalization", "Phase boundaries (normalized DNA)",
        "Positivity method", "Configuration"
      ),
      value = c(
        nrow(manifest), length(unique(manifest$replicate_index)),
        paste(unique(manifest$condition[order(manifest$condition_index)]), collapse = ", "),
        analysis$config$dna_channel, analysis$config$target_channel,
        paste(populations, collapse = ", "),
        paste0("G1/2N mapped to ", analysis$config$dna_2n_value,
               " using the configured ", analysis$config$g1_anchor, " anchor"),
        if (is.data.frame(gates) && nrow(gates)) {
          paste(sprintf("%s [%s–%s]", gates$gate, gates$xmin, gates$xmax), collapse = ", ")
        } else "Not available",
        analysis$config$ph3_positivity_method,
        if (is.character(analysis$provenance$config_path) &&
            length(analysis$provenance$config_path) == 1L &&
            !is.na(analysis$provenance$config_path) &&
            nzchar(analysis$provenance$config_path)) {
          basename(analysis$provenance$config_path)
        } else "Completed analysis artifact"
      ), stringsAsFactors = FALSE
    ),
    samples = manifest[, intersect(c(
      "replicate", "technical_replicate", "condition", "prefix", "is_reference"
    ), names(manifest)), drop = FALSE],
    inputs = input
  )
}

#' Arrange pH3 pseudocolor panels into a replicate x condition grid
#'
#' Simpler than facs_report_edu_responsive_groups: pH3 replicates have no
#' `reference`/`is_reference`/`model_group` pairing concept, so this groups
#' purely by biological replicate, in configured condition order.
facs_report_ph3_responsive_groups <- function(
    analysis, panels, max_condition_columns = 4L, min_panel_width = 220L
) {
  manifest <- analysis$sample_manifest
  prefixes <- as.character(manifest$prefix)
  if (is.null(names(panels)) || !identical(names(panels), prefixes)) {
    stop("pH3 responsive overview panels must be named in exact manifest-prefix order.", call. = FALSE)
  }
  condition_order <- unique(as.character(
    manifest$condition[order(manifest$condition_index, seq_len(nrow(manifest)))]
  ))
  replicate_index_order <- unique(manifest$replicate_index[order(manifest$replicate_index)])
  groups <- lapply(seq_along(replicate_index_order), function(group_index) {
    replicate_index <- replicate_index_order[[group_index]]
    rows <- manifest[manifest$replicate_index == replicate_index, , drop = FALSE]
    replicate <- unique(as.character(rows$replicate))
    if (length(replicate) != 1L) {
      stop("Each pH3 biological replicate index must have exactly one replicate label.", call. = FALSE)
    }
    cards <- lapply(condition_order, function(condition) {
      found <- rows[rows$condition == condition, , drop = FALSE]
      if (!nrow(found)) {
        return(list(
          condition = condition, prefix = NA_character_, status = "missing",
          plot = ggplot2::ggplot() +
            ggplot2::annotate("text", .5, .5, label = "MISSING", fontface = "bold",
                              colour = "#9b2c2c", size = 5) +
            ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) + ggplot2::theme_void()
        ))
      }
      prefix <- as.character(found$prefix[[1L]])
      list(condition = condition, prefix = prefix, status = "available",
           plot = panels[[prefix]])
    })
    list(
      group_index = group_index, replicate = replicate[[1L]],
      replicate_index = replicate_index, cards = cards
    )
  })
  list(
    groups = groups, condition_order = condition_order,
    max_condition_columns = max_condition_columns, min_panel_width = min_panel_width
  )
}

#' Single Cells / G1 / pH3+ gating cards, one three-panel card per sample
#'
#' Mirrors facs_report_edu_gating_cards's shape and CSS conventions, adapted
#' for pH3: the Single Cells panel shades Single Cells over all acquisition
#' events on raw DNA-A versus DNA-H (mandatory, same convention as EdU/
#' CDC45 reports); the G1 and pH3+ panels shade those populations over
#' Single Cells on normalized DNA versus target signal. Unlike EdU, event
#' identities are not required to be unique/matched between all-events and
#' gated exports (this dataset's legacy exports carry no event_index), so
#' overlay points are layered on top of the full Single Cells context
#' rather than excluded from it -- visually equivalent, since the colored
#' overlay points are drawn on top regardless.
facs_report_ph3_gating_cards <- function(analysis, all_events = NULL,
                                         dna_height_channel = NULL,
                                         max_points = 3000L) {
  if (!inherits(analysis, "facs_analysis") ||
      !identical(analysis$config$plot_type, "ph3") ||
      length(analysis$normalized_data) != nrow(analysis$sample_manifest) ||
      !is.numeric(max_points) || length(max_points) != 1L ||
      is.na(max_points) || !is.finite(max_points) || max_points < 1 ||
      max_points != floor(max_points)) {
    stop("pH3 gating cards require one valid pH3 analysis and a positive point limit.",
         call. = FALSE)
  }
  max_points <- as.integer(max_points)
  prefixes <- as.character(analysis$sample_manifest$prefix)
  if (!is.list(all_events) ||
      !identical(names(all_events), prefixes)) {
    stop("pH3 DNA-A/DNA-H gating displays require explicit all-events exports in exact sample order.",
         call. = FALSE)
  }
  if (!is.character(dna_height_channel) || length(dna_height_channel) != 1L ||
      is.na(dna_height_channel) || !nzchar(dna_height_channel)) {
    stop("pH3 DNA-A/DNA-H gating displays require an explicit DNA-height channel.",
         call. = FALSE)
  }
  retain_points <- function(data) {
    if (nrow(data) <= max_points) return(data)
    data[unique(as.integer(round(seq(1, nrow(data), length.out = max_points)))),
         , drop = FALSE]
  }
  robust_axis_limits <- function(values) {
    values <- values[is.finite(values)]
    if (length(values) < 2L) {
      stop("DNA-A/DNA-H display limits require at least two finite values.", call. = FALSE)
    }
    limits <- as.numeric(stats::quantile(values, c(0.001, 0.999), names = FALSE, type = 7))
    span <- diff(limits)
    if (!all(is.finite(limits)) || span <= 0) {
      stop("DNA-A/DNA-H display limits are degenerate.", call. = FALSE)
    }
    limits + c(-1, 1) * span * 0.03
  }
  target_name <- analysis$config$target_name
  overlay_plot <- function(parent, child, title, child_colour, identity) {
    parent <- parent[is.finite(parent$target_norm) & parent$target_norm > 0, , drop = FALSE]
    child <- child[is.finite(child$target_norm) & child$target_norm > 0, , drop = FALSE]
    ggplot2::ggplot() +
      ggplot2::geom_point(
        data = retain_points(parent),
        ggplot2::aes(x = dna_norm, y = target_norm),
        colour = "#4D4D4D", shape = 1, size = 0.34, stroke = 0.25, alpha = 0.65
      ) +
      ggplot2::geom_point(
        data = retain_points(child),
        ggplot2::aes(x = dna_norm, y = target_norm),
        colour = child_colour, shape = 16, size = 0.42, stroke = 0, alpha = 0.78
      ) +
      ggplot2::labs(title = title, subtitle = identity, x = "DNA", y = target_name) +
      ggplot2::scale_y_log10() +
      ggplot2::theme_classic(base_size = 8) +
      ggplot2::theme(aspect.ratio = 1, plot.title = ggplot2::element_text(face = "bold"))
  }
  manifest <- analysis$sample_manifest
  lapply(seq_len(nrow(manifest)), function(i) {
    sample <- analysis$normalized_data[[i]]
    prefix <- as.character(manifest$prefix[[i]])
    condition <- as.character(manifest$condition[[i]])
    identity <- paste0(condition, " | Sample: ", prefix)
    single <- sample$data
    g1 <- sample$g1
    positive <- sample$ph3_positive
    if (!is.data.frame(single) || !is.data.frame(g1) || !is.data.frame(positive) ||
        !nrow(single) || !nrow(g1) || !nrow(positive)) {
      stop("Validated Single Cells, G1, and ", target_name, "-positive inputs are required for gating card ",
           prefix, ".", call. = FALSE)
    }
    all_data <- all_events[[prefix]]
    dna_channel <- analysis$config$dna_channel
    required_all_events <- c(dna_channel, dna_height_channel)
    if (!is.data.frame(all_data) ||
        any(!required_all_events %in% names(all_data)) || !nrow(all_data)) {
      stop("Validated all-events DNA-A/DNA-H coordinates are unavailable for ", prefix, ".",
           call. = FALSE)
    }
    single_x_limits <- robust_axis_limits(all_data[[dna_channel]])
    single_y_limits <- robust_axis_limits(all_data[[dna_height_channel]])
    single_plot <- ggplot2::ggplot() +
      ggplot2::geom_point(
        data = retain_points(all_data),
        ggplot2::aes(x = .data[[dna_channel]], y = .data[[dna_height_channel]]),
        colour = "#4D4D4D", shape = 1, size = 0.34, stroke = 0.25, alpha = 0.65
      ) +
      ggplot2::geom_point(
        data = retain_points(single),
        ggplot2::aes(x = .data[[dna_channel]], y = .data[[dna_height_channel]]),
        colour = "#0072B2", shape = 16, size = 0.42, stroke = 0, alpha = 0.78
      ) +
      ggplot2::labs(title = "Single Cells", subtitle = identity, x = "DNA-A", y = "DNA-H") +
      ggplot2::coord_cartesian(xlim = single_x_limits, ylim = single_y_limits) +
      ggplot2::theme_classic(base_size = 8) +
      ggplot2::theme(aspect.ratio = 1, plot.title = ggplot2::element_text(face = "bold"))
    g1_plot <- overlay_plot(single, g1, "G1", "#0072B2", identity)
    positive_plot <- overlay_plot(single, positive, paste0(target_name, "+"), "#D55E00", identity)
    list(
      prefix = prefix, condition = condition,
      plot = cowplot::plot_grid(single_plot, g1_plot, positive_plot, nrow = 1L, align = "hv"),
      single_cells_n = nrow(single), g1_n = nrow(g1), ph3_positive_n = nrow(positive),
      displayed_point_limit = max_points
    )
  })
}

#' Raw versus G1-anchored-normalized DNA content, one card per sample
#'
#' The only real "normalization" step in pH3-legacy mode: raw DNA content
#' is anchored to the G1 median and mapped to the configured dna_2n_value
#' (2N). There is no background-subtracted target signal to compare in
#' this mode (ph3_positivity_method: flowjo_legacy_v1 reuses the FlowJo
#' pH3+ gate directly rather than fitting a background model), so unlike
#' EdU's background-correction cards, this compares the DNA axis itself
#' before and after normalization. Configured phase boundaries are drawn
#' as dashed reference lines on the normalized panel.
facs_report_ph3_normalization_cards <- function(analysis, max_points = 3000L) {
  if (!inherits(analysis, "facs_analysis") ||
      !identical(analysis$config$plot_type, "ph3") ||
      length(analysis$normalized_data) != nrow(analysis$sample_manifest)) {
    stop("pH3 normalization cards require one valid pH3 analysis.", call. = FALSE)
  }
  retain_points <- function(data) {
    if (nrow(data) <= max_points) return(data)
    data[unique(as.integer(round(seq(1, nrow(data), length.out = max_points)))),
         , drop = FALSE]
  }
  gates <- analysis$quantitation$ph3$gates
  manifest <- analysis$sample_manifest
  dna_channel <- analysis$config$dna_channel
  lapply(seq_len(nrow(manifest)), function(i) {
    sample <- analysis$normalized_data[[i]]$data
    prefix <- as.character(manifest$prefix[[i]])
    condition <- as.character(manifest$condition[[i]])
    identity <- paste0(condition, " | Sample: ", prefix)
    if (!is.data.frame(sample) || !nrow(sample) ||
        any(!c(dna_channel, "dna_norm") %in% names(sample))) {
      stop("Single Cells raw and normalized DNA content are required for ", prefix, ".",
           call. = FALSE)
    }
    raw_plot <- ggplot2::ggplot(retain_points(sample), ggplot2::aes(x = .data[[dna_channel]])) +
      ggplot2::geom_histogram(bins = 80L, fill = "#4D4D4D", colour = NA) +
      ggplot2::labs(title = "Raw DNA content", subtitle = identity, x = "DNA (raw)", y = "Events") +
      ggplot2::theme_classic(base_size = 8) +
      ggplot2::theme(aspect.ratio = 0.8, plot.title = ggplot2::element_text(face = "bold"))
    normalized_plot <- ggplot2::ggplot(retain_points(sample), ggplot2::aes(x = dna_norm))
    if (is.data.frame(gates) && nrow(gates)) {
      # Contiguous phase regions (each gate's xmax is the next gate's xmin)
      # share edges, so plain boundary lines are indistinguishable from one
      # another. Shade each region with a distinct, low-alpha color from
      # the same qualitative palette used for the pseudocolor phase-gate
      # overlay (add_phase_gates_to_plot()) instead, drawn first so the
      # histogram sits on top of it.
      shading <- gates
      shading$fill_color <- facspseudocolor:::facs_named_palette("colorblind", nrow(shading))[
        shading$gate_index %||% seq_len(nrow(shading))
      ]
      normalized_plot <- normalized_plot + ggplot2::geom_rect(
        data = shading,
        ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf,
                    fill = fill_color),
        inherit.aes = FALSE, alpha = 0.16, show.legend = FALSE
      ) + ggplot2::scale_fill_identity()
    }
    normalized_plot <- normalized_plot +
      ggplot2::geom_histogram(bins = 80L, fill = "#0072B2", colour = NA) +
      ggplot2::labs(
        title = paste0("Normalized DNA content (2N = ", analysis$config$dna_2n_value, ")"),
        subtitle = identity, x = "DNA (normalized)", y = "Events"
      ) +
      ggplot2::theme_classic(base_size = 8) +
      ggplot2::theme(aspect.ratio = 0.8, plot.title = ggplot2::element_text(face = "bold"))
    list(
      prefix = prefix, condition = condition,
      plot = cowplot::plot_grid(raw_plot, normalized_plot, nrow = 1L, align = "hv"),
      single_cells_n = nrow(sample)
    )
  })
}

# POI Standard v1 report helpers ------------------------------------------
#
# These mirror the EdU/pH3 Standard report helpers above, adapted for
# plot_type: poi. Unlike EdU/pH3, POI mode has only one per-sample
# population (`complete`/Single Cells -- there is no separate G1 or
# positive-population file), and its background correction is a real
# per-replicate regression against a declared `reference` (background
# control) acquisition, closely mirroring EdU's own background-fit shape
# (`sample$data` carries the same `baseline`/`target_bgsub`/`target_norm`
# fields EdU uses). POI replicates DO have `is_reference`/
# `reference_condition`, so cards mark the configured background control.
#
# Known, deliberately NOT auto-"fixed" scientific caveat (see
# poi_current_implementation_audit_5de77a9.md /
# poi_method_and_output_decision_memo_5de77a9.md in the owner's analysis
# workspace, POI-D06): the `background_quantile` cutoff line is computed in
# background-divided units but drawn, by default, on whatever
# `pseudocolor_signal` is configured -- for the default
# `background_subtracted` display this is a real unit mismatch, not a
# validated positivity boundary. This report does not attempt a "corrected"
# line (there is no valid constant-offset equivalent: the true conversion is
# DNA-position-dependent, not a horizontal line, and would itself be a new
# unapproved scientific method). Instead it surfaces `show_cutoff_line`
# (new, additive `facs_config_keys()` entry, default `TRUE` so existing
# configs are unaffected) so a `background_subtracted`-display config can
# explicitly hide the mismatched line, and the report explains why whenever
# it is hidden for that reason.

facs_report_poi_overview <- function(analysis) {
  facspseudocolor:::validate_analysis_object(analysis)
  if (!identical(analysis$config$plot_type, "poi")) {
    stop("The POI report overview requires a POI analysis.", call. = FALSE)
  }
  manifest <- analysis$sample_manifest
  input <- analysis$input_report
  target_name <- analysis$config$target_name
  signal <- analysis$config$pseudocolor_signal
  signal_label <- if (identical(signal, "background_subtracted")) {
    "background-subtracted (raw − predicted background)"
  } else {
    "background-divided (1000 × raw / predicted background)"
  }
  reference_labels <- unique(as.character(
    manifest$condition[manifest$is_reference %in% TRUE]
  ))
  list(
    summary = data.frame(
      item = c(
        "Samples / acquisitions", "Biological replicates", "Conditions",
        "DNA channel", paste(target_name, "channel"), "Background control",
        "DNA alignment", "Canonical signal", "Background cutoff line",
        "Configuration"
      ),
      value = c(
        nrow(manifest), length(unique(manifest$replicate_index)),
        paste(unique(manifest$condition[order(manifest$condition_index)]), collapse = ", "),
        analysis$config$dna_channel, analysis$config$target_channel,
        if (length(reference_labels)) paste(reference_labels, collapse = ", ") else "Not configured",
        paste0(analysis$config$poi_dna_align, " 2N peak detection, mapped to ",
               analysis$config$dna_2n_value),
        signal_label,
        if (isTRUE(analysis$config$show_cutoff_line)) {
          "Shown (background-divided units -- see Methods)"
        } else {
          "Hidden: not valid on the background-subtracted display (see Methods)"
        },
        if (is.character(analysis$provenance$config_path) &&
            length(analysis$provenance$config_path) == 1L &&
            !is.na(analysis$provenance$config_path) &&
            nzchar(analysis$provenance$config_path)) {
          basename(analysis$provenance$config_path)
        } else "Completed analysis artifact"
      ), stringsAsFactors = FALSE
    ),
    samples = manifest[, intersect(c(
      "replicate", "technical_replicate", "condition", "prefix", "is_reference"
    ), names(manifest)), drop = FALSE],
    inputs = input
  )
}

#' Arrange POI pseudocolor panels into a replicate x condition grid
#'
#' Simpler than facs_report_edu_responsive_groups (no model_group/paired-
#' acquisition or shared-legend logic), but -- unlike
#' facs_report_ph3_responsive_groups -- marks the configured background-
#' control (`is_reference`) card, since POI replicates carry a real
#' reference concept.
facs_report_poi_responsive_groups <- function(
    analysis, panels, max_condition_columns = 4L, min_panel_width = 220L
) {
  manifest <- analysis$sample_manifest
  prefixes <- as.character(manifest$prefix)
  # plot_pseudocolor_panels() omits the reference/background-control
  # acquisition from `panels` by default (config `show_reference_panel`;
  # per the POI audit this setting also changes what gets quantified, so
  # this report does not silently flip it) -- so `panels` is normally a
  # proper subset of the manifest prefixes, not an exact match. Still
  # require whatever names ARE present to be a subsequence in manifest
  # order, to catch a genuine mismatch (wrong analysis object, stale cache).
  if (is.null(names(panels)) || !all(names(panels) %in% prefixes) ||
      !identical(names(panels), prefixes[prefixes %in% names(panels)])) {
    stop("POI responsive overview panels must be named in manifest-prefix order.", call. = FALSE)
  }
  condition_order <- unique(as.character(
    manifest$condition[order(manifest$condition_index, seq_len(nrow(manifest)))]
  ))
  replicate_index_order <- unique(manifest$replicate_index[order(manifest$replicate_index)])
  groups <- lapply(seq_along(replicate_index_order), function(group_index) {
    replicate_index <- replicate_index_order[[group_index]]
    rows <- manifest[manifest$replicate_index == replicate_index, , drop = FALSE]
    replicate <- unique(as.character(rows$replicate))
    if (length(replicate) != 1L) {
      stop("Each POI biological replicate index must have exactly one replicate label.", call. = FALSE)
    }
    cards <- lapply(condition_order, function(condition) {
      found <- rows[rows$condition == condition, , drop = FALSE]
      if (!nrow(found)) {
        return(list(
          condition = condition, prefix = NA_character_, status = "missing",
          is_reference = FALSE,
          plot = ggplot2::ggplot() +
            ggplot2::annotate("text", .5, .5, label = "MISSING", fontface = "bold",
                              colour = "#9b2c2c", size = 5) +
            ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) + ggplot2::theme_void()
        ))
      }
      prefix <- as.character(found$prefix[[1L]])
      is_reference <- found$is_reference[[1L]] %in% TRUE
      if (!prefix %in% names(panels)) {
        # Intentionally omitted (the reference/background-control
        # acquisition, when show_reference_panel is not enabled) -- not an
        # error, so this is explicitly labeled rather than styled as MISSING.
        return(list(
          condition = condition, prefix = prefix, status = "not_plotted",
          is_reference = is_reference,
          plot = ggplot2::ggplot() +
            ggplot2::annotate("text", .5, .58, label = "Not plotted", fontface = "bold",
                              colour = "#46565c", size = 4.5) +
            ggplot2::annotate("text", .5, .42,
                              label = "background control\n(show_reference_panel: false)",
                              colour = "#46565c", size = 3) +
            ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) + ggplot2::theme_void()
        ))
      }
      list(condition = condition, prefix = prefix, status = "available",
           is_reference = is_reference, plot = panels[[prefix]])
    })
    is_ref <- vapply(cards, `[[`, logical(1), "is_reference")
    cards <- c(cards[is_ref], cards[!is_ref])
    list(
      group_index = group_index, replicate = replicate[[1L]],
      replicate_index = replicate_index, cards = cards
    )
  })
  list(
    groups = groups, condition_order = condition_order,
    max_condition_columns = max_condition_columns, min_panel_width = min_panel_width
  )
}

#' Single Cells gating card, one per sample
#'
#' POI has only one population per sample (no G1/positive files), so unlike
#' the EdU/pH3 gating cards this is a single panel: Single Cells shaded over
#' all acquisition events on raw DNA-A versus DNA-H (mandatory, same
#' convention as the EdU/CDC45/pH3 reports).
facs_report_poi_gating_cards <- function(analysis, all_events = NULL,
                                         dna_height_channel = NULL,
                                         max_points = 3000L) {
  if (!inherits(analysis, "facs_analysis") ||
      !identical(analysis$config$plot_type, "poi") ||
      length(analysis$normalized_data) != nrow(analysis$sample_manifest) ||
      !is.numeric(max_points) || length(max_points) != 1L ||
      is.na(max_points) || !is.finite(max_points) || max_points < 1 ||
      max_points != floor(max_points)) {
    stop("POI gating cards require one valid POI analysis and a positive point limit.",
         call. = FALSE)
  }
  max_points <- as.integer(max_points)
  prefixes <- as.character(analysis$sample_manifest$prefix)
  if (!is.list(all_events) || !identical(names(all_events), prefixes)) {
    stop("POI DNA-A/DNA-H gating displays require explicit all-events exports in exact sample order.",
         call. = FALSE)
  }
  if (!is.character(dna_height_channel) || length(dna_height_channel) != 1L ||
      is.na(dna_height_channel) || !nzchar(dna_height_channel)) {
    stop("POI DNA-A/DNA-H gating displays require an explicit DNA-height channel.",
         call. = FALSE)
  }
  retain_points <- function(data) {
    if (nrow(data) <= max_points) return(data)
    data[unique(as.integer(round(seq(1, nrow(data), length.out = max_points)))),
         , drop = FALSE]
  }
  robust_axis_limits <- function(values) {
    values <- values[is.finite(values)]
    if (length(values) < 2L) {
      stop("DNA-A/DNA-H display limits require at least two finite values.", call. = FALSE)
    }
    limits <- as.numeric(stats::quantile(values, c(0.001, 0.999), names = FALSE, type = 7))
    span <- diff(limits)
    if (!all(is.finite(limits)) || span <= 0) {
      stop("DNA-A/DNA-H display limits are degenerate.", call. = FALSE)
    }
    limits + c(-1, 1) * span * 0.03
  }
  manifest <- analysis$sample_manifest
  dna_channel <- analysis$config$dna_channel
  lapply(seq_len(nrow(manifest)), function(i) {
    sample <- analysis$normalized_data[[i]]
    prefix <- as.character(manifest$prefix[[i]])
    condition <- as.character(manifest$condition[[i]])
    identity <- paste0(condition, " | Sample: ", prefix)
    single <- sample$data
    if (!is.data.frame(single) || !nrow(single)) {
      stop("Validated Single Cells input is required for gating card ", prefix, ".", call. = FALSE)
    }
    all_data <- all_events[[prefix]]
    required_all_events <- c(dna_channel, dna_height_channel)
    if (!is.data.frame(all_data) ||
        any(!required_all_events %in% names(all_data)) || !nrow(all_data)) {
      stop("Validated all-events DNA-A/DNA-H coordinates are unavailable for ", prefix, ".", call. = FALSE)
    }
    single_x_limits <- robust_axis_limits(all_data[[dna_channel]])
    single_y_limits <- robust_axis_limits(all_data[[dna_height_channel]])
    single_plot <- ggplot2::ggplot() +
      ggplot2::geom_point(
        data = retain_points(all_data),
        ggplot2::aes(x = .data[[dna_channel]], y = .data[[dna_height_channel]]),
        colour = "#4D4D4D", shape = 1, size = 0.34, stroke = 0.25, alpha = 0.65
      ) +
      ggplot2::geom_point(
        data = retain_points(single),
        ggplot2::aes(x = .data[[dna_channel]], y = .data[[dna_height_channel]]),
        colour = "#0072B2", shape = 16, size = 0.42, stroke = 0, alpha = 0.78
      ) +
      ggplot2::labs(title = "Single Cells", subtitle = identity, x = "DNA-A", y = "DNA-H") +
      ggplot2::coord_cartesian(xlim = single_x_limits, ylim = single_y_limits) +
      ggplot2::theme_classic(base_size = 8) +
      ggplot2::theme(aspect.ratio = 1, plot.title = ggplot2::element_text(face = "bold"))
    list(prefix = prefix, condition = condition, plot = single_plot, single_cells_n = nrow(single))
  })
}

#' Raw-plus-fitted-background versus background-corrected signal cards
#'
#' Mirrors facs_report_edu_gating_cards's background-comparison panel pair
#' exactly (same `baseline`/`target_bgsub` fields, same layout), relabeled
#' for POI's target and corrected-signal terminology. The left plot always
#' shows raw target with the fitted background line; the right plot shows
#' whichever `pseudocolor_signal` is configured -- background-subtracted
#' (default, `target_bgsub` + display offset) or background-divided
#' (`1000 * target_raw / background_fitted`), explicitly labeled either way.
facs_report_poi_background_cards <- function(analysis, max_points = 3000L) {
  if (!inherits(analysis, "facs_analysis") ||
      !identical(analysis$config$plot_type, "poi") ||
      length(analysis$normalized_data) != nrow(analysis$sample_manifest)) {
    stop("POI background-correction cards require one valid POI analysis.", call. = FALSE)
  }
  retain_points <- function(data) {
    if (nrow(data) <= max_points) return(data)
    data[unique(as.integer(round(seq(1, nrow(data), length.out = max_points)))),
         , drop = FALSE]
  }
  manifest <- analysis$sample_manifest
  dna_channel <- analysis$config$dna_channel
  target_name <- analysis$config$target_name
  dna_2n_value <- analysis$config$dna_2n_value
  divided <- !identical(analysis$config$pseudocolor_signal, "background_subtracted")
  lapply(seq_len(nrow(manifest)), function(i) {
    sample <- analysis$normalized_data[[i]]
    single <- sample$data
    prefix <- as.character(manifest$prefix[[i]])
    condition <- as.character(manifest$condition[[i]])
    identity <- paste0(condition, " | Sample: ", prefix)
    if (!is.data.frame(single) ||
        any(!c("dna_norm", "target_raw", "baseline") %in% names(single))) {
      stop("Fitted background coordinates are required for ", prefix, ".", call. = FALSE)
    }
    background_curve <- single[
      is.finite(single$dna_norm) & is.finite(single$baseline) & single$baseline > 0,
      c("dna_norm", "baseline"), drop = FALSE
    ]
    background_curve <- background_curve[order(background_curve$dna_norm), , drop = FALSE]
    if (nrow(background_curve) < 2L) {
      stop("Fitted background coordinates are unavailable for ", prefix, ".", call. = FALSE)
    }
    raw_plot <- ggplot2::ggplot() +
      ggplot2::geom_point(
        data = retain_points(single[is.finite(single$dna_norm) & is.finite(single$target_raw) &
                                     single$target_raw > 0, , drop = FALSE]),
        ggplot2::aes(x = dna_norm, y = target_raw),
        colour = "#000000", size = 0.25, stroke = 0, alpha = 0.45
      ) +
      ggplot2::geom_line(
        data = background_curve, ggplot2::aes(x = dna_norm, y = baseline),
        colour = "#542788", linewidth = 0.55, linetype = "dashed"
      ) +
      ggplot2::labs(title = paste("Raw", target_name, "+ fitted background"), subtitle = identity,
                    x = "DNA", y = target_name) +
      ggplot2::scale_y_log10() +
      ggplot2::theme_classic(base_size = 8) +
      ggplot2::theme(aspect.ratio = 1, plot.title = ggplot2::element_text(face = "bold"))
    if (divided) {
      corrected <- single[is.finite(single$dna_norm) & is.finite(single$baseline) &
                            single$baseline > 0 & is.finite(single$target_raw), , drop = FALSE]
      corrected$displayed_signal <- dna_2n_value * corrected$target_raw / corrected$baseline
      corrected_title <- paste("Background-divided", target_name)
    } else {
      corrected <- single[is.finite(single$dna_norm) & is.finite(single$target_bgsub), , drop = FALSE]
      corrected$displayed_signal <- corrected$target_bgsub
      corrected_title <- paste("Background-subtracted", target_name)
    }
    if (nrow(corrected) < 2L) {
      stop("Too few finite corrected-signal events for ", prefix, ".", call. = FALSE)
    }
    corrected_plot <- ggplot2::ggplot(
      retain_points(corrected), ggplot2::aes(x = dna_norm, y = displayed_signal)
    ) +
      ggplot2::geom_point(colour = "#4D4D4D", size = 0.25, stroke = 0, alpha = 0.55) +
      ggplot2::labs(title = corrected_title, subtitle = identity, x = "DNA", y = target_name) +
      ggplot2::theme_classic(base_size = 8) +
      ggplot2::theme(aspect.ratio = 1, plot.title = ggplot2::element_text(face = "bold"))
    if (!divided) {
      corrected_plot <- corrected_plot +
        ggplot2::scale_y_continuous() +
        ggplot2::labs(caption = "Finite negative values are real (background-subtracted); a symmetric linear axis is used, not log10.")
    } else {
      corrected_plot <- corrected_plot + ggplot2::scale_y_log10()
    }
    list(
      prefix = prefix, condition = condition,
      plot = cowplot::plot_grid(raw_plot, corrected_plot, nrow = 1L, align = "hv"),
      single_cells_n = nrow(single)
    )
  })
}
