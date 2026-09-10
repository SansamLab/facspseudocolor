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
    figure_height = 2.15 + 0.15 * rows
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
    "print(report_plot)",
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
  population_names <- c(
    complete = "FlowJo Single Cells", g1 = "FlowJo G1",
    edu_positive = "FlowJo EdU Positive"
  )
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
        "DNA channel", "EdU channel", "Population sources",
        "DNA normalization", "EdU correction", "Intensity reference",
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
          "Independent EdU-negative baseline slope fitted for each acquisition ",
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
    panel + ggplot2::labs(y = "EdU signal (display scale)") +
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
        ggplot2::labs(colour = "Within-panel relative density") +
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
      shared_legend_title = "Within-panel relative density",
      legend_source_prefix = legend_source_prefix,
      panel_legends_visible = FALSE,
      y_axis_title = "EdU signal (display scale)",
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
      y = "EdU signal (display scale)"
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
        ggplot2::labs(colour = "Within-panel relative density") +
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
