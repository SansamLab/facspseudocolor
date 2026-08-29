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

ph3_configuration_digest <- function(config) {
  value <- as.list(config)
  location_only <- c(
    "data_dir", "output_pdf", "output_png", "ph3_export_operation_dirs"
  )
  value <- value[setdiff(names(value), location_only)]
  listify_vectors <- function(x) {
    if (is.list(x)) return(lapply(x, listify_vectors))
    if (length(x) > 1L) return(lapply(as.list(x), listify_vectors))
    x
  }
  value <- listify_vectors(value)
  digest <- as.character(openssl::sha256(charToRaw(
    paste0(ph3_canonical_json(value), "\n")
  )))
  attributes(digest) <- NULL
  digest
}

ph3_classification_fail <- function(acquisition_id, reason, detail) {
  stop(
    "PH3 event classification failed [", reason, "] for ", acquisition_id,
    ": ", detail, ".", call. = FALSE
  )
}

ph3_analysis_id <- function(config_digest, export_manifests) {
  operations <- lapply(export_manifests, function(x) list(
    export_operation_id = x$export_operation_id, manifest_sha256 = x$sha256
  ))
  order_key <- vapply(operations, `[[`, character(1), "export_operation_id")
  operations <- operations[order(order_key)]
  digest <- as.character(openssl::sha256(charToRaw(paste0(
    ph3_canonical_json(list(
      config_sha256 = config_digest, export_operations = operations
    )), "\n"
  ))))
  attributes(digest) <- NULL
  paste0("ph3-analysis-sha256:", digest)
}

ph3_classification_context <- function(
    prefix, containment, export_manifests
) {
  rows <- containment[containment$prefix == prefix, , drop = FALSE]
  if (nrow(rows) != 2L ||
      !setequal(rows$child_population_key, c("g1", "ph3_positive")) ||
      any(rows$containment_status != "validated") ||
      length(unique(rows$acquisition_id)) != 1L ||
      length(unique(rows$parent_row_count)) != 1L ||
      length(unique(rows$export_operation_id)) != 1L ||
      length(unique(rows$manifest_digest)) != 1L ||
      length(unique(rows$identity_method_id)) != 1L ||
      length(unique(rows$identity_method_version)) != 1L ||
      length(unique(rows$identity_source)) != 1L ||
      length(unique(rows$containment_method_id)) != 1L) {
    ph3_classification_fail(
      prefix, "invalid_containment_context",
      "exact validated G1 and pH3-positive containment is required"
    )
  }
  acquisition_id <- rows$acquisition_id[[1L]]
  operations <- Filter(function(x) {
    identical(x$export_operation_id, rows$export_operation_id[[1L]])
  }, export_manifests)
  if (length(operations) != 1L ||
      !identical(operations[[1L]]$sha256, rows$manifest_digest[[1L]])) {
    ph3_classification_fail(
      acquisition_id, "manifest_linkage_mismatch",
      "containment must link to exactly one verified export manifest"
    )
  }
  acquisition <- Filter(function(x) identical(x$acquisition_id, acquisition_id),
                        operations[[1L]]$manifest$acquisitions)
  if (length(acquisition) != 1L || !identical(acquisition[[1L]]$prefix, prefix)) {
    ph3_classification_fail(
      acquisition_id, "acquisition_linkage_mismatch",
      "the verified manifest acquisition and configured prefix disagree"
    )
  }
  list(rows = rows, operation = operations[[1L]], acquisition = acquisition[[1L]])
}

build_ph3_event_classification <- function(
    normalized, validated_inputs, config, manifest_row, containment,
    export_manifests
) {
  prefix <- manifest_row$prefix[[1L]]
  context <- ph3_classification_context(prefix, containment, export_manifests)
  acquisition_id <- context$acquisition$acquisition_id
  parent <- normalized$data
  g1 <- normalized$g1
  positive <- normalized$ph3_positive
  required_identity <- c(
    "acquisition_id", "event_index", "event_identity", "identity_source",
    "identity_method_id", "identity_method_version", "duplicate_occurrence",
    "export_profile", "export_operation_id", "export_manifest_digest"
  )
  validated_parent <- validated_inputs$complete
  validated_g1 <- validated_inputs$g1
  validated_positive <- validated_inputs$ph3_positive
  if (any(!required_identity %in% names(parent)) ||
      any(!required_identity %in% names(g1)) ||
      any(!required_identity %in% names(positive)) ||
      any(!required_identity %in% names(validated_parent)) ||
      any(!required_identity %in% names(validated_g1)) ||
      any(!required_identity %in% names(validated_positive)) ||
      any(!c("dna_norm", "target_raw", config$dna_channel) %in% names(parent)) ||
      anyNA(parent$event_identity) || any(!nzchar(parent$event_identity)) ||
      anyDuplicated(parent$event_identity) ||
      nrow(parent) != context$rows$parent_row_count[[1L]] ||
      !identical(parent$event_identity, validated_parent$event_identity) ||
      any(parent$acquisition_id != acquisition_id) ||
      any(parent$export_operation_id != context$rows$export_operation_id[[1L]]) ||
      any(parent$export_manifest_digest !=
          context$operation$manifest$manifest_binding$digest)) {
    ph3_classification_fail(
      acquisition_id, "parent_identity_inconsistency",
      "validated parent identity, row count, or immutable binding changed"
    )
  }
  validate_direct_identity <- function(table, population) {
    canonical_index <- !is.na(table$event_index) &
      grepl("^(0|[1-9][0-9]*)$", table$event_index)
    if (any(!canonical_index) ||
        any(table$event_identity != paste0(
          acquisition_id, ":event_index:", table$event_index
        )) ||
        any(table$identity_source != context$rows$identity_source[[1L]]) ||
        any(table$identity_method_id != context$rows$identity_method_id[[1L]]) ||
        any(table$identity_method_version !=
              context$rows$identity_method_version[[1L]]) ||
        any(table$export_profile != context$operation$manifest$profile) ||
        anyNA(table$duplicate_occurrence) ||
        any(table$duplicate_occurrence != "1")) {
      ph3_classification_fail(
        acquisition_id, paste0(population, "_direct_identity_mismatch"),
        "canonical direct event identity or its immutable provenance changed"
      )
    }
  }
  validate_direct_identity(parent, "parent")
  child_membership <- function(child, validated_child, population) {
    evidence <- context$rows[
      context$rows$child_population_key == population, , drop = FALSE
    ]
    if (anyNA(child$event_identity) || anyDuplicated(child$event_identity) ||
        any(!child$event_identity %in% parent$event_identity) ||
        !identical(child$event_identity, validated_child$event_identity) ||
        nrow(child) != evidence$child_row_count[[1L]] ||
        length(unique(child$event_identity)) !=
          evidence$child_unique_identity_count[[1L]] ||
        nrow(child) != evidence$matched_child_count[[1L]] ||
        any(child$acquisition_id != acquisition_id) ||
        any(child$export_operation_id !=
              context$rows$export_operation_id[[1L]]) ||
        any(child$export_manifest_digest !=
              context$operation$manifest$manifest_binding$digest) ||
        any(child$identity_method_id != context$rows$identity_method_id[[1L]]) ||
        any(child$identity_source != context$rows$identity_source[[1L]])) {
      ph3_classification_fail(
        acquisition_id, paste0(population, "_identity_inconsistency"),
        "child membership no longer matches validated direct identity"
      )
    }
    validate_direct_identity(child, population)
    parent$event_identity %in% child$event_identity
  }
  g1_member <- child_membership(g1, validated_g1, "g1")
  ph3_positive_member <- child_membership(
    positive, validated_positive, "ph3_positive"
  )

  gates <- ph3_gate_table(config)
  bounds <- c(gates$xmin[[1L]], gates$xmax)
  if (length(bounds) != 6L || any(!is.finite(bounds)) ||
      any(diff(bounds) <= 0)) {
    ph3_classification_fail(
      acquisition_id, "invalid_configured_boundaries",
      "b0 through b5 must be six finite strictly increasing boundaries"
    )
  }
  dna <- parent$dna_norm
  dna_finite <- is.finite(dna)
  identity_valid <- rep(TRUE, nrow(parent))
  eligible <- identity_valid & dna_finite & dna >= bounds[[1L]] &
    dna <= bounds[[6L]]
  sub_four <- eligible & dna < bounds[[5L]]
  four_n <- eligible & dna >= bounds[[5L]]
  phase <- as.character(assign_ph3_dna_phase(dna, gates))
  phase_assigned <- eligible & phase != "Unassigned"
  exclusion_reason <- ifelse(
    !identity_valid, "identity_invalid",
    ifelse(!dna_finite, "dna_nonfinite",
      ifelse(dna < bounds[[1L]], "below_b0",
        ifelse(dna > bounds[[6L]], "above_b5", "none")))
  )
  unassigned_reason <- ifelse(
    eligible & !phase_assigned, "configured_phase_gap", "none"
  )
  config_digest <- ph3_configuration_digest(config)
  positivity_method_id <- context$operation$manifest$approval$positivity_method_id
  geometry_status <- context$operation$manifest$geometry_overlay_status
  classification <- data.frame(
    classification_schema_version = "ph3-event-classification-1.0.0",
    analysis_id = ph3_analysis_id(config_digest, export_manifests),
    acquisition_id = acquisition_id,
    sample_id = context$acquisition$sample_id,
    prefix = prefix,
    parent_row_number = seq_len(nrow(parent)),
    event_index = parent$event_index,
    event_identity = parent$event_identity,
    event_identity_base = parent$event_identity,
    identity_occurrence = parent$duplicate_occurrence,
    identity_source = parent$identity_source,
    identity_method_id = parent$identity_method_id,
    identity_method_version = parent$identity_method_version,
    identity_valid = identity_valid,
    containment_method_id = context$rows$containment_method_id[[1L]],
    g1_member = g1_member,
    ph3_positive_member = ph3_positive_member,
    g1_containment_status = context$rows$containment_status[
      match("g1", context$rows$child_population_key)
    ],
    ph3_containment_status = context$rows$containment_status[
      match("ph3_positive", context$rows$child_population_key)
    ],
    dna_raw = parent[[config$dna_channel]],
    dna_norm = dna,
    dna_finite = dna_finite,
    b0 = bounds[[1L]], b1 = bounds[[2L]], b2 = bounds[[3L]],
    b3 = bounds[[4L]], b4 = bounds[[5L]], b5 = bounds[[6L]],
    eligible_2to4n = eligible,
    sub_4n_member = sub_four,
    four_n_member = four_n,
    configured_phase_id = phase,
    configured_phase_assigned = phase_assigned,
    eligibility_exclusion_reason = exclusion_reason,
    unassigned_reason = unassigned_reason,
    target_raw = parent$target_raw,
    target_display = parent$target_raw,
    target_display_finite = is.finite(parent$target_raw),
    display_transform_id = "raw_identity_v1",
    geometry_linkage_status = geometry_status,
    input_manifest_key = context$rows$manifest_digest[[1L]],
    input_manifest_reference = context$rows$manifest_reference[[1L]],
    config_digest = config_digest,
    export_profile = parent$export_profile,
    export_operation_id = context$rows$export_operation_id[[1L]],
    positivity_method_id = positivity_method_id,
    eligibility_method_id = "identity_validated_finite_dna_b0_b5_inclusive_v1",
    interval_method_id = "configured_shared_boundaries_left_closed_v1",
    four_n_method_id = "configured_b4_b5_closed_v1",
    sub_four_n_method_id =
      "configured_b0_b4_left_closed_right_open_v1",
    output_schema_version = "ph3-1.0.0",
    stringsAsFactors = FALSE
  )
  exclusion_levels <- c(
    "none", "identity_invalid", "dna_nonfinite", "below_b0", "above_b5"
  )
  if (nrow(classification) != nrow(parent) ||
      any(!classification$eligibility_exclusion_reason %in% exclusion_levels) ||
      sum(classification$eligible_2to4n) +
        sum(classification$eligibility_exclusion_reason != "none") != nrow(parent) ||
      any(classification$sub_4n_member & classification$four_n_member) ||
      sum(classification$sub_4n_member) + sum(classification$four_n_member) !=
        sum(classification$eligible_2to4n) ||
      sum(classification$configured_phase_assigned) +
        sum(classification$eligible_2to4n &
              classification$unassigned_reason != "none") !=
        sum(classification$eligible_2to4n) ||
      !identical(classification$g1_member, g1_member) ||
      !identical(classification$ph3_positive_member, ph3_positive_member)) {
    ph3_classification_fail(
      acquisition_id, "classification_reconciliation_failure",
      "row, exclusion, partition, phase, or membership counts do not reconcile"
    )
  }
  classification
}

#' Quantify user-gated pH3-positive events
#'
#' Every percentage uses all Single Cell events from the same sample as its
#' denominator. PH3-positive events outside the explicit DNA gates, including
#' nonfinite DNA values, are retained as `Unassigned`. This legacy quantitation
#' is unavailable for production direct-identity analyses.
#'
#' @param analysis A PH3-mode `facs_analysis` using the explicit legacy
#'   count-only input profile.
#' @return The updated analysis with `quantitation$ph3`.
#' @export
quantify_ph3 <- function(analysis) {
  validate_analysis_object(analysis)
  if (!identical(analysis$config$plot_type, "ph3")) {
    stop("`quantify_ph3()` requires a PH3-mode analysis.", call. = FALSE)
  }
  if (!identical(analysis$config$ph3_input_profile,
                 "legacy_count_only_unverified_v1")) {
    stop(
      "`quantify_ph3()` requires explicit `legacy_count_only_unverified_v1` input; production PH3 biological quantitation is withheld.",
      call. = FALSE
    )
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
