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
  ph3_positive_events <- if (is.null(ph3_positive_events)) {
    # Pilot profile has no imported pH3-positive population.  This empty table
    # is a structural placeholder only; positivity is computed later.
    events[FALSE, , drop = FALSE]
  } else {
    read_facs_sample(
      ph3_positive_events, dna_channel, target_channel,
      paste(sample_id, "pH3 positive"), 0L
    )
  }
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

ph3_metrics_fail <- function(acquisition_id, reason, detail) {
  stop(
    "PH3 acquisition metrics failed [", reason, "] for ", acquisition_id,
    ": ", detail, ".", call. = FALSE
  )
}

ph3_classification_scalar <- function(classification, field, acquisition_id) {
  value <- unique(classification[[field]])
  if (length(value) != 1L || is.na(value) ||
      (is.character(value) && !nzchar(value))) {
    ph3_metrics_fail(
      acquisition_id, "classification_provenance_inconsistency",
      paste0("`", field, "` must contain one nonmissing value")
    )
  }
  value[[1L]]
}

validate_ph3_metric_classification <- function(
    classification, manifest_row, config, export_manifests, containment
) {
  required <- c(
    "classification_schema_version", "analysis_id", "acquisition_id",
    "sample_id", "prefix", "event_identity", "identity_valid",
    "ph3_positive_member", "g1_containment_status",
    "ph3_containment_status", "dna_norm", "dna_finite", "b0", "b1",
    "b2", "b3", "b4", "b5", "eligible_2to4n", "sub_4n_member",
    "four_n_member", "configured_phase_id", "configured_phase_assigned",
    "eligibility_exclusion_reason", "unassigned_reason",
    "input_manifest_key", "config_digest", "export_operation_id",
    "positivity_method_id", "eligibility_method_id", "interval_method_id",
    "four_n_method_id", "sub_four_n_method_id", "containment_method_id",
    "output_schema_version"
  )
  acquisition_id <- if (is.data.frame(classification) &&
                          "acquisition_id" %in% names(classification) &&
                          nrow(classification)) {
    as.character(classification$acquisition_id[[1L]])
  } else {
    as.character(manifest_row$prefix[[1L]])
  }
  if (!is.data.frame(classification) || any(!required %in% names(classification))) {
    ph3_metrics_fail(
      acquisition_id, "missing_or_malformed_classification",
      "the retained Slice 3 classification is absent or lacks required fields"
    )
  }
  logical_fields <- c(
    "identity_valid", "ph3_positive_member", "dna_finite",
    "eligible_2to4n", "sub_4n_member", "four_n_member",
    "configured_phase_assigned"
  )
  if (!nrow(classification) || any(vapply(
      classification[logical_fields], function(x) !is.logical(x) || anyNA(x),
      logical(1)
    )) || anyNA(classification$event_identity) ||
      any(!nzchar(classification$event_identity)) ||
      anyDuplicated(classification$event_identity)) {
    ph3_metrics_fail(
      acquisition_id, "missing_or_malformed_classification",
      "classification rows, identities, or logical membership fields are invalid"
    )
  }
  scalar_fields <- c(
    "classification_schema_version", "analysis_id", "acquisition_id",
    "sample_id", "prefix", "g1_containment_status",
    "ph3_containment_status", paste0("b", 0:5), "input_manifest_key",
    "config_digest", "export_operation_id", "positivity_method_id",
    "eligibility_method_id", "interval_method_id", "four_n_method_id",
    "sub_four_n_method_id", "containment_method_id", "output_schema_version"
  )
  values <- stats::setNames(lapply(
    scalar_fields,
    function(field) ph3_classification_scalar(classification, field, acquisition_id)
  ), scalar_fields)
  bounds <- as.numeric(unlist(values[paste0("b", 0:5)]))
  active_gates <- ph3_gate_table(config)
  active_bounds <- c(active_gates$xmin[[1L]], active_gates$xmax)
  phases <- c("G1", "Early S", "Mid S", "Late S", "G2/M")
  exclusion_levels <- c(
    "none", "identity_invalid", "dna_nonfinite", "below_b0", "above_b5"
  )
  eligible <- classification$eligible_2to4n
  active_config_digest <- ph3_configuration_digest(config)
  active_operations <- Filter(function(operation) {
    identical(operation$export_operation_id, values$export_operation_id)
  }, export_manifests)
  active_acquisitions <- if (length(active_operations) == 1L) {
    Filter(function(acquisition) {
      identical(acquisition$acquisition_id, values$acquisition_id)
    }, active_operations[[1L]]$manifest$acquisitions)
  } else {
    list()
  }
  containment_fields <- c(
    "acquisition_id", "prefix", "child_population_key", "containment_status",
    "containment_method_id", "export_operation_id", "manifest_digest"
  )
  active_containment <- if (is.data.frame(containment) &&
                              all(containment_fields %in% names(containment))) {
    containment[
      containment$acquisition_id == values$acquisition_id, , drop = FALSE
    ]
  } else {
    data.frame()
  }
  computed_cutoff <- identical(values$positivity_method_id,
                               "ph3_raw_4n_density_cutoff_v1")
  configured_computed_cutoff <- identical(config$ph3_positivity_method,
                                           "ph3_raw_4n_density_cutoff_v1")
  computed_membership_valid <- TRUE
  if (computed_cutoff) {
    source_method_id <- if (length(active_operations) == 1L) {
      active_operations[[1L]]$manifest$approval$positivity_method_id
    } else NA_character_
    required_computed <- c("flowjo_ph3_positive_member", "positivity_call_status",
                           "positivity_call_reason_code", "raw_4n_cutoff", "target_raw")
    cutoff_failure <- all(classification$positivity_call_status == "unavailable_cutoff_failure")
    computed_membership_valid <- all(required_computed %in% names(classification)) &&
      is.logical(classification$flowjo_ph3_positive_member) &&
      !anyNA(classification$flowjo_ph3_positive_member) &&
      "flowjo_positivity_method_id" %in% names(classification) &&
      identical(ph3_classification_scalar(classification,
                                          "flowjo_positivity_method_id",
                                          acquisition_id),
                source_method_id) &&
      is.character(classification$positivity_call_status) &&
      is.character(classification$positivity_call_reason_code) &&
      is.numeric(classification$raw_4n_cutoff) &&
      length(unique(classification$raw_4n_cutoff)) == 1L &&
      (cutoff_failure && all(is.na(classification$raw_4n_cutoff)) &&
       all(!classification$ph3_positive_member) &&
       all(!is.na(classification$positivity_call_reason_code)) ||
       !cutoff_failure && is.finite(classification$raw_4n_cutoff[[1L]]) &&
       identical(classification$ph3_positive_member,
                 is.finite(classification$target_raw) &
                   classification$target_raw > classification$raw_4n_cutoff[[1L]]) &&
       identical(classification$positivity_call_status,
                 ifelse(is.finite(classification$target_raw), "called",
                        "unavailable_nonfinite_raw_signal")) &&
       identical(classification$positivity_call_reason_code,
                 ifelse(is.finite(classification$target_raw), NA_character_,
                        "nonfinite_raw_signal")))
  }
  if (length(active_operations) != 1L ||
      length(active_acquisitions) != 1L ||
      !identical(active_acquisitions[[1L]]$sample_id, values$sample_id) ||
      !identical(active_acquisitions[[1L]]$prefix, values$prefix) ||
      !identical(active_operations[[1L]]$sha256, values$input_manifest_key) ||
      (!computed_cutoff && !identical(active_operations[[1L]]$manifest$approval$positivity_method_id,
                                      values$positivity_method_id)) ||
      !identical(active_config_digest, values$config_digest) ||
      !identical(ph3_analysis_id(active_config_digest, export_manifests),
                 values$analysis_id) ||
      !identical(values$containment_method_id,
                 "exact_direct_identity_multiset_containment") ||
      nrow(active_containment) != 2L ||
      !setequal(active_containment$child_population_key,
                c("g1", "ph3_positive")) ||
      any(active_containment$containment_status != "validated") ||
      any(active_containment$containment_method_id !=
            values$containment_method_id) ||
      any(active_containment$prefix != values$prefix) ||
      any(active_containment$export_operation_id !=
            values$export_operation_id) ||
      any(active_containment$manifest_digest != values$input_manifest_key) ||
      !computed_membership_valid ||
      !identical(computed_cutoff, configured_computed_cutoff)) {
    ph3_metrics_fail(
      acquisition_id, "active_provenance_mismatch",
      "classification configuration, manifest, operation, positivity, or analysis binding does not match the active analysis"
    )
  }
  expected_eligible <- classification$identity_valid &
    classification$dna_finite & classification$dna_norm >= bounds[[1L]] &
    classification$dna_norm <= bounds[[6L]]
  expected_sub_four <- expected_eligible &
    classification$dna_norm < bounds[[5L]]
  expected_four_n <- expected_eligible &
    classification$dna_norm >= bounds[[5L]]
  phase_gates <- data.frame(
    gate = phases, gate_index = seq_along(phases),
    xmin = bounds[1:5], xmax = bounds[2:6], stringsAsFactors = FALSE
  )
  expected_phase <- as.character(assign_ph3_dna_phase(
    classification$dna_norm, phase_gates
  ))
  expected_phase_assigned <- expected_eligible & expected_phase != "Unassigned"
  expected_exclusion <- ifelse(
    !classification$identity_valid, "identity_invalid",
    ifelse(!classification$dna_finite, "dna_nonfinite",
      ifelse(classification$dna_norm < bounds[[1L]], "below_b0",
        ifelse(classification$dna_norm > bounds[[6L]], "above_b5", "none")))
  )
  expected_unassigned <- ifelse(
    expected_eligible & !expected_phase_assigned, "configured_phase_gap", "none"
  )
  if (!identical(values$classification_schema_version,
                 "ph3-event-classification-1.0.0") ||
      !identical(values$output_schema_version, "ph3-1.0.0") ||
      !identical(values$eligibility_method_id,
                 "identity_validated_finite_dna_b0_b5_inclusive_v1") ||
      !identical(values$interval_method_id,
                 "configured_shared_boundaries_left_closed_v1") ||
      !identical(values$four_n_method_id,
                 "configured_b4_b5_closed_v1") ||
      !identical(values$sub_four_n_method_id,
                 "configured_b0_b4_left_closed_right_open_v1") ||
      !identical(values$g1_containment_status, "validated") ||
      !identical(values$ph3_containment_status, "validated") ||
      !identical(values$prefix, as.character(manifest_row$prefix[[1L]])) ||
      any(!is.finite(bounds)) || any(diff(bounds) <= 0) ||
      !identical(bounds, active_bounds) ||
      !all(classification$identity_valid) ||
      !identical(classification$dna_finite,
                 is.finite(classification$dna_norm)) ||
      any(!classification$eligibility_exclusion_reason %in% exclusion_levels) ||
      any(!classification$unassigned_reason %in%
            c("none", "configured_phase_gap")) ||
      !identical(classification$eligible_2to4n, expected_eligible) ||
      !identical(classification$sub_4n_member, expected_sub_four) ||
      !identical(classification$four_n_member, expected_four_n) ||
      !identical(as.character(classification$configured_phase_id),
                 expected_phase) ||
      !identical(classification$configured_phase_assigned,
                 expected_phase_assigned) ||
      !identical(classification$eligibility_exclusion_reason,
                 expected_exclusion) ||
      !identical(classification$unassigned_reason, expected_unassigned) ||
      sum(eligible) +
        sum(classification$eligibility_exclusion_reason != "none") !=
        nrow(classification) ||
      any(eligible !=
        (classification$eligibility_exclusion_reason == "none")) ||
      any(classification$sub_4n_member & classification$four_n_member) ||
      any((classification$sub_4n_member | classification$four_n_member) &
        !eligible) ||
      sum(classification$sub_4n_member) +
        sum(classification$four_n_member) != sum(eligible) ||
      any(classification$configured_phase_assigned !=
        (eligible & classification$configured_phase_id %in% phases)) ||
      any(eligible & !classification$configured_phase_assigned &
        classification$unassigned_reason == "none") ||
      any((!eligible | classification$configured_phase_assigned) &
        classification$unassigned_reason != "none")) {
    ph3_metrics_fail(
      acquisition_id, "classification_reconciliation_failure",
      "Slice 3 provenance, eligibility, partition, phase, or reason fields do not reconcile"
    )
  }
  list(acquisition_id = acquisition_id, values = values, bounds = bounds,
       phases = phases)
}

ph3_percentage_result <- function(numerator, denominator) {
  numerator <- as.integer(numerator)
  denominator <- as.integer(denominator)
  if (denominator == 0L) {
    return(list(value = NA_real_, status = "undefined_zero_denominator"))
  }
  value <- 100 * numerator / denominator
  if (!is.finite(value)) {
    stop("PH3 percentage calculation produced a nonfinite value.", call. = FALSE)
  }
  list(value = value, status = "ok")
}

ph3_metric_common <- function(manifest_row, context) {
  values <- context$values
  data.frame(
    analysis_id = values$analysis_id,
    acquisition_id = values$acquisition_id,
    sample_id = values$sample_id,
    prefix = as.character(manifest_row$prefix[[1L]]),
    condition = as.character(manifest_row$condition[[1L]]),
    condition_index = as.integer(manifest_row$condition_index[[1L]]),
    replicate = as.character(manifest_row$replicate[[1L]]),
    replicate_index = as.integer(manifest_row$replicate_index[[1L]]),
    technical_replicate = as.character(manifest_row$technical_replicate[[1L]]),
    aggregation_level = "acquisition",
    classification_schema_version = values$classification_schema_version,
    output_schema_version = values$output_schema_version,
    positivity_method_id = values$positivity_method_id,
    eligibility_method_id = values$eligibility_method_id,
    interval_method_id = values$interval_method_id,
    four_n_method_id = values$four_n_method_id,
    sub_four_n_method_id = values$sub_four_n_method_id,
    containment_method_id = values$containment_method_id,
    config_digest = values$config_digest,
    export_operation_id = values$export_operation_id,
    input_manifest_key = values$input_manifest_key,
    stringsAsFactors = FALSE
  )
}

derive_ph3_acquisition_tables <- function(
    classification, manifest_row, config, export_manifests, containment
) {
  context <- validate_ph3_metric_classification(
    classification, manifest_row, config, export_manifests, containment
  )
  common <- ph3_metric_common(manifest_row, context)
  cutoff_unavailable <- "positivity_call_status" %in% names(classification) &&
    all(classification$positivity_call_status == "unavailable_cutoff_failure")
  b <- context$bounds
  eligible <- classification$eligible_2to4n
  positive <- classification$ph3_positive_member
  sub_four <- classification$sub_4n_member
  four_n <- classification$four_n_member
  a_n <- as.integer(sum(eligible))
  s_n <- as.integer(sum(sub_four))
  f_n <- as.integer(sum(four_n))
  pa_n <- as.integer(sum(positive & eligible))
  ps_n <- as.integer(sum(positive & sub_four))
  pf_n <- as.integer(sum(positive & four_n))
  specifications <- list(
    list("ph3_2to4n_positivity_percent", pa_n, a_n, "eligible_2to4n",
         b[[1L]], b[[6L]], TRUE, TRUE),
    list("ph3_within_4n_positivity_percent", pf_n, f_n, "eligible_4n",
         b[[5L]], b[[6L]], TRUE, TRUE),
    list("ph3_4n_positive_prevalence_within_2to4n_percent", pf_n, a_n,
         "eligible_4n_positive_within_2to4n", b[[5L]], b[[6L]], TRUE, TRUE),
    list("ph3_within_sub_4n_positivity_percent", ps_n, s_n,
         "eligible_sub_4n", b[[1L]], b[[5L]], TRUE, FALSE),
    list("ph3_sub_4n_positive_prevalence_within_2to4n_percent", ps_n, a_n,
         "eligible_sub_4n_positive_within_2to4n", b[[1L]], b[[5L]], TRUE, FALSE)
  )
  metric_rows <- lapply(specifications, function(specification) {
    result <- if (cutoff_unavailable) {
      list(value = NA_real_, status = "unavailable_cutoff_failure")
    } else ph3_percentage_result(specification[[2L]], specification[[3L]])
    cbind(
      common,
      data.frame(
        metric_id = specification[[1L]], value_percent = result$value,
        numerator_count = as.integer(specification[[2L]]),
        denominator_count = as.integer(specification[[3L]]),
        population_id = specification[[4L]],
        interval_lower = as.numeric(specification[[5L]]),
        interval_upper = as.numeric(specification[[6L]]),
        lower_inclusive = specification[[7L]],
        upper_inclusive = specification[[8L]], result_status = result$status,
        stringsAsFactors = FALSE
      )
    )
  })
  metrics <- do.call(rbind, metric_rows)
  rownames(metrics) <- NULL

  phase_bounds <- data.frame(
    phase_id = context$phases, phase_index = seq_along(context$phases),
    interval_lower = b[1:5], interval_upper = b[2:6],
    lower_inclusive = TRUE,
    upper_inclusive = c(FALSE, FALSE, FALSE, FALSE, TRUE),
    stringsAsFactors = FALSE
  )
  phase_rows <- lapply(seq_len(nrow(phase_bounds)), function(i) {
    phase <- phase_bounds$phase_id[[i]]
    numerator <- as.integer(sum(
      eligible & positive & classification$configured_phase_id == phase
    ))
    result <- if (cutoff_unavailable) {
      list(value = NA_real_, status = "unavailable_cutoff_failure")
    } else ph3_percentage_result(numerator, a_n)
    cbind(common, phase_bounds[i, , drop = FALSE], data.frame(
      metric_id = "ph3_phase_positive_prevalence_within_2to4n_percent",
      ph3_positive_count = numerator, eligible_2to4n_count = a_n,
      value_percent = result$value, result_status = result$status,
      stringsAsFactors = FALSE
    ))
  })
  phase_prevalence <- do.call(rbind, phase_rows)
  rownames(phase_prevalence) <- NULL

  qc_specifications <- list(
    list("rows", "imported_classification_rows", rep(TRUE, nrow(classification))),
    list("identity", "identity_valid", classification$identity_valid),
    list("eligibility", "eligible_2to4n", eligible),
    list("eligibility_exclusion", "identity_invalid",
         classification$eligibility_exclusion_reason == "identity_invalid"),
    list("eligibility_exclusion", "dna_nonfinite",
         classification$eligibility_exclusion_reason == "dna_nonfinite"),
    list("eligibility_exclusion", "below_b0",
         classification$eligibility_exclusion_reason == "below_b0"),
    list("eligibility_exclusion", "above_b5",
         classification$eligibility_exclusion_reason == "above_b5"),
    list("partition", "sub_4n", sub_four),
    list("partition", "4n", four_n),
    list("configured_phase", "configured_phase_assigned",
         classification$configured_phase_assigned)
  )
  for (phase in context$phases) {
    qc_specifications[[length(qc_specifications) + 1L]] <- list(
      "configured_phase", phase,
      eligible & classification$configured_phase_id == phase
    )
  }
  unassigned_reasons <- sort(unique(
    classification$unassigned_reason[
      eligible & classification$unassigned_reason != "none"
    ]
  ))
  if (!length(unassigned_reasons)) unassigned_reasons <- "none"
  for (reason in unassigned_reasons) {
    selected <- eligible & !classification$configured_phase_assigned
    if (!identical(reason, "none")) {
      selected <- selected & classification$unassigned_reason == reason
    }
    qc_specifications[[length(qc_specifications) + 1L]] <- list(
      "eligible_configured_phase_unassigned", reason, selected
    )
  }
  qc_rows <- lapply(qc_specifications, function(specification) {
    selected <- specification[[3L]]
    cbind(common, data.frame(
      qc_dimension = specification[[1L]], reason = specification[[2L]],
      event_count = as.integer(sum(selected)),
      ph3_positive_count = as.integer(sum(selected & positive)),
      imported_classification_count = as.integer(nrow(classification)),
      identity_valid_count = as.integer(sum(classification$identity_valid)),
      eligible_2to4n_count = a_n, eligible_sub_4n_count = s_n,
      eligible_4n_count = f_n,
      configured_phase_assigned_count = as.integer(sum(
        classification$configured_phase_assigned
      )),
      eligibility_excluded_count = as.integer(nrow(classification) - a_n),
      eligible_configured_phase_unassigned_count = as.integer(sum(
        eligible & !classification$configured_phase_assigned
      )), eligibility_reconciliation_difference = as.integer(
        nrow(classification) - a_n -
          sum(classification$eligibility_exclusion_reason != "none")
      ), partition_reconciliation_difference = as.integer(a_n - s_n - f_n),
      configured_phase_reconciliation_difference = as.integer(
        a_n - sum(classification$configured_phase_assigned) -
          sum(eligible & !classification$configured_phase_assigned)
      ), qc_status = "ok", stringsAsFactors = FALSE
    ))
  })
  eligibility_qc <- do.call(rbind, qc_rows)
  rownames(eligibility_qc) <- NULL

  if (!config_scalar_number(config$ph3_boundary_sensitivity_fraction) ||
      !config_scalar_number(config$dna_2n_value)) {
    ph3_metrics_fail(
      context$acquisition_id, "invalid_boundary_sensitivity_interval",
      "delta source parameters must be present finite scalar numbers"
    )
  }
  delta <- config$ph3_boundary_sensitivity_fraction * config$dna_2n_value
  variants <- data.frame(
    perturbation_id = c("primary", "lower_outward", "lower_inward", "upper_inward"),
    perturbation_status = c("primary", rep("perturbation", 3L)),
    interval_lower = c(b[[5L]], b[[5L]] - delta, b[[5L]] + delta, b[[5L]]),
    interval_upper = c(b[[6L]], b[[6L]], b[[6L]], b[[6L]] - delta),
    lower_inclusive = TRUE, upper_inclusive = TRUE,
    stringsAsFactors = FALSE
  )
  if (!is.finite(delta) || delta <= 0 ||
      any(!is.finite(variants$interval_lower)) ||
      any(!is.finite(variants$interval_upper)) ||
      any(variants$interval_lower >= variants$interval_upper)) {
    ph3_metrics_fail(
      context$acquisition_id, "invalid_boundary_sensitivity_interval",
      "delta and all four closed perturbation intervals must be finite, positive, and noninverted"
    )
  }
  sensitivity_rows <- lapply(seq_len(nrow(variants)), function(i) {
    region <- eligible & classification$dna_norm >= variants$interval_lower[[i]] &
      classification$dna_norm <= variants$interval_upper[[i]]
    regional_n <- as.integer(sum(region))
    regional_positive_n <- as.integer(sum(region & positive))
    within <- ph3_percentage_result(regional_positive_n, regional_n)
    prevalence <- ph3_percentage_result(regional_positive_n, a_n)
    cbind(common, variants[i, , drop = FALSE], data.frame(
      delta = as.numeric(delta),
      ph3_boundary_sensitivity_fraction =
        as.numeric(config$ph3_boundary_sensitivity_fraction),
      dna_2n_value = as.numeric(config$dna_2n_value),
      eligible_ph3_positive_regional_count = regional_positive_n,
      eligible_regional_count = regional_n, eligible_2to4n_count = a_n,
      within_region_positivity_percent = within$value,
      within_region_result_status = within$status,
      positive_prevalence_within_2to4n_percent = prevalence$value,
      prevalence_result_status = prevalence$status,
      stringsAsFactors = FALSE
    ))
  })
  sensitivity <- do.call(rbind, sensitivity_rows)
  rownames(sensitivity) <- NULL
  primary_metric <- metrics[
    metrics$metric_id == "ph3_within_4n_positivity_percent", , drop = FALSE
  ]
  primary_prevalence <- metrics[
    metrics$metric_id ==
      "ph3_4n_positive_prevalence_within_2to4n_percent", , drop = FALSE
  ]
  if (nrow(primary_metric) != 1L || nrow(primary_prevalence) != 1L ||
      sensitivity$eligible_ph3_positive_regional_count[[1L]] !=
        primary_metric$numerator_count[[1L]] ||
      sensitivity$eligible_regional_count[[1L]] !=
        primary_metric$denominator_count[[1L]] ||
      !identical(sensitivity$within_region_positivity_percent[[1L]],
                 primary_metric$value_percent[[1L]]) ||
      !identical(sensitivity$positive_prevalence_within_2to4n_percent[[1L]],
                 primary_prevalence$value_percent[[1L]])) {
    ph3_metrics_fail(
      context$acquisition_id, "primary_sensitivity_invariance_failure",
      "the primary sensitivity row does not reproduce the approved 4N metrics"
    )
  }
  list(
    ph3_metrics_acquisition = metrics,
    ph3_phase_prevalence = phase_prevalence,
    ph3_event_eligibility_qc = eligibility_qc,
    ph3_4n_boundary_sensitivity_qc = sensitivity
  )
}

ph3_aggregation_fail <- function(reason, detail) {
  stop(
    "PH3 replicate/condition aggregation failed [", reason, "]: ", detail,
    ".", call. = FALSE
  )
}

ph3_canonical_string_array <- function(value) {
  ph3_canonical_json(as.list(sort(unique(as.character(value)))))
}

ph3_require_one_group_value <- function(data, field) {
  value <- unique(data[[field]])
  if (length(value) != 1L || is.na(value) ||
      (is.character(value) && !nzchar(value))) {
    ph3_aggregation_fail(
      "mixed_or_missing_provenance",
      paste0("`", field, "` must be one nonmissing value per aggregate group")
    )
  }
  value[[1L]]
}

ph3_validate_slice4_aggregation_inputs <- function(
    quantitation, manifest, analysis_provenance
) {
  table_names <- c(
    "ph3_metrics_acquisition", "ph3_phase_prevalence",
    "ph3_event_eligibility_qc", "ph3_4n_boundary_sensitivity_qc"
  )
  if (!is.data.frame(manifest) || !nrow(manifest) ||
      any(!c("prefix", "condition", "condition_index", "replicate",
             "replicate_index", "technical_replicate") %in% names(manifest)) ||
      anyNA(manifest[c("prefix", "condition", "condition_index", "replicate",
                       "replicate_index", "technical_replicate")]) ||
      any(!vapply(manifest[c("prefix", "condition", "replicate",
                             "technical_replicate")], function(value) {
        all(nzchar(as.character(value)))
      }, logical(1))) ||
      anyDuplicated(manifest$prefix) ||
      anyDuplicated(manifest[c("condition", "replicate", "technical_replicate")])) {
    ph3_aggregation_fail(
      "invalid_explicit_manifest",
      "the active manifest must declare unique acquisition prefixes and condition/replicate/technical membership"
    )
  }
  if (any(!table_names %in% names(quantitation)) ||
      any(!vapply(quantitation[table_names], is.data.frame, logical(1)))) {
    ph3_aggregation_fail(
      "missing_slice4_table", "all four authoritative Slice 4 tables are required"
    )
  }
  common <- c(
    "analysis_id", "acquisition_id", "sample_id", "prefix", "condition",
    "condition_index", "replicate", "replicate_index",
    "technical_replicate", "aggregation_level",
    "classification_schema_version", "output_schema_version",
    "positivity_method_id", "eligibility_method_id", "interval_method_id",
    "four_n_method_id", "sub_four_n_method_id", "containment_method_id",
    "config_digest", "export_operation_id", "input_manifest_key"
  )
  expected_rows <- c(5L, 5L, 16L, 4L)
  source_acquisitions <- NULL
  source_mapping <- NULL
  condition_mapping <- unique(manifest[c("condition", "condition_index")])
  replicate_mapping <- unique(manifest[c(
    "condition", "replicate", "replicate_index"
  )])
  if (anyDuplicated(condition_mapping$condition) ||
      anyDuplicated(condition_mapping$condition_index) ||
      anyDuplicated(replicate_mapping[c("condition", "replicate")]) ||
      anyDuplicated(replicate_mapping[c("condition", "replicate_index")])) {
    ph3_aggregation_fail(
      "unstable_manifest_group_index",
      "conditions require global index bijection and replicates require condition-scoped index bijection"
    )
  }
  operations <- analysis_provenance$ph3_export_manifests
  containment <- analysis_provenance$ph3_containment
  if (!is.list(operations) || !length(operations) ||
      !is.data.frame(containment)) {
    ph3_aggregation_fail(
      "missing_active_analysis_provenance",
      "active export-manifest and containment provenance are required"
    )
  }
  active_rows <- list()
  active_index <- 0L
  for (operation in operations) {
    acquisitions <- operation$manifest$acquisitions
    if (!is.list(acquisitions) || !length(acquisitions) ||
        !is.character(operation$sha256) || length(operation$sha256) != 1L ||
        is.na(operation$sha256) || !nzchar(operation$sha256) ||
        !is.character(operation$export_operation_id) ||
        length(operation$export_operation_id) != 1L ||
        is.na(operation$export_operation_id) ||
        !nzchar(operation$export_operation_id)) {
      ph3_aggregation_fail(
        "missing_active_analysis_provenance",
        "active export operations lack acquisition or digest identity"
      )
    }
    for (acquisition in acquisitions) {
      required <- c("acquisition_id", "sample_id", "prefix")
      required_values <- if (all(required %in% names(acquisition))) {
        vapply(acquisition[required], as.character, character(1))
      } else {
        character()
      }
      if (length(required_values) != length(required) || anyNA(required_values) ||
          any(!nzchar(required_values))) {
        ph3_aggregation_fail(
          "missing_active_analysis_provenance",
          "an active manifest acquisition lacks required identity"
        )
      }
      active_index <- active_index + 1L
      active_rows[[active_index]] <- data.frame(
        acquisition_id = acquisition$acquisition_id,
        sample_id = acquisition$sample_id, prefix = acquisition$prefix,
        export_operation_id = operation$export_operation_id,
        input_manifest_key = operation$sha256, stringsAsFactors = FALSE
      )
    }
  }
  active_mapping <- do.call(rbind, active_rows)
  if (nrow(active_mapping) != nrow(manifest) ||
      anyDuplicated(active_mapping$acquisition_id) ||
      anyDuplicated(active_mapping$prefix) ||
      !setequal(active_mapping$prefix, manifest$prefix)) {
    ph3_aggregation_fail(
      "active_provenance_membership_mismatch",
      "active manifest acquisitions do not match the explicit sample manifest"
    )
  }
  containment_required <- c(
    "acquisition_id", "prefix", "child_population_key", "containment_status",
    "export_operation_id", "manifest_digest"
  )
  if (any(!containment_required %in% names(containment))) {
    ph3_aggregation_fail(
      "missing_active_analysis_provenance",
      "active containment provenance lacks required identity fields"
    )
  }
  if (nrow(containment) != 2L * nrow(active_mapping) ||
      !setequal(unique(containment$acquisition_id), active_mapping$acquisition_id)) {
    ph3_aggregation_fail(
      "active_provenance_membership_mismatch",
      "active containment acquisition membership is not exact"
    )
  }
  for (i in seq_along(table_names)) {
    name <- table_names[[i]]
    data <- quantitation[[name]]
    string_common <- c(
      "analysis_id", "acquisition_id", "sample_id", "prefix", "condition",
      "replicate", "technical_replicate", "aggregation_level",
      "classification_schema_version", "output_schema_version",
      "positivity_method_id", "eligibility_method_id", "interval_method_id",
      "four_n_method_id", "sub_four_n_method_id", "containment_method_id",
      "config_digest", "export_operation_id", "input_manifest_key"
    )
    if (any(!common %in% names(data)) || nrow(data) != nrow(manifest) * expected_rows[[i]] ||
        anyNA(data[common]) || any(data$aggregation_level != "acquisition") ||
        any(data$output_schema_version != "ph3-1.0.0") ||
        any(!vapply(data[string_common], function(value) {
          all(nzchar(as.character(value)))
        }, logical(1)))) {
      ph3_aggregation_fail(
        "malformed_or_incomplete_slice4_table",
        paste0("`", name, "` has missing fields, rows, values, or schema")
      )
    }
    keys <- if (identical(name, "ph3_metrics_acquisition")) {
      c("acquisition_id", "metric_id")
    } else if (identical(name, "ph3_phase_prevalence")) {
      c("acquisition_id", "phase_id")
    } else if (identical(name, "ph3_event_eligibility_qc")) {
      c("acquisition_id", "qc_dimension", "reason")
    } else {
      c("acquisition_id", "perturbation_id")
    }
    if (any(!keys %in% names(data)) || anyDuplicated(data[keys])) {
      ph3_aggregation_fail(
        "duplicated_or_malformed_slice4_row",
        paste0("`", name, "` does not have unique required rows")
      )
    }
    current_acquisitions <- sort(unique(data$acquisition_id))
    if (is.null(source_acquisitions)) source_acquisitions <- current_acquisitions
    if (!identical(current_acquisitions, source_acquisitions)) {
      ph3_aggregation_fail(
        "unreconciled_acquisition_membership",
        "Slice 4 tables do not contain the same acquisitions"
      )
    }
    current_mapping <- unique(data[c("prefix", "acquisition_id")])
    current_mapping <- current_mapping[order(current_mapping$prefix), , drop = FALSE]
    rownames(current_mapping) <- NULL
    if (is.null(source_mapping)) source_mapping <- current_mapping
    if (!identical(current_mapping, source_mapping)) {
      ph3_aggregation_fail(
        "unreconciled_acquisition_membership",
        "Slice 4 tables disagree on the manifest-prefix/acquisition mapping"
      )
    }
    for (j in seq_len(nrow(manifest))) {
      rows <- data[data$prefix == manifest$prefix[[j]], , drop = FALSE]
      if (nrow(rows) != expected_rows[[i]] ||
          length(unique(rows$acquisition_id)) != 1L ||
          any(rows$condition != manifest$condition[[j]]) ||
          any(rows$condition_index != manifest$condition_index[[j]]) ||
          any(rows$replicate != manifest$replicate[[j]]) ||
          any(rows$replicate_index != manifest$replicate_index[[j]]) ||
          any(rows$technical_replicate != manifest$technical_replicate[[j]])) {
        ph3_aggregation_fail(
          "manifest_membership_mismatch",
          paste0("`", name, "` does not exactly match manifest prefix `",
                 manifest$prefix[[j]], "`")
        )
      }
    }
  }
  metric_order <- c(
    "ph3_2to4n_positivity_percent", "ph3_within_4n_positivity_percent",
    "ph3_4n_positive_prevalence_within_2to4n_percent",
    "ph3_within_sub_4n_positivity_percent",
    "ph3_sub_4n_positive_prevalence_within_2to4n_percent"
  )
  phase_order <- c("G1", "Early S", "Mid S", "Late S", "G2/M")
  metrics <- quantitation$ph3_metrics_acquisition
  phases <- quantitation$ph3_phase_prevalence
  for (acquisition in source_acquisitions) {
    if (!identical(metrics$metric_id[metrics$acquisition_id == acquisition],
                   metric_order) ||
        !identical(phases$phase_id[phases$acquisition_id == acquisition],
                   phase_order) ||
        !identical(phases$phase_index[phases$acquisition_id == acquisition],
                   seq_along(phase_order))) {
      ph3_aggregation_fail(
        "missing_extra_or_reordered_required_row",
        "each acquisition must contain the exact fixed metric and phase rows"
      )
    }
    sensitivity_ids <- quantitation$ph3_4n_boundary_sensitivity_qc$perturbation_id[
      quantitation$ph3_4n_boundary_sensitivity_qc$acquisition_id == acquisition
    ]
    if (!identical(sensitivity_ids,
                   c("primary", "lower_outward", "lower_inward", "upper_inward"))) {
      ph3_aggregation_fail(
        "missing_extra_or_reordered_required_row",
        "each acquisition must retain the exact four Slice 4 sensitivity rows"
      )
    }
    qc <- quantitation$ph3_event_eligibility_qc[
      quantitation$ph3_event_eligibility_qc$acquisition_id == acquisition,
      , drop = FALSE
    ]
    expected_qc_dimension <- c(
      "rows", "identity", "eligibility", rep("eligibility_exclusion", 4L),
      rep("partition", 2L), rep("configured_phase", 6L),
      "eligible_configured_phase_unassigned"
    )
    expected_qc_reason <- c(
      "imported_classification_rows", "identity_valid", "eligible_2to4n",
      "identity_invalid", "dna_nonfinite", "below_b0", "above_b5",
      "sub_4n", "4n", "configured_phase_assigned", "G1", "Early S",
      "Mid S", "Late S", "G2/M"
    )
    if (!identical(qc$qc_dimension, expected_qc_dimension) ||
        !identical(qc$reason[seq_along(expected_qc_reason)], expected_qc_reason) ||
        !qc$reason[[16L]] %in% c("none", "configured_phase_gap")) {
      ph3_aggregation_fail(
        "invalid_slice4_eligibility_qc_sequence",
        "each acquisition must retain the exact Slice 4 eligibility QC rows"
      )
    }
    active <- active_mapping[
      active_mapping$acquisition_id == acquisition, , drop = FALSE
    ]
    containment_rows <- containment[
      containment$acquisition_id == acquisition, , drop = FALSE
    ]
    if (nrow(active) != 1L ||
        nrow(containment_rows) != 2L ||
        !setequal(containment_rows$child_population_key, c("g1", "ph3_positive")) ||
        any(containment_rows$containment_status != "validated") ||
        any(containment_rows$prefix != active$prefix[[1L]]) ||
        any(containment_rows$export_operation_id != active$export_operation_id[[1L]]) ||
        any(containment_rows$manifest_digest != active$input_manifest_key[[1L]])) {
      ph3_aggregation_fail(
        "active_provenance_membership_mismatch",
        paste0("active containment provenance does not bind acquisition `",
               acquisition, "`")
      )
    }
    source_identity <- unique(metrics[
      metrics$acquisition_id == acquisition,
      c("acquisition_id", "sample_id", "prefix", "export_operation_id",
        "input_manifest_key"), drop = FALSE
    ])
    rownames(source_identity) <- NULL
    rownames(active) <- NULL
    if (!identical(source_identity, active)) {
      ph3_aggregation_fail(
        "active_provenance_membership_mismatch",
        paste0("Slice 4 identity does not match active provenance for `",
               acquisition, "`")
      )
    }
    for (field in common) {
      reference <- unique(metrics[[field]][metrics$acquisition_id == acquisition])
      if (length(reference) != 1L) {
        ph3_aggregation_fail(
          "mixed_or_missing_provenance",
          paste0("acquisition `", acquisition, "` has inconsistent `", field, "`")
        )
      }
      for (name in table_names) {
        candidate <- unique(quantitation[[name]][[field]][
          quantitation[[name]]$acquisition_id == acquisition
        ])
        if (!identical(candidate, reference)) {
          ph3_aggregation_fail(
            "mixed_or_missing_provenance",
            paste0("Slice 4 tables disagree on acquisition `", acquisition,
                   "` field `", field, "`")
          )
        }
      }
    }
  }
  provenance <- c(
    "analysis_id", "classification_schema_version", "output_schema_version",
    "positivity_method_id", "eligibility_method_id", "interval_method_id",
    "four_n_method_id", "sub_four_n_method_id", "containment_method_id",
    "config_digest"
  )
  for (field in provenance) {
    reference <- ph3_require_one_group_value(metrics, field)
    for (name in table_names) {
      if (!identical(ph3_require_one_group_value(quantitation[[name]], field),
                     reference)) {
        ph3_aggregation_fail(
          "mixed_or_missing_provenance",
          paste0("Slice 4 tables disagree on `", field, "`")
        )
      }
    }
  }
  exact_methods <- c(
    classification_schema_version = "ph3-event-classification-1.0.0",
    output_schema_version = "ph3-1.0.0",
    eligibility_method_id =
      "identity_validated_finite_dna_b0_b5_inclusive_v1",
    interval_method_id = "configured_shared_boundaries_left_closed_v1",
    four_n_method_id = "configured_b4_b5_closed_v1",
    sub_four_n_method_id =
      "configured_b0_b4_left_closed_right_open_v1",
    containment_method_id = "exact_direct_identity_multiset_containment"
  )
  for (field in names(exact_methods)) {
    if (!identical(ph3_require_one_group_value(metrics, field),
                   exact_methods[[field]])) {
      ph3_aggregation_fail(
        "mixed_or_missing_provenance",
        paste0("Slice 4 `", field, "` is not the approved exact identifier")
      )
    }
  }
  if (!ph3_require_one_group_value(metrics, "positivity_method_id") %in% c(
    "flowjo_owner_approved_positive_population_v1",
    "ph3_raw_4n_density_cutoff_v1"
  )) {
    ph3_aggregation_fail(
      "mixed_or_missing_provenance",
      "Slice 4 `positivity_method_id` is not an approved exact identifier"
    )
  }
  invisible(list(metric_order = metric_order, phase_order = phase_order))
}

ph3_aggregate_source_group <- function(data) {
  if (!is.numeric(data$value_percent) || any(is.nan(data$value_percent)) ||
      any(is.infinite(data$value_percent)) ||
      any(is.finite(data$value_percent) & data$result_status != "ok") ||
      any(is.na(data$value_percent) &
            !data$result_status %in% c("undefined_zero_denominator", "unavailable_cutoff_failure"))) {
    ph3_aggregation_fail(
      "invalid_source_value_status", "Slice 4 values and statuses disagree"
    )
  }
  finite <- is.finite(data$value_percent)
  total_n <- as.integer(nrow(data))
  finite_n <- as.integer(sum(finite))
  undefined_n <- as.integer(total_n - finite_n)
  cutoff_failure <- any(data$result_status == "unavailable_cutoff_failure")
  list(
    value = if (finite_n) mean(data$value_percent[finite]) else NA_real_,
    total_n = total_n, finite_n = finite_n, undefined_n = undefined_n,
    status = if (cutoff_failure) "unavailable_cutoff_failure" else if (!finite_n) "undefined_no_finite_values" else if (undefined_n) {
      "ok_partial_undefined"
    } else {
      "ok"
    }
  )
}

derive_ph3_replicate_condition_tables <- function(
    quantitation, manifest, analysis_provenance
) {
  orders <- ph3_validate_slice4_aggregation_inputs(
    quantitation, manifest, analysis_provenance
  )
  metrics <- quantitation$ph3_metrics_acquisition
  phases <- quantitation$ph3_phase_prevalence
  provenance <- c(
    "classification_schema_version", "output_schema_version",
    "positivity_method_id", "eligibility_method_id", "interval_method_id",
    "four_n_method_id", "sub_four_n_method_id", "containment_method_id",
    "config_digest"
  )
  replicate_groups <- unique(manifest[c(
    "condition", "condition_index", "replicate", "replicate_index"
  )])
  replicate_groups <- replicate_groups[order(
    replicate_groups$condition_index, replicate_groups$replicate_index
  ), , drop = FALSE]

  aggregate_replicates <- function(source, identity_field, identity_order,
                                   descriptor_fields) {
    rows <- list()
    row_index <- 0L
    for (i in seq_len(nrow(replicate_groups))) {
      group <- replicate_groups[i, , drop = FALSE]
      for (identity in identity_order) {
        selected <- source[
          source$condition == group$condition &
            source$condition_index == group$condition_index &
            source$replicate == group$replicate &
            source$replicate_index == group$replicate_index &
            source[[identity_field]] == identity, , drop = FALSE
        ]
        result <- ph3_aggregate_source_group(selected)
        fixed <- c(descriptor_fields, provenance)
        fixed_values <- lapply(fixed, function(field) {
          ph3_require_one_group_value(selected, field)
        })
        names(fixed_values) <- fixed
        row_index <- row_index + 1L
        rows[[row_index]] <- c(
          list(
            analysis_id = ph3_require_one_group_value(selected, "analysis_id"),
            condition = group$condition[[1L]],
            condition_index = as.integer(group$condition_index[[1L]]),
            replicate = group$replicate[[1L]],
            replicate_index = as.integer(group$replicate_index[[1L]]),
            aggregation_level = "biological_replicate"
          ), fixed_values,
          list(
            value_percent = result$value,
            technical_acquisition_count = result$total_n,
            finite_technical_acquisition_count = result$finite_n,
            undefined_technical_acquisition_count = result$undefined_n,
            result_status = result$status,
            source_acquisition_ids =
              ph3_canonical_string_array(selected$acquisition_id),
            source_export_operation_ids =
              ph3_canonical_string_array(selected$export_operation_id),
            source_input_manifest_keys =
              ph3_canonical_string_array(selected$input_manifest_key),
            aggregation_method_id =
              "unweighted_technical_percentage_mean_within_biorep_condition_v1"
          )
        )
      }
    }
    as.data.frame(do.call(rbind, lapply(rows, as.data.frame)),
                  stringsAsFactors = FALSE)
  }

  metric_replicates <- aggregate_replicates(
    metrics, "metric_id", orders$metric_order,
    c("metric_id", "population_id", "interval_lower", "interval_upper",
      "lower_inclusive", "upper_inclusive")
  )
  phase_replicates <- aggregate_replicates(
    phases, "phase_id", orders$phase_order,
    c("metric_id", "phase_id", "phase_index", "interval_lower",
      "interval_upper", "lower_inclusive", "upper_inclusive")
  )

  aggregate_conditions <- function(source, identity_field, identity_order,
                                   descriptor_fields) {
    conditions <- unique(manifest[c("condition", "condition_index")])
    conditions <- conditions[order(conditions$condition_index), , drop = FALSE]
    rows <- list()
    row_index <- 0L
    for (i in seq_len(nrow(conditions))) {
      group <- conditions[i, , drop = FALSE]
      for (identity in identity_order) {
        selected <- source[
          source$condition == group$condition &
            source$condition_index == group$condition_index &
            source[[identity_field]] == identity, , drop = FALSE
        ]
        finite <- is.finite(selected$value_percent)
        total_n <- as.integer(nrow(selected))
        finite_n <- as.integer(sum(finite))
        undefined_n <- as.integer(total_n - finite_n)
        mean_value <- if (finite_n) mean(selected$value_percent[finite]) else NA_real_
        sd_value <- if (finite_n >= 2L) stats::sd(selected$value_percent[finite]) else NA_real_
        sem_value <- if (finite_n >= 2L) sd_value / sqrt(finite_n) else NA_real_
        cutoff_failure <- any(selected$result_status == "unavailable_cutoff_failure")
        status <- if (cutoff_failure) {
          "available_partial_unavailable_cutoff_failure"
        } else if (!finite_n) {
          "undefined_no_finite_values"
        } else if (undefined_n) {
          "ok_partial_undefined"
        } else if (finite_n == 1L && total_n == 1L) {
          "ok_single_biological_replicate"
        } else {
          "ok"
        }
        fixed <- c(descriptor_fields, provenance)
        fixed_values <- lapply(fixed, function(field) {
          ph3_require_one_group_value(selected, field)
        })
        names(fixed_values) <- fixed
        row_index <- row_index + 1L
        rows[[row_index]] <- c(
          list(
            analysis_id = ph3_require_one_group_value(selected, "analysis_id"),
            condition = group$condition[[1L]],
            condition_index = as.integer(group$condition_index[[1L]]),
            aggregation_level = "condition"
          ), fixed_values,
          list(
            mean_percent = mean_value, sd_percent = sd_value,
            sem_percent = sem_value,
            biological_replicate_count = total_n,
            finite_biological_replicate_count = finite_n,
            undefined_biological_replicate_count = undefined_n,
            result_status = status,
            source_replicate_ids =
              ph3_canonical_string_array(selected$replicate),
            aggregation_method_id =
              "unweighted_biological_replicate_mean_within_condition_v1"
          )
        )
      }
    }
    as.data.frame(do.call(rbind, lapply(rows, as.data.frame)),
                  stringsAsFactors = FALSE)
  }

  metric_conditions <- aggregate_conditions(
    metric_replicates, "metric_id", orders$metric_order,
    c("metric_id", "population_id", "interval_lower", "interval_upper",
      "lower_inclusive", "upper_inclusive")
  )
  phase_conditions <- aggregate_conditions(
    phase_replicates, "phase_id", orders$phase_order,
    c("metric_id", "phase_id", "phase_index", "interval_lower",
      "interval_upper", "lower_inclusive", "upper_inclusive")
  )

  metric_replicate_columns <- c(
    "analysis_id", "condition", "condition_index", "replicate",
    "replicate_index", "aggregation_level", "metric_id", "population_id",
    "interval_lower", "interval_upper", "lower_inclusive", "upper_inclusive",
    "value_percent", "technical_acquisition_count",
    "finite_technical_acquisition_count", "undefined_technical_acquisition_count",
    "result_status", "source_acquisition_ids", "source_export_operation_ids",
    "source_input_manifest_keys", "classification_schema_version",
    "output_schema_version", "aggregation_method_id", "positivity_method_id",
    "eligibility_method_id", "interval_method_id", "four_n_method_id",
    "sub_four_n_method_id", "containment_method_id", "config_digest"
  )
  phase_replicate_columns <- append(metric_replicate_columns, "phase_id", 7L)
  phase_replicate_columns <- append(phase_replicate_columns, "phase_index", 8L)
  phase_replicate_columns <- phase_replicate_columns[
    phase_replicate_columns != "population_id"
  ]
  metric_condition_columns <- c(
    "analysis_id", "condition", "condition_index", "aggregation_level",
    "metric_id", "population_id", "interval_lower", "interval_upper",
    "lower_inclusive", "upper_inclusive", "mean_percent", "sd_percent",
    "sem_percent", "biological_replicate_count",
    "finite_biological_replicate_count", "undefined_biological_replicate_count",
    "result_status", "source_replicate_ids", "classification_schema_version",
    "output_schema_version", "aggregation_method_id", "positivity_method_id",
    "eligibility_method_id", "interval_method_id", "four_n_method_id",
    "sub_four_n_method_id", "containment_method_id", "config_digest"
  )
  phase_condition_columns <- append(metric_condition_columns, "phase_id", 5L)
  phase_condition_columns <- append(phase_condition_columns, "phase_index", 6L)
  phase_condition_columns <- phase_condition_columns[
    phase_condition_columns != "population_id"
  ]
  tables <- list(
    ph3_metrics_biological_replicate =
      metric_replicates[metric_replicate_columns],
    ph3_phase_prevalence_biological_replicate =
      phase_replicates[phase_replicate_columns],
    ph3_metrics_condition_summary =
      metric_conditions[metric_condition_columns],
    ph3_phase_prevalence_condition_summary =
      phase_conditions[phase_condition_columns]
  )
  lapply(tables, function(data) {
    rownames(data) <- NULL
    data
  })
}

quantify_ph3_production_acquisitions <- function(analysis) {
  if (!identical(analysis$config$plot_type, "ph3") ||
      !identical(analysis$config$ph3_input_profile,
                 "production_direct_identity_v1")) {
    stop("Production PH3 acquisition metrics require the production PH3 profile.",
         call. = FALSE)
  }
  tables <- lapply(seq_len(nrow(analysis$sample_manifest)), function(i) {
    classification <-
      analysis$normalized_data[[i]]$ph3_event_classification
    derive_ph3_acquisition_tables(
      classification, analysis$sample_manifest[i, , drop = FALSE],
      analysis$config, analysis$provenance$ph3_export_manifests,
      analysis$provenance$ph3_containment
    )
  })
  names_to_attach <- names(tables[[1L]])
  for (name in names_to_attach) {
    analysis$quantitation[[name]] <- do.call(
      rbind, lapply(tables, `[[`, name)
    )
    rownames(analysis$quantitation[[name]]) <- NULL
  }
  aggregate_tables <- derive_ph3_replicate_condition_tables(
    analysis$quantitation, analysis$sample_manifest, analysis$provenance
  )
  for (name in names(aggregate_tables)) {
    analysis$quantitation[[name]] <- aggregate_tables[[name]]
  }
  analysis
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
