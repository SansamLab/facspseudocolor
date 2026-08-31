# Validated pH3 output-contract model ---------------------------------------

ph3_output_contract_fail <- function(reason, detail) {
  stop(
    "PH3 output-contract validation failed [", reason, "]: ", detail, ".",
    call. = FALSE
  )
}

ph3_output_contract_exact_names <- function(value, expected, context) {
  if (!is.list(value) || is.data.frame(value) ||
      !identical(names(value), expected)) {
    ph3_output_contract_fail(
      "invalid_schema",
      paste0("`", context, "` must contain exactly these ordered keys: ",
             paste(expected, collapse = ", "))
    )
  }
}

ph3_output_contract_key <- function(value) {
  config_scalar_string(value) &&
    grepl("^[A-Za-z0-9][A-Za-z0-9._-]*$", value)
}

ph3_output_contract_string_vector <- function(value) {
  is.character(value) && length(value) > 0L && !anyNA(value) &&
    all(nzchar(value)) && !anyDuplicated(value)
}

validate_ph3_output_contract_config <- function(contract, replicates) {
  ph3_output_contract_exact_names(
    contract,
    c("schema_version", "experiment_id", "conditions", "reference",
      "comparisons", "geometry"),
    "ph3_output_contract"
  )
  if (!identical(contract$schema_version,
                 "ph3-output-contract-config-1.0.0")) {
    ph3_output_contract_fail(
      "unsupported_schema_version",
      "`ph3_output_contract.schema_version` must be `ph3-output-contract-config-1.0.0`"
    )
  }
  if (!ph3_output_contract_key(contract$experiment_id)) {
    ph3_output_contract_fail(
      "invalid_experiment_id",
      "`ph3_output_contract.experiment_id` must be one stable filesystem-safe key"
    )
  }
  if (!is.list(replicates) || !length(replicates)) {
    ph3_output_contract_fail(
      "missing_replicate_sets", "production pH3 requires explicit replicates"
    )
  }
  replicate_ids <- vapply(replicates, function(replicate) {
    value <- replicate$id
    if (!ph3_output_contract_key(value)) NA_character_ else value
  }, character(1))
  if (anyNA(replicate_ids) || anyDuplicated(replicate_ids)) {
    ph3_output_contract_fail(
      "invalid_replicate_set_identity",
      "every pH3 replicate must declare one unique stable `id`"
    )
  }

  conditions <- contract$conditions
  if (!is.list(conditions) || !length(conditions)) {
    ph3_output_contract_fail(
      "missing_conditions", "`ph3_output_contract.conditions` cannot be empty"
    )
  }
  condition_rows <- lapply(seq_along(conditions), function(i) {
    condition <- conditions[[i]]
    ph3_output_contract_exact_names(
      condition, c("id", "label", "role"),
      paste0("ph3_output_contract.conditions[[", i, "]]"))
    if (!ph3_output_contract_key(condition$id) ||
        !config_scalar_string(condition$label) ||
        !config_scalar_string(condition$role) ||
        !condition$role %in% c("control", "treatment")) {
      ph3_output_contract_fail(
        "invalid_condition",
        "each condition requires a stable `id`, nonempty `label`, and role `control` or `treatment`"
      )
    }
    data.frame(
      condition_id = condition$id, condition_label = condition$label,
      role = condition$role, condition_index = as.integer(i),
      stringsAsFactors = FALSE
    )
  })
  condition_table <- do.call(rbind, condition_rows)
  if (anyDuplicated(condition_table$condition_id) ||
      anyDuplicated(condition_table$condition_label)) {
    ph3_output_contract_fail(
      "duplicate_condition_identity", "condition IDs and labels must be unique"
    )
  }
  configured_condition_labels <- unique(vapply(
    replicates[[1L]]$samples, `[[`, character(1), "label"
  ))
  if (!identical(condition_table$condition_label,
                 configured_condition_labels)) {
    ph3_output_contract_fail(
      "condition_manifest_mismatch",
      "condition records must match the first replicate's explicit condition labels and order"
    )
  }

  reference <- contract$reference
  if (!is.list(reference) || is.data.frame(reference) ||
      !config_scalar_string(reference$status) ||
      !reference$status %in% c("not_configured", "configured")) {
    ph3_output_contract_fail(
      "invalid_reference_state",
      "reference status must be exactly `not_configured` or `configured`"
    )
  }
  if (identical(reference$status, "not_configured")) {
    ph3_output_contract_exact_names(
      reference, "status", "ph3_output_contract.reference"
    )
  } else {
    ph3_output_contract_exact_names(
      reference, c("status", "condition_id"),
      "ph3_output_contract.reference"
    )
    if (!ph3_output_contract_key(reference$condition_id) ||
        !reference$condition_id %in% condition_table$condition_id) {
      ph3_output_contract_fail(
        "unknown_reference_condition",
        "configured reference must name one declared condition ID"
      )
    }
  }

  comparisons <- contract$comparisons
  if (!is.list(comparisons) || is.data.frame(comparisons) ||
      !is.null(names(comparisons))) {
    ph3_output_contract_fail(
      "invalid_comparison_schema",
      "`ph3_output_contract.comparisons` must be an explicitly supplied unnamed list"
    )
  }
  comparison_rows <- list()
  for (i in seq_along(comparisons)) {
    comparison <- comparisons[[i]]
    ph3_output_contract_exact_names(
      comparison,
      c("id", "control_condition_id", "treatment_condition_id",
        "outcome_ids"),
      paste0("ph3_output_contract.comparisons[[", i, "]]"))
    if (!ph3_output_contract_key(comparison$id) ||
        !ph3_output_contract_key(comparison$control_condition_id) ||
        !ph3_output_contract_key(comparison$treatment_condition_id) ||
        !ph3_output_contract_string_vector(comparison$outcome_ids) ||
        !all(comparison$outcome_ids %in% c("A", "B", "C", "D"))) {
      ph3_output_contract_fail(
        "invalid_comparison",
        "each comparison requires stable condition keys and unique outcome IDs drawn from A-D"
      )
    }
    control <- condition_table[
      condition_table$condition_id == comparison$control_condition_id,
      , drop = FALSE
    ]
    treatment <- condition_table[
      condition_table$condition_id == comparison$treatment_condition_id,
      , drop = FALSE
    ]
    if (nrow(control) != 1L || nrow(treatment) != 1L) {
      ph3_output_contract_fail(
        "unknown_comparison_condition",
        "comparison conditions must both be declared"
      )
    }
    if (identical(comparison$control_condition_id,
                  comparison$treatment_condition_id)) {
      ph3_output_contract_fail(
        "self_comparison", "control and treatment conditions must differ"
      )
    }
    if (!identical(control$role, "control") ||
        !identical(treatment$role, "treatment")) {
      ph3_output_contract_fail(
        "comparison_direction_mismatch",
        "comparison direction must agree with control/treatment metadata"
      )
    }
    comparison_rows[[i]] <- data.frame(
      comparison_id = comparison$id,
      control_condition_id = comparison$control_condition_id,
      treatment_condition_id = comparison$treatment_condition_id,
      outcome_id = comparison$outcome_ids,
      comparison_index = as.integer(i),
      outcome_index = as.integer(match(comparison$outcome_ids,
                                       c("A", "B", "C", "D"))),
      stringsAsFactors = FALSE
    )
  }
  comparison_table <- if (length(comparison_rows)) {
    do.call(rbind, comparison_rows)
  } else {
    data.frame(
      comparison_id = character(), control_condition_id = character(),
      treatment_condition_id = character(), outcome_id = character(),
      comparison_index = integer(), outcome_index = integer(),
      stringsAsFactors = FALSE
    )
  }
  if (anyDuplicated(vapply(comparisons, function(x) x$id, character(1))) ||
      anyDuplicated(comparison_table[c(
        "control_condition_id", "treatment_condition_id", "outcome_id"
      )])) {
    ph3_output_contract_fail(
      "duplicate_comparison", "comparison IDs and directed outcome pairs must be unique"
    )
  }

  geometry <- contract$geometry
  ph3_output_contract_exact_names(
    geometry, c("verified", "computed"), "ph3_output_contract.geometry"
  )
  verified <- geometry$verified
  if (!is.list(verified) || !config_scalar_string(verified$status) ||
      !verified$status %in% c("not_configured", "geometry_set_index_declared")) {
    ph3_output_contract_fail(
      "invalid_verified_geometry_state",
      "verified geometry status must be `not_configured` or `geometry_set_index_declared`"
    )
  }
  if (identical(verified$status, "not_configured")) {
    ph3_output_contract_exact_names(
      verified, "status", "ph3_output_contract.geometry.verified"
    )
  } else {
    ph3_output_contract_exact_names(
      verified, c("status", "geometry_set_index_id"),
      "ph3_output_contract.geometry.verified"
    )
    if (!ph3_output_contract_key(verified$geometry_set_index_id)) {
      ph3_output_contract_fail(
        "missing_geometry_set_index",
        "declared verified geometry requires one stable geometry-set index ID"
      )
    }
  }
  ph3_output_contract_exact_names(
    geometry$computed, "status", "ph3_output_contract.geometry.computed"
  )
  if (!identical(geometry$computed$status, "not_computed")) {
    ph3_output_contract_fail(
      "invalid_computed_geometry_state",
      "Slice 1 requires computed geometry status `not_computed`"
    )
  }

  attr(contract, "condition_table") <- condition_table
  attr(contract, "comparison_table") <- comparison_table
  attr(contract, "replicate_set_ids") <- replicate_ids
  contract
}

ph3_output_model_schema <- function() {
  list(
    schema_version = "ph3-output-contract-model-1.0.0",
    experiment_columns = c("experiment_id", "analysis_id"),
    experiment_types = c(
      experiment_id = "character", analysis_id = "character"
    ),
    outcome_columns = c(
      "outcome_id", "outcome_name", "label", "population_id",
      "denominator_population_id", "value_kind", "unit",
      "source_metric_id", "signal_basis_applicability",
      "reference_applicability"
    ),
    outcome_types = c(
      outcome_id = "character", outcome_name = "character",
      label = "character", population_id = "character",
      denominator_population_id = "character", value_kind = "character",
      unit = "character", source_metric_id = "character",
      signal_basis_applicability = "character",
      reference_applicability = "character"
    ),
    outcome_ids = c("A", "B", "C", "D"),
    condition_columns = c(
      "condition_id", "condition_label", "role", "condition_index"
    ),
    condition_types = c(
      condition_id = "character", condition_label = "character",
      role = "character", condition_index = "integer"
    ),
    replicate_set_columns = c(
      "experiment_id", "replicate_set_id", "replicate_label",
      "manifest_replicate_index"
    ),
    replicate_set_types = c(
      experiment_id = "character", replicate_set_id = "character",
      replicate_label = "character", manifest_replicate_index = "integer"
    ),
    sample_columns = c(
      "experiment_id", "replicate_set_id", "sample_id", "condition_id",
      "acquisition_count", "source_acquisition_ids"
    ),
    sample_types = c(
      experiment_id = "character", replicate_set_id = "character",
      sample_id = "character", condition_id = "character",
      acquisition_count = "integer", source_acquisition_ids = "character"
    ),
    correction_columns = c(
      "experiment_id", "replicate_set_id", "status", "signal_basis",
      "reason_code", "reason_detail"
    ),
    correction_types = c(
      experiment_id = "character", replicate_set_id = "character",
      status = "character", signal_basis = "character",
      reason_code = "character", reason_detail = "character"
    ),
    correction_status = c("not_computed", "selected"),
    signal_basis = c("individual_corrected", "pooled_corrected", "raw"),
    correction_reason_codes = "background_correction_not_computed_slice_1",
    reference_columns = c(
      "experiment_id", "replicate_set_id", "status", "condition_id",
      "sample_id", "reason_code", "reason_detail"
    ),
    reference_types = c(
      experiment_id = "character", replicate_set_id = "character",
      status = "character", condition_id = "character",
      sample_id = "character", reason_code = "character",
      reason_detail = "character"
    ),
    reference_status = c("configured", "not_configured"),
    reference_reason_codes = "reference_not_configured",
    comparison_columns = c(
      "comparison_id", "control_condition_id", "treatment_condition_id",
      "outcome_id", "comparison_index", "outcome_index"
    ),
    comparison_types = c(
      comparison_id = "character", control_condition_id = "character",
      treatment_condition_id = "character", outcome_id = "character",
      comparison_index = "integer", outcome_index = "integer"
    ),
    geometry_columns = c(
      "geometry_class", "geometry_id", "required_provenance_label", "status",
      "source_id", "reason_code", "reason_detail"
    ),
    geometry_types = c(
      geometry_class = "character", geometry_id = "character",
      required_provenance_label = "character", status = "character",
      source_id = "character", reason_code = "character",
      reason_detail = "character"
    ),
    geometry_status = c(
      "not_configured", "geometry_set_index_declared", "not_computed"
    ),
    geometry_reason_codes = c(
      "verified_geometry_not_configured",
      "computed_geometry_not_generated_slice_1"
    ),
    artifact_columns = c(
      "artifact_id", "scope", "format", "media_type", "path_template",
      "required"
    ),
    artifact_types = c(
      artifact_id = "character", scope = "character", format = "character",
      media_type = "character", path_template = "character",
      required = "logical"
    )
  )
}

ph3_validate_model_table <- function(data, columns, types, context) {
  actual_types <- vapply(data, typeof, character(1))
  if (!is.data.frame(data) || !identical(names(data), columns) ||
      !identical(actual_types, types)) {
    ph3_output_contract_fail(
      "internal_schema_mismatch",
      paste0("`", context, "` changed ordered columns or scalar types")
    )
  }
  invisible(data)
}

ph3_owner_confirmed_outcomes <- function() {
  data.frame(
    outcome_id = c("A", "B", "C", "D"),
    outcome_name = c(
      "ph3_4n_positive_prevalence",
      "ph3_below_4n_positive_prevalence",
      "ph3_4n_positive_signal",
      "ph3_below_4n_positive_signal"
    ),
    label = c(
      "4N pH3-positive prevalence",
      "Below-4N pH3-positive prevalence",
      "pH3 signal in 4N pH3-positive cells",
      "pH3 signal in below-4N pH3-positive cells"
    ),
    population_id = c(
      "ph3_positive_4n", "ph3_positive_below_4n",
      "ph3_positive_4n", "ph3_positive_below_4n"
    ),
    denominator_population_id = c(
      "eligible_2to4n", "eligible_2to4n", NA_character_, NA_character_
    ),
    value_kind = c(
      "prevalence_percent", "prevalence_percent",
      "population_signal_median", "population_signal_median"
    ),
    unit = c("percent", "percent", "source_signal_unit", "source_signal_unit"),
    source_metric_id = c(
      "ph3_4n_positive_prevalence_within_2to4n_percent",
      "ph3_sub_4n_positive_prevalence_within_2to4n_percent",
      NA_character_, NA_character_
    ),
    signal_basis_applicability = c(
      "not_applicable", "not_applicable", "required", "required"
    ),
    reference_applicability = c(
      "not_applicable", "not_applicable", "optional", "optional"
    ),
    stringsAsFactors = FALSE
  )
}

ph3_required_output_registry <- function() {
  roots <- data.frame(
    artifact_id = c(
      "report", "sample_results", "condition_summary", "regression_qc",
      "provenance_manifest"
    ),
    scope = "package", format = c("html", "csv", "csv", "csv", "json"),
    media_type = c(
      "text/html", "text/csv", "text/csv", "text/csv",
      "application/json"
    ),
    path_template = c(
      "ph3_report.html", "sample_results.csv", "condition_summary.csv",
      "regression_qc.csv", "provenance_manifest.json"
    ),
    required = TRUE, stringsAsFactors = FALSE
  )
  condition_names <- c(
    "outcome_A_4n_ph3_positive_prevalence",
    "outcome_B_below4n_ph3_positive_prevalence",
    "outcome_C_4n_ph3_positive_signal",
    "outcome_D_below4n_ph3_positive_signal"
  )
  condition <- do.call(rbind, lapply(condition_names, function(name) {
    data.frame(
      artifact_id = paste0("condition_plot_", name, "_", c("pdf", "png")),
      scope = "condition", format = c("pdf", "png"),
      media_type = c("application/pdf", "image/png"),
      path_template = paste0("plots/condition/", name, c(".pdf", ".png")),
      required = TRUE, stringsAsFactors = FALSE
    )
  }))
  sample_names <- c(
    "dna_vs_ph3_pseudocolor", "gating_overview", "dna_content_distribution"
  )
  sample <- do.call(rbind, lapply(sample_names, function(name) {
    data.frame(
      artifact_id = paste0("sample_plot_", name, "_", c("pdf", "png")),
      scope = "sample", format = c("pdf", "png"),
      media_type = c("application/pdf", "image/png"),
      path_template = paste0(
        "plots/samples/{sample_id}/", name, c(".pdf", ".png")
      ),
      required = TRUE, stringsAsFactors = FALSE
    )
  }))
  registry <- rbind(roots, condition, sample)
  rownames(registry) <- NULL
  registry
}

ph3_validate_output_source_tables <- function(analysis) {
  required <- c(
    "ph3_metrics_acquisition", "ph3_phase_prevalence",
    "ph3_event_eligibility_qc", "ph3_4n_boundary_sensitivity_qc",
    "ph3_metrics_biological_replicate",
    "ph3_phase_prevalence_biological_replicate",
    "ph3_metrics_condition_summary",
    "ph3_phase_prevalence_condition_summary"
  )
  if (!identical(names(analysis$quantitation), required) ||
      !all(vapply(analysis$quantitation, is.data.frame, logical(1)))) {
    ph3_output_contract_fail(
      "missing_or_extra_source_table",
      "the exact authoritative Slice 4 and Slice 5 table set is required"
    )
  }
  classifications <- lapply(
    analysis$normalized_data, `[[`, "ph3_event_classification"
  )
  if (length(classifications) != nrow(analysis$sample_manifest) ||
      any(!vapply(classifications, is.data.frame, logical(1)))) {
    ph3_output_contract_fail(
      "missing_classification",
      "one authoritative event-classification table is required per manifest acquisition"
    )
  }
  expected_slice4 <- lapply(seq_len(nrow(analysis$sample_manifest)), function(i) {
    derive_ph3_acquisition_tables(
      classifications[[i]], analysis$sample_manifest[i, , drop = FALSE],
      analysis$config, analysis$provenance$ph3_export_manifests,
      analysis$provenance$ph3_containment
    )
  })
  for (name in required[1:4]) {
    expected <- do.call(rbind, lapply(expected_slice4, `[[`, name))
    rownames(expected) <- NULL
    if (!identical(analysis$quantitation[[name]], expected)) {
      ph3_output_contract_fail(
        "slice4_source_mismatch",
        paste0("`", name, "` is not byte-for-byte derived from active membership")
      )
    }
  }
  expected_slice5 <- derive_ph3_replicate_condition_tables(
    analysis$quantitation, analysis$sample_manifest, analysis$provenance
  )
  for (name in required[5:8]) {
    if (!identical(analysis$quantitation[[name]], expected_slice5[[name]])) {
      ph3_output_contract_fail(
        "slice5_source_mismatch",
        paste0("`", name, "` is not byte-for-byte derived from Slice 4")
      )
    }
  }
  classifications
}

build_ph3_output_model <- function(analysis) {
  if (!inherits(analysis, "facs_analysis") ||
      !identical(analysis$config$plot_type, "ph3") ||
      !identical(analysis$config$ph3_input_profile,
                 "production_direct_identity_v1")) {
    ph3_output_contract_fail(
      "invalid_analysis_profile",
      "the model requires one production direct-identity pH3 analysis"
    )
  }
  contract <- validate_ph3_output_contract_config(
    analysis$config$ph3_output_contract, analysis$config$replicates
  )
  classifications <- ph3_validate_output_source_tables(analysis)
  conditions <- attr(contract, "condition_table")
  comparisons <- attr(contract, "comparison_table")
  replicate_ids <- attr(contract, "replicate_set_ids")
  manifest <- analysis$sample_manifest
  expected_manifest <- make_sample_manifest(replicates = analysis$config$replicates)
  expected_manifest$replicate_set_id <- replicate_ids[
    expected_manifest$replicate_index
  ]
  manifest_identity <- c(
    "replicate_set_id", "replicate", "replicate_index",
    "technical_replicate", "condition", "condition_index", "prefix"
  )
  if (!all(manifest_identity %in% names(manifest)) ||
      !identical(manifest[manifest_identity],
                 expected_manifest[manifest_identity])) {
    ph3_output_contract_fail(
      "replicate_set_manifest_mismatch",
      "the active manifest must retain the explicit configured replicate-set mapping"
    )
  }
  manifest_replicate_mapping <- unique(manifest[c(
    "replicate_set_id", "replicate_index"
  )])
  rownames(manifest_replicate_mapping) <- NULL
  if (length(replicate_ids) != length(unique(manifest$replicate_index)) ||
      any(!manifest$replicate_index %in% seq_along(replicate_ids)) ||
      !identical(
        manifest_replicate_mapping,
        data.frame(
          replicate_set_id = replicate_ids,
          replicate_index = as.integer(seq_along(replicate_ids)),
          stringsAsFactors = FALSE
        )
      )) {
    ph3_output_contract_fail(
      "replicate_set_manifest_mismatch",
      "stable replicate-set IDs do not map exactly to the active manifest"
    )
  }
  condition_ids <- stats::setNames(
    conditions$condition_id, conditions$condition_label
  )
  if (any(!manifest$condition %in% names(condition_ids))) {
    ph3_output_contract_fail(
      "condition_manifest_mismatch",
      "active manifest conditions do not map to declared condition IDs"
    )
  }
  acquisition_mapping <- unique(analysis$quantitation$ph3_metrics_acquisition[c(
    "acquisition_id", "sample_id", "prefix", "condition", "replicate_index"
  )])
  if (nrow(acquisition_mapping) != nrow(manifest) ||
      anyDuplicated(acquisition_mapping$acquisition_id) ||
      anyDuplicated(acquisition_mapping$prefix)) {
    ph3_output_contract_fail(
      "unstable_sample_identity",
      "active provenance must map every acquisition and prefix exactly once"
    )
  }
  acquisition_mapping$experiment_id <- contract$experiment_id
  acquisition_mapping$replicate_set_id <- manifest$replicate_set_id[
    match(acquisition_mapping$prefix, manifest$prefix)
  ]
  acquisition_mapping$condition_id <- unname(
    condition_ids[acquisition_mapping$condition]
  )
  sample_groups <- unique(acquisition_mapping[c(
    "sample_id", "experiment_id", "replicate_set_id", "condition_id"
  )])
  if (anyDuplicated(sample_groups$sample_id)) {
    ph3_output_contract_fail(
      "unstable_sample_identity",
      "one explicit sample ID cannot span replicate sets or conditions"
    )
  }
  samples <- do.call(rbind, lapply(seq_len(nrow(sample_groups)), function(i) {
    group <- sample_groups[i, , drop = FALSE]
    rows <- acquisition_mapping[
      acquisition_mapping$sample_id == group$sample_id, , drop = FALSE
    ]
    data.frame(
      experiment_id = group$experiment_id,
      replicate_set_id = group$replicate_set_id,
      sample_id = group$sample_id, condition_id = group$condition_id,
      acquisition_count = as.integer(nrow(rows)),
      source_acquisition_ids = ph3_canonical_string_array(rows$acquisition_id),
      stringsAsFactors = FALSE
    )
  }))
  rownames(samples) <- NULL
  source_sample_mapping <- samples[c(
    "experiment_id", "replicate_set_id", "sample_id", "condition_id"
  )]
  replicate_sets <- data.frame(
    experiment_id = contract$experiment_id,
    replicate_set_id = replicate_ids,
    replicate_label = vapply(
      analysis$config$replicates, `[[`, character(1), "label"
    ),
    manifest_replicate_index = as.integer(seq_along(replicate_ids)),
    stringsAsFactors = FALSE
  )
  correction <- data.frame(
    experiment_id = contract$experiment_id,
    replicate_set_id = replicate_ids,
    status = "not_computed", signal_basis = NA_character_,
    reason_code = "background_correction_not_computed_slice_1",
    reason_detail = NA_character_, stringsAsFactors = FALSE
  )
  reference_condition <- if (identical(contract$reference$status,
                                        "configured")) {
    contract$reference$condition_id
  } else NA_character_
  reference <- do.call(rbind, lapply(replicate_ids, function(replicate_id) {
    if (identical(contract$reference$status, "not_configured")) {
      sample_id <- NA_character_
      reason_code <- "reference_not_configured"
    } else {
      matched <- samples[
        samples$replicate_set_id == replicate_id &
          samples$condition_id == reference_condition, , drop = FALSE
      ]
      if (nrow(matched) != 1L) {
        ph3_output_contract_fail(
          "reference_resolution_failure",
          paste0("reference condition must resolve to exactly one sample in `",
                 replicate_id, "`")
        )
      }
      sample_id <- matched$sample_id[[1L]]
      reason_code <- NA_character_
    }
    data.frame(
      experiment_id = contract$experiment_id,
      replicate_set_id = replicate_id, status = contract$reference$status,
      condition_id = reference_condition, sample_id = sample_id,
      reason_code = reason_code, reason_detail = NA_character_,
      stringsAsFactors = FALSE
    )
  }))
  rownames(reference) <- NULL
  verified <- contract$geometry$verified
  verified_source <- if (identical(verified$status,
                                    "geometry_set_index_declared")) {
    verified$geometry_set_index_id
  } else NA_character_
  geometry <- data.frame(
    geometry_class = c("verified", "computed"),
    geometry_id = c("source_gate_geometry", "analysis_region_geometry"),
    required_provenance_label = c("VERIFIED", "COMPUTED"),
    status = c(verified$status, "not_computed"),
    source_id = c(verified_source, NA_character_),
    reason_code = c(
      if (identical(verified$status, "not_configured")) {
        "verified_geometry_not_configured"
      } else NA_character_,
      "computed_geometry_not_generated_slice_1"
    ),
    reason_detail = NA_character_, stringsAsFactors = FALSE
  )
  schema <- ph3_output_model_schema()
  outcomes <- ph3_owner_confirmed_outcomes()
  registry <- ph3_required_output_registry()
  ph3_validate_model_table(
    outcomes, schema$outcome_columns, schema$outcome_types, "outcomes"
  )
  ph3_validate_model_table(
    correction, schema$correction_columns, schema$correction_types,
    "correction"
  )
  ph3_validate_model_table(
    reference, schema$reference_columns, schema$reference_types, "reference"
  )
  ph3_validate_model_table(
    comparisons, schema$comparison_columns, schema$comparison_types,
    "comparisons"
  )
  ph3_validate_model_table(
    geometry, schema$geometry_columns, schema$geometry_types, "geometry"
  )
  ph3_validate_model_table(
    registry, schema$artifact_columns, schema$artifact_types,
    "output_registry"
  )
  analysis_ids <- unique(
    analysis$quantitation$ph3_metrics_acquisition$analysis_id
  )
  if (!identical(outcomes$outcome_id, schema$outcome_ids) ||
      length(analysis_ids) != 1L || is.na(analysis_ids) ||
      !nzchar(analysis_ids)) {
    ph3_output_contract_fail(
      "mixed_experiment_scope",
      "one model requires the exact A-D registry and one active analysis identity"
    )
  }
  experiment <- data.frame(
    experiment_id = contract$experiment_id,
    analysis_id = analysis_ids, stringsAsFactors = FALSE
  )
  ph3_validate_model_table(
    experiment, schema$experiment_columns, schema$experiment_types,
    "experiment"
  )
  ph3_validate_model_table(
    conditions, schema$condition_columns, schema$condition_types,
    "conditions"
  )
  ph3_validate_model_table(
    replicate_sets, schema$replicate_set_columns,
    schema$replicate_set_types, "replicate_sets"
  )
  ph3_validate_model_table(
    samples, schema$sample_columns, schema$sample_types, "samples"
  )
  structure(list(
    schema = schema,
    experiment = experiment,
    outcomes = outcomes, conditions = conditions,
    replicate_sets = replicate_sets, samples = samples,
    correction = correction, reference = reference,
    comparisons = comparisons, geometry = geometry,
    output_registry = registry,
    source = list(
      sample_manifest = manifest,
      sample_mapping = source_sample_mapping,
      event_classifications = classifications,
      quantitation = analysis$quantitation
    )
  ), class = "ph3_output_model")
}
