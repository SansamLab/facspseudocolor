# SYNTHETIC aggregation arithmetic only. No value in this file is an
# experimental observation or biological conclusion.

synthetic_slice5_inputs <- function() {
  manifest <- data.frame(
    prefix = paste0("SYNTHETIC_prefix_", 1:6),
    condition = c(rep("SYNTHETIC A", 4), rep("SYNTHETIC B", 2)),
    condition_index = c(rep(1L, 4), rep(2L, 2)),
    replicate = c(rep("SYNTHETIC R1", 3), "SYNTHETIC R2",
                  "SYNTHETIC R3", "SYNTHETIC R4"),
    replicate_index = c(1L, 1L, 1L, 2L, 3L, 4L),
    technical_replicate = c("1", "2", "3", "1", "1", "1"),
    stringsAsFactors = FALSE
  )
  acquisition_values <- c(10, NA, 30, 40, 0, NA)
  common <- function(i, n) data.frame(
    analysis_id = "SYNTHETIC-analysis",
    acquisition_id = paste0("SYNTHETIC-A", i),
    sample_id = paste0("SYNTHETIC-sample-", i),
    prefix = manifest$prefix[[i]], condition = manifest$condition[[i]],
    condition_index = manifest$condition_index[[i]],
    replicate = manifest$replicate[[i]],
    replicate_index = manifest$replicate_index[[i]],
    technical_replicate = manifest$technical_replicate[[i]],
    aggregation_level = "acquisition",
    classification_schema_version = "ph3-event-classification-1.0.0",
    output_schema_version = "ph3-1.0.0",
    positivity_method_id = "flowjo_owner_approved_positive_population_v1",
    eligibility_method_id =
      "identity_validated_finite_dna_b0_b5_inclusive_v1",
    interval_method_id = "configured_shared_boundaries_left_closed_v1",
    four_n_method_id = "configured_b4_b5_closed_v1",
    sub_four_n_method_id =
      "configured_b0_b4_left_closed_right_open_v1",
    containment_method_id = "exact_direct_identity_multiset_containment",
    config_digest = "SYNTHETIC-config",
    export_operation_id = paste0("SYNTHETIC-operation-", i),
    input_manifest_key = paste0("SYNTHETIC-manifest-", i),
    stringsAsFactors = FALSE
  )[rep(1L, n), , drop = FALSE]
  metric_order <- c(
    "ph3_2to4n_positivity_percent", "ph3_within_4n_positivity_percent",
    "ph3_4n_positive_prevalence_within_2to4n_percent",
    "ph3_within_sub_4n_positivity_percent",
    "ph3_sub_4n_positive_prevalence_within_2to4n_percent"
  )
  phase_order <- c("G1", "Early S", "Mid S", "Late S", "G2/M")
  metric_tables <- lapply(seq_len(nrow(manifest)), function(i) {
    cbind(common(i, 5L), data.frame(
      metric_id = metric_order, value_percent = rep(acquisition_values[[i]], 5L),
      numerator_count = 1:5, denominator_count = rep(10L, 5L),
      population_id = paste0("SYNTHETIC-population-", 1:5),
      interval_lower = 0:4, interval_upper = 1:5,
      lower_inclusive = TRUE, upper_inclusive = c(TRUE, TRUE, TRUE, FALSE, FALSE),
      result_status = if (is.na(acquisition_values[[i]])) {
        "undefined_zero_denominator"
      } else "ok", stringsAsFactors = FALSE
    ))
  })
  phase_tables <- lapply(seq_len(nrow(manifest)), function(i) {
    cbind(common(i, 5L), data.frame(
      phase_id = phase_order, phase_index = 1:5,
      interval_lower = 0:4, interval_upper = 1:5,
      lower_inclusive = TRUE,
      upper_inclusive = c(FALSE, FALSE, FALSE, FALSE, TRUE),
      metric_id = "ph3_phase_positive_prevalence_within_2to4n_percent",
      ph3_positive_count = 1:5, eligible_2to4n_count = rep(10L, 5L),
      value_percent = rep(acquisition_values[[i]], 5L),
      result_status = if (is.na(acquisition_values[[i]])) {
        "undefined_zero_denominator"
      } else "ok", stringsAsFactors = FALSE
    ))
  })
  qc_tables <- lapply(seq_len(nrow(manifest)), function(i) {
    cbind(common(i, 16L), data.frame(
      qc_dimension = c(
        "rows", "identity", "eligibility", rep("eligibility_exclusion", 4L),
        rep("partition", 2L), rep("configured_phase", 6L),
        "eligible_configured_phase_unassigned"
      ),
      reason = c(
        "imported_classification_rows", "identity_valid", "eligible_2to4n",
        "identity_invalid", "dna_nonfinite", "below_b0", "above_b5",
        "sub_4n", "4n", "configured_phase_assigned", "G1", "Early S",
        "Mid S", "Late S", "G2/M", "none"
      ), stringsAsFactors = FALSE
    ))
  })
  sensitivity_tables <- lapply(seq_len(nrow(manifest)), function(i) {
    cbind(common(i, 4L), data.frame(
      perturbation_id = c("primary", "lower_outward", "lower_inward",
                          "upper_inward"), stringsAsFactors = FALSE
    ))
  })
  acquisitions <- lapply(seq_len(nrow(manifest)), function(i) list(
    acquisition_id = paste0("SYNTHETIC-A", i),
    sample_id = paste0("SYNTHETIC-sample-", i), prefix = manifest$prefix[[i]]
  ))
  containment <- do.call(rbind, lapply(seq_len(nrow(manifest)), function(i) {
    data.frame(
      acquisition_id = paste0("SYNTHETIC-A", i), prefix = manifest$prefix[[i]],
      child_population_key = c("g1", "ph3_positive"),
      containment_status = "validated",
      export_operation_id = paste0("SYNTHETIC-operation-", i),
      manifest_digest = paste0("SYNTHETIC-manifest-", i),
      stringsAsFactors = FALSE
    )
  }))
  provenance <- list(
    ph3_export_manifests = lapply(seq_len(nrow(manifest)), function(i) list(
      sha256 = paste0("SYNTHETIC-manifest-", i),
      export_operation_id = paste0("SYNTHETIC-operation-", i),
      manifest = list(acquisitions = list(acquisitions[[i]]))
    )),
    ph3_containment = containment
  )
  list(manifest = manifest, provenance = provenance, quantitation = list(
    ph3_metrics_acquisition = do.call(rbind, metric_tables),
    ph3_phase_prevalence = do.call(rbind, phase_tables),
    ph3_event_eligibility_qc = do.call(rbind, qc_tables),
    ph3_4n_boundary_sensitivity_qc = do.call(rbind, sensitivity_tables)
  ))
}

test_that("SYNTHETIC technical acquisitions cannot create pseudoreplicates", {
  inputs <- synthetic_slice5_inputs()
  original <- inputs$quantitation
  result <- derive_ph3_replicate_condition_tables(
    inputs$quantitation, inputs$manifest, inputs$provenance
  )
  expect_identical(inputs$quantitation, original)
  expect_identical(names(result), c(
    "ph3_metrics_biological_replicate",
    "ph3_phase_prevalence_biological_replicate",
    "ph3_metrics_condition_summary",
    "ph3_phase_prevalence_condition_summary"
  ))
  replicate_metrics <- result$ph3_metrics_biological_replicate
  first <- replicate_metrics[
    replicate_metrics$condition == "SYNTHETIC A" &
      replicate_metrics$replicate == "SYNTHETIC R1" &
      replicate_metrics$metric_id == "ph3_2to4n_positivity_percent", , drop = FALSE
  ]
  expect_identical(first$value_percent, 20)
  expect_identical(first$technical_acquisition_count, 3L)
  expect_identical(first$finite_technical_acquisition_count, 2L)
  expect_identical(first$undefined_technical_acquisition_count, 1L)
  expect_identical(first$result_status, "ok_partial_undefined")
  expect_identical(first$source_acquisition_ids,
                   '["SYNTHETIC-A1","SYNTHETIC-A2","SYNTHETIC-A3"]')
  condition <- result$ph3_metrics_condition_summary
  first_condition <- condition[
    condition$condition == "SYNTHETIC A" &
      condition$metric_id == "ph3_2to4n_positivity_percent", , drop = FALSE
  ]
  expect_identical(first_condition$mean_percent, 30)
  expect_equal(first_condition$sd_percent, sqrt(200))
  expect_equal(first_condition$sem_percent, 10)
  expect_identical(first_condition$biological_replicate_count, 2L)
  expect_identical(first_condition$result_status, "ok")
})

test_that("SYNTHETIC zero, single-replicate, and undefined statuses are exact", {
  inputs <- synthetic_slice5_inputs()
  result <- derive_ph3_replicate_condition_tables(
    inputs$quantitation, inputs$manifest, inputs$provenance
  )
  condition <- result$ph3_metrics_condition_summary
  b <- condition[condition$condition == "SYNTHETIC B", , drop = FALSE]
  expect_true(all(b$mean_percent == 0))
  expect_true(all(is.na(b$sd_percent)))
  expect_true(all(is.na(b$sem_percent)))
  expect_true(all(b$biological_replicate_count == 2L))
  expect_true(all(b$finite_biological_replicate_count == 1L))
  expect_true(all(b$undefined_biological_replicate_count == 1L))
  expect_true(all(b$result_status == "ok_partial_undefined"))

  keep <- inputs$manifest$replicate != "SYNTHETIC R4"
  one_manifest <- inputs$manifest[keep, , drop = FALSE]
  one_quantitation <- lapply(inputs$quantitation, function(data) {
    data[data$prefix %in% one_manifest$prefix, , drop = FALSE]
  })
  one_provenance <- inputs$provenance
  one_provenance$ph3_export_manifests <- one_provenance$ph3_export_manifests[1:5]
  one_provenance$ph3_containment <- one_provenance$ph3_containment[
    one_provenance$ph3_containment$prefix %in% one_manifest$prefix, , drop = FALSE
  ]
  one <- derive_ph3_replicate_condition_tables(
    one_quantitation, one_manifest, one_provenance
  )
  one_b <- one$ph3_metrics_condition_summary[
    one$ph3_metrics_condition_summary$condition == "SYNTHETIC B", , drop = FALSE
  ]
  expect_true(all(one_b$mean_percent == 0))
  expect_true(all(one_b$result_status == "ok_single_biological_replicate"))

  undefined_manifest <- inputs$manifest[
    inputs$manifest$condition != "SYNTHETIC B" |
      inputs$manifest$replicate == "SYNTHETIC R4", , drop = FALSE
  ]
  undefined_quantitation <- lapply(inputs$quantitation, function(data) {
    data[data$prefix %in% undefined_manifest$prefix, , drop = FALSE]
  })
  undefined <- derive_ph3_replicate_condition_tables(
    undefined_quantitation, undefined_manifest,
    list(
      ph3_export_manifests = inputs$provenance$ph3_export_manifests[c(1:4, 6)],
      ph3_containment = inputs$provenance$ph3_containment[
        inputs$provenance$ph3_containment$prefix %in% undefined_manifest$prefix,
        , drop = FALSE
      ]
    )
  )
  undefined_b <- undefined$ph3_metrics_condition_summary[
    undefined$ph3_metrics_condition_summary$condition == "SYNTHETIC B",
    , drop = FALSE
  ]
  expect_true(all(is.na(undefined_b$mean_percent)))
  expect_true(all(undefined_b$result_status == "undefined_no_finite_values"))
})

test_that("SYNTHETIC Slice 5 schemas and ordering are exact", {
  inputs <- synthetic_slice5_inputs()
  result <- derive_ph3_replicate_condition_tables(
    inputs$quantitation, inputs$manifest, inputs$provenance
  )
  expect_identical(names(result$ph3_metrics_biological_replicate), c(
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
  ))
  expect_identical(names(result$ph3_phase_prevalence_condition_summary), c(
    "analysis_id", "condition", "condition_index", "aggregation_level",
    "metric_id", "phase_id", "phase_index", "interval_lower",
    "interval_upper", "lower_inclusive", "upper_inclusive", "mean_percent",
    "sd_percent", "sem_percent", "biological_replicate_count",
    "finite_biological_replicate_count", "undefined_biological_replicate_count",
    "result_status", "source_replicate_ids", "classification_schema_version",
    "output_schema_version", "aggregation_method_id", "positivity_method_id",
    "eligibility_method_id", "interval_method_id", "four_n_method_id",
    "sub_four_n_method_id", "containment_method_id", "config_digest"
  ))
  expect_identical(
    result$ph3_metrics_biological_replicate$metric_id[1:5],
    c("ph3_2to4n_positivity_percent", "ph3_within_4n_positivity_percent",
      "ph3_4n_positive_prevalence_within_2to4n_percent",
      "ph3_within_sub_4n_positivity_percent",
      "ph3_sub_4n_positive_prevalence_within_2to4n_percent")
  )
  expect_true(all(result$ph3_metrics_biological_replicate$output_schema_version ==
                    "ph3-1.0.0"))
})

test_that("SYNTHETIC missing duplicated mixed or remapped inputs fail closed", {
  inputs <- synthetic_slice5_inputs()
  missing <- inputs$quantitation
  missing$ph3_metrics_acquisition <- missing$ph3_metrics_acquisition[-1, ]
  expect_error(
    derive_ph3_replicate_condition_tables(
      missing, inputs$manifest, inputs$provenance
    ),
    "malformed_or_incomplete_slice4_table"
  )
  duplicated <- inputs$quantitation
  duplicated$ph3_phase_prevalence$phase_id[[2L]] <-
    duplicated$ph3_phase_prevalence$phase_id[[1L]]
  expect_error(
    derive_ph3_replicate_condition_tables(
      duplicated, inputs$manifest, inputs$provenance
    ),
    "duplicated_or_malformed_slice4_row"
  )
  mixed <- inputs$quantitation
  mixed$ph3_metrics_acquisition$analysis_id[[1L]] <- "SYNTHETIC-other"
  expect_error(
    derive_ph3_replicate_condition_tables(
      mixed, inputs$manifest, inputs$provenance
    ),
    "mixed_or_missing_provenance"
  )
  remapped <- inputs$quantitation
  remapped$ph3_phase_prevalence$condition[[1L]] <- "SYNTHETIC wrong"
  expect_error(
    derive_ph3_replicate_condition_tables(
      remapped, inputs$manifest, inputs$provenance
    ),
    "manifest_membership_mismatch"
  )
})

test_that("SYNTHETIC active provenance rejects coordinated identity mutation", {
  inputs <- synthetic_slice5_inputs()
  mutations <- list(
    acquisition_id = "SYNTHETIC-mutated-acquisition",
    sample_id = "SYNTHETIC-mutated-sample",
    export_operation_id = "SYNTHETIC-mutated-operation",
    input_manifest_key = "SYNTHETIC-mutated-manifest"
  )
  for (field in names(mutations)) {
    changed <- inputs$quantitation
    for (name in names(changed)) {
      selected <- changed[[name]]$acquisition_id == "SYNTHETIC-A1"
      changed[[name]][[field]][selected] <- mutations[[field]]
    }
    expect_error(
      derive_ph3_replicate_condition_tables(
        changed, inputs$manifest, inputs$provenance
      ),
      "active_provenance_membership_mismatch", info = field
    )
  }

  mixed <- inputs$quantitation
  mixed$ph3_phase_prevalence$sample_id[
    mixed$ph3_phase_prevalence$acquisition_id == "SYNTHETIC-A1"
  ] <- "SYNTHETIC-mutated-sample"
  expect_error(
    derive_ph3_replicate_condition_tables(
      mixed, inputs$manifest, inputs$provenance
    ),
    "mixed_or_missing_provenance"
  )
})

test_that("SYNTHETIC condition and condition-scoped replicate indices are bijective", {
  inputs <- synthetic_slice5_inputs()
  split_condition <- inputs$manifest
  split_condition$condition_index[[2L]] <- 99L
  expect_error(
    derive_ph3_replicate_condition_tables(
      inputs$quantitation, split_condition, inputs$provenance
    ),
    "unstable_manifest_group_index"
  )
  shared_condition_index <- inputs$manifest
  shared_condition_index$condition_index[
    shared_condition_index$condition == "SYNTHETIC B"
  ] <- 1L
  expect_error(
    derive_ph3_replicate_condition_tables(
      inputs$quantitation, shared_condition_index, inputs$provenance
    ),
    "unstable_manifest_group_index"
  )
  shared_replicate_index <- inputs$manifest
  shared_replicate_index$replicate_index[
    shared_replicate_index$replicate == "SYNTHETIC R2"
  ] <- 1L
  expect_error(
    derive_ph3_replicate_condition_tables(
      inputs$quantitation, shared_replicate_index, inputs$provenance
    ),
    "unstable_manifest_group_index"
  )
  split_replicate <- inputs$manifest
  split_replicate$replicate_index[[2L]] <- 99L
  expect_error(
    derive_ph3_replicate_condition_tables(
      inputs$quantitation, split_replicate, inputs$provenance
    ),
    "unstable_manifest_group_index"
  )
})

test_that("SYNTHETIC replicate identity may be reused across conditions", {
  inputs <- synthetic_slice5_inputs()
  reused_manifest <- inputs$manifest
  selected_manifest <- reused_manifest$prefix == "SYNTHETIC_prefix_5"
  reused_manifest$replicate[selected_manifest] <- "SYNTHETIC R1"
  reused_manifest$replicate_index[selected_manifest] <- 1L
  reused_quantitation <- inputs$quantitation
  for (name in names(reused_quantitation)) {
    selected <- reused_quantitation[[name]]$prefix == "SYNTHETIC_prefix_5"
    reused_quantitation[[name]]$replicate[selected] <- "SYNTHETIC R1"
    reused_quantitation[[name]]$replicate_index[selected] <- 1L
  }
  result <- derive_ph3_replicate_condition_tables(
    reused_quantitation, reused_manifest, inputs$provenance
  )
  condition_b <- result$ph3_metrics_biological_replicate[
    result$ph3_metrics_biological_replicate$condition == "SYNTHETIC B" &
      result$ph3_metrics_biological_replicate$replicate == "SYNTHETIC R1",
    , drop = FALSE
  ]
  expect_identical(nrow(condition_b), 5L)
  expect_true(all(condition_b$replicate_index == 1L))
})

test_that("SYNTHETIC stale containment acquisition rows fail closed", {
  inputs <- synthetic_slice5_inputs()
  stale <- inputs$provenance
  stale$ph3_containment <- rbind(
    stale$ph3_containment,
    data.frame(
      acquisition_id = "SYNTHETIC-stale", prefix = "SYNTHETIC-stale-prefix",
      child_population_key = "g1", containment_status = "validated",
      export_operation_id = "SYNTHETIC-stale-operation",
      manifest_digest = "SYNTHETIC-stale-manifest",
      stringsAsFactors = FALSE
    )
  )
  expect_error(
    derive_ph3_replicate_condition_tables(
      inputs$quantitation, inputs$manifest, stale
    ),
    "active_provenance_membership_mismatch"
  )
})

test_that("SYNTHETIC empty operation and manifest provenance fail closed", {
  inputs <- synthetic_slice5_inputs()
  for (field in c("export_operation_id", "input_manifest_key")) {
    changed <- inputs$quantitation
    for (name in names(changed)) changed[[name]][[field]] <- ""
    expect_error(
      derive_ph3_replicate_condition_tables(
        changed, inputs$manifest, inputs$provenance
      ),
      "malformed_or_incomplete_slice4_table", info = field
    )
  }

  empty_operation <- inputs
  empty_operation$provenance$ph3_export_manifests[[1L]]$export_operation_id <- ""
  empty_operation$provenance$ph3_containment$export_operation_id[
    empty_operation$provenance$ph3_containment$acquisition_id == "SYNTHETIC-A1"
  ] <- ""
  for (name in names(empty_operation$quantitation)) {
    selected <- empty_operation$quantitation[[name]]$acquisition_id ==
      "SYNTHETIC-A1"
    empty_operation$quantitation[[name]]$export_operation_id[selected] <- ""
  }
  expect_error(
    derive_ph3_replicate_condition_tables(
      empty_operation$quantitation, empty_operation$manifest,
      empty_operation$provenance
    ),
    "missing_active_analysis_provenance"
  )

  empty_digest <- inputs
  empty_digest$provenance$ph3_export_manifests[[1L]]$sha256 <- ""
  empty_digest$provenance$ph3_containment$manifest_digest[
    empty_digest$provenance$ph3_containment$acquisition_id == "SYNTHETIC-A1"
  ] <- ""
  for (name in names(empty_digest$quantitation)) {
    selected <- empty_digest$quantitation[[name]]$acquisition_id ==
      "SYNTHETIC-A1"
    empty_digest$quantitation[[name]]$input_manifest_key[selected] <- ""
  }
  expect_error(
    derive_ph3_replicate_condition_tables(
      empty_digest$quantitation, empty_digest$manifest,
      empty_digest$provenance
    ),
    "missing_active_analysis_provenance"
  )
})

test_that("SYNTHETIC eligibility QC must retain the exact Slice 4 sequence", {
  inputs <- synthetic_slice5_inputs()
  replaced <- inputs$quantitation
  replaced$ph3_event_eligibility_qc$reason[[1L]] <- "SYNTHETIC-replacement"
  expect_error(
    derive_ph3_replicate_condition_tables(
      replaced, inputs$manifest, inputs$provenance
    ),
    "invalid_slice4_eligibility_qc_sequence"
  )
  missing <- inputs$quantitation
  missing$ph3_event_eligibility_qc <-
    missing$ph3_event_eligibility_qc[-1L, , drop = FALSE]
  expect_error(
    derive_ph3_replicate_condition_tables(
      missing, inputs$manifest, inputs$provenance
    ),
    "malformed_or_incomplete_slice4_table"
  )
  extra <- inputs$quantitation
  extra$ph3_event_eligibility_qc <- rbind(
    extra$ph3_event_eligibility_qc,
    extra$ph3_event_eligibility_qc[1L, , drop = FALSE]
  )
  expect_error(
    derive_ph3_replicate_condition_tables(
      extra, inputs$manifest, inputs$provenance
    ),
    "malformed_or_incomplete_slice4_table"
  )
})
