# Every constructed value in this file is explicitly SYNTHETIC and test-only.
# These fixtures exercise contract identity and validation, not biology.

synthetic_ph3_output_model_config <- function(reference = TRUE) {
  comparisons <- list(list(
    id = "SYNTHETIC-treatment-vs-control",
    control_condition_id = "SYNTHETIC-control",
    treatment_condition_id = "SYNTHETIC-treatment",
    outcome_ids = c("A", "B", "C", "D")
  ))
  list(
    plot_type = "ph3", ph3_input_profile = "production_direct_identity_v1",
    ph3_boundary_sensitivity_fraction = 0.1, dna_2n_value = 10,
    g1_x_range = c(0, 3),
    s_phase_bins = list(
      early = c(3, 6), mid = c(6, 9), late = c(9, 12)
    ),
    g2m_x_range = c(12, 15),
    replicates = list(
      list(
        id = "SYNTHETIC-replicate-set-1", label = "SYNTHETIC R1",
        samples = list(
          list(label = "SYNTHETIC Control", prefix = "SYNTHETIC_r1_control"),
          list(label = "SYNTHETIC Treatment", prefix = "SYNTHETIC_r1_treatment")
        )
      ),
      list(
        id = "SYNTHETIC-replicate-set-2", label = "SYNTHETIC R2",
        samples = list(
          list(label = "SYNTHETIC Control", prefix = "SYNTHETIC_r2_control"),
          list(label = "SYNTHETIC Treatment", prefix = "SYNTHETIC_r2_treatment")
        )
      )
    ),
    ph3_output_contract = list(
      schema_version = "ph3-output-contract-config-1.0.0",
      experiment_id = "SYNTHETIC-experiment",
      conditions = list(
        list(
          id = "SYNTHETIC-control", label = "SYNTHETIC Control",
          role = "control"
        ),
        list(
          id = "SYNTHETIC-treatment", label = "SYNTHETIC Treatment",
          role = "treatment"
        )
      ),
      reference = if (isTRUE(reference)) {
        list(status = "configured", condition_id = "SYNTHETIC-control")
      } else list(status = "not_configured"),
      comparisons = comparisons,
      geometry = list(
        verified = list(status = "not_configured"),
        computed = list(status = "not_computed")
      )
    )
  )
}

synthetic_ph3_output_model_analysis <- function(reference = TRUE) {
  config <- synthetic_ph3_output_model_config(reference)
  manifest <- do.call(rbind, lapply(seq_along(config$replicates), function(i) {
    replicate <- config$replicates[[i]]
    do.call(rbind, lapply(seq_along(replicate$samples), function(j) {
      sample <- replicate$samples[[j]]
      data.frame(
        replicate_set_id = replicate$id,
        replicate = replicate$label, replicate_index = as.integer(i),
        technical_replicate = "1", condition = sample$label,
        condition_index = as.integer(j), prefix = sample$prefix,
        stringsAsFactors = FALSE
      )
    }))
  }))
  rownames(manifest) <- NULL
  operations <- lapply(seq_len(nrow(manifest)), function(i) {
    acquisition_id <- paste0("SYNTHETIC-acquisition-", i)
    sample_id <- paste0("SYNTHETIC-sample-", i)
    operation_id <- paste0("SYNTHETIC-operation-", i)
    digest <- paste(rep(letters[[i]], 64L), collapse = "")
    list(
      sha256 = digest, export_operation_id = operation_id,
      manifest = list(
        approval = list(
          positivity_method_id =
            "flowjo_owner_approved_positive_population_v1"
        ),
        acquisitions = list(list(
          acquisition_id = acquisition_id, sample_id = sample_id,
          prefix = manifest$prefix[[i]]
        ))
      )
    )
  })
  containment <- do.call(rbind, lapply(seq_len(nrow(manifest)), function(i) {
    operation <- operations[[i]]
    acquisition <- operation$manifest$acquisitions[[1L]]
    data.frame(
      acquisition_id = acquisition$acquisition_id,
      prefix = acquisition$prefix,
      child_population_key = c("g1", "ph3_positive"),
      containment_status = "validated",
      containment_method_id = "exact_direct_identity_multiset_containment",
      export_operation_id = operation$export_operation_id,
      manifest_digest = operation$sha256, stringsAsFactors = FALSE
    )
  }))
  config_digest <- ph3_configuration_digest(config)
  analysis_id <- ph3_analysis_id(config_digest, operations)
  classifications <- lapply(seq_len(nrow(manifest)), function(i) {
    operation <- operations[[i]]
    acquisition <- operation$manifest$acquisitions[[1L]]
    dna <- c(0, 2, 4, 6, 8, 10, 12, 14)
    positive <- c(TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE)
    phases <- c(
      "G1", "G1", "Early S", "Mid S", "Mid S", "Late S", "G2/M", "G2/M"
    )
    data.frame(
      classification_schema_version = "ph3-event-classification-1.0.0",
      analysis_id = analysis_id,
      acquisition_id = acquisition$acquisition_id,
      sample_id = acquisition$sample_id, prefix = acquisition$prefix,
      event_identity = paste0(
        acquisition$acquisition_id, ":event_index:", seq_along(dna) - 1L
      ),
      identity_valid = TRUE, ph3_positive_member = positive,
      g1_containment_status = "validated",
      ph3_containment_status = "validated", dna_norm = dna,
      dna_raw = dna, target_raw = 2 + 3 * dna,
      dna_finite = TRUE, b0 = 0, b1 = 3, b2 = 6, b3 = 9, b4 = 12, b5 = 15,
      eligible_2to4n = TRUE, sub_4n_member = dna < 12,
      four_n_member = dna >= 12, configured_phase_id = phases,
      configured_phase_assigned = TRUE,
      eligibility_exclusion_reason = "none", unassigned_reason = "none",
      input_manifest_key = operation$sha256, config_digest = config_digest,
      export_operation_id = operation$export_operation_id,
      positivity_method_id = "flowjo_owner_approved_positive_population_v1",
      eligibility_method_id =
        "identity_validated_finite_dna_b0_b5_inclusive_v1",
      interval_method_id = "configured_shared_boundaries_left_closed_v1",
      four_n_method_id = "configured_b4_b5_closed_v1",
      sub_four_n_method_id =
        "configured_b0_b4_left_closed_right_open_v1",
      containment_method_id = "exact_direct_identity_multiset_containment",
      output_schema_version = "ph3-1.0.0", stringsAsFactors = FALSE
    )
  })
  acquisition_tables <- lapply(seq_len(nrow(manifest)), function(i) {
    derive_ph3_acquisition_tables(
      classifications[[i]], manifest[i, , drop = FALSE], config,
      operations, containment
    )
  })
  quantitation <- lapply(names(acquisition_tables[[1L]]), function(name) {
    value <- do.call(rbind, lapply(acquisition_tables, `[[`, name))
    rownames(value) <- NULL
    value
  })
  names(quantitation) <- names(acquisition_tables[[1L]])
  aggregate <- derive_ph3_replicate_condition_tables(
    quantitation, manifest,
    list(ph3_export_manifests = operations, ph3_containment = containment)
  )
  quantitation <- c(quantitation, aggregate)
  normalized <- lapply(classifications, function(classification) {
    list(ph3_event_classification = classification)
  })
  names(normalized) <- manifest$prefix
  structure(list(
    config = config, sample_manifest = manifest, normalized_data = normalized,
    quantitation = quantitation,
    provenance = list(
      ph3_export_manifests = operations, ph3_containment = containment
    )
  ), class = "facs_analysis")
}

test_that("SYNTHETIC output model has exact schemas and owner-confirmed order", {
  model <- build_ph3_output_model(synthetic_ph3_output_model_analysis())
  expect_s3_class(model, "ph3_output_model")
  expect_identical(names(model), c(
    "schema", "experiment", "outcomes", "conditions", "replicate_sets",
    "samples", "correction", "reference", "comparisons", "geometry",
    "output_registry", "source"
  ))
  expect_identical(model$schema$schema_version,
                   "ph3-output-contract-model-1.0.0")
  expect_identical(model$outcomes$outcome_id, c("A", "B", "C", "D"))
  expect_identical(names(model$outcomes), model$schema$outcome_columns)
  expect_identical(names(model$experiment), model$schema$experiment_columns)
  expect_identical(names(model$conditions), model$schema$condition_columns)
  expect_identical(names(model$replicate_sets),
                   model$schema$replicate_set_columns)
  expect_identical(names(model$samples), model$schema$sample_columns)
  expect_identical(names(model$correction), model$schema$correction_columns)
  expect_identical(names(model$reference), model$schema$reference_columns)
  expect_identical(names(model$comparisons), model$schema$comparison_columns)
  expect_identical(names(model$geometry), model$schema$geometry_columns)
  expect_identical(names(model$output_registry), model$schema$artifact_columns)
  expect_identical(model$replicate_sets$replicate_set_id, c(
    "SYNTHETIC-replicate-set-1", "SYNTHETIC-replicate-set-2"
  ))
  expect_true(all(model$correction$status == "not_computed"))
  expect_true(all(is.na(model$correction$signal_basis)))
  expect_identical(
    model$geometry$required_provenance_label, c("VERIFIED", "COMPUTED")
  )
  expect_identical(vapply(model$outcomes, typeof, character(1)),
                   model$schema$outcome_types)
  expect_identical(vapply(model$comparisons, typeof, character(1)),
                   model$schema$comparison_types)
  expect_true(all(model$output_registry$required))
  expect_identical(model$output_registry$path_template[1:5], c(
    "ph3_report.html", "sample_results.csv", "condition_summary.csv",
    "regression_qc.csv", "provenance_manifest.json"
  ))
})

test_that("SYNTHETIC membership and approved A/B values remain byte-identical", {
  analysis <- synthetic_ph3_output_model_analysis()
  membership_before <- serialize(lapply(
    analysis$normalized_data,
    function(x) x$ph3_event_classification[c(
      "event_identity", "ph3_positive_member", "eligible_2to4n",
      "sub_4n_member", "four_n_member"
    )]
  ), NULL, version = 3L)
  metric <- analysis$quantitation$ph3_metrics_acquisition
  ab <- metric$metric_id %in% c(
    "ph3_4n_positive_prevalence_within_2to4n_percent",
    "ph3_sub_4n_positive_prevalence_within_2to4n_percent"
  )
  ab_before <- serialize(metric[ab, , drop = FALSE], NULL, version = 3L)
  model <- build_ph3_output_model(analysis)
  membership_after <- serialize(lapply(
    model$source$event_classifications,
    function(x) x[c(
      "event_identity", "ph3_positive_member", "eligible_2to4n",
      "sub_4n_member", "four_n_member"
    )]
  ), NULL, version = 3L)
  source_metric <- model$source$quantitation$ph3_metrics_acquisition
  source_ab <- source_metric$metric_id %in% c(
    "ph3_4n_positive_prevalence_within_2to4n_percent",
    "ph3_sub_4n_positive_prevalence_within_2to4n_percent"
  )
  expect_identical(membership_after, membership_before)
  expect_identical(
    serialize(source_metric[source_ab, , drop = FALSE], NULL, version = 3L),
    ab_before
  )
  expect_identical(model$source$sample_manifest, analysis$sample_manifest)
  expect_identical(model$source$quantitation, analysis$quantitation)
})

test_that("SYNTHETIC reference and comparison state is explicit and exact", {
  model <- build_ph3_output_model(synthetic_ph3_output_model_analysis())
  expect_true(all(model$reference$status == "configured"))
  expect_identical(model$reference$condition_id,
                   rep("SYNTHETIC-control", 2L))
  expect_identical(model$comparisons$outcome_id, c("A", "B", "C", "D"))

  no_reference <- build_ph3_output_model(
    synthetic_ph3_output_model_analysis(reference = FALSE)
  )
  expect_true(all(no_reference$reference$status == "not_configured"))
  expect_true(all(is.na(no_reference$reference$sample_id)))
  expect_true(all(no_reference$reference$reason_code ==
                    "reference_not_configured"))
})

test_that("SYNTHETIC missing extra mixed and unstable contract state fails closed", {
  base <- synthetic_ph3_output_model_config()$ph3_output_contract
  missing <- base
  missing$geometry <- NULL
  expect_error(
    validate_ph3_output_contract_config(
      missing, synthetic_ph3_output_model_config()$replicates
    ), "invalid_schema"
  )
  extra <- base
  extra$undocumented <- TRUE
  expect_error(
    validate_ph3_output_contract_config(
      extra, synthetic_ph3_output_model_config()$replicates
    ), "invalid_schema"
  )
  wrong_direction <- base
  wrong_direction$comparisons[[1L]]$control_condition_id <-
    "SYNTHETIC-treatment"
  wrong_direction$comparisons[[1L]]$treatment_condition_id <-
    "SYNTHETIC-control"
  expect_error(
    validate_ph3_output_contract_config(
      wrong_direction, synthetic_ph3_output_model_config()$replicates
    ), "comparison_direction_mismatch"
  )
  duplicated <- base
  duplicated$comparisons <- c(duplicated$comparisons, duplicated$comparisons)
  expect_error(
    validate_ph3_output_contract_config(
      duplicated, synthetic_ph3_output_model_config()$replicates
    ), "duplicate_comparison"
  )
  unstable_config <- synthetic_ph3_output_model_config()
  unstable_config$replicates[[2L]]$id <- unstable_config$replicates[[1L]]$id
  expect_error(
    validate_ph3_output_contract_config(
      unstable_config$ph3_output_contract, unstable_config$replicates
    ), "invalid_replicate_set_identity"
  )
  unsafe_geometry <- base
  unsafe_geometry$geometry$verified <- list(
    status = "geometry_set_index_declared",
    geometry_set_index_id = "../SYNTHETIC-unsafe-path"
  )
  expect_error(
    validate_ph3_output_contract_config(
      unsafe_geometry, synthetic_ph3_output_model_config()$replicates
    ), "missing_geometry_set_index"
  )

  mixed <- synthetic_ph3_output_model_analysis()
  mixed$quantitation$ph3_metrics_acquisition$analysis_id[[1L]] <-
    "SYNTHETIC-mixed-analysis"
  expect_error(build_ph3_output_model(mixed), "slice4_source_mismatch")

  remapped <- synthetic_ph3_output_model_analysis()
  remapped$sample_manifest$replicate_set_id[
    remapped$sample_manifest$replicate_index == 2L
  ] <- "SYNTHETIC-remapped-replicate"
  expect_error(build_ph3_output_model(remapped),
               "replicate_set_manifest_mismatch")
})

test_that("SYNTHETIC validated output model is the Slice 2 boundary", {
  validated <- build_ph3_output_model(synthetic_ph3_output_model_analysis())
  expect_s3_class(validated, "ph3_output_model")
  result <- apply_ph3_background_regression(validated)
  expect_true(all(result$correction$status == "selected"))
  expect_true(all(result$correction$signal_basis == "raw"))
  expect_true(all(result$background_regression$set_decisions$signal_basis ==
                    "raw"))
  expect_identical(result$source$event_classifications,
                   validated$source$event_classifications)
  expect_identical(result$source$quantitation, validated$source$quantitation)
})
