# SYNTHETIC operational arithmetic only. These in-memory rows are not
# experimental observations and encode no biological conclusion.

synthetic_slice4_manifest <- function() {
  data.frame(
    replicate = "SYNTHETIC Replicate 1", replicate_index = 1L,
    technical_replicate = "1", condition = "SYNTHETIC condition",
    condition_index = 1L, prefix = "SYNTHETIC_acquisition",
    stringsAsFactors = FALSE
  )
}

synthetic_slice4_config <- function() {
  list(
    plot_type = "ph3", ph3_input_profile = "production_direct_identity_v1",
    ph3_boundary_sensitivity_fraction = 0.1, dna_2n_value = 10,
    g1_x_range = c(0, 3),
    s_phase_bins = list(
      early = c(3, 6), mid = c(6, 9), late = c(9, 12)
    ),
    g2m_x_range = c(12, 15)
  )
}

synthetic_slice4_export_manifests <- function() {
  digest <- paste(rep("a", 64), collapse = "")
  list(list(
    sha256 = digest, export_operation_id = "SYNTHETIC-operation",
    manifest = list(
      approval = list(
        positivity_method_id = "flowjo_owner_approved_positive_population_v1"
      ),
      acquisitions = list(list(
        acquisition_id = "SYNTHETIC-A", sample_id = "SYNTHETIC-A.fcs",
        prefix = "SYNTHETIC_acquisition"
      ))
    )
  ))
}

synthetic_slice4_containment <- function() {
  data.frame(
    acquisition_id = "SYNTHETIC-A",
    prefix = "SYNTHETIC_acquisition",
    child_population_key = c("g1", "ph3_positive"),
    containment_status = "validated",
    containment_method_id = "exact_direct_identity_multiset_containment",
    export_operation_id = "SYNTHETIC-operation",
    manifest_digest = paste(rep("a", 64), collapse = ""),
    stringsAsFactors = FALSE
  )
}

synthetic_slice4_classification <- function(
    dna = c(0, 2, 4, 6, 8, 10, 12, 14),
    positive = c(TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE),
    phases = c("G1", "G1", "Early S", "Mid S", "Mid S", "Late S",
               "G2/M", "G2/M"),
    eligible = rep(TRUE, length(dna)),
    exclusion = rep("none", length(dna)),
    unassigned = rep("none", length(dna))
) {
  stopifnot(length(dna) == length(positive), length(dna) == length(phases))
  n <- length(dna)
  b <- c(0, 3, 6, 9, 12, 15)
  assigned <- eligible & phases %in% c("G1", "Early S", "Mid S", "Late S", "G2/M")
  config <- synthetic_slice4_config()
  manifests <- synthetic_slice4_export_manifests()
  data.frame(
    classification_schema_version = "ph3-event-classification-1.0.0",
    analysis_id = ph3_analysis_id(ph3_configuration_digest(config), manifests),
    acquisition_id = "SYNTHETIC-A",
    sample_id = "SYNTHETIC-A.fcs", prefix = "SYNTHETIC_acquisition",
    event_identity = paste0("SYNTHETIC-A:event_index:", seq_len(n) - 1L),
    identity_valid = TRUE, ph3_positive_member = positive,
    g1_containment_status = "validated",
    ph3_containment_status = "validated", dna_norm = dna,
    dna_finite = is.finite(dna), b0 = b[[1L]], b1 = b[[2L]],
    b2 = b[[3L]], b3 = b[[4L]], b4 = b[[5L]], b5 = b[[6L]],
    eligible_2to4n = eligible,
    sub_4n_member = eligible & dna < b[[5L]],
    four_n_member = eligible & dna >= b[[5L]],
    configured_phase_id = phases, configured_phase_assigned = assigned,
    eligibility_exclusion_reason = exclusion, unassigned_reason = unassigned,
    input_manifest_key = manifests[[1L]]$sha256,
    config_digest = ph3_configuration_digest(config),
    export_operation_id = "SYNTHETIC-operation",
    positivity_method_id = "flowjo_owner_approved_positive_population_v1",
    eligibility_method_id = "identity_validated_finite_dna_b0_b5_inclusive_v1",
    interval_method_id = "configured_shared_boundaries_left_closed_v1",
    four_n_method_id = "configured_b4_b5_closed_v1",
    sub_four_n_method_id = "configured_b0_b4_left_closed_right_open_v1",
    containment_method_id = "exact_direct_identity_multiset_containment",
    output_schema_version = "ph3-1.0.0", stringsAsFactors = FALSE
  )
}

derive_synthetic_slice4 <- function(
    classification, config = synthetic_slice4_config()
) {
  derive_ph3_acquisition_tables(
    classification, synthetic_slice4_manifest(), config,
    synthetic_slice4_export_manifests(), synthetic_slice4_containment()
  )
}

test_that("SYNTHETIC five acquisition metrics use the approved denominators", {
  classification <- synthetic_slice4_classification()
  result <- derive_synthetic_slice4(classification)$ph3_metrics_acquisition
  expect_identical(result$metric_id, c(
    "ph3_2to4n_positivity_percent",
    "ph3_within_4n_positivity_percent",
    "ph3_4n_positive_prevalence_within_2to4n_percent",
    "ph3_within_sub_4n_positivity_percent",
    "ph3_sub_4n_positive_prevalence_within_2to4n_percent"
  ))
  expect_identical(result$numerator_count, c(4L, 0L, 0L, 4L, 4L))
  expect_identical(result$denominator_count, c(8L, 2L, 8L, 6L, 8L))
  expect_equal(result$value_percent, c(50, 0, 0, 100 * 4 / 6, 50))
  expect_true(all(result$result_status == "ok"))
  expect_false(result$upper_inclusive[[4L]])
  expect_false(result$upper_inclusive[[5L]])
})

test_that("SYNTHETIC zero and missing denominator semantics are exact", {
  zero_positive <- synthetic_slice4_classification(positive = rep(FALSE, 8))
  result <- derive_synthetic_slice4(zero_positive)$ph3_metrics_acquisition
  expect_true(all(result$value_percent == 0))
  expect_true(all(result$result_status == "ok"))

  no_four <- synthetic_slice4_classification(dna = 0:7, positive = rep(FALSE, 8),
    phases = c(rep("G1", 3), rep("Early S", 3), rep("Mid S", 2)))
  result <- derive_synthetic_slice4(no_four)$ph3_metrics_acquisition
  expect_true(is.na(result$value_percent[[2L]]))
  expect_identical(result$result_status[[2L]], "undefined_zero_denominator")

  no_sub <- synthetic_slice4_classification(dna = 12:14, positive = rep(FALSE, 3),
    phases = rep("G2/M", 3))
  result <- derive_synthetic_slice4(no_sub)$ph3_metrics_acquisition
  expect_true(is.na(result$value_percent[[4L]]))
  expect_identical(result$result_status[[4L]], "undefined_zero_denominator")

  no_eligible <- synthetic_slice4_classification(
    dna = c(-1, 16), positive = c(TRUE, TRUE), phases = c("Unassigned", "Unassigned"),
    eligible = c(FALSE, FALSE), exclusion = c("below_b0", "above_b5")
  )
  result <- derive_synthetic_slice4(no_eligible)$ph3_metrics_acquisition
  expect_true(all(is.na(result$value_percent)))
  expect_true(all(result$result_status == "undefined_zero_denominator"))
  expect_true(all(result$numerator_count == 0L))
})

test_that("SYNTHETIC positives outside eligibility never enter metrics", {
  classification <- synthetic_slice4_classification(
    dna = c(-1, 0, 12, 15, 16), positive = rep(TRUE, 5),
    phases = c("Unassigned", "G1", "G2/M", "G2/M", "Unassigned"),
    eligible = c(FALSE, TRUE, TRUE, TRUE, FALSE),
    exclusion = c("below_b0", "none", "none", "none", "above_b5")
  )
  result <- derive_synthetic_slice4(classification)$ph3_metrics_acquisition
  expect_identical(result$numerator_count, c(3L, 2L, 2L, 1L, 1L))
  expect_identical(result$denominator_count, c(3L, 2L, 3L, 1L, 3L))
})

test_that("SYNTHETIC phase prevalence and eligibility QC reconcile", {
  classification <- synthetic_slice4_classification()
  tables <- derive_synthetic_slice4(classification)
  phase <- tables$ph3_phase_prevalence
  expect_identical(phase$phase_id, c("G1", "Early S", "Mid S", "Late S", "G2/M"))
  expect_identical(phase$ph3_positive_count, c(1L, 1L, 1L, 1L, 0L))
  expect_true(all(phase$eligible_2to4n_count == 8L))
  expect_equal(phase$value_percent, c(12.5, 12.5, 12.5, 12.5, 0))
  qc <- tables$ph3_event_eligibility_qc
  expect_identical(qc$event_count[qc$reason == "imported_classification_rows"], 8L)
  expect_identical(qc$event_count[qc$reason == "eligible_2to4n"], 8L)
  expect_identical(qc$event_count[qc$reason == "sub_4n"], 6L)
  expect_identical(qc$event_count[qc$reason == "4n"], 2L)
  expect_identical(qc$event_count[qc$reason == "identity_valid"], 8L)
  expect_identical(
    qc$event_count[qc$reason == "configured_phase_assigned"], 8L
  )
  expect_identical(qc$event_count[
    qc$qc_dimension == "eligible_configured_phase_unassigned" &
      qc$reason == "none"
  ], 0L)
})

test_that("SYNTHETIC sensitivity has exactly four closed in-span variants", {
  classification <- synthetic_slice4_classification(
    dna = c(11, 12, 13, 14, 15), positive = c(TRUE, TRUE, TRUE, FALSE, TRUE),
    phases = c("Late S", rep("G2/M", 4))
  )
  tables <- derive_synthetic_slice4(classification)
  sensitivity <- tables$ph3_4n_boundary_sensitivity_qc
  expect_identical(sensitivity$perturbation_id,
                   c("primary", "lower_outward", "lower_inward", "upper_inward"))
  expect_identical(sensitivity$interval_lower, c(12, 11, 13, 12))
  expect_identical(sensitivity$interval_upper, c(15, 15, 15, 14))
  expect_true(all(sensitivity$lower_inclusive & sensitivity$upper_inclusive))
  expect_false(any(grepl("upper.*outward", sensitivity$perturbation_id)))
  expect_identical(sensitivity$eligible_regional_count, c(4L, 5L, 3L, 3L))
  primary_within <- tables$ph3_metrics_acquisition[
    tables$ph3_metrics_acquisition$metric_id ==
      "ph3_within_4n_positivity_percent", , drop = FALSE
  ]
  expect_identical(sensitivity$within_region_positivity_percent[[1L]],
                   primary_within$value_percent[[1L]])
})

test_that("SYNTHETIC invalid sensitivity and classification fail closed", {
  config <- synthetic_slice4_config()
  config$ph3_boundary_sensitivity_fraction <- 4
  invalid_delta <- synthetic_slice4_classification()
  invalid_delta$config_digest <- ph3_configuration_digest(config)
  invalid_delta$analysis_id <- ph3_analysis_id(
    invalid_delta$config_digest[[1L]], synthetic_slice4_export_manifests()
  )
  expect_error(
    derive_synthetic_slice4(invalid_delta, config),
    "invalid_boundary_sensitivity_interval"
  )
  expect_error(
    derive_synthetic_slice4(NULL),
    "missing_or_malformed_classification"
  )
  duplicated <- synthetic_slice4_classification()
  duplicated$event_identity[[2L]] <- duplicated$event_identity[[1L]]
  expect_error(
    derive_synthetic_slice4(duplicated), "missing_or_malformed_classification"
  )
  malformed <- synthetic_slice4_classification()
  malformed$four_n_member[[1L]] <- TRUE
  expect_error(
    derive_synthetic_slice4(malformed), "classification_reconciliation_failure"
  )
})

test_that("SYNTHETIC active provenance and Slice 3 predicates fail closed", {
  provenance <- synthetic_slice4_classification()
  provenance$config_digest <- "SYNTHETIC-wrong-config"
  expect_error(derive_synthetic_slice4(provenance), "active_provenance_mismatch")

  manifest <- synthetic_slice4_classification()
  manifest$input_manifest_key <- paste(rep("b", 64), collapse = "")
  expect_error(derive_synthetic_slice4(manifest), "active_provenance_mismatch")

  cases <- list(
    identity_valid = function(x) { x$identity_valid[[1L]] <- FALSE; x },
    dna_finite = function(x) { x$dna_finite[[1L]] <- FALSE; x },
    eligible = function(x) { x$eligible_2to4n[[1L]] <- FALSE; x },
    sub_four = function(x) { x$sub_4n_member[[1L]] <- FALSE; x },
    phase = function(x) { x$configured_phase_id[[1L]] <- "Early S"; x },
    phase_assigned = function(x) {
      x$configured_phase_assigned[[1L]] <- FALSE; x
    },
    exclusion = function(x) {
      x$eligibility_exclusion_reason[[1L]] <- "below_b0"; x
    },
    unassigned = function(x) {
      x$unassigned_reason[[1L]] <- "SYNTHETIC_UNAPPROVED_REASON"; x
    }
  )
  for (name in names(cases)) {
    expect_error(
      derive_synthetic_slice4(cases[[name]](synthetic_slice4_classification())),
      "classification_reconciliation_failure", info = name
    )
  }
  expect_error(
    derive_ph3_acquisition_tables(
      synthetic_slice4_classification(), synthetic_slice4_manifest(),
      synthetic_slice4_config(), list(), synthetic_slice4_containment()
    ), "active_provenance_mismatch"
  )

  identity_tamper <- synthetic_slice4_classification()
  identity_tamper$identity_valid[[1L]] <- FALSE
  identity_tamper$eligible_2to4n[[1L]] <- FALSE
  identity_tamper$sub_4n_member[[1L]] <- FALSE
  identity_tamper$configured_phase_assigned[[1L]] <- FALSE
  identity_tamper$eligibility_exclusion_reason[[1L]] <- "identity_invalid"
  expect_error(
    derive_synthetic_slice4(identity_tamper),
    "classification_reconciliation_failure"
  )

  finite_tamper <- synthetic_slice4_classification()
  finite_tamper$dna_finite[[1L]] <- FALSE
  finite_tamper$eligible_2to4n[[1L]] <- FALSE
  finite_tamper$sub_4n_member[[1L]] <- FALSE
  finite_tamper$configured_phase_assigned[[1L]] <- FALSE
  finite_tamper$eligibility_exclusion_reason[[1L]] <- "dna_nonfinite"
  expect_error(
    derive_synthetic_slice4(finite_tamper),
    "classification_reconciliation_failure"
  )

  containment <- synthetic_slice4_containment()
  containment$containment_method_id <- "SYNTHETIC_WRONG_CONTAINMENT"
  expect_error(
    derive_ph3_acquisition_tables(
      synthetic_slice4_classification(), synthetic_slice4_manifest(),
      synthetic_slice4_config(), synthetic_slice4_export_manifests(),
      containment
    ), "active_provenance_mismatch"
  )

  containment_mutations <- list(
    prefix = "SYNTHETIC_wrong_prefix",
    export_operation_id = "SYNTHETIC-wrong-operation",
    manifest_digest = paste(rep("c", 64), collapse = "")
  )
  for (field in names(containment_mutations)) {
    changed_containment <- synthetic_slice4_containment()
    changed_containment[[field]] <- containment_mutations[[field]]
    expect_error(
      derive_ph3_acquisition_tables(
        synthetic_slice4_classification(), synthetic_slice4_manifest(),
        synthetic_slice4_config(), synthetic_slice4_export_manifests(),
        changed_containment
      ), "active_provenance_mismatch", info = field
    )
  }
  missing_field <- synthetic_slice4_containment()
  missing_field$manifest_digest <- NULL
  expect_error(
    derive_ph3_acquisition_tables(
      synthetic_slice4_classification(), synthetic_slice4_manifest(),
      synthetic_slice4_config(), synthetic_slice4_export_manifests(),
      missing_field
    ), "active_provenance_mismatch"
  )
})

test_that("SYNTHETIC production attachment exposes only Slice 4 tables", {
  analysis <- structure(list(
    config = synthetic_slice4_config(),
    sample_manifest = synthetic_slice4_manifest(),
    normalized_data = list(SYNTHETIC_acquisition = list(
      ph3_event_classification = synthetic_slice4_classification()
    )), quantitation = list(), models = list(), warnings = character(),
    provenance = list(
      ph3_export_manifests = synthetic_slice4_export_manifests(),
      ph3_containment = synthetic_slice4_containment()
    ), input_report = data.frame()
  ), class = "facs_analysis")
  result <- quantify_ph3_production_acquisitions(analysis)
  expect_identical(names(result$quantitation), c(
    "ph3_metrics_acquisition", "ph3_phase_prevalence",
    "ph3_event_eligibility_qc", "ph3_4n_boundary_sensitivity_qc"
  ))
  expect_false(any(grepl(
    "biological|aggregate|csv|plot|report|geometry",
    names(result$quantitation), ignore.case = TRUE
  )))
  expect_identical(result$normalized_data[[1L]]$ph3_event_classification,
                   analysis$normalized_data[[1L]]$ph3_event_classification)
})
