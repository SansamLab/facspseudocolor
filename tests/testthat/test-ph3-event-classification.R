# SYNTHETIC operational fixtures only. These values encode classification
# arithmetic and contract failures, not experimental observations or biology.

synthetic_classification_config <- function() {
  config <- minimal_config("ph3")
  config$ph3_input_profile <- "production_direct_identity_v1"
  config$ph3_export_operation_dirs <- "/SYNTHETIC/operation"
  config$g1_x_range <- c(1000, 1200)
  config$s_phase_bins <- list(
    early = c(1200, 1400), mid = c(1400, 1600), late = c(1600, 1800)
  )
  config$g2m_x_range <- c(1800, 2000)
  config$dna_channel <- "DNA"
  config$target_channel <- "Target"
  config
}

synthetic_classification_inputs <- function(
    dna = c(999, 1000, 1000.1, 1799.9, 1800, 1999.9, 2000, 2001,
            NA, NaN, Inf, -Inf),
    target = seq_along(dna)
) {
  ids <- paste0("SYNTHETIC-A:event_index:", seq_along(dna) - 1L)
  bound <- paste(rep("a", 64), collapse = "")
  operation_id <- "SYNTHETIC-OP"
  table <- data.frame(
    acquisition_id = "SYNTHETIC-A",
    event_index = as.character(seq_along(dna) - 1L),
    event_identity = ids,
    identity_source = "flowkit_get_gate_events_index",
    identity_method_id = "flowkit_source_event_index_scoped_v1",
    identity_method_version = "1.0.0",
    duplicate_occurrence = "1",
    export_profile = "production_direct_identity_v1",
    export_operation_id = operation_id,
    export_manifest_digest = bound,
    DNA = dna,
    Target = target,
    target_raw = target,
    dna_norm = dna,
    stringsAsFactors = FALSE
  )
  g1_ids <- ids[c(2, 3)]
  positive_ids <- ids[c(2, 5, 9, 11)]
  subset_ids <- function(selected) table[match(selected, ids), , drop = FALSE]
  normalized <- list(
    data = table,
    g1 = subset_ids(g1_ids),
    ph3_positive = subset_ids(positive_ids)
  )
  validated_inputs <- list(
    complete = table, g1 = normalized$g1,
    ph3_positive = normalized$ph3_positive
  )
  containment <- do.call(rbind, lapply(c("g1", "ph3_positive"), function(key) {
    child <- normalized[[key]]
    data.frame(
      acquisition_id = "SYNTHETIC-A", prefix = "reference",
      child_population_key = key, containment_status = "validated",
      parent_row_count = nrow(table), export_operation_id = operation_id,
      child_row_count = nrow(child),
      child_unique_identity_count = length(unique(child$event_identity)),
      matched_child_count = nrow(child),
      manifest_digest = bound,
      identity_method_id = "flowkit_source_event_index_scoped_v1",
      identity_method_version = "1.0.0",
      identity_source = "flowkit_get_gate_events_index",
      containment_method_id = "exact_direct_identity_multiset_containment",
      manifest_reference = "/SYNTHETIC/export-manifest.json",
      stringsAsFactors = FALSE
    )
  }))
  export_manifests <- list(list(
    sha256 = bound, export_operation_id = operation_id,
    manifest = list(
      profile = "production_direct_identity_v1",
      manifest_binding = list(digest = bound),
      approval = list(positivity_method_id = "SYNTHETIC_POSITIVITY_METHOD"),
      geometry_overlay_status = "not_requested",
      acquisitions = list(list(
        acquisition_id = "SYNTHETIC-A", sample_id = "SYNTHETIC-A.fcs",
        prefix = "reference"
      ))
    )
  ))
  list(
    normalized = normalized, validated_inputs = validated_inputs,
    containment = containment,
    export_manifests = export_manifests,
    manifest_row = data.frame(prefix = "reference", stringsAsFactors = FALSE),
    ids = ids
  )
}

test_that("SYNTHETIC classification preserves identity and exact boundaries", {
  fixture <- synthetic_classification_inputs()
  result <- build_ph3_event_classification(
    fixture$normalized, fixture$validated_inputs,
    synthetic_classification_config(),
    fixture$manifest_row, fixture$containment, fixture$export_manifests
  )

  expect_identical(nrow(result), nrow(fixture$normalized$data))
  expect_identical(result$event_identity, fixture$ids)
  expect_identical(result$parent_row_number, seq_along(fixture$ids))
  expect_identical(which(result$eligible_2to4n), 2:7)
  expect_identical(which(result$sub_4n_member), 2:4)
  expect_identical(which(result$four_n_member), 5:7)
  expect_false(any(result$sub_4n_member & result$four_n_member))
  expect_identical(
    sum(result$sub_4n_member) + sum(result$four_n_member),
    sum(result$eligible_2to4n)
  )
  expect_identical(
    result$eligibility_exclusion_reason,
    c("below_b0", rep("none", 6), "above_b5", rep("dna_nonfinite", 4))
  )
  expect_true(all(result$unassigned_reason == "none"))
  expect_identical(result$configured_phase_id[[5]], "G2/M")
  expect_true(all(result$output_schema_version == "ph3-1.0.0"))
})

test_that("SYNTHETIC target finiteness never changes eligibility or positivity", {
  fixture <- synthetic_classification_inputs(
    dna = c(1000, 1800, 2000), target = c(NA, Inf, -Inf)
  )
  fixture$normalized$g1 <- fixture$normalized$data[1:2, , drop = FALSE]
  fixture$normalized$ph3_positive <- fixture$normalized$data[c(1, 3), , drop = FALSE]
  fixture$validated_inputs$g1 <- fixture$normalized$g1
  fixture$validated_inputs$ph3_positive <- fixture$normalized$ph3_positive
  fixture$containment$parent_row_count <- 3L
  fixture$containment$child_row_count <- c(2L, 2L)
  fixture$containment$child_unique_identity_count <- c(2L, 2L)
  fixture$containment$matched_child_count <- c(2L, 2L)
  result <- build_ph3_event_classification(
    fixture$normalized, fixture$validated_inputs,
    synthetic_classification_config(),
    fixture$manifest_row, fixture$containment, fixture$export_manifests
  )
  expect_true(all(result$eligible_2to4n))
  expect_identical(result$ph3_positive_member, c(TRUE, FALSE, TRUE))
  expect_false(any(result$target_display_finite))
})

test_that("SYNTHETIC configured phase gaps remain explicitly Unassigned", {
  fixture <- synthetic_classification_inputs(dna = 1500, target = 1)
  fixture$normalized$g1 <- fixture$normalized$data
  fixture$normalized$ph3_positive <- fixture$normalized$data
  fixture$validated_inputs$g1 <- fixture$normalized$g1
  fixture$validated_inputs$ph3_positive <- fixture$normalized$ph3_positive
  fixture$containment$parent_row_count <- 1L
  fixture$containment$child_row_count <- 1L
  fixture$containment$child_unique_identity_count <- 1L
  fixture$containment$matched_child_count <- 1L
  config <- synthetic_classification_config()
  config$s_phase_bins$mid <- c(1400, 1450)
  config$s_phase_bins$late <- c(1550, 1800)
  result <- build_ph3_event_classification(
    fixture$normalized, fixture$validated_inputs, config, fixture$manifest_row,
    fixture$containment, fixture$export_manifests
  )
  expect_true(result$eligible_2to4n)
  expect_false(result$configured_phase_assigned)
  expect_identical(result$configured_phase_id, "Unassigned")
  expect_identical(result$unassigned_reason, "configured_phase_gap")
})

test_that("SYNTHETIC identity and containment inconsistencies fail closed", {
  fixture <- synthetic_classification_inputs(dna = c(1000, 1800), target = 1:2)
  fixture$normalized$g1 <- fixture$normalized$data[1, , drop = FALSE]
  fixture$normalized$ph3_positive <- fixture$normalized$data[2, , drop = FALSE]
  fixture$validated_inputs$g1 <- fixture$normalized$g1
  fixture$validated_inputs$ph3_positive <- fixture$normalized$ph3_positive
  fixture$containment$parent_row_count <- 2L
  fixture$containment$child_row_count <- 1L
  fixture$containment$child_unique_identity_count <- 1L
  fixture$containment$matched_child_count <- 1L
  altered <- fixture
  altered$normalized$ph3_positive$event_identity <- "SYNTHETIC-A:event_index:99"
  expect_error(
    build_ph3_event_classification(
      altered$normalized, altered$validated_inputs,
      synthetic_classification_config(),
      altered$manifest_row, altered$containment, altered$export_manifests
    ),
    "ph3_positive_identity_inconsistency"
  )
  fixture$containment$containment_status[[1L]] <- "unverified"
  expect_error(
    build_ph3_event_classification(
      fixture$normalized, fixture$validated_inputs,
      synthetic_classification_config(),
      fixture$manifest_row, fixture$containment, fixture$export_manifests
    ),
    "invalid_containment_context"
  )
})

test_that("SYNTHETIC validated child loss and direct provenance mutation fail", {
  fixture <- synthetic_classification_inputs()
  dropped <- fixture
  dropped$normalized$ph3_positive <- dropped$normalized$ph3_positive[-1, , drop = FALSE]
  expect_error(
    build_ph3_event_classification(
      dropped$normalized, dropped$validated_inputs,
      synthetic_classification_config(), dropped$manifest_row,
      dropped$containment, dropped$export_manifests
    ),
    "ph3_positive_identity_inconsistency"
  )

  mutated <- fixture
  mutated$normalized$g1$event_index[[1L]] <- "01"
  mutated$validated_inputs$g1 <- mutated$normalized$g1
  expect_error(
    build_ph3_event_classification(
      mutated$normalized, mutated$validated_inputs,
      synthetic_classification_config(), mutated$manifest_row,
      mutated$containment, mutated$export_manifests
    ),
    "g1_direct_identity_mismatch"
  )

  provenance <- fixture
  provenance$normalized$ph3_positive$identity_source[[1L]] <-
    "SYNTHETIC_MUTATED_SOURCE"
  provenance$validated_inputs$ph3_positive <-
    provenance$normalized$ph3_positive
  expect_error(
    build_ph3_event_classification(
      provenance$normalized, provenance$validated_inputs,
      synthetic_classification_config(), provenance$manifest_row,
      provenance$containment, provenance$export_manifests
    ),
    "ph3_positive_identity_inconsistency"
  )

  profile <- fixture
  profile$normalized$g1$export_profile[[1L]] <- "SYNTHETIC_MUTATED_PROFILE"
  profile$validated_inputs$g1 <- profile$normalized$g1
  expect_error(
    build_ph3_event_classification(
      profile$normalized, profile$validated_inputs,
      synthetic_classification_config(), profile$manifest_row,
      profile$containment, profile$export_manifests
    ),
    "g1_direct_identity_mismatch"
  )
})

test_that("SYNTHETIC configuration identity excludes location spelling", {
  first <- synthetic_classification_config()
  second <- first
  first$data_dir <- "C:/SYNTHETIC/data"
  first$output_pdf <- "C:/SYNTHETIC/results/out.pdf"
  first$output_png <- "C:/SYNTHETIC/results/out.png"
  first$ph3_export_operation_dirs <- "C:/SYNTHETIC/operation"
  second$data_dir <- "/SYNTHETIC/data"
  second$output_pdf <- "/SYNTHETIC/results/out.pdf"
  second$output_png <- "/SYNTHETIC/results/out.png"
  second$ph3_export_operation_dirs <- "/SYNTHETIC/operation"
  expect_identical(ph3_configuration_digest(first),
                   ph3_configuration_digest(second))
  scientific_changes <- list(
    g2m_x_range = c(1800, 2100),
    g1_x_range = c(900, 1200),
    dna_channel = "SYNTHETIC_OTHER_DNA"
  )
  for (field in names(scientific_changes)) {
    changed <- second
    changed[[field]] <- scientific_changes[[field]]
    expect_false(
      identical(ph3_configuration_digest(first),
                ph3_configuration_digest(changed)),
      info = field
    )
  }
})
