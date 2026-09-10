test_that("high-level EdU analysis returns a structured non-saving result", {
  config <- read_facs_config(system.file(
    "config", "config_edu.yml", package = "facspseudocolor"
  ))
  expect_warning(
    result <- analyze_facs_experiment(config),
    "aliases.*deprecated"
  )

  expect_s3_class(result, "facs_analysis")
  expect_equal(nrow(result$sample_manifest), 6)
  expect_length(result$models, 6)
  expect_length(result$normalized_data, 6)
  expect_identical(names(result$normalized_data), result$sample_manifest$prefix)
  expect_equal(
    vapply(result$normalized_data, function(x) {
      stats::median(x$data$target_norm, na.rm = TRUE)
    }, numeric(1)),
    c(rep1_NT = 1752.82011390766, rep1_2h = 1504.51043111236,
      rep1_4h = 1368.51065076182, rep2_NT = 2349.03387094802,
      rep2_2h = 1643.1981974707, rep2_4h = 1472.14509202877),
    tolerance = 1e-8
  )
  expect_identical(names(result$quantitation$by_signal),
                   "background_subtracted")
  expect_true(all(c("phase_medians", "whole_medians") %in%
                    names(result$quantitation$by_signal$background_subtracted)))
  expect_true("phase_percentages" %in% names(result$quantitation))
  expect_identical(result$provenance$package, "facspseudocolor")
})

test_that("high-level POI analysis preserves current normalized results", {
  config <- read_facs_config(system.file(
    "config", "config_poi.yml", package = "facspseudocolor"
  ))
  result <- analyze_facs_experiment(config)

  expect_s3_class(result, "facs_analysis")
  expect_length(result$models, 2)
  expect_equal(
    vapply(result$normalized_data, function(x) {
      stats::median(x$data$target_norm, na.rm = TRUE)
    }, numeric(1)),
    c(rep1_bg = 956.742596862670, rep1_NT = 4991.65220845589,
      rep1_2h = 1638.55373071234, rep1_4h = 1826.48029208043,
      rep2_bg = 968.992295650424, rep2_NT = 4098.33453848040,
      rep2_2h = 1648.79195919279, rep2_4h = 1689.28981136966),
    tolerance = 1e-8
  )
})

test_that("analysis print method summarizes without exposing event data", {
  object <- structure(list(
    config = list(plot_type = "poi"),
    sample_manifest = data.frame(replicate_index = c(1, 1)),
    models = list(list()), warnings = character()
  ), class = "facs_analysis")
  expect_output(print(object), "mode:       poi")
  expect_output(print(object), "samples:    2")
})

test_that("model-only EdU manifests reject synthetic or incomplete provenance", {
  directory <- withr::local_tempdir(pattern = "SYNTHETIC_model_gate_")
  populations <- c(single_cells = "complete", g1 = "g1",
                   edu_positive = "edu_positive")
  paths <- file.path(directory, paste0("SYNTHETIC_A_", names(populations), ".csv"))
  for (path in paths) writeLines(c("event_identity,event_index,DNA content,EdU",
                                   "SYNTHETIC_A:event_index:0,0,1,2"), path)
  input_report <- data.frame(
    prefix = "SYNTHETIC_A", population = unname(populations), path = paths,
    exists = TRUE, stringsAsFactors = FALSE
  )
  analysis <- structure(list(
    config = list(plot_type = "edu", g1_source = "model"),
    sample_manifest = data.frame(prefix = "SYNTHETIC_A"),
    input_report = input_report
  ), class = "facs_analysis")
  make_manifest <- function(status = "EXPERIMENTAL MODEL-DERIVED NON-PRODUCTION")
    list(schema_version = 3L, status_label = status,
         predicted_parent_authorization = "SYNTHETIC-NOT-APPROVED",
         g1_source = "model", single_cells_source = "model",
         edu_positive_source = "model", flowjo_workspace_gate_used = FALSE,
         feature_scopes = list(), approved_mapping = list(sha256 = "SYNTHETIC"),
         artifacts = list(), acquisitions = list())
  write_manifest <- function(value) {
    path <- tempfile("SYNTHETIC-manifest-", tmpdir = directory, fileext = ".json")
    jsonlite::write_json(value, path, auto_unbox = TRUE)
    path
  }
  expect_error(validate_edu_model_gate_manifest(
    analysis, write_manifest(make_manifest())), "authorization")
  expect_error(validate_edu_model_gate_manifest(
    analysis, write_manifest(make_manifest(status = "production"))),
    "schema, authorization")
})

test_that("hardened synthetic model memberships validate and tampering fails", {
  parent <- data.frame(event_identity = paste0("SYNTHETIC:event_index:", 0:4))
  g1 <- parent[c(1, 3, 5), , drop = FALSE]
  edu <- parent[c(2, 4), , drop = FALSE]
  expect_true(facspseudocolor:::validate_model_parent_membership(parent, g1, edu))

  outside <- g1
  outside$event_identity[[2]] <- "SYNTHETIC:event_index:99"
  expect_error(facspseudocolor:::validate_model_parent_membership(
    parent, outside, edu), "order-preserving subset")

  reordered <- g1[c(3, 1, 2), , drop = FALSE]
  expect_error(facspseudocolor:::validate_model_parent_membership(
    parent, reordered, edu), "order-preserving subset")
})

test_that("empty model QC records remain empty", {
  expect_identical(
    facspseudocolor:::format_model_qc_warnings("SYNTHETIC", character()),
    character()
  )
  expect_identical(
    facspseudocolor:::format_model_qc_warnings("SYNTHETIC", "LOW_EVENT_SUPPORT"),
    "SYNTHETIC:LOW_EVENT_SUPPORT"
  )
})

test_that("valid hardened synthetic manifest header passes and tampering fails", {
  mapping_hash <- paste(rep("a", 64), collapse = "")
  manifest <- list(
    schema_version = 3L,
    status_label = "EXPERIMENTAL MODEL-DERIVED NON-PRODUCTION",
    predicted_parent_authorization =
      "OWNER_AUTHORIZED_PREDICTED_PARENT_EXPERIMENTAL_BRANCH_ONLY",
    single_cells_source = "model", g1_source = "model",
    edu_positive_source = "model", flowjo_workspace_gate_used = FALSE,
    approved_mapping = list(sha256 = mapping_hash),
    feature_scopes = list(single_cells = "complete_source_frame",
      g1 = "predicted_single_cells_parent",
      edu_positive = "complete_source_frame_then_predicted_single_cells_parent")
  )
  expect_true(facspseudocolor:::validate_model_manifest_header(
    manifest, mapping_hash))
  tampered <- manifest
  tampered$feature_scopes$g1 <- "complete_source_frame"
  expect_error(facspseudocolor:::validate_model_manifest_header(
    tampered, mapping_hash), "feature scopes")
  tampered <- manifest
  tampered$flowjo_workspace_gate_used <- TRUE
  expect_error(facspseudocolor:::validate_model_manifest_header(
    tampered, mapping_hash), "sources")
})
