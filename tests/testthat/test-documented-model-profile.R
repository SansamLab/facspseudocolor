test_that("documented model profile rejects incomplete provenance", {
  directory <- withr::local_tempdir(pattern = "SYNTHETIC_documented_model_")
  populations <- c(single_cells = "complete", g1 = "g1",
                   edu_positive = "edu_positive")
  paths <- file.path(directory, paste0("SYNTHETIC_A_", names(populations), ".csv"))
  for (path in paths) {
    writeLines(c("event_identity,event_index,DNA content,EdU",
                 "SYNTHETIC_A:event_index:0,0,1,2"), path)
  }
  analysis <- structure(list(
    config = list(
      plot_type = "edu",
      gating = list(mode = "model_experimental",
                    profile = "documented_edu_g1_v1")
    ),
    sample_manifest = data.frame(prefix = "SYNTHETIC_A"),
    input_report = data.frame(
      prefix = "SYNTHETIC_A", population = unname(populations), path = paths,
      exists = TRUE, stringsAsFactors = FALSE
    )
  ), class = "facs_analysis")
  manifest <- list(
    schema_version = 3L,
    status_label = "EXPERIMENTAL MODEL-DERIVED NON-PRODUCTION",
    predicted_parent_authorization = "SYNTHETIC-NOT-APPROVED",
    g1_source = "model", single_cells_source = "model",
    edu_positive_source = "model", flowjo_workspace_gate_used = FALSE,
    feature_scopes = list(), approved_mapping = list(sha256 = "SYNTHETIC"),
    artifacts = list(), acquisitions = list()
  )
  path <- file.path(directory, "manifest.json")
  jsonlite::write_json(manifest, path, auto_unbox = TRUE)
  expect_error(validate_edu_model_gate_manifest(analysis, path),
               "authorization")
})

test_that("documented model membership checks are profile-specific", {
  parent <- data.frame(event_identity = paste0("SYNTHETIC:event_index:", 0:4))
  g1 <- parent[c(1, 3, 5), , drop = FALSE]
  edu <- parent[c(2, 4), , drop = FALSE]
  expect_true(facspseudocolor:::validate_documented_model_parent_membership(
    parent, g1, edu
  ))
  outside <- g1
  outside$event_identity[[2]] <- "SYNTHETIC:event_index:99"
  expect_error(
    facspseudocolor:::validate_documented_model_parent_membership(
      parent, outside, edu
    ),
    "order-preserving subset"
  )
})

test_that("documented model QC and header checks remain hardened", {
  expect_identical(
    facspseudocolor:::format_documented_model_qc_warnings(
      "SYNTHETIC", character()
    ),
    character()
  )
  mapping_hash <- paste(rep("a", 64), collapse = "")
  manifest <- list(
    schema_version = 3L,
    status_label = "EXPERIMENTAL MODEL-DERIVED NON-PRODUCTION",
    predicted_parent_authorization =
      "OWNER_AUTHORIZED_PREDICTED_PARENT_EXPERIMENTAL_BRANCH_ONLY",
    single_cells_source = "model", g1_source = "model",
    edu_positive_source = "model", flowjo_workspace_gate_used = FALSE,
    approved_mapping = list(sha256 = mapping_hash),
    feature_scopes = list(
      single_cells = "complete_source_frame",
      g1 = "predicted_single_cells_parent",
      edu_positive = "complete_source_frame_then_predicted_single_cells_parent"
    )
  )
  expect_true(facspseudocolor:::validate_documented_model_manifest_header(
    manifest, mapping_hash
  ))
  manifest$feature_scopes$g1 <- "complete_source_frame"
  expect_error(
    facspseudocolor:::validate_documented_model_manifest_header(
      manifest, mapping_hash
    ),
    "feature scopes"
  )
})
