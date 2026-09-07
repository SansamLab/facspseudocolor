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

test_that("model-derived EdU manifests are hash-bound and fail closed", {
  directory <- withr::local_tempdir(pattern = "SYNTHETIC_model_gate_")
  populations <- c(single_cells = "complete", g1 = "g1",
                   edu_positive = "edu_positive")
  paths <- file.path(directory, paste0("SYNTHETIC_A_", names(populations), ".csv"))
  for (path in paths) writeLines(c("DNA content,EdU", "1,2"), path)
  input_report <- data.frame(
    prefix = "SYNTHETIC_A", population = unname(populations), path = paths,
    exists = TRUE, stringsAsFactors = FALSE
  )
  analysis <- structure(list(
    config = list(plot_type = "edu"),
    sample_manifest = data.frame(prefix = "SYNTHETIC_A"),
    input_report = input_report
  ), class = "facs_analysis")
  outputs <- setNames(lapply(paths, function(path) list(
    path = basename(path),
    sha256 = paste0(as.character(openssl::sha256(file(path))), collapse = "")
  )), names(populations))
  make_manifest <- function(status = "EXPERIMENTAL MODEL-DERIVED NON-PRODUCTION",
                            prefix = "SYNTHETIC_A", records = outputs) {
    artifact_path <- file.path(directory, "SYNTHETIC-rule.json")
    if (!file.exists(artifact_path)) jsonlite::write_json(list(
      sha256 = "SYNTHETIC-semantic", schema_version = "fd-feature-v2",
      threshold = 0.5
    ), artifact_path, auto_unbox = TRUE)
    artifact <- list(path = artifact_path,
                     byte_sha256 = paste0(
                       as.character(openssl::sha256(file(artifact_path))),
                       collapse = ""
                     ),
                     semantic_sha256 = "SYNTHETIC-semantic")
    list(status_label = status, feature_schema_version = "fd-feature-v2",
         artifacts = list(single_cells = artifact, g1 = artifact),
         acquisitions = list(list(prefix = prefix, outputs = records)))
  }
  write_manifest <- function(value) {
    path <- tempfile("SYNTHETIC-manifest-", tmpdir = directory, fileext = ".json")
    jsonlite::write_json(value, path, auto_unbox = TRUE)
    path
  }
  expect_equal(nrow(validate_edu_model_gate_manifest(
    analysis, write_manifest(make_manifest()))), 3L)
  expect_error(validate_edu_model_gate_manifest(
    analysis, write_manifest(make_manifest(status = "production"))),
    "non-production status")
  expect_error(validate_edu_model_gate_manifest(
    analysis, write_manifest(make_manifest(prefix = "SYNTHETIC_OTHER"))),
    "prefixes do not exactly match")
  wrong_path <- outputs
  wrong_path$single_cells$path <- basename(paths[[2]])
  expect_error(validate_edu_model_gate_manifest(
    analysis, write_manifest(make_manifest(records = wrong_path))),
    "input path differs")
  wrong_hash <- outputs
  wrong_hash$g1$sha256 <- paste(rep("0", 64), collapse = "")
  expect_error(validate_edu_model_gate_manifest(
    analysis, write_manifest(make_manifest(records = wrong_hash))), "SHA-256")
  traversal <- outputs
  traversal$single_cells$path <- "../SYNTHETIC_A_single_cells.csv"
  expect_error(validate_edu_model_gate_manifest(
    analysis, write_manifest(make_manifest(records = traversal))),
    "path-free basename")
})
