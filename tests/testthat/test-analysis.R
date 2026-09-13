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
  writeLines(c("acquisition_id,event_identity,event_index,DNA content,EdU",
               "SYNTHETIC_A,SYNTHETIC_A:event_index:1,1,1,2",
               "SYNTHETIC_A,SYNTHETIC_A:event_index:2,2,2,3"), paths[[1]])
  writeLines(c("acquisition_id,event_identity,event_index,DNA content,EdU",
               "SYNTHETIC_A,SYNTHETIC_A:event_index:1,1,1,2"),
             paths[[2]])
  writeLines(c("acquisition_id,event_identity,event_index,DNA content,EdU",
               "SYNTHETIC_A,SYNTHETIC_A:event_index:2,2,2,3"),
             paths[[3]])
  source_fcs <- file.path(directory, "SYNTHETIC_A.fcs")
  writeLines("SYNTHETIC", source_fcs)
  single_path <- file.path(directory, "SYNTHETIC-single.json")
  g1_path <- file.path(directory, "SYNTHETIC-g1.json")
  edu_path <- file.path(directory, "SYNTHETIC-edu.portable.json")
  jsonlite::write_json(list(sha256 = "single-semantic", schema_version = "fd-feature-v2",
                            threshold = 0.205, target = "single_cells_reference"),
                       single_path, auto_unbox = TRUE)
  writeLines(paste0(
    '{"sha256":"g1-semantic","schema_version":"fd-feature-v2",',
    '"threshold":0.25499999999999995,"target":"g1_reference"}'
  ), g1_path)
  expect_match(readLines(g1_path), '"threshold":0.25499999999999995',
               fixed = TRUE)
  writeLines("SYNTHETIC", edu_path)
  artifact_hashes <- c(single_cells = "single-byte", g1 = "g1-byte",
                       edu_positive = "9e6ed466687efe5f7523dea633e9e9244381c0e613c86f0944f861c06047fdf4")
  artifacts <- list(
    single_cells = list(path = single_path, byte_sha256 = artifact_hashes[[1]],
                        semantic_sha256 = "single-semantic"),
    g1 = list(path = g1_path, byte_sha256 = artifact_hashes[[2]],
              semantic_sha256 = "g1-semantic"),
    edu_positive = list(path = edu_path, byte_sha256 = artifact_hashes[[3]]))
  roles <- list(fsc_area = "FSC-A", ssc_area = "SSC-A", dna_area = "FL2-A",
                dna_height = "FL2-H", edu = "FL4-A")
  acquisition_config <- list(prefix = "SYNTHETIC_A", acquisition_id = "SYNTHETIC_A",
                             fcs_path = source_fcs, fcs_sha256 = "source-hash",
                             channel_roles = roles)
  input_report <- data.frame(
    prefix = "SYNTHETIC_A", population = unname(populations), path = paths,
    exists = TRUE, stringsAsFactors = FALSE
  )
  analysis <- structure(list(
    config = list(plot_type = "edu", gating = list(
      mode = "model_experimental", profile = "single_g1_edu_frozen_v1",
      output_dir = directory, artifacts = artifacts,
      acquisitions = list(acquisition_config))),
    sample_manifest = data.frame(prefix = "SYNTHETIC_A"),
    input_report = input_report
  ), class = "facs_analysis")
  all_events_path <- file.path(directory, "SYNTHETIC_A_all_events.csv")
  writeLines(c(paste(
                 "acquisition_id,event_identity,event_index,DNA content,EdU",
                 "model_single_cells_probability,model_g1_probability",
                 "model_edu_positive_probability", sep = ","),
               "SYNTHETIC_A,SYNTHETIC_A:event_index:0,0,1,1,0.1,0.9,0.1",
               "SYNTHETIC_A,SYNTHETIC_A:event_index:1,1,1,2,0.9,0.9,0.1",
               "SYNTHETIC_A,SYNTHETIC_A:event_index:2,2,2,3,0.9,0.9,0.9"),
             all_events_path)
  outputs <- setNames(lapply(paths, function(path) list(
    path = basename(path),
    sha256 = paste0(as.character(openssl::sha256(file(path))), collapse = ""),
    rows = length(readLines(path)) - 1L
  )), names(populations))
  outputs$all_events <- list(
    path = basename(all_events_path),
    sha256 = paste0(as.character(openssl::sha256(file(all_events_path))),
                    collapse = ""),
    rows = 3L
  )
  make_manifest <- function(records = outputs) {
    list(status_label = "EXPERIMENTAL MODEL-DERIVED NON-PRODUCTION",
         gating_mode = "model_experimental", profile = "single_g1_edu_frozen_v1",
         feature_schema_version = "fd-feature-v2",
         frozen_thresholds = list(single_cells = 0.205, g1 = 0.255,
                                  edu_positive = 0.385),
         g1_dna_minimum = list(
           method = "median_positive_dna_area_of_preliminary_edu_negative_g1",
           fraction = 0.35),
         edu_predictor = list(
           format = "PORTABLE_HGB_V1",
           implementation = "embedded dependency-free reference",
           derived_from_reference_sha256 = "2edf1bbd529159e8752d01731e87f19efbc63b514409b1d246d7a4ffc266a2c3"),
         artifacts = artifacts,
         acquisitions = list(list(prefix = "SYNTHETIC_A", acquisition_id = "SYNTHETIC_A",
           source_fcs = source_fcs, source_fcs_sha256 = "source-hash",
           workspace_used = FALSE,
           channel_roles = roles[c("dna_area", "dna_height", "edu",
                                   "fsc_area", "ssc_area")],
           audit = list(model_single_cells_count = 2L, model_g1_count = 1L,
                        model_edu_positive_count = 1L, g1_subset_single_cells = TRUE,
                        all_event_count = 3L,
                        edu_positive_subset_single_cells = TRUE,
                        g1_excludes_edu_positive = TRUE,
                        preliminary_g1_count = 1L,
                        g1_dna_positive_candidate_count = 1L,
                        g1_dna_2n_center_raw = 1,
                        g1_dna_minimum_fraction = 0.35,
                        g1_dna_minimum_raw = 0.35,
                        g1_dna_minimum_excluded_count = 0L), outputs = records)))
  }
  write_manifest <- function(value) {
    path <- tempfile("SYNTHETIC-manifest-", tmpdir = directory, fileext = ".json")
    jsonlite::write_json(value, path, auto_unbox = TRUE)
    path
  }
  real_hash <- facs_sha256_file
  source_observed <- "source-hash"
  fake_hash <- function(path) {
    normalized <- normalizePath(path, mustWork = TRUE)
    if (normalized == normalizePath(single_path)) return("single-byte")
    if (normalized == normalizePath(g1_path)) return("g1-byte")
    if (normalized == normalizePath(edu_path)) return(artifact_hashes[[3]])
    if (normalized == normalizePath(source_fcs)) return(source_observed)
    real_hash(path)
  }
  testthat::local_mocked_bindings(
    facs_sha256_file = fake_hash, .package = "facspseudocolor"
  )
  expect_equal(nrow(validate_edu_model_gate_manifest(
    analysis, write_manifest(make_manifest()))), 3L)
  without_all <- make_manifest(outputs[names(outputs) != "all_events"])
  expect_error(validate_edu_model_gate_manifest(
    analysis, write_manifest(without_all)), "outputs are invalid")
  reordered <- make_manifest()
  reordered$frozen_thresholds <- list(
    edu_positive = 0.385, g1 = 0.255, single_cells = 0.205
  )
  expect_equal(nrow(validate_edu_model_gate_manifest(
    analysis, write_manifest(reordered))), 3L)
  reordered_predictor <- make_manifest()
  reordered_predictor$edu_predictor <- reordered_predictor$edu_predictor[
    c("derived_from_reference_sha256", "format", "implementation")
  ]
  expect_equal(nrow(validate_edu_model_gate_manifest(
    analysis, write_manifest(reordered_predictor))), 3L)
  with_all <- make_manifest()
  validated_all <- validate_edu_model_gate_manifest(
    analysis, write_manifest(with_all)
  )
  expect_identical(names(attr(validated_all, "all_events")), "SYNTHETIC_A")
  expect_identical(nrow(attr(validated_all, "all_events")[[1L]]), 3L)
  bad_all_hash <- with_all
  bad_all_hash$acquisitions[[1L]]$outputs$all_events$sha256 <-
    paste(rep("0", 64L), collapse = "")
  expect_error(validate_edu_model_gate_manifest(
    analysis, write_manifest(bad_all_hash)), "all-events path or SHA-256")
  bad_all_count <- with_all
  bad_all_count$acquisitions[[1L]]$audit$all_event_count <- 2L
  expect_error(validate_edu_model_gate_manifest(
    analysis, write_manifest(bad_all_count)), "identity/count audit")
  duplicated_predictor <- make_manifest()
  duplicated_predictor$edu_predictor <- structure(
    c(duplicated_predictor$edu_predictor, list("changed")),
    names = c(names(duplicated_predictor$edu_predictor), "format")
  )
  expect_error(validate_edu_model_gate_manifest(
    analysis, write_manifest(duplicated_predictor)),
    "portable EdU predictor provenance differs")
  duplicated <- make_manifest()
  duplicated$frozen_thresholds <- structure(
    list(0.205, 0.255, 0.385, 0.5),
    names = c("single_cells", "g1", "edu_positive", "single_cells")
  )
  expect_error(validate_edu_model_gate_manifest(
    analysis, write_manifest(duplicated)), "frozen thresholds differ")
  changed_threshold <- make_manifest()
  changed_threshold$frozen_thresholds$single_cells <- 0.206
  expect_error(validate_edu_model_gate_manifest(
    analysis, write_manifest(changed_threshold)), "frozen thresholds differ")
  bad <- make_manifest()
  bad$acquisitions[[1]]$source_fcs_sha256 <- "changed"
  expect_error(validate_edu_model_gate_manifest(analysis, write_manifest(bad)),
               "acquisition mapping differs")
  source_observed <- "changed-source-bytes"
  expect_error(validate_edu_model_gate_manifest(
    analysis, write_manifest(make_manifest())), "source FCS SHA-256 differs")
  source_observed <- "source-hash"
  writeLines(c("acquisition_id,event_identity,event_index,DNA content,EdU",
               "SYNTHETIC_A,SYNTHETIC_A:event_index:3,3,1,2"), paths[[2]])
  changed <- make_manifest()
  changed$acquisitions[[1]]$outputs$g1$sha256 <- real_hash(paths[[2]])
  expect_error(validate_edu_model_gate_manifest(analysis, write_manifest(changed)),
               "not contained")
  writeLines(c("acquisition_id,event_identity,event_index,DNA content,EdU",
               "SYNTHETIC_A,SYNTHETIC_A:event_index:2,2,2,3"), paths[[2]])
  overlapping <- make_manifest()
  overlapping$acquisitions[[1]]$outputs$g1$sha256 <- real_hash(paths[[2]])
  expect_error(validate_edu_model_gate_manifest(
    analysis, write_manifest(overlapping)), "child identities")
  writeLines(c("acquisition_id,event_identity,event_index,DNA content,EdU",
               "SYNTHETIC_A,noncanonical,1,1,2"), paths[[2]])
  noncanonical <- make_manifest()
  noncanonical$acquisitions[[1]]$outputs$g1$sha256 <- real_hash(paths[[2]])
  expect_error(validate_edu_model_gate_manifest(
    analysis, write_manifest(noncanonical)), "unique direct event identities")
})

test_that("model manifest rejects thresholds before accepting artifacts", {
  directory <- withr::local_tempdir(pattern = "SYNTHETIC_model_contract_")
  analysis <- structure(list(
    config = list(plot_type = "edu", gating = list(
      mode = "model_experimental", profile = "single_g1_edu_frozen_v1",
      output_dir = directory)),
    sample_manifest = data.frame(prefix = "SYNTHETIC_A"),
    input_report = data.frame()
  ), class = "facs_analysis")
  manifest <- list(
    status_label = "EXPERIMENTAL MODEL-DERIVED NON-PRODUCTION",
    gating_mode = "model_experimental", profile = "single_g1_edu_frozen_v1",
    g1_dna_minimum = list(
      method = "median_positive_dna_area_of_preliminary_edu_negative_g1",
      fraction = 0.35),
    frozen_thresholds = list(single_cells = 0.205, g1 = 0.255,
                             edu_positive = 0.5))
  path <- file.path(directory, "model-gating-manifest.json")
  jsonlite::write_json(manifest, path, auto_unbox = TRUE)
  expect_error(validate_edu_model_gate_manifest(analysis, path),
               "frozen thresholds differ")
})

test_that("model CSV hashing and parsing use the same bytes", {
  path <- withr::local_tempfile(pattern = "SYNTHETIC_model_csv_", fileext = ".csv")
  original <- charToRaw("event_identity,event_index\noriginal,1\n")
  writeBin(original, path)
  testthat::local_mocked_bindings(
    facs_read_file_bytes = function(file_path, size) {
      writeLines(c("event_identity,event_index", "replacement,2"), file_path)
      original
    },
    .package = "facspseudocolor"
  )
  observed <- facspseudocolor:::facs_read_hashed_csv(path)
  expect_identical(observed$sha256,
                   paste0(openssl::sha256(original), collapse = ""))
  expect_identical(observed$data$event_identity, "original")
  expect_identical(observed$data$event_index, 1L)
})
