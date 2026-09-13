test_that("configuration validation applies centralized mode defaults", {
  edu <- validate_facs_config(minimal_config("edu"))
  poi <- validate_facs_config(minimal_config("poi"))

  expect_s3_class(edu, "facs_config")
  expect_identical(edu$suffixes$edu_positive, "_edu_positive.csv")
  expect_identical(poi$suffixes, list(complete = "_single_cells.csv"))
  expect_equal(edu$dna_2n_value, 1000)
  expect_identical(edu$layout, "cowplot")
  expect_identical(edu$pseudocolor_signal, "background_subtracted")
  expect_identical(edu$background_subtracted_offset, "auto")
  expect_identical(edu$gating, list(mode = "flowjo"))
})

test_that("experimental model gating is explicit and has no FlowJo fallback", {
  config <- minimal_config("edu")
  config$data_dir <- "SYNTHETIC-model-output"
  single_artifact <- list(
    path = "SYNTHETIC-single-artifact",
    byte_sha256 = "71b459d8fb00930f31a6d289a21f587226fd2d6be4b31ebc782b7bae7839a376",
    semantic_sha256 = "f9a53e9fd78d2f39cf5980b8f470c3c44e7a187fa942cbcc0324c563fa68a31a")
  g1_artifact <- list(
    path = "SYNTHETIC-g1-artifact",
    byte_sha256 = "d653d388d15af3bd22185cb9c8e1addea132e0b5d1d18e60e37c65a7548163b7",
    semantic_sha256 = "32723ccc768ed517affa040f2ce62a125dad465399d869f471bc331549f47f6a")
  config$gating <- list(
    mode = "model_experimental", profile = "single_g1_edu_frozen_v1",
    status_label = "EXPERIMENTAL MODEL-DERIVED NON-PRODUCTION",
    output_dir = "SYNTHETIC-model-output",
    artifacts = list(single_cells = single_artifact, g1 = g1_artifact,
                     edu_positive = list(path = "SYNTHETIC-model.portable.json",
                                         byte_sha256 = "9e6ed466687efe5f7523dea633e9e9244381c0e613c86f0944f861c06047fdf4")),
    acquisitions = lapply(config$replicates[[1]]$samples, function(sample) list(
      prefix = sample$prefix, acquisition_id = paste0("SYNTHETIC-", sample$prefix),
      fcs_path = paste0("/SYNTHETIC/", sample$prefix, ".fcs"),
      fcs_sha256 = paste0(rep("0", 64), collapse = ""),
      channel_roles = list(fsc_area = "FSC-A", ssc_area = "SSC-A",
                           dna_area = "FL2-A", dna_height = "FL2-H", edu = "FL4-A")
    ))
  )
  expect_identical(validate_facs_config(config)$gating$mode, "model_experimental")
  config$gating$profile <- "automatic"
  expect_error(validate_facs_config(config), "approved frozen profile")
  config$gating$profile <- "single_g1_edu_frozen_v1"
  config$gating$artifacts$g1$byte_sha256 <- paste0(rep("0", 64), collapse = "")
  expect_error(validate_facs_config(config), "differs from frozen profile for g1")
  config$gating$artifacts$g1$byte_sha256 <-
    "d653d388d15af3bd22185cb9c8e1addea132e0b5d1d18e60e37c65a7548163b7"
  config$data_dir <- "SYNTHETIC-flowjo-fallback"
  expect_error(validate_facs_config(config), "no fallback input")
})

test_that("pseudocolor signal and offset settings are validated", {
  config <- minimal_config("edu")
  config$pseudocolor_signal <- "background_subtracted"
  config$background_subtracted_offset <- 10000
  result <- validate_facs_config(config)
  expect_equal(result$background_subtracted_offset, 10000)

  config$background_subtracted_offset <- "nearest"
  expect_error(validate_facs_config(config), "background_subtracted_offset")
  config$background_subtracted_offset <- -1
  expect_error(validate_facs_config(config), "background_subtracted_offset")
})

test_that("EdU report settings preserve DNA-H and interactive output choices", {
  config <- minimal_config("edu")
  config$report <- list(
    all_events_dir = "all_events_csv",
    dna_height_channel = "FL2-H",
    show_apex_comparison_toggle = TRUE,
    show_phase_gate_toggle = TRUE,
    embed_pseudocolor_pdf_downloads = TRUE
  )
  result <- validate_facs_config(config)
  expect_identical(result$report, config$report)

  config$report$dna_height_channel <- ""
  expect_error(validate_facs_config(config), "dna_height_channel")
  config$report$dna_height_channel <- "FL2-H"
  config$report$unknown <- TRUE
  expect_error(validate_facs_config(config), "Unknown EdU `report` setting")

  poi <- minimal_config("poi")
  poi$report <- result$report
  expect_error(validate_facs_config(poi), "supported only for EdU")
})

test_that("unknown and missing settings are reported together", {
  config <- minimal_config("poi")
  config$dna_channel <- NULL
  config$misspelled_setting <- TRUE

  error <- tryCatch(validate_facs_config(config), error = identity)
  expect_s3_class(error, "error")
  expect_match(conditionMessage(error), "misspelled_setting")
  expect_match(conditionMessage(error), "dna_channel")
})

test_that("configuration validates references, prefixes, and ranges", {
  config <- minimal_config("poi")
  config$replicates[[1]]$reference <- "Not present"
  config$x_limits <- c(2000, 1000)
  expect_error(validate_facs_config(config), "exactly one matching.*x_limits")

  config <- minimal_config("poi")
  config$replicates[[1]]$samples[[2]]$prefix <- "reference"
  expect_error(validate_facs_config(config), "Duplicate sample prefix")
})

test_that("EdU permits an explicit matched reference in every replicate pair", {
  config <- minimal_config("edu")
  config$replicates[[1]]$reference <- "Reference"
  validated <- validate_facs_config(config)
  manifest <- build_sample_manifest(validated)

  expect_identical(manifest$prefix[manifest$is_reference], "reference")

  second <- config$replicates[[1]]
  second$label <- "Replicate 2"
  second$samples <- lapply(second$samples, function(sample) {
    sample$prefix <- paste0("second_", sample$prefix)
    sample
  })
  config$replicates[[2]] <- second
  expect_silent(validate_facs_config(config))

  config$replicates[[2]]$reference <- NULL
  expect_error(
    validate_facs_config(config),
    "EdU replicate `reference` samples must be absent.*exactly one"
  )

  config$replicates[[2]]$reference <- "Missing"
  expect_error(
    validate_facs_config(config),
    "EdU replicate `reference` samples must be absent.*exactly one"
  )

  single_invalid <- minimal_config("edu")
  single_invalid$replicates[[1]]$reference <- "Missing"
  expect_error(
    validate_facs_config(single_invalid),
    "each declared reference must name exactly one matching sample"
  )

  all_invalid <- config
  all_invalid$replicates[[1]]$reference <- "Missing"
  all_invalid$replicates[[2]]$reference <- "Also missing"
  expect_error(
    validate_facs_config(all_invalid),
    "each declared reference must name exactly one matching sample"
  )
})

test_that("flat configurations cannot express replicate reference structure", {
  config <- minimal_config("poi")
  config$samples <- config$replicates[[1]]$samples
  config$replicates <- NULL
  expect_error(validate_facs_config(config), "flat `samples` list.*reference")
})

test_that("configuration reader uses only the explicit supplied path", {
  path <- withr::local_tempfile(fileext = ".yml")
  yaml::write_yaml(minimal_config("poi"), path)
  config <- read_facs_config(path)

  expect_s3_class(config, "facs_config")
  expect_identical(attr(config, "config_path"), normalizePath(path))
  expect_identical(attr(config, "config_dir"), dirname(normalizePath(path)))
  expect_error(read_facs_config(paste0(path, "-missing")), "not found")
})

test_that("included EdU and POI configurations validate", {
  edu <- read_facs_config(system.file(
    "config", "config_edu.yml", package = "facspseudocolor"
  ))
  poi <- read_facs_config(system.file(
    "config", "config_poi.yml", package = "facspseudocolor"
  ))
  expect_identical(edu$plot_type, "edu")
  expect_identical(poi$plot_type, "poi")
})

test_that("plotgardener is optional and reports actionable installation help", {
  expect_invisible(require_plotgardener(TRUE))
  expect_error(
    require_plotgardener(FALSE),
    "BiocManager::install\\('plotgardener'\\)"
  )
})
