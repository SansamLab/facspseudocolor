test_that("configuration validation applies centralized mode defaults", {
  edu <- validate_facs_config(minimal_config("edu"))
  poi <- validate_facs_config(minimal_config("poi"))

  expect_s3_class(edu, "facs_config")
  expect_identical(edu$suffixes$edu_positive, "_edu_positive.csv")
  expect_identical(poi$suffixes, list(complete = "_single_cells.csv"))
  expect_equal(edu$dna_2n_value, 1000)
  expect_identical(edu$layout, "plotgardener")
})

test_that("unknown and missing settings are reported together", {
  config <- minimal_config()
  config$dna_channel <- NULL
  config$misspelled_setting <- TRUE

  error <- tryCatch(validate_facs_config(config), error = identity)
  expect_s3_class(error, "error")
  expect_match(conditionMessage(error), "misspelled_setting")
  expect_match(conditionMessage(error), "dna_channel")
})

test_that("configuration validates references, prefixes, and ranges", {
  config <- minimal_config()
  config$replicates[[1]]$reference <- "Not present"
  config$x_limits <- c(2000, 1000)
  expect_error(validate_facs_config(config), "exactly one matching.*x_limits")

  config <- minimal_config()
  config$replicates[[1]]$samples[[2]]$prefix <- "reference"
  expect_error(validate_facs_config(config), "Duplicate sample prefix")
})

test_that("one-replicate configurations must still name a reference", {
  config <- minimal_config()
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
