test_that("public manifest builder preserves configured order", {
  config <- validate_facs_config(minimal_config("poi"))
  manifest <- build_sample_manifest(config)
  expect_identical(manifest$condition, c("Reference", "Treatment"))
  expect_identical(manifest$prefix, c("reference", "treatment"))
  expect_identical(which(manifest$is_reference), 1L)
})

test_that("sample reader accepts CSVs and data frames without changing rows", {
  events <- data.frame(DNA = c(1, 2, NA), Target = c(10, 20, 30), Extra = 4:6)
  from_memory <- read_facs_sample(events, "DNA", "Target", "memory sample")
  expect_equal(from_memory, events, ignore_attr = TRUE)
  expect_equal(attr(from_memory, "facs_input")$nonfinite_event_n, 1)

  path <- withr::local_tempfile(fileext = ".csv")
  utils::write.csv(events, path, row.names = FALSE)
  from_csv <- read_facs_sample(path, "DNA", "Target", "CSV sample")
  expect_equal(from_csv, events, ignore_attr = TRUE)
  expect_identical(attr(from_csv, "facs_input")$input_type, "csv")
})

test_that("sample reader rejects missing and nonnumeric channels", {
  expect_error(read_facs_sample(data.frame(DNA = 1:3), "DNA", "Target"),
               "Missing required channel")
  expect_error(read_facs_sample(
    data.frame(DNA = 1:3, Target = letters[1:3]), "DNA", "Target"),
    "must be numeric"
  )
  expect_error(read_facs_sample(
    data.frame(DNA = c(NA_real_, NA_real_), Target = c(1, 2)), "DNA", "Target"),
    "Too few finite events"
  )
})

test_that("input validation reports every missing required file", {
  config <- validate_facs_config(minimal_config("poi"))
  directory <- withr::local_tempdir()
  error <- tryCatch(validate_facs_inputs(config, directory), error = identity)
  expect_s3_class(error, "error")
  expect_match(conditionMessage(error), "reference_single_cells.csv")
  expect_match(conditionMessage(error), "treatment_single_cells.csv")
})

test_that("input validation checks included example files", {
  edu <- read_facs_config(system.file(
    "config", "config_edu.yml", package = "facspseudocolor"
  ))
  poi <- read_facs_config(system.file(
    "config", "config_poi.yml", package = "facspseudocolor"
  ))
  edu_files <- validate_facs_inputs(edu, example_data_dir("example"))
  poi_files <- validate_facs_inputs(poi, example_data_dir("example_poi"))

  expect_true(all(edu_files$exists))
  expect_true(all(poi_files$exists))
  expect_true(all(edu_files$finite_event_n > 0))
  expect_true(all(poi_files$finite_event_n > 0))
})
