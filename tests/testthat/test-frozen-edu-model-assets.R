test_that("frozen EdU profile bundles its approved model collection", {
  root <- test_path("..", "..", "inst", "models")
  collection <- file.path(root, "frozen_edu_models_v1")
  manifest <- jsonlite::fromJSON(
    file.path(collection, "MANIFEST.json"), simplifyVector = FALSE
  )
  expect_identical(manifest$schema, "FROZEN_EDU_MODEL_COLLECTION_V1")
  expect_identical(manifest$profile, "single_g1_edu_frozen_v1")
  expect_identical(
    manifest$status_label, "EXPERIMENTAL MODEL-DERIVED NON-PRODUCTION"
  )
  expect_identical(manifest$protected_data_present, FALSE)
  expect_setequal(
    list.files(collection),
    c("MANIFEST.json", "single_cells_frozen_development_rule.json",
      "g1_frozen_development_rule.json")
  )
  expect_identical(
    names(manifest$models), c("single_cells", "g1", "edu_positive")
  )
  expected <- c(
    single_cells = "71b459d8fb00930f31a6d289a21f587226fd2d6be4b31ebc782b7bae7839a376",
    g1 = "d653d388d15af3bd22185cb9c8e1addea132e0b5d1d18e60e37c65a7548163b7",
    edu_positive = "9e6ed466687efe5f7523dea633e9e9244381c0e613c86f0944f861c06047fdf4"
  )
  for (target in names(expected)) {
    record <- manifest$models[[target]]
    path <- normalizePath(file.path(collection, record$path), mustWork = TRUE)
    expect_true(startsWith(
      path,
      paste0(normalizePath(root, mustWork = TRUE), .Platform$file.sep)
    ))
    observed <- paste0(as.character(openssl::sha256(file(path))), collapse = "")
    expect_identical(record$byte_sha256, unname(expected[[target]]))
    expect_identical(observed, unname(expected[[target]]))
  }
  expect_identical(
    manifest$models$single_cells$semantic_sha256,
    "f9a53e9fd78d2f39cf5980b8f470c3c44e7a187fa942cbcc0324c563fa68a31a"
  )
  expect_identical(
    manifest$models$g1$semantic_sha256,
    "32723ccc768ed517affa040f2ce62a125dad465399d869f471bc331549f47f6a"
  )
  expect_identical(
    manifest$g1_definition$exclude_edu_positive_before_background, TRUE
  )
  expect_identical(
    manifest$g1_dna_minimum$method,
    "median_positive_dna_area_of_preliminary_edu_negative_g1"
  )
  expect_equal(
    manifest$g1_dna_minimum$fraction, 0.35
  )
})
