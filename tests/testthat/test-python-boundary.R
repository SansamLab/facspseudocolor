test_that("installed package contains no Python launcher", {
  namespace <- asNamespace("facspseudocolor")
  expect_false(exists("prepare_flowjo_csvs", envir = namespace, inherits = FALSE))
  expect_false(exists("prepare_flowjo_csvs_external", envir = namespace,
                      inherits = FALSE))
  expect_identical(system.file("python", package = "facspseudocolor"), "")
})
