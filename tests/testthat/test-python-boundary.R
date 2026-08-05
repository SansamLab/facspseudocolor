test_that("installed package contains no Python launcher", {
  namespace <- asNamespace("facspseudocolor")
  expect_false(exists("prepare_flowjo_csvs", envir = namespace, inherits = FALSE))
  expect_false(exists("prepare_flowjo_csvs_external", envir = namespace,
                      inherits = FALSE))

  package_root <- getNamespaceInfo(namespace, "path")
  if (file.exists(file.path(package_root, ".Rbuildignore"))) {
    # pkgload resolves system.file() against the source checkout. In that
    # context, verify the rule that keeps the repository-only tools out of the
    # built package; the installed-package branch below verifies the result.
    build_ignore <- readLines(file.path(package_root, ".Rbuildignore"),
                              warn = FALSE)
    expect_true("^python$" %in% build_ignore)
    expect_false(dir.exists(file.path(package_root, "inst", "python")))
  } else {
    expect_identical(system.file("python", package = "facspseudocolor"), "")
  }
})
