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

test_that("repository orchestration prohibits sequential identity fallback", {
  namespace <- asNamespace("facspseudocolor")
  package_root <- getNamespaceInfo(namespace, "path")
  source_checkout <- file.exists(file.path(package_root, ".Rbuildignore"))
  orchestration_path <- test_path(
    "..", "..", "tools", "flowjo-orchestration.R"
  )
  if (source_checkout) {
    expect_true(file.exists(orchestration_path))
    orchestration <- readLines(orchestration_path, warn = FALSE)
    text <- paste(orchestration, collapse = "\n")
    expect_match(text, "Sequential identity fallback is prohibited", fixed = TRUE)
    expect_false(grepl("seq_len(sum(selected)) - 1L", text, fixed = TRUE))
    expect_match(text, '"event_identity"', fixed = TRUE)
    expect_match(text, 'colClasses = stats::setNames', fixed = TRUE)
    expect_match(text, "Legacy or ambiguous FlowJo exports cannot be consumed",
                 fixed = TRUE)
    expect_match(text, '"export_manifest_reference"', fixed = TRUE)
    expect_match(text, '"export_manifest_digest"', fixed = TRUE)
    expect_false(grepl("utils::write.csv(output, output_file", text, fixed = TRUE))
    expect_match(text, '"--verify-operation"', fixed = TRUE)
    expect_match(text, "Finalized manifest or consumed population artifact verification failed",
                 fixed = TRUE)
  } else {
    expect_false(file.exists(orchestration_path))
    expect_identical(
      system.file("tools", "flowjo-orchestration.R", package = "facspseudocolor"),
      ""
    )
  }
})
