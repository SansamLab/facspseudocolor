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
    expect_match(text, "prepare_edu_g1_inputs_external", fixed = TRUE)
    expect_match(text, "Documented model profile requires its existing exported manifest",
                 fixed = TRUE)
  } else {
    expect_false(file.exists(orchestration_path))
    expect_identical(
      system.file("tools", "flowjo-orchestration.R", package = "facspseudocolor"),
      ""
    )
  }
})

test_that("repository orchestration dispatches validated FlowJo configs canonically", {
  orchestration_path <- test_path("..", "..", "tools", "flowjo-orchestration.R")
  skip_if_not(file.exists(orchestration_path))
  orchestration <- new.env(parent = globalenv())
  sys.source(orchestration_path, envir = orchestration)
  called <- FALSE
  orchestration$prepare_flowjo_csvs_external <- function(config, verbose = TRUE) {
    called <<- TRUE
    invisible("flowjo")
  }
  config <- structure(list(
    plot_type = "edu", gating = list(mode = "flowjo")
  ), class = c("facs_config", "list"))
  expect_identical(
    orchestration$prepare_edu_g1_inputs_external(config, verbose = FALSE),
    "flowjo"
  )
  expect_true(called)
})

test_that("repository orchestration reuses a documented-profile manifest", {
  orchestration_path <- test_path("..", "..", "tools", "flowjo-orchestration.R")
  skip_if_not(file.exists(orchestration_path))
  orchestration <- new.env(parent = globalenv())
  sys.source(orchestration_path, envir = orchestration)
  manifest <- tempfile(fileext = ".json")
  writeLines("{}", manifest)
  config <- structure(list(
    plot_type = "edu",
    gating = list(mode = "model_experimental", profile = "documented_edu_g1_v1",
                  manifest = manifest)
  ), class = c("facs_config", "list"))
  expect_identical(
    orchestration$prepare_edu_g1_inputs_external(config, verbose = FALSE),
    normalizePath(manifest)
  )
})

test_that("repository orchestration invokes the frozen profile caller", {
  orchestration_path <- test_path("..", "..", "tools", "flowjo-orchestration.R")
  skip_if_not(file.exists(orchestration_path))
  orchestration <- new.env(parent = globalenv())
  sys.source(orchestration_path, envir = orchestration)
  output_dir <- tempfile("frozen-output-")
  dir.create(output_dir)
  config_path <- tempfile(fileext = ".yml")
  caller <- tempfile(fileext = ".py")
  writeLines("plot_type: edu", config_path)
  writeLines("# test caller", caller)
  orchestration$system2 <- function(command, args) {
    writeLines("{}", file.path(output_dir, "model-gating-manifest.json"))
    0L
  }
  config <- structure(list(
    plot_type = "edu", flowjo = list(python = Sys.which("python3")),
    gating = list(mode = "model_experimental", profile = "single_g1_edu_frozen_v1",
                  output_dir = output_dir)
  ), class = c("facs_config", "list"))
  attr(config, "config_path") <- config_path
  expect_identical(
    orchestration$prepare_edu_g1_inputs_external(
      config, frozen_model_caller = caller, verbose = FALSE
    ),
    normalizePath(file.path(output_dir, "model-gating-manifest.json"))
  )
})
