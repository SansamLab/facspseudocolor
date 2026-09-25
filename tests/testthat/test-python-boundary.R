test_that("installed package contains inert Python runtime but no launcher", {
  namespace <- asNamespace("facspseudocolor")
  expect_false(exists("prepare_flowjo_csvs", envir = namespace, inherits = FALSE))
  expect_false(exists("prepare_flowjo_csvs_external", envir = namespace,
                      inherits = FALSE))

  package_root <- getNamespaceInfo(namespace, "path")
  if (file.exists(file.path(package_root, ".Rbuildignore"))) {
    # pkgload resolves system.file() against the source checkout. In that
    # context, verify the rule that keeps repository-only launchers out of the
    # built package; the installed-package branch below verifies the runtime.
    build_ignore <- readLines(file.path(package_root, ".Rbuildignore"),
                              warn = FALSE)
    expect_true("^python$" %in% build_ignore)
    expect_true(file.exists(file.path(package_root, "inst", "python",
                                      "export_flowjo_populations.py")))
  } else {
    runtime <- system.file("python", package = "facspseudocolor")
    expect_true(nzchar(runtime))
    expect_true(file.exists(file.path(runtime, "export_flowjo_populations.py")))
    expect_true(file.exists(file.path(runtime, "export_contract.py")))
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

flowjo_orchestration_fixture <- function(plot_type, extra_flowjo = list()) {
  root <- tempfile("flowjo-orchestration-")
  dir.create(root)
  source_dir <- file.path(root, "source")
  dir.create(source_dir)
  writeLines("<wsp/>", file.path(source_dir, "workspace.wsp"))
  data_dir <- file.path(root, "data")
  dir.create(data_dir)
  exporter <- file.path(root, "export_flowjo_populations.py")
  writeLines("# stub", exporter)
  writeLines("# stub", file.path(root, "export_contract.py"))
  python_stub <- file.path(root, "python3")
  writeLines("# stub", python_stub)

  flowjo <- utils::modifyList(
    list(source_dir = source_dir, workspace = "workspace.wsp",
         dna_source_channel = "FL2-A", target_source_channel = "FL4-A",
         python = python_stub),
    extra_flowjo
  )
  config <- structure(list(
    plot_type = plot_type, data_dir = data_dir,
    suffixes = list(complete = "_single_cells.csv"),
    flowjo = flowjo,
    replicates = list(list(
      label = "Replicate 1", flowjo = list(),
      samples = list(list(label = "Sample 1", prefix = "sample1"))
    ))
  ), class = c("facs_config", "list"))
  list(config = config, root = root, exporter = exporter, data_dir = data_dir)
}

test_that("non-pH3 FlowJo orchestration uses the legacy profile without the production contract", {
  orchestration_path <- test_path("..", "..", "tools", "flowjo-orchestration.R")
  skip_if_not(file.exists(orchestration_path))
  orchestration <- new.env(parent = globalenv())
  sys.source(orchestration_path, envir = orchestration)

  fixture <- flowjo_orchestration_fixture("edu")
  export_dir <- file.path(fixture$data_dir, ".flowjo_export", "Replicate_1")
  captured_args <- NULL
  orchestration$system2 <- function(command, args, ...) {
    if ("--verify-operation" %in% args) return(0L)
    if (identical(args[[1]], "-c")) return(0L)
    captured_args <<- args
    dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)
    writeLines("gate_name,exported_count", file.path(export_dir, "population_counts.csv"))
    identity_columns <- c(
      "acquisition_id", "event_index", "event_identity", "identity_source",
      "identity_method_id", "identity_method_version", "duplicate_occurrence",
      "export_profile", "export_operation_id", "export_manifest_digest",
      "export_manifest_reference"
    )
    row <- setNames(as.list(rep("", length(identity_columns))), identity_columns)
    row$export_profile <- "legacy_count_only_unverified_v1"
    row$export_operation_id <- "LEGACY-UNVERIFIED"
    utils::write.csv(
      as.data.frame(row, stringsAsFactors = FALSE),
      file.path(export_dir, "sample1_single_cells.csv"), row.names = FALSE
    )
    0L
  }

  expect_silent(orchestration$prepare_flowjo_csvs_external(
    fixture$config, exporter = fixture$exporter, verbose = FALSE
  ))
  expect_true("legacy_count_only_unverified_v1" %in% captured_args)
  expect_false("production_direct_identity_v1" %in% captured_args)
  expect_false(any(grepl("--contract-metadata", captured_args, fixed = TRUE)))
  expect_false(any(grepl("--export-operation-id", captured_args, fixed = TRUE)))
  expect_false("--direct-index-semantics-verified" %in% captured_args)
})

test_that("pH3 FlowJo orchestration still requires and uses the production contract", {
  orchestration_path <- test_path("..", "..", "tools", "flowjo-orchestration.R")
  skip_if_not(file.exists(orchestration_path))
  orchestration <- new.env(parent = globalenv())
  sys.source(orchestration_path, envir = orchestration)

  fixture <- flowjo_orchestration_fixture("ph3")
  expect_error(
    orchestration$prepare_flowjo_csvs_external(
      fixture$config, exporter = fixture$exporter, verbose = FALSE
    ),
    "Missing FlowJo setting"
  )

  fixture <- flowjo_orchestration_fixture("ph3", extra_flowjo = list(
    contract_metadata = "contract.json", export_operation_id = "SYNTHETIC-OP",
    direct_index_semantics_verified = TRUE
  ))
  writeLines("{}", file.path(fixture$config$flowjo$source_dir, "contract.json"))
  export_dir <- file.path(fixture$data_dir, ".flowjo_export", "Replicate_1")
  captured_args <- NULL
  orchestration$system2 <- function(command, args, ...) {
    if ("--verify-operation" %in% args) return(0L)
    if (identical(args[[1]], "-c")) return(0L)
    captured_args <<- args
    dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)
    writeLines("gate_name,exported_count", file.path(export_dir, "population_counts.csv"))
    identity_columns <- c(
      "acquisition_id", "event_index", "event_identity", "identity_source",
      "identity_method_id", "identity_method_version", "duplicate_occurrence",
      "export_profile", "export_operation_id", "export_manifest_digest",
      "export_manifest_reference"
    )
    row <- setNames(as.list(rep("", length(identity_columns))), identity_columns)
    row$export_profile <- "production_direct_identity_v1"
    row$export_operation_id <- "SYNTHETIC-OP"
    utils::write.csv(
      as.data.frame(row, stringsAsFactors = FALSE),
      file.path(export_dir, "sample1_single_cells.csv"), row.names = FALSE
    )
    0L
  }
  expect_silent(orchestration$prepare_flowjo_csvs_external(
    fixture$config, exporter = fixture$exporter, verbose = FALSE
  ))
  expect_true("production_direct_identity_v1" %in% captured_args)
  expect_true(any(grepl("--contract-metadata", captured_args, fixed = TRUE)))
  expect_true("--direct-index-semantics-verified" %in% captured_args)
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
