synthetic_ph3_population <- function(
    key, identities, acquisition = "SYNTHETIC-A",
    sample = "SYNTHETIC-A.fcs", gate = paste0("root/Single Cells/", key),
    parent = "root/Single Cells"
) {
  events <- data.frame(event_identity = identities, stringsAsFactors = FALSE)
  list(
    events = events,
    population = list(
      acquisition_id = acquisition, sample_id = sample, population_key = key,
      gate_path = gate, parent_population_path = parent
    ),
    artifact = list(sha256 = paste(rep("a", 64), collapse = "")),
    path = paste0("/SYNTHETIC/", key, ".csv"), duplicate_rows = 0L,
    duplicate_bases = 0L, intentionally_empty = !length(identities)
  )
}

synthetic_ph3_operation <- function() {
  list(
    full_digest = paste(rep("b", 64), collapse = ""),
    manifest_path = "/SYNTHETIC/export-manifest.json",
    manifest = list(
      profile = "production_direct_identity_v1",
      export_operation_id = "SYNTHETIC-OP",
      identity_method = list(
        id = "flowkit_source_event_index_scoped_v1", version = "1.0.0",
        source = "flowkit_get_gate_events_index"
      )
    )
  )
}

test_that("direct G1 and pH3 containment is independent and may overlap", {
  parent <- synthetic_ph3_population(
    "complete", paste0("SYNTHETIC-A:event_index:", 1:4),
    gate = "root/Single Cells", parent = "root"
  )
  g1 <- synthetic_ph3_population(
    "g1", paste0("SYNTHETIC-A:event_index:", c(1, 2))
  )
  positive <- synthetic_ph3_population(
    "ph3_positive", paste0("SYNTHETIC-A:event_index:", c(2, 4))
  )
  g1_result <- ph3_containment_row(synthetic_ph3_operation(), "SYNTHETIC_A",
                                   parent, g1)
  positive_result <- ph3_containment_row(
    synthetic_ph3_operation(), "SYNTHETIC_A", parent, positive
  )
  expect_identical(g1_result$containment_status, "validated")
  expect_identical(positive_result$containment_status, "validated")
  expect_identical(g1_result$child_row_count, 2L)
  expect_identical(positive_result$matched_child_count, 2L)
  expect_identical(g1_result$unmatched_child_count, 0L)
  expect_identical(g1_result$excess_occurrence_count, 0L)
})

test_that("cross-acquisition and unmatched identities fail with exact scope", {
  parent <- synthetic_ph3_population(
    "complete", "SYNTHETIC-A:event_index:1", gate = "root/Single Cells",
    parent = "root"
  )
  wrong <- synthetic_ph3_population(
    "g1", "SYNTHETIC-B:event_index:1", acquisition = "SYNTHETIC-B",
    sample = "SYNTHETIC-B.fcs"
  )
  expect_error(
    ph3_containment_row(synthetic_ph3_operation(), "SYNTHETIC_A", parent, wrong),
    "parent_child_acquisition_mismatch.*SYNTHETIC-B.*g1"
  )
  unmatched <- synthetic_ph3_population("ph3_positive",
                                        "SYNTHETIC-A:event_index:9")
  expect_error(
    ph3_containment_row(synthetic_ph3_operation(), "SYNTHETIC_A", parent,
                        unmatched),
    "child_not_contained_in_parent.*SYNTHETIC-A.*ph3_positive"
  )
})

test_that("child less than or equal to parent succeeds and excess fails", {
  parent <- synthetic_ph3_population(
    "complete", paste0("SYNTHETIC-A:event_index:", 1:2),
    gate = "root/Single Cells", parent = "root"
  )
  equal <- synthetic_ph3_population(
    "g1", paste0("SYNTHETIC-A:event_index:", 1:2)
  )
  expect_identical(
    ph3_containment_row(synthetic_ph3_operation(), "SYNTHETIC_A", parent,
                        equal)$child_row_count,
    2L
  )
  duplicate <- synthetic_ph3_population(
    "g1", c("SYNTHETIC-A:event_index:1", "SYNTHETIC-A:event_index:1")
  )
  duplicate$duplicate_rows <- 2L
  duplicate$duplicate_bases <- 1L
  expect_error(
    ph3_containment_row(synthetic_ph3_operation(), "SYNTHETIC_A", parent,
                        duplicate),
    "excess occurrences"
  )
})

test_that("structured containment preserves empty proof and immutable references", {
  parent <- synthetic_ph3_population(
    "complete", "SYNTHETIC-A:event_index:1", gate = "root/Single Cells",
    parent = "root"
  )
  empty <- synthetic_ph3_population("ph3_positive", character())
  result <- ph3_containment_row(synthetic_ph3_operation(), "SYNTHETIC_A",
                                parent, empty)
  expect_true(result$intentionally_empty)
  expect_identical(result$intentionally_empty_proof_status, "verified")
  expect_identical(result$containment_reason_code,
                   "validated_intentional_empty")
  expect_identical(result$containment_schema_version,
                   "ph3-input-containment-1.0.0")
  expect_match(result$manifest_reference, "SYNTHETIC")
})

test_that("PH3 input profile is explicit and legacy is isolated", {
  config <- minimal_config("ph3")
  config$g1_x_range <- c(750, 1250)
  config$s_phase_bins <- list(early = c(1250, 1450), mid = c(1450, 1650),
                              late = c(1650, 1800))
  config$g2m_x_range <- c(1800, 2400)
  config$ph3_input_profile <- NULL
  expect_error(validate_facs_config(config), "requires explicit")
  config$ph3_input_profile <- "production_direct_identity_v1"
  expect_error(validate_facs_config(config), "ph3_export_operation_dirs")
  legacy <- config
  legacy$ph3_input_profile <- "legacy_count_only_unverified_v1"
  expect_identical(validate_facs_config(legacy)$ph3_input_profile,
                   "legacy_count_only_unverified_v1")
})

test_that("legacy input warns once and yields one unverified row per child", {
  directory <- withr::local_tempdir(pattern = "SYNTHETIC_legacy_")
  config <- minimal_config("ph3")
  config$g1_x_range <- c(750, 1250)
  config$s_phase_bins <- list(early = c(1250, 1450), mid = c(1450, 1650),
                              late = c(1650, 1800))
  config$g2m_x_range <- c(1800, 2400)
  config <- validate_facs_config(config)
  for (prefix in c("reference", "treatment")) {
    complete <- data.frame(DNA = 1:10, Target = 101:110)
    g1 <- complete[1:2, , drop = FALSE]
    positive <- complete[1, , drop = FALSE]
    utils::write.csv(complete, file.path(directory, paste0(prefix,
      "_single_cells.csv")), row.names = FALSE)
    utils::write.csv(g1, file.path(directory, paste0(prefix, "_g1.csv")),
                     row.names = FALSE)
    utils::write.csv(positive, file.path(directory, paste0(prefix,
      "_ph3_positive.csv")), row.names = FALSE)
  }
  warnings <- character()
  report <- withCallingHandlers(
    validate_facs_inputs(config, directory),
    warning = function(warning) {
      warnings <<- c(warnings, conditionMessage(warning))
      invokeRestart("muffleWarning")
    }
  )
  containment <- attr(report, "ph3_containment")
  expect_length(warnings, 1L)
  expect_match(warnings[[1L]], "LEGACY COUNT-ONLY PH3 INPUT")
  expect_equal(nrow(containment), 4L)
  expect_setequal(containment$child_population_key, c("g1", "ph3_positive"))
  expect_true(all(containment$containment_status == "unverified"))
  expect_true(all(is.na(containment$identity_method_id)))
  expect_false(any(containment$containment_status == "validated"))
})

test_that("audited production namespace functions contain no unsafe APIs", {
  namespace <- asNamespace("facspseudocolor")
  binding_names <- ls(namespace, all.names = TRUE)
  bindings <- mget(binding_names, envir = namespace, inherits = FALSE)
  local_functions <- bindings[vapply(
    bindings,
    function(value) is.function(value) && !is.primitive(value),
    logical(1)
  )]
  sources <- vapply(local_functions, function(fun) {
    paste(deparse(body(fun), width.cutoff = 500L), collapse = "\n")
  }, character(1))
  all_text <- paste(sources, collapse = "\n")
  for (name in names(sources)) {
    expect_false(grepl(
      "(?s)event_identity\\s*(?:<-|=)\\s*.{0,300}(?:seq_len|seq_along|as[.]numeric)\\s*\\(",
      sources[[name]], perl = TRUE
    ), info = name)
    expect_false(grepl(
      "as[.]numeric\\s*\\(\\s*[^)]*event_identity",
      sources[[name]], perl = TRUE
    ), info = name)
  }
  expect_match(all_text, "duplicate_direct_identity", fixed = TRUE)
  expect_match(all_text, "manifest_digest_mismatch", fixed = TRUE)
  audited <- sources[grepl(
    paste0(
      "^(ph3_|normalize_ph3$|assign_ph3_dna_phase$|quantify_ph3$|",
      "plot_ph3_|analysis_settings$|ph3_normalization_inputs$|",
      "new_facs_analysis$|analyze_facs_experiment$|print[.]facs_analysis$|",
      "facs_config_|config_|validate_facs_config$|read_facs_config$|",
      "build_sample_manifest$|resolve_facs_directory$|",
      "required_population_keys$|facs_input_files$|read_facs_sample$|",
      "validate_facs_inputs$|validate_ph3_export_operations$)"
    ),
    names(sources)
  )]
  expect_gt(length(audited), 0L)
  forbidden <- paste(
    c(
      "https?://", "download[.]file\\s*\\(", "url\\s*\\(",
      "browseURL\\s*\\(", "socketConnection\\s*\\(", "curl::", "RCurl::",
      "httr::", "httr2::", "upload", "telemetry", "processx::", "callr::",
      "system\\s*\\(", "system2\\s*\\(", "shell\\s*\\(", "pipe\\s*\\(",
      "unlink\\s*\\(", "file[.]remove\\s*\\(", "file[.]rename\\s*\\(",
      "file[.]copy\\s*\\(", "Sys[.]setenv\\s*\\(",
      "Sys[.]unsetenv\\s*\\(", "dyn[.]load\\s*\\(", "source\\s*\\(",
      "eval\\s*\\(\\s*parse"
    ),
    collapse = "|"
  )
  for (name in names(audited)) {
    expect_false(grepl(forbidden, audited[[name]],
                       perl = TRUE, ignore.case = TRUE), info = name)
  }
})

test_that("R canonical JSON matches the fixed Slice 1 canonical form", {
  expect_identical(ph3_canonical_json(list(b = 2L, a = 1L)),
                   "{\"a\":1,\"b\":2}")
  expect_identical(
    ph3_canonical_json(list(z = list(TRUE, FALSE), a = "SYNTHETIC")),
    "{\"a\":\"SYNTHETIC\",\"z\":[true,false]}"
  )
})

write_synthetic_ph3_operation <- function(
    root, empty_positive = FALSE,
    prefixes = c("reference", "treatment"), operation_id = "SYNTHETIC-OP"
) {
  binding <- paste(rep("0", 64), collapse = "")
  identity_columns <- c(
    "acquisition_id", "event_index", "event_identity", "identity_source",
    "identity_method_id", "identity_method_version", "duplicate_occurrence",
    "export_profile", "export_operation_id", "export_manifest_digest",
    "export_manifest_reference"
  )
  acquisitions <- lapply(prefixes, function(prefix) list(
    acquisition_id = paste0("SYNTHETIC-", prefix),
    sample_id = paste0("SYNTHETIC-", prefix, ".fcs"), prefix = prefix,
    source_fcs_reference = paste0("SYNTHETIC-", prefix, ".fcs"),
    source_fcs_sha256 = paste(rep("1", 64), collapse = "")
  ))
  populations <- list()
  artifacts <- list()
  counts <- list()
  for (acquisition in acquisitions) {
    parent_ids <- paste0(acquisition$acquisition_id, ":event_index:", 1:10)
    identities <- list(complete = parent_ids, g1 = parent_ids[1:2],
                       ph3_positive = if (empty_positive) character() else parent_ids[2])
    gates <- list(complete = "root/Single Cells",
                  g1 = "root/Single Cells/G1",
                  ph3_positive = "root/Single Cells/pH3 Positive")
    for (key in names(identities)) {
      ids <- identities[[key]]
      index <- sub(".*:event_index:", "", ids)
      events <- data.frame(
        acquisition_id = rep(acquisition$acquisition_id, length(ids)),
        event_index = index, event_identity = ids,
        identity_source = rep("flowkit_get_gate_events_index", length(ids)),
        identity_method_id = rep("flowkit_source_event_index_scoped_v1", length(ids)),
        identity_method_version = rep("1.0.0", length(ids)),
        duplicate_occurrence = rep("1", length(ids)),
        export_profile = rep("production_direct_identity_v1", length(ids)),
        export_operation_id = rep(operation_id, length(ids)),
        export_manifest_digest = rep(binding, length(ids)),
        export_manifest_reference = rep(
          "export-manifest.json + export-manifest.sha256", length(ids)
        ),
        sample_id = rep(acquisition$sample_id, length(ids)),
        raw__DNA = seq_along(ids) + 100, raw__Target = seq_along(ids) + 200,
        scaled__DNA = seq_along(ids) + 100,
        scaled__Target = seq_along(ids) + 200,
        check.names = FALSE, stringsAsFactors = FALSE
      )
      path <- paste0(acquisition$prefix, "_", key, ".csv")
      utils::write.csv(events, file.path(root, path), row.names = FALSE,
                       quote = TRUE)
      digest <- ph3_sha256_file(file.path(root, path))
      gate_name <- c(complete = "Single Cells", g1 = "G1",
                     ph3_positive = "pH3 Positive")[[key]]
      populations[[length(populations) + 1L]] <- list(
        population_key = key, gate_name = gate_name, gate_path = gates[[key]],
        gate_type = "SYNTHETIC_GATE",
        parent_population_path = if (key == "complete") "root" else gates$complete,
        acquisition_id = acquisition$acquisition_id,
        sample_id = acquisition$sample_id, prefix = acquisition$prefix,
        channels = as.list(c("DNA", "Target")),
        gate_channels = as.list(c("DNA", "Target")), row_count = length(ids),
        identity_field = "event_identity",
        identity_method_id = "flowkit_source_event_index_scoped_v1",
        unique_identity_count = length(ids),
        duplicate_base_combination_count = 0L, duplicate_row_count = 0L,
        intentionally_empty = !length(ids), export_operation_id = operation_id,
        artifact_path = path, artifact_sha256 = digest
      )
      artifacts[[length(artifacts) + 1L]] <- list(
        role = "population_events", path = path, sha256 = digest,
        byte_size = unname(file.info(file.path(root, path))$size),
        row_count = length(ids), columns = as.list(names(events)),
        identity_columns = as.list(c("acquisition_id", "event_index",
                                     "event_identity")),
        identity_method_id = "flowkit_source_event_index_scoped_v1",
        intentionally_empty = !length(ids), export_operation_id = operation_id,
        acquisition_id = acquisition$acquisition_id,
        sample_id = acquisition$sample_id, population_key = key,
        gate_path = gates[[key]], channels = as.list(c("DNA", "Target"))
      )
      counts[[length(counts) + 1L]] <- data.frame(
        sample_id = acquisition$sample_id, gate_name = gate_name,
        exported_count = length(ids), flowjo_saved_count = length(ids),
        count_difference = 0L, stringsAsFactors = FALSE
      )
    }
  }
  count_report <- do.call(rbind, counts)
  count_path <- file.path(root, "population_counts.csv")
  utils::write.csv(count_report, count_path, row.names = FALSE, quote = TRUE)
  artifacts[[length(artifacts) + 1L]] <- list(
    role = "population_count_report", path = "population_counts.csv",
    sha256 = ph3_sha256_file(count_path),
    byte_size = unname(file.info(count_path)$size), row_count = nrow(count_report),
    columns = as.list(names(count_report)), identity_columns = list(),
    identity_method_id = NULL, intentionally_empty = FALSE,
    export_operation_id = operation_id, acquisition_id = NULL, sample_id = NULL,
    population_key = NULL, gate_path = NULL, channels = list()
  )
  manifest <- list(
    manifest_schema = list(name = "facspseudocolor-flowjo-export", version = "1.0.0"),
    export_operation_id = operation_id, created_at = "SYNTHETIC-TIME",
    status = "complete", profile = "production_direct_identity_v1",
    legacy_warning = NULL,
    export_profile = list(id = "ph3_export_identity_provenance_v1", version = "1.0.0"),
    identity_method = list(id = "flowkit_source_event_index_scoped_v1",
      version = "1.0.0", source = "flowkit_get_gate_events_index",
      semantics = "conditional_on_pinned_environment_verification",
      verified = TRUE),
    workspace = list(filename = "SYNTHETIC.wsp", local_path = "/SYNTHETIC.wsp",
                     sha256 = paste(rep("2", 64), collapse = "")),
    software = list(exporter = "SYNTHETIC", exporter_version = "1.0.0",
      source_commit = "SYNTHETIC", flowjo_version = "SYNTHETIC",
      flowkit_version = "1.3.2", supported_flowkit_version = "1.3.2",
      python_version = "SYNTHETIC"),
    transformation_context = "SYNTHETIC", compensation_context = "SYNTHETIC",
    approval = list(gate_owner = "SYNTHETIC OWNER", approver = "SYNTHETIC APPROVER",
      approval_date = "2026-08-24", approval_record = "SYNTHETIC/approval.md",
      positivity_method_id = "flowjo_owner_approved_positive_population_v1",
      positivity_method_version = "1"),
    acquisitions = acquisitions,
    requested_populations = as.list(c("complete", "g1", "ph3_positive")),
    populations = populations, artifacts = artifacts,
    geometry_overlay_status = "not_requested",
    manifest_binding = list(algorithm = "sha256",
      canonical_object = "manifest_binding_object_v1",
      excludes = as.list(c("created_at", "populations", "artifacts",
                           "full_manifest_digest")), digest = binding),
    manifest_digest = list(algorithm = "sha256",
                           target = "canonical UTF-8 export-manifest.json bytes",
                           storage = "export-manifest.sha256 sidecar")
  )
  binding <- ph3_manifest_binding_digest(manifest)
  manifest$manifest_binding$digest <- binding
  for (i in seq_along(populations)) {
    path <- file.path(root, populations[[i]]$artifact_path)
    events <- utils::read.csv(path, check.names = FALSE, colClasses = "character")
    events$export_manifest_digest <- rep(binding, nrow(events))
    utils::write.csv(events, path, row.names = FALSE, quote = TRUE)
    digest <- ph3_sha256_file(path)
    manifest$populations[[i]]$artifact_sha256 <- digest
    manifest$artifacts[[i]]$sha256 <- digest
    manifest$artifacts[[i]]$byte_size <- unname(file.info(path)$size)
  }
  manifest_text <- paste0(ph3_canonical_json(manifest), "\n")
  writeBin(charToRaw(manifest_text), file.path(root, "export-manifest.json"))
  full_digest <- ph3_sha256_file(file.path(root, "export-manifest.json"))
  writeLines(paste0(full_digest, "  export-manifest.json"),
             file.path(root, "export-manifest.sha256"), useBytes = TRUE)
  invisible(manifest)
}

synthetic_production_config <- function(operation_dir) {
  operation_dir <- vapply(operation_dir, normalizePath, character(1),
                          winslash = "/", mustWork = FALSE, USE.NAMES = FALSE)
  config <- minimal_config("ph3")
  config$data_dir <- dirname(operation_dir[[1L]])
  config$ph3_input_profile <- "production_direct_identity_v1"
  config$ph3_export_operation_dirs <- operation_dir
  config$dna_channel <- "raw__DNA"
  config$target_channel <- "raw__Target"
  config$g1_x_range <- c(750, 1250)
  config$s_phase_bins <- list(early = c(1250, 1450), mid = c(1450, 1650),
                              late = c(1650, 1800))
  config$g2m_x_range <- c(1800, 2400)
  validate_facs_config(config)
}

rewrite_synthetic_manifest <- function(root, manifest) {
  manifest_text <- paste0(ph3_canonical_json(manifest), "\n")
  writeBin(charToRaw(manifest_text), file.path(root, "export-manifest.json"))
  digest <- ph3_sha256_file(file.path(root, "export-manifest.json"))
  writeLines(paste0(digest, "  export-manifest.json"),
             file.path(root, "export-manifest.sha256"), useBytes = TRUE)
}

rewrite_synthetic_population <- function(root, manifest, index, events) {
  path <- file.path(root, manifest$populations[[index]]$artifact_path)
  utils::write.csv(events, path, row.names = FALSE, quote = TRUE)
  digest <- ph3_sha256_file(path)
  manifest$populations[[index]]$artifact_sha256 <- digest
  manifest$artifacts[[index]]$sha256 <- digest
  manifest$artifacts[[index]]$byte_size <- unname(file.info(path)$size)
  rewrite_synthetic_manifest(root, manifest)
  manifest
}

rewrite_synthetic_count_report <- function(root, manifest, report) {
  path <- file.path(root, "population_counts.csv")
  utils::write.csv(report, path, row.names = FALSE, quote = TRUE)
  index <- which(vapply(manifest$artifacts, function(artifact) {
    identical(artifact$role, "population_count_report")
  }, logical(1)))[[1L]]
  manifest$artifacts[[index]]$sha256 <- ph3_sha256_file(path)
  manifest$artifacts[[index]]$byte_size <- unname(file.info(path)$size)
  manifest$artifacts[[index]]$row_count <- nrow(report)
  rewrite_synthetic_manifest(root, manifest)
  manifest
}

test_that("full production validator retains verified in-memory events", {
  operation <- withr::local_tempdir(pattern = "SYNTHETIC_operation_")
  write_synthetic_ph3_operation(operation)
  config <- synthetic_production_config(operation)
  report <- validate_facs_inputs(config)
  expect_equal(nrow(attr(report, "ph3_containment")), 4L)
  retained <- attr(report, "ph3_verified_events")$reference$complete
  writeLines("SYNTHETIC MUTATION", report$path[report$prefix == "reference" &
                                               report$population == "complete"])
  consumed <- ph3_normalization_inputs(
    report, "reference", "production_direct_identity_v1"
  )
  expect_identical(consumed$complete, retained)
})

test_that("production validator rejects altered bindings, artifacts, and counts", {
  operation <- withr::local_tempdir(pattern = "SYNTHETIC_operation_")
  manifest <- write_synthetic_ph3_operation(operation)
  config <- synthetic_production_config(operation)
  manifest$manifest_binding$digest <- paste(rep("f", 64), collapse = "")
  manifest_text <- paste0(ph3_canonical_json(manifest), "\n")
  writeBin(charToRaw(manifest_text), file.path(operation, "export-manifest.json"))
  full <- ph3_sha256_file(file.path(operation, "export-manifest.json"))
  writeLines(paste0(full, "  export-manifest.json"),
             file.path(operation, "export-manifest.sha256"), useBytes = TRUE)
  expect_error(validate_facs_inputs(config), "manifest_binding_digest_mismatch")

  operation2 <- withr::local_tempdir(pattern = "SYNTHETIC_operation_")
  write_synthetic_ph3_operation(operation2)
  config2 <- synthetic_production_config(operation2)
  writeLines("SYNTHETIC MUTATION", file.path(operation2, "reference_g1.csv"))
  expect_error(validate_facs_inputs(config2), "artifact_ledger_mismatch")

  operation3 <- withr::local_tempdir(pattern = "SYNTHETIC_operation_")
  manifest3 <- write_synthetic_ph3_operation(operation3)
  config3 <- synthetic_production_config(operation3)
  counts <- utils::read.csv(file.path(operation3, "population_counts.csv"))
  counts$exported_count[[1L]] <- counts$exported_count[[1L]] + 1L
  counts$flowjo_saved_count[[1L]] <- counts$flowjo_saved_count[[1L]] + 1L
  rewrite_synthetic_count_report(operation3, manifest3, counts)
  expect_error(validate_facs_inputs(config3), "population_count_mismatch")

  operation4 <- withr::local_tempdir(pattern = "SYNTHETIC_operation_")
  write_synthetic_ph3_operation(operation4)
  writeLines("SYNTHETIC MUTATION",
             file.path(operation4, "population_counts.csv"))
  expect_error(validate_facs_inputs(synthetic_production_config(operation4)),
               "count_report_mismatch")
})

test_that("detached digest format and extra acquisition ledgers fail closed", {
  operation <- withr::local_tempdir(pattern = "SYNTHETIC_operation_")
  write_synthetic_ph3_operation(operation)
  config <- synthetic_production_config(operation)
  write(file = file.path(operation, "export-manifest.sha256"),
        "extra", append = TRUE)
  expect_error(validate_facs_inputs(config), "invalid_detached_digest_format")
})

test_that("exact manifest identity, digest, binding, and acquisition metadata is required", {
  cases <- list(
    created_at = list(field = "created_at", reason = "manifest_schema_or_profile_mismatch"),
    identity_semantics = list(field = "identity_semantics", reason = "identity_method_mismatch"),
    digest_target = list(field = "digest_target", reason = "missing_manifest_provenance"),
    binding_object = list(field = "binding_object", reason = "missing_manifest_binding"),
    source_fcs_hash = list(field = "source_fcs_hash", reason = "acquisition_schema_mismatch")
  )
  for (case in cases) {
    operation <- withr::local_tempdir(pattern = "SYNTHETIC_operation_")
    manifest <- write_synthetic_ph3_operation(operation)
    if (case$field == "created_at") manifest$created_at <- NULL
    if (case$field == "identity_semantics") {
      manifest$identity_method$semantics <- "SYNTHETIC-WRONG"
    }
    if (case$field == "digest_target") {
      manifest$manifest_digest$target <- "SYNTHETIC-WRONG"
    }
    if (case$field == "binding_object") {
      manifest$manifest_binding$canonical_object <- "SYNTHETIC-WRONG"
    }
    if (case$field == "source_fcs_hash") {
      manifest$acquisitions[[1L]]$source_fcs_sha256 <- "SYNTHETIC-WRONG"
    }
    rewrite_synthetic_manifest(operation, manifest)
    expect_error(validate_facs_inputs(synthetic_production_config(operation)),
                 case$reason, info = case$field)
  }
})

test_that("source FCS references follow confined Slice 1 Path semantics", {
  bad_references <- c(
    "/SYNTHETIC/input.fcs", "../SYNTHETIC.fcs",
    "SYNTHETIC/../input.fcs",
    "..\\SYNTHETIC.fcs", "SYNTHETIC\\..\\input.fcs",
    "C:\\SYNTHETIC\\input.fcs", "C:SYNTHETIC.fcs",
    "\\\\SYNTHETIC\\input.fcs", "."
  )
  for (reference in bad_references) {
    operation <- withr::local_tempdir(pattern = "SYNTHETIC_operation_")
    manifest <- write_synthetic_ph3_operation(operation)
    manifest$acquisitions[[1L]]$source_fcs_reference <- reference
    rewrite_synthetic_manifest(operation, manifest)
    expect_error(validate_facs_inputs(synthetic_production_config(operation)),
                 "acquisition_schema_mismatch", info = reference)
  }
  harmless_references <- c(
    "SYNTHETIC//input.fcs", "SYNTHETIC/./input.fcs",
    "SYNTHETIC\\.\\input.fcs", "~", "SYNTHETIC"
  )
  expect_true(all(vapply(
    harmless_references, ph3_confined_relative_reference, logical(1)
  )))
})

test_that("required manifest strings reject whitespace-only values", {
  cases <- list(
    approval = list(reason = "missing_approval_metadata", mutate = function(x) {
      x$approval$approver <- "   "
      x
    }),
    provenance = list(reason = "missing_manifest_provenance", mutate = function(x) {
      x$transformation_context <- "\t"
      x
    }),
    acquisition = list(reason = "acquisition_schema_mismatch", mutate = function(x) {
      x$acquisitions[[1L]]$sample_id <- "   "
      x
    }),
    population = list(reason = "population_schema_mismatch", mutate = function(x) {
      x$populations[[1L]]$gate_name <- "   "
      x
    }),
    artifact = list(reason = "artifact_schema_mismatch", mutate = function(x) {
      x$artifacts[[1L]]$gate_path <- "   "
      x
    })
  )
  for (case in cases) {
    operation <- withr::local_tempdir(pattern = "SYNTHETIC_operation_")
    manifest <- case$mutate(write_synthetic_ph3_operation(operation))
    rewrite_synthetic_manifest(operation, manifest)
    expect_error(validate_facs_inputs(synthetic_production_config(operation)),
                 case$reason)
  }
})

test_that("distinct operation directories require unique operation IDs", {
  first <- withr::local_tempdir(pattern = "SYNTHETIC_operation_")
  second <- withr::local_tempdir(pattern = "SYNTHETIC_operation_")
  write_synthetic_ph3_operation(
    first, prefixes = "reference", operation_id = "SYNTHETIC-DUPLICATE-OP"
  )
  write_synthetic_ph3_operation(
    second, prefixes = "treatment", operation_id = "SYNTHETIC-DUPLICATE-OP"
  )
  expect_error(
    validate_facs_inputs(synthetic_production_config(c(first, second))),
    "duplicate_export_operation_id"
  )
})

test_that("unexpected artifact roles and missing nested fields fail with stable reasons", {
  extra <- withr::local_tempdir(pattern = "SYNTHETIC_operation_")
  manifest <- write_synthetic_ph3_operation(extra)
  unexpected <- manifest$artifacts[[length(manifest$artifacts)]]
  unexpected$role <- "SYNTHETIC_UNEXPECTED"
  manifest$artifacts[[length(manifest$artifacts) + 1L]] <- unexpected
  rewrite_synthetic_manifest(extra, manifest)
  expect_error(validate_facs_inputs(synthetic_production_config(extra)),
               "artifact_ledger_coverage_mismatch")

  missing <- withr::local_tempdir(pattern = "SYNTHETIC_operation_")
  manifest <- write_synthetic_ph3_operation(missing)
  manifest$artifacts[[1L]]$gate_path <- NULL
  rewrite_synthetic_manifest(missing, manifest)
  expect_error(validate_facs_inputs(synthetic_production_config(missing)),
               "artifact_schema_mismatch")
})

test_that("relative operation directories cannot escape data_dir", {
  base <- withr::local_tempdir(pattern = "SYNTHETIC_data_")
  config <- synthetic_production_config(base)
  config$ph3_export_operation_dirs <- "../SYNTHETIC_ESCAPE"
  expect_error(ph3_resolve_operation_dirs(config, base), "operation_path_escape")

  config$ph3_export_operation_dirs <- "SYNTHETIC/./operation"
  normalized_base <- normalizePath(base, winslash = "/", mustWork = TRUE)
  resolved <- ph3_resolve_operation_dirs(config, base)
  expect_identical(
    resolved,
    normalizePath(file.path(normalized_base, "SYNTHETIC", "operation"),
                  winslash = "/", mustWork = FALSE)
  )
  expect_false(any(grepl("\\\\", resolved, fixed = TRUE)))
})

test_that("fractional counts, duplicate identity columns, and linkage fail", {
  fractional <- withr::local_tempdir(pattern = "SYNTHETIC_operation_")
  manifest <- write_synthetic_ph3_operation(fractional)
  manifest$artifacts[[1L]]$row_count <- 10.5
  rewrite_synthetic_manifest(fractional, manifest)
  expect_error(validate_facs_inputs(synthetic_production_config(fractional)),
               "artifact_schema_mismatch")

  duplicate <- withr::local_tempdir(pattern = "SYNTHETIC_operation_")
  manifest <- write_synthetic_ph3_operation(duplicate)
  manifest$artifacts[[1L]]$identity_columns <- as.list(c(
    "acquisition_id", "event_index", "event_identity", "event_identity"
  ))
  rewrite_synthetic_manifest(duplicate, manifest)
  expect_error(validate_facs_inputs(synthetic_production_config(duplicate)),
               "artifact_ledger_mismatch")

  reused <- withr::local_tempdir(pattern = "SYNTHETIC_operation_")
  manifest <- write_synthetic_ph3_operation(reused)
  manifest$populations[[2L]]$artifact_path <-
    manifest$populations[[1L]]$artifact_path
  manifest$artifacts[[2L]]$path <- manifest$artifacts[[1L]]$path
  rewrite_synthetic_manifest(reused, manifest)
  expect_error(validate_facs_inputs(synthetic_production_config(reused)),
               "artifact_ledger_coverage_mismatch")

  linkage <- withr::local_tempdir(pattern = "SYNTHETIC_operation_")
  manifest <- write_synthetic_ph3_operation(linkage)
  manifest$populations[[1L]]$export_operation_id <- "SYNTHETIC-WRONG"
  rewrite_synthetic_manifest(linkage, manifest)
  expect_error(validate_facs_inputs(synthetic_production_config(linkage)),
               "acquisition_operation_linkage_mismatch")
})

test_that("missing count report cannot prove an empty or nonempty population", {
  operation <- withr::local_tempdir(pattern = "SYNTHETIC_operation_")
  write_synthetic_ph3_operation(operation)
  unlink(file.path(operation, "population_counts.csv"))
  expect_error(validate_facs_inputs(synthetic_production_config(operation)),
               "missing_count_report")
})

test_that("missing population artifacts and wrong required-channel schemas fail", {
  missing <- withr::local_tempdir(pattern = "SYNTHETIC_operation_")
  manifest <- write_synthetic_ph3_operation(missing)
  unlink(file.path(missing, manifest$populations[[3L]]$artifact_path))
  expect_error(validate_facs_inputs(synthetic_production_config(missing)),
               "missing_population_artifact")

  schema <- withr::local_tempdir(pattern = "SYNTHETIC_operation_")
  manifest <- write_synthetic_ph3_operation(schema)
  path <- file.path(schema, manifest$populations[[1L]]$artifact_path)
  events <- utils::read.csv(path, check.names = FALSE)
  events$raw__Target <- NULL
  manifest$artifacts[[1L]]$columns <- as.list(names(events))
  rewrite_synthetic_population(schema, manifest, 1L, events)
  expect_error(validate_facs_inputs(synthetic_production_config(schema)),
               "configured_channel_schema_mismatch")
})

test_that("required-channel finite-event minima are exactly 10, 2, and 0", {
  complete <- withr::local_tempdir(pattern = "SYNTHETIC_operation_")
  manifest <- write_synthetic_ph3_operation(complete)
  path <- file.path(complete, manifest$populations[[1L]]$artifact_path)
  events <- utils::read.csv(path, check.names = FALSE)
  events$raw__DNA[[1L]] <- NA_real_
  rewrite_synthetic_population(complete, manifest, 1L, events)
  expect_error(validate_facs_inputs(synthetic_production_config(complete)),
               "insufficient_finite_events.*require at least 10")

  g1 <- withr::local_tempdir(pattern = "SYNTHETIC_operation_")
  manifest <- write_synthetic_ph3_operation(g1)
  path <- file.path(g1, manifest$populations[[2L]]$artifact_path)
  events <- utils::read.csv(path, check.names = FALSE)
  events$raw__DNA[[1L]] <- NA_real_
  rewrite_synthetic_population(g1, manifest, 2L, events)
  expect_error(validate_facs_inputs(synthetic_production_config(g1)),
               "insufficient_finite_events.*require at least 2")

  zero <- withr::local_tempdir(pattern = "SYNTHETIC_operation_")
  write_synthetic_ph3_operation(zero, empty_positive = TRUE)
  expect_silent(validate_facs_inputs(synthetic_production_config(zero)))
})

test_that("production validator proves intentional empty and rejects missing proof", {
  valid <- withr::local_tempdir(pattern = "SYNTHETIC_operation_")
  write_synthetic_ph3_operation(valid, empty_positive = TRUE)
  result <- validate_facs_inputs(synthetic_production_config(valid))
  empty_rows <- attr(result, "ph3_containment")$child_population_key ==
    "ph3_positive"
  expect_true(all(attr(result, "ph3_containment")$intentionally_empty[empty_rows]))
  expect_true(all(attr(result, "ph3_containment")$
                  intentionally_empty_proof_status[empty_rows] == "verified"))

  unproven <- withr::local_tempdir(pattern = "SYNTHETIC_operation_")
  manifest <- write_synthetic_ph3_operation(unproven, empty_positive = TRUE)
  positive <- which(vapply(manifest$populations, function(x) {
    identical(x$population_key, "ph3_positive")
  }, logical(1)))[[1L]]
  manifest$populations[[positive]]$intentionally_empty <- FALSE
  manifest$artifacts[[positive]]$intentionally_empty <- FALSE
  rewrite_synthetic_manifest(unproven, manifest)
  expect_error(validate_facs_inputs(synthetic_production_config(unproven)),
               "unproven_empty_population")

  missing_gate <- withr::local_tempdir(pattern = "SYNTHETIC_operation_")
  manifest <- write_synthetic_ph3_operation(missing_gate, empty_positive = TRUE)
  positive <- which(vapply(manifest$populations, function(x) {
    identical(x$population_key, "ph3_positive")
  }, logical(1)))[[1L]]
  manifest$populations[[positive]]$gate_path <- ""
  manifest$artifacts[[positive]]$gate_path <- ""
  rewrite_synthetic_manifest(missing_gate, manifest)
  expect_error(validate_facs_inputs(synthetic_production_config(missing_gate)),
               "population_schema_mismatch|unproven_empty_population")
})

test_that("duplicate and missing or mixed direct row identity fail production", {
  duplicate <- withr::local_tempdir(pattern = "SYNTHETIC_operation_")
  manifest <- write_synthetic_ph3_operation(duplicate)
  path <- file.path(duplicate, manifest$populations[[1L]]$artifact_path)
  events <- utils::read.csv(path, check.names = FALSE, colClasses = "character")
  events$event_identity[[2L]] <- events$event_identity[[1L]]
  events$event_index[[2L]] <- events$event_index[[1L]]
  manifest <- rewrite_synthetic_population(duplicate, manifest, 1L, events)
  expect_error(validate_facs_inputs(synthetic_production_config(duplicate)),
               "duplicate_direct_identity")

  mixed <- withr::local_tempdir(pattern = "SYNTHETIC_operation_")
  manifest <- write_synthetic_ph3_operation(mixed)
  path <- file.path(mixed, manifest$populations[[1L]]$artifact_path)
  events <- utils::read.csv(path, check.names = FALSE, colClasses = "character")
  events$export_operation_id[[1L]] <- "SYNTHETIC-WRONG"
  rewrite_synthetic_population(mixed, manifest, 1L, events)
  expect_error(validate_facs_inputs(synthetic_production_config(mixed)),
               "row_binding_mismatch")

  missing <- withr::local_tempdir(pattern = "SYNTHETIC_operation_")
  manifest <- write_synthetic_ph3_operation(missing)
  path <- file.path(missing, manifest$populations[[1L]]$artifact_path)
  events <- utils::read.csv(path, check.names = FALSE)
  events$event_identity <- NULL
  manifest$artifacts[[1L]]$columns <- as.list(names(events))
  rewrite_synthetic_population(missing, manifest, 1L, events)
  expect_error(validate_facs_inputs(synthetic_production_config(missing)),
               "missing_identity_fields")

  for (field in c("export_profile", "export_operation_id",
                  "export_manifest_digest", "sample_id")) {
    absent <- withr::local_tempdir(pattern = "SYNTHETIC_operation_")
    manifest <- write_synthetic_ph3_operation(absent)
    path <- file.path(absent, manifest$populations[[1L]]$artifact_path)
    events <- utils::read.csv(path, check.names = FALSE)
    events[[field]] <- NULL
    manifest$artifacts[[1L]]$columns <- as.list(names(events))
    rewrite_synthetic_population(absent, manifest, 1L, events)
    expect_error(validate_facs_inputs(synthetic_production_config(absent)),
                 "missing_identity_fields", info = field)
  }

  noncanonical <- withr::local_tempdir(pattern = "SYNTHETIC_operation_")
  manifest <- write_synthetic_ph3_operation(noncanonical)
  path <- file.path(noncanonical, manifest$populations[[1L]]$artifact_path)
  events <- utils::read.csv(path, check.names = FALSE, colClasses = "character")
  events$event_index[[1L]] <- "01"
  events$event_identity[[1L]] <- paste0(events$acquisition_id[[1L]],
                                        ":event_index:01")
  rewrite_synthetic_population(noncanonical, manifest, 1L, events)
  expect_error(validate_facs_inputs(synthetic_production_config(noncanonical)),
               "direct_identity_contract_mismatch")
})

test_that("structured production results and provenance are retained exactly", {
  operation <- withr::local_tempdir(pattern = "SYNTHETIC_operation_")
  write_synthetic_ph3_operation(operation)
  config <- synthetic_production_config(operation)
  report <- validate_facs_inputs(config)
  containment <- attr(report, "ph3_containment")
  expect_equal(nrow(containment), 4L)
  expect_true(all(containment$containment_status == "validated"))
  expect_true(all(containment$containment_reason_code == "contained"))
  expect_true(all(containment$unmatched_child_count == 0L))
  expect_true(all(containment$child_duplicate_row_count == 0L))
  reference_g1 <- containment[containment$prefix == "reference" &
                                containment$child_population_key == "g1", ]
  expect_identical(reference_g1$parent_row_count, 10L)
  expect_identical(reference_g1$parent_unique_identity_count, 10L)
  expect_identical(reference_g1$child_row_count, 2L)
  expect_identical(reference_g1$child_unique_identity_count, 2L)
  expect_identical(reference_g1$matched_child_count, 2L)
  expect_identical(reference_g1$excess_occurrence_count, 0L)
  analysis <- new_facs_analysis(
    config, build_sample_manifest(config), report,
    normalized_data = list(), models = list()
  )
  expect_identical(analysis$provenance$ph3_containment, containment)
  expect_length(analysis$provenance$ph3_export_manifests, 1L)
  expect_null(attr(analysis$input_report, "ph3_verified_events"))
})

test_that("production high-level analysis retains Slice 4 and appends Slice 5 tables", {
  operation <- withr::local_tempdir(pattern = "SYNTHETIC_operation_")
  write_synthetic_ph3_operation(operation)
  analysis <- analyze_facs_experiment(synthetic_production_config(operation))
  expect_s3_class(analysis, "facs_analysis")
  expect_null(attr(analysis$input_report, "ph3_verified_events"))
  expect_equal(nrow(analysis$provenance$ph3_containment), 4L)
  expect_length(analysis$provenance$ph3_export_manifests, 1L)
  expect_length(analysis$normalized_data, 2L)
  expect_true(all(vapply(
    analysis$normalized_data,
    function(x) is.data.frame(x$ph3_event_classification), logical(1)
  )))
  expect_true(all(vapply(
    analysis$normalized_data,
    function(x) nrow(x$ph3_event_classification) == nrow(x$data), logical(1)
  )))
  expect_true(all(vapply(
    analysis$normalized_data,
    function(x) identical(
      x$ph3_event_classification$event_identity, x$data$event_identity
    ), logical(1)
  )))
  approved_names <- c(
    "ph3_metrics_acquisition", "ph3_phase_prevalence",
    "ph3_event_eligibility_qc", "ph3_4n_boundary_sensitivity_qc",
    "ph3_metrics_biological_replicate",
    "ph3_phase_prevalence_biological_replicate",
    "ph3_metrics_condition_summary",
    "ph3_phase_prevalence_condition_summary"
  )
  expect_identical(names(analysis$quantitation), approved_names)
  expect_true(all(vapply(
    analysis$quantitation, is.data.frame, logical(1)
  )))
  expect_identical(nrow(analysis$quantitation$ph3_metrics_acquisition), 10L)
  expect_identical(nrow(analysis$quantitation$ph3_phase_prevalence), 10L)
  expect_identical(
    nrow(analysis$quantitation$ph3_4n_boundary_sensitivity_qc), 8L
  )
  expect_identical(nrow(analysis$quantitation$ph3_event_eligibility_qc), 32L)
  expect_true(all(vapply(
    analysis$quantitation[approved_names[1:4]],
    function(x) identical(unique(x$aggregation_level), "acquisition"),
    logical(1)
  )))
  expect_identical(
    nrow(analysis$quantitation$ph3_metrics_biological_replicate), 10L
  )
  expect_identical(
    nrow(analysis$quantitation$ph3_phase_prevalence_biological_replicate), 10L
  )
  expect_identical(nrow(analysis$quantitation$ph3_metrics_condition_summary), 10L)
  expect_identical(
    nrow(analysis$quantitation$ph3_phase_prevalence_condition_summary), 10L
  )
  expect_true(all(vapply(
    analysis$quantitation[approved_names[5:6]],
    function(x) identical(unique(x$aggregation_level), "biological_replicate"),
    logical(1)
  )))
  expect_true(all(vapply(
    analysis$quantitation[approved_names[7:8]],
    function(x) identical(unique(x$aggregation_level), "condition"),
    logical(1)
  )))
  expected_acquisition_ids <- c("SYNTHETIC-reference", "SYNTHETIC-treatment")
  expected_rows_per_acquisition <- c(
    ph3_metrics_acquisition = 5L,
    ph3_phase_prevalence = 5L,
    ph3_event_eligibility_qc = 16L,
    ph3_4n_boundary_sensitivity_qc = 4L
  )
  for (name in approved_names[1:4]) {
    table_value <- analysis$quantitation[[name]]
    expect_identical(
      sort(unique(table_value$acquisition_id)), expected_acquisition_ids,
      info = name
    )
    acquisition_rows <- table(table_value$acquisition_id)
    expect_identical(
      as.integer(acquisition_rows[expected_acquisition_ids]),
      rep(expected_rows_per_acquisition[[name]], 2L), info = name
    )
  }
  expect_false("ph3" %in% names(analysis$quantitation))
  expect_false(any(grepl(
    "aggregate|condition_stat|csv|plot|report|geometry|compat",
    names(analysis$quantitation), ignore.case = TRUE
  )))
  quantitation_before_rejected_call <- analysis$quantitation
  expect_error(
    quantify_ph3(analysis),
    "legacy_count_only_unverified_v1"
  )
  expect_identical(analysis$quantitation, quantitation_before_rejected_call)
})

test_that("EdU and POI result provenance remains unchanged", {
  for (mode in c("edu", "poi")) {
    config <- validate_facs_config(minimal_config(mode))
    report <- data.frame(path = character(), exists = logical())
    analysis <- new_facs_analysis(
      config, build_sample_manifest(config), report,
      normalized_data = list(), models = list()
    )
    expect_false(any(startsWith(names(analysis$provenance), "ph3_")),
                 info = mode)
  }
})

test_that("missing provenance and malformed channel/linkage ledgers fail", {
  provenance <- withr::local_tempdir(pattern = "SYNTHETIC_operation_")
  manifest <- write_synthetic_ph3_operation(provenance)
  manifest$transformation_context <- NULL
  manifest$manifest_binding$digest <- ph3_manifest_binding_digest(manifest)
  rewrite_synthetic_manifest(provenance, manifest)
  expect_error(validate_facs_inputs(synthetic_production_config(provenance)),
               "manifest_schema_or_profile_mismatch|missing_manifest_provenance")

  channels <- withr::local_tempdir(pattern = "SYNTHETIC_operation_")
  manifest <- write_synthetic_ph3_operation(channels)
  manifest$populations[[1L]]$channels <- as.list(c("DNA", "DNA"))
  manifest$artifacts[[1L]]$channels <- as.list(c("DNA", "DNA"))
  rewrite_synthetic_manifest(channels, manifest)
  expect_error(validate_facs_inputs(synthetic_production_config(channels)),
               "configured_channel_schema_mismatch")

  acquisition <- withr::local_tempdir(pattern = "SYNTHETIC_operation_")
  manifest <- write_synthetic_ph3_operation(acquisition)
  manifest$artifacts[[1L]]$acquisition_id <- "SYNTHETIC-WRONG"
  rewrite_synthetic_manifest(acquisition, manifest)
  expect_error(validate_facs_inputs(synthetic_production_config(acquisition)),
               "population_ledger_coverage_mismatch|artifact_ledger_mismatch")
})
