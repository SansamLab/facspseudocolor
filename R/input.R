# Input manifests, reading, and validation -----------------------------------

#' Build the sample manifest for an analysis
#'
#' Converts the configured samples or replicates into one row per sample while
#' preserving the configured order.
#'
#' @param config A validated `facs_config` object or a configuration list.
#'
#' @return A data frame describing sample, condition, replicate, prefix, and
#'   reference status.
#' @export
build_sample_manifest <- function(config) {
  if (!inherits(config, "facs_config")) config <- validate_facs_config(config)
  manifest <- make_sample_manifest(
    samples = config$samples, replicates = config$replicates
  )
  if (identical(config$plot_type, "ph3") &&
      identical(config$ph3_input_profile,
                "production_direct_identity_v1")) {
    replicate_set_ids <- vapply(
      config$replicates, `[[`, character(1), "id"
    )
    manifest$replicate_set_id <- replicate_set_ids[manifest$replicate_index]
  }
  manifest
}

resolve_facs_directory <- function(config, data_dir = NULL) {
  directory <- if (is.null(data_dir)) config$data_dir else data_dir
  if (!config_scalar_string(directory)) {
    stop("`data_dir` must be one explicit, nonempty directory path.", call. = FALSE)
  }
  directory <- path.expand(directory)
  if (!grepl("^(/|[A-Za-z]:[/\\\\])", directory)) {
    config_dir <- attr(config, "config_dir")
    if (is.null(data_dir) && config_scalar_string(config_dir)) {
      directory <- file.path(config_dir, directory)
    } else if (is.null(data_dir)) {
      stop(
        "Relative `data_dir` requires a configuration read from a file; otherwise supply an explicit `data_dir`.",
        call. = FALSE
      )
    }
  }
  normalizePath(directory, mustWork = FALSE)
}

required_population_keys <- function(config, manifest) {
  if (identical(config$plot_type, "poi")) return("complete")
  if (identical(config$plot_type, "ph3")) {
    return(c("complete", "g1", "ph3_positive"))
  }

  list(
    common = c("complete", "g1"),
    positive_rows = seq_len(nrow(manifest))
  )
}

facs_input_files <- function(config, data_dir = NULL) {
  manifest <- build_sample_manifest(config)
  directory <- resolve_facs_directory(config, data_dir)
  rows <- list()
  add_row <- function(i, key, required = TRUE) {
    suffix <- config$suffixes[[key]]
    data.frame(
      replicate = manifest$replicate[[i]],
      replicate_index = manifest$replicate_index[[i]],
      condition = manifest$condition[[i]],
      prefix = manifest$prefix[[i]],
      population = key,
      required = required,
      path = file.path(directory, paste0(manifest$prefix[[i]], suffix)),
      stringsAsFactors = FALSE
    )
  }

  if (identical(config$plot_type, "poi")) {
    rows <- lapply(seq_len(nrow(manifest)), add_row, key = "complete")
  } else if (identical(config$plot_type, "ph3")) {
    for (i in seq_len(nrow(manifest))) {
      rows[[length(rows) + 1L]] <- add_row(i, "complete")
      rows[[length(rows) + 1L]] <- add_row(i, "g1")
      rows[[length(rows) + 1L]] <- add_row(i, "ph3_positive")
    }
  } else {
    for (i in seq_len(nrow(manifest))) {
      rows[[length(rows) + 1L]] <- add_row(i, "complete")
      rows[[length(rows) + 1L]] <- add_row(i, "g1")
      positive_required <- i %in% required_population_keys(config, manifest)$positive_rows
      rows[[length(rows) + 1L]] <- add_row(i, "edu_positive", positive_required)
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Read and validate one event-level FACS table
#'
#' Reads a CSV path or validates an in-memory data frame. The function does not
#' filter, repair, transform, or substitute event values.
#'
#' @param input An explicit CSV path or an in-memory data frame.
#' @param dna_channel Exact name of the DNA channel column.
#' @param target_channel Exact name of the target channel column.
#' @param sample_id Label used in validation errors and provenance.
#' @param minimum_finite_events Minimum number of rows in which both required
#'   channels are finite.
#'
#' @return The original event table with validation and provenance attributes.
#' @export
read_facs_sample <- function(
    input, dna_channel, target_channel, sample_id = "sample",
    minimum_finite_events = 2L
) {
  if (!config_scalar_string(dna_channel) || !config_scalar_string(target_channel)) {
    stop("`dna_channel` and `target_channel` must be explicit nonempty strings.",
         call. = FALSE)
  }
  if (identical(dna_channel, target_channel)) {
    stop("DNA and target channels must be different columns.", call. = FALSE)
  }
  if (!config_scalar_number(as.numeric(minimum_finite_events)) ||
      minimum_finite_events < 0 || minimum_finite_events %% 1 != 0) {
    stop("`minimum_finite_events` must be a nonnegative integer.", call. = FALSE)
  }

  source_path <- NULL
  if (is.data.frame(input)) {
    events <- input
    input_type <- "data.frame"
  } else if (config_scalar_string(input)) {
    if (!file.exists(input)) {
      stop("Input CSV not found for ", sample_id, ": ", input, call. = FALSE)
    }
    source_path <- normalizePath(input, mustWork = TRUE)
    events <- utils::read.csv(source_path, check.names = FALSE)
    input_type <- "csv"
  } else {
    stop("`input` must be an explicit CSV path or a data frame.", call. = FALSE)
  }

  missing_columns <- setdiff(c(dna_channel, target_channel), names(events))
  if (length(missing_columns)) {
    stop(
      "Missing required channel(s) for ", sample_id, ": ",
      paste(missing_columns, collapse = ", "), ". Available columns: ",
      paste(names(events), collapse = ", "), call. = FALSE
    )
  }
  valid_required_type <- function(value) {
    is.numeric(value) || (nrow(events) == 0L && is.logical(value))
  }
  nonnumeric <- c(dna_channel, target_channel)[
    !vapply(events[c(dna_channel, target_channel)],
            valid_required_type, logical(1))
  ]
  if (length(nonnumeric)) {
    stop("Required channel(s) must be numeric for ", sample_id, ": ",
         paste(nonnumeric, collapse = ", "), call. = FALSE)
  }

  finite <- is.finite(events[[dna_channel]]) & is.finite(events[[target_channel]])
  finite_n <- sum(finite)
  if (finite_n < minimum_finite_events) {
    stop(
      "Too few finite events for ", sample_id, ": found ", finite_n,
      ", require at least ", minimum_finite_events, ".", call. = FALSE
    )
  }

  attr(events, "facs_input") <- list(
    sample_id = sample_id,
    input_type = input_type,
    source_path = source_path,
    event_n = nrow(events),
    finite_event_n = finite_n,
    nonfinite_event_n = nrow(events) - finite_n,
    dna_channel = dna_channel,
    target_channel = target_channel
  )
  events
}

#' Validate all required input CSVs for an analysis
#'
#' Resolves the configured sample files, checks that every required file exists,
#' and validates required channel columns and numeric types. Event values are
#' never filtered or repaired.
#'
#' @param config A validated `facs_config` object or configuration list.
#' @param data_dir Optional explicit data directory overriding `config$data_dir`.
#'
#' @return A data frame with one row per expected population file, including
#'   existence and event-count diagnostics.
#' @export
validate_facs_inputs <- function(config, data_dir = NULL) {
  if (!inherits(config, "facs_config")) config <- validate_facs_config(config)
  if (identical(config$plot_type, "ph3") &&
      identical(config$ph3_input_profile, "production_direct_identity_v1")) {
    return(validate_ph3_export_operations(config, data_dir))
  }
  files <- facs_input_files(config, data_dir)
  files$exists <- file.exists(files$path)

  missing <- files$required & !files$exists
  if (any(missing)) {
    details <- paste0(
      files$replicate[missing], " / ", files$condition[missing], " / ",
      files$population[missing], ": ", files$path[missing]
    )
    stop(
      paste0("Missing required FACS input file(s):\n- ",
             paste(details, collapse = "\n- ")), call. = FALSE
    )
  }

  files$event_n <- NA_integer_
  files$finite_event_n <- NA_integer_
  files$nonfinite_event_n <- NA_integer_
  files$subset_membership_validated <- NA
  for (i in which(files$exists)) {
    minimum <- if (files$population[[i]] == "complete") {
      10L
    } else if (files$population[[i]] == "ph3_positive") {
      0L
    } else {
      2L
    }
    events <- read_facs_sample(
      files$path[[i]], config$dna_channel, config$target_channel,
      sample_id = paste(files$replicate[[i]], files$condition[[i]],
                        files$population[[i]], sep = " / "),
      minimum_finite_events = minimum
    )
    info <- attr(events, "facs_input")
    files$event_n[[i]] <- info$event_n
    files$finite_event_n[[i]] <- info$finite_event_n
    files$nonfinite_event_n[[i]] <- info$nonfinite_event_n
  }
  if (identical(config$plot_type, "ph3")) {
    warning(
      "LEGACY COUNT-ONLY PH3 INPUT: event identity and population containment are unverified.",
      call. = FALSE
    )
    grouped <- split(files, files$prefix)
    count_errors <- vapply(grouped, function(x) {
      complete_n <- x$event_n[x$population == "complete"]
      positive_n <- x$event_n[x$population == "ph3_positive"]
      length(complete_n) == 1L && length(positive_n) == 1L &&
        is.finite(complete_n) && is.finite(positive_n) &&
        positive_n > complete_n
    }, logical(1))
    if (any(count_errors)) {
      stop(
        "pH3-positive event count exceeds the Single Cell count for: ",
        paste(names(grouped)[count_errors], collapse = ", "),
        ". Verify that the exported pH3 gate is nested within Single Cells.",
        call. = FALSE
      )
    }
    for (prefix in names(grouped)) {
      group <- grouped[[prefix]]
      complete_path <- group$path[group$population == "complete"]
      positive_path <- group$path[group$population == "ph3_positive"]
      complete <- utils::read.csv(complete_path, check.names = FALSE)
      positive <- utils::read.csv(positive_path, check.names = FALSE)
      if ("event_index" %in% names(complete) &&
          "event_index" %in% names(positive)) {
        invalid_index <- function(value) {
          anyNA(value) || any(!nzchar(trimws(as.character(value)))) ||
            anyDuplicated(value)
        }
        if (invalid_index(complete$event_index) ||
            invalid_index(positive$event_index)) {
          stop(
            "Invalid, missing, or duplicate event_index values for legacy PH3 sample ",
            prefix, ".", call. = FALSE
          )
        }
        outside <- setdiff(positive$event_index, complete$event_index)
        if (length(outside)) {
          stop(
            "The legacy pH3-positive population contains event_index values outside ",
            "the Single Cell population for ", prefix,
            ". Verify the FlowJo gate hierarchy.", call. = FALSE
          )
        }
      }
    }
    child <- files$population %in% c("g1", "ph3_positive")
    files$subset_membership_validated[child] <- FALSE
    legacy_rows <- lapply(names(grouped), function(prefix) {
      group <- grouped[[prefix]]
      parent <- group[group$population == "complete", , drop = FALSE]
      do.call(rbind, lapply(c("g1", "ph3_positive"), function(key) {
        child_row <- group[group$population == key, , drop = FALSE]
        data.frame(
          containment_schema_version = "ph3-input-containment-1.0.0",
          acquisition_id = NA_character_, sample_id = NA_character_,
          prefix = prefix, parent_population_key = "complete",
          parent_population_path = parent$path[[1L]],
          child_population_key = key,
          child_population_path = child_row$path[[1L]],
          export_profile = "legacy_count_only_unverified_v1",
          export_operation_id = NA_character_, manifest_digest = NA_character_,
          manifest_reference = NA_character_, identity_method_id = NA_character_,
          identity_method_version = NA_character_, identity_source = NA_character_,
          parent_row_count = parent$event_n[[1L]],
          parent_unique_identity_count = NA_integer_,
          child_row_count = child_row$event_n[[1L]],
          child_unique_identity_count = NA_integer_,
          parent_duplicated_identity_count = NA_integer_,
          parent_duplicate_row_count = NA_integer_,
          child_duplicated_identity_count = NA_integer_,
          child_duplicate_row_count = NA_integer_, matched_child_count = NA_integer_,
          unmatched_child_count = NA_integer_, excess_occurrence_count = NA_integer_,
          intentionally_empty = identical(child_row$event_n[[1L]], 0L),
          intentionally_empty_proof_status = "unverified_legacy",
          containment_method_id = "unverified_count_only_legacy",
          containment_method_version = "1.0.0", containment_status = "unverified",
          containment_reason_code = "legacy_identity_unavailable",
          validation_severity = "warning", parent_artifact_path = parent$path[[1L]],
          parent_artifact_sha256 = NA_character_,
          child_artifact_path = child_row$path[[1L]],
          child_artifact_sha256 = NA_character_, stringsAsFactors = FALSE
        )
      }))
    })
    attr(files, "ph3_containment") <- do.call(rbind, legacy_rows)
  }
  files
}

ph3_required_identity_columns <- function() {
  c(
    "acquisition_id", "event_index", "event_identity", "identity_source",
    "identity_method_id", "identity_method_version", "duplicate_occurrence",
    "export_profile", "export_operation_id", "export_manifest_digest",
    "export_manifest_reference"
  )
}

ph3_sha256_file <- function(path) {
  connection <- file(path, open = "rb")
  on.exit(close(connection), add = TRUE)
  digest <- as.character(openssl::sha256(connection))
  attributes(digest) <- NULL
  digest
}

ph3_canonical_json <- function(value) {
  if (is.null(value)) return("null")
  if (is.list(value)) {
    if (!is.null(names(value))) {
      keys <- sort(names(value))
      members <- vapply(keys, function(key) {
        paste0(ph3_canonical_json(key), ":", ph3_canonical_json(value[[key]]))
      }, character(1))
      return(paste0("{", paste(members, collapse = ","), "}"))
    }
    members <- vapply(value, ph3_canonical_json, character(1))
    return(paste0("[", paste(members, collapse = ","), "]"))
  }
  if (length(value) != 1L || is.na(value) ||
      (is.numeric(value) && !is.finite(value))) {
    stop("Manifest canonicalization encountered an invalid scalar.", call. = FALSE)
  }
  if (is.logical(value)) return(if (value) "true" else "false")
  as.character(jsonlite::toJSON(value, auto_unbox = TRUE, digits = NA,
                                null = "null", na = "null"))
}

ph3_manifest_binding_digest <- function(manifest) {
  fields <- c(
    "manifest_schema", "export_operation_id", "profile", "export_profile",
    "identity_method", "workspace", "software", "transformation_context",
    "compensation_context", "approval", "acquisitions",
    "requested_populations"
  )
  object <- manifest[fields]
  bytes <- charToRaw(paste0(ph3_canonical_json(object), "\n"))
  digest <- as.character(openssl::sha256(bytes))
  attributes(digest) <- NULL
  digest
}

ph3_nonnegative_integer <- function(value) {
  is.numeric(value) && length(value) == 1L && is.finite(value) &&
    value >= 0 && value %% 1 == 0
}

ph3_trimmed_nonblank_string <- function(value) {
  config_scalar_string(value) && nzchar(trimws(value))
}

ph3_exact_object <- function(value, fields) {
  is.list(value) && !is.data.frame(value) && !is.null(names(value)) &&
    !anyDuplicated(names(value)) && setequal(names(value), fields)
}

ph3_string_array <- function(value, allow_empty = FALSE) {
  is.list(value) && is.null(names(value)) &&
    (allow_empty || length(value) > 0L) &&
    all(vapply(value, ph3_trimmed_nonblank_string, logical(1)))
}

ph3_record_array <- function(value, allow_empty = FALSE) {
  is.list(value) && is.null(names(value)) &&
    (allow_empty || length(value) > 0L)
}

ph3_lower_sha256 <- function(value) {
  config_scalar_string(value) && grepl("^[0-9a-f]{64}$", value)
}

ph3_confined_relative_reference <- function(value) {
  if (!ph3_trimmed_nonblank_string(value)) return(FALSE)
  portable <- gsub("\\\\", "/", value)
  if (startsWith(portable, "/") || grepl("^[A-Za-z]:", portable)) {
    return(FALSE)
  }
  components <- strsplit(portable, "/", fixed = TRUE)[[1L]]
  if (any(components == "..")) return(FALSE)
  components <- components[!components %in% c("", ".")]
  length(components) > 0L
}

ph3_fail <- function(acquisition, population, reason, detail) {
  stop(
    "PH3 input validation failed [", reason, "] for acquisition '",
    acquisition, "', population '", population, "': ", detail,
    call. = FALSE
  )
}

ph3_resolve_operation_dirs <- function(config, data_dir = NULL) {
  directories <- unlist(config$ph3_export_operation_dirs, use.names = FALSE)
  base <- normalizePath(resolve_facs_directory(config, data_dir),
                        winslash = "/", mustWork = FALSE)
  resolved <- vapply(directories, function(path) {
    expanded <- path.expand(path)
    portable <- gsub("\\\\", "/", expanded)
    relative <- !startsWith(portable, "/") &&
      !grepl("^[A-Za-z]:/", portable)
    if (relative) {
      components <- strsplit(portable, "/", fixed = TRUE)[[1L]]
      if (any(components == "..")) {
        ph3_fail("unknown", "manifest", "operation_path_escape",
                 "relative export-operation directories cannot contain parent traversal")
      }
      components <- components[!components %in% c("", ".")]
      expanded <- if (length(components)) {
        file.path(base, paste(components, collapse = .Platform$file.sep))
      } else {
        base
      }
    }
    path <- normalizePath(expanded, winslash = "/", mustWork = FALSE)
    if (relative && !identical(path, base) &&
        !startsWith(path, paste0(base, "/"))) {
      ph3_fail("unknown", "manifest", "operation_path_escape",
               "relative export-operation directories must remain beneath data_dir")
    }
    path
  }, character(1), USE.NAMES = FALSE)
  if (anyDuplicated(resolved)) {
    ph3_fail("unknown", "manifest", "duplicate_operation_directory",
             "each explicit export-operation directory must be unique")
  }
  resolved
}

ph3_operation_member <- function(
    operation, relative, acquisition, population, reason
) {
  if (!ph3_trimmed_nonblank_string(relative) ||
      grepl("^(/|[A-Za-z]:[/\\\\])", path.expand(relative))) {
    ph3_fail(acquisition, population, reason,
             "ledger member path must be one nonempty relative path")
  }
  path <- normalizePath(file.path(operation$directory, relative),
                        winslash = "/", mustWork = FALSE)
  root <- paste0(operation$directory, "/")
  if (!startsWith(path, root) || !file.exists(path) ||
      !utils::file_test("-f", path)) {
    ph3_fail(acquisition, population, reason,
             "ledger member is absent, is not a regular file, or escapes the operation")
  }
  path
}

ph3_validate_manifest_structure <- function(manifest) {
  required <- c(
    "manifest_schema", "export_operation_id", "created_at", "status", "profile",
    "export_profile", "identity_method", "workspace", "software", "approval",
    "acquisitions", "requested_populations", "populations", "artifacts",
    "geometry_overlay_status", "manifest_binding", "manifest_digest",
    "transformation_context", "compensation_context"
  )
  if (!is.list(manifest) || is.null(names(manifest)) || anyDuplicated(names(manifest))) {
    ph3_fail("unknown", "manifest", "manifest_schema_or_profile_mismatch",
             "finalized manifest must be one named JSON object")
  }
  missing <- setdiff(required, names(manifest))
  if (length(missing) || !setequal(names(manifest), c(required, "legacy_warning")) ||
      !is.null(manifest$legacy_warning) ||
      !ph3_trimmed_nonblank_string(manifest$export_operation_id) ||
      grepl("[/\\\\]", manifest$export_operation_id) ||
      !ph3_trimmed_nonblank_string(manifest$created_at) ||
      !identical(manifest$status, "complete") ||
      !identical(manifest$profile, "production_direct_identity_v1") ||
      !ph3_exact_object(manifest$manifest_schema, c("name", "version")) ||
      !identical(manifest$manifest_schema$name, "facspseudocolor-flowjo-export") ||
      !identical(manifest$manifest_schema$version, "1.0.0") ||
      !ph3_record_array(manifest$acquisitions) ||
      !ph3_record_array(manifest$populations) ||
      !ph3_record_array(manifest$artifacts)) {
    ph3_fail("unknown", "manifest", "manifest_schema_or_profile_mismatch",
             paste(c("missing", missing), collapse = " "))
  }

  if (!ph3_exact_object(manifest$export_profile, c("id", "version")) ||
      !identical(manifest$export_profile$id,
                 "ph3_export_identity_provenance_v1") ||
      !identical(manifest$export_profile$version, "1.0.0") ||
      !ph3_exact_object(manifest$identity_method,
                        c("id", "version", "source", "semantics", "verified")) ||
      !identical(manifest$identity_method$id,
                 "flowkit_source_event_index_scoped_v1") ||
      !identical(manifest$identity_method$version, "1.0.0") ||
      !identical(manifest$identity_method$source,
                 "flowkit_get_gate_events_index") ||
      !identical(manifest$identity_method$semantics,
                 "conditional_on_pinned_environment_verification") ||
      !isTRUE(manifest$identity_method$verified)) {
    ph3_fail("unknown", "manifest", "identity_method_mismatch",
             "production direct identity metadata does not match Slice 1 exactly")
  }

  workspace_fields <- c("filename", "local_path", "sha256")
  software_fields <- c(
    "exporter", "exporter_version", "source_commit", "flowjo_version",
    "flowkit_version", "supported_flowkit_version", "python_version"
  )
  digest_fields <- c("algorithm", "target", "storage")
  if (!ph3_exact_object(manifest$workspace, workspace_fields) ||
      any(!vapply(manifest$workspace[c("filename", "local_path")],
                  ph3_trimmed_nonblank_string, logical(1))) ||
      !ph3_lower_sha256(manifest$workspace$sha256) ||
      !ph3_exact_object(manifest$software, software_fields) ||
      any(!vapply(manifest$software, ph3_trimmed_nonblank_string, logical(1))) ||
      !identical(manifest$software$flowkit_version,
                 manifest$software$supported_flowkit_version) ||
      !ph3_trimmed_nonblank_string(manifest$transformation_context) ||
      !ph3_trimmed_nonblank_string(manifest$compensation_context) ||
      !ph3_exact_object(manifest$manifest_digest, digest_fields) ||
      !identical(manifest$manifest_digest$algorithm, "sha256") ||
      !identical(manifest$manifest_digest$target,
                 "canonical UTF-8 export-manifest.json bytes") ||
      !identical(manifest$manifest_digest$storage,
                 "export-manifest.sha256 sidecar") ||
      !identical(manifest$geometry_overlay_status, "not_requested")) {
    ph3_fail("unknown", "manifest", "missing_manifest_provenance",
             "workspace/software/transform/compensation/digest/geometry metadata is not exact")
  }

  approval_fields <- c(
    "gate_owner", "approver", "approval_date", "approval_record",
    "positivity_method_id", "positivity_method_version"
  )
  if (!ph3_exact_object(manifest$approval, approval_fields) ||
      any(!vapply(manifest$approval, ph3_trimmed_nonblank_string, logical(1)))) {
    ph3_fail("unknown", "manifest", "missing_approval_metadata",
             "production gate/positivity approval metadata is incomplete")
  }

  binding_fields <- c("algorithm", "canonical_object", "excludes", "digest")
  expected_excludes <- c("created_at", "populations", "artifacts",
                         "full_manifest_digest")
  if (!ph3_exact_object(manifest$manifest_binding, binding_fields) ||
      !identical(manifest$manifest_binding$algorithm, "sha256") ||
      !identical(manifest$manifest_binding$canonical_object,
                 "manifest_binding_object_v1") ||
      !ph3_string_array(manifest$manifest_binding$excludes) ||
      !identical(unlist(manifest$manifest_binding$excludes, use.names = FALSE),
                 expected_excludes) ||
      !ph3_lower_sha256(manifest$manifest_binding$digest)) {
    ph3_fail("unknown", "manifest", "missing_manifest_binding",
             "manifest binding metadata does not match Slice 1 exactly")
  }

  requested <- unlist(manifest$requested_populations, use.names = FALSE)
  if (!ph3_string_array(manifest$requested_populations) ||
      !identical(requested, c("complete", "g1", "ph3_positive"))) {
    ph3_fail("unknown", "manifest", "requested_population_mismatch",
             "production PH3 requires ordered complete, g1, and ph3_positive populations")
  }

  acquisition_fields <- c(
    "acquisition_id", "sample_id", "prefix", "source_fcs_reference",
    "source_fcs_sha256"
  )
  acquisition_ok <- vapply(manifest$acquisitions, function(acquisition) {
    ph3_exact_object(acquisition, acquisition_fields) &&
      all(vapply(acquisition, ph3_trimmed_nonblank_string, logical(1))) &&
      grepl("^[A-Za-z0-9][A-Za-z0-9._-]*$", acquisition$prefix) &&
      !acquisition$prefix %in% c(".", "..") &&
      ph3_confined_relative_reference(acquisition$source_fcs_reference) &&
      ph3_lower_sha256(acquisition$source_fcs_sha256)
  }, logical(1))
  acquisition_ids <- if (all(acquisition_ok)) {
    vapply(manifest$acquisitions, `[[`, character(1), "acquisition_id")
  } else character()
  sample_ids <- if (all(acquisition_ok)) {
    vapply(manifest$acquisitions, `[[`, character(1), "sample_id")
  } else character()
  prefixes <- if (all(acquisition_ok)) {
    vapply(manifest$acquisitions, `[[`, character(1), "prefix")
  } else character()
  if (!all(acquisition_ok) || anyDuplicated(acquisition_ids) ||
      anyDuplicated(sample_ids) || anyDuplicated(prefixes)) {
    bad <- if (any(!acquisition_ok)) which(!acquisition_ok)[[1L]] else 1L
    acquisition <- manifest$acquisitions[[bad]]
    acquisition_id <- if (is.list(acquisition) &&
                              ph3_trimmed_nonblank_string(acquisition$acquisition_id)) {
      acquisition$acquisition_id
    } else "unknown"
    ph3_fail(acquisition_id, "manifest", "acquisition_schema_mismatch",
             "acquisition fields, source FCS hashes, IDs, samples, and prefixes must be exact and unique")
  }

  population_fields <- c(
    "population_key", "gate_name", "gate_path", "gate_type",
    "parent_population_path", "acquisition_id", "sample_id", "prefix",
    "channels", "gate_channels", "row_count", "identity_field",
    "identity_method_id", "unique_identity_count",
    "duplicate_base_combination_count", "duplicate_row_count",
    "intentionally_empty", "export_operation_id", "artifact_path",
    "artifact_sha256"
  )
  population_ok <- vapply(manifest$populations, function(population) {
    text_fields <- c(
      "population_key", "gate_name", "gate_path", "gate_type",
      "parent_population_path", "acquisition_id", "sample_id", "prefix",
      "identity_field", "identity_method_id", "export_operation_id",
      "artifact_path", "artifact_sha256"
    )
    ph3_exact_object(population, population_fields) &&
      all(vapply(population[text_fields], ph3_trimmed_nonblank_string,
                 logical(1))) &&
      ph3_string_array(population$channels) &&
      ph3_string_array(population$gate_channels) &&
      ph3_nonnegative_integer(population$row_count) &&
      ph3_nonnegative_integer(population$unique_identity_count) &&
      ph3_nonnegative_integer(population$duplicate_base_combination_count) &&
      ph3_nonnegative_integer(population$duplicate_row_count) &&
      (identical(population$intentionally_empty, TRUE) ||
       identical(population$intentionally_empty, FALSE)) &&
      ph3_lower_sha256(population$artifact_sha256)
  }, logical(1))
  if (!all(population_ok)) {
    population <- manifest$populations[[which(!population_ok)[[1L]]]]
    acquisition_id <- if (is.list(population) &&
                              ph3_trimmed_nonblank_string(population$acquisition_id)) {
      population$acquisition_id
    } else "unknown"
    population_key <- if (is.list(population) &&
                              ph3_trimmed_nonblank_string(population$population_key)) {
      population$population_key
    } else "manifest"
    ph3_fail(acquisition_id, population_key, "population_schema_mismatch",
             "one or more production population ledger entries are incomplete or malformed")
  }

  artifact_fields <- c(
    "role", "path", "sha256", "byte_size", "row_count", "columns",
    "identity_columns", "identity_method_id", "intentionally_empty",
    "export_operation_id", "acquisition_id", "sample_id", "population_key",
    "gate_path", "channels"
  )
  artifact_ok <- vapply(manifest$artifacts, function(artifact) {
    base_ok <- ph3_exact_object(artifact, artifact_fields) &&
      ph3_trimmed_nonblank_string(artifact$role) &&
      ph3_trimmed_nonblank_string(artifact$path) &&
      ph3_lower_sha256(artifact$sha256) &&
      ph3_nonnegative_integer(artifact$byte_size) &&
      ph3_nonnegative_integer(artifact$row_count) &&
      ph3_string_array(artifact$columns) &&
      ph3_string_array(artifact$identity_columns, allow_empty = TRUE) &&
      ph3_string_array(artifact$channels, allow_empty = TRUE) &&
      (identical(artifact$intentionally_empty, TRUE) ||
       identical(artifact$intentionally_empty, FALSE)) &&
      ph3_trimmed_nonblank_string(artifact$export_operation_id)
    if (!base_ok) return(FALSE)
    if (identical(artifact$role, "population_events")) {
      linked <- c("identity_method_id", "acquisition_id", "sample_id",
                  "population_key", "gate_path")
      return(all(vapply(artifact[linked], ph3_trimmed_nonblank_string,
                        logical(1))) &&
               length(artifact$identity_columns) > 0L &&
               length(artifact$channels) > 0L)
    }
    TRUE
  }, logical(1))
  if (!all(artifact_ok)) {
    artifact <- manifest$artifacts[[which(!artifact_ok)[[1L]]]]
    acquisition_id <- if (is.list(artifact) &&
                              ph3_trimmed_nonblank_string(artifact$acquisition_id)) {
      artifact$acquisition_id
    } else "unknown"
    population_key <- if (is.list(artifact) &&
                              ph3_trimmed_nonblank_string(artifact$population_key)) {
      artifact$population_key
    } else "manifest"
    ph3_fail(acquisition_id, population_key, "artifact_schema_mismatch",
             "one or more artifact ledger entries are incomplete or malformed")
  }
  roles <- vapply(manifest$artifacts, `[[`, character(1), "role")
  artifact_paths <- vapply(manifest$artifacts, `[[`, character(1), "path")
  population_paths <- vapply(manifest$populations, `[[`, character(1),
                             "artifact_path")
  if (sum(roles == "population_events") != length(manifest$acquisitions) * 3L ||
      sum(roles == "population_count_report") != 1L ||
      any(!roles %in% c("population_events", "population_count_report")) ||
      anyDuplicated(artifact_paths) || anyDuplicated(population_paths)) {
    ph3_fail("unknown", "manifest", "artifact_ledger_coverage_mismatch",
             "base Slice 1 manifest must contain only exact population artifacts and one count report")
  }
  invisible(manifest)
}

ph3_read_verified_operation <- function(operation_dir) {
  manifest_path <- file.path(operation_dir, "export-manifest.json")
  digest_path <- file.path(operation_dir, "export-manifest.sha256")
  if (!file.exists(manifest_path) || !file.exists(digest_path)) {
    ph3_fail("unknown", "manifest", "missing_finalized_manifest",
             paste0("required manifest/digest pair is absent in ", operation_dir))
  }
  operation_root <- normalizePath(operation_dir, winslash = "/", mustWork = TRUE)
  manifest_path <- normalizePath(manifest_path, winslash = "/", mustWork = TRUE)
  digest_path <- normalizePath(digest_path, winslash = "/", mustWork = TRUE)
  root_prefix <- paste0(operation_root, "/")
  if (!startsWith(manifest_path, root_prefix) ||
      !startsWith(digest_path, root_prefix) ||
      !utils::file_test("-f", manifest_path) ||
      !utils::file_test("-f", digest_path)) {
    ph3_fail("unknown", "manifest", "manifest_path_escape",
             "manifest/digest must be regular files inside the explicit operation")
  }
  digest_lines <- readLines(digest_path, warn = FALSE)
  if (length(digest_lines) != 1L ||
      !grepl("^[0-9a-f]{64}  export-manifest[.]json$", digest_lines)) {
    ph3_fail("unknown", "manifest", "invalid_detached_digest_format",
             "detached digest must contain exactly one canonical line")
  }
  digest_tokens <- strsplit(digest_lines, "  ", fixed = TRUE)[[1L]]
  if (!identical(digest_tokens[[1L]], ph3_sha256_file(manifest_path))) {
    ph3_fail("unknown", "manifest", "manifest_digest_mismatch",
             "detached digest does not verify the exact manifest bytes")
  }
  manifest <- tryCatch(
    jsonlite::fromJSON(manifest_path, simplifyVector = FALSE),
    error = function(error) {
      ph3_fail("unknown", "manifest", "invalid_manifest_json",
               conditionMessage(error))
    }
  )
  canonical_manifest <- tryCatch(
    charToRaw(paste0(ph3_canonical_json(manifest), "\n")),
    error = function(error) {
      ph3_fail("unknown", "manifest", "invalid_manifest_json",
               conditionMessage(error))
    }
  )
  connection <- file(manifest_path, open = "rb")
  on.exit(close(connection), add = TRUE)
  actual_manifest <- readBin(connection, what = "raw", n = file.info(manifest_path)$size)
  if (!identical(actual_manifest, canonical_manifest)) {
    ph3_fail("unknown", "manifest", "noncanonical_manifest_bytes",
             "manifest is not exact sorted compact UTF-8 JSON with one trailing newline")
  }
  ph3_validate_manifest_structure(manifest)
  binding <- manifest$manifest_binding$digest
  if (!identical(binding, ph3_manifest_binding_digest(manifest))) {
    ph3_fail("unknown", "manifest", "manifest_binding_digest_mismatch",
             "manifest_binding_object_v1 does not match its SHA-256")
  }
  if (!identical(digest_tokens[[1L]], ph3_sha256_file(manifest_path))) {
    ph3_fail("unknown", "manifest", "manifest_changed_during_validation",
             "manifest bytes changed between digest verification and consumption")
  }
  list(
    directory = operation_root, manifest_path = manifest_path,
    digest_path = digest_path,
    full_digest = tolower(digest_tokens[[1L]]), manifest = manifest
  )
}

ph3_verify_count_report <- function(operation) {
  records <- Filter(function(x) identical(x$role, "population_count_report"),
                    operation$manifest$artifacts)
  if (length(records) != 1L) {
    ph3_fail("unknown", "population_counts", "missing_count_report",
             "manifest must ledger exactly one population count report")
  }
  record <- records[[1L]]
  path <- ph3_operation_member(
    operation, record$path, "unknown", "population_counts",
    "missing_count_report"
  )
  if (!identical(as.numeric(file.info(path)$size),
                 as.numeric(record$byte_size)) ||
      !identical(ph3_sha256_file(path), record$sha256)) {
    ph3_fail("unknown", "population_counts", "count_report_mismatch",
             "count report size or hash disagrees before parsing")
  }
  report <- tryCatch(
    utils::read.csv(
      path, check.names = FALSE, stringsAsFactors = FALSE,
      colClasses = c(sample_id = "character", gate_name = "character")
    ),
    error = function(error) {
      ph3_fail("unknown", "population_counts", "count_report_mismatch",
               conditionMessage(error))
    }
  )
  required <- c("sample_id", "gate_name", "exported_count",
                "flowjo_saved_count", "count_difference")
  if (!identical(unlist(record$columns, use.names = FALSE), required) ||
      !identical(names(report), required) ||
      !identical(nrow(report), record$row_count) ||
      !identical(as.numeric(file.info(path)$size), as.numeric(record$byte_size)) ||
      !identical(ph3_sha256_file(path), record$sha256) ||
      !identical(record$export_operation_id,
                 operation$manifest$export_operation_id) ||
      length(record$identity_columns) != 0L ||
      !is.null(record$identity_method_id) ||
      !identical(record$intentionally_empty, FALSE) ||
      !is.null(record$acquisition_id) || !is.null(record$sample_id) ||
      !is.null(record$population_key) || !is.null(record$gate_path) ||
      length(record$channels) != 0L) {
    ph3_fail("unknown", "population_counts", "count_report_mismatch",
             "count report schema, count, size, hash, or operation linkage disagrees")
  }
  count_values <- c(report$exported_count, report$flowjo_saved_count,
                    report$count_difference)
  if (any(!vapply(count_values, ph3_nonnegative_integer, logical(1))) ||
      any(report$count_difference != 0) ||
      any(report$exported_count - report$flowjo_saved_count !=
          report$count_difference)) {
    ph3_fail("unknown", "population_counts", "count_report_mismatch",
             "FlowJo saved and exported counts do not reconcile")
  }
  list(report = report, record = record, path = path)
}

ph3_verify_population_artifact <- function(
    operation, prefix, population_key, count_report
) {
  manifest <- operation$manifest
  populations <- Filter(function(x) {
    identical(x$prefix, prefix) && identical(x$population_key, population_key)
  }, manifest$populations)
  if (length(populations) != 1L) {
    ph3_fail(prefix, population_key, "ambiguous_population_mapping",
             "manifest must contain exactly one population ledger entry")
  }
  population <- populations[[1L]]
  artifacts <- Filter(function(x) {
    identical(x$role, "population_events") &&
      identical(x$population_key, population_key) &&
      identical(x$acquisition_id, population$acquisition_id) &&
      identical(x$sample_id, population$sample_id)
  }, manifest$artifacts)
  if (length(artifacts) != 1L) {
    ph3_fail(population$acquisition_id, population_key, "missing_or_ambiguous_artifact",
             "manifest must contain exactly one linked population artifact")
  }
  artifact <- artifacts[[1L]]
  path <- ph3_operation_member(
    operation, artifact$path, population$acquisition_id, population_key,
    "missing_population_artifact"
  )
  if (!identical(as.numeric(file.info(path)$size),
                 as.numeric(artifact$byte_size)) ||
      !identical(ph3_sha256_file(path), artifact$sha256)) {
    ph3_fail(population$acquisition_id, population_key,
             "artifact_ledger_mismatch",
             "population artifact size or hash disagrees before parsing")
  }
  expected_text <- unique(c(
    "export_operation_id", "export_manifest_digest", "acquisition_id",
    "sample_id", "export_profile", ph3_required_identity_columns()
  ))
  header <- tryCatch(
    names(utils::read.csv(path, check.names = FALSE, nrows = 0L)),
    error = function(error) {
      ph3_fail(population$acquisition_id, population_key,
               "population_artifact_unreadable", conditionMessage(error))
    }
  )
  missing <- setdiff(expected_text, header)
  if (length(missing)) {
    ph3_fail(population$acquisition_id, population_key, "missing_identity_fields",
             paste(missing, collapse = ", "))
  }
  events <- tryCatch(
    utils::read.csv(
      path, check.names = FALSE,
      colClasses = stats::setNames(rep("character", length(expected_text)),
                                   expected_text)
    ),
    error = function(error) {
      ph3_fail(population$acquisition_id, population_key,
               "population_artifact_unreadable", conditionMessage(error))
    }
  )
  identity_columns <- unlist(artifact$identity_columns, use.names = FALSE)
  artifact_columns <- unlist(artifact$columns, use.names = FALSE)
  population_channels <- unlist(population$channels, use.names = FALSE)
  artifact_channels <- unlist(artifact$channels, use.names = FALSE)
  if (anyDuplicated(artifact_columns) ||
      !identical(names(events), artifact_columns) ||
      !identical(nrow(events), artifact$row_count) ||
      !identical(as.numeric(file.info(path)$size), as.numeric(artifact$byte_size)) ||
      !identical(ph3_sha256_file(path), artifact$sha256) ||
      !identical(population$artifact_path, artifact$path) ||
      !identical(population$artifact_sha256, artifact$sha256) ||
      !identical(artifact$export_operation_id, manifest$export_operation_id) ||
      !identical(artifact$gate_path, population$gate_path) ||
      !identical(artifact$acquisition_id, population$acquisition_id) ||
      !identical(artifact$sample_id, population$sample_id) ||
      !identical(artifact$population_key, population$population_key) ||
      !identical(artifact_channels, population_channels) ||
      !identical(artifact$identity_method_id, manifest$identity_method$id) ||
      !identical(artifact$intentionally_empty,
                 population$intentionally_empty) ||
      !identical(identity_columns,
                 c("acquisition_id", "event_index", "event_identity")) ||
      !identical(population$row_count, nrow(events))) {
    ph3_fail(population$acquisition_id, population_key, "artifact_ledger_mismatch",
             "schema, count, size, hash, or population linkage disagrees")
  }
  expected <- list(
    export_operation_id = manifest$export_operation_id,
    export_manifest_digest = manifest$manifest_binding$digest,
    acquisition_id = population$acquisition_id,
    sample_id = population$sample_id,
    export_profile = manifest$profile,
    identity_method_id = manifest$identity_method$id,
    identity_method_version = manifest$identity_method$version,
    identity_source = manifest$identity_method$source
  )
  for (field in names(expected)) {
    if (nrow(events) && any(is.na(events[[field]]) | !nzchar(events[[field]]) |
                           events[[field]] != expected[[field]])) {
      ph3_fail(population$acquisition_id, population_key, "row_binding_mismatch",
               paste0("field '", field, "' is blank, mixed, or inconsistent"))
    }
  }
  identities <- events$event_identity
  if (nrow(events) && (anyNA(identities) || any(!nzchar(identities)))) {
    ph3_fail(population$acquisition_id, population_key, "blank_event_identity",
             "direct event identities must be nonblank character data")
  }
  duplicate_rows <- sum(duplicated(identities) |
                        duplicated(identities, fromLast = TRUE))
  duplicate_bases <- length(unique(identities[duplicated(identities)]))
  if (duplicate_rows > 0L) {
    ph3_fail(population$acquisition_id, population_key,
             "duplicate_direct_identity",
             "direct identities must be unique within each population")
  }
  if (nrow(events) &&
      (any(is.na(events$event_index) |
           !grepl("^(0|[1-9][0-9]*)$", events$event_index)) ||
       any(events$event_identity != paste0(
         events$acquisition_id, ":event_index:", events$event_index
       )) ||
       anyNA(events$duplicate_occurrence) ||
       any(events$duplicate_occurrence != "1") ||
       anyNA(events$export_manifest_reference) ||
       any(events$export_manifest_reference !=
           "export-manifest.json + export-manifest.sha256"))) {
    ph3_fail(population$acquisition_id, population_key,
             "direct_identity_contract_mismatch",
             "event identity, occurrence, or manifest reference is inconsistent")
  }
  if (!identical(population$identity_field, "event_identity") ||
      !identical(population$identity_method_id, manifest$identity_method$id) ||
      !identical(population$unique_identity_count, length(unique(identities))) ||
      (!is.null(population$duplicate_base_combination_count) &&
       (!ph3_nonnegative_integer(population$duplicate_base_combination_count) ||
        !identical(population$duplicate_base_combination_count, duplicate_bases))) ||
      !identical(population$duplicate_row_count, duplicate_rows)) {
    ph3_fail(population$acquisition_id, population_key,
             "population_identity_ledger_mismatch",
             "identity field/method or unique/duplicate counts disagree")
  }
  count_rows <- count_report$sample_id == population$sample_id &
    count_report$gate_name == population$gate_name
  if (sum(count_rows) != 1L ||
      !ph3_nonnegative_integer(count_report$exported_count[count_rows]) ||
      !ph3_nonnegative_integer(count_report$flowjo_saved_count[count_rows]) ||
      !identical(count_report$exported_count[count_rows], nrow(events)) ||
      !identical(count_report$flowjo_saved_count[count_rows], nrow(events))) {
    ph3_fail(population$acquisition_id, population_key,
             "population_count_mismatch",
             "population artifact and count report do not agree exactly")
  }
  intentionally_empty <- isTRUE(artifact$intentionally_empty) &&
    isTRUE(population$intentionally_empty)
  if (nrow(events) > 0L && intentionally_empty) {
    ph3_fail(population$acquisition_id, population_key,
             "artifact_ledger_mismatch",
             "a nonempty population cannot be marked intentionally empty")
  }
  if (nrow(events) == 0L) {
    proof <- intentionally_empty && identical(artifact$row_count, 0L) &&
      identical(population$row_count, 0L) &&
      is.character(population$gate_path) && nzchar(population$gate_path) &&
      is.character(population$gate_type) && nzchar(population$gate_type) &&
      length(manifest$approval) > 0L &&
      identical(count_report$exported_count[count_rows], 0L) &&
      identical(count_report$flowjo_saved_count[count_rows], 0L)
    if (!proof) {
      ph3_fail(population$acquisition_id, population_key, "unproven_empty_population",
               "zero rows lack exact gate, artifact, count, schema, and approval proof")
    }
  }
  list(events = events, population = population, artifact = artifact, path = path,
       duplicate_rows = duplicate_rows, duplicate_bases = duplicate_bases,
       intentionally_empty = nrow(events) == 0L && intentionally_empty)
}

ph3_containment_row <- function(operation, prefix, parent, child) {
  acquisition <- child$population$acquisition_id
  if (!identical(parent$population$acquisition_id, acquisition) ||
      !identical(parent$population$sample_id, child$population$sample_id) ||
      !identical(child$population$parent_population_path,
                 parent$population$gate_path) ||
      !startsWith(child$population$gate_path,
                  paste0(parent$population$gate_path, "/"))) {
    ph3_fail(acquisition, child$population$population_key,
             "parent_child_acquisition_mismatch",
             "parent and child do not identify the same acquisition/sample")
  }
  parent_ids <- parent$events$event_identity
  child_ids <- child$events$event_identity
  unmatched <- child_ids[!child_ids %in% parent_ids]
  excess <- sum(pmax(table(factor(child_ids, levels = union(parent_ids, child_ids))) -
                     table(factor(parent_ids, levels = union(parent_ids, child_ids))), 0L))
  if (length(unmatched) || excess > 0L) {
    ph3_fail(acquisition, child$population$population_key,
             "child_not_contained_in_parent",
             paste0(length(unmatched), " unmatched identities; ", excess,
                    " excess occurrences"))
  }
  data.frame(
    containment_schema_version = "ph3-input-containment-1.0.0",
    acquisition_id = acquisition, sample_id = child$population$sample_id,
    prefix = prefix, parent_population_key = parent$population$population_key,
    parent_population_path = parent$population$gate_path,
    child_population_key = child$population$population_key,
    child_population_path = child$population$gate_path,
    export_profile = operation$manifest$profile,
    export_operation_id = operation$manifest$export_operation_id,
    manifest_digest = operation$full_digest,
    manifest_reference = operation$manifest_path,
    identity_method_id = operation$manifest$identity_method$id,
    identity_method_version = operation$manifest$identity_method$version,
    identity_source = operation$manifest$identity_method$source,
    parent_row_count = nrow(parent$events),
    parent_unique_identity_count = length(unique(parent$events$event_identity)),
    child_row_count = nrow(child$events),
    child_unique_identity_count = length(unique(child$events$event_identity)),
    parent_duplicated_identity_count = parent$duplicate_bases,
    parent_duplicate_row_count = parent$duplicate_rows,
    child_duplicated_identity_count = child$duplicate_bases,
    child_duplicate_row_count = child$duplicate_rows,
    matched_child_count = length(child_ids), unmatched_child_count = length(unmatched),
    excess_occurrence_count = excess,
    intentionally_empty = child$intentionally_empty,
    intentionally_empty_proof_status = if (child$intentionally_empty) "verified" else "not_applicable",
    containment_method_id = "exact_direct_identity_multiset_containment",
    containment_method_version = "1.0.0", containment_status = "validated",
    containment_reason_code = if (child$intentionally_empty) "validated_intentional_empty" else "contained",
    validation_severity = "none", parent_artifact_path = parent$path,
    parent_artifact_sha256 = parent$artifact$sha256,
    child_artifact_path = child$path,
    child_artifact_sha256 = child$artifact$sha256,
    stringsAsFactors = FALSE
  )
}

validate_ph3_export_operations <- function(config, data_dir = NULL) {
  manifest <- build_sample_manifest(config)
  operations <- lapply(ph3_resolve_operation_dirs(config, data_dir),
                       ph3_read_verified_operation)
  operation_ids <- vapply(operations, function(operation) {
    operation$manifest$export_operation_id
  }, character(1))
  if (anyDuplicated(operation_ids)) {
    ph3_fail("unknown", "manifest", "duplicate_export_operation_id",
             "distinct operation directories must have unique export_operation_id values")
  }
  count_reports <- lapply(operations, ph3_verify_count_report)
  configured_prefixes <- manifest$prefix
  operation_prefixes <- unlist(lapply(operations, function(operation) {
    vapply(operation$manifest$acquisitions, `[[`, character(1), "prefix")
  }), use.names = FALSE)
  if (anyDuplicated(operation_prefixes) ||
      !setequal(operation_prefixes, configured_prefixes)) {
    ph3_fail("unknown", "manifest", "acquisition_ledger_coverage_mismatch",
             "explicit operations must cover configured prefixes exactly once with no extras")
  }
  for (j in seq_along(operations)) {
    operation <- operations[[j]]
    acquisitions <- operation$manifest$acquisitions
    expected <- do.call(rbind, lapply(acquisitions, function(acquisition) {
      expand.grid(
        acquisition_id = acquisition$acquisition_id,
        sample_id = acquisition$sample_id,
        population_key = c("complete", "g1", "ph3_positive"),
        stringsAsFactors = FALSE
      )
    }))
    population_keys <- vapply(operation$manifest$populations, function(x) {
      paste(x$acquisition_id, x$sample_id, x$population_key, sep = "\r")
    }, character(1))
    artifact_entries <- Filter(function(x) identical(x$role, "population_events"),
                               operation$manifest$artifacts)
    artifact_keys <- vapply(artifact_entries, function(x) {
      paste(x$acquisition_id, x$sample_id, x$population_key, sep = "\r")
    }, character(1))
    expected_keys <- paste(expected$acquisition_id, expected$sample_id,
                           expected$population_key, sep = "\r")
    if (anyDuplicated(population_keys) || anyDuplicated(artifact_keys) ||
        !setequal(population_keys, expected_keys) ||
        !setequal(artifact_keys, expected_keys)) {
      ph3_fail("unknown", "manifest", "population_ledger_coverage_mismatch",
               "population and artifact ledgers must exactly cover every acquisition")
    }
    count_report <- count_reports[[j]]$report
    gate_lookup <- stats::setNames(
      vapply(operation$manifest$populations, `[[`, character(1), "gate_name"),
      population_keys
    )
    expected_count_keys <- paste(expected$sample_id,
      gate_lookup[expected_keys], sep = "\r")
    observed_count_keys <- paste(count_report$sample_id,
                                 count_report$gate_name, sep = "\r")
    if (anyDuplicated(observed_count_keys) ||
        !setequal(observed_count_keys, expected_count_keys)) {
      ph3_fail("unknown", "population_counts", "count_report_coverage_mismatch",
               "count report must exactly cover every acquisition/population")
    }
  }
  population_rows <- list()
  containment_rows <- list()
  verified_events <- list()
  for (i in seq_len(nrow(manifest))) {
    prefix <- manifest$prefix[[i]]
    hits <- vapply(operations, function(operation) {
      any(vapply(operation$manifest$acquisitions, function(x) identical(x$prefix, prefix), logical(1)))
    }, logical(1))
    if (sum(hits) != 1L) {
      ph3_fail(prefix, "manifest", "ambiguous_acquisition_mapping",
               "prefix must map to exactly one explicit export operation")
    }
    operation <- operations[[which(hits)]]
    count_report <- count_reports[[which(hits)]]$report
    acquisition_rows <- Filter(function(x) identical(x$prefix, prefix),
                               operation$manifest$acquisitions)
    if (length(acquisition_rows) != 1L) {
      ph3_fail(prefix, "manifest", "ambiguous_acquisition_mapping",
               "operation acquisition ledger must contain the prefix exactly once")
    }
    populations <- lapply(c("complete", "g1", "ph3_positive"), function(key) {
      ph3_verify_population_artifact(operation, prefix, key, count_report)
    })
    names(populations) <- c("complete", "g1", "ph3_positive")
    acquisition <- acquisition_rows[[1L]]
    if (any(!vapply(populations, function(x) {
      identical(x$population$acquisition_id, acquisition$acquisition_id) &&
        identical(x$population$sample_id, acquisition$sample_id) &&
        identical(x$population$export_operation_id,
                  operation$manifest$export_operation_id)
    }, logical(1)))) {
      ph3_fail(acquisition$acquisition_id, "manifest",
               "acquisition_operation_linkage_mismatch",
               "population ledgers disagree with acquisition or operation")
    }
    containment_rows[[length(containment_rows) + 1L]] <-
      ph3_containment_row(operation, prefix, populations$complete, populations$g1)
    containment_rows[[length(containment_rows) + 1L]] <-
      ph3_containment_row(operation, prefix, populations$complete, populations$ph3_positive)
    for (key in names(populations)) {
      item <- populations[[key]]
      minimum <- if (key == "complete") 10L else if (key == "g1") 2L else 0L
      missing_channels <- setdiff(c(config$dna_channel, config$target_channel),
                                  names(item$events))
      declared_channels <- unlist(item$population$channels, use.names = FALSE)
      gate_channels <- unlist(item$population$gate_channels, use.names = FALSE)
      measurement_columns <- unlist(lapply(declared_channels, function(channel) {
        c(paste0("raw__", channel), paste0("scaled__", channel))
      }), use.names = FALSE)
      invalid_declared_channels <- !length(declared_channels) ||
        anyNA(declared_channels) || any(!nzchar(declared_channels)) ||
        anyDuplicated(declared_channels) || !length(gate_channels) ||
        anyNA(gate_channels) || any(!nzchar(gate_channels)) ||
        anyDuplicated(gate_channels) || !all(gate_channels %in% declared_channels) ||
        !all(measurement_columns %in% names(item$events)) ||
        !all(c(config$dna_channel, config$target_channel) %in% measurement_columns)
      valid_channel_type <- function(value) {
        is.numeric(value) || (nrow(item$events) == 0L && is.logical(value))
      }
      if (invalid_declared_channels || length(missing_channels) ||
          any(!vapply(item$events[c(config$dna_channel, config$target_channel)],
                      valid_channel_type, logical(1)))) {
        ph3_fail(item$population$acquisition_id, key,
                 "configured_channel_schema_mismatch",
                 "configured DNA/target columns are missing or nonnumeric")
      }
      finite_n <- sum(is.finite(item$events[[config$dna_channel]]) &
                      is.finite(item$events[[config$target_channel]]))
      if (finite_n < minimum) {
        ph3_fail(item$population$acquisition_id, key,
                 "insufficient_finite_events",
                 paste0("found ", finite_n, ", require at least ", minimum))
      }
      population_rows[[length(population_rows) + 1L]] <- data.frame(
        replicate = manifest$replicate[[i]],
        replicate_index = manifest$replicate_index[[i]],
        condition = manifest$condition[[i]], prefix = prefix,
        population = key, required = TRUE, path = item$path, exists = TRUE,
        event_n = nrow(item$events),
        finite_event_n = finite_n,
        nonfinite_event_n = sum(!is.finite(item$events[[config$dna_channel]]) |
                                !is.finite(item$events[[config$target_channel]])),
        subset_membership_validated = key != "complete", stringsAsFactors = FALSE
      )
    }
    verified_events[[prefix]] <- lapply(populations, `[[`, "events")
  }
  report <- do.call(rbind, population_rows)
  rownames(report) <- NULL
  attr(report, "ph3_containment") <- do.call(rbind, containment_rows)
  attr(report, "ph3_export_manifests") <- lapply(operations, function(x) list(
    reference = x$manifest_path, detached_digest_reference = x$digest_path,
    sha256 = x$full_digest, export_operation_id = x$manifest$export_operation_id,
    manifest = x$manifest
  ))
  attr(report, "ph3_verified_events") <- verified_events
  report
}
