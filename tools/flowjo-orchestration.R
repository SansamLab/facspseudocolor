# Optional repository-level FlowJo orchestration.
#
# This file is intentionally outside R/. It may launch the unchanged Python
# exporter, but it is not installed with or required by the facspseudocolor R
# package. Users who already have event-level CSV files never source this file.

flowjo_or <- function(value, default) if (is.null(value)) default else value

flowjo_safe_name <- function(value) {
  value <- gsub("[^A-Za-z0-9._-]+", "_", flowjo_or(value, "replicate"))
  value <- sub("^_+", "", sub("_+$", "", value))
  if (nzchar(value)) value else "replicate"
}

flowjo_export_column <- function(column_names, channel, source = "raw") {
  exact <- paste0(source, "__", channel)
  if (exact %in% column_names) return(exact)
  hits <- column_names[startsWith(column_names, paste0(exact, " "))]
  if (length(hits) == 1L) return(hits[[1]])
  hits <- column_names[startsWith(column_names, exact)]
  if (length(hits) == 1L) return(hits[[1]])
  stop(
    "Could not find a unique '", source, "__", channel,
    "' column in the FlowJo export. Available: ",
    paste(column_names, collapse = ", "), call. = FALSE
  )
}

flowjo_sample_id <- function(sample_ids, fcs) {
  target <- tolower(sub("\\.fcs$", "", fcs, ignore.case = TRUE))
  normalized <- tolower(sub("\\.fcs$", "", sample_ids, ignore.case = TRUE))
  hits <- which(normalized == target | tolower(sample_ids) == tolower(fcs))
  if (!length(hits)) {
    stop(
      "No FlowJo sample matched FCS file '", fcs,
      "'. Exported sample ids: ", paste(sample_ids, collapse = ", "),
      call. = FALSE
    )
  }
  hits[[1]]
}

prepare_flowjo_csvs_external <- function(
    config,
    exporter = "python/export_flowjo_populations.py",
    verbose = TRUE
) {
  if (!inherits(config, "facs_config")) {
    config <- facspseudocolor::validate_facs_config(config)
  }
  replicates <- config$replicates
  if (is.null(replicates)) {
    stop("FlowJo orchestration currently requires a `replicates` configuration.",
         call. = FALSE)
  }
  data_dir <- config$data_dir
  if (!grepl("^/", data_dir)) {
    config_dir <- attr(config, "config_dir")
    if (is.null(config_dir)) {
      stop("Relative data paths require a file-backed configuration.",
           call. = FALSE)
    }
    data_dir <- file.path(config_dir, data_dir)
  }
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

  for (replicate in replicates) {
    flowjo <- utils::modifyList(
      flowjo_or(config$flowjo, list()), flowjo_or(replicate$flowjo, list())
    )
    if (!length(flowjo)) {
      stop("A top-level or per-replicate `flowjo` block is required.",
           call. = FALSE)
    }
    required <- c("source_dir", "workspace", "dna_source_channel",
                  "target_source_channel")
    missing <- required[vapply(flowjo[required], function(x) {
      is.null(x) || !is.character(x) || length(x) != 1L || !nzchar(x)
    }, logical(1))]
    if (length(missing)) {
      stop("Missing FlowJo setting(s) for ", replicate$label, ": ",
           paste(missing, collapse = ", "), call. = FALSE)
    }

    source_dir <- flowjo$source_dir
    workspace <- file.path(source_dir, flowjo$workspace)
    python <- flowjo_or(flowjo$python, Sys.which("python3"))
    rebuild <- isTRUE(flowjo_or(flowjo$rebuild, TRUE))
    population_map <- flowjo$populations
    if (is.null(population_map)) {
      population_map <- list(
        complete = flowjo_or(flowjo$population, "Single Cells")
      )
    }
    population_map <- population_map[
      intersect(names(population_map), names(config$suffixes))
    ]
    if (!length(population_map)) {
      stop("No configured FlowJo populations match the required file suffixes.",
           call. = FALSE)
    }
    if (!nzchar(python) || !file.exists(python)) {
      stop("Python interpreter not found: ", python, call. = FALSE)
    }
    if (!file.exists(exporter)) {
      stop("Python exporter not found: ", exporter, call. = FALSE)
    }
    if (!file.exists(workspace)) {
      stop("FlowJo workspace not found: ", workspace, call. = FALSE)
    }
    dependency_status <- suppressWarnings(system2(
      python, c("-c", shQuote("import flowkit, pandas, lxml")),
      stdout = FALSE, stderr = FALSE
    ))
    if (dependency_status != 0) {
      stop(
        "The selected Python environment is missing flowkit, pandas, or lxml: ",
        python, call. = FALSE
      )
    }

    prefixes <- vapply(replicate$samples, `[[`, character(1), "prefix")
    expected <- unlist(lapply(names(population_map), function(key) {
      file.path(data_dir, paste0(prefixes, config$suffixes[[key]]))
    }))
    if (!rebuild && all(file.exists(expected))) {
      if (verbose) message("FlowJo CSVs already exist for ", replicate$label)
      next
    }

    export_dir <- file.path(
      data_dir, ".flowjo_export", flowjo_safe_name(replicate$label)
    )
    dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)
    populations <- unique(unlist(population_map, use.names = FALSE))
    if (verbose) message("Running unchanged FlowJo exporter for ", replicate$label)
    status <- system2(python, c(
      shQuote(exporter), shQuote(workspace),
      "--fcs-dir", shQuote(source_dir),
      "--output-dir", shQuote(export_dir),
      "--populations", vapply(populations, shQuote, character(1))
    ))
    if (status != 0) {
      stop("FlowJo export failed for ", replicate$label, call. = FALSE)
    }

    for (key in names(population_map)) {
      population <- population_map[[key]]
      export_csv <- file.path(
        export_dir,
        paste0(gsub("[^A-Za-z0-9._-]+", "_", population), ".csv")
      )
      if (!file.exists(export_csv)) {
        stop("Expected FlowJo export was not created: ", export_csv,
             call. = FALSE)
      }
      combined <- utils::read.csv(export_csv, check.names = FALSE)
      if (!"sample_id" %in% names(combined) || !nrow(combined)) {
        stop("FlowJo export is empty or lacks sample_id: ", export_csv,
             call. = FALSE)
      }
      dna_column <- flowjo_export_column(
        names(combined), flowjo$dna_source_channel
      )
      target_column <- flowjo_export_column(
        names(combined), flowjo$target_source_channel
      )
      sample_ids <- unique(combined$sample_id)
      for (sample in replicate$samples) {
        index <- flowjo_sample_id(sample_ids, sample$fcs)
        selected <- combined$sample_id == sample_ids[[index]]
        output <- data.frame(
          combined[[dna_column]][selected],
          combined[[target_column]][selected]
        )
        names(output) <- c(config$dna_channel, config$target_channel)
        output_file <- file.path(
          data_dir, paste0(sample$prefix, config$suffixes[[key]])
        )
        utils::write.csv(output, output_file, row.names = FALSE)
        if (verbose) message("Wrote ", output_file)
      }
    }
  }
  invisible(TRUE)
}
