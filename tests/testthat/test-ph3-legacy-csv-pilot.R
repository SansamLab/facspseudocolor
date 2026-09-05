# Every constructed value in this file is explicitly SYNTHETIC and test-only.

synthetic_ph3_legacy_pilot_config <- function() {
  list(
    plot_type = "ph3", ph3_input_profile = "legacy_csv_pilot_v1",
    ph3_positivity_method = "ph3_raw_4n_density_cutoff_v1", data_dir = "SYNTHETIC",
    dna_channel = "SYNTHETIC DNA", target_channel = "SYNTHETIC pH3", target_name = "SYNTHETIC pH3",
    suffixes = list(complete = "_single.csv", g1 = "_g1.csv"), dna_2n_value = 1000,
    g1_x_range = c(900, 1150), s_phase_bins = list(early = c(1150, 1350), mid = c(1350, 1600), late = c(1600, 1850)),
    g2m_x_range = c(1850, 2300), ph3_boundary_sensitivity_fraction = .05,
    output_pdf = "SYNTHETIC.pdf", output_png = "SYNTHETIC.png", flowjo = list(rebuild = FALSE),
    ph3_pilot = list(control_label = "Untreated", provenance_label = "PILOT / LIMITED-PROVENANCE — NOT PUBLICATION-GRADE"),
    replicates = lapply(1:3, function(i) list(id = paste0("SYNTHETIC-set-", i), label = paste("SYNTHETIC", i), samples = list(
      list(label = "Untreated", prefix = paste0("SYNTHETIC-", i, "-control")),
      list(label = "2 h NMPP1", prefix = paste0("SYNTHETIC-", i, "-treatment")))))
  )
}

# Constructed values below are unmistakably SYNTHETIC.  Each control has a
# deliberately different raw-pH3 scale, so applying another clone's cutoff to
# its treatment would change the expected positivity result.
synthetic_ph3_legacy_pilot_data <- function(four_n_raw) {
  data.frame(
    # Distinct finite DNA coordinates keep the display-only SYNTHETIC KDE
    # well-defined while remaining entirely inside the synthetic 4N region.
    dna_norm = rep(seq(1900, 2100, length.out = 20L), length.out = length(four_n_raw)),
    target_raw = as.numeric(four_n_raw),
    stringsAsFactors = FALSE
  )
}

# This fixture has a deliberately unique, dominant negative peak and a
# separated positive peak in every SYNTHETIC sample.  The 160/40 split also
# keeps both membership groups on identical repeated DNA coordinates, so the
# approved background-fit coverage rule can be evaluated without extrapolation.
synthetic_ph3_legacy_pilot_signal <- function(negative_center, positive_center) {
  c(
    rep(negative_center + c(-8, -4, 0, 4, 8), c(10L, 30L, 80L, 30L, 10L)),
    rep(positive_center + c(-8, 0, 8), c(5L, 30L, 5L))
  )
}

synthetic_ph3_legacy_pilot_analysis <- function(cutoff_failure_clone = NULL) {
  config <- validate_facs_config(synthetic_ph3_legacy_pilot_config())
  manifest <- build_sample_manifest(config)
  control_raw <- list(
    synthetic_ph3_legacy_pilot_signal(115, 275),
    synthetic_ph3_legacy_pilot_signal(355, 815),
    synthetic_ph3_legacy_pilot_signal(915, 2200)
  )
  treatment_raw <- list(
    synthetic_ph3_legacy_pilot_signal(120, 250),
    synthetic_ph3_legacy_pilot_signal(360, 750),
    synthetic_ph3_legacy_pilot_signal(960, 2000)
  )
  normalized_data <- lapply(seq_len(nrow(manifest)), function(i) {
    clone <- manifest$replicate_index[[i]]
    raw <- if (identical(manifest$condition[[i]], "Untreated")) {
      control_raw[[clone]]
    } else {
      treatment_raw[[clone]]
    }
    if (identical(manifest$condition[[i]], "Untreated") &&
        identical(clone, cutoff_failure_clone)) {
      # Explicitly SYNTHETIC zero-range 4N control: it must fail closed.
      raw <- rep(100, length(raw))
    }
    list(data = synthetic_ph3_legacy_pilot_data(raw))
  })
  names(normalized_data) <- manifest$prefix
  structure(list(
    config = config,
    sample_manifest = manifest,
    normalized_data = normalized_data,
    warnings = character()
  ), class = "facs_analysis")
}

test_that("SYNTHETIC legacy pilot accepts only complete and G1 inputs", {
  config <- validate_facs_config(synthetic_ph3_legacy_pilot_config())
  expect_identical(required_population_keys(config, build_sample_manifest(config)), c("complete", "g1"))
  # This in-memory SYNTHETIC config has no config-file directory attribute.
  # Supplying an explicit directory exercises the same input inventory without
  # relaxing production handling of relative configured data directories.
  expect_identical(sort(unique(facs_input_files(config, data_dir = tempdir())$population)), c("complete", "g1"))
})

test_that("SYNTHETIC reordered normalized-data tables fail before pilot event access", {
  analysis <- synthetic_ph3_legacy_pilot_analysis()
  analysis$normalized_data <- analysis$normalized_data[c(2L, 1L, 3L, 4L, 5L, 6L)]
  expect_error(
    ph3_build_legacy_csv_pilot(analysis),
    "normalized_data_prefix_alignment"
  )
})

test_that("SYNTHETIC cutoff failure makes every matched-clone outcome unavailable and cannot enter summaries", {
  result <- ph3_build_legacy_csv_pilot(
    synthetic_ph3_legacy_pilot_analysis(cutoff_failure_clone = 1L)
  )
  pilot <- result$ph3_legacy_pilot
  failed_cutoff <- pilot$cutoff_records[pilot$cutoff_records$replicate_index == 1L, , drop = FALSE]
  failed_values <- pilot$biological_replicate_values[
    pilot$biological_replicate_values$replicate_index == 1L, , drop = FALSE
  ]

  expect_identical(failed_cutoff$cutoff_status, "unavailable_cutoff_failure")
  expect_identical(failed_cutoff$reason_code, "unavailable_sample_background_correction")
  expect_true(all(is.na(failed_values$value)))
  expect_true(all(failed_values$value_status == "unavailable"))
  expect_true(all(failed_values$reason_code == "unavailable_cutoff_failure"))
  expect_true(all(is.na(failed_values$cutoff_corrected_signal)))

  # A display target is all-or-nothing within the explicit clone group.  The
  # failed cutoff cannot create a partial visual target from its treatment.
  failed_offsets <- pilot$display_offset_qc[pilot$display_offset_qc$replicate_index == 1L, , drop = FALSE]
  expect_true(all(failed_offsets$shared_target_status == "unavailable"))
  expect_true(all(failed_offsets$shared_target_reason_code ==
                    "configured_group_member_display_input_unavailable"))
  expect_true(all(failed_offsets$display_offset_status == "unavailable"))
  report <- plot_ph3_legacy_pilot_report(result)
  expect_true(all(vapply(report$pseudocolor_panels[failed_offsets$prefix], function(panel) {
    !any(vapply(panel$layers, function(layer) inherits(layer$geom, "GeomPoint") ||
      inherits(layer$geom, "GeomHline"), logical(1)))
  }, logical(1))))

  # The failed clone contributes neither zero nor any other value to means.
  # Compare within each condition: the Untreated and treatment values are
  # distinct outcomes and must never be pooled merely for this assertion.
  for (summary_row in seq_len(nrow(pilot$condition_summary))) {
    reported <- pilot$condition_summary[summary_row, , drop = FALSE]
    surviving <- pilot$biological_replicate_values[
      pilot$biological_replicate_values$outcome_id == reported$outcome_id[[1L]] &
        pilot$biological_replicate_values$condition_index == reported$condition_index[[1L]] &
        pilot$biological_replicate_values$replicate_index != 1L,
      "value", drop = TRUE
    ]
    finite <- is.finite(surviving)
    expect_identical(reported$finite_biological_replicate_count[[1L]], sum(finite))
    if (any(finite)) {
      expect_equal(reported$mean_value[[1L]], mean(surviving[finite]))
    } else {
      expect_true(is.na(reported$mean_value[[1L]]))
      expect_identical(reported$summary_status[[1L]], "unavailable")
    }
  }
})

test_that("SYNTHETIC unexpected cutoff errors stop rather than becoming unavailable results", {
  testthat::local_mocked_bindings(
    ph3_raw_4n_density_cutoff = function(...) {
      stop("SYNTHETIC unexpected internal failure", call. = FALSE)
    },
    .package = "facspseudocolor"
  )
  expect_error(
    ph3_build_legacy_csv_pilot(synthetic_ph3_legacy_pilot_analysis()),
    "SYNTHETIC unexpected internal failure"
  )
})

test_that("SYNTHETIC legacy pilot requires its prominent provenance label and no rebuild", {
  bad_label <- synthetic_ph3_legacy_pilot_config()
  bad_label$ph3_pilot$provenance_label <- "SYNTHETIC hidden"
  expect_error(validate_facs_config(bad_label), "ph3_pilot")
  bad_rebuild <- synthetic_ph3_legacy_pilot_config()
  bad_rebuild$flowjo$rebuild <- TRUE
  expect_error(validate_facs_config(bad_rebuild), "rebuild")
  bad_method <- synthetic_ph3_legacy_pilot_config()
  bad_method$ph3_positivity_method <- "flowjo_legacy_v1"
  expect_error(validate_facs_config(bad_method), "raw_4n_density")
})

test_that("SYNTHETIC legacy pilot derives final corrected cutoff only from matched untreated clone", {
  result <- ph3_build_legacy_csv_pilot(synthetic_ph3_legacy_pilot_analysis())
  pilot <- result$ph3_legacy_pilot
  cutoffs <- pilot$cutoff_records
  values <- pilot$biological_replicate_values

  expect_identical(cutoffs$control_prefix, c(
    "SYNTHETIC-1-control", "SYNTHETIC-2-control", "SYNTHETIC-3-control"
  ))
  control_rows <- match(cutoffs$control_prefix, result$sample_manifest$prefix)
  expect_false(anyNA(control_rows))
  expect_identical(
    result$sample_manifest$replicate_index[control_rows],
    cutoffs$replicate_index
  )
  expect_true(all(result$sample_manifest$condition[control_rows] == "Untreated"))
  expect_true(all(cutoffs$cutoff_status == "available"))
  expect_true(all(is.finite(cutoffs$cutoff_corrected_signal)))
  expect_true(all(pilot$provisional_background_qc$correction_status == "available"))

  for (clone in seq_len(3L)) {
    clone_values <- values[values$replicate_index == clone, , drop = FALSE]
    expected_cutoff <- cutoffs$cutoff_corrected_signal[[clone]]
    expect_identical(
      unique(clone_values$cutoff_corrected_signal), expected_cutoff,
      info = paste("SYNTHETIC clone", clone, "must use its own untreated corrected cutoff")
    )
  }

  # Background correction has no required cross-clone cutoff ordering. The
  # direct prefix, clone, condition, and value propagation checks above are
  # the auditable matched-control contract.
})

test_that("SYNTHETIC legacy pilot exposes one corrected pH3-versus-DNA panel per configured sample", {
  analysis <- ph3_build_legacy_csv_pilot(synthetic_ph3_legacy_pilot_analysis())
  result <- plot_ph3_legacy_pilot_report(analysis)
  metadata <- result$pseudocolor_metadata
  manifest <- analysis$sample_manifest

  expect_identical(metadata$prefix, manifest$prefix)
  expect_identical(names(result$pseudocolor_panels), manifest$prefix)
  expect_identical(nrow(metadata), 6L)
  expect_true(all(metadata$signal_basis == "background_subtracted"))
  expect_true(all(metadata$provenance_label ==
                    "PILOT / LIMITED-PROVENANCE — NOT PUBLICATION-GRADE"))
  expect_identical(metadata$cutoff_control_prefix,
                   analysis$ph3_legacy_pilot$cutoff_records$control_prefix[
                     match(metadata$replicate_index,
                           analysis$ph3_legacy_pilot$cutoff_records$replicate_index)
                   ])
  expect_identical(metadata$cutoff_corrected_signal,
                   analysis$ph3_legacy_pilot$cutoff_records$cutoff_corrected_signal[
                     match(metadata$replicate_index,
                           analysis$ph3_legacy_pilot$cutoff_records$replicate_index)
                   ])
  offset_qc <- result$display_offset_qc
  expect_identical(offset_qc$prefix, manifest$prefix)
  expect_true(all(offset_qc$display_offset_status == "available"))
  expect_true(all(is.finite(offset_qc$display_offset)))
  expect_identical(metadata$display_offset, offset_qc$display_offset)
  expect_identical(metadata$displayed_cutoff_signal, offset_qc$displayed_cutoff_signal)
  log_qc <- result$log_display_qc
  expect_identical(log_qc$prefix, manifest$prefix)
  expect_true(all(log_qc$log_display_status == "available"))
  expect_true(all(log_qc$positive_domain_display_event_count > 0L))
  expect_true(all(log_qc$nonpositive_display_event_count == 0L))
  expect_identical(metadata$log_display_status, log_qc$log_display_status)
  expect_true(all(vapply(result$pseudocolor_panels, inherits, logical(1), "ggplot")))
  # The visual window is 1.6N--4.4N, while the three analytic-region lines
  # remain independent display guides rather than plotting limits.
  expect_true(all(vapply(result$pseudocolor_panels, function(panel) {
    identical(panel$coordinates$limits$x, c(800, 2200))
  }, logical(1))))
  expect_true(all(vapply(result$pseudocolor_panels, function(panel) {
    any(vapply(panel$layers, function(layer) inherits(layer$geom, "GeomHline"), logical(1)))
  }, logical(1))))
  expect_true(all(vapply(result$pseudocolor_panels, function(panel) {
    grepl("PILOT / LIMITED-PROVENANCE", panel$labels$caption, fixed = TRUE) &&
      grepl("Background-subtracted", panel$labels$subtitle, fixed = TRUE) &&
      grepl("matched Untreated", panel$labels$subtitle, fixed = TRUE)
  }, logical(1))))
  expect_true(all(vapply(result$pseudocolor_panels, function(panel) {
    grepl("display offset; log10 display", panel$labels$subtitle, fixed = TRUE) &&
      grepl("display offset restored; log10 display", panel$labels$y, fixed = TRUE)
  }, logical(1))))
  expect_true(all(vapply(result$pseudocolor_panels, function(panel) {
    identical(panel$scales$get_scales("y")$trans$name, "log-10")
  }, logical(1))))
  # The KDE must use the log10 displayed-signal coordinate used by the panel,
  # rather than linear fluorescence followed by a log-scale transformation.
  # This is a SYNTHETIC display-only assertion; event values and analysis
  # outputs are checked separately below.
  expect_true(all(vapply(result$pseudocolor_panels, function(panel) {
    panel_data <- panel$data
    expected_density <- compute_point_density(
      panel_data$dna_norm, log10(panel_data$displayed_signal), n = 300L,
      bandwidth_multiplier = 0.5
    )
    isTRUE(all.equal(panel_data$density, expected_density, tolerance = 1e-12))
  }, logical(1))))
  # Log10 cannot transform the former +/-Inf rectangle extent.  The analytic
  # region shading must instead use the finite, positive panel limits derived
  # from the already filtered display events.
  expect_true(all(vapply(result$pseudocolor_panels, function(panel) {
    rect_layers <- Filter(function(layer) inherits(layer$geom, "GeomRect"), panel$layers)
    length(rect_layers) == 1L &&
      all(is.finite(rect_layers[[1L]]$data$ymin)) &&
      all(is.finite(rect_layers[[1L]]$data$ymax)) &&
      all(rect_layers[[1L]]$data$ymin > 0) &&
      all(rect_layers[[1L]]$data$ymax > rect_layers[[1L]]$data$ymin)
  }, logical(1))))
  # The report's physical figure dimensions work with this square data
  # viewport. This protects the display contract without changing the
  # displayed analytic regions, events, or cutoff calculations.
  expect_true(all(vapply(result$pseudocolor_panels, function(panel) {
    identical(panel$theme$aspect.ratio, 1)
  }, logical(1))))
  expect_true(all(vapply(result$pseudocolor_panels, function(panel) {
    identical(panel$theme$legend.title$size, 6) &&
      identical(panel$theme$legend.text$size, 5)
  }, logical(1))))
})

test_that("SYNTHETIC shared clone display target equally weights sample medians without changing analysis values", {
  analysis <- ph3_build_legacy_csv_pilot(synthetic_ph3_legacy_pilot_analysis())
  before_values <- analysis$ph3_legacy_pilot$biological_replicate_values
  before_cutoffs <- analysis$ph3_legacy_pilot$cutoff_records
  before_events <- analysis$ph3_legacy_pilot$event_corrections
  result <- plot_ph3_legacy_pilot_report(analysis)
  qc <- result$display_offset_qc

  for (i in seq_len(nrow(qc))) {
    prefix <- qc$prefix[[i]]
    events <- before_events[[prefix]]
    final_negative <- ph3_pilot_regions(analysis$config, events)$eligible &
      events$corrected_signal <= qc$cutoff_corrected_signal[[i]]
    expect_identical(qc$eligible_final_negative_event_count[[i]], sum(final_negative))
    expect_equal(qc$raw_negative_median[[i]], stats::median(events$target_raw[final_negative]))
    expect_equal(qc$corrected_negative_median[[i]], stats::median(events$corrected_signal[final_negative]))
    group <- qc$replicate_index == qc$replicate_index[[i]]
    expected_shared_target <- stats::median(qc$raw_negative_median[group])
    expect_identical(qc$shared_target_status[[i]], "available")
    expect_equal(qc$shared_raw_negative_median[[i]], expected_shared_target)
    expect_equal(qc$display_offset[[i]], expected_shared_target - qc$corrected_negative_median[[i]])
    expect_equal(qc$corrected_negative_median[[i]] + qc$display_offset[[i]], expected_shared_target)
    expect_equal(qc$displayed_cutoff_signal[[i]], qc$cutoff_corrected_signal[[i]] + qc$display_offset[[i]])
    panel <- result$pseudocolor_panels[[prefix]]
    hlines <- Filter(function(layer) inherits(layer$geom, "GeomHline"), panel$layers)
    expect_identical(length(hlines), 1L)
    # geom_hline() stores an explicitly supplied intercept in the layer's
    # one-row data, rather than in aes_params.  This verifies the final
    # display-scale cutoff without inspecting a presentation-only aesthetic.
    expect_identical(nrow(hlines[[1L]]$data), 1L)
    expect_equal(hlines[[1L]]$data$yintercept[[1L]], qc$displayed_cutoff_signal[[i]])
  }

  expect_identical(analysis$ph3_legacy_pilot$biological_replicate_values, before_values)
  expect_identical(analysis$ph3_legacy_pilot$cutoff_records, before_cutoffs)
  expect_identical(analysis$ph3_legacy_pilot$event_corrections, before_events)
})

test_that("SYNTHETIC shared clone display target fails every panel in an incomplete group", {
  analysis <- ph3_build_legacy_csv_pilot(synthetic_ph3_legacy_pilot_analysis())
  before_values <- analysis$ph3_legacy_pilot$biological_replicate_values
  before_cutoffs <- analysis$ph3_legacy_pilot$cutoff_records
  corrections <- analysis$ph3_legacy_pilot$event_corrections
  failed_prefix <- analysis$sample_manifest$prefix[[1L]]
  corrections[[failed_prefix]]$corrected_signal[[1L]] <- NA_real_
  qc <- ph3_pilot_display_offset_qc(analysis, analysis$ph3_legacy_pilot$cutoff_records, corrections)

  failed_group <- qc$replicate_index == 1L
  intact_group <- qc$replicate_index == 2L
  expect_true(all(qc$shared_target_status[failed_group] == "unavailable"))
  expect_true(all(qc$shared_target_reason_code[failed_group] ==
                    "configured_group_member_display_input_unavailable"))
  expect_true(all(qc$display_offset_status[failed_group] == "unavailable"))
  expect_true(all(qc$display_offset_reason_code[failed_group] ==
                    "clone_group_shared_raw_negative_median_unavailable"))
  expect_true(all(is.na(qc$display_offset[failed_group])))
  expect_true(all(qc$shared_target_status[intact_group] == "available"))
  expect_identical(analysis$ph3_legacy_pilot$biological_replicate_values, before_values)
  expect_identical(analysis$ph3_legacy_pilot$cutoff_records, before_cutoffs)
})

test_that("SYNTHETIC log10 display QC excludes only nonpositive restored values", {
  analysis <- ph3_build_legacy_csv_pilot(synthetic_ph3_legacy_pilot_analysis())
  before_values <- analysis$ph3_legacy_pilot$biological_replicate_values
  prefix <- analysis$sample_manifest$prefix[[1L]]
  corrections <- analysis$ph3_legacy_pilot$event_corrections
  offset <- analysis$ph3_legacy_pilot$display_offset_qc[
    analysis$ph3_legacy_pilot$display_offset_qc$prefix == prefix, "display_offset", drop = TRUE
  ]
  # This is an explicit SYNTHETIC display-only mutation: it makes a subset of
  # visual-window restored values nonpositive without altering the analysis
  # object, its corrected cutoff, or any stored outcome.
  corrections[[prefix]]$corrected_signal[seq_len(10L)] <- -offset - 1
  qc <- ph3_pilot_log_display_qc(
    analysis, analysis$ph3_legacy_pilot$display_offset_qc, corrections
  )
  row <- qc[qc$prefix == prefix, , drop = FALSE]
  expect_identical(row$log_display_status, "available_nonpositive_display_values_excluded")
  expect_identical(row$log_display_reason_code, "nonpositive_display_values_excluded")
  expect_true(row$positive_domain_display_event_count > 0L)
  expect_true(row$nonpositive_display_event_count > 0L)
  expect_identical(analysis$ph3_legacy_pilot$biological_replicate_values, before_values)

  corrections[[prefix]]$corrected_signal <- -offset - 1
  zero_qc <- ph3_pilot_log_display_qc(
    analysis, analysis$ph3_legacy_pilot$display_offset_qc, corrections
  )
  zero_row <- zero_qc[zero_qc$prefix == prefix, , drop = FALSE]
  expect_identical(zero_row$log_display_status, "unavailable")
  expect_identical(zero_row$log_display_reason_code, "no_positive_domain_display_events")
  expect_identical(zero_row$positive_domain_display_event_count, 0L)
})

test_that("SYNTHETIC zero-positive log display is explicit and has no cutoff overlay", {
  analysis <- ph3_build_legacy_csv_pilot(synthetic_ph3_legacy_pilot_analysis())
  prefix <- analysis$sample_manifest$prefix[[1L]]
  i <- match(prefix, analysis$ph3_legacy_pilot$log_display_qc$prefix)
  analysis$ph3_legacy_pilot$log_display_qc$log_display_status[[i]] <- "unavailable"
  analysis$ph3_legacy_pilot$log_display_qc$log_display_reason_code[[i]] <-
    "no_positive_domain_display_events"
  analysis$ph3_legacy_pilot$log_display_qc$positive_domain_display_event_count[[i]] <- 0L
  analysis$ph3_legacy_pilot$log_display_qc$nonpositive_display_event_count[[i]] <- 200L
  metadata_i <- match(prefix, analysis$ph3_legacy_pilot$pseudocolor_metadata$prefix)
  analysis$ph3_legacy_pilot$pseudocolor_metadata$log_display_status[[metadata_i]] <- "unavailable"
  analysis$ph3_legacy_pilot$pseudocolor_metadata$log_display_reason_code[[metadata_i]] <-
    "no_positive_domain_display_events"
  analysis$ph3_legacy_pilot$pseudocolor_metadata$positive_domain_display_event_count[[metadata_i]] <- 0L
  analysis$ph3_legacy_pilot$pseudocolor_metadata$nonpositive_display_event_count[[metadata_i]] <- 200L

  result <- plot_ph3_legacy_pilot_report(analysis)
  panel <- result$pseudocolor_panels[[prefix]]
  expect_match(panel$labels$subtitle, "LOG10 DISPLAY", ignore.case = TRUE)
  expect_false(any(vapply(panel$layers, function(layer) inherits(layer$geom, "GeomPoint") ||
    inherits(layer$geom, "GeomHline"), logical(1))))
})

test_that("SYNTHETIC unavailable display offset renders no substituted signal", {
  analysis <- ph3_build_legacy_csv_pilot(synthetic_ph3_legacy_pilot_analysis())
  prefix <- analysis$sample_manifest$prefix[[1L]]
  i <- match(prefix, analysis$ph3_legacy_pilot$display_offset_qc$prefix)
  analysis$ph3_legacy_pilot$display_offset_qc$display_offset_status[[i]] <- "unavailable"
  analysis$ph3_legacy_pilot$display_offset_qc$display_offset_reason_code[[i]] <- "no_eligible_final_negative_events"
  analysis$ph3_legacy_pilot$display_offset_qc$display_offset[[i]] <- NA_real_
  analysis$ph3_legacy_pilot$display_offset_qc$displayed_cutoff_signal[[i]] <- NA_real_
  analysis$ph3_legacy_pilot$pseudocolor_metadata$display_offset_status[[i]] <- "unavailable"
  analysis$ph3_legacy_pilot$pseudocolor_metadata$display_offset_reason_code[[i]] <- "no_eligible_final_negative_events"
  analysis$ph3_legacy_pilot$pseudocolor_metadata$display_offset[[i]] <- NA_real_
  analysis$ph3_legacy_pilot$pseudocolor_metadata$displayed_cutoff_signal[[i]] <- NA_real_

  result <- plot_ph3_legacy_pilot_report(analysis)
  panel <- result$pseudocolor_panels[[prefix]]
  expect_match(panel$labels$subtitle, "Corrected display offset unavailable", fixed = TRUE)
  expect_match(panel$labels$caption, "No corrected event display is fabricated", fixed = TRUE)
  expect_false(any(vapply(panel$layers, function(layer) inherits(layer$geom, "GeomPoint") ||
    inherits(layer$geom, "GeomHline"), logical(1))))
})

test_that("SYNTHETIC unavailable cutoff panel is explicit and has no cutoff overlay", {
  analysis <- ph3_build_legacy_csv_pilot(synthetic_ph3_legacy_pilot_analysis())
  # Constructed display metadata is changed together with its SYNTHETIC
  # authoritative cutoff row, so this checks the unavailable-display branch
  # without representing any experimental observation.
  analysis$ph3_legacy_pilot$cutoff_records$cutoff_corrected_signal[[1L]] <- NA_real_
  analysis$ph3_legacy_pilot$cutoff_records$cutoff_status[[1L]] <- "unavailable_cutoff_failure"
  analysis$ph3_legacy_pilot$cutoff_records$reason_code[[1L]] <- "zero_or_invalid_control_signal_range"
  first_clone <- analysis$ph3_legacy_pilot$pseudocolor_metadata$replicate_index == 1L
  analysis$ph3_legacy_pilot$pseudocolor_metadata$cutoff_corrected_signal[first_clone] <- NA_real_
  analysis$ph3_legacy_pilot$pseudocolor_metadata$cutoff_status[first_clone] <- "unavailable_cutoff_failure"
  analysis$ph3_legacy_pilot$pseudocolor_metadata$cutoff_reason_code[first_clone] <- "zero_or_invalid_control_signal_range"

  result <- plot_ph3_legacy_pilot_report(analysis)
  panel <- result$pseudocolor_panels[["SYNTHETIC-1-control"]]
  expect_match(panel$labels$subtitle, "Corrected display offset unavailable", fixed = TRUE)
  expect_false(any(vapply(panel$layers, function(layer) inherits(layer$geom, "GeomHline"), logical(1))))
})

test_that("SYNTHETIC provisional raw cutoff is audit-only and correction changes final scale", {
  analysis <- ph3_build_legacy_csv_pilot(synthetic_ph3_legacy_pilot_analysis())
  pilot <- analysis$ph3_legacy_pilot
  qc <- pilot$provisional_background_qc
  final <- pilot$cutoff_records
  expect_true(all(qc$provisional_cutoff_status == "available"))
  expect_true(all(is.finite(qc$provisional_cutoff_raw_signal)))
  expect_true(all(final$cutoff_basis == "background_subtracted_4n_pH3"))
  expect_false(identical(final$cutoff_corrected_signal,
                         qc$provisional_cutoff_raw_signal[match(final$control_prefix, qc$prefix)]))
  corrected <- pilot$event_corrections[["SYNTHETIC-1-control"]]
  expect_true(all(is.finite(corrected$predicted_background)))
  expect_true(any(corrected$corrected_signal != corrected$target_raw))
  expect_true(all(c("provisional_negative_member", "provisional_positive_member") %in% names(corrected)))
})

test_that("SYNTHETIC invalid background fit fails an entire clone closed", {
  model <- synthetic_ph3_legacy_pilot_analysis()
  # This unmistakably SYNTHETIC mutation removes DNA variation only from the
  # matched clone's provisional-negative fit; no experimental data are used.
  control <- model$sample_manifest$prefix[[1L]]
  model$normalized_data[[control]]$data$dna_norm <- 2000
  result <- ph3_build_legacy_csv_pilot(model)
  pilot <- result$ph3_legacy_pilot
  failed <- pilot$biological_replicate_values[pilot$biological_replicate_values$replicate_index == 1L, , drop = FALSE]
  expect_true(all(is.na(failed$value)))
  expect_true(all(failed$reason_code == "unavailable_cutoff_failure"))
  expect_true(all(pilot$condition_summary$finite_biological_replicate_count <= 2L))
})

test_that("SYNTHETIC unavailable correction renders a placeholder without numerical leakage", {
  model <- synthetic_ph3_legacy_pilot_analysis()
  # This unmistakably SYNTHETIC mutation invalidates one sample-level fit only.
  # It must never result in raw or invented corrected events in that panel.
  failed_prefix <- model$sample_manifest$prefix[[1L]]
  model$normalized_data[[failed_prefix]]$data$dna_norm <- 2000
  analysis <- ph3_build_legacy_csv_pilot(model)
  result <- plot_ph3_legacy_pilot_report(analysis)
  failed_panel <- result$pseudocolor_panels[[failed_prefix]]
  valid_prefix <- model$sample_manifest$prefix[[3L]]
  valid_panel <- result$pseudocolor_panels[[valid_prefix]]
  failed_qc <- result$provisional_background_qc[
    result$provisional_background_qc$prefix == failed_prefix, , drop = FALSE
  ]
  expect_identical(failed_qc$correction_status, "unavailable")
  expect_match(failed_panel$labels$subtitle, "Background correction unavailable", fixed = TRUE)
  expect_match(failed_panel$labels$caption, "No corrected event display is fabricated", fixed = TRUE)
  expect_false(any(vapply(failed_panel$layers, function(layer) inherits(layer$geom, "GeomPoint") ||
    inherits(layer$geom, "GeomHline"), logical(1))))
  expect_true(any(vapply(failed_panel$layers, function(layer) inherits(layer$geom, "GeomText"), logical(1))))
  expect_true(any(vapply(valid_panel$layers, function(layer) inherits(layer$geom, "GeomPoint"), logical(1))))
  failed_values <- result$biological_replicate_values[
    result$biological_replicate_values$replicate_index == 1L, , drop = FALSE
  ]
  expect_true(all(is.na(failed_values$value)))
  expect_true(all(result$condition_summary$finite_biological_replicate_count <= 2L))
})
