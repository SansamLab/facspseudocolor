# =============================================================================
# pseudocolor_helpers.R
#
# Reusable helpers for DNA-content-normalized pseudocolor plots of a single
# channel (an EdU incorporation signal or a protein-of-interest signal) versus
# DNA content, laid out as a replicate (row) x condition (column) panel grid.
#
# Two processing modes, chosen with `plot_type`:
#
#   "edu"  Reference EdU-negative regression baseline normalization
#          (reproduces the main Figure 1 pipeline). For each replicate, the
#          designated reference sample is used to (1) fit a positive-minimum
#          boundary that separates EdU+ from EdU- events, (2) linearly regress
#          the EdU-negative signal on normalized DNA content to get a slope,
#          and (3) build a per-sample baseline line that shares that slope but
#          is forced through the sample's own G1 anchor point. Every sample's
#          signal is then divided by this baseline (and scaled to dna_2n_value),
#          so EdU-negative cells sit at the baseline and EdU+ cells rise above.
#
#   "poi"  Protein-of-interest normalization. DNA content is normalized the
#          same way, but the target channel is either divided by its own G1
#          median (normalize_target = TRUE) or plotted raw
#          (normalize_target = FALSE). No reference sample is required.
#
# Packages are called with explicit namespaces, so sourcing this file does not
# attach packages or alter the caller's search path. Required packages:
# ggplot2, MASS, scales, cowplot, grid.
# =============================================================================


# ---------------------------------------------------------------------------
# Parameter helpers
# ---------------------------------------------------------------------------

# Quarto wraps complex YAML params in a list(value = ...). This unwraps them so
# the same code works whether a value came from Quarto params or a plain list.
unwrap_param <- function(value) {
  if (is.list(value) && identical(names(value), "value")) {
    value$value
  } else {
    value
  }
}

# Read a named setting with a fallback default.
get_setting <- function(settings, name, default = NULL) {
  value <- settings[[name]]
  if (is.null(value)) default else value
}


# ---------------------------------------------------------------------------
# Color palettes
# ---------------------------------------------------------------------------

# Project-standard "refined" density palette used for Figure 1
# (light blue -> blue -> deep blue -> cyan -> yellow -> dark red).
refined_density_palette <- function() {
  c("#D7EBFF", "#69B7FF", "#2255D6", "#12BED0", "#F2D13D", "#7F0000")
}

# FlowJo-style rainbow density palette (low -> high density).
flowjo_palette <- function() {
  c(
    "#503AFF", "#3480FF", "#27D7FC", "#38FFA1", "#74FF4D", "#B2FF33",
    "#F8FF33", "#FFEA4D", "#FFBF4D", "#FF844D", "#FF4A51"
  )
}

# Resolve a palette from a name ("refined" / "flowjo") or a custom hex vector.
resolve_palette <- function(palette) {
  if (is.character(palette) && length(palette) == 1) {
    switch(
      palette,
      refined = refined_density_palette(),
      flowjo  = flowjo_palette(),
      stop(paste0("Unknown palette name: ", palette))
    )
  } else if (is.character(palette) && length(palette) >= 2) {
    palette
  } else {
    stop("palette must be 'refined', 'flowjo', or a vector of >= 2 hex colors.")
  }
}


# ---------------------------------------------------------------------------
# Density estimation and color scaling
# ---------------------------------------------------------------------------

# 2D kernel density estimate evaluated at each point (x, y). Returns one density
# value per input row (NA for non-finite inputs). Faithful to MASS::kde2d.
compute_point_density <- function(x, y, n = 300, bandwidth_multiplier = 0.5) {
  valid <- is.finite(x) & is.finite(y)
  density_values <- rep(NA_real_, length(x))

  x_valid <- x[valid]
  y_valid <- y[valid]

  if (length(x_valid) < 10) {
    stop("Too few valid events for KDE density calculation.")
  }

  default_h <- c(MASS::bandwidth.nrd(x_valid), MASS::bandwidth.nrd(y_valid))

  if (any(!is.finite(default_h)) || any(default_h <= 0)) {
    stop("Could not calculate a valid KDE bandwidth.")
  }

  dens <- MASS::kde2d(
    x = x_valid, y = y_valid, n = n, h = bandwidth_multiplier * default_h
  )

  ix <- findInterval(x_valid, dens$x, all.inside = TRUE)
  iy <- findInterval(y_valid, dens$y, all.inside = TRUE)

  density_values[valid] <- dens$z[cbind(ix, iy)]
  density_values
}

# Map raw density onto a clipped, gamma-adjusted 0-1 color scale so the palette
# spreads across the informative density range rather than a few hot pixels.
prepare_density_color <- function(
    density,
    lower_clip_quantile = 0.05,
    upper_clip_quantile = 0.90,
    gamma = 0.35
) {
  valid_density <- density[is.finite(density)]

  if (length(valid_density) < 2) {
    stop("Too few finite density values for color scaling.")
  }

  if (
    lower_clip_quantile < 0 || upper_clip_quantile > 1 ||
    lower_clip_quantile >= upper_clip_quantile
  ) {
    stop("Density clipping quantiles must satisfy 0 <= lower < upper <= 1.")
  }

  if (!is.finite(gamma) || gamma <= 0) {
    stop("Density gamma must be a positive finite number.")
  }

  lower_limit <- stats::quantile(
    valid_density, probs = lower_clip_quantile, na.rm = TRUE, names = FALSE
  )
  upper_limit <- stats::quantile(
    valid_density, probs = upper_clip_quantile, na.rm = TRUE, names = FALSE
  )

  if (!is.finite(lower_limit) || !is.finite(upper_limit) ||
      upper_limit <= lower_limit) {
    stop("Invalid density clipping limits.")
  }

  density_clipped <- pmin(pmax(density, lower_limit), upper_limit)
  density_scaled <- scales::rescale(
    density_clipped, from = c(lower_limit, upper_limit), to = c(0, 1)
  )
  pmin(pmax(density_scaled^gamma, 0), 1)
}


# ---------------------------------------------------------------------------
# Sample manifest (replicates x conditions)
# ---------------------------------------------------------------------------

# Build a tidy manifest from either a flat `samples` list (single row) or a
# nested `replicates` list (one row per biological replicate). Each sample needs
# a `label` (condition name) and a `prefix` (file-name stem). For "edu" mode a
# `reference` condition should be named on each replicate (matches one sample).
make_sample_manifest <- function(
    samples = NULL,
    replicates = NULL,
    default_replicate_label = ""
) {
  samples <- unwrap_param(samples)
  replicates <- unwrap_param(replicates)

  if (!is.null(samples) && !is.null(replicates)) {
    stop("Specify either samples or replicates, not both.")
  }
  if (is.null(samples) && is.null(replicates)) {
    stop("No samples or biological replicates were supplied.")
  }

  validate_samples <- function(
      sample_list, replicate_label, replicate_index,
      reference_condition = NA_character_) {
    sample_list <- unwrap_param(sample_list)

    if (!is.list(sample_list) || length(sample_list) == 0) {
      stop(paste("No samples were supplied for", replicate_label))
    }

    rows <- lapply(seq_along(sample_list), function(sample_index) {
      sample <- sample_list[[sample_index]]

      if (is.null(sample$label) || is.null(sample$prefix) ||
          !nzchar(sample$label) || !nzchar(sample$prefix)) {
        stop(paste(
          "Every sample must have a nonempty label and prefix in",
          replicate_label
        ))
      }

      data.frame(
        replicate = replicate_label,
        replicate_index = replicate_index,
        reference_condition = reference_condition,
        is_reference = !is.na(reference_condition) &&
          sample$label == reference_condition,
        condition = sample$label,
        condition_index = sample_index,
        prefix = sample$prefix,
        stringsAsFactors = FALSE
      )
    })

    do.call(rbind, rows)
  }

  if (!is.null(samples)) {
    return(validate_samples(samples, default_replicate_label, 1))
  }

  replicate_rows <- lapply(seq_along(replicates), function(replicate_index) {
    replicate <- replicates[[replicate_index]]
    replicate_label <- replicate$label

    if (is.null(replicate_label) || !nzchar(replicate_label)) {
      replicate_label <- paste("Replicate", replicate_index)
    }

    reference_condition <- replicate$reference
    if (is.null(reference_condition) || !nzchar(reference_condition)) {
      reference_condition <- NA_character_
    }

    validate_samples(
      replicate$samples, replicate_label, replicate_index, reference_condition
    )
  })

  manifest <- do.call(rbind, replicate_rows)
  rownames(manifest) <- NULL

  condition_counts <- table(manifest$replicate_index)
  if (length(unique(condition_counts)) != 1) {
    stop("Every biological replicate must contain the same number of conditions.")
  }

  conditions_by_replicate <- split(manifest$condition, manifest$replicate_index)
  first_conditions <- conditions_by_replicate[[1]]
  matching <- vapply(
    conditions_by_replicate, identical, logical(1), y = first_conditions
  )
  if (!all(matching)) {
    stop("Replicates must contain the same conditions in the same order.")
  }

  manifest
}


# ---------------------------------------------------------------------------
# G1 anchor for the target channel
# ---------------------------------------------------------------------------

# The value the baseline line is forced through at the 2N (G1) position.
# "median" = median of the G1-gated target signal (Figure 1 default).
# "mode"   = peak of the kernel density of the G1-gated target signal.
g1_target_anchor <- function(values, method = "median") {
  values <- values[is.finite(values)]
  if (length(values) < 2) stop("Too few finite G1 target values.")

  if (method == "median") {
    stats::median(values)
  } else if (method == "mode") {
    dens <- stats::density(values, n = 512)
    dens$x[which.max(dens$y)]
  } else {
    stop(paste0("Unknown g1_anchor method: ", method))
  }
}


# ---------------------------------------------------------------------------
# Read one sample and normalize
# ---------------------------------------------------------------------------

# Reads a sample's "all cells" table (and, for the G1-based methods, its
# G1-gated table) and normalizes:
#   * DNA content -> 2N at dna_2n_value, 4N at 2x. The 2N anchor is the sample's
#     G1-gated DNA median (g1 methods) or the replicate background sample's DNA
#     density peak carried in `background_model` (poi method).
#   * target signal, depending on `normalization_method`:
#       "g1_median"  -> raw / g1_target_anchor * dna_2n_value   (generic poi)
#                       or raw, if normalize_target is FALSE
#       "reference_negative_regression" -> raw / baseline * dna_2n_value, with
#                       baseline = g1_anchor + baseline_slope *
#                       (dna_norm - dna_2n_value)                (edu mode)
#       "background_reference_regression" -> dna_2n_value * raw / background_fit,
#                       with background_fit = intercept + slope * dna_norm
#                       (floored positive), from the replicate's background
#                       control sample                           (poi mode)
#
# Files expected at:
#   <data_dir>/<prefix><complete_suffix>   (all single cells; every method)
#   <data_dir>/<prefix><g1_suffix>         (G1-gated cells; g1 methods only)
read_and_normalize_sample <- function(
    prefix,
    condition_label,
    data_dir,
    file_suffixes,
    settings,
    baseline_slope = NULL,
    background_model = NULL,
    normalization_method = NULL
) {
  dna_channel <- settings$dna_channel
  target_channel <- settings$target_channel
  dna_2n_value <- get_setting(settings, "dna_2n_value", 1000)
  normalize_target <- isTRUE(get_setting(settings, "normalize_target", TRUE))
  g1_anchor_method <- get_setting(settings, "g1_anchor", "median")

  if (is.null(normalization_method)) {
    normalization_method <- get_setting(
      settings, "normalization_method", "g1_median"
    )
  }

  uses_g1 <- normalization_method %in%
    c("g1_median", "reference_negative_regression")

  complete_file <- file.path(data_dir, paste0(prefix, file_suffixes$complete))
  required_files <- complete_file
  g1_file <- NA_character_
  if (uses_g1) {
    g1_file <- file.path(data_dir, paste0(prefix, file_suffixes$g1))
    required_files <- c(complete_file, g1_file)
  }

  missing_files <- required_files[!file.exists(required_files)]
  if (length(missing_files) > 0) {
    stop(paste0(
      "Missing files for ", condition_label, ":\n",
      paste(missing_files, collapse = "\n")
    ))
  }

  dat <- utils::read.csv(complete_file, check.names = FALSE)
  required_channels <- c(dna_channel, target_channel)
  missing_channels <- setdiff(required_channels, names(dat))
  if (length(missing_channels) > 0) {
    stop(paste0(
      "Missing channels in the all-cells table for ", condition_label, ": ",
      paste(missing_channels, collapse = ", "),
      "\nAvailable columns: ", paste(names(dat), collapse = ", ")
    ))
  }
  dat$target_raw <- dat[[target_channel]]

  g1_median_dna <- NA_real_
  g1_anchor_target <- NA_real_
  cutoff <- NA_real_
  edu_positive_data <- NULL

  if (uses_g1) {
    g1 <- utils::read.csv(g1_file, check.names = FALSE)
    missing_g1 <- setdiff(required_channels, names(g1))
    if (length(missing_g1) > 0) {
      stop(paste0(
        "Missing channels in the G1 table for ", condition_label, ": ",
        paste(missing_g1, collapse = ", ")
      ))
    }
    g1_median_dna <- stats::median(g1[[dna_channel]], na.rm = TRUE)
    g1_anchor_target <- g1_target_anchor(g1[[target_channel]], g1_anchor_method)
    if (!is.finite(g1_median_dna) || g1_median_dna == 0) {
      stop(paste("Invalid G1 DNA median for", condition_label))
    }
    if (!is.finite(g1_anchor_target) || g1_anchor_target == 0) {
      stop(paste("Invalid G1 target anchor for", condition_label))
    }
    dat$dna_norm <- dat[[dna_channel]] / g1_median_dna * dna_2n_value

    if (normalization_method == "g1_median") {
      if (normalize_target) {
        dat$baseline <- rep(g1_anchor_target, nrow(dat))
        dat$target_norm <- dat$target_raw / dat$baseline * dna_2n_value
      } else {
        dat$baseline <- NA_real_
        dat$target_norm <- dat$target_raw
      }
    } else {  # reference_negative_regression
      edu <- NULL
      if (!is.null(file_suffixes$edu_positive)) {
        edu_file <- file.path(data_dir, paste0(prefix, file_suffixes$edu_positive))
        if (file.exists(edu_file)) {
          edu <- utils::read.csv(edu_file, check.names = FALSE)
        }
      }
      normalized <- normalize_edu(
        events = dat,
        g1_events = g1,
        dna_channel = dna_channel,
        target_channel = target_channel,
        baseline_slope = baseline_slope,
        dna_2n_value = dna_2n_value,
        g1_anchor = g1_anchor_method,
        edu_positive_events = edu,
        sample_id = condition_label
      )
      dat <- normalized$data
      edu_positive_data <- normalized$edu_positive
      g1_median_dna <- normalized$g1_median_dna
      g1_anchor_target <- normalized$g1_anchor_target
    }
  } else if (normalization_method == "background_reference_regression") {
    normalized <- normalize_poi(
      events = dat,
      background_model = background_model,
      dna_channel = dna_channel,
      target_channel = target_channel,
      dna_2n_value = dna_2n_value,
      dna_align = get_setting(settings, "poi_dna_align", "per_sample"),
      peak_failure = get_setting(settings, "poi_peak_failure", "error"),
      sample_id = condition_label
    )
    dat <- normalized$data
    cutoff <- normalized$cutoff
  } else {
    stop(paste0("Unknown normalization method: ", normalization_method))
  }

  # Background-subtracted signal (raw - fitted baseline). Used by the
  # quantitation medians when quant_signal = "background_subtracted".
  dat$target_bgsub <- dat$target_raw - dat$baseline

  list(
    data = dat,
    edu_positive = edu_positive_data,
    g1_median_dna = g1_median_dna,
    g1_anchor_target = g1_anchor_target,
    normalization_method = normalization_method,
    baseline_slope = if (is.null(baseline_slope)) NA_real_ else baseline_slope,
    cutoff = cutoff,
    complete_file = complete_file,
    g1_file = g1_file
  )
}


# ---------------------------------------------------------------------------
# POI mode: background-control normalization
# ---------------------------------------------------------------------------

# Estimate the lower (2N) major DNA density peak, as used to place G1 at 2N in
# the MTBP analysis. Ported from the Figure 1 MTBP pipeline.
estimate_lower_dna_peak <- function(x, adjust = 0.75) {
  x <- x[is.finite(x)]
  if (length(x) < 10) stop("Too few finite DNA values to find a 2N peak.")
  density_fit <- stats::density(x, n = 2048, adjust = adjust)
  maxima <- which(diff(sign(diff(density_fit$y))) == -2) + 1
  maxima <- maxima[
    density_fit$x[maxima] >= stats::quantile(x, 0.02) &
    density_fit$x[maxima] <= stats::quantile(x, 0.80)
  ]
  if (!length(maxima)) stop("No DNA density peak could be identified.")
  height_cutoff <- 0.15 * max(density_fit$y[maxima])
  candidates <- maxima[density_fit$y[maxima] >= height_cutoff]
  density_fit$x[min(candidates)]
}

# Fit one replicate's background-control model from its background/control
# sample: estimate the 2N DNA peak, regress raw target on normalized DNA, and
# take the requested quantile of the corrected background as the cutoff line.
fit_background_model <- function(
    prefix,
    condition_label,
    replicate_label,
    data_dir,
    file_suffixes,
    settings
) {
  dna_channel <- settings$dna_channel
  target_channel <- settings$target_channel
  dna_2n_value <- get_setting(settings, "dna_2n_value", 1000)
  background_quantile <- get_setting(settings, "background_quantile", 0.95)

  complete_file <- file.path(data_dir, paste0(prefix, file_suffixes$complete))
  if (!file.exists(complete_file)) {
    stop(paste0("Missing background file for ", condition_label, ":\n",
                complete_file))
  }
  bg <- utils::read.csv(complete_file, check.names = FALSE)
  if (!all(c(dna_channel, target_channel) %in% names(bg))) {
    stop(paste("Background file is missing required channels for",
               condition_label))
  }

  dna_peak <- estimate_lower_dna_peak(bg[[dna_channel]])
  bg$dna_norm <- bg[[dna_channel]] / dna_peak * dna_2n_value
  bg$target_raw <- bg[[target_channel]]

  fit <- stats::lm(target_raw ~ dna_norm, data = bg)
  coefs <- stats::coef(fit)
  if (any(!is.finite(coefs))) {
    stop(paste("Invalid background regression for", replicate_label))
  }

  predicted <- coefs[[1]] + coefs[[2]] * bg$dna_norm
  positive_floor <- max(
    .Machine$double.eps,
    stats::quantile(predicted[predicted > 0], 0.001, na.rm = TRUE)
  )
  predicted <- pmax(predicted, positive_floor)
  corrected <- dna_2n_value * bg$target_raw / predicted
  cutoff <- unname(stats::quantile(corrected, background_quantile, na.rm = TRUE))

  bg$background_fitted <- predicted

  diagnostic_plot <- ggplot2::ggplot(bg, ggplot2::aes(x = dna_norm, y = target_raw)) +
    ggplot2::geom_point(size = 0.25, alpha = 0.18, color = "#2952A3") +
    ggplot2::geom_line(
      ggplot2::aes(y = background_fitted),
      color = "black", linewidth = 0.7
    ) +
    ggplot2::scale_x_continuous(breaks = c(dna_2n_value, 2 * dna_2n_value),
                                labels = c("2N", "4N")) +
    ggplot2::labs(
      title = replicate_label, subtitle = paste("Background:", condition_label),
      x = "Normalized DNA content", y = paste("Raw", target_channel)
    ) +
    ggplot2::theme_classic(base_size = 10)

  list(
    replicate = replicate_label,
    background_condition = condition_label,
    dna_peak = dna_peak,
    intercept = unname(coefs[[1]]),
    slope = unname(coefs[[2]]),
    floor = positive_floor,
    r_squared = summary(fit)$r.squared,
    background_event_n = nrow(bg),
    cutoff = cutoff,
    diagnostic_plot = diagnostic_plot
  )
}

# Fit a background model for every replicate; returns a list keyed by
# replicate_index. The background sample is the replicate's `is_reference` row.
fit_replicate_background_models <- function(
    sample_manifest, data_dir, file_suffixes, settings) {
  replicate_groups <- split(sample_manifest, sample_manifest$replicate_index)

  lapply(replicate_groups, function(replicate_data) {
    background_rows <- replicate_data[replicate_data$is_reference, , drop = FALSE]
    if (nrow(background_rows) != 1) {
      stop(paste(
        "Exactly one background sample must be selected for",
        replicate_data$replicate[[1]],
        "- set `reference:` on each replicate to the background control label."
      ))
    }
    fit_background_model(
      prefix = background_rows$prefix[[1]],
      condition_label = background_rows$condition[[1]],
      replicate_label = background_rows$replicate[[1]],
      data_dir = data_dir, file_suffixes = file_suffixes, settings = settings
    )
  })
}

# Pool normalized target values across a set of samples to derive shared y-axis
# limits, so every panel is directly comparable (used by poi mode). Returns a
# length-2 numeric vector of quantile-based limits.
collect_common_y_limits <- function(
    sample_manifest,
    data_dir,
    file_suffixes,
    settings,
    baseline_slope_list = NULL,
    background_model_list = NULL,
    lower = 0.001,
    upper = 0.999
) {
  values <- unlist(lapply(seq_len(nrow(sample_manifest)), function(i) {
    sample_data <- read_and_normalize_sample(
      prefix = sample_manifest$prefix[[i]],
      condition_label = sample_manifest$condition[[i]],
      data_dir = data_dir, file_suffixes = file_suffixes, settings = settings,
      baseline_slope = if (is.null(baseline_slope_list)) NULL else
        baseline_slope_list[[i]],
      background_model = if (is.null(background_model_list)) NULL else
        background_model_list[[i]]
    )
    v <- sample_data$data$target_norm
    v[is.finite(v) & v > 0]
  }), use.names = FALSE)

  if (length(values) < 2) stop("Too few values to set common y limits.")
  stats::quantile(values, c(lower, upper), names = FALSE)
}


# ---------------------------------------------------------------------------
# EdU mode: reference EdU-negative regression baseline
# ---------------------------------------------------------------------------

# Fit a positive-minimum boundary that traces the lower edge of the EdU+ cloud
# across DNA content, by binning DNA and taking the minimum signal per bin.
fit_positive_minimum_boundary <- function(
    edu,
    x_column = "dna_norm",
    y_column,
    fit_x_range = c(1000, 2000),
    boundary_bins = 20,
    minimum_events_per_bin = 20,
    condition_label = "reference sample"
) {
  missing_columns <- setdiff(c(x_column, y_column), names(edu))
  if (length(missing_columns) > 0) {
    stop(paste0(
      "Missing columns for positive-boundary fitting in ", condition_label,
      ": ", paste(missing_columns, collapse = ", ")
    ))
  }

  x <- edu[[x_column]]
  y <- edu[[y_column]]
  valid <- is.finite(x) & is.finite(y) &
    x >= fit_x_range[[1]] & x <= fit_x_range[[2]]

  fit_data <- data.frame(boundary_x = x[valid], boundary_y = y[valid])
  if (nrow(fit_data) < 2 * minimum_events_per_bin) {
    stop(paste("Too few EdU+ reference events for boundary fitting in",
               condition_label))
  }

  breaks <- seq(fit_x_range[[1]], fit_x_range[[2]], length.out = boundary_bins + 1)
  fit_data$bin <- pmin(
    findInterval(fit_data$boundary_x, breaks, rightmost.closed = TRUE),
    boundary_bins
  )

  bin_counts <- table(fit_data$bin)
  usable_bins <- as.integer(names(bin_counts)[bin_counts >= minimum_events_per_bin])

  boundary_points <- do.call(rbind, lapply(usable_bins, function(bin_number) {
    bin_data <- fit_data[fit_data$bin == bin_number, , drop = FALSE]
    minimum_row <- which.min(bin_data$boundary_y)
    data.frame(
      bin = bin_number,
      boundary_x = bin_data$boundary_x[[minimum_row]],
      boundary_y = bin_data$boundary_y[[minimum_row]]
    )
  }))

  if (is.null(boundary_points) || nrow(boundary_points) < 2) {
    stop(paste("Fewer than two DNA bins had enough events for boundary fitting in",
               condition_label))
  }

  boundary_fit <- stats::lm(boundary_y ~ boundary_x, data = boundary_points)
  if (any(!is.finite(stats::coef(boundary_fit)))) {
    stop(paste("Invalid positive-boundary fit for", condition_label))
  }

  list(fit = boundary_fit, points = boundary_points)
}

# For one replicate's reference sample: fit the positive boundary, classify the
# EdU-negative events, and regress the negative signal on normalized DNA to get
# the shared baseline slope. Returns the slope plus a diagnostic plot.
fit_reference_negative_model <- function(
    prefix,
    condition_label,
    replicate_label,
    data_dir,
    file_suffixes,
    settings
) {
  dna_channel <- settings$dna_channel
  target_channel <- settings$target_channel
  dna_2n_value <- get_setting(settings, "dna_2n_value", 1000)
  fit_x_range <- get_setting(settings, "baseline_fit_x_range", c(1000, 2000))
  boundary_bins <- get_setting(settings, "baseline_boundary_bins", 20)
  min_per_bin <- get_setting(settings, "baseline_minimum_events_per_bin", 20)
  min_negative <- get_setting(settings, "baseline_minimum_negative_events", 100)

  # Reference sample DNA-normalized (target not baseline-corrected yet).
  reference_data <- read_and_normalize_sample(
    prefix = prefix, condition_label = condition_label, data_dir = data_dir,
    file_suffixes = file_suffixes, settings = settings,
    normalization_method = "g1_median"
  )
  complete <- reference_data$data

  # Read the EdU-positive events and DNA-normalize with the same G1 median.
  edu_file <- file.path(data_dir, paste0(prefix, file_suffixes$edu_positive))
  if (!file.exists(edu_file)) {
    stop(paste0("Missing EdU-positive file for ", condition_label, ":\n", edu_file))
  }
  edu <- utils::read.csv(edu_file, check.names = FALSE)
  if (!all(c(dna_channel, target_channel) %in% names(edu))) {
    stop(paste("EdU-positive file is missing required channels for",
               condition_label))
  }
  edu$dna_norm <- edu[[dna_channel]] / reference_data$g1_median_dna * dna_2n_value
  edu$target_raw <- edu[[target_channel]]

  boundary_result <- fit_positive_minimum_boundary(
    edu = edu, x_column = "dna_norm", y_column = "target_raw",
    fit_x_range = fit_x_range, boundary_bins = boundary_bins,
    minimum_events_per_bin = min_per_bin,
    condition_label = paste(replicate_label, condition_label)
  )

  complete$positive_boundary <- stats::predict(
    boundary_result$fit,
    newdata = data.frame(boundary_x = complete$dna_norm)
  )

  is_negative <-
    is.finite(complete$dna_norm) & is.finite(complete$target_raw) &
    is.finite(complete$positive_boundary) &
    complete$dna_norm >= fit_x_range[[1]] & complete$dna_norm <= fit_x_range[[2]] &
    complete$target_raw < complete$positive_boundary

  negative_events <- complete[is_negative, , drop = FALSE]
  if (nrow(negative_events) < min_negative) {
    stop(paste("Too few inferred EdU-negative events in reference sample",
               condition_label, "for", replicate_label))
  }

  negative_fit <- stats::lm(
    target_raw ~ dna_norm,
    data = data.frame(target_raw = negative_events$target_raw,
                      dna_norm = negative_events$dna_norm)
  )
  reference_slope <- unname(stats::coef(negative_fit)[[2]])
  if (!is.finite(reference_slope)) {
    stop(paste("Invalid EdU-negative regression slope for", replicate_label))
  }

  # Diagnostic: grey = all cells, blue = inferred EdU-, green = EdU+,
  # red points/line = positive-minimum boundary, black line = negative fit
  # forced through the sample's G1 anchor at 2N.
  g1_anchor <- reference_data$g1_anchor_target
  line_x <- seq(fit_x_range[[1]], fit_x_range[[2]], length.out = 200)
  boundary_line <- data.frame(
    dna_norm = line_x,
    target_raw = stats::predict(boundary_result$fit,
                                newdata = data.frame(boundary_x = line_x))
  )
  # Baseline actually used downstream: shares the slope, anchored at G1 point.
  baseline_line <- data.frame(
    dna_norm = line_x,
    target_raw = g1_anchor + reference_slope * (line_x - dna_2n_value)
  )

  diagnostic_plot <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = complete,
      ggplot2::aes(x = dna_norm, y = target_raw),
      color = "grey75", size = 0.25, alpha = 0.5
    ) +
    ggplot2::geom_point(
      data = negative_events,
      ggplot2::aes(x = dna_norm, y = target_raw),
      color = "#3480FF", size = 0.3, alpha = 0.6
    ) +
    ggplot2::geom_point(
      data = edu,
      ggplot2::aes(x = dna_norm, y = target_raw),
      color = "#18BFA5", size = 0.3, alpha = 0.55
    ) +
    ggplot2::geom_point(
      data = boundary_result$points,
      ggplot2::aes(x = boundary_x, y = boundary_y),
      color = "#C00000", size = 1.6
    ) +
    ggplot2::geom_line(
      data = boundary_line, ggplot2::aes(x = dna_norm, y = target_raw),
      color = "#C00000", linewidth = 0.7
    ) +
    ggplot2::geom_line(
      data = baseline_line, ggplot2::aes(x = dna_norm, y = target_raw),
      color = "black", linewidth = 0.8
    ) +
    ggplot2::geom_point(
      data = data.frame(dna_norm = dna_2n_value, target_raw = g1_anchor),
      ggplot2::aes(x = dna_norm, y = target_raw),
      shape = 21, size = 2.4, stroke = 0.7, fill = "white", color = "black"
    ) +
    ggplot2::coord_cartesian(xlim = fit_x_range) +
    ggplot2::scale_x_continuous(breaks = c(dna_2n_value, 2 * dna_2n_value),
                                labels = c("2N", "4N")) +
    ggplot2::scale_y_continuous(labels = scales::label_comma()) +
    ggplot2::labs(
      title = replicate_label, subtitle = paste("Reference:", condition_label),
      x = "Normalized DNA content",
      y = paste("Raw", target_channel)
    ) +
    ggplot2::theme_classic(base_size = 10)

  list(
    replicate = replicate_label,
    reference_condition = condition_label,
    slope = reference_slope,
    negative_intercept = unname(stats::coef(negative_fit)[[1]]),
    negative_r_squared = summary(negative_fit)$r.squared,
    negative_event_n = nrow(negative_events),
    g1_anchor_target = g1_anchor,
    diagnostic_plot = diagnostic_plot
  )
}

# Fit a reference model for every replicate; returns a list keyed by
# replicate_index.
fit_replicate_reference_models <- function(
    sample_manifest, data_dir, file_suffixes, settings) {
  replicate_groups <- split(sample_manifest, sample_manifest$replicate_index)

  lapply(replicate_groups, function(replicate_data) {
    reference_rows <- replicate_data[replicate_data$is_reference, , drop = FALSE]
    if (nrow(reference_rows) != 1) {
      stop(paste(
        "Exactly one reference sample must be selected for",
        replicate_data$replicate[[1]],
        "- set `reference:` on each replicate to a matching condition label."
      ))
    }
    fit_reference_negative_model(
      prefix = reference_rows$prefix[[1]],
      condition_label = reference_rows$condition[[1]],
      replicate_label = reference_rows$replicate[[1]],
      data_dir = data_dir, file_suffixes = file_suffixes, settings = settings
    )
  })
}


# ---------------------------------------------------------------------------
# EdU apex line (highest EdU density at the S-phase arc apex)
# ---------------------------------------------------------------------------

# Within a DNA-content window (mid-S by default), find the highest-density EdU+
# signal by taking the mode of the log10 signal density. This marks the apex of
# the EdU+ arc. Ported from the Figure 1 calculate_mid_density helper.
calculate_edu_apex_density <- function(
    edu,
    x_column = "dna_norm",
    y_column = "target_norm",
    mid_x_range = c(1400, 1600),
    density_adjust = 1,
    minimum_events = 10,
    condition_label = "sample"
) {
  x <- edu[[x_column]]
  y <- edu[[y_column]]
  keep <- is.finite(x) & is.finite(y) &
    x >= mid_x_range[[1]] & x <= mid_x_range[[2]] & y > 0
  y_win <- y[keep]

  if (length(y_win) < minimum_events) {
    warning(paste("Too few EdU+ events in the apex window for", condition_label))
    return(list(y = NA_real_, n = length(y_win)))
  }

  dens <- stats::density(log10(y_win), n = 512, adjust = density_adjust)
  list(y = 10^dens$x[which.max(dens$y)], n = length(y_win))
}


# Build a decoration-free panel (data area only): points + the 2N line + any
# horizontal reference lines, with the data mapped exactly to the panel edges
# (expand = FALSE). Used by the plotgardener layout so every panel's plotting
# area is an identical fixed size.
build_naked_panel <- function(
    data, palette, settings, x_limits, y_limits, y_log10, dna_2n_value, hlines
) {
  p <- ggplot2::ggplot(
    data, ggplot2::aes(x = dna_norm, y = target_norm, color = density_color)
  ) +
    ggplot2::geom_point(
      shape = 16, size = get_setting(settings, "point_size", 0.3),
      stroke = 0, alpha = 1
    ) +
    ggplot2::scale_color_gradientn(
      colors = palette, limits = c(0, 1), oob = scales::squish, guide = "none"
    ) +
    ggplot2::geom_vline(
      xintercept = dna_2n_value, linetype = "dashed",
      color = "grey50", linewidth = 0.35
    )
  for (h in hlines) {
    p <- p + ggplot2::geom_hline(
      yintercept = h$y, color = h$color, linewidth = h$linewidth,
      linetype = h$linetype %||% "solid"
    )
  }
  if (y_log10) {
    p <- p + ggplot2::scale_y_log10(expand = ggplot2::expansion(0, 0))
  } else {
    p <- p + ggplot2::scale_y_continuous(expand = ggplot2::expansion(0, 0))
  }
  p +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(0, 0)) +
    ggplot2::coord_cartesian(
      xlim = x_limits, ylim = y_limits, expand = FALSE, clip = "off"
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "none",
                   plot.margin = ggplot2::margin(0, 0, 0, 0))
}


# ---------------------------------------------------------------------------
# Single pseudocolor panel
# ---------------------------------------------------------------------------

# Builds one pseudocolor density plot of the target signal (y) vs normalized
# DNA content (x). `baseline_slope` is used only in EdU mode. Returns both a
# decorated `plot` (cowplot layout / HTML preview) and a `panel` with no
# decorations (plotgardener layout), plus the horizontal-line specs.
make_pseudocolor_plot <- function(
    prefix,
    condition_label,
    data_dir,
    file_suffixes,
    settings,
    baseline_slope = NULL,
    background_model = NULL,
    palette = refined_density_palette()
) {
  sample_data <- read_and_normalize_sample(
    prefix = prefix, condition_label = condition_label, data_dir = data_dir,
    file_suffixes = file_suffixes, settings = settings,
    baseline_slope = baseline_slope, background_model = background_model
  )
  make_pseudocolor_plot_from_data(
    sample_data = sample_data,
    condition_label = condition_label,
    settings = settings,
    palette = palette
  )
}

# Build a panel from an already-normalized sample without reading input files.
make_pseudocolor_plot_from_data <- function(
    sample_data,
    condition_label,
    settings,
    palette = refined_density_palette()
) {
  dat <- sample_data$data

  dna_2n_value <- get_setting(settings, "dna_2n_value", 1000)
  x_limits <- get_setting(settings, "x_limits", c(700, 2250))
  y_limits <- get_setting(settings, "y_limits", NULL)
  y_log10 <- isTRUE(get_setting(settings, "y_log10", TRUE))

  if (is.null(y_limits)) {
    finite_y <- dat$target_norm[is.finite(dat$target_norm) & dat$target_norm > 0]
    if (length(finite_y) < 10) {
      stop(paste("Too few positive target values for", condition_label))
    }
    y_limits <- stats::quantile(
      finite_y,
      c(get_setting(settings, "y_limit_lower_quantile", 0.001),
        get_setting(settings, "y_limit_upper_quantile", 0.999)),
      names = FALSE
    )
  }

  keep <-
    is.finite(dat$dna_norm) & is.finite(dat$target_norm) &
    dat$dna_norm >= x_limits[1] & dat$dna_norm <= x_limits[2] &
    dat$target_norm >= y_limits[1] & dat$target_norm <= y_limits[2]
  if (y_log10) keep <- keep & dat$target_norm > 0

  dat_display <- dat[keep, , drop = FALSE]
  if (nrow(dat_display) < 10) {
    stop(paste("Too few events in the display region for", condition_label))
  }

  density_y <- if (y_log10) log10(dat_display$target_norm) else dat_display$target_norm
  dat_display$density <- compute_point_density(
    x = dat_display$dna_norm, y = density_y,
    n = get_setting(settings, "density_grid_n", 300),
    bandwidth_multiplier = get_setting(settings, "density_bandwidth", 0.5)
  )
  dat_display$density_color <- prepare_density_color(
    density = dat_display$density,
    lower_clip_quantile = get_setting(settings, "density_lower_clip", 0.05),
    upper_clip_quantile = get_setting(settings, "density_upper_clip", 0.90),
    gamma = get_setting(settings, "density_gamma", 0.35)
  )
  dat_display <- dat_display[order(dat_display$density_color), , drop = FALSE]

  # Horizontal reference lines (computed once, used by both layouts).
  hlines <- list()
  if (isTRUE(get_setting(settings, "show_baseline_line", TRUE)) &&
      sample_data$normalization_method == "reference_negative_regression") {
    hlines[[length(hlines) + 1]] <- list(
      y = dna_2n_value, color = "grey40", linewidth = 0.35, linetype = "dotted")
  }
  if (isTRUE(get_setting(settings, "show_cutoff_line", TRUE)) &&
      sample_data$normalization_method == "background_reference_regression" &&
      is.finite(sample_data$cutoff)) {
    hlines[[length(hlines) + 1]] <- list(
      y = sample_data$cutoff,
      color = get_setting(settings, "cutoff_color", "#C00000"),
      linewidth = get_setting(settings, "cutoff_linewidth", 0.65),
      linetype = "solid")
  }
  if (isTRUE(get_setting(settings, "show_edu_apex_line", FALSE)) &&
      sample_data$normalization_method == "reference_negative_regression" &&
      !is.null(sample_data$edu_positive)) {
    apex <- calculate_edu_apex_density(
      edu = sample_data$edu_positive,
      mid_x_range = get_setting(settings, "edu_apex_x_range",
                                c(1.4, 1.6) * dna_2n_value),
      density_adjust = get_setting(settings, "edu_apex_density_adjust", 1),
      condition_label = condition_label)
    if (is.finite(apex$y)) {
      hlines[[length(hlines) + 1]] <- list(
        y = apex$y,
        color = get_setting(settings, "edu_apex_color", "#C00000"),
        linewidth = get_setting(settings, "edu_apex_linewidth", 0.65),
        linetype = "solid")
    }
  }

  x_breaks <- get_setting(settings, "x_breaks", c(dna_2n_value, 2 * dna_2n_value))
  x_labels <- get_setting(settings, "x_labels", c("2N", "4N"))

  # Decorated single-panel plot (cowplot layout and HTML diagnostics).
  plot <- ggplot2::ggplot(
    dat_display,
    ggplot2::aes(x = dna_norm, y = target_norm, color = density_color)
  ) +
    ggplot2::geom_point(
      shape = 16, size = get_setting(settings, "point_size", 0.3),
      stroke = 0, alpha = 1) +
    ggplot2::scale_color_gradientn(
      colors = palette, limits = c(0, 1),
      oob = scales::squish, name = "Relative density") +
    ggplot2::scale_x_continuous(breaks = x_breaks, labels = x_labels,
                                minor_breaks = NULL) +
    ggplot2::geom_vline(xintercept = dna_2n_value, linetype = "dashed",
                        color = "grey50", linewidth = 0.35)
  for (h in hlines) {
    plot <- plot + ggplot2::geom_hline(
      yintercept = h$y, color = h$color, linewidth = h$linewidth,
      linetype = h$linetype %||% "solid")
  }
  plot <- plot +
    ggplot2::coord_cartesian(xlim = x_limits, ylim = y_limits, clip = "off") +
    ggplot2::labs(
      title = condition_label,
      x = get_setting(settings, "x_axis_title", "DNA content"),
      y = get_setting(settings, "y_axis_title", "Target intensity")) +
    ggplot2::theme_classic(base_size = get_setting(settings, "base_font_size", 10)) +
    ggplot2::theme(
      aspect.ratio = 1, legend.position = "none",
      plot.title = ggplot2::element_text(
        face = "bold", hjust = 0.5,
        size = get_setting(settings, "title_font_size", 12)),
      plot.margin = ggplot2::margin(t = 6, r = 10, b = 6, l = 16))
  if (y_log10) {
    plot <- plot +
      ggplot2::scale_y_log10(
        breaks = get_setting(settings, "y_breaks", 10^(2:5)),
        labels = get_setting(settings, "y_labels", scales::label_log()),
        minor_breaks = NULL) +
      ggplot2::annotation_logticks(
        sides = "l", outside = TRUE, short = grid::unit(0.08, "cm"),
        mid = grid::unit(0.14, "cm"), long = grid::unit(0.20, "cm"))
  } else {
    plot <- plot + ggplot2::scale_y_continuous(labels = scales::label_comma())
  }

  # Decoration-free panel for exact fixed-size placement (plotgardener).
  panel <- build_naked_panel(
    data = dat_display, palette = palette, settings = settings,
    x_limits = x_limits, y_limits = y_limits, y_log10 = y_log10,
    dna_2n_value = dna_2n_value, hlines = hlines)

  list(plot = plot, panel = panel, data = dat_display,
       sample_data = sample_data, y_limits = y_limits, hlines = hlines)
}


# ---------------------------------------------------------------------------
# Panel grid (replicate rows x condition columns)
# ---------------------------------------------------------------------------

# Strip per-panel axis titles so a single shared x title / y title can be used,
# keeping the y title only on the first column of each row.
format_panel_plots <- function(plot_list, panel_columns) {
  if (length(panel_columns) != 1 || !is.finite(panel_columns) ||
      panel_columns < 1) {
    stop("panel_columns must be one positive number.")
  }

  lapply(seq_along(plot_list), function(i) {
    current <- plot_list[[i]] + ggplot2::labs(subtitle = NULL, x = NULL)
    is_first_column <- ((i - 1) %% panel_columns) == 0
    if (!is_first_column) current <- current + ggplot2::labs(y = NULL)
    current
  })
}

# Assemble the formatted panels into a grid with an optional bold row label per
# replicate and a single shared x-axis title along the bottom.
make_panel_grid <- function(
    plot_list,
    panel_columns,
    shared_x_title = "DNA content",
    shared_x_size = 10,
    row_labels = NULL,
    row_label_size = 10
) {
  number_of_rows <- ceiling(length(plot_list) / panel_columns)

  if (!is.null(row_labels) && length(row_labels) != number_of_rows) {
    stop("row_labels must contain one label for each row of plots.")
  }

  if (is.null(row_labels)) {
    core_grid <- cowplot::plot_grid(
      plotlist = plot_list, ncol = panel_columns, align = "hv", axis = "tblr"
    )
  } else {
    labeled_rows <- lapply(seq_len(number_of_rows), function(row_index) {
      first_plot <- (row_index - 1) * panel_columns + 1
      last_plot <- min(row_index * panel_columns, length(plot_list))

      row_grid <- cowplot::plot_grid(
        plotlist = plot_list[first_plot:last_plot],
        ncol = panel_columns, align = "hv", axis = "tblr"
      )

      label_panel <- cowplot::ggdraw() +
        cowplot::draw_label(
          row_labels[[row_index]], angle = 90,
          fontface = "bold", size = row_label_size
        )

      cowplot::plot_grid(
        label_panel, row_grid, ncol = 2, rel_widths = c(0.045, 0.955)
      )
    })

    core_grid <- cowplot::plot_grid(
      plotlist = labeled_rows, ncol = 1, align = "v"
    )
  }

  cowplot::ggdraw() +
    cowplot::draw_plot(core_grid, x = 0, y = 0.10, width = 1, height = 0.90) +
    cowplot::draw_label(
      shared_x_title, x = 0.5, y = 0.025, hjust = 0.5, vjust = 0,
      size = shared_x_size
    )
}


# ---------------------------------------------------------------------------
# plotGardener layout: fixed-size panels on an auto-sized canvas
# ---------------------------------------------------------------------------

# Compute the page geometry (all values in inches). Every panel's plotting area
# is exactly `panel_size` x `panel_size`; the canvas grows with the grid.
facs_page_layout <- function(n_cols, n_rows, opts = list()) {
  o <- function(name, default) if (is.null(opts[[name]])) default else opts[[name]]
  panel <- o("panel_size", 1.0)
  gap_x <- o("panel_gap_x", 0.10)
  gap_y <- o("panel_gap_y", 0.14)
  # Default left margin fits log tick marks + y title + row labels (no numbers).
  # If y_axis_numbers is on, widen to ~1.05 to fit the decade labels.
  y_axis_numbers <- isTRUE(o("y_axis_numbers", FALSE))
  margin_left <- o("margin_left", if (y_axis_numbers) 1.05 else 0.62)
  margin_top <- o("margin_top", 0.34)
  margin_right <- o("margin_right", 0.14)
  margin_bottom <- o("margin_bottom", 0.62)

  grid_w <- n_cols * panel + (n_cols - 1) * gap_x
  grid_h <- n_rows * panel + (n_rows - 1) * gap_y

  list(
    panel = panel, gap_x = gap_x, gap_y = gap_y,
    margin_left = margin_left, margin_top = margin_top,
    margin_right = margin_right, margin_bottom = margin_bottom,
    grid_w = grid_w, grid_h = grid_h,
    width = margin_left + grid_w + margin_right,
    height = margin_top + grid_h + margin_bottom,
    panel_x = function(col) margin_left + (col - 1) * (panel + gap_x),
    panel_y = function(row) margin_top + (row - 1) * (panel + gap_y),
    axis_fontsize = o("axis_fontsize", 7),
    title_fontsize = o("title_fontsize", 8),
    col_title_fontsize = o("col_title_fontsize", 8),
    row_label_fontsize = o("row_label_fontsize", 9),
    line_color = o("panel_border_color", "black"),
    line_width = o("panel_border_width", 0.5),
    y_axis_numbers = y_axis_numbers
  )
}

# Draw the whole figure with plotGardener onto the active graphics device.
# `panels` is a row-major list of decoration-free ggplots (one per grid cell).
assemble_facs_plotgardener <- function(
    panels, n_cols, n_rows, condition_labels, replicate_labels,
    x_limits, y_limits, y_log10 = TRUE,
    x_breaks = c(1000, 2000), x_labels = c("2N", "4N"),
    y_breaks = 10^(0:6), y_axis_title = "Signal",
    x_axis_title = "DNA content", layout = NULL
) {
  if (!requireNamespace("plotgardener", quietly = TRUE)) {
    stop("The plotgardener package is required for layout: plotgardener.")
  }
  if (is.null(layout)) layout <- facs_page_layout(n_cols, n_rows)
  L <- layout
  panel <- L$panel

  plotgardener::pageCreate(
    width = L$width, height = L$height, default.units = "inches",
    showGuides = FALSE, xgrid = 0, ygrid = 0
  )

  # Panels.
  for (r in seq_len(n_rows)) {
    for (c in seq_len(n_cols)) {
      idx <- (r - 1) * n_cols + c
      if (idx > length(panels) || is.null(panels[[idx]])) next
      px <- L$panel_x(c)
      py <- L$panel_y(r)
      plotgardener::plotGG(
        plot = panels[[idx]], x = px, y = py,
        width = panel, height = panel,
        just = c("left", "top"), default.units = "inches"
      )
      plotgardener::plotRect(
        x = px, y = py, width = panel, height = panel,
        just = c("left", "top"), default.units = "inches",
        fill = NA, linecolor = L$line_color, lwd = L$line_width
      )
    }
  }

  # Column titles (condition) above the top row.
  for (c in seq_len(n_cols)) {
    plotgardener::plotText(
      label = condition_labels[[c]],
      x = L$panel_x(c) + panel / 2, y = L$margin_top - 0.05,
      just = c("center", "bottom"), default.units = "inches",
      fontface = "bold", fontsize = L$col_title_fontsize
    )
  }

  # Row labels (replicate) rotated at the far left.
  if (!is.null(replicate_labels)) {
    for (r in seq_len(n_rows)) {
      plotgardener::plotText(
        label = replicate_labels[[r]],
        x = 0.14, y = L$panel_y(r) + panel / 2, rot = 90,
        just = c("center", "center"), default.units = "inches",
        fontface = "bold", fontsize = L$row_label_fontsize
      )
    }
  }

  # X axis (2N / 4N) under the bottom row, per column.
  y_bottom <- L$panel_y(n_rows) + panel
  x_span <- diff(x_limits)
  for (c in seq_len(n_cols)) {
    px <- L$panel_x(c)
    for (k in seq_along(x_breaks)) {
      frac <- (x_breaks[[k]] - x_limits[[1]]) / x_span
      if (frac < -1e-9 || frac > 1 + 1e-9) next
      xk <- px + frac * panel
      plotgardener::plotSegments(
        x0 = xk, y0 = y_bottom, x1 = xk, y1 = y_bottom + 0.04,
        default.units = "inches", lwd = L$line_width
      )
      plotgardener::plotText(
        label = x_labels[[k]], x = xk, y = y_bottom + 0.06,
        just = c("center", "top"), default.units = "inches",
        fontsize = L$axis_fontsize
      )
    }
  }
  plotgardener::plotText(
    label = x_axis_title, x = L$margin_left + L$grid_w / 2,
    y = y_bottom + 0.24, just = c("center", "top"),
    default.units = "inches", fontsize = L$title_fontsize
  )

  # Y axis at the left of each row.
  if (y_log10) {
    y_lo <- log10(y_limits[[1]]); y_hi <- log10(y_limits[[2]])
    if (isTRUE(L$y_axis_numbers)) {
      # Numbered decade ticks.
      y_ticks <- y_breaks[y_breaks >= y_limits[[1]] & y_breaks <= y_limits[[2]]]
      y_tick_labels <- scales::label_comma()(y_ticks)
      for (r in seq_len(n_rows)) {
        py <- L$panel_y(r)
        for (k in seq_along(y_ticks)) {
          yk <- py + panel * (1 - (log10(y_ticks[[k]]) - y_lo) / (y_hi - y_lo))
          plotgardener::plotSegments(
            x0 = L$margin_left - 0.04, y0 = yk, x1 = L$margin_left, y1 = yk,
            default.units = "inches", lwd = L$line_width)
          plotgardener::plotText(
            label = y_tick_labels[[k]], x = L$margin_left - 0.06, y = yk,
            just = c("right", "center"), default.units = "inches",
            fontsize = L$axis_fontsize)
        }
      }
    } else {
      # Classic log tick marks with no numbers: long ticks at decades (1),
      # medium at 5, short at 2-4 and 6-9 (as in ggplot2::annotation_logticks).
      decades <- seq(floor(y_lo), ceiling(y_hi))
      long_len <- 0.075; mid_len <- 0.050; short_len <- 0.030
      for (r in seq_len(n_rows)) {
        py <- L$panel_y(r)
        for (k in decades) {
          for (j in 1:9) {
            v <- j * 10^k
            if (v < y_limits[[1]] || v > y_limits[[2]]) next
            len <- if (j == 1) long_len else if (j == 5) mid_len else short_len
            yk <- py + panel * (1 - (log10(v) - y_lo) / (y_hi - y_lo))
            plotgardener::plotSegments(
              x0 = L$margin_left - len, y0 = yk, x1 = L$margin_left, y1 = yk,
              default.units = "inches", lwd = L$line_width)
          }
        }
      }
    }
  }

  # Y-axis title, rotated, at the far left.
  plotgardener::plotText(
    label = y_axis_title, x = 0.30,
    y = L$margin_top + L$grid_h / 2, rot = 90,
    just = c("center", "center"), default.units = "inches",
    fontsize = L$title_fontsize
  )

  invisible(L)
}


# ---------------------------------------------------------------------------
