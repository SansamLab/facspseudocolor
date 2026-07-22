current_settings_from_config <- function(cfg) {
  list(
    dna_channel = cfg$dna_channel,
    target_channel = cfg$target_channel,
    normalization_method = if (cfg$plot_type == "edu") {
      "reference_negative_regression"
    } else {
      "background_reference_regression"
    },
    normalize_target = cfg$normalize_target,
    y_log10 = cfg$y_log10,
    dna_2n_value = cfg$dna_2n_value,
    g1_anchor = cfg$g1_anchor,
    x_limits = as.numeric(unlist(cfg$x_limits)),
    y_limits = if (is.null(cfg$y_limits)) NULL else as.numeric(unlist(cfg$y_limits)),
    point_size = cfg$point_size,
    density_bandwidth = cfg$density_bandwidth,
    density_lower_clip = cfg$density_lower_clip,
    density_upper_clip = cfg$density_upper_clip,
    density_gamma = cfg$density_gamma,
    baseline_fit_x_range = as.numeric(unlist(cfg$baseline_fit_x_range)),
    baseline_boundary_bins = cfg$baseline_boundary_bins,
    baseline_minimum_events_per_bin = cfg$baseline_minimum_events_per_bin,
    baseline_minimum_negative_events = cfg$baseline_minimum_negative_events,
    background_quantile = cfg$background_quantile,
    poi_dna_align = cfg$poi_dna_align,
    y_limit_lower_quantile = cfg$y_limit_lower_quantile,
    y_limit_upper_quantile = cfg$y_limit_upper_quantile,
    show_edu_apex_line = cfg$show_edu_apex_line,
    edu_apex_x_range = as.numeric(unlist(cfg$edu_apex_x_range)),
    edu_apex_density_adjust = cfg$edu_apex_density_adjust
  )
}

minimal_config <- function(plot_type = "edu") {
  samples <- list(
    list(label = "Reference", prefix = "reference"),
    list(label = "Treatment", prefix = "treatment")
  )
  list(
    plot_type = plot_type,
    data_dir = "explicit/data/path",
    dna_channel = "DNA",
    target_channel = "Target",
    target_name = "Target",
    output_pdf = "results/panels.pdf",
    output_png = "results/panels.png",
    replicates = list(list(
      label = "Replicate 1", reference = "Reference", samples = samples
    ))
  )
}

read_example_config <- function(filename) {
  path <- system.file("config", filename, package = "facspseudocolor")
  if (!nzchar(path)) stop("Installed example config not found: ", filename)
  yaml::read_yaml(path)
}

example_data_dir <- function(name) {
  path <- system.file("extdata", name, package = "facspseudocolor")
  if (!nzchar(path)) stop("Installed example data not found: ", name)
  path
}
