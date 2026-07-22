# Rebuild README gallery images from the public POI example only.
devtools::load_all(".")

analysis <- analyze_facs_experiment("config_poi.yml")
style <- resolve_facs_appearance(
  analysis,
  list(phase_lineplot = FALSE)
)
bundle <- build_facs_figure_bundle(analysis, appearance = style)
directory <- "docs/gallery"
dir.create(directory, recursive = TRUE, showWarnings = FALSE)

save_page <- function(path, page, width, height, resolution = 160) {
  grDevices::png(path, width = width, height = height, units = "in",
                 res = resolution, bg = "white")
  tryCatch(print(page), finally = grDevices::dev.off())
}

save_page(file.path(directory, "pseudocolor.png"), bundle$pseudocolor,
          bundle$pseudocolor$page_layout$width,
          bundle$pseudocolor$page_layout$height)

ggplot2::ggsave(
  file.path(directory, "quantitation.png"),
  bundle$quantitation$phase_median,
  width = 7.2, height = 4.5, dpi = 160, bg = "white"
)
ggplot2::ggsave(
  file.path(directory, "cell-cycle.png"),
  bundle$quantitation$phase_percent,
  width = 7.2, height = 4.5, dpi = 160, bg = "white"
)

diagnostic <- plot_facs_gate_assignments(analysis, appearance = style, seed = 1)
ggplot2::ggsave(
  file.path(directory, "diagnostics.png"), diagnostic,
  width = 9, height = 5.5, dpi = 160, bg = "white"
)

complete_plots <- c(
  bundle$pseudocolor$decorated_plots[seq_len(min(4, length(bundle$pseudocolor$decorated_plots)))],
  bundle$quantitation[c("phase_median", "whole_median")]
)
complete <- assemble_facs_panels(
  complete_plots, ncol = 2,
  layout_options = list(panel_size = 2.8, panel_gap_x = 0.25, panel_gap_y = 0.25,
                        margin_left = 0.2, margin_top = 0.2,
                        margin_right = 0.2, margin_bottom = 0.2)
)
save_page(file.path(directory, "complete.png"), complete,
          complete$layout$width, complete$layout$height)
