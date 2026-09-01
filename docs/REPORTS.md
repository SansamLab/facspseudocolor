# Reusable reports and RDS artifacts

Every EdU analysis automatically calculates the seven approved canonical tables:
historical five-gate assigned composition, whole-Single-Cells composition,
six-gate assigned composition, Early/Mid/Late-S regional positivity, overall
G1-through-G2/M positivity, computed-positive regional intensity, and whole
computed-positive-population intensity. Acquisition-level versions are retained
before the unchanged unweighted technical-acquisition aggregation.

Canonical intensities use unoffset background-subtracted EdU. Event-level
pseudocolor axes explicitly identify background-subtracted EdU with the display
offset; the offset is recorded in table and analysis provenance and is not used
for classification or quantitative medians. The legacy background-divided signal
is calculated only when explicitly selected with `quant_signal: normalized`.
Optional reference ratios are calculated afterward from the selected signal;
with the defaults, this means background subtraction followed by normalization
to the explicitly selected experimental reference condition. This is distinct
from a POI background control and is optional for EdU.

For pH3 analyses, the reports instead display all supported pH3 quantities:
overall pH3-positive percentage, percentages of all Single Cell events that are
pH3-positive in each DNA phase, unassigned pH3-positive events, the FlowJo gate
diagnostic, and G2/M boundary sensitivity.

## Choose one input

Every report requires exactly one of:

- `config`: validate CSVs and perform the analysis; or
- `analysis_rds`: reuse a previously completed analysis.

Supplying both or neither is an error.

## Templates

| Template | Purpose |
|---|---|
| `facs_configurator.qmd` | Interactive RStudio GUI for discovering, ordering, typing, validating, and saving samples. |
| `facs_complete.qmd` | Pseudocolor, quantitation, tables, and provenance. |
| `facs_pseudocolor.qmd` | Editable signal-versus-DNA panels. |
| `facs_quantitation.qmd` | Phase and whole-population signal summaries. |
| `facs_cell_cycle.qmd` | Cell-cycle phase percentages and gates. |
| `facs_diagnostics.qmd` | Gate assignments, fits, input checks, and warnings. |
| `facs_ph3_4n.qmd` | Exact FlowJo pH3 gate intersected with configured G2/M DNA, pseudocolor, and percentage. |
| `facs_ph3_output_contract.qmd` | Canonical four-panel pH3 condition report for a completed production output-contract analysis. |

## Interactive configuration in RStudio

Set the RStudio working directory to the folder containing the exported CSVs,
run `facspseudocolor::set_facs_configurator_directory()` in the RStudio Console,
then open `inst/quarto/facs_configurator.qmd` and click **Run Document**. The
Shiny-backed document discovers FlowJo `.wsp` files recursively, reads their
FCS sample names, populations, and detector channels; sets the overall EdU,
POI, or pH3 analysis type; then uses three explicit stages to create sample
names, assign one or more workspace FCS acquisitions to each sample, identify
each file's biological and technical replicate, and define each sample's
mode-specific role and order. Technical acquisitions are processed independently
and averaged after quantitation within their biological replicate. It selects channels, validates
the result; and downloads a standard YAML configuration containing the FlowJo
export block. Population CSVs are generated later by the repository's external
FlowJo orchestration step rather than serving as GUI inputs. Save the YAML with
the analysis to preserve the exact setup.

## Example gallery

All previews below are generated from the public POI example data by
`tools/build-report-gallery.R`.

| Template | Preview |
|---|---|
| Complete | ![Complete report preview](gallery/complete.png) |
| Pseudocolor | ![Pseudocolor preview](gallery/pseudocolor.png) |
| Quantitation | ![Quantitation preview](gallery/quantitation.png) |
| Cell cycle | ![Cell-cycle preview](gallery/cell-cycle.png) |
| Diagnostics | ![Diagnostics preview](gallery/diagnostics.png) |

Installed templates are available through:

```r
system.file("quarto", "facs_complete.qmd", package = "facspseudocolor")
```

## Appearance options

Template appearance parameters default to `null`. Null values are ignored, so
package defaults remain authoritative. An optional appearance YAML can provide
report-specific overrides without changing the analysis RDS.

Condition colors use a dynamically sized viridis palette by default. The
package determines the number and ordering of conditions from the analysis;
templates never hard-code sample labels.

```yaml
condition_palette: "colorblind"
y_limits: [600, 30000]
condition_colors:
  "Untreated": "#808080"
```

The package discovers condition labels and panel counts from the analysis.
Named condition colors may be partial; unmatched conditions receive stable
automatic colors. Unknown labels are rejected.

## EdU CSV tables

`save_facs_results(..., output_csv_dir = "results/tables")` writes the seven
canonical aggregate tables and their acquisition-level forms. Deprecated
`phase_percentages`, `phase_medians`, and `whole_medians` aliases are omitted
unless `include_deprecated_csv = TRUE` is explicitly requested. Background-
divided and legacy Figure 1 populations are not part of this automatic export.

## Durable artifacts

`save_facs_analysis()` writes the scientific source of truth, including event
data, fits, all quantitation, configuration, warnings, and provenance.

`save_facs_figure_bundle()` writes editable individual ggplots, summarized
plot-ready tables, resolved appearance, and Plotgardener layout metadata. It
does not duplicate event-level data.

```r
analysis <- read_facs_analysis("results/analysis.rds")
bundle <- build_facs_figure_bundle(analysis)
save_facs_figure_bundle(bundle, "results/figures.rds")

panel <- get_facs_panel(bundle, "rep1_NT")
panel + ggplot2::theme_classic(base_size = 8)
```

Multipanel bundles are assembled with Plotgardener. Individual components stay
as editable ggplot objects for later publication layouts.

The gallery currently uses the public POI example. A pH3 gallery example will
require a suitably licensed real pH3 dataset; production examples are not
silently generated or simulated.

The focused pH3 report additionally requires exactly one
`gate_geometry_csv`, created separately from the FlowJo workspace. Its optional
gate fill, opacity, label color, label size, and label precision parameters
default to `null`, preserving package presentation defaults until overridden.
The optional `gate_top_inset_fraction` controls only how far the visible upper
outline sits below the y-axis ceiling.

## pH3 output-contract report

`facs_ph3_output_contract.qmd` consumes one completed production pH3 analysis
and shows only the four owner-confirmed condition-level outcomes: 4N and
below-4N pH3-positive prevalence, followed by the corresponding pH3 signal
outcomes. It retains biological-replicate points, condition summaries, and
concise availability/basis QC. Its plot layer does not recalculate pH3
results, write CSV/RDS/JSON artifacts, show inferential statistics, or
substitute legacy pH3 plots. Supply exactly one `config` or `analysis_rds`;
`config` runs the configured analysis, while `analysis_rds` renders an already
completed analysis without repeating it.
