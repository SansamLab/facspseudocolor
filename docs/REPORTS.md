# Reusable reports and RDS artifacts

## Reader-facing report convention

Reader-facing pH3, EdU, and future POI reports show their approved figures,
concise explanatory text, warnings, and availability states; they do not print
raw or source tables. The completed analysis RDS and approved machine-readable
exports retain the full source tables, QC records, and provenance needed for
audit. Every reader-facing report ends with a **Software used** section that
lists the package, R, Quarto, and plotting/rendering software versions. Future
POI report templates must follow the same convention.

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
| `facs_configurator.qmd` | **Unfinished/experimental.** Guided FlowJo-gated EdU configuration prototype; not part of the stable report workflow. |
| `facs_complete.qmd` | Pseudocolor, quantitation, tables, and provenance. |
| `facs_pseudocolor.qmd` | Editable signal-versus-DNA panels. |
| `facs_quantitation.qmd` | Phase and whole-population signal summaries. |
| `facs_cell_cycle.qmd` | Cell-cycle phase percentages and gates. |
| `facs_diagnostics.qmd` | Gate assignments, fits, input checks, and warnings. |
| `facs_ph3_4n.qmd` | Exact FlowJo pH3 gate intersected with configured G2/M DNA, pseudocolor, and percentage. |
| `facs_ph3_output_contract.qmd` | Canonical four-panel pH3 condition report for a completed production output-contract analysis. |
| `facs_edu_pseudocolor_output_contract.qmd` | Stable legacy-compatible EdU report retained at its established filename for 0.1-line render commands. |
| `facs_edu_standard_v2.qmd` | Current enhanced EdU report with overview galleries, gating and correction diagnostics, quantitative cards, interactive apex and phase-gate controls, and PDF downloads. |
| `facs_edu_model_gated.qmd` | Compact experimental model-gating report emphasizing model-derived panels and gating provenance. This remains separate from the optimized full EdU report. |

The report-version policy and complete EdU series are documented in
[`EDU_REPORT_VERSIONS.md`](EDU_REPORT_VERSIONS.md) and recorded for tools in
`inst/quarto/report-catalog.yml`. New designs receive new filenames; established
report templates are not replaced.

For `facs_edu_standard_v2.qmd`, `embed_pseudocolor_pdf_downloads: true` embeds a
vector PDF for every displayed overview card and adds a download button below
that card. Each PDF uses a 3-by-3-inch page. Enabled apex-line and phase-gate
controls also add the corresponding plain, apex, phase, and combined variants
to the all-plots download. The PDFs are embedded in the
self-contained HTML; no companion download directory is created.
Because the PDFs preserve plotted event coordinates at higher fidelity than the
overview PNGs, the report displays a sharing warning when downloads are enabled.
Rendering fails closed if any PDF exceeds 8 MiB or the combined embedded PDFs
exceed 96 MiB when phase-gate variants are enabled (64 MiB otherwise).

## Interactive configuration in RStudio

Set the RStudio working directory to the experiment folder,
run `facspseudocolor::set_facs_configurator_directory()` in the RStudio Console,
then open `inst/quarto/facs_configurator.qmd` and click **Run Document**. The
Shiny-backed document discovers FlowJo `.wsp` files recursively, reads their
FCS acquisition names and detector channels, and creates a FlowJo-gated EdU
configuration. It proposes one row per acquisition, but validation and download
remain disabled until the user explicitly confirms the DNA-H detector and every
included acquisition, condition, and replicate assignment against the
experiment record. It saves those report settings in `config.yml` and never
substitutes another detector when DNA-H is not recognized. Technical
acquisitions are processed independently and averaged after quantitation within
their biological replicate.

This first version configures analysis and report rendering; it does not run or
validate a FlowJo export operation. The required `flowjo_gated_csv` and
`all_events_csv` directories must be prepared separately before the generated
report command is run.

## All-events Single Cells gating card (`flowjo_all_events_dir`)

Passing `-P flowjo_all_events_dir=<dir> -P flowjo_dna_height_channel=<channel>`
to a FlowJo-mode EdU render (or the equivalent `report:` config keys) makes
the Single Cells gating card show real DNA-A-versus-DNA-H doublet
discrimination instead of the config's DNA/target axes. Omitting it is valid
and renders a labeled fallback banner instead of failing -- but the two CSV
inputs have exact, otherwise-undocumented column requirements
(`facs_report_edu_gating_cards()` in `inst/quarto/_report-setup.R`):

- Each `<prefix>_all_events.csv` under `flowjo_all_events_dir` needs a
  `sample_id` column (must exactly equal the configured `fcs:` basename), an
  `event_index` column (as text), the configured `dna_channel` column, and the
  `dna_height_channel` column.
- The `complete`/Single-Cells per-population CSV (in the ordinary
  `data_dir`/`flowjo_gated_csv`) must *also* carry `event_index` and the
  `dna_height_channel` value per event, in addition to the DNA/target columns
  every population CSV already needs -- because accepted Single Cells events
  are overlaid directly onto the all-events background by matching
  `event_index`. The `g1` and `edu_positive`/`s_phase` CSVs do not need
  `dna_height_channel`.

A missing column or an empty population now raises an error that names the
exact population and column(s) involved (for example, `Single Cells missing
column(s): DNA height`) rather than the generic "lacks required channels or
events" message from earlier package versions.

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
