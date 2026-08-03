# Normalized pseudocolor plots: signal vs DNA content

A reusable R package with a Quarto front end for making publication-style **pseudocolor
(2D kernel-density) plots** of a single channel (EdU incorporation, a
protein-of-interest signal, or pH3) versus DNA content, laid out as a
**replicate × condition** panel grid. DNA content is normalized per sample so
the G1 population sits at **2N**. By default, the observed G2/M signal is
centered at **1.9× the G1 position** to accommodate typical fluorescence
compression; this is configurable.

This is a generalized distillation of the pseudocolor-plotting code used to
build Figure 1 of the MTBP/EdU project, reproducing that figure's EdU
normalization and default color scheme.

## Three processing modes (`plot_type`)

**`edu`** — fits background independently for every acquisition. The package
uses that acquisition's EdU-positive gate to fit a boundary separating EdU+
from EdU− events, regresses the readily identifiable EdU-negative population
against normalized DNA content, and anchors the fitted baseline at that
acquisition's G1 signal. The default plots and quantitation use raw EdU minus
this fitted baseline. No separate EdU background-control sample is required.

**`poi`** — protein-of-interest, reproducing the Figure 1 CB/Total MTBP
pipeline. Each replicate has an untagged **background control** sample (named
via `reference:`). The tool finds the 2N DNA density peak in that control, fits a
line of background signal versus DNA content, and corrects every cell as
`1000 × signal / background_fit(DNA)`. Cells at the autofluorescence background
sit at ~1000; a red **95th-percentile cutoff line** is drawn on each panel, and
panels share one y-axis for comparability. The background panel is hidden by
default. Needs only the `<prefix>_single_cells.csv` file per sample (no G1 or
EdU-positive files).

**`ph3`** — mitotic-marker analysis using a user-drawn FlowJo pH3-positive
polygon. Each sample supplies Single Cells, G1, and pH3-positive populations.
The G1 gate normalizes DNA; the package reports overall pH3 positivity and the
percentage of all Single Cells that are pH3-positive in G1, Early S, Mid S,
Late S, G2/M, or Unassigned. PH3 phase boundaries are explicit and required.
See [`docs/PH3_MODE.md`](docs/PH3_MODE.md).

A ready-to-run POI example (Total MTBP) is included:

```bash
quarto render pseudocolor_plots.qmd -P config:config_poi.yml
```

## Project components

- `config.yml` — **the one file you edit.** A commented parameters file holding
  every setting: data folder, channel names, mode, palette, replicates, outputs.
  Set up for the EdU example.
- `config_poi.yml` — a ready-to-run POI (Total MTBP) example config.
- `examples/config_ph3.yml` — a PH3 configuration template with explicit
  placeholder paths and phase boundaries.
- `inst/quarto/facs_ph3_4n.qmd` — focused exact-FlowJo-gate pH3 report.
- `pseudocolor_plots.qmd` — a thin report that calls the installed package.
- `inst/quarto/` — focused complete, pseudocolor, quantitation, cell-cycle,
  diagnostics, exact-gate pH3, and interactive configurator templates. The
  configurator supports guided entry or an Excel table prepopulated from the
  selected FlowJo workspace.
- `examples/appearance/` — optional presentation-only YAML examples.
- `R/` — the `facspseudocolor` package implementation.
- `python/export_flowjo_populations.py` — the unchanged optional event exporter.
- `python/export_flowjo_gate_geometry.py` — separate optional exact gate
  geometry extractor.
- `tools/flowjo-orchestration.R` — optional repository-level Python launcher;
  it is not installed with the R package.
- `inst/extdata/example/` — EdU example data (2 replicates × 3 conditions).
- `inst/extdata/example_poi/` — Total MTBP example data (2 replicates, background + 3).
- `results/` — where rendered figures are written.

## Installation

- [Quarto](https://quarto.org)
- R with these required packages: `ggplot2`, `MASS`, `scales`, `cowplot`, and
  `yaml`. The default `cowplot` layout does not require Bioconductor packages.

```r
install.packages(c("ggplot2", "MASS", "scales", "cowplot", "yaml"))
```

`plotgardener` is optional. Install it only when fixed-size publication panels
or Plotgardener-backed figure bundles are needed:

```r
install.packages("BiocManager")
BiocManager::install("plotgardener")
```

Install the package from the project directory:

```bash
R CMD INSTALL .
```

Python is not required to install, load, test, or use the package with existing
CSV files.

## Unit tests

The unit and integration suite checks configuration and input validation,
normalization, quantitation, plotting, the Python boundary, and numerical
regressions for all three modes. Run it from the project directory:

```bash
Rscript -e 'devtools::test()'
R CMD build .
R CMD check --no-manual facspseudocolor_*.tar.gz
```

The tests require `testthat` and `withr`:

```r
install.packages(c("testthat", "withr"))
```

## Input data format

Per sample you provide CSV files in your data folder. Which files are needed
depends on the mode:

| File | Contents | Needed for |
|------|----------|-----------|
| `<prefix>_single_cells.csv` | all single cells | all modes |
| `<prefix>_g1.csv` | the G1-gated single cells (used to normalize) | `edu` and `ph3` modes |
| `<prefix>_edu_positive.csv` | the EdU-positive cells | `edu` mode, every acquisition |
| `<prefix>_ph3_positive.csv` | user-gated pH3-positive cells | `ph3` mode |

`poi` mode needs **only** `<prefix>_single_cells.csv` (no G1 or EdU-positive
files); it subtracts a fitted background from the matched background control.

`ph3` mode needs Single Cells, G1, and pH3-positive files for every sample.
The pH3 gate is defined in FlowJo; the package never guesses a positivity
threshold.

Every file must contain a column for the **DNA channel** and a column for the
**target channel**, named exactly as they appear in your export. For example:

```
Propidium Iodide PE-A,Alexa-647 APC-A
404416.9,85233.4
...
```

Export these from FlowJo (or your gating tool). The `<prefix>` is the file-name
stem you assign in the config, e.g. prefix `rep1_NT` reads
`rep1_NT_single_cells.csv` (and, in `edu` mode, `rep1_NT_g1.csv`).

### Building CSVs from raw FCS + a FlowJo workspace (optional)

If you only have FCS files and a FlowJo `.wsp` (no CSVs yet), the repository-level
orchestration script can build the per-sample CSVs. Add a `flowjo:` block to the
config; when explicitly enabled, it reconstructs saved gates with FlowKit (via
`python/export_flowjo_populations.py`) and writes each `<prefix>_single_cells.csv`
with the DNA and target detector values, then plots. See `config_fig3_2.yml`
for a worked example:

```yaml
flowjo:
  source_dir: "../path/to/fcs_and_wsp_folder"
  workspace: "workspace.wsp"
  python: "../.venv/bin/python"        # python with flowkit/pandas/lxml
  population: "Single Cells"
  dna_source_channel: "FL2-A"          # raw FCS detector for DNA
  target_source_channel: "FL1-A"       # raw FCS detector for the target
  rebuild: true                        # false = skip once CSVs exist

# each sample additionally names its source FCS file:
#   - { label: "siCtrl", prefix: "sictrl", fcs: "... siCtrl.fcs" }
```

This optional step needs **python3 with `flowkit`, `pandas`, and `lxml`** and
must be requested explicitly:

```bash
quarto render pseudocolor_plots.qmd -P run_flowjo_export:true
```

The installed R package never launches Python. See
`docs/PYTHON_INTERFACE.md` for the exact boundary.

## Quick start

Render the example as-is:

```bash
quarto render pseudocolor_plots.qmd
```

This reads `config.yml` and produces `results/pseudocolor_panels.pdf` and
`.png`, plus an HTML report (which echoes the exact settings used).

## Reusable reports and editable artifacts

Every analysis calculates all supported phase medians, whole-population
medians, and phase percentages for the configured signal, which is
background-subtracted by default. Reference normalization, when requested, is
applied afterward to the background-subtracted medians.

To configure an experiment with a GUI, set RStudio's working directory to the
folder containing the exported CSVs, run
`facspseudocolor::set_facs_configurator_directory()`, then open
`inst/quarto/facs_configurator.qmd` and click **Run Document**. It starts from a
FlowJo `.wsp`, discovers its FCS samples, populations, and channels, and can set the overall
analysis type, assign sample roles and replicates, set order, validate the
configuration, and download the resulting YAML file.

Analysis RDS files are the durable scientific source of truth. Figure-bundle
RDS files contain editable individual ggplots, plot-ready summaries, and
Plotgardener layout metadata without duplicating event-level data.

See [`docs/REPORTS.md`](docs/REPORTS.md) for the report catalog, appearance
overrides, RDS workflow, and examples.

The saving layer refuses to overwrite existing results. To deliberately replace
the configured PDF and PNG:

```bash
quarto render pseudocolor_plots.qmd -P overwrite:true
```

The same analysis can be run directly from R without plotting or saving:

```r
library(facspseudocolor)
analysis <- analyze_facs_experiment("config.yml")
analysis <- quantify_cell_cycle(analysis)
plots <- plot_pseudocolor_panels(analysis)
```

## Using it with your own data

1. Put your CSVs in a folder (e.g. `data/mydata/`).
2. Open **`config.yml`** and edit (it is fully commented):
   - `data_dir` — your data folder.
   - `dna_channel`, `target_channel` — the exact column names in your CSVs.
   - `target_name` — label for the y-axis (e.g. `"MTBP"`, `"γH2AX"`).
   - `plot_type` — `"edu"` or `"poi"`.
   - `replicates` — list your biological replicates, each with the same set of
     conditions in the same order (`label` = panel title, `prefix` = file stem).
3. Render with `quarto render pseudocolor_plots.qmd`.

### Multiple configurations

Keep a separate settings file per target/experiment by copying `config.yml`,
then point the document at it:

```bash
quarto render pseudocolor_plots.qmd -P config:config_MTBP.yml
```

### One replicate only

Use a `replicates` list with one entry. This keeps the scientifically required
reference condition explicit:

```yaml
replicates:
  - label: "Experiment"
    reference: "NT"
    samples:
      - { label: "NT", prefix: "NT" }
      - { label: "2h", prefix: "2h" }
```

## Key options (in `config.yml`)

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `layout` | `"cowplot"` | `"cowplot"` (installation-safe aligned grid) or optional `"plotgardener"` (fixed 1×1-inch panels, identical plot areas, auto-sized canvas). |
| `layout_options.panel_size` | `1.0` | (plotgardener) size of each panel's plotting area, in inches. |
| `plot_type` | `"edu"` | `"edu"` (EdU-negative regression baseline), `"poi"` (background-control), or `"ph3"` (FlowJo-defined pH3-positive gate). |
| `palette` | `"refined"` | `"refined"` (Figure 1 default), `"flowjo"`, or a custom hex list. |
| `g1_anchor` | `"median"` | (edu) G1 anchor the baseline passes through: `"median"` (Figure 1) or `"mode"`. |
| `background_quantile` | `0.95` | (poi) percentile of corrected background used for the cutoff line. |
| `show_reference_panel` | mode-based | Show the reference/control as its own panel. Default: `true` for edu, `false` for poi. |
| `y_limits` | auto | Fixed shared y-axis for all panels, e.g. `[700, 6000]`. If unset, limits are automatic (see below). |
| `y_limit_lower_quantile` | `0.001` | Lower quantile of the pooled signal for the automatic shared y-axis. |
| `y_limit_upper_quantile` | `0.999` | Upper quantile of the pooled signal for the automatic shared y-axis. |
| `poi_dna_align` | `per_sample` | (poi) align each sample by its own 2N/G1 peak, or `shared_background` for the old behavior. |
| `poi_peak_failure` | `error` | Stop if a per-sample POI DNA peak cannot be detected. `use_background` explicitly enables the legacy fallback. |
| `show_phase_gates` | `false` | Overlay the G1 / Early S / Mid S / Late S / G2/M gates on every panel. EdU uses 2D gates; POI and pH3 use DNA-content bands. |
| `ph3_boundary_sensitivity_fraction` | `0.05` | (pH3) Diagnostic displacement of each G2/M DNA boundary, expressed as a fraction of the configured 2N value. The configured gate remains the primary result. |
| `show_gate_assignment` | `false` | Diagnostic scatter (per sample) coloring each event by the phase gate it falls in. |
| `quantify_phase_median` | `false` | Legacy switch controlling display in the original report; values are always calculated. |
| `quantify_whole_median` | `false` | Legacy switch controlling display in the original report; values are always calculated. |
| `quantify_phase_percent` | `false` | Legacy switch controlling display in the original report; values are always calculated. |
| `quant_signal` | `background_subtracted` | Signal used for the quantitation medians: `background_subtracted` (raw − fitted baseline, then ratio to reference; matches Figure 1) or `normalized` (the divided signal shown in panels). Percentages unaffected. |
| `quantify_reference_normalized` | `false` | Also make median bar plots normalized within each replicate to the reference (reference = 1). |
| `quant_reference_condition` | reference | Condition to normalize to; defaults to the replicate reference (e.g. `siCtrl NT`). |
| `quant_show_points` | `true` | Show the individual replicate points on the quantitation bar plots. |
| `quant_phase_lineplot` | `true` | Also draw a line-plot of the per-phase quantitation (phase on x, one line per condition) alongside the bar chart. |
| `bar_colors` | (palette) | Named map of condition/sample name → bar color, e.g. `{"NT": "#999999", "2h Auxin": "#3480FF"}`. |
| `y_log10` | `true` | Log10 y-axis (typical for IF/EdU intensity). |
| `dna_2n_value` | `1000` | Value the G1/2N population maps to; 4N falls at 2×. |
| `x_limits` | `[700, 2250]` | DNA-axis display range (in normalized units). |
| `density_gamma` | `0.35` | Lower = more contrast in the density coloring. |
| `point_size` | `0.3` | Point size in the scatter. |

The EdU baseline anchor defaults to the **G1 median** of the target channel,
matching Figure 1. Set `g1_anchor: "mode"` to force the line through the peak of
the G1 target density instead.

## Public API

- Configuration: `read_facs_config()`, `validate_facs_config()`
- Inputs: `build_sample_manifest()`, `read_facs_sample()`, `validate_facs_inputs()`
- Normalization: `normalize_edu()`, `normalize_poi()`, `normalize_ph3()`
- Analysis: `analyze_facs_experiment()`, `quantify_cell_cycle()`, `quantify_ph3()`
- Plots: `plot_pseudocolor_panels()`, `plot_facs_quantitation()`,
  `plot_ph3_overall()`, `plot_ph3_phase()`, `plot_ph3_diagnostic()`,
  `plot_ph3_boundary_sensitivity()`, `read_ph3_gate_geometry()`,
  `plot_ph3_4n_gate_panels()`, `build_ph3_4n_figure_bundle()`
- Saving: `save_facs_results()`

Low-level density, regression, gate, and layout helpers remain internal.

## Reproducibility notes

- `save_facs_results()` writes only explicitly requested destinations and never
  derives an RDS filename. Set `output_rds` explicitly when the complete analysis
  object should be serialized.
- The 2D KDE is deterministic; the document also sets `set.seed(1)`.
- All inputs and options live in `config.yml`. The rendered HTML echoes the full
  config used and prints `sessionInfo()`, so each report fully documents how its
  figure was made. Commit `config.yml` alongside your results for a complete
  record.

See `docs/CONFIGURATION.md` for the complete configuration contract and
`docs/MIGRATION.md` for changes from the original sourced-script workflow.
