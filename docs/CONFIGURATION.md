# Configuration reference

`read_facs_config()` reads one explicit YAML file, rejects unknown top-level
keys, applies centralized defaults, and validates the complete structure before
experimental data are processed. Relative paths are resolved from the YAML
file's directory.

## Required settings

| Setting | Type | Meaning |
|---|---|---|
| `plot_type` | `edu`, `poi`, or `ph3` | Selects the scientific analysis method. |
| `data_dir` | path | Directory containing event-level CSV files. |
| `dna_channel` | string | Exact numeric DNA column name. |
| `target_channel` | string | Exact numeric EdU, POI, or pH3 column name. |
| `target_name` | string | Human-readable target label. |
| `output_pdf` | path | Explicit main-figure PDF destination. |
| `output_png` | path | Explicit main-figure PNG destination. |
| `replicates` | list | Biological/technical acquisitions, conditions, and unique prefixes. |

Every sample must contain `label` and `prefix`; optional FlowJo orchestration
also requires `fcs`. EdU acquisitions do not have a `reference`: their
background is fitted independently from the EdU-negative cells in each
acquisition. POI requires one matched background-control `reference` for each
biological/technical replicate pair.

PH3 replicates contain `label` and `samples` but no `reference`. PH3 requires
explicit contiguous `g1_x_range`, Early/Mid/Late `s_phase_bins`, and
`g2m_x_range` values.

There are no package defaults for these scientific phase boundaries. The
values in `examples/config_ph3.yml` are visible starting examples only and are
never inserted into a configuration automatically.

The pH3-positive population must be drawn by the user in FlowJo. The R package
does not infer a pH3 cutoff. Its primary phase calls use the configured DNA
ranges, while `ph3_boundary_sensitivity_fraction` produces a diagnostic showing
how the G2/M percentage changes when both G2/M boundaries move by that fraction
of the configured 2N value.

## Input suffixes

The `suffixes` map defaults by mode:

```yaml
# EdU
suffixes:
  complete: "_single_cells.csv"
  g1: "_g1.csv"
  edu_positive: "_edu_positive.csv"

# POI
suffixes:
  complete: "_single_cells.csv"

# PH3
suffixes:
  complete: "_single_cells.csv"
  g1: "_g1.csv"
  ph3_positive: "_ph3_positive.csv"
```

Changing a suffix changes the exact required filename; the package does not
search for alternatives.

## Scientific normalization

| Setting | Default | Applies to |
|---|---:|---|
| `dna_2n_value` | `1000` | All modes |
| `g1_anchor` | `median` | EdU and PH3 |
| `baseline_fit_x_range` | `[1000, 2000]` | EdU |
| `baseline_boundary_bins` | `20` | EdU |
| `baseline_minimum_events_per_bin` | `20` | EdU |
| `baseline_minimum_negative_events` | `100` | EdU |
| `background_quantile` | `0.95` | POI |
| `poi_dna_align` | `per_sample` | POI |
| `poi_peak_failure` | `error` | POI |
| `ph3_boundary_sensitivity_fraction` | `0.05` | PH3 diagnostic only |

`poi_peak_failure: use_background` is an explicit request for the legacy
background-peak fallback. It is never selected silently.

## Phase gates and quantitation

| Setting | Default | Meaning |
|---|---:|---|
| `show_phase_gates` | `false` | Overlay configured gates on panels. |
| `show_gate_assignment` | `false` | Show deterministic gate-assignment diagnostics. |
| `quantify_phase_median` | `false` | Legacy report-display switch; values are always calculated. |
| `quantify_whole_median` | `false` | Legacy report-display switch; values are always calculated. |
| `quantify_phase_percent` | `false` | Legacy report-display switch; values are always calculated. |
| `quant_signal` | `background_subtracted` | Signal used for median quantitation. |
| `quantify_reference_normalized` | `false` | Add within-replicate reference ratios. |
| `quant_reference_condition` | none | Explicit displayed ratio reference. |
| `quant_error_bar` | `sd` | `sd`, `sem`, or `none`. |

Custom geometry may be supplied through `s_phase_bins`, `g1_x_range`,
`g2m_x_range`, `negative_y_range`, and `s_phase_y_range`. Ranges must contain
two increasing finite numbers.

`quant_reference_condition` is not a background control. It is an optional
experimental condition used only after background subtraction, so reported
medians can be expressed relative to that condition. Technical acquisitions
are processed separately before their summaries are averaged.

## Display and density settings

Presentation defaults include `palette: refined`, `layout: plotgardener`,
`y_log10: true`, `x_limits: [700, 2250]`, and quantile-based shared y limits.
See the fully commented `config.yml` for density clipping, point size, layout,
gate styling, and quantitation styling options.

`pseudocolor_signal` selects the event-level y value used in pseudocolor
panels. The default, `background_subtracted`, displays
`target_raw - baseline` plus an offset. Set it to `normalized` only to request
the legacy baseline-divided display. `background_subtracted_offset` may be a nonnegative number or `auto`
(the default). Automatic mode takes the median positive raw G1 target anchor
across displayed samples and rounds it to the nearest power of ten; for
example, a typical raw background near 6,000 produces an offset of 10,000.
The offset is display-only and never changes background-subtracted
quantitation. Optional reference normalization is applied afterward to the
background-subtracted medians, within each biological replicate.

## Optional FlowJo block

The package accepts but never executes `flowjo`. The repository-only launcher
uses `source_dir`, `workspace`, `python`, `dna_source_channel`,
`target_source_channel`, `populations`, and `rebuild`. See
`PYTHON_INTERFACE.md`.

## Output policy

The package never invents missing output destinations and refuses to overwrite
existing files unless explicitly authorized. `output_pdf` and `output_png` are
used by the Quarto front end. An analysis RDS requires an explicit `output_rds`
argument to `save_facs_results()`.
