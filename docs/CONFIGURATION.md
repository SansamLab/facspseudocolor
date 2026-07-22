# Configuration reference

`read_facs_config()` reads one explicit YAML file, rejects unknown top-level
keys, applies centralized defaults, and validates the complete structure before
experimental data are processed. Relative paths are resolved from the YAML
file's directory.

## Required settings

| Setting | Type | Meaning |
|---|---|---|
| `plot_type` | `edu` or `poi` | Selects the scientific normalization method. |
| `data_dir` | path | Directory containing event-level CSV files. |
| `dna_channel` | string | Exact numeric DNA column name. |
| `target_channel` | string | Exact numeric EdU or POI column name. |
| `target_name` | string | Human-readable target label. |
| `output_pdf` | path | Explicit main-figure PDF destination. |
| `output_png` | path | Explicit main-figure PNG destination. |
| `replicates` | list | Replicates, references, conditions, and unique prefixes. |

Each replicate must contain `label`, `reference`, and `samples`. Every sample
must contain `label` and `prefix`; optional FlowJo orchestration also requires
`fcs`. Replicates must contain the same conditions in the same order, and each
`reference` must match exactly one condition.

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
```

Changing a suffix changes the exact required filename; the package does not
search for alternatives.

## Scientific normalization

| Setting | Default | Applies to |
|---|---:|---|
| `dna_2n_value` | `1000` | Both modes |
| `g1_anchor` | `median` | EdU |
| `baseline_fit_x_range` | `[1000, 2000]` | EdU |
| `baseline_boundary_bins` | `20` | EdU |
| `baseline_minimum_events_per_bin` | `20` | EdU |
| `baseline_minimum_negative_events` | `100` | EdU |
| `background_quantile` | `0.95` | POI |
| `poi_dna_align` | `per_sample` | POI |
| `poi_peak_failure` | `error` | POI |

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

## Display and density settings

Presentation defaults include `palette: refined`, `layout: plotgardener`,
`y_log10: true`, `x_limits: [700, 2250]`, and quantile-based shared y limits.
See the fully commented `config.yml` for density clipping, point size, layout,
gate styling, and quantitation styling options.

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
