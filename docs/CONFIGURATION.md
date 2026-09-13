# Configuration reference

## EdU population sources

Use the `gating` block shown below for new EdU configurations. Omitting it
retains historical FlowJo behavior. Legacy paired `g1_source` and
`edu_positive_source` fields remain accepted only for backward compatibility.
Paired `model` values require `g1_model_manifest` and normalize to the
`documented_edu_g1_v1` profile. The package verifies that profile's model-derived
Single Cells, G1, and EdU Positive declarations; the immutable eight-acquisition
Figure 6 mapping; FCS and consumed-snapshot digests; exact panel/channel roles;
artifact and metadata/configuration hashes; feature scopes; thresholds; row
counts; containment; and event-identity digests before normalization. The
default caller reads no workspace and never fits, tunes, or downloads a model.
The only nonfatal model QC codes are `LOW_EVENT_SUPPORT`,
`EXTREME_PREDICTED_FRACTION`, and `ZERO_PREDICTED_POSITIVES`; they are surfaced
in analysis provenance, warnings, and the model-only report without changing a
threshold. Unknown warnings, mapping/digest/runtime/schema mismatches,
incomplete identities, unsupported panels, collisions, and containment failures
stop analysis.

Setting both legacy sources to `flowjo` selects FlowJo compatibility behavior.
Legacy source fields cannot be mixed with a modern `gating` block.

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
also requires `fcs`. EdU background is fitted independently from the
EdU-negative cells in each acquisition. EdU configurations may either omit
replicate `reference` values everywhere or name exactly one explicitly matched
reference sample in every biological/technical replicate pair. Such a reference
supports explicitly requested reference-normalized outputs; it does not supply
the EdU background model. POI requires one matched background-control
`reference` for each biological/technical replicate pair.

PH3 replicates contain `label` and `samples` but no `reference`. PH3 requires
explicit contiguous `g1_x_range`, Early/Mid/Late `s_phase_bins`, and
`g2m_x_range` values.

PH3 also requires an explicit `ph3_input_profile`. Production uses
`production_direct_identity_v1` plus one or more explicit
`ph3_export_operation_dirs`. Each directory must be a finalized Slice 1 export
operation. Relative operation directories are resolved only beneath `data_dir`;
the package does not search for a nearby operation. The historical profile
`legacy_count_only_unverified_v1` is an explicit opt-in, warns once during
analysis, and never claims identity or containment validation. Its structured
rows are keyed by configured prefix and child population; acquisition,
identity, manifest, and hash fields remain explicitly unavailable rather than
being inferred from filenames or row order.

In Slice 2, production PH3 analysis returns validated normalized events and
containment/provenance only; its quantitation list remains empty until the
approved production metric slices are implemented. The explicit legacy profile
retains historical pH3 quantitation and conditional `event_index` safety checks,
without reporting identity containment as validated.

There are no package defaults for these scientific phase boundaries. The
values in `examples/config_ph3.yml` are visible starting examples only and are
never inserted into a configuration automatically.

`ph3_positivity_method` is explicit. `flowjo_legacy_v1` preserves the
historical FlowJo-positive membership. `ph3_raw_4n_density_cutoff_v1` uses the
matched untreated control's retained 4N raw-pH3 density, selecting the first
right-side local minimum after its unique dominant peak; the frozen cutoff is
then applied to every matched treatment across DNA regions. The FlowJo
pH3-positive export remains provenance during the direct-identity transition,
not the numerical authority for the computed method. Its primary phase calls
use the configured DNA ranges, while `ph3_boundary_sensitivity_fraction`
produces a diagnostic showing how the G2/M percentage changes when both G2/M
boundaries move by that fraction of the configured 2N value.

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

## Optional EdU report block

An EdU configuration may preserve settings needed by the optimized HTML report:

```yaml
report:
  all_events_dir: all_events_csv
  dna_height_channel: FL2-H
  show_apex_comparison_toggle: true
  show_phase_gate_toggle: true
  embed_pseudocolor_pdf_downloads: true
```

`all_events_dir` is resolved relative to the configuration file. The report
uses `dna_height_channel` only for the DNA-H-versus-DNA-A all-events view. These
settings do not create or validate FlowJo exports.

## Output policy

The package never invents missing output destinations and refuses to overwrite
existing files unless explicitly authorized. `output_pdf` and `output_png` are
used by the Quarto front end. An analysis RDS requires an explicit `output_rds`
argument to `save_facs_results()`.
# EdU gating source

EdU configurations select one explicit, fail-closed gating source:

```yaml
gating:
  mode: "flowjo"
```

or the owner-approved experimental FCS-only profile:

```yaml
gating:
  mode: "model_experimental"
  profile: "single_g1_edu_frozen_v1"
  status_label: "EXPERIMENTAL MODEL-DERIVED NON-PRODUCTION"
  output_dir: "/absolute/new/model-output-directory"
  artifacts: # exact paths and frozen hashes for Single Cells, G1, and EdU+
  acquisitions: # exact FCS paths, hashes, identities, and channel roles
```

Model mode uses no FlowJo workspace and never falls back to FlowJo or existing
CSVs. It calls Single Cells, G1, and EdU+ at frozen thresholds 0.205, 0.255, and
0.385. For EdU experiments, the final G1 population excludes every model-called
EdU+ event before the G1 table is used to calculate the normalization anchor or
background. It then calculates the median positive raw DNA-A value among those
preliminary G1 candidates and requires final G1 events to have DNA-A at least
35% of that acquisition-specific center. The operation records the center,
minimum, and excluded-event count in its unified hash/count/containment manifest. The output directory
must not already exist. This profile is **EXPERIMENTAL MODEL-DERIVED
NON-PRODUCTION**: the G1 and EdU models were originally evaluated within
expert-gated Single Cells, so their composition after model-derived Single
Cells is not validated biological ground truth.

The EdU artifact is the exact SHA-256-pinned `PORTABLE_HGB_V1` JSON export of
`EDU_POSITIVE_MODEL002`. It retains the original model, ordered features, and
0.385 threshold while removing pickle loading and the exact Python 3.10.12 /
scikit-learn 1.7.2 runtime requirement. NumPy, pandas, SciPy, FlowKit, and
PyYAML remain pipeline dependencies.
The generated manifest records `derived_from_reference_sha256` for the external
validated predictor source; that digest describes the implementation lineage
and is not presented as a hash of the embedded predictor code.

The separately preserved documented profile is selected as follows:

```yaml
gating:
  mode: "model_experimental"
  profile: "documented_edu_g1_v1"
  status_label: "EXPERIMENTAL MODEL-DERIVED NON-PRODUCTION"
  manifest: "/absolute/path/to/model-gating-manifest.json"
```

It validates the provenance-bound documented model operation and its 0.29 G1
threshold. It does not apply the frozen profile's 0.255 G1 threshold or 35%
acquisition-relative DNA minimum. The two profiles are intentionally distinct.
