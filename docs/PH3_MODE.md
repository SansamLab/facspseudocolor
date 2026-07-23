# PH3 mode

PH3 mode measures a user-defined FlowJo pH3-positive population. The R package
does not estimate, regress, move, or reinterpret the pH3 positivity boundary.

## Required populations

Each configured sample requires:

```text
<prefix>_single_cells.csv
<prefix>_g1.csv
<prefix>_ph3_positive.csv
```

- `Single Cells` supplies the denominator.
- `G1` supplies the per-sample DNA normalization anchor.
- `pH3 Positive` is a user-drawn FlowJo polygon nested within Single Cells.

An empty pH3-positive CSV is valid and produces zero percent positive. When
`event_index` is available, validation confirms that pH3-positive events are a
subset of Single Cells. Without event indices, the package still checks that
the positive count does not exceed the Single Cell count.

## DNA normalization

For each sample:

```text
dna_norm = raw DNA / G1 DNA anchor * dna_2n_value
```

The G1 anchor is the configured `median` or `mode` of the sample's FlowJo G1
population. Target signal is retained as supplied; PH3 mode performs no target
normalization or background-threshold fitting.

## Explicit phase boundaries

PH3 mode requires explicit, contiguous ranges for G1, Early S, Mid S, Late S,
and G2/M. No scientific phase thresholds are silently defaulted. The final
G2/M upper boundary is inclusive; other upper boundaries are exclusive.

The ranges shown in the example YAML are starting examples, not package
defaults. They must be reviewed against the normalized-DNA distributions for
the experiment.

Events outside the configured ranges or with nonfinite normalized DNA are
reported as `Unassigned`.

## Denominators

Overall positivity:

```text
100 * pH3-positive events / all Single Cell events
```

Each phase:

```text
100 * pH3-positive events assigned to the phase / all Single Cell events
```

The five phase percentages plus the Unassigned percentage equal the overall
pH3-positive percentage, apart from numerical rounding.

## Outputs

The analysis object stores:

- Single Cell, pH3-positive, assigned, and unassigned counts.
- Overall pH3-positive percentage.
- Phase-specific counts and percentages.
- Explicit phase-gate definitions.
- G2/M boundary-sensitivity diagnostics.
- Per-sample G1 anchors and DNA normalization factors.
- Input paths, configuration, and software provenance.

Editable plots include overall positivity, grouped and stacked phase
percentages, pH3/DNA assignment diagnostics, and G2/M boundary sensitivity.

## G2/M clearance diagnostic

`ph3_boundary_sensitivity_fraction` defines a diagnostic shift expressed as a
fraction of `dna_2n_value`. At the default `0.05` with 2N at `1000`, the lower
and upper G2/M boundaries are moved by 50 normalized DNA units. These
alternative counts are diagnostic only; the configured boundaries remain the
reported scientific result.
