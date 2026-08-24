# Migration from the original Quarto workflow

## What changed

The reusable analysis is now the `facspseudocolor` R package. The active Quarto
report is a thin client that reads YAML, calls package functions, displays
results, and optionally saves explicitly configured outputs.

The scientific calculations for the included EdU and POI examples are protected
by numerical regression tests. Their fitted models, normalized medians, event
counts, and ordering did not change during the package conversion.

## Installation

Install the package before rendering the report:

```bash
R CMD INSTALL .
quarto render pseudocolor_plots.qmd
```

## Paths and example data

Example data moved from `data/` to `inst/extdata/`, the standard location for
files installed with an R package. Paths in file-backed YAML configurations are
resolved relative to the YAML file rather than an implicit working directory.

## Saving and overwrite behavior

The original report silently overwrote configured outputs and derived several
RDS and quantitation filenames. The package no longer invents output filenames.

- PDF and PNG paths remain explicit in YAML.
- Existing files cause an error unless `overwrite = TRUE` or Quarto parameter
  `overwrite:true` is supplied.
- Saving an analysis RDS requires an explicit `output_rds` argument to
  `save_facs_results()`.
- Quantitation plots are returned as ggplot objects and displayed by Quarto;
  save them only to destinations explicitly chosen by the caller.
- EdU canonical CSV tables require an explicit `output_csv_dir`. Deprecated
  aliases require the additional explicit `include_deprecated_csv = TRUE`.

## EdU output schema 2

EdU analyses now add seven scientifically explicit `edu_*` tables and matching
`*_acquisition` tables. The old names retain exactly their former populations:

- `phase_percentages`: historical five-gate assigned composition;
- `phase_medians`: computed-positive Early/Mid/Late-S intensity;
- `whole_medians`: whole computed-positive-population intensity.

These aliases are deprecated for one major-release compatibility window and
emit one warning per analysis. They are never redirected to whole-Single-Cells
or other scientifically different values. Existing saved artifacts are not
rewritten.

## POI peak-detection failure

The old code silently substituted the background-control DNA peak when
per-sample peak detection failed. The new default is:

```yaml
poi_peak_failure: "error"
```

This stops with an actionable error. To deliberately reproduce the former
fallback for a particular analysis, set:

```yaml
poi_peak_failure: "use_background"
```

## One-replicate experiments

Use a `replicates` list containing one entry. In POI mode, name its background
control with `reference`; EdU mode instead fits background from the EdU-negative
population of every acquisition. The old flat `samples` example has been
removed from the recommended workflow.

## Python preprocessing

Python is no longer reachable from the installed package. The unchanged
exporter remains at `python/export_flowjo_populations.py`; optional launching is
implemented in the repository-only `tools/flowjo-orchestration.R`. Enable it in
Quarto only when raw FlowJo/FCS preprocessing is required:

```bash
quarto render pseudocolor_plots.qmd -P run_flowjo_export:true
```
