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

Use a `replicates` list containing one entry and name its `reference`. The old
flat `samples` example could not encode the reference required by EdU or POI
normalization and has been removed from the recommended workflow.

## Python preprocessing

Python is no longer reachable from the installed package. The unchanged
exporter remains at `python/export_flowjo_populations.py`; optional launching is
implemented in the repository-only `tools/flowjo-orchestration.R`. Enable it in
Quarto only when raw FlowJo/FCS preprocessing is required:

```bash
quarto render pseudocolor_plots.qmd -P run_flowjo_export:true
```
