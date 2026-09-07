# EdU reference validation for the Figure 6 CSV-only pilot

**Date:** 2026-09-07
**Status:** IMPLEMENTED, INDEPENDENTLY REVIEWED, AND VERIFIED; Figure 6
real-data report model validated, with HTML rendering pending a Quarto-enabled
environment.

## Purpose

The documented Figure 6 EdU CSV-only pilot explicitly names `Untreated` as a
reference in each of its two biological-replicate runs. The prior canonical
validator rejected every EdU replicate reference even though the package's
reference-normalized output path already requires an explicit configured
reference. The scientific owner approved this narrow correction rather than
removing the documented references from the pilot configuration.

## Change

For EdU only, validation now permits either of these complete configurations:

1. no replicate `reference` values in any biological/technical replicate pair;
   or
2. exactly one configured reference matching a sample label in **every** pair.

Mixed presence, a missing reference label, or more than one matching reference
remain invalid. In particular, the validator distinguishes an absent reference
from a declared label that matches no sample, including when every replicate
declares an invalid label. PH3 remains prohibited from declaring replicate
references; POI continues to require exactly one matched background-control
reference per pair. EdU backgrounds remain independently fit from each
acquisition's EdU-negative events. This update changes no gate, threshold,
normalization, or reference selection: it makes the explicitly configured
reference available to the already-defined reference-normalized output pathway.

The Figure 6 pilot README now accurately distinguishes the external temporary
build/install/HTML-render outputs from its retained, non-written configured
PDF/PNG paths. Those paths resolve under the pilot directory if an explicit
future artifact writer is used; they do not target the historical experiment.

## Files

Modified:

- `R/config.R`
- `tests/testthat/test-config.R`
- `docs/CONFIGURATION.md`
- `test_reports/edu_figure6_nmpp1_csv_pilot_2026-09-07/README.md`

Generated:

- this implementation record.

Experimental inputs read or modified: **none**. The Figure 6 CSV/FCS/WSP files
were not opened. The setup task had previously inspected only the source YAML
and names; this correction did not inspect it again.

Sample mappings, exclusions, gates, thresholds, normalization, statistics, and
biological claims changed: **none**.

## Verification and review

After static diff review, the temporary execution exception was used for:

```sh
Rscript -e 'devtools::test(filter = "config", stop_on_failure = TRUE)'
```

The initial independent scientific review found that an invalid declared label
could be mistaken for a fully absent reference when no group matched it. The
validator and focused synthetic tests were corrected to reject both one-group
and every-group invalid labels. The focused test was rerun:

```sh
Rscript -e 'devtools::test(filter = "config", stop_on_failure = TRUE)'
```

Result: **PASS** — 30 expectations, zero failures, warnings, or skips. The
full package suite then passed **1,753** expectations with zero failures,
warnings, or skips, and `R CMD check --no-manual` returned **Status: OK**.

Independent scientific-integrity re-review: clean. It confirmed that an
explicitly declared invalid reference cannot be mistaken for an absent
reference, that all-absent EdU references remain permitted, and that mixed,
missing, or multiply matched references fail closed. Independent
artifact/security re-review: clean. It found no unsafe I/O, network activity,
telemetry, secrets, fallbacks, or output-containment regression.

The approved Figure 6 CSV-only pilot then read its declared existing CSV
exports and constructed all eight named EdU report panels with status
`available`. It did not modify any experimental input or write into the
historical experiment directory. A known compatibility deprecation warning for
legacy EdU table aliases was emitted; no new warning or fallback was observed.
The sandbox lacks Quarto, so no HTML was rendered there.
