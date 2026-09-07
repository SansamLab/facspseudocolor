# EdU report Slice 1 implementation record

**Date:** 2026-09-07
**Status:** IMPLEMENTED, INDEPENDENTLY REVIEWED, AND VERIFIED (except Quarto render availability)

## Scope

Implements only the owner-authorized every-sample DNA-versus-EdU pseudocolor
report model from `edu_report_slice1_readiness_2026-09-07.md`.

For every completed EdU acquisition, the background-fit model now retains the
exact inferred-negative row membership and a SHA-256 fingerprint of the
corresponding DNA/raw-signal values. The display report uses that proof to
calculate `median(raw negative) - median(corrected negative)` for that exact
sample only. It applies the result only to a temporary display coordinate.

## New surface

- `build_edu_pseudocolor_output_contract()` returns one named ggplot per
  manifest prefix plus display-offset and panel-QC tables.
- `facs_edu_pseudocolor_output_contract.qmd` is a thin Quarto front end for a
  completed analysis or an explicit configuration.

If the model identity, retained membership proof, raw/corrected values, offset,
or display coordinate cannot be validated, that panel is explicitly suppressed
with its reason. No raw, global, rounded, pooled, or substitute offset is used.

## Deliberately unchanged

No changes were made to sample mapping, event inclusion, background fitting,
positivity, DNA normalization, phase regions, thresholds, corrected medians,
reference normalization, statistics, canonical EdU tables, or quantitative
exports. This slice does not add density summaries, aggregate figures, report
manifests, or new QC criteria.

## Inputs and artifacts

Experimental inputs read or modified: none. Generated experimental outputs:
none. Synthetic tests use only unmistakably SYNTHETIC in-memory events.

## Verification and reviews

Executed under the temporary owner-approved execution exception:

```sh
Rscript -e 'devtools::test(filter = "edu-pseudocolor-output-contract", stop_on_failure = TRUE)'
Rscript -e 'devtools::test(filter = "edu-output-contract", stop_on_failure = TRUE)'
Rscript -e 'devtools::test(stop_on_failure = TRUE)'
```

Results: focused Slice 1 test PASS (14 assertions); existing EdU output-contract
test PASS (89 assertions); full package suite PASS (1,747 assertions; no
failures, warnings, or skips). No experimental input was read.

An external temporary-directory synthetic render was prepared using only the
SYNTHETIC test fixture. Package build and install passed, but rendering could
not start because the current execution environment has no `quarto` executable
on its PATH. No HTML report or repository artifact was produced. This is an
environment limitation, not a package/render error; perform the synthetic
render on a workstation with Quarto installed before accepting a real-data
render.

Independent scientific-integrity review: clean. It confirmed that each offset
uses exact retained fit membership, cross-sample reassignment fails visibly,
and numerical science is unchanged. Independent artifact/security review:
clean. It found no network, telemetry, unsafe I/O, hard-coded experimental
paths, generated repository artifacts, or fallback path. Both reviewers noted
that concurrent pH3 changes and `.claude/` must remain outside an EdU-only
commit.
