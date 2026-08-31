# pH3 output contract Slice 3 — signal outcomes

**Date:** 2026-08-31
**Status:** IMPLEMENTED — corrected after focused local verification failure; local re-verification required
**Canonical source:** `facs_pseudocolor_workflow` on branch
`feature/ph3-output-contract-signal-outcomes`, based on `main` at
`174bf4f`.

## Purpose and bounded scope

This slice implements only the owner-confirmed quantitative signal outcomes:

- **C:** median pH3 signal among pH3-positive 4N events;
- **D:** median pH3 signal among pH3-positive below-4N events.

It consumes the validated Slice 1 output model and Slice 2 event-specific
background result. It makes no new membership, threshold, gating, DNA-region,
normalization, or reference-selection decisions.

For every technical acquisition and each of C/D, the implementation retains one
acquisition median. A biological sample's direct C/D value is the equal-weight
arithmetic mean of its available technical-acquisition medians. An acquisition
with no qualifying pH3-positive events, or with a nonfinite required signal, is
explicitly unavailable. It is not silently dropped: the sample has a partial
coverage status when other acquisition medians remain, and is unavailable when
none remain.

The selected signal basis is required to be one exact analysis-wide value:
`individual_corrected`, `pooled_corrected`, or `raw`. A mixed basis fails before
calculations occur. The retained Slice 2 correction and reference records are
included as `basis_qc` and `reference_qc` beside the new value tables.

When a reference condition is configured, the Slice 1 pre-resolved reference
sample is used exactly once per biological-replicate set. Ratios are population
specific. A missing, unavailable, zero, or nonfinite reference value makes only
the affected ratio unavailable; no substitute, pseudocount, or alternate
reference is used. With no configured reference, ratios are `not_applicable`.

## Explicitly excluded

This slice does **not** implement condition summaries, paired statistics,
figures, pseudocolor display offsets, Quarto, CSV/RDS exports, report writing,
geometry, or any EdU/POI behavior. Existing Slice 1/2 tables, event identities,
pH3 positivity membership, 4N/below-4N membership, and prevalence outputs are
retained unchanged.

## Files read

- `AGENTS.md`
- owner-confirmed pH3 output contract and decision-recovery records under
  `/Users/sansamc/Documents/Codex/2026-08-30/`
- `R/analysis.R`, `R/ph3_outputs.R`, and `R/ph3_background.R`
- focused pH3 output-model and background-regression tests.

Experimental `.fcs`, `.wsp`, `.prism`, CSV/RDS, historical reports, and
generated figures were not read or modified.

## Files modified and generated

Modified:

- `R/analysis.R`
- `R/ph3_outputs.R` (retains the immutable sample-to-condition mapping already
  established by the validated output model)

Generated source/test/documentation:

- `R/ph3_signal.R`
- `tests/testthat/test-ph3-signal-outcomes.R`
- this implementation record.

The pre-existing untracked output-contract audit documents are outside this
slice and remain untouched.

## Scientific-change declaration

Sample mappings: none. Exclusions: none. Gates: none. Thresholds: none. DNA
normalization: none. Background fitting method: none; the existing Slice 2
result is consumed. Statistics: none. Biological claims: none.

## Tests prepared, not run

The new tests use only explicitly labelled `SYNTHETIC` in-memory values. They
cover equal technical-acquisition weighting, population-specific C/D values,
optional reference behavior, zero-reference unavailability, partial technical
coverage, unchanged source event signals, and rejection of mixed bases.

Codex did not execute R, Python, package loading, tests, checks, builds,
renders, analyses, or benchmarks. Verification is **NOT RUN — local execution
required**.

## Local verification correction

The user-reported first focused run of `ph3-signal-outcomes` failed in all five
nominal SYNTHETIC cases before outcome calculation with
`sample_provenance_mismatch`. The fixture's event and configured sample values
were the same. The validator compared those keyed tables with `identical()`
after `unique()` but without normalizing their row names; `unique()` retained
the first event-row names while the configured table retained its own row
names. This was an implementation defect in the provenance comparison, not a
scientific mismatch or a fixture defect.

The correction sorts both keyed sample-provenance comparisons by `sample_id`
and clears only their row names before exact comparison. This applies both to
the `(replicate_set_id, sample_id)` check and the subsequent
`(experiment_id, replicate_set_id, sample_id)` check, which otherwise retained
the same event-derived row-name difference. The strict required values, keyed
sample/replicate-set matching, duplicate-sample rejection, acquisition mapping
check, subsequent experiment/replicate-set provenance check, and immutable
condition mapping check are unchanged. No gates, memberships, signals,
background basis, reference rule, or scientific method changed.

Codex did not execute the corrected project code or any tests. The corrected
tree is **NOT RUN — local execution required**.

## Independent review and resolution

The independent scientific-integrity review found two important issues in the
initial implementation: an unavailable technical acquisition in the reference
sample incorrectly made its otherwise valid partial reference value
unavailable, and a tampered configured reference record could degrade to a
local unavailable ratio rather than fail as an invalid configuration. The
implementation now derives each reference value from its already calculated
sample-level direct value, propagates partial coverage to the ratio, and
validates configured reference sample, condition, and replicate-set identity
against the authoritative sample table before calculation.

The independent artifact/security review found three important fail-closed
validation gaps: incomplete experiment/condition provenance validation,
unvalidated event-signal status values, and the same malformed-reference issue.
The implementation now requires exact experiment/sample/condition/replicate-set
mapping against the immutable source mapping retained in the Slice 1 model,
validates the authoritative event status/reason combinations, and rejects
malformed reference records. Focused SYNTHETIC negative tests cover those
conditions. Both reviewers made no edits and executed no project code.

For the local-verification correction above, the renewed independent
scientific-integrity review was clean. The renewed artifact/security review
initially identified the parallel row-name issue in the detailed
experiment/replicate-set comparison; it was corrected with the same
value-preserving normalization. A final artifact/security re-review is required
before human re-verification. Neither reviewer executed project code.

From the canonical repository, with R >= 4.2 and the declared dependencies,
run:

```sh
Rscript -e 'devtools::test(filter = "ph3-signal-outcomes", stop_on_failure = TRUE)'
Rscript -e 'devtools::test(filter = "ph3-output-model|ph3-background-regression", stop_on_failure = TRUE)'
Rscript -e 'devtools::test(filter = "^ph3$", stop_on_failure = TRUE)'
Rscript -e 'devtools::test(stop_on_failure = TRUE)'
git diff --check
git status --short --untracked-files=all
```

Expected success: every test command has zero failures/errors; `git diff
--check` has no output; and the only new Slice 3 paths are the three listed
above plus the modified `R/analysis.R` and `R/ph3_outputs.R`. Stop on a changed event identity or
membership, any mixed analytical basis, any silent unavailable acquisition,
reference substitution, unexpected warning, or unexpected file. Report only
the command, PASS/FAIL, failed test names, short relevant excerpt, and
unexpected files.

## Assumptions and remaining uncertainty

Slice 1's validated config/reference mapping and Slice 2's analysis-wide basis
decision remain authoritative. Condition-level summaries and paired tests are
intentionally deferred; therefore this slice does not combine C/D values across
conditions or basis strata.
