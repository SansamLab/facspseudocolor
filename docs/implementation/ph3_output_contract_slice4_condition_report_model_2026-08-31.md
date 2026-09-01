# pH3 output contract Slice 4 — condition report model

**Date:** 2026-08-31
**Status:** IMPLEMENTED — corrected local verification required
**Canonical source:** `facs_pseudocolor_workflow`, branch
`feature/ph3-output-contract-signal-outcomes`, based on Slice 3 commit
`2ff1622`.

## Purpose and bounded scope

This slice creates a pure, validated, report-ready condition model for the
four owner-confirmed pH3 outcomes:

- A: 4N pH3-positive prevalence within Analysis singlets (2N–4N);
- B: below-4N pH3-positive prevalence within Analysis singlets (2N–4N);
- C: pH3 signal in 4N pH3-positive cells; and
- D: pH3 signal in below-4N pH3-positive cells.

The model retains one explicit biological-replicate value per condition and
outcome, condition summaries, signal-basis strata for C/D, and concise
per-condition/outcome QC flags. It consumes the authoritative Slice 1–3 model
only. It does not recalculate event membership, prevalence, background fits,
or acquisition/sample signal medians.

For A/B, the model consumes the existing validated equal-weight
technical-acquisition biological-replicate values. For C/D, it consumes the
Slice 3 sample-level value: the reference ratio when a reference is configured,
otherwise the direct median. It records partial and unavailable values without
substitution.

Condition means give each finite biological-replicate value equal weight. SD is
recorded only for three or more finite biological replicates, matching the
owner-confirmed display rule. C/D values on more than one signal basis never
receive a combined condition mean; their overall condition row is explicitly
`unavailable_mixed_signal_basis` and the individual basis strata remain
available for a later plotting/report slice.

## Excluded work

This slice does not implement plotting, Quarto rendering, display offsets or
positive-domain safeguards, CSV/RDS/JSON writing, geometry generation,
inferential statistics, comparisons, or EdU/POI work. It does not alter the
approved pH3-positive, 4N, below-4N, or eligible-population definitions.

## Files read

- workspace `AGENTS.md`;
- owner-confirmed pH3 output/QC/report contract and implementation-decision
  recovery record (external, nonportable scientific-owner consultation
  references; no local path is required by this repository);
- `R/analysis.R`, `R/ph3_outputs.R`, `R/ph3_background.R`,
  `R/ph3_signal.R`, and the prior Slice 1–3 implementation records;
- focused synthetic pH3 output-model and signal-outcome tests.

Experimental `.fcs`, `.wsp`, `.prism`, CSV/RDS, historical reports, and
generated figures were not read or modified.

## Files modified and generated

Modified:

- `R/analysis.R` — attaches the completed condition report model only after
  the existing validated Slice 1–3 sequence succeeds.
- `R/ph3_background.R` and `R/ph3_signal.R` — restore the owner-confirmed
  per-biological-replicate-set fallback after independent review found the
  inherited global fallback; this is valid individual fit, otherwise valid
  pooled set fit, otherwise raw signal for that set only.
- `tests/testthat/test-ph3-background-regression.R` — updates the existing
  SYNTHETIC assertions from the superseded global-fallback expectation to the
  owner-confirmed per-set behavior.

Generated source/test/documentation:

- `R/ph3_report_model.R`;
- `tests/testthat/test-ph3-condition-report-model.R`;
- this record.

The three pre-existing untracked audit documents remain outside this slice and
untouched.

## Scientific-change declaration

Sample mappings: none. Exclusions: none. Gates: none. Thresholds: none. DNA
normalization: none. Positivity and DNA-region membership: none. Statistics:
none. Biological claims: none. New behavior is deterministic packaging and
validation of already approved values for a later report/plot layer. The
inherited global fallback was corrected to the already approved
per-replicate-set method; no new background method was introduced.

## Tests prepared, not run

Focused tests use explicitly labelled `SYNTHETIC` in-memory values only. They
cover the exact A–D ordering, reference-ratio versus direct-median mode,
equal biological-replicate accounting, SD suppression below three finite
replicates, explicit partial prevalence state, mixed-basis no-combined-summary
with retained strata, and fail-closed rejection of tampered signal-basis or
prevalence provenance. Existing background tests also cover per-set pooled and
raw fallbacks alongside unaffected sets.

Codex did not execute R, Python, package loading, tests, checks, builds,
renders, analyses, or benchmarks. Verification is **NOT RUN — local execution
required**.

## User-reported verification correction

The initial focused local run of `ph3-condition-report-model` reported five
test failures. Static review determined that all five were stale SYNTHETIC test
expectations rather than production-model failures:

- the expected biological-replicate table had eight rows, but the approved
  schema requires one row per condition, biological-replicate set, and outcome
  (2 conditions × 2 sets × 4 outcomes = 16 rows);
- the mixed-basis fixture targeted the nonexistent
  `SYNTHETIC-replicate-set-2`, while the fixture's configured second set is
  `SYNTHETIC-set-2`. It therefore left all C/D values on one basis and could
  not exercise the required no-combined-summary behavior.

The correction changes only `tests/testthat/test-ph3-condition-report-model.R`:
it requires the 16-row condition×set×outcome table and targets the configured
second SYNTHETIC replicate set. Production R source is unchanged. The test
continues to require C/D condition summaries to be unavailable with `NA` means
when bases differ, while retaining the eight condition×outcome×basis strata.
Codex did not run the corrected test; local verification remains required.

Independent read-only scientific-integrity re-review was clean: the corrected
row count and ordering match the one-value-per-condition×replicate-set×outcome
contract, and the corrected fixture now exercises the required mixed-basis
prohibition without changing any scientific method. Independent read-only
artifact/security re-review was also clean: no unsafe I/O, network activity,
secrets, generated artifacts, hidden fallback, or machine-path dependency was
introduced. Neither reviewer ran project code or accessed experimental inputs.

## User-reported background-regression correction

The subsequent affected-regression command reported six errors in
`test-ph3-background-regression.R`. Each stopped at
`ph3_background_fit_row()` while materializing pooled-fit QC, with
`data.frame()` reporting a mix of scalar and zero-length fields. The failure
was not a fit-validity, coverage, membership, or outcome-calculation result.

Static diagnosis found that the new per-replicate-set decision code assigned
pooled fits inside an `lapply()` callback. That callback-local assignment did
not populate the surrounding `pooled_fits` list used later to construct the
pooled QC table. Consequently, the QC-row helper received an empty fit entry.

The minimal production correction replaces that callback with an explicit,
deterministic loop over the validated replicate-set order. Every iteration now
stores exactly one complete fit object before the decision and QC tables are
constructed: a valid or invalid pooled fit after an individual failure, or the
existing explicit scalar-`NA` `not_attempted` fit when all individual fits are
valid. It does not change individual-fit rules, the per-set
individual→pooled→raw selection policy, fit inputs, coverage policy, event
signals, counts, provenance, or report-model logic. Tests were not weakened or
modified for this correction.

Codex did not run the corrected code or any project tests. The corrected tree
is **NOT RUN — local execution required**.

Independent read-only scientific-integrity review initially flagged the
zero-length QC construction as blocking. After the explicit-loop correction,
re-review was clean: every replicate-set entry is populated before QC-row
construction, and the approved individual→pooled→raw policy remains confined
to that replicate set. Independent read-only artifact/security review was
clean with no blocking, important, or advisory findings. Neither reviewer ran
project code or accessed experimental inputs.

## Independent review and resolution

An independent scientific-integrity reviewer initially found one blocking
inherited issue: the existing background fallback was global rather than
per-replicate-set, which prevented the required mixed-basis signal reporting.
It also found that the initial A/B conversion obscured an explicit partial
status. Both were corrected without weakening validation: fallback decisions
and event/signal validation are now set-scoped; report values retain a stable
partial status/reason; and focused SYNTHETIC tests prove mixed-basis condition
summaries remain unavailable while their basis strata remain available. A final
scientific re-review was clean. No fabricated data,
changed sample mapping, gate, threshold, normalization, statistic, or
biological claim was introduced.

An independent artifact/security reviewer found one advisory: this record
initially embedded a machine-specific absolute path for external owner records.
The reference is now explicitly described as an external, nonportable
consultation reference; no runtime code uses that location. The final artifact
review was otherwise clean: no network/telemetry, secrets, shell execution,
unsafe path, generated artifact, hidden fallback, or destructive I/O was
introduced. Reviewers made no edits and executed no project code.

From the canonical repository, with R >= 4.2 and declared dependencies, run:

```sh
Rscript -e 'devtools::test(filter = "ph3-condition-report-model", stop_on_failure = TRUE)'
Rscript -e 'devtools::test(filter = "ph3-signal-outcomes|ph3-background-regression|ph3-output-model", stop_on_failure = TRUE)'
Rscript -e 'devtools::test(filter = "^ph3$", stop_on_failure = TRUE)'
Rscript -e 'devtools::test(stop_on_failure = TRUE)'
git diff --check
git status --short --untracked-files=all
```

Expected success: every test command has zero failures/errors; `git diff
--check` has no output; no generated repository artifacts appear; and only the
Slice 4 paths above plus modified `R/analysis.R`, `R/ph3_background.R`,
`R/ph3_signal.R`, and `tests/testthat/test-ph3-background-regression.R`, as
well as the three existing audit documents, are present. Stop on any changed event membership, prevalence,
signal median, reference value, signal basis, unexpected warning, test failure,
or unexpected file. Report only command, PASS/FAIL, failed test names, a short
relevant excerpt, and unexpected files.

## Assumptions and unresolved uncertainty

The merged Slice 1–3 output model, background basis, and signal outcome
schemas are authoritative. Signal-basis strata are preserved so a later
plotting/report slice can show them without manufacturing a cross-basis
summary. Statistics, file exports, display transforms, and geometry remain
intentionally deferred.
