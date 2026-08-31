# pH3 output contract Slice 2 implementation handoff

**Date:** 2026-08-31
**Status:** CORRECTED — FOCUSED TEST USER-REPORTED PASS
**Scope:** Audit Slice 2 only — regression validity and fallback engine

## Canonical source and baseline

Canonical source directory:

`/Users/sansamc/OMRF Dropbox/Chris Sansam/Synchd/Organized Dropbox/Projects and Data/Projects/Cowork_Projects/2026_Analyze_FACS_Data/facs_pseudocolor_workflow`

Branch at implementation start: `feature/ph3-output-contract-audit`
Starting `HEAD`: `2f5f709fdf4871c3617e4754e2d535eae26c0236`
Starting description: `Add pH3 output contract model`

The following pre-existing untracked audit records were protected and were not
modified or treated as implementation outputs:

- `docs/implementation/edu_output_contract_gap_audit_2026-08-30.md`
- `docs/implementation/ph3_output_contract_gap_audit_2026-08-30.md`
- `docs/implementation/ph3_slice6_output_contract_audit_2026-08-29.md`

## Authorized boundary implemented

Slice 2 consumes one validated Slice 1 `ph3_output_model` and adds an in-memory
event-specific pH3 background-correction result. For each sample, it fits raw
pH3 signal versus raw DNA among events with
`eligible_2to4n & !ph3_positive_member`. It requires at least 100 negative
events, finite fit inputs, finite nonzero DNA variation, a full-rank
intercept-plus-slope design, finite fitted coefficients and required
predictions, and negative-event DNA coverage spanning the inclusive eligible
pH3-positive DNA range. No event is silently removed to make a model valid.

Tier selection follows the later owner-confirmed experiment-wide decision in
the recovery record:

1. Use individual-corrected signal throughout the experiment only when every
   sample has a valid individual fit.
2. If any individual fit fails, fit a separate pooled-negative model for every
   biological-replicate set. Use pooled-corrected signal throughout the
   experiment only when every pooled fit is valid.
3. If any pooled fit fails, use raw pH3 signal throughout the experiment.

The result retains every individual and pooled diagnostic, separate truthful
individual-trigger and pooled-failure reasons, one common selected basis,
per-set decision rows, and one event-signal row per source event. Nonfinite
analytical event signal is retained as explicitly unavailable rather than
omitted or substituted. A pooled fit that is not attempted has explicit
`not_attempted` status and typed missing diagnostic values rather than
plausible-looking zero counts.

The Slice 1 source classifications and quantitation tables are not mutated.
The engine validates their schemas and exact experiment, replicate-set, sample,
acquisition, prefix, event, and membership mappings before fitting.

## Owner-confirmed decisions applied

- Raw DNA and raw pH3 are the regression inputs.
- Eligible pH3-negative events are the only fit observations.
- Coverage equality at either DNA-range boundary is valid; any positive-range
  gap is invalid because it would require extrapolation.
- An empty eligible positive population makes coverage not applicable; it does
  not manufacture a C/D value or invalidate an otherwise valid fit.
- With one intercept and one finite DNA predictor, exact design-rank deficiency
  implies zero DNA variation. The ordered validator reports the more specific
  `zero_dna_variation` reason first and retains a separate numerical rank guard.
- Nonfinite positive pH3 is not a regression input. The affected event signal
  is retained with explicit unavailable status for later C/D handling.
- One common analytical tier applies throughout the experiment. Mixed
  individual-corrected, pooled-corrected, and raw bases are prohibited.
- Regression failure does not alter or invalidate A/B prevalence outputs.

## Strict non-goals preserved

No C/D sample value, acquisition reduction, reference ratio, condition
summary, statistic, Holm adjustment, display offset, positive-domain display
safeguard, plot, report, Quarto source/render, CSV/JSON/RDS export, writer,
geometry generation, legacy path, EdU behavior, or POI behavior was added or
changed.

## Files read

Governance and decisions:

- workspace `AGENTS.md`;
- `docs/implementation/ph3_output_contract_gap_audit_2026-08-30.md`;
- `docs/implementation/ph3_output_contract_slice1_implementation_2026-08-30.md`;
- `/Users/sansamc/Documents/Codex/2026-08-30/ph3-output-contract-decisions/outputs/2026-08-30_final_ph3_output_qc_report_contract.md`;
- `/Users/sansamc/Documents/Codex/2026-08-30/ph3-implementation-decision-recovery/outputs/2026-08-30_ph3_implementation_decision_recovery_record.md`.

Canonical source and tests:

- `R/analysis.R`, `R/ph3.R`, and `R/ph3_outputs.R`;
- `DESCRIPTION` and `NAMESPACE`;
- `tests/testthat/helper-load.R`;
- `tests/testthat/helper-current-workflow.R`;
- `tests/testthat/test-ph3-output-model.R`;
- `tests/testthat/test-ph3-event-classification.R`.

Read-only Git status, log, and diff metadata were also inspected. No
experimental or generated analysis artifact was opened.

## Files modified

- `R/analysis.R`
- `tests/testthat/test-ph3-output-model.R`

## Files generated

Source-controlled implementation artifacts created by this task:

- `R/ph3_background.R`
- `tests/testthat/test-ph3-background-regression.R`
- `docs/implementation/ph3_output_contract_slice2_implementation_2026-08-31.md`

No report, plot, CSV, JSON, RDS, cache, binary, package archive, temporary
inspection dump, or experimental output was generated.

## Scientific changes and protected inputs

Experimental inputs read: **none**.
Experimental inputs modified: **none**.

Sample mappings changed: **none**.
Exclusions changed: **none**.
Gates changed: **none**.
pH3-positive membership changed: **none**.
Eligible, 4N, or below-4N membership changed: **none**.
Thresholds/boundaries changed: **none**.
DNA normalization changed: **none**.
Background correction changed: **owner-confirmed Slice 2 method added**.
Statistics changed: **none**.
Biological claims changed: **none**.
A/B outputs changed: **none; exact source objects are retained unchanged**.

## Focused deterministic tests prepared

The new in-memory tests use only values explicitly labeled `SYNTHETIC`. They
cover 99 versus 100 negative events, nonfinite fit and prediction inputs,
zero DNA variation and rank-validation order, deterministic direct validation
of nonfinite fitted coefficients and nonfinite required predictions,
inclusive coverage equality, individual and pooled coverage gaps,
experiment-wide pooled and raw propagation, typed inapplicable pooled
diagnostics, empty 4N and below-4N positive compartments, explicit nonfinite
positive-signal unavailability, exact A/B and membership preservation,
fail-closed mapping/membership validation, and a genuine Slice 1 builder
boundary.

## Tests and checks

R, Python, project code, tests, package loading/build/check, Quarto, analyses,
renders, snapshots, benchmarks, and installs run by Codex: **none**.

Runtime verification: **NOT RUN — local execution required**.

`git diff --check` was run as a text-format check and reported no whitespace
errors. This is not runtime or scientific verification.

### User-reported focused-test failure and correction

Report date: **2026-08-31**.

The user ran:

```sh
Rscript -e 'devtools::test(filter = "ph3-background-regression", stop_on_failure = TRUE)'
```

Reported result: **FAIL** at
`tests/testthat/test-ph3-background-regression.R:210`. The SYNTHETIC test used
finite response values alternating between `-1e308` and `1e308` and assumed
that `lm.fit()` would necessarily return a fit error, nonfinite coefficient, or
nonfinite prediction. On the user's R platform the fit remained finite, so the
platform-dependent expected-reason assertion failed.

Root cause: the test attempted to induce a numerical failure indirectly through
floating-point behavior that is not guaranteed across R/BLAS platforms. The
production requirements to reject nonfinite fitted coefficients and required
predictions were correct and were not weakened.

Correction: finite-solution validation is now routed through the small internal
`ph3_validate_background_solution()` helper used by the production fit path.
The focused test directly supplies a nonfinite coefficient and, separately,
finite maximum-magnitude coefficients whose required prediction is
deterministically nonfinite. The platform-dependent `lm.fit()` overflow
assumption was removed.

Corrected-tree runtime verification status:
Codex did not execute the failed command or any replacement test command.

### User-reported corrected-tree focused verification

Report date: **2026-08-31**.

The user reran:

```sh
Rscript -e 'devtools::test(filter = "ph3-background-regression", stop_on_failure = TRUE)'
```

User-reported result: **PASS**. No failed test name or error excerpt was
reported. Unexpected generated workspace files were not reported. Codex did
not execute or independently reproduce this result.

The directly affected multi-context regression command and complete package
test suite remain **NOT RUN — local execution required** unless separately
reported by the user or a laboratory member.

### User-reported full-suite safety-guard failure and correction

Report date: **2026-08-31**.

The user ran the complete package suite. It reached the `ph3-input-containment`
context and failed its existing production-namespace safety guard at
`test-ph3-input-containment.R:189`. The guard detected text matching the
prohibited pattern that could construct `event_identity` from a sequential
counter.

Root cause: `ph3_background_event_table()` correctly copied the validated
`data$event_identity`, but created a separate `event_row` counter with
`as.integer(seq_len(nrow(data)))` immediately afterward inside the same
`data.frame()` call. The strict lexical guard conservatively matched both
expressions as one potential identity construction. No sequential identity
fallback was implemented or reachable: `event_identity` continued to come
directly from the validated Slice 1 classification.

Correction: `event_row` is now computed before the `data.frame()` call and
passed by name. This is behavior-preserving: it retains the same integer row
counter while making the separation from canonical event identity explicit.
The safety assertion was not weakened or changed. Direct canonical identity,
exact mapping, and fail-closed validation remain unchanged.

Corrected-tree runtime verification status: **NOT RUN — local execution
required**. The user-reported focused pass predates this source-only correction;
the focused and complete test commands below must be rerun.

Independent rereview of this correction:

- scientific-integrity reviewer: **APPROVED — no blocking, important, or
  advisory findings**; and
- artifact/security reviewer: **APPROVED — no blocking, important, or advisory
  findings**.

Both reviewers confirmed that the production finite-coefficient and
finite-prediction requirements remain fail closed, the replacement tests are
deterministic, and no membership, tier, A/B, I/O, network, or unrelated change
was introduced.

## Independent reviews and resolution

The first read-only scientific-integrity and artifact/security reviews found:

- a blocking mismatch between the older per-set Slice 2 wording and the later
  owner-confirmed experiment-wide tier amendment;
- important fabricated zero counts for pooled fits that were not attempted;
- important ambiguity that paired one sample ID with a pooled-fit reason;
- important need for stronger fail-closed validation and a genuine Slice 1
  boundary test; and
- focused-test gaps for pooled coverage and empty positive compartments.

The first rereview then found two important remaining boundary gaps: a complete
technical-acquisition classification could be deleted or duplicated without an
exact acquisition-set check, and positivity/eligibility membership could be
changed while stale Slice 1 acquisition counts remained present.

Resolution:

- tier selection now follows the later owner-confirmed experiment-wide rule;
- unattempted pooled numerical diagnostics are typed missing values;
- individual and pooled triggers have separate IDs and reason codes;
- the Slice 1 boundary validates exact schemas and identity/membership mappings;
- a genuine `build_ph3_output_model()` integration test and the missing focused
  cases were added; and
- an orchestration review caught and corrected validation of the established
  colon-bearing `ph3-analysis-sha256:<digest>` analysis ID.
- the boundary now requires exactly one classification table for every
  authoritative acquisition, with exact set and cardinality equality; and
- each acquisition's eligible, positive-eligible, sub-4N, 4N,
  positive-sub-4N, and positive-4N event counts are reconciled against all five
  authoritative Slice 1 metric numerator/denominator records before fitting.

Final scientific-integrity rereview: **APPROVED — no blocking, important, or
advisory findings**.

Final artifact/security rereview: **APPROVED — no blocking, important, or
advisory findings**.

## Exact local verification instructions

Working directory for every command:

`/Users/sansamc/OMRF Dropbox/Chris Sansam/Synchd/Organized Dropbox/Projects and Data/Projects/Cowork_Projects/2026_Analyze_FACS_Data/facs_pseudocolor_workflow`

Prerequisites:

- R 4.2 or newer;
- package dependencies in `DESCRIPTION` already installed;
- `devtools` and `testthat` installed;
- no production input substitution or synthesis; and
- preferably a disposable clean working copy if local package tooling is
  configured to create artifacts.

Run the required focused test:

```sh
Rscript -e 'devtools::test(filter = "ph3-background-regression", stop_on_failure = TRUE)'
```

Expected success: all `ph3-background-regression` tests complete with zero
failures and errors. Expected workspace artifacts: none. Test-only temporary
files, if any, must remain under the R session temporary directory.

Run the directly affected pH3 regression set:

```sh
Rscript -e 'devtools::test(filter = "ph3-output-model|ph3-background-regression|ph3-acquisition-metrics|ph3-event-classification|ph3-replicate-condition-aggregation", stop_on_failure = TRUE)'
```

Expected success: every selected test completes with zero failures and errors.
Expected workspace artifacts: none.

Run the complete package test suite:

```sh
Rscript -e 'devtools::test(stop_on_failure = TRUE)'
```

Expected success: the complete suite finishes with zero failures and errors.
Expected workspace artifacts: none.

Stop and report `FAIL` for any failed test, error, unexpected warning suggesting
changed scientific behavior, basis mixture, schema/provenance/membership/A/B
mismatch, nonfinite result without its explicit reason, unexpected network
access, or unexpected generated workspace file. Do not weaken validation,
alter membership, or substitute data to force success. Report only the exact
command, `PASS` or `FAIL`, failed test names, a short relevant error excerpt,
and unexpected generated files; do not reproduce large successful logs.

## Assumptions, warnings, and unresolved uncertainty

- The later owner-confirmed recovery amendment controls the older Slice 2
  per-set wording because the recovery record was explicitly required input
  and defines one common experiment-wide analytical basis.
- Regression uses retained raw DNA and raw pH3 fields; normalized DNA remains
  immutable membership context only.
- Runtime behavior remains unverified until a laboratory member runs the local
  commands above.
- No owner-level scientific decision remains unresolved within this bounded
  Slice 2 scope.
