# pH3 output contract Slice 1 implementation handoff

**Date:** 2026-08-30
**Status:** CORRECTED — ALL REQUIRED LOCAL VERIFICATION COMMANDS USER-REPORTED PASS
**Scope:** Audit Slice 1 only — validated in-memory pH3 output/configuration
model

## Canonical source and baseline

Canonical source directory:

`/Users/sansamc/OMRF Dropbox/Chris Sansam/Synchd/Organized Dropbox/Projects and Data/Projects/Cowork_Projects/2026_Analyze_FACS_Data/facs_pseudocolor_workflow`

Branch at implementation start: `feature/ph3-output-contract-audit`
Baseline/current starting `HEAD`: `ead9e13e18c8f6243e9622161e0c9b475a42da0b`
Baseline description: `Merge pull request #21 from SansamLab/feature/ph3-replicate-condition-aggregation`

The pre-existing untracked audit files listed by `git status` at implementation
start were not modified or treated as implementation outputs:

- `docs/implementation/edu_output_contract_gap_audit_2026-08-30.md`
- `docs/implementation/ph3_output_contract_gap_audit_2026-08-30.md`
- `docs/implementation/ph3_slice6_output_contract_audit_2026-08-29.md`

## Exact authorized scope implemented

The production direct-identity pH3 configuration now requires an explicit,
versioned `ph3_output_contract` declaration with:

- one stable experiment ID;
- stable biological-replicate-set IDs that are separate from replicate labels,
  technical-replicate labels, model groups, and file order;
- exact condition IDs, labels, and designated `control`/`treatment` roles;
- explicit configured or not-configured reference state;
- an explicit, possibly empty comparison list whose directed conditions and
  applicable A–D outcomes are validated;
- declared or not-configured verified-geometry state; and
- explicit `not_computed` computed-geometry state for this slice.

One `ph3_output_model` is attached to a completed production pH3
`facs_analysis`. The model contains explicit ordered column/type/enum/reason
schemas; the four owner-confirmed A–D outcome definitions; experiment,
condition, replicate-set, sample, correction, reference, comparison, and
geometry state; the owner-confirmed required package/plot registry; and the
unchanged authoritative manifest, classifications, and Slice 4–5 tables.

Validation fails closed for missing/extra/reordered configuration fields,
unsupported versions, unsafe or duplicate stable keys, condition/manifest
mismatches, unknown or self comparisons, role-inconsistent directions,
duplicate comparison/outcome pairs, invalid reference state, invalid geometry
state, mixed analysis provenance, altered classification-derived tables,
altered Slice 5 tables, unstable sample identity, and remapped replicate-set
membership.

Correction basis is intentionally recorded as `not_computed` with stable reason
`background_correction_not_computed_slice_1`; no correction result is invented.
Computed geometry is intentionally recorded as `not_computed` with stable
reason `computed_geometry_not_generated_slice_1`; no geometry is generated or
represented as verified. Missing declared verified geometry is explicit and
does not suppress the retained numerical inputs.

## Owner-confirmed decisions applied

- The outcome registry order is exactly A, B, C, D.
- A is 4N pH3-positive prevalence with denominator `eligible_2to4n`.
- B is below-4N pH3-positive prevalence with denominator `eligible_2to4n`.
- C and D are signal outcomes whose analytical basis is required but remains
  not computed in Slice 1.
- Reference state is explicit; when configured, a reference condition must
  resolve to exactly one source sample ID in each biological-replicate set.
- Condition roles do not trigger statistics. Only separate explicit comparison
  records identify authorized control/treatment/outcome directions.
- The required root filenames, condition plot PDF/PNG names, and sample plot
  PDF/PNG path templates follow D0-11 exactly.
- Verified and computed geometry are separate states. A declared geometry-set
  index is not read or promoted to verified evidence in this slice.
- One model contains one stable experiment identity and exact active analysis
  provenance.

## Strict non-goals preserved

No background-regression or fallback engine, C/D value calculation, reference
ratio, condition summary, statistical test, Holm adjustment, display offset,
positive-domain safeguard, plot, Quarto source/render, CSV/JSON/RDS writer,
report export, geometry generation, legacy-field reinterpretation, experimental
input access, Git commit, push, PR, or merge was added or performed.

## Files read

Governance and source-of-truth records:

- workspace `AGENTS.md`;
- `docs/implementation/ph3_output_contract_gap_audit_2026-08-30.md`;
- `/Users/sansamc/Documents/Codex/2026-08-30/ph3-output-contract-decisions/outputs/2026-08-30_final_ph3_output_qc_report_contract.md`; and
- `/Users/sansamc/Documents/Codex/2026-08-30/ph3-implementation-decision-recovery/outputs/2026-08-30_ph3_implementation_decision_recovery_record.md`.

Canonical source/config/package files:

- `R/config.R`;
- `R/analysis.R`;
- `R/input.R`;
- `R/ph3.R`;
- `R/phase_quantitation.R`;
- `R/pseudocolor_helpers.R`;
- `examples/config_ph3.yml`;
- `DESCRIPTION`; and
- `NAMESPACE`.

Canonical tests and helpers:

- `tests/testthat/helper-load.R`;
- `tests/testthat/helper-current-workflow.R`;
- `tests/testthat/test-config.R`;
- `tests/testthat/test-ph3.R`;
- `tests/testthat/test-ph3-acquisition-metrics.R`;
- `tests/testthat/test-ph3-event-classification.R`;
- `tests/testthat/test-ph3-input-containment.R`;
- `tests/testthat/test-ph3-replicate-condition-aggregation.R`; and
- `tests/testthat/test-edu-output-contract.R` for established contract-test
  architecture only.

Read-only Git metadata/status/diff output was also inspected. No experimental or
generated analysis artifact was opened.

## Files modified

- `R/config.R`
- `R/analysis.R`
- `R/input.R`
- `examples/config_ph3.yml`
- `tests/testthat/helper-current-workflow.R`
- `tests/testthat/test-ph3-acquisition-metrics.R`
- `tests/testthat/test-ph3-event-classification.R`
- `tests/testthat/test-ph3-input-containment.R`

## Files generated

Source-controlled implementation artifacts created by this task:

- `R/ph3_outputs.R`
- `tests/testthat/test-ph3-output-model.R`
- `docs/implementation/ph3_output_contract_slice1_implementation_2026-08-30.md`

No report, figure, CSV, JSON, RDS, cache, binary, temporary inspection dump, or
experimental output was generated.

## Scientific changes and protected inputs

Experimental inputs read: **none**. No `.fcs`, `.wsp`, `.prism`, experimental
CSV/export, generated report, figure, or prior analysis result was accessed or
modified.

Experimental inputs modified: **none**.

Sample mappings changed: **none**. Stable replicate-set IDs and condition IDs
are newly required configuration identity; active manifest membership is
validated unchanged.

Exclusions changed: **none**.
Gates changed: **none**.
Thresholds/boundaries changed: **none**.
Normalization changed: **none**.
Statistics changed: **none**.
Biological claims changed: **none**.
Approved pH3-positive, eligible, 4N, and below-4N membership changed: **none**.
Approved A/B values changed: **none**; the new builder requires byte-identical
reproduction from active membership and retains the source tables unchanged.

## Tests and checks

Project code/tests/checks run by Codex: **none**.

Codex did not execute R, Python, project code, tests, package loading/build/check,
Quarto, analyses, renders, snapshots, benchmarks, or installs. Only read-only
source/Git inspection and text edits were performed. `git diff --check` was used
as a text-format check and reported no whitespace errors; it is not scientific
or runtime verification.

### Reported local verification failure and correction

Correction date: **2026-08-31**.

The user reported that this required command failed before any test could load:

```sh
Rscript -e 'devtools::test(filter = "ph3-output-model", stop_on_failure = TRUE)'
```

Reported failure: `devtools::test()` stopped inside `load_all()`/parsing with
`unexpected '}'` at `R/ph3_outputs.R:17`, at the closing brace following the
expected-schema error message.

Root cause: `ph3_output_contract_exact_names()` opened a
`ph3_output_contract_fail(...)` call at line 13 but omitted the call's closing
`)` after the `paste0(...)` argument. The parser therefore reached the closing
brace while the function call was still open.

Correction: inserted the single missing `)` immediately after the completed
`paste0(...)` argument. No validation, schema, assertion, scientific method, or
test expectation was weakened, removed, or otherwise changed.

Corrected-tree runtime verification status:
**NOT RUN — local execution required**. Codex performed static source and diff
inspection only and did not rerun the failed command.

### Reported focused-test identity failure and correction

Correction date: **2026-08-31**.

After the parse correction, the user reran:

```sh
Rscript -e 'devtools::test(filter = "ph3-output-model", stop_on_failure = TRUE)'
```

The package and focused test context loaded successfully. The reported result
was `[FAIL 3 | WARN 0 | SKIP 0 | PASS 8]`. All three errors stopped at the same
fail-closed guard with
`PH3 output-contract validation failed [replicate_set_manifest_mismatch]: stable replicate-set IDs do not map exactly to the active manifest`
from the reported pre-correction location `R/ph3_outputs.R:579`.

Root cause: the compared replicate-set IDs and integer replicate indices were
equal, but `unique(manifest[...])` retained nonsequential source row names
(`1`, `3`) while the newly constructed expected two-row data frame used
canonical row names (`1`, `2`). Because `identical()` compares attributes as
well as values and types, the irrelevant row-name difference caused valid
SYNTHETIC mappings to fail.

Correction: assign canonical row names to the derived two-column
`manifest_replicate_mapping` before the existing `identical()` comparison. The
guard still requires exact column order, values, scalar types, and row order;
no mapping, validation, or assertion was removed or weakened.

Corrected-tree runtime verification status after this correction:
**NOT RUN — local execution required**. Codex performed static source and diff
inspection only and did not rerun the focused test.

### User-reported corrected-tree focused verification

Report date: **2026-08-31**.

The user reran:

```sh
Rscript -e 'devtools::test(filter = "ph3-output-model", stop_on_failure = TRUE)'
```

User-reported result: **PASS**. No failed test name or error excerpt was
reported. Unexpected generated workspace files were not reported. Codex did
not execute or independently reproduce this result.

At the time of this focused-test report, the standalone example-configuration
validation had not yet been run; its subsequent **PASS** is recorded below.

### User-reported broader corrected-tree verification

Report date: **2026-08-31**.

The user reported **PASS** for the directly affected regression set:

```sh
Rscript -e 'devtools::test(filter = "config|ph3-output-model|ph3-acquisition-metrics|ph3-event-classification|ph3-input-containment", stop_on_failure = TRUE)'
```

The user also reported **PASS** for the complete package test suite:

```sh
Rscript -e 'devtools::test(stop_on_failure = TRUE)'
```

No failed test names or error excerpts were reported. Unexpected generated
workspace files were not reported. Codex did not execute or independently
reproduce either result.

### User-reported standalone example-configuration verification

Report date: **2026-08-31**.

The user reported **PASS** for the standalone example-configuration validation:

```sh
Rscript -e 'devtools::load_all(); cfg <- read_facs_config("examples/config_ph3.yml"); stopifnot(inherits(cfg, "facs_config"))'
```

No error excerpt was reported. Unexpected generated workspace files were not
reported. Codex did not execute or independently reproduce this result. All
required local verification commands are now user-reported **PASS**.

### Required local verification

Working directory for every command:

`/Users/sansamc/OMRF Dropbox/Chris Sansam/Synchd/Organized Dropbox/Projects and Data/Projects/Cowork_Projects/2026_Analyze_FACS_Data/facs_pseudocolor_workflow`

Prerequisites:

- R 4.2 or newer;
- package dependencies in `DESCRIPTION` already installed;
- `devtools` and `testthat` installed;
- no command is permitted to substitute or synthesize production inputs; and
- run in a disposable/clean working copy if local package tooling is configured
  to create artifacts.

Run the required focused test:

```sh
Rscript -e 'devtools::test(filter = "ph3-output-model", stop_on_failure = TRUE)'
```

Expected success: the `ph3-output-model` tests complete with no failed tests or
errors. Expected workspace files: none. Test-only temporary files, if any, must
remain under the R session temporary directory.

Then run the directly affected focused regression set:

```sh
Rscript -e 'devtools::test(filter = "config|ph3-output-model|ph3-acquisition-metrics|ph3-event-classification|ph3-input-containment", stop_on_failure = TRUE)'
```

Expected success: all selected tests complete with no failed tests or errors.
Expected workspace files: none.

Validate the edited example configuration without analyzing data:

```sh
Rscript -e 'devtools::load_all(); cfg <- read_facs_config("examples/config_ph3.yml"); stopifnot(inherits(cfg, "facs_config"))'
```

Expected success: exit status 0 and a validated `facs_config`; no input path is
read and no analysis is performed. Expected workspace files: none.

Finally, run the complete package test suite:

```sh
Rscript -e 'devtools::test(stop_on_failure = TRUE)'
```

Expected success: the full suite completes with no failed tests or errors.
Expected workspace files: none.

Stop and report `FAIL` if any command has a failed test, error, unexpected
warning suggesting changed scientific behavior, source/classification/table
mismatch, reference-resolution failure, unstable identity failure, unexpected
network access, or unexpected generated workspace file. Do not weaken a
validation or alter membership/configuration to force success. Report only the
exact command, `PASS` or `FAIL`, failed test names, a relevant short error
excerpt, and any unexpected generated files; do not reproduce a large
successful log.

## Independent reviews

Scientific-integrity reviewer: `/root/scientific_review` (read-only; did not
author or edit the implementation).

- Final classification: **no blocking, important, or advisory findings**.
- Confirmed no fabricated/substituted observations, experimental-input access,
  changed membership, exclusions, gates, thresholds, normalization, statistics,
  biological claims, or A/B values.
- Confirmed the exact source-table reproduction checks, explicit reference and
  comparison direction, stable manifest-bound replicate-set identity, and
  truthful `not_computed` states.
- During review, the remapped-identity test was clarified to mutate only the
  in-memory SYNTHETIC manifest ID and expect
  `replicate_set_manifest_mismatch`; the reviewer confirmed the resolved test
  reaches the intended guard without changing the configuration digest or
  scientific source tables.

Artifact/security reviewer: `/root/artifact_security_review` (read-only; did
not author or edit the implementation).

- Final classification: **no blocking or unresolved important findings**.
- The initial important finding was that a declared geometry-set reference
  accepted any nonempty string. Resolution: the schema now requires the exact
  key `geometry_set_index_id`, validates it as a stable filesystem-safe ID, and
  includes a SYNTHETIC traversal-like rejection case. The reviewer confirmed
  the finding resolved with no new issue.
- One advisory remains: the inert registry template retains raw `{sample_id}`.
  It is not expanded, read, or written in Slice 1. Collision-safe encoding and
  output-root containment preflight remain explicitly required before a later
  writer/publisher may consume it.
- Confirmed no network/telemetry/credential access, shell/process invocation,
  file write/delete, renderer, hidden fallback, generated artifact, or unrelated
  production behavior was introduced.

Both reviewers inspected the latest complete diff after resolutions. Neither
reviewer executed project code or tests. At review time, runtime verification
was **NOT RUN — local execution required**; subsequent user-reported results are
recorded above.

Renewed correction reviews on **2026-08-31**:

- `/root/scientific_review` statically confirmed that the delimiter structure at
  `R/ph3_outputs.R:10-19` is now balanced, the correction only closes the
  existing validation-failure call, and no membership, A/B value, mapping,
  gate, threshold, normalization, statistic, or biological claim changed.
- `/root/artifact_security_review` statically confirmed that the correction
  introduces no validation-boundary, I/O, write/delete, network, process, path,
  fallback, or security change. The prior inert `{sample_id}` template advisory
  remains unchanged and deferred.
- Renewed-review result: **no blocking, important, or new advisory findings**.
  Neither reviewer executed project code or tests; corrected-tree verification
  was **NOT RUN — local execution required** at review time. Subsequent
  user-reported results are recorded above.

Renewed identity-guard correction reviews on **2026-08-31**:

- `/root/scientific_review` confirmed that `unique()` retained noncanonical
  source row names, that `identical()` compared those attributes, and that
  canonicalizing only the derived mapping's row names is the minimal correction.
  Exact ordered columns, replicate-set IDs, integer indices, and row order remain
  required. No scientific-integrity finding remains.
- `/root/artifact_security_review` confirmed no validation weakening, fallback,
  I/O, path, network, process, security, or scope change. Its documentation-only
  advisory about the shifted source line was resolved by labeling line 579 as
  the reported pre-correction location.
- Renewed-review result: **no blocking, important, or new advisory findings**.
  The existing inert `{sample_id}` template advisory remains deferred. Neither
  reviewer executed project code or tests; corrected-tree verification was
  **NOT RUN — local execution required** at review time. Subsequent
  user-reported results are recorded above.

## Assumptions, warnings, and remaining uncertainty

- Stable experiment, condition, replicate-set, comparison, outcome, artifact,
  and geometry-declaration key spellings are bounded architect-owned interface
  decisions authorized for this slice.
- Sample identity is taken only from the active immutable export-manifest
  `sample_id`; technical acquisitions may share a sample identity only when the
  source manifests explicitly say so.
- The sample plot registry intentionally retains `{sample_id}` as a path
  template. Collision-safe filesystem encoding and path publication remain
  later-slice responsibilities; the scientific sample ID is not rewritten.
- A declared stable geometry-set index ID is configuration state only. Paths,
  contents, digests, bundle bindings, verified status, and geometry are not read or
  validated until the separately authorized geometry slice.
- `not_computed` correction and computed-geometry states are intentional, stable
  Slice 1 states, not analytical failures or substitutes.
- All required local verification commands were subsequently user-reported
  **PASS**. Codex did not independently reproduce them, and unexpected generated
  workspace files were not reported.
