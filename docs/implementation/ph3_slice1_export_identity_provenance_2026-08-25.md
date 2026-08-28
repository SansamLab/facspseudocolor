# pH3 Slice 1 export identity and provenance implementation

**Status:** IMPLEMENTED — LOCALLY VERIFIED BY USER

## Scope and provenance

- Canonical source: `/Users/sansamc/OMRF Dropbox/Chris Sansam/Synchd/Organized Dropbox/Projects and Data/Projects/Cowork_Projects/2026_Analyze_FACS_Data/facs_pseudocolor_workflow/`
- Branch: `feature/ph3-export-identity-provenance`
- Starting commit: `5de77a9e6f45c2b1f6f7675a5a37aa67a312a2bb`
- Starting state: no tracked modifications; the protected pre-existing untracked
  `PROMPT_WHOLE_DATASET_ORACLE_AUDIT.md` was not opened or changed.
- Authorities: workspace-root
  `ph3_scientific_owner_approval_2026-08-24.md`,
  `ph3_implementation_plan_and_test_spec_5de77a9.md`, and
  `ph3_current_implementation_audit_5de77a9.md`.

This slice changes only upstream export identity, immutable local provenance,
and optional geometry linkage. It does not implement Slice 2, pH3 metrics,
denominators, eligibility, aggregation, result schemas, plots, reports, gates,
thresholds, normalization, statistics, sample mappings, exclusions, or
biological interpretations.

## Files

Read: the three authority records above;
`python/export_flowjo_populations.py`,
`python/export_flowjo_gate_geometry.py`, `tools/flowjo-orchestration.R`,
`tests/testthat/test-python-boundary.R`, `docs/PYTHON_INTERFACE.md`, `README.md`,
and `NEWS.md`. Repository file names under `python/`, `tests/`, and `docs/` were
listed only to locate the established boundary tests and documentation.

Modified: `python/export_flowjo_populations.py`,
`python/export_flowjo_gate_geometry.py`, `tools/flowjo-orchestration.R`,
`tests/testthat/test-python-boundary.R`, `docs/PYTHON_INTERFACE.md`, `README.md`,
and `NEWS.md`.

Generated as source/test/documentation (not analysis output):
`python/export_contract.py`, `tests/python/test_export_contract.py`,
`tests/python/test_flowkit_source_index_semantics.py`, and this report. No
experimental or generated historical input was read, modified, copied, or
derived. No cache, binary, package archive, rendered report, or analysis output
was generated.

`python/export_contract.py` is the one narrowly required additional production
file: it gives both exporters one canonical hashing, atomic-finalization,
artifact-ledger, and geometry-supplement implementation instead of duplicating
security-sensitive contract logic.

## Identity contract

Production direct identity is
`<acquisition_id>:event_index:<FlowKit get_gate_events index>`. Population CSVs
carry `acquisition_id`, character `event_index`, `event_identity`,
`identity_source`, method ID/version, `duplicate_occurrence`, operation ID,
`export_manifest_digest`, and the full-manifest sidecar reference. The row
digest is SHA-256 of canonical `manifest_binding_object_v1`, which covers the
immutable operation/method/provenance/acquisition/request scope while excluding
time and artifact/full-manifest hashes. The detached `export-manifest.sha256`
separately covers the complete finalized manifest and artifact ledger. This
distinction provides a truthful non-self-referential row binding. Duplicate direct indices within one
sample/population are fatal. The R boundary forces identity columns to
character and the former `seq_len(...)` replacement identity was removed.
There is no row-order, file-order, measurement, DNA/target, or nearby-file
fallback. Composite identity remains blocked and unimplemented because exact
pre-coercion emitted channel tokens are not guaranteed by this path.

Direct semantics are **verified for the user-tested FlowKit 1.3.2 environment**.
Static inspection alone showed only that the FlowKit frame index is passed
through, so production creation continues to require an exact
operator-supplied supported FlowKit version and the explicit
`--direct-index-semantics-verified` attestation. On 2026-08-28, the user ran the
pinned-environment SYNTHETIC parent/child/acquisition test successfully against
a manually created two-acquisition FlowJo workspace. The test compares gated
raw measurements with source-acquisition rows at every reported index. This
result does not establish semantics for other FlowKit versions or environments,
which require their own verification before attestation.

## Manifest and fail-closed behavior

One population operation writes `export-manifest.json` plus detached
`export-manifest.sha256`. Canonical bytes are sorted-key, compact UTF-8 JSON
with one trailing newline; the digest is not included in the hashed document.
The schema, operation ID, timezone timestamp, profiles/methods, workspace path
and hash, software versions, transformation/compensation context, explicit
approval metadata, acquisition-to-prefix records, population hierarchy/counts,
and artifact ledger are recorded. Population artifacts are one file per
acquisition/requested population. Their ledger entries include path, SHA-256,
size, row count, schema, identity method, empty status, acquisition, sample,
population, gate path, channels, and operation linkage. Hashes are checked
before atomic manifest finalization. Existing completed operations and nonempty
operation directories are not overwritten. Failed operations have no finalized
manifest/digest pair.

Production requires explicit owner, approver, date, local approval record,
positivity method, transform and compensation contexts, unique acquisition and
sample IDs, prefixes, FCS references/hashes, exact supported FlowKit version,
and the environment-specific verified direct-index attestation. Values are
never inferred.
Production loads only metadata-declared FCS references confined below the
explicit FCS directory, requires exact workspace-sample coverage, and verifies
every FCS SHA-256. Metadata keys are allowlisted so credentials cannot be
serialized.
Header-only zero-event population output is intentionally empty only when the
gate exists for its acquisition and its zero count/schema are recorded; an
absent population is fatal in production.

The explicit `legacy_count_only_unverified_v1` profile records unavailable
fields and a prominent warning. It cannot claim production identity,
containment, or verified geometry. Its identity fields are blank, and R
orchestration refuses legacy or ambiguous profiles as production.
R consumes and returns the immutable operation artifacts directly. It no longer
creates or silently overwrites unledgered per-sample CSV copies. At that return
boundary it invokes the local verifier, which rechecks the finalized full
manifest digest and each consumed artifact's hash, size, row count, ordered
schema, plus row operation ID, binding digest, acquisition/sample IDs, and
profile. Any mutation is fatal.

## Geometry linkage

No requested geometry is explicitly `not_requested`. Requested geometry must
extend the completed population operation in the same directory. Its immutable
supplement links the population-manifest digest and rechecks operation ID,
workspace hash, exact population key, acquisition/sample coverage, gate path,
and ordered gate channels. It never relies on filenames. A mismatch is fatal;
failed finalization removes the geometry CSV. R orchestration writes one linked
operation artifact at a time and does not rewrite it. Geometry remains
nonquantitative.
The geometry artifact ledger and supplement carry an identical structured
multi-acquisition linkage scope, reconciled against all geometry rows, with one
acquisition/sample/full-gate-path/ordered-channel entry per scope member.

## SYNTHETIC tests added

- `tests/python/test_export_contract.py`: deterministic identity/acquisition
  scoping, exact FCS/hash binding, missing-index failure/no sequential fallback,
  blank legacy identity, mixed-acquisition missing-gate rejection, gate-present
  zero ledgers, complete production schema, exact geometry coverage,
  conditional production,
  canonical digesting/nonfinite rejection, schema and confined-path validation,
  metadata allowlisting, mutation detection, overwrite refusal, and operation-ID
  geometry rejection, consumption-boundary row/ledger reconciliation, and
  structured multi-acquisition geometry scope. Exact expected outcomes are assertions; failures are
  blocking because they protect identity and immutable provenance.
- `tests/python/test_flowkit_source_index_semantics.py`: user-supplied,
  unmistakably SYNTHETIC two-acquisition fixture; checks raw/transformed index
  equality, source-row measurement equality, uniqueness, and G1/pH3-positive
  containment in Single Cells. Failure
  is blocking for production direct identity.
- `tests/testthat/test-python-boundary.R`: statically confirms the R boundary
  preserves identity columns as character and contains no former sequence
  fallback. Failure is blocking.

No fixture encodes a biological conclusion. Agents did not execute these tests;
the user-run results are recorded below.

## Local verification

Working directory for every command is the absolute canonical directory
recorded above.
Prerequisites: Python 3 used for production; the exact candidate pinned FlowKit
version plus pandas/lxml/numpy; R >= 4.2; declared package dependencies;
`testthat`, `withr`, and `devtools`; Quarto only if the laboratory later elects
to render (no render is required here).

1. Contract/manifest/linkage tests:
   `python3 -m unittest -v tests/python/test_export_contract.py`
   Expected: all tests `OK`; only controlled OS temporary directories are
   created and removed. Stop on any failure or residual file.
2. Pinned direct-index proof:
   set `SYNTHETIC_FLOWJO_WORKSPACE`, `SYNTHETIC_FLOWJO_FCS_DIR`,
   `SYNTHETIC_SINGLE_CELLS_GATE`, `SYNTHETIC_G1_GATE`, and
   `SYNTHETIC_PH3_POSITIVE_GATE` to an explicitly labeled SYNTHETIC fixture,
   then run
   `python3 -m unittest -v tests/python/test_flowkit_source_index_semantics.py`.
   Expected: one test `OK`; no project file is generated. Stop if the fixture is
   not synthetic, has fewer than two acquisitions, gate names are ambiguous, an
   index differs across raw/transformed exports, identities duplicate, or child
   containment fails. Do not attest production support after any failure.
3. R boundary test:
   `Rscript -e 'devtools::test(filter = "python-boundary")'`
   Expected: the focused file passes; ordinary testthat temporary files only.
4. Full package tests:
   `Rscript -e 'devtools::test()'`
   Expected: zero failures/errors; ordinary session-temporary files only.
5. Temporary build/check because R orchestration and package tests/docs changed:
   `ph3_check_dir=$(mktemp -d)` then
   `(cd "$ph3_check_dir" && R CMD build "/Users/sansamc/OMRF Dropbox/Chris Sansam/Synchd/Organized Dropbox/Projects and Data/Projects/Cowork_Projects/2026_Analyze_FACS_Data/facs_pseudocolor_workflow")` then
   `R CMD check --no-manual "$ph3_check_dir"/facspseudocolor_*.tar.gz --output="$ph3_check_dir"`.
   Expected: build succeeds and check reports no errors, warnings, or notes;
   expected generated files are only the tarball and check directory under
   `/tmp`. Stop on any error/warning/note or generated project-tree artifact.
6. Repository inspection:
   `git status --short`.
   Expected: only the files listed in this report plus the protected pre-existing
   untracked audit prompt. Stop for caches, temporary files, archives, rendered
   outputs, or any unrelated modification.

Return only command, `PASS`/`FAIL`, failed test names, a relevant error excerpt,
check status, and unexpected generated files; do not return large successful
logs.

### User-run verification results

On 2026-08-28, the user reported all required local verification complete:

- `python3 -m unittest -v tests/python/test_export_contract.py`: **PASS**,
  18/18.
- `python3 -m unittest -v tests/python/test_create_synthetic_flowjo_index_fixture.py`:
  **PASS**, 7/7.
- `python3 -m unittest -v tests/python/test_flowkit_source_index_semantics.py`:
  **PASS**, 1/1, using the manually created two-acquisition SYNTHETIC FlowJo
  workspace with FlowKit 1.3.2 and the documented environment variables.
- `Rscript -e 'devtools::test(stop_on_failure = TRUE)'`: **PASS**, 397 tests,
  0 failures, warnings, or skips.
- Temporary `R CMD build` followed by `R CMD check --no-manual`: **PASS**,
  final `Status: OK`.
- `git diff --check`: **PASS**, clean.

These commands and results were executed and reported by the user. No agent ran
R, Python, FlowKit, tests, package build/check, or project analysis code.

**Verification:** COMPLETE — USER-RUN LOCAL VERIFICATION PASSED

## SYNTHETIC FlowJo source-index fixture preparation

**Status:** SOURCE PREPARED AND VERIFIED — USER-RUN LOCAL VERIFICATION PASSED

The deterministic local generator is
`tools/create_synthetic_flowjo_index_fixture.py`; its manual generation,
FlowJo gating, and verification guide is
`tests/fixtures/SYNTHETIC_flowjo_index_semantics/README.md`. FlowIO was selected
because the pinned FlowKit environment already supplies it and its installed
`create_fcs` API writes FCS 3.1 float data. The import is isolated inside the
test-only writer and failure gives an explicit pinned-environment prerequisite;
no production dependency was added.

The generator writes exactly two explicitly named SYNTHETIC acquisitions of
360 deterministic events each. Construction uses no random sampling: each has
40 deliberately off-diagonal debris/non-single-cell points and 320 intended
single-cell points divided into separated G1-like/later-DNA and low/high
synthetic-pH3 signal groups. Acquisition 2 has a small declared arithmetic
offset. Four synthetic-only channel names and metadata declare the generator,
construction, acquisition, nonexperimental status, absence of biological
interpretation, and truthful lack of compensation/transformation. Suggested
manual gate regions are fixture mechanics only, never biological thresholds.

The generator requires a user-supplied resolved destination containing
`SYNTHETIC`, refuses filesystem-root and symbolic-link targets, refuses existing
named outputs by default, writes through same-directory temporary files, and
offers a narrowly scoped `--overwrite` for only its two filenames. It performs
no network, telemetry, upload, shell, recursive deletion, experimental-file
discovery, or production analysis behavior. Generated `.fcs` and `.wsp` files
remain local and untracked; narrow ignore rules cover only binaries placed
directly in this fixture recipe directory.

Both acquisitions are staged and synchronized before either final filename is
installed. A second-acquisition failure removes only the generator's exact
staged targets; overwrite mode temporarily preserves and restores either exact
pre-existing named fixture file if installation fails.

First local command, from the canonical repository directory in the pinned
FlowKit/FlowIO Python environment:

```bash
python3 tools/create_synthetic_flowjo_index_fixture.py \
  --output-dir /private/tmp/SYNTHETIC_ph3_flowkit_index_fixture
```

Expected: exactly `SYNTHETIC_INDEX_ACQUISITION_1.fcs` and
`SYNTHETIC_INDEX_ACQUISITION_2.fcs`, followed by the manual FlowJo workspace
steps and exact environment-variable verification command in the fixture
README. Stop on any missing/extra output, prerequisite error, ambiguous or
empty gate, index/measurement/containment failure, experimental-file selection,
or unexpected repository artifact. The focused pure-helper/safety test is:

```bash
python3 -m unittest -v tests/python/test_create_synthetic_flowjo_index_fixture.py
```

Expected: seven tests report `ok` and the result is `OK`; only an automatically
removed OS temporary directory is used. No agent ran the generator, FlowJo,
FlowKit verification, or tests, and no agent generated an `.fcs` or `.wsp`.
The user subsequently reported the helper test passed 7/7 and the manual
two-acquisition SYNTHETIC FlowJo/FlowKit verification passed 1/1.

Fixture-preparation provenance: canonical source and branch remain those above,
at starting commit `5de77a9e6f45c2b1f6f7675a5a37aa67a312a2bb`, with all pre-existing Slice 1
changes preserved. Experimental inputs read: none. Experimental inputs
modified: none. Sample mappings, exclusions, experimental gates, thresholds,
normalization, statistics, or biological claims changed: none. The protected
audit prompt and generated Python cache were not opened or modified; the cache
must be removed before commit after separate authorization. Direct-index
semantics are **VERIFIED for the user-tested FlowKit 1.3.2 environment** by the
successful user-run test against the manually created two-acquisition SYNTHETIC
workspace. Other FlowKit versions or environments remain unverified until the
same test passes there.

Fixture-review resolution added fail-closed checks for exactly the two
generator filenames/sample IDs and their embedded SYNTHETIC provenance, direct
child gate paths, and nonempty strict containment that excludes the declared
debris group. Artifact review additionally required path-free safe acquisition
prefixes, explicit operation-directory confinement, exclusive unpredictable
temporary CSV files with atomic no-clobber publication and cleanup in both exporters,
operation-level fixture staging/rollback, and a root-only ignore rule for the
protected unrelated prompt. These changes alter no scientific method.

Final artifact re-review further required that a failed overwrite rollback
never delete its only recoverable original: each backup is now removed from the
rollback ledger only after successful restoration, and any retained backup path
is surfaced in the raised error for manual recovery. Immutable CSV, manifest,
and digest publication now uses same-filesystem atomic hard-link creation that
fails when a concurrent exporter has already claimed the final path; it never
replaces the concurrent file. Focused tests cover retained-backup reporting and
no-clobber preservation. Agents did not execute these tests; the user later
reported the complete focused Python verification passed.

### Local FlowKit bookkeeping-column correction

The user's first local source-index run reached FlowKit but failed before the
measurement comparison: `Workspace.get_gate_events()` adds a `sample_id`
bookkeeping column, while `Sample.as_dataframe(source="raw")` contains detector
measurements only. Static inspection of the installed FlowKit implementation
confirms that `get_gate_events()` unconditionally inserts `sample_id` after
selecting the source-indexed event rows.

The focused proof now requests the source acquisition with documented simple
PnN column labels, requires exactly the four declared SYNTHETIC detector
columns, requires `sample_id` to be the sole gated-table bookkeeping column,
and validates that every bookkeeping value equals the current workspace sample
ID. It then compares only detector measurements—without renumbering or resetting
the index—against the source acquisition at the exact gated indices. Raw versus
transformed index equality, uniqueness, acquisition scoping, direct-child gate
paths, and strict containment checks remain unchanged. This is a test-harness
correction, not a weakening of the source-index proof or a production behavior
change. Agents did not run the corrected test. On 2026-08-28, the user ran it
locally against the manually created two-acquisition SYNTHETIC FlowJo workspace:
one test reported `ok` and the final result was `OK`. The pinned-environment
source-index verification is therefore **PASS — user-run local verification**.
The user reported FlowKit version `1.3.2`, completing the environment-version
provenance for this verification.

### Rollback failure-injection correction

The locally reported failures in
`test_restore_failure_retains_and_reports_recoverable_backup` were isolated to
the test's failure-injection boundary rather than the rollback preservation
rule. First, the test patched `replace` on Python's shared `os` module object;
the generator now routes only its exact atomic install, backup, and restore
moves through a narrow internal `replace_file` wrapper, and the test patches
only that wrapper. Second, the user rerun showed that the injected failure was
not reached on macOS: `TemporaryDirectory` presented a `/var/...` path while
the generator's required destination validation canonicalized it to
`/private/var/...`, so exact target-path equality was false. The test now uses
the same `require_synthetic_destination()` result as `write_fixture()` before
constructing its expected paths. Its saved real wrapper still performs every
non-injected move.

The test continues to require both a simulated second-install failure and a
simulated first-backup restoration failure. It still asserts that the first
original's bytes remain at the exact backup path reported by the raised error,
that the second original is restored, and that the failed final path is not
present. No rollback check was removed or weakened, and production behavior
still preserves and reports every backup whose restoration fails. Agents did
not run the corrected test. On 2026-08-28, the user reran the single focused
test locally after the macOS path-alias correction; one test reported `ok` and
the final result was `OK`. Focused status is therefore **PASS — user-run local
verification**.

## Review and unresolved uncertainty

Initial independent reviews found blocking/important gaps in FCS/hash binding,
per-acquisition missing/empty handling, legacy identity separation, the digest
field, source-row proof, geometry completeness/cleanup, path/schema validation,
metadata allowlisting, nonfinite JSON, and sanitized-name collisions. These
were resolved by exact metadata-selected hash-verified FCS loading; fatal
per-acquisition absence and gate-present zero ledgers; blank legacy identity
plus R rejection; distinct canonical binding and full-manifest digests;
immutable per-acquisition artifacts with no R-side copies; source-row comparison;
exact operation-complete geometry linkage and cleanup; path/schema validators;
metadata allowlists; `allow_nan = FALSE`; and collision checks. Independent
re-review then identified a required row-level digest and unledgered R-created
copies. Those were resolved with the separately canonicalized immutable
manifest-binding digest, per-acquisition exporter artifacts, exact artifact
ledger linkage, and removal of all R-side population CSV writes. Focused tests
were expanded for FCS hashes, mixed missing gates, intentional zero exports,
legacy blanks, production schema, and exact geometry coverage. Final independent
scientific-integrity and artifact/security re-reviews found no remaining
blocking, important, or advisory findings.
The final important findings—verification at the actual R consumption boundary
and explicit multi-acquisition geometry artifact linkage—were resolved with a
local fail-closed verification mode invoked by R and a row-reconciled structured
geometry linkage scope in both artifact ledger and supplement.

Assumptions: FlowKit exposes `__version__`; configured FlowJo gate paths returned
by `find_matching_gate_paths()` have the same semantics in population and
geometry calls; operator metadata is authoritative only when explicitly
provided. Direct FlowKit index semantics are verified only for the user-tested
FlowKit 1.3.2 environment; another version or environment requires the same
pinned SYNTHETIC verification before attestation. Missing or ambiguous
provenance, mapping, gates, identity, artifacts, or linkage remains fatal in
production. Composite identity remains blocked and unimplemented because exact
pre-coercion emitted channel tokens are not guaranteed by this path.

Slice 2 and every downstream pH3 metric remain unimplemented.
