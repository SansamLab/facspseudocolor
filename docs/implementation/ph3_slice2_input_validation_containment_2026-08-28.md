# pH3 Slice 2 input validation and containment implementation

**Status:** CORRECTED AFTER USER-REPORTED FAIL — NOT RUN AFTER CORRECTION,
LOCAL RERUN REQUIRED

## Scope and provenance

Canonical source: `facs_pseudocolor_workflow/`, branch
`feature/ph3-input-validation-containment`, exact baseline
`3fe0f2816023538987a441a9b01622096a001031`.

This slice adds only fail-closed R consumption of the merged Slice 1 export
contract and independent exact-identity containment of FlowJo G1 and
pH3-positive populations within the same approved Single Cells acquisition.
It does not add eligibility, DNA-region classification, biological metrics,
boundary sensitivity, aggregation, canonical biological output schemas,
reports, or quantitative plots.

Experimental inputs read or modified: none. Protected FCS, WSP, Prism,
experimental CSV, historical generated, cache, and archive files were not
opened, hashed, copied, or derived. Python exporters were inspected only as
source and were not modified.

## Implemented contract

PH3 configuration now requires an explicit production or legacy input profile.
Production requires explicit export-operation directories. R verifies the
detached full-manifest SHA-256; exact manifest, export-profile, identity-method,
binding, digest, workspace, software, approval, acquisition, population, and
artifact subcontracts; exact PH3 population/count coverage; confined regular-
file paths; artifact SHA-256, byte size, ordered schema, row count, channels,
ledger linkage; and row-level operation/manifest/acquisition/sample/profile/
method bindings. Relative operation directories cannot escape `data_dir`.
Distinct operation directories must have unique operation IDs. Recorded source
FCS references follow the merged Slice 1 `Path` contract: absolute, UNC,
drive-qualified, and parent-traversal variants are rejected, harmless empty or
`.` components are normalized, and at least one component must remain. The
referenced experimental files are not opened. All required production metadata
strings must contain a non-whitespace character.

For each acquisition, G1 and pH3-positive are checked independently against
Single Cells by exact character `event_identity`. No coercion, rounding,
reformatting, row-order identity, numeric reconstruction, duplicate suffixing,
or nearby-file fallback exists. Source-index tokens must use the exact Slice 1
canonical nonnegative-integer character form. Direct duplicates, cross-acquisition linkage,
unmatched children, excess occurrences, and ambiguous mappings are fatal.
Children may overlap. Zero-row children require matching intentional-empty
flags, zero counts in the immutable count report, exact schema/artifact
verification, gate metadata, and approval metadata. Required configured
channels retain the existing finite-event minima of 10 Single Cells, 2 G1,
and 0 pH3-positive events.

One deterministic row per acquisition and child is retained as
`ph3-input-containment-1.0.0`, attached to the input report and analysis
provenance together with immutable nested manifest references. Explicit legacy
mode retains historical count-only numerical behavior, warns once per
analysis, and reports one prefix-by-child `unverified_count_only_legacy` row;
unavailable acquisition/identity/manifest proof remains `NA` and cannot claim
production containment. Verified production event frames are retained only
transiently between verification and normalization, then stripped before the
analysis object is returned. EdU and POI result provenance is unchanged.
Production analysis stops at normalized events and containment/provenance: it
does not call or attach the existing legacy `quantify_ph3()` biological
summaries. The exported `quantify_ph3()` API also fails closed unless the
analysis input profile is exactly `legacy_count_only_unverified_v1`, preventing
manual production-profile bypass. Explicit legacy analysis retains its
historical quantitation and conditional `event_index` safety failures, but those
checks do not upgrade legacy containment above `unverified`.

## Files

Modified: `DESCRIPTION`, `R/config.R`, `R/input.R`, `R/analysis.R`, `R/ph3.R`,
`tests/testthat/helper-current-workflow.R`, `tests/testthat/test-ph3.R`,
`examples/config_ph3.yml`,
`docs/CONFIGURATION.md`, `docs/PYTHON_INTERFACE.md`, `README.md`, and `NEWS.md`.

Implementation source read: the complete handoff and workspace `AGENTS.md`;
the three authoritative owner/plan/audit records; the Slice 1 implementation
record; complete `python/export_contract.py`,
`python/export_flowjo_populations.py`, and `tools/flowjo-orchestration.R`;
the modified R/configuration/documentation files; and the relevant existing R
and Python tests. The complete user-supplied 2026-08-28 local-verification log
was read from its explicit Codex attachment path. No experimental input
directory or experimental file was opened.

Generated as source documentation/tests only:
`tests/testthat/test-ph3-input-containment.R` and this implementation record.
No analysis output, report, image, cache, archive, or binary was generated.

The following pre-existing ignored residue was supplied as an exact inventory
for handoff. It is outside the Slice 2 diff and was not opened, read, modified,
or generated by the implementer:

- `.R-library/`
- `PROMPT_WHOLE_DATASET_ORACLE_AUDIT.md`
- `facspseudocolor.Rcheck/`
- `facspseudocolor_0.1.0.9000.tar.gz`
- `inst/.DS_Store`
- `inst/quarto/facs_configurator.html`
- `inst/quarto/facs_configurator_files/`
- `local_facs_assistant/`
- `python/__pycache__/`
- `tests/python/__pycache__/`
- `tools/__pycache__/`

Sample mappings, exclusions, experimental gates, thresholds, normalization,
statistics, and biological claims changed: none. Input acceptance is stricter
for production PH3 by the approved identity/provenance contract. Legacy PH3
normalization and numerical semantics are unchanged.

## Local verification

Agents did not execute R, Python, tests, package commands, renders, analyses,
or benchmarks. **NOT RUN — local execution required.**

Working directory:
`/Users/sansamc/OMRF Dropbox/Chris Sansam/Synchd/Organized Dropbox/Projects and Data/Projects/Cowork_Projects/2026_Analyze_FACS_Data/facs_pseudocolor_workflow`

Prerequisites: R >= 4.2; declared package dependencies including `jsonlite` and
`openssl`; `testthat`, `withr`, and `devtools`. No experimental input is needed
for the focused tests.

1. Run `Rscript -e 'devtools::test(filter = "ph3-input-containment", stop_on_failure = TRUE)'`.
   Expected: focused deterministic SYNTHETIC tests pass with zero failures or
   errors; production retains normalized data and containment with an empty
   quantitation list; only automatically removed session-temporary files are
   created.
2. Run `Rscript -e 'devtools::test(filter = "python-boundary", stop_on_failure = TRUE)'`.
   Expected: the static Python/R boundary tests pass with zero failures or
   errors, no external process launch, and no generated files.
3. Run `Rscript -e 'devtools::test(stop_on_failure = TRUE)'`.
   Expected: the full suite passes with zero failures or errors and no
   project-tree outputs. Explicit legacy PH3 still produces its historical
   quantitation, while invalid duplicate/missing/outside legacy event indices
   fail without claiming validated containment.
4. Set
   `ph3_repo="/Users/sansamc/OMRF Dropbox/Chris Sansam/Synchd/Organized Dropbox/Projects and Data/Projects/Cowork_Projects/2026_Analyze_FACS_Data/facs_pseudocolor_workflow"` and
   `ph3_check_dir=$(mktemp -d)`, then run
   `(cd "$ph3_check_dir" && R CMD build "$ph3_repo")`, then
   `R CMD check --no-manual "$ph3_check_dir"/facspseudocolor_*.tar.gz --output="$ph3_check_dir"`.
   Expected: `Status: OK`; only the temporary tarball/check directory is
   generated outside the repository.
5. Run `git diff --check`, `git status --short`, and
   `git status --short --ignored`.
   Expected: no whitespace errors; the non-ignored change is limited to the
   files listed above; and ignored items match the exact pre-existing residue
   inventory above. Those ignored paths are outside the Slice 2 diff and must
   not be mistaken for outputs of this implementation.

Stop on any test/check failure or warning, any acceptance of a missing/mutated
manifest or artifact, any legacy result labeled validated, any experimental
fixture selection, or any unexpected generated project file. Report only the
command, `PASS` or `FAIL`, failed test names, a relevant error excerpt, and
unexpected generated files; do not return large successful logs.

## User local verification and corrective follow-up

The user supplied the complete log from a local run on 2026-08-28. The focused
`ph3-input-containment` command failed with 20 failures/errors and 29 passing
expectations. The full test run and package-check test run also failed. Package
build completed successfully; `R CMD check --no-manual` ended with **1 ERROR,
1 NOTE**, with its test tail reporting 27 failures, 1 warning, and 408 passing
expectations. The NOTE was the unqualified `file_test` reference. The final
`git diff --check` in the supplied command sequence apparently produced no
output. These are preserved as failed pre-correction results and are not treated
as verification of the corrected tree.

The supplied failures identified exact-type mismatches in containment
expectations; reason-code ordering in scoped error regexes; an incomplete legacy
configuration fixture; overbroad collapsed-source reconstruction regexes;
`openssl` digest objects leaking their `hash`/`sha256` classes into SYNTHETIC
manifest values; scalar assignment into a zero-row fixture; nonexistent
relative operation paths containing lexical parent traversal; unqualified
`file_test`; and legacy high-level test capture that did not retain the final
quantified analysis object. The implementer corrected those items without
weakening production validation: SHA helpers now return plain character
scalars, zero-row bindings use length-zero vectors, relative operation paths
reject lexical `..` while preserving harmless `.` components,
`utils::file_test` is explicit, identity-reconstruction scans target actual
assignments/coercions, the legacy profile/ranges are complete, and the legacy
analysis value is captured independently of the warning expectation. The
SYNTHETIC acquisition, population, artifact, and count-report lists were also
inspected for duplicated named fields; no duplicated `acquisition_id` or similar
field remains.

After these corrections, agents did not rerun R, Python, focused/full tests,
package build/check, renders, analyses, or benchmarks. Post-correction status is
**NOT RUN — local rerun required**; the user-provided failed results above were
the latest executed evidence until the second stopped run recorded below.

The user supplied a second stopped-run log later on 2026-08-28. Its focused
`ph3-input-containment` run reported 18 failures/errors and 38 passing
expectations; all but the macOS path-spelling expectation arose because the
SYNTHETIC production helper retained `minimal_config()`'s relative `data_dir`
while callers intentionally omitted a `data_dir` override. The separate
`python-boundary` run passed all 15 expectations. The full suite passed the
16 `analysis` expectations and was manually interrupted during
`appearance-and-artifacts` after 5 expectations, so it did not complete. The
package build failed before reading the repository because the handoff command
expanded `$PWD` only after changing into the temporary directory; it looked for
the temporary directory's nonexistent `DESCRIPTION`. Consequently, no source
tarball existed and `R CMD check` did not run. `git diff --check` again produced
no reported output, the tracked/untracked status was unchanged, and the ignored
inventory exactly matched the pre-existing list above.

The follow-up correction is confined to SYNTHETIC test setup and this handoff:
`synthetic_production_config()` now normalizes every explicit operation
directory, sets an absolute `data_dir` to the first operation's parent, and
continues to support multiple explicit operation directories. The harmless-dot
path expectation now normalizes the existing base before appending its
nonexistent child, avoiding macOS `/var` versus `/private/var` spelling. The
build instruction now records the canonical repository in `ph3_repo` before
changing directories. Production path-resolution behavior was not weakened or
otherwise changed.

After this follow-up correction, agents again did not run project code, R,
Python, focused/full tests, package build/check, renders, analyses, or
benchmarks. In a third user run on 2026-08-28, the corrected focused
`ph3-input-containment` command passed all 117 expectations with zero failures,
warnings, or skips. The full suite then reported **2 failures/errors, 0 warnings,
and 517 passing expectations**. Both failures were confined to the legacy
`test-ph3.R` assertion: assigning the return value of `expect_warning()` stored
the warning expectation result instead of the validated input report, so its
`ph3_containment` attribute appeared `NULL`. The test now captures and muffles
the single expected warning independently while retaining the actual report;
production behavior is unchanged.

After this test-only correction, agents did not rerun project code or tests.
The focused 117/117 PASS is preserved as executed evidence for the preceding
tree, while the full-suite result remains FAIL for that tree. Current status is
**NOT RUN after the final correction — local rerun required**. The earlier
`python-boundary` result remains 15 passing expectations; package build/check
was not included in the third log.

The user subsequently reported that the targeted corrected legacy `ph3` test
context completed with zero failures, errors, or warnings. The complete full
suite was then reported as passing. Package build also succeeded. `R CMD check`
completed all checks successfully until installed-package tests, where the
static Slice 2 source scan tried to open `../../R/input.R` from the installed
test location. The check ended with **Status: 1 ERROR** and a test summary of
**1 failure, 1 warning, and 510 passing expectations**; the warning and failure
both arose from that same unavailable source-file connection.

The test-only correction removes the source-tree path dependency without a skip
or bypass. It obtains the package namespace's own bindings with
`inherits = FALSE`, deparses local function bodies, applies identity-
reconstruction checks separately to every local package function, and applies
the full unsafe-API assertions to the production functions corresponding to the
audited Slice 2 R sources. Imported functions are not scanned. Production code and
validation behavior are unchanged.

After this installed-package test correction, agents did not rerun project code,
tests, package build/check, renders, analyses, or benchmarks. The user-reported
full-suite PASS and build success are preserved for the preceding tree; current
package-check status is **NOT RUN after the final correction — local rerun
required**.

The user then rebuilt `facspseudocolor_0.1.0.9000.tar.gz` from the final
corrected tree and ran `R CMD check --no-manual` under R 4.5.0 on
`aarch64-apple-darwin20`. Installation, loading, dependencies, code checks,
documentation, examples, and installed-package tests all completed successfully.
Final package-check result: **Status: OK** with no errors, warnings, or notes.
The tarball and `facspseudocolor.Rcheck` directory were generated only beneath
the reported temporary check directory, outside the repository.

## Assumptions and unresolved uncertainty

The merged Slice 1 schema, ordered identity columns, exact method metadata,
canonical source-index form, and digest metadata remain authoritative. SHA-256
is checked over exact bytes; R does not rewrite or lossily replace the finalized
manifest. Direct-index semantics remain conditional on the explicit verified
FlowKit version recorded by Slice 1, and R requires recorded installed and
supported versions to agree. Absolute operation directories are permitted only
because they are explicit configuration; relative directories must remain
beneath `data_dir`. Legacy files do not contain trustworthy acquisition/sample
identity, so their structured rows carry configured prefixes and explicit `NA`
for unavailable proof instead of inventing identifiers.

## Independent review and resolution

The initial independent reviews found that production analysis still called
legacy biological quantitation; legacy conditional event-index safety had been
removed; operation IDs were not unique across directories; source FCS
references and whitespace-only metadata were too permissive; and the static
unsafe-API test and local verification handoff were incomplete. The implementer
resolved these findings by withholding `quantify_ph3()` only for production;
restoring legacy missing/duplicate/outside-Single-Cells index failures without
claiming containment proof; requiring unique operation IDs; validating confined
portable relative FCS references and trimmed-nonblank metadata; scanning all
changed production R files for network, telemetry, shell, broad mutation, and
destructive APIs; and adding the static Python-boundary verification command.

A subsequent scientific-integrity re-review identified that callers could still
invoke exported `quantify_ph3()` manually on a production-profile analysis. The
implementer resolved that blocking public-API bypass with a direct legacy-profile
guard and a focused regression assertion that production quantitation remains
empty after the rejected call.

The artifact/security re-review found that the first R-side source-FCS reference
check rejected harmless components accepted by the merged Slice 1 `Path`
contract and that the pre-existing ignored residue needed an explicit handoff
inventory. The implementer aligned the check and SYNTHETIC cases to normalize
empty/`.` components while still rejecting parent traversal and all
absolute/UNC/drive forms, and recorded the supplied residue inventory without
opening or modifying it.

After those resolutions, the final independent scientific-integrity and
artifact/security re-reviews were clean: neither reviewer reported any blocking,
important, or advisory finding. The implementer did not self-certify either
independent review role.

The focused deterministic `SYNTHETIC` matrix prepared for local execution
covers valid independent/overlapping children; less/equal/excess/unmatched and
cross-acquisition identities; exact canonical JSON; full and binding digests;
exact method/digest/acquisition subcontracts; confined paths; extra and missing
ledger fields; hashes, sizes, schemas, channels, row bindings, canonical event
indices, duplicates, counts, 10/2/0 finite-event minima, intentional empty
proof, missing gates/artifacts, transient verified-event consumption,
structured provenance, production-withheld quantitation, unique multi-operation
IDs, confined source references, whitespace-only metadata rejection, explicit
per-child legacy results and legacy index safety, static unsafe-API absence,
rejection of manual production `quantify_ph3()` calls, and unchanged EdU/POI
provenance shape.

Static check performed by the implementer: `git diff --check` — **PASS** (no
output). This is not project-code or test execution. All R/Python tests,
package build/check, renders, analyses, and benchmarks remain **NOT RUN — local
execution required**.
