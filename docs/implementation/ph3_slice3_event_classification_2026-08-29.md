# pH3 Slice 3 eligibility and shared event classification

**Status:** VERIFIED — USER-REPORTED LOCAL PASS; WINDOWS CI PENDING

## Scope and provenance

Canonical source:
`/Users/sansamc/OMRF Dropbox/Chris Sansam/Synchd/Organized Dropbox/Projects and Data/Projects/Cowork_Projects/2026_Analyze_FACS_Data/facs_pseudocolor_workflow`

Branch: `feature/ph3-event-classification`. Starting commit:
`4fba05aad354c768a99b987ef3aed333cdf30647` with no tracked or non-ignored
untracked changes. The exact pre-existing ignored residue inventory carried
forward from the Slice 2 record was `.R-library/`,
`PROMPT_WHOLE_DATASET_ORACLE_AUDIT.md`, `facspseudocolor.Rcheck/`,
`facspseudocolor_0.1.0.9000.tar.gz`, `inst/.DS_Store`,
`inst/quarto/facs_configurator.html`, `inst/quarto/facs_configurator_files/`,
`local_facs_assistant/`, `python/__pycache__/`,
`tests/python/__pycache__/`, and `tools/__pycache__/`. This inventory is cited
from the authoritative Slice 2 implementation record; those paths were not
opened, read, modified, or generated during Slice 3. This slice consumes only Slice 2-validated production identity,
containment, event tables, and immutable manifest linkage. It creates one
classification per acquisition and retains it only at
`analysis$normalized_data[[acquisition]]$ph3_event_classification`.

Production biological metrics, result tables, CSVs, sensitivity calculations,
aggregation, plots, reports, and geometry rendering remain withheld. Explicit
legacy quantitation remains on its unchanged count-only path. EdU and POI code
paths are unchanged.

Experimental inputs read or modified: none. No FCS, WSP, Prism, exported CSV,
RDS, figure, report, archive, or historical analysis artifact was opened,
hashed, copied, modified, or generated.

## Implemented classification

For every validated Single Cells row, the retained table preserves the parent
row number and exact direct `event_identity`, identity source/method/occurrence,
acquisition and prefix, immutable operation/manifest linkage, and validated G1
and pH3-positive memberships. Membership is joined by exact identity; parent
row order is never used as identity or as a membership fallback. Any changed
parent count, duplicate/missing/out-of-parent identity, acquisition, operation,
manifest binding, identity method/source, containment status, or manifest link
is fatal with a stable reason code. Each current child is also compared to the
original Slice 2-validated child table and its population-specific row, unique-
identity, and matched counts; dropping even one valid child after containment
is therefore fatal. Parent and child `event_index` values must remain canonical
nonnegative decimal tokens, `event_identity` must remain exactly
`<acquisition_id>:event_index:<event_index>`, and identity source, method, and
occurrence must still match the verified manifest and containment evidence.
For the approved direct-index method, `event_identity_base` is the exact direct
identity token itself: there is no separately suffixed composite base, and
`identity_occurrence` is the exporter-proven direct occurrence value `"1"`.

The existing G1-anchor DNA normalization is unchanged. Classification then:

1. marks normalized-DNA finiteness;
2. assigns inclusive 2to4N eligibility `[b0,b5]`;
3. partitions eligible rows into sub-4N `[b0,b4)` and 4N `[b4,b5]`;
4. uses the existing configured phase assignment and boundary ownership;
5. records one eligibility reason (`none`, `identity_invalid`,
   `dna_nonfinite`, `below_b0`, or `above_b5`) and one phase-Unassigned reason
   (`none` or `configured_phase_gap`);
6. validates exact row, exclusion, partition, phase, and membership
   reconciliation before retaining the table.

Production inputs have already passed direct identity validation, so
`identity_invalid` is structurally available but any later identity mutation
fails closed rather than being silently classified. `b0`, `b4`, and `b5`
ownership is exact. Raw target intensity is retained as raw display metadata,
with `target_display_finite`; target finiteness does not affect eligibility or
FlowJo-positive membership. Geometry remains `not_requested`, and the declared
display transform is the raw identity transform only.

The configuration identity digest is SHA-256 over canonical JSON of the
validated analysis configuration after explicitly excluding the location-only
keys `data_dir`, `output_pdf`, `output_png`, and
`ph3_export_operation_dirs`. Their exact locations remain separately preserved
in the existing validated input report, manifest references, and analysis
provenance; they do not change scientific configuration identity when only
Windows/macOS/Linux path spelling changes. The analysis ID is SHA-256 over that digest plus the
sorted immutable export-operation IDs and full-manifest SHA-256 values. These
keys are deterministic across supported operating systems and do not use local
row order or filenames as scientific identity.

Method identifiers retained are the approved positivity ID from the manifest,
`identity_validated_finite_dna_b0_b5_inclusive_v1`,
`configured_b0_b4_left_closed_right_open_v1`,
`configured_b4_b5_closed_v1`, the shared configured-boundary interval ID,
containment method ID, classification schema
`ph3-event-classification-1.0.0`, and withheld future output schema
`ph3-1.0.0`.

## Files

Read before implementation: workspace `AGENTS.md`; the complete Slice 3
handoff; `ph3_scientific_owner_approval_2026-08-24.md`;
`ph3_implementation_plan_and_test_spec_5de77a9.md`; Slice 1 and Slice 2
implementation records; `R/ph3.R`, `R/analysis.R`, focused portions of
`R/config.R` and `R/input.R`; `tests/testthat/helper-current-workflow.R`,
`tests/testthat/test-ph3.R`, and
`tests/testthat/test-ph3-input-containment.R`.

Modified: `R/ph3.R`, `R/analysis.R`, and
`tests/testthat/test-ph3-input-containment.R`.

Generated as source/test documentation only:
`tests/testthat/test-ph3-event-classification.R` and this implementation
record. No generated analysis artifact, cache, archive, binary, rendered
document, or package output was created.

Sample mappings, experimental gates, thresholds, G1 normalization, statistics,
or biological claims changed: none. The owner-approved eligibility,
exclusion-reason, sub-4N/4N, and configured-phase classifications are newly
implemented exactly as approved. No numerical biological result is produced.

## SYNTHETIC tests prepared

The new focused test uses only in-memory values explicitly labeled SYNTHETIC
and asserts exact identity/row preservation; `b0`, just above `b0`, just below
`b4`, exact `b4`, just below `b5`, exact `b5`; below/above-span; `NA`, `NaN`,
positive and negative infinity; exclusion exclusivity; sub-4N/4N exclusivity
and reconciliation; configured phase ownership and explicit gap handling;
target-finiteness independence; exact FlowJo-positive membership; fatal
identity/containment inconsistencies; fatal loss of a validated child row;
canonical direct-index and identity-provenance mutation rejection; and equal
configuration identity across Windows/POSIX location spelling while a
scientific boundary change alters the digest. The existing Slice 2 high-level
SYNTHETIC production test now also requires retained classification with exact
parent identity order while confirming biological quantitation is empty and
manual production quantitation remains rejected.

## Local verification contract

Agents did not execute R, Python, project code, tests, builds, package checks,
renders, analyses, snapshots, or benchmarks. **NOT RUN — local execution
required.**

Working directory for all repository commands:
`/Users/sansamc/OMRF Dropbox/Chris Sansam/Synchd/Organized Dropbox/Projects and Data/Projects/Cowork_Projects/2026_Analyze_FACS_Data/facs_pseudocolor_workflow`

Prerequisites: R >= 4.2; declared package dependencies including `openssl` and
`jsonlite`; `testthat`, `withr`, and `devtools`; supported Python only for the
repository boundary checks. Dependencies must already be installed. No
experimental input is needed or permitted.

1. Slice 3: `Rscript -e 'devtools::test(filter = "ph3-event-classification", stop_on_failure = TRUE)'`.
2. Slice 2 regression: `Rscript -e 'devtools::test(filter = "ph3-input-containment", stop_on_failure = TRUE)'`.
3. Existing pH3 regression: `Rscript -e 'devtools::test(filter = "^ph3$", stop_on_failure = TRUE)'`.
4. Export/R boundary: `Rscript -e 'devtools::test(filter = "python-boundary", stop_on_failure = TRUE)'`.
5. Repository-only Python contract: `python3 -m unittest -v tests/python/test_export_contract.py`.
6. Full R suite: `Rscript -e 'devtools::test(stop_on_failure = TRUE)'`.
7. Create a temporary directory outside the repository with
   `ph3_check_dir=$(mktemp -d)`, then run
   `(cd "$ph3_check_dir" && R CMD build "/Users/sansamc/OMRF Dropbox/Chris Sansam/Synchd/Organized Dropbox/Projects and Data/Projects/Cowork_Projects/2026_Analyze_FACS_Data/facs_pseudocolor_workflow")`, followed by
   `R CMD check --no-manual "$ph3_check_dir"/facspseudocolor_*.tar.gz --output="$ph3_check_dir"`.
8. Run `git diff --check`, `git status --short --untracked-files=all`, and
   `git status --short --ignored` from the canonical repository.

Expected success: focused and full tests have zero failures/errors and only
explicitly asserted legacy warnings; Python reports `OK`; package check reports
`Status: OK` with zero errors, warnings, or notes; the tarball and check
directory exist only beneath the temporary directory; repository status shows
only the files listed above and no cache, report, image, RDS, tarball, check
directory, or other unexpected artifact.

Stop on any identity/provenance/containment mismatch; lost, duplicated, or
reordered parent identity; changed membership; eligibility/partition/phase
reconciliation difference; changed legacy number; nonempty production
biological quantitation; EdU/POI regression; network attempt; warning beyond an
asserted legacy warning; test/check failure; check note requiring review; or
unexpected repository artifact. Report only the command, `PASS` or `FAIL`,
failed test names, a short relevant error excerpt, and unexpected generated
files—not a large successful log.

### First user-run local verification and correction

The user supplied a local verification transcript after the initial reviewed
implementation. The focused `ph3-event-classification` run reported **FAIL 6,
WARN 0, PASS 1**. Five cases stopped prematurely with
`parent_identity_inconsistency`; the sixth showed unequal configuration digests
for location-only Windows/POSIX spelling. The root causes were confined to the
new SYNTHETIC fixture and configuration-identity helper:

- the fixture copied the normalized list names (`data`, `g1`,
  `ph3_positive`) into `validated_inputs`, while the actual Slice 1/2 contract
  and production plumbing use `complete`, `g1`, `ph3_positive`; therefore the
  fixture's validated parent lookup was absent before intended assertions;
- assigning `NULL` through a multi-name list subset did not reliably implement
  the intended exact omission contract, leaving location-only values able to
  affect canonical configuration bytes.

The correction constructs the fixture with the actual `complete` parent key.
It replaces the ambiguous list mutation with an explicit positive selection:
`value[setdiff(names(value), location_only)]`. Exactly the four documented
location-only fields are excluded; all other configuration fields remain in
canonical identity, while export operation IDs and full manifest digests remain
separately bound into `analysis_id`. Existing negative fixtures continue to
mutate their validated comparison table only when the intended assertion is
canonical direct provenance, ensuring child loss reaches
`ph3_positive_identity_inconsistency`, noncanonical index/profile reaches
`*_direct_identity_mismatch`, and invalid containment reaches
`invalid_containment_context`.

Other results in the supplied transcript were:

- Slice 2 `ph3-input-containment`: **PASS 493**, zero failures or warnings;
- existing `ph3`: interrupted after 39 passing expectations, **NOT
  COMPLETED**;
- `python-boundary`: **PASS 15**, zero failures or warnings;
- Python export contract: **PASS 18**;
- full R suite: interrupted during the first context, **NOT COMPLETED**;
- package build and check: **NOT RUN**.

Agents did not execute any verification after these corrections. The corrected
tree is **NOT RUN — local execution required**. The commands and stop criteria
above remain the required rerun contract; interrupted commands must be rerun to
completion, and build/check must be performed only after all focused and full
tests complete successfully.

### Second user-run local verification and correction

The user supplied a second transcript. The focused
`ph3-event-classification` run reported **FAIL 2, WARN 0, PASS 23**. The intended
noncanonical-direct-index case instead produced `g1_identity_inconsistency`,
and changing a scientific boundary did not change the configuration digest.
The Slice 2 run was manually interrupted after 419 passing expectations and is
**NOT COMPLETED**. The existing pH3 run was manually interrupted after 20
passing expectations and is **NOT COMPLETED**. `python-boundary` reported
**PASS 15**, and the Python export contract reported **PASS 18**. The full R
suite was manually interrupted and is **NOT COMPLETED**. Package build/check
were **NOT RUN**.

Static diagnosis found two exact causes:

- `attributes(value) <- NULL` removed the names from the configuration list.
  The subsequent name-based positive selection therefore selected an empty
  list, erasing both location and scientific semantics from the digest.
- the negative fixture changed both `event_index` and `event_identity` to a new
  identity absent from the parent. Exact child containment correctly failed
  before canonical direct-index validation, so the observed stable reason was
  scientifically appropriate for that compound mutation but not the intended
  isolated assertion.

The corrected digest begins with `as.list(config)`, which removes the class
without removing list names, then uses exact `setdiff()` selection of only the
four documented location-only keys. Tests now require digest changes for
multiple scientific fields (`g2m_x_range`, `g1_x_range`, and `dna_channel`) in
addition to requiring equality across location spelling. The direct-index
negative case now changes only `event_index` to the noncanonical token `"01"`,
leaves the approved parent-linked `event_identity` unchanged, and synchronizes
the validated child. Containment remains valid and the isolated canonical
relation failure must reach `g1_direct_identity_mismatch`.

Agents did not execute any verification after this second correction. The
current corrected tree remains **NOT RUN — local execution required**. Rerun
the exact commands above, allowing long-running focused, regression, and full
suite commands to finish unless they report an actual failure or other stated
stop condition. A progress counter or temporarily quiet context is not a
failure and is not a reason to interrupt. Build/check remain last, after every
required test completes successfully.

### Final user-run local verification

After the second correction, the user reported that the focused Slice 3
classification test passed. The user then reported that every remaining
required command completed successfully: the Slice 2 containment regression,
existing pH3 regression, Python boundary tests, Python export-contract tests,
the full R test suite, package build, and `R CMD check --no-manual`. The package
check result was reported as `Status: OK`, and `git diff --check` produced no
output.

The final repository inventory contained only the five intended Slice 3 source,
test, and implementation-record paths listed above. The separately reported
ignored paths matched the inherited pre-existing residue inventory; no
unexpected repository output was reported. Package build/check outputs were
created only in the user-selected temporary directory outside the repository.

These commands were executed by the user, not by an agent. The user reported
their outcome as pass/OK without supplying per-test assertion counts in the
final concise result. GitHub Actions, including the Windows package-check job,
remain pending until this reviewed change is committed and pushed.

## Assumptions and unresolved uncertainty

The merged Slice 1/2 direct-identity schema, FlowKit-version attestation,
containment rows, and manifest hashes remain authoritative. Exact configured
phase contiguity is still enforced upstream; explicit `configured_phase_gap`
handling is retained as fail-visible defense. The approved classification is
now the intended sole future source for Slice 4 metric predicates, but Slice 4
does not yet exist. Local execution remains required; the completed independent
reviews are recorded below.

## Independent-review resolution

The first scientific-integrity review found a blocking post-containment child-
loss gap and an important direct-identity provenance gap. Resolution passes the
original Slice 2-validated parent and child tables separately into
classification, compares exact identity vectors and population-specific
containment counts, and rejects canonical-index, identity-source, method,
version, occurrence, profile, acquisition, operation, or binding mutation.
Focused SYNTHETIC loss and mutation cases were added. The exact semantics of
the direct identity base and occurrence are documented above.

The first artifact/security review found that configuration and analysis keys
could vary with machine-specific input/output path spelling. Resolution defines
the explicit location-independent configuration identity contract above,
retains locations through their existing provenance owners, and adds a
Windows/POSIX spelling regression assertion plus a scientific-boundary control.
No network, telemetry, unsafe execution, destructive I/O, fallback, or generated
artifact behavior was added. Both independent read-only re-reviews of that
pre-user-verification-failure tree found no remaining blocking, important, or
advisory findings.

After the user-reported fixture and portable-digest failures were corrected as
recorded above, renewed independent scientific-integrity and artifact/security
reviews inspected that first corrected, pre-second-rerun diff and relevant
surrounding code. Both reviews were clean, with no blocking, important, or
advisory findings.

After the second user-reported failures and their focused corrections, renewed
independent scientific-integrity and artifact/security reviews inspected the
second corrected complete diff and relevant surrounding code. Both renewed
reviews were clean, with no blocking, important, or advisory findings. These
review outcomes do not change execution status: the current corrected tree
remains **NOT RUN — local execution required**.
