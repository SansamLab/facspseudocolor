# pH3 Slice 5 biological-replicate and condition aggregation

**Status:** IMPLEMENTED, INDEPENDENTLY REVIEWED, AND USER-VERIFIED

The scientific owner resolved every stop-gate item in
`ph3_slice5_aggregation_scientific_owner_decision_2026-08-29.md` on 2026-08-29.
That confirmed supplement is authoritative and remains unchanged. It limits
numerical aggregation to the Slice 4 metric and phase-prevalence tables,
defines the exact four output tables and column order, confirms unweighted
finite-value means at both levels, defines counts/statuses/SD/SEM/provenance,
requires manifest-exact completeness and fail-closed validation, and fixes
deterministic ordering. Slice 4 eligibility and sensitivity QC remain attached
unchanged and acquisition-only.

## Scope and provenance

Canonical source:
`/Users/sansamc/OMRF Dropbox/Chris Sansam/Synchd/Organized Dropbox/Projects and Data/Projects/Cowork_Projects/2026_Analyze_FACS_Data/facs_pseudocolor_workflow`

Exact baseline: `bb48c0347822570476c6d335b13b36c5c195f2ff`
(`Merge pull request #20 from SansamLab/feature/ph3-acquisition-metrics-sensitivity`).
Working branch: `feature/ph3-replicate-condition-aggregation`. The tracked and
non-ignored tree was clean before this record was created. No commit, push,
pull request, merge, or branch deletion occurred.

Experimental inputs read or modified: none. Pre-existing ignored residue was
not opened, modified, deleted, staged, or inventoried anew.

## Initial contract extraction and stop decision (resolved)

The confirmed owner record approves unweighted arithmetic means of finite
technical-acquisition percentages within each explicit biological-replicate /
condition group; biological-replicate means are equal-weight independent points
in condition summaries; acquisition values and `technical_n`, aggregation
method, and heterogeneity QC are retained; event pooling and denominator
weighting are prohibited; and an all-nonfinite technical group produces a
missing value with structured status. The plan additionally specifies the
technical example `{10, NA, 30}` as mean 20 with total/finite/missing counts and
min/max/SD/range, and the two-bioreplicate example `{20, 40}` as condition mean
30. Approved identifiers are
`aggregation_method_id = "unweighted_technical_percentage_mean_within_biorep_condition_v1"`
and `output_schema_version = "ph3-1.0.0"`.

The authoritative Slice 4 numerical inputs are attached directly below
`analysis$quantitation` as `ph3_metrics_acquisition`,
`ph3_phase_prevalence`, `ph3_event_eligibility_qc`, and
`ph3_4n_boundary_sensitivity_qc`. Their rows retain explicit acquisition,
sample, prefix, condition, replicate, technical-replicate, schema, method,
configuration, operation, and manifest provenance. Slice 4 expressly attaches
no aggregate or condition output.

The required 15-point gate resolves as follows:

1. Aggregation source: all four Slice 4 tables are authoritative inputs, but
   the records do not say which percentage-bearing tables are aggregated.
2. Biological-replicate output: `ph3_metrics_biological_replicate` is approved
   by name; its exact columns and whether phase/sensitivity results share or
   receive separate tables are unspecified.
3. Condition output: no exact condition-level table name or schema is approved.
4. Experimental units: explicit technical acquisitions within an explicit
   biological-replicate/condition group; biological replicates at condition
   level.
5. Multiple technical acquisitions: unweighted arithmetic mean of finite
   percentages, never pooled or denominator weighted.
6. Ratio handling: percentages are averaged; counts are not approved for
   recomputation of the aggregate value. The disposition of aggregate count
   columns is unspecified.
7. Phase prevalence: the records do not explicitly state whether Slice 4 phase
   rows use the same aggregation rule or remain acquisition-only.
8. Eligibility and sensitivity QC: acquisition values are retained, but the
   records do not specify whether either table is aggregated or summarized.
9. Condition summaries: equal-weight biological-replicate points and an
   arithmetic-mean example are approved; exact required descriptive fields and
   undefined-value accounting are unspecified.
10. Inferential statistics: none are approved for Slice 5.
11. Missing/excluded/zero-denominator behavior: partial and all-missing
    technical-value behavior is partly specified; missing acquisition rows must
    fail. Excluded-replicate representation and condition-level undefined
    behavior are unspecified.
12. Ordering: deterministic ordering is required, but the exact key precedence
    for aggregate tables is unspecified.
13. Minimum replicates: no minimum biological-replicate rule is approved.
14. Provenance: the method/schema IDs above and complete Slice 4 provenance are
    required, but the exact aggregate-row provenance columns and lossless
    technical-value representation remain unspecified.
15. Attachment: Slice 4 tables are at `analysis$quantitation`; the exact names
    and attachment positions for all Slice 5 biological-replicate and condition
    outputs are not fully specified. CSVs, plots, reports, geometry,
    compatibility aliases, inferential statistics, and Slice 6+ work are
    prohibited.

Sources: `ph3_scientific_owner_approval_2026-08-24.md` (PH3-D09 and stable
schema section), `ph3_implementation_plan_and_test_spec_5de77a9.md` (Slice 5
and focused deterministic test matrix), and final Slice 3 and Slice 4
implementation records.

## Narrow scientific-owner decision memo (resolved)

Before implementation, approve one exact Slice 5 schema answering:

1. Which percentage rows are aggregated: primary metrics only; primary plus
   phase prevalence; or primary, phase prevalence, and both sensitivity
   percentages. State whether eligibility QC remains acquisition-only.
2. The exact biological-replicate and condition table names, column order,
   unique keys, and attachment names beneath `analysis$quantitation`, including
   whether phase/sensitivity aggregates use separate tables.
3. Whether numerator/denominator/event-count fields in aggregate rows are
   omitted, retained only as lossless source-value lists, or summarized in an
   explicitly named non-estimand form. Summed-count recomputation would conflict
   with the approved unweighted-percentage estimand unless separately approved.
4. Exact biological-replicate heterogeneity fields and definitions, including
   SD behavior for one finite value, range representation, and serialization of
   individual acquisition IDs/statuses/values.
5. Exact condition descriptive fields and definitions (at minimum whether the
   approved arithmetic mean is accompanied by SD, range, counts, or other
   summaries), plus partial/all-undefined biological-replicate handling.
6. Missing and excluded biological-replicate representation, minimum replicate
   requirements, and whether a one-bioreplicate condition is valid descriptive
   output or a structured undefined result.
7. Exact deterministic ordering keys for every new table.

Scientifically meaningful alternatives affect which estimands exist, whether
undefined acquisitions or replicates influence summaries, whether QC is
misrepresented as biology, and whether downstream code can mistake summed
counts for the approved unweighted mean. No conventional default was selected.

## Files and changes

Read: workspace `AGENTS.md`; complete Slice 5 handoff; complete scientific-owner
approval; complete implementation plan/test specification; final Slice 3 and
Slice 4 records; focused `R/ph3.R`, `R/analysis.R`,
`tests/testthat/test-ph3-acquisition-metrics.R`,
`tests/testthat/test-ph3-input-containment.R`, and
`tests/testthat/helper-current-workflow.R`.

Modified: none.

Generated as source documentation only: this implementation record.

Sample mappings, exclusions, gates, thresholds, normalization, aggregation
implementation, statistics, and biological claims changed: none. Tests
prepared: none, because their expected schemas/formulas would encode decisions
that are not approved.

## Verification and review status

**NOT RUN — local execution required.** No R, Python, project code, tests,
builds, checks, renders, analyses, snapshots, or benchmarks were executed.
Independent scientific-integrity and artifact/security reviews were not
requested because production implementation stopped at the mandatory contract
gate. After the owner supplies the decisions above and implementation/reviews
are complete, the required ordered human sequence is: focused Slice 5; Slice 4;
Slice 3; Slice 2; legacy pH3; Python boundary; repository-only Python export
contract; full R suite; external-temporary `R CMD build` and
`R CMD check --no-manual`; then `git diff --check` and exact
tracked/untracked/ignored inventory. Exact executable commands must be finalized
against the approved test filename and implemented files; none should be run
against this intentionally stopped documentation-only state.

At the initial stop, the merged Slice 1–4 schemas, manifest mappings,
provenance, and attachment locations were treated as authoritative. The seven
questions above were subsequently resolved by the confirmed supplement.

## Implemented contract after owner confirmation

The confirmed supplement resolves the historical stop memo above. Slice 5 now
validates all four attached Slice 4 tables against the active explicit manifest
and against one another. Missing, extra, duplicated, reordered, remapped,
mixed-analysis, mixed-schema, mixed-method, or provenance-inconsistent rows
fail closed. Aggregation uses only `value_percent` from
`ph3_metrics_acquisition` and `ph3_phase_prevalence`. The eligibility and
sensitivity QC tables remain unchanged and acquisition-only.

Exactly four tables are appended after the four Slice 4 tables:

- `ph3_metrics_biological_replicate`;
- `ph3_phase_prevalence_biological_replicate`;
- `ph3_metrics_condition_summary`;
- `ph3_phase_prevalence_condition_summary`.

Within explicit condition × replicate manifest groups, finite technical values
receive an unweighted arithmetic mean. Within condition, finite biological-
replicate means receive an unweighted arithmetic mean. Acquisition counts are
never pooled or copied into aggregate tables. Condition SD is sample SD; SEM is
SD divided by the square root of the finite biological-replicate count; both
are missing below two finite replicates. Exact total/finite/undefined counts and
the approved `ok`, `ok_partial_undefined`,
`undefined_no_finite_values`, and `ok_single_biological_replicate` statuses are
implemented. Zero remains numeric zero. Provenance arrays are canonical JSON
arrays of unique lexicographically sorted strings. Fixed metric/phase ordering
and the exact approved column order are enforced.

No event classification, positivity, eligibility, phase membership, 4N/sub-4N
membership, acquisition formula, sensitivity predicate, sample mapping,
exclusion, gate, threshold, normalization, inferential statistic, biological
claim, CSV/RDS export, plot, report, geometry, GUI, compatibility alias, legacy
pH3 behavior, EdU behavior, or POI behavior was changed.

## Final implementation file inventory

Additional authoritative file read unchanged:
`docs/implementation/ph3_slice5_aggregation_scientific_owner_decision_2026-08-29.md`.

Modified: `R/ph3.R`,
`tests/testthat/test-ph3-input-containment.R`, and this implementation record.

Generated as source/test material:
`tests/testthat/test-ph3-replicate-condition-aggregation.R` and this record.
The scientific-owner decision record was supplied externally and preserved
unchanged. No generated analysis output, cache, archive, binary, report, image,
RDS, package tree, or temporary artifact was created.

The focused in-memory `SYNTHETIC` tests cover unequal-denominator technical
values 10/undefined/30 yielding 20 rather than a pooled result; prevention of
technical pseudoreplication; equal-replicate condition mean, sample SD, and SEM;
numeric zero; partial undefined, all undefined, and single-replicate statuses;
canonical provenance arrays; exact schemas and ordering; preservation of Slice
4 inputs; and fatal missing, duplicated, mixed-analysis, and remapped rows. The
high-level Slice 2 test now requires the exact ordered eight-table attachment
while preserving containment, classification, acquisition rows, legacy
rejection, EdU, and POI assertions.

Agent execution status remains **NOT RUN — local execution required**. No R,
Python, project code, test, build, check, render, analysis, snapshot, or
benchmark was executed.

## Ordered human verification commands

Working directory for commands 1–8 and 10:
`/Users/sansamc/OMRF Dropbox/Chris Sansam/Synchd/Organized Dropbox/Projects and Data/Projects/Cowork_Projects/2026_Analyze_FACS_Data/facs_pseudocolor_workflow`.

Prerequisites: R >= 4.2 with declared package dependencies plus `testthat`,
`withr`, and `devtools`; supported Python with the existing repository export-
contract dependencies. All dependencies must already be installed. No
experimental input is needed or permitted.

1. Slice 2 containment/high-level regression: `Rscript -e 'devtools::test(filter = "ph3-input-containment", stop_on_failure = TRUE)'`.
2. Slice 4: `Rscript -e 'devtools::test(filter = "ph3-acquisition-metrics", stop_on_failure = TRUE)'`.
3. Slice 5: `Rscript -e 'devtools::test(filter = "ph3-replicate-condition-aggregation", stop_on_failure = TRUE)'`.
4. Slice 3: `Rscript -e 'devtools::test(filter = "ph3-event-classification", stop_on_failure = TRUE)'`.
5. Legacy pH3: `Rscript -e 'devtools::test(filter = "^ph3$", stop_on_failure = TRUE)'`.
6. Python boundary: `Rscript -e 'devtools::test(filter = "python-boundary", stop_on_failure = TRUE)'`.
7. Repository-only Python export contract:
   `python3 -m unittest -v tests/python/test_export_contract.py`.
8. Complete R suite: `Rscript -e 'devtools::test(stop_on_failure = TRUE)'`.
9. External temporary build/check: set
   `ph3_repo="/Users/sansamc/OMRF Dropbox/Chris Sansam/Synchd/Organized Dropbox/Projects and Data/Projects/Cowork_Projects/2026_Analyze_FACS_Data/facs_pseudocolor_workflow"` and
   `ph3_check_dir=$(mktemp -d)`; run
   `(cd "$ph3_check_dir" && R CMD build "$ph3_repo")`, then
   `R CMD check --no-manual "$ph3_check_dir"/facspseudocolor_*.tar.gz --output="$ph3_check_dir"`.
10. Repository inspection: `git diff --check`,
    `git status --short --untracked-files=all`, and
    `git status --short --ignored`.

Expected success: each focused/full test has zero failures/errors and only
explicitly asserted warnings; Python reports `OK`; package check reports
`Status: OK` with zero errors, warnings, or notes; and `git diff --check` has no
output. Commands 1–8 and 10 must generate no repository file. Command 9 may
generate only the package tarball and check directory beneath
`$ph3_check_dir`, never in the repository. Stop on any changed Slice 1–4,
legacy, EdU, or POI result; aggregation/schema/provenance mismatch; unexpected
warning; test/check failure or reviewable note; network attempt; experimental-
input access; or unexpected repository artifact. Report only command, `PASS`
or `FAIL`, failed names, a short relevant error excerpt, and unexpected files.

There is no unresolved scientific uncertainty after the confirmed supplement;
the sole implementation assumption is that the merged Slice 1–4 tables and
active explicit manifest remain authoritative at the immediate aggregation
boundary.

## First independent-review findings

The scientific-integrity review classified two findings as blocking. First,
cross-table agreement alone did not bind acquisition/sample/operation/manifest
identity back to the active analysis provenance. Second, grouping occurred
without first proving one stable `condition_index` per condition and one stable
`replicate_index` per condition × replicate. The artifact/security review
classified one finding as important: eligibility QC accepted any sixteen unique
dimension/reason pairs rather than the exact authoritative Slice 4 sequence.

The sole implementer is correcting all three findings with active provenance
binding, pre-grouping index validation, exact QC sequence validation, and
focused coordinated-mutation, split-index, and QC replacement/missing/extra
tests. No aggregation formula or approved schema is changing. Execution remains
**NOT RUN — local execution required**.

### First-review corrections

All three findings are resolved in the corrected diff. Aggregation now receives
the active analysis provenance and reconstructs the authoritative acquisition
map only from its retained export operations: acquisition ID, sample ID,
prefix, export-operation ID, and full-manifest digest. Every Slice 4 table must
agree per acquisition, and the metric identity must exactly match that active
map. The two active containment rows must also bind the same acquisition,
prefix, operation, and manifest digest with validated G1 and pH3-positive
statuses. Coordinated mutations across all four Slice 4 tables therefore cannot
manufacture a new internally consistent identity.

Before any grouping, the manifest must now provide exactly one
`condition_index` per condition and exactly one `replicate_index` per condition
× replicate. Split indices fail with a stable reason instead of creating extra
aggregate groups.

Eligibility QC validation now requires, in exact Slice 4 order, the imported-
row, identity-valid, eligibility, four exclusion, two partition, assigned-phase,
five configured-phase, and one eligible-Unassigned rows. The final reason is
limited to `none` or `configured_phase_gap`. The SYNTHETIC fixture now uses
those authoritative categories, and focused tests cover replacement, missing,
and extra QC rows.

Additional focused cases coordinate mutations of acquisition ID, sample ID,
export-operation ID, or manifest digest across every Slice 4 table while
leaving active provenance authoritative; each must fail. A mixed-table sample
mutation, split condition index, and split replicate index must also fail.

Files modified/generated remain exactly those in the final implementation
inventory above. The owner decision record remains unchanged. No experimental
input or ignored residue was accessed. No scientific method, output schema, or
aggregation formula changed during review resolution. Verification remains
**NOT RUN — local execution required**, and the ordered commands and stop
conditions above remain current. Both independent re-reviews are pending.

## Second independent-review findings

The re-review found two remaining fail-closed gaps. Index validation proved one
index per label but did not prove the reverse mapping, so distinct conditions
or replicates could share an ordering index and create ordering ties. The
repository manifest constructs replicate labels and indices globally, so both
condition and replicate label/index pairs require global bijections. In
addition, active containment was validated per expected acquisition without
rejecting stale rows for an unknown acquisition.

The sole implementer is adding global condition↔condition-index and
replicate↔replicate-index bijections before grouping, plus exact containment
acquisition membership and exact total/two-child row counts. Focused tests will
cover two conditions sharing an index, two replicate identities sharing an
index, one replicate label receiving inconsistent global indices, and a stale
unknown-acquisition containment row. No scientific formula or output schema is
changing. Execution remains **NOT RUN — local execution required**.

### Second-review corrections

Both findings are resolved. Before ordering or grouping, condition labels and
condition indices must now form a global bijection, and replicate labels and
replicate indices must form the repository manifest's global bijection. This
rejects label splits, reverse-index collisions, and all ordering ties that could
otherwise depend on manifest appearance. Focused SYNTHETIC cases cover two
conditions sharing one index, two replicate identities sharing one index, and
one replicate label assigned inconsistent indices.

Containment now must contain exactly the same acquisition-ID set as the active
export-manifest map and exactly two rows per active acquisition before any
per-acquisition validation. A stale unknown-acquisition row is fatal and has a
focused regression case. The existing per-acquisition validation continues to
require the exact G1/pH3-positive children and matching prefix, operation,
manifest digest, and validated status.

Files, scientific methods, schemas, and ordered human verification commands
remain unchanged. The owner decision record remains untouched. No experimental
input or ignored residue was accessed. Verification remains **NOT RUN — local
execution required**. Both independent re-reviews of this second corrected diff
are pending.

## Final artifact re-review finding

The artifact re-review found that global replicate↔replicate-index uniqueness
was stricter than the approved experimental unit. Replicates are identified
within condition, so a valid manifest may reuse the same replicate label and
index in another condition. The global condition↔condition-index bijection
remains required, but both directions of the replicate mapping must be scoped
to condition. The sole implementer is narrowing that validation and adding a
positive cross-condition reuse case while retaining within-condition split and
collision failures. No aggregation formula or output schema is changing, and
verification remains **NOT RUN — local execution required**.

The same review also found that active operation identifiers/digests and Slice
4 string provenance were checked for presence and agreement but not uniformly
for nonempty values. The correction will require non-`NA`, nonempty operation,
manifest, and other required string provenance values and will add coordinated
empty-value mutation cases.

### Final artifact-review corrections

Both findings are resolved. Conditions retain a global two-way label/index
mapping. Replicates now use the approved condition-scoped two-way mapping:
within a condition, one replicate label has one index and one index has one
replicate label. Reusing the same replicate label/index pair in another
condition is explicitly accepted and covered by a positive SYNTHETIC test;
within-condition label splits and reverse collisions remain fatal. Ordering is
therefore deterministic without imposing an unapproved cross-condition
replicate identity.

Active export-operation IDs and full-manifest digests must now be scalar,
non-`NA`, and nonempty. Active acquisition IDs, sample IDs, and prefixes receive
the same check. Every required string provenance field in all four Slice 4
tables must be non-`NA` and nonempty in addition to matching authoritative
provenance. Focused cases cover coordinated empty source operation/digest
fields and fully coordinated empty active/source/containment operation or
digest values; each fails closed.

The owner decision record, file inventory, aggregation formulas, schemas, and
ordered verification commands remain unchanged. No experimental input or
ignored residue was accessed. Verification remains **NOT RUN — local execution
required**. The corrected complete diff is ready for both final independent
re-reviews.

## Final independent-review outcomes

The final scientific-integrity review inspected the complete corrected diff
and relevant surrounding code and reported no blocking, important, or advisory
findings. It found the approved unweighted technical-acquisition and biological-
replicate estimands preserved, with no pseudoreplication, pooling, count
weighting, mapping alteration, silent source loss, zero/undefined error,
unsupported statistic, gate/threshold/normalization change, or legacy/EdU/POI
regression introduced.

The final artifact/security review inspected the same complete corrected diff
and reported no blocking, important, or advisory findings. It found no network
call, telemetry, upload, credential, unsafe execution, destructive I/O,
machine-specific production path, fallback, partial continuation, warning
suppression, generated artifact, or unrelated edit. Provenance binding,
condition-scoped replicate identity, deterministic ordering, exact QC rows,
and nonempty identifiers were clean.

Across review rounds, the sole implementer resolved active acquisition/sample/
operation/manifest and containment binding; stable condition and replicate
indices; the exact Slice 4 eligibility-QC sequence; reverse-index uniqueness;
rejection of stale containment acquisitions; condition-scoped rather than
global replicate identity; and nonempty active/source provenance. Each
correction added focused SYNTHETIC regression coverage without changing the
approved formulas or schemas.

Final non-ignored inventory relative to baseline:

- modified: `R/ph3.R`;
- modified: `tests/testthat/test-ph3-acquisition-metrics.R`;
- modified: `tests/testthat/test-ph3-input-containment.R`;
- untracked scientific-owner input preserved unchanged:
  `docs/implementation/ph3_slice5_aggregation_scientific_owner_decision_2026-08-29.md`;
- untracked implementation record:
  `docs/implementation/ph3_slice5_replicate_condition_aggregation_2026-08-29.md`;
- untracked focused test:
  `tests/testthat/test-ph3-replicate-condition-aggregation.R`.

No other tracked or non-ignored file is changed. Experimental inputs and
pre-existing ignored residue were not accessed or modified. No commit, push,
pull request, merge, or branch deletion occurred.

Final agent verification status is **NOT RUN — local execution required**. No
R, Python, project code, tests, packages, builds, checks, renders, analyses,
snapshots, or benchmarks were executed. The exact ordered human verification
commands, prerequisites, expected results, permitted external-temporary build
outputs, and stop conditions under “Ordered human verification commands” above
are the final handoff sequence.

## First user-run verification and test-only correction

The user reported focused Slice 5
`ph3-replicate-condition-aggregation`: **PASS**. The subsequent Slice 4
`ph3-acquisition-metrics` regression reported **FAIL 2, WARN 0, PASS 59**, at
the stale attachment assertions then located at
`tests/testthat/test-ph3-acquisition-metrics.R` lines 348 and 352. No later
command in the ordered sequence is recorded as run.

Static inspection found no production defect. The Slice 4 test still required
exactly the four pre-Slice-5 table names and rejected any biological table,
which directly contradicted the confirmed additive Slice 5 attachment. The
test-only correction retitles that case and now requires exactly the four
retained Slice 4 tables followed by the four owner-approved Slice 5 tables. It
compares every retained Slice 4 table exactly with independently derived Slice
4 output; requires all eight objects to be data frames; requires exact 5/5/16/4
Slice 4 and 5/5/5/5 Slice 5 row shapes; binds the aggregate tables to the exact
SYNTHETIC condition, replicate, acquisition provenance, aggregation levels,
and source-replicate identity; preserves the retained classification assertion;
and continues to reject the legacy `ph3` attachment plus CSV, RDS, plot,
report, geometry, GUI, compatibility/alias, and Slice 6+ outputs. No unrelated
assertion was weakened.

The corrected tree is **NOT RUN — local execution required**. The final
non-ignored inventory is the inventory above with
`tests/testthat/test-ph3-acquisition-metrics.R` now included as a modified
test-only correction. Human verification must restart with the revised ordered
sequence above: Slice 4 first, then Slice 5 and every remaining command in
order. The earlier Slice 5 pass applies only to the pre-correction tree and does
not verify this corrected tree. No production source, owner decision, formula,
schema, experimental input, or ignored residue changed during this correction.

### Renewed reviews after the test-only correction

The renewed scientific-integrity review inspected the complete corrected diff,
including the Slice 4 regression correction, and reported no blocking,
important, or advisory findings. It confirmed that the test now protects exact
Slice 4 retention and the approved Slice 5 attachment without changing or
weakening any scientific formula, mapping, eligibility, gate, threshold,
normalization, aggregation unit, statistic, or legacy/EdU/POI behavior.

The renewed artifact/security review inspected the same complete corrected
diff and reported no blocking, important, or advisory findings. It confirmed
that the correction is test/documentation-only, preserves strict object and
out-of-scope assertions, and adds no network behavior, telemetry, unsafe I/O,
fallback, generated artifact, machine-specific production path, or unrelated
edit.

Final corrected-tree non-ignored inventory relative to baseline:

- modified: `R/ph3.R`;
- modified: `tests/testthat/test-ph3-acquisition-metrics.R`;
- modified: `tests/testthat/test-ph3-input-containment.R`;
- untracked owner decision preserved unchanged:
  `docs/implementation/ph3_slice5_aggregation_scientific_owner_decision_2026-08-29.md`;
- untracked implementation record:
  `docs/implementation/ph3_slice5_replicate_condition_aggregation_2026-08-29.md`;
- untracked focused Slice 5 test:
  `tests/testthat/test-ph3-replicate-condition-aggregation.R`.

No other tracked or non-ignored file is changed. No commit, push, pull request,
merge, or branch deletion occurred. Experimental inputs and ignored residue
were not accessed or modified. The corrected tree remains **NOT RUN — local
execution required**. At that correction stage, the human sequence began with Slice 4
`ph3-acquisition-metrics`, then Slice 5
`ph3-replicate-condition-aggregation`, followed by the remaining commands under
“Ordered human verification commands” in exact order.

## Second user-run verification and fixture correction

The user reported the corrected-tree results:

- Slice 4 `ph3-acquisition-metrics`: **PASS**;
- Slice 5 `ph3-replicate-condition-aggregation`: **PASS**;
- Slice 3 `ph3-event-classification`: **PASS**;
- Slice 2 `ph3-input-containment`: **FAIL 1, WARN 0, PASS 516**.

The exact failure was in
`tests/testthat/test-ph3-input-containment.R` at the high-level test beginning
near line 841, with stable reason `mixed_or_missing_provenance`: the Slice 4
`positivity_method_id` was not the approved exact identifier. The reported
backtrace entered replicate/condition aggregation validation through
`R/ph3.R` near the then-current lines 1476, 1250, and 1213.

Static inspection of the complete fixture and canonical approval/provenance
path found no production defect. The Slice 1 manifest validator intentionally
requires complete nonblank owner approval metadata but does not reinterpret an
arbitrary identifier as scientifically approved. Slice 3 retains the active
manifest identifier, Slice 4 binds to that same manifest, and Slice 5 correctly
requires the confirmed production identifier
`flowjo_owner_approved_positive_population_v1`. The high-level SYNTHETIC
operation fixture alone still declared the stale placeholder
`SYNTHETIC_METHOD`.

The smallest correction changes only that fixture's positivity method ID to
the confirmed canonical value. Positivity membership, semantics, mappings,
gates, thresholds, classification, acquisition formulas, aggregation formulas,
and production validation are unchanged. Arbitrary, missing, mixed, mutated,
or non-approved positivity IDs remain fatal; no validator was relaxed.

The corrected tree is again **NOT RUN — local execution required**. The final
non-ignored inventory remains the six paths listed above; the existing modified
`tests/testthat/test-ph3-input-containment.R` now also contains this one-line
canonical fixture correction. Verification must restart with Slice 2
`ph3-input-containment`, then run Slice 4, Slice 5, Slice 3, and every remaining
command in the revised ordered sequence above. No command was executed by an
agent, and no commit, production-source correction, owner-record edit,
experimental-input access, or ignored-residue access occurred.

### Renewed reviews after the second human correction

The renewed scientific-integrity review inspected the complete corrected diff
and reported no blocking, important, or advisory findings. It confirmed that
the fixture now names the exact approved positivity method while production
continues to reject arbitrary, mixed, missing, mutated, or non-approved method
IDs. Positivity membership and semantics, mappings, gates, thresholds,
classification, acquisition metrics, aggregation formulas, experimental units,
statistics, and legacy/EdU/POI behavior remain unchanged.

The renewed artifact/security review inspected the same complete corrected
diff and reported no blocking, important, or advisory findings. It confirmed
that the correction is confined to the canonical SYNTHETIC fixture and this
record, with no network behavior, telemetry, unsafe or destructive I/O,
fallback, generated artifact, production-path change, or unrelated edit.

Final corrected-tree non-ignored inventory relative to baseline:

- modified: `R/ph3.R`;
- modified: `tests/testthat/test-ph3-acquisition-metrics.R`;
- modified: `tests/testthat/test-ph3-input-containment.R`;
- untracked owner decision preserved unchanged:
  `docs/implementation/ph3_slice5_aggregation_scientific_owner_decision_2026-08-29.md`;
- untracked implementation record:
  `docs/implementation/ph3_slice5_replicate_condition_aggregation_2026-08-29.md`;
- untracked focused Slice 5 test:
  `tests/testthat/test-ph3-replicate-condition-aggregation.R`.

No other tracked or non-ignored file is changed. At the time of that renewed
review, the corrected tree was **NOT RUN — local execution required** and the
requested rerun order began with Slice 2 containment. That historical status
was subsequently superseded by the final user-run verification below. No
commit, push, pull request, merge, branch deletion, experimental-input access,
or ignored-residue access occurred.

## Final user-run verification

The user reported the final corrected-tree results exactly as follows:

- Slice 2 `ph3-input-containment`: **PASS**;
- Slice 4 `ph3-acquisition-metrics`: **PASS**;
- Slice 5 `ph3-replicate-condition-aggregation`: **PASS**;
- Slice 3 `ph3-event-classification`: **PASS**;
- legacy `ph3`: **PASS**;
- `python-boundary`: **PASS**;
- `tests/python/test_export_contract.py`: **OK**;
- complete R suite: **PASS**;
- external-temporary `R CMD build` and `R CMD check --no-manual`: user
  reported **OK / Status: OK**;
- `git diff --check`: no output.

The user supplied the exact final non-ignored inventory:

- modified: `R/ph3.R`;
- modified: `tests/testthat/test-ph3-acquisition-metrics.R`;
- modified: `tests/testthat/test-ph3-input-containment.R`;
- untracked owner decision preserved unchanged:
  `docs/implementation/ph3_slice5_aggregation_scientific_owner_decision_2026-08-29.md`;
- untracked implementation record:
  `docs/implementation/ph3_slice5_replicate_condition_aggregation_2026-08-29.md`;
- untracked focused Slice 5 test:
  `tests/testthat/test-ph3-replicate-condition-aggregation.R`.

The user reported that the ignored inventory was unchanged from the inherited
pre-existing inventory. Build/check outputs were confined to the external
temporary location. No unexpected non-ignored repository output was reported.

This final verification was performed by the user, not by an agent. Agents ran
no R, Python, project code, tests, package operations, builds, checks, renders,
analyses, snapshots, or benchmarks. The successful user-run results supersede
the historical pending/`NOT RUN` state for the current corrected tree while
preserving the earlier failure and correction chronology above. The final
assumption is that the user-reported concise results and exact inventories
faithfully describe the verified local tree. The implementation remains local
and uncommitted; no push, pull request, merge, or branch deletion occurred.
