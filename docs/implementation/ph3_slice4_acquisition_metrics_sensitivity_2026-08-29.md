# pH3 Slice 4 acquisition metrics and boundary-sensitivity QC

**Status:** IMPLEMENTED, INDEPENDENTLY REVIEWED, AND USER-VERIFIED

## Scope and provenance

Canonical source:
`/Users/sansamc/OMRF Dropbox/Chris Sansam/Synchd/Organized Dropbox/Projects and Data/Projects/Cowork_Projects/2026_Analyze_FACS_Data/facs_pseudocolor_workflow`

Branch: `feature/ph3-acquisition-metrics-sensitivity`. Exact starting commit:
`ac6cdf71fb9999128fa82c92b41699113dc5689a`, with a clean working tree.

This slice consumes only the retained authoritative Slice 3 table at
`analysis$normalized_data[[acquisition]]$ph3_event_classification`. It does not
read or independently classify raw or normalized event tables, derive pH3
positivity, or fall back to legacy calculations. Production fails closed when
the classification is absent, empty, duplicated, malformed, provenance-
inconsistent, containment-unvalidated, or unreconciled. Before calculation,
the retained classification's configuration digest, analysis ID, export
operation, full-manifest digest, positivity method, acquisition ID, sample ID,
and prefix must match the active validated configuration and immutable export
manifest references retained in analysis provenance.

No experimental `.fcs`, `.wsp`, `.prism`, CSV, RDS, figure, report, archive, or
historical analysis artifact was opened, hashed, copied, modified, or
generated. No external data transmission, telemetry, upload, or network call
was added.

## Implemented tables and formulas

Production analysis attaches exactly these acquisition-level tables beneath
`analysis$quantitation`:

- `ph3_metrics_acquisition`: five normative long-form rows per acquisition.
  `ph3_2to4n_positivity_percent = 100*|P∩A|/|A|`;
  `ph3_within_4n_positivity_percent = 100*|P∩F|/|F|`;
  `ph3_4n_positive_prevalence_within_2to4n_percent = 100*|P∩F|/|A|`;
  `ph3_within_sub_4n_positivity_percent = 100*|P∩S|/|S|`; and
  `ph3_sub_4n_positive_prevalence_within_2to4n_percent = 100*|P∩S|/|A|`.
- `ph3_phase_prevalence`: one row each for G1, Early S, Mid S, Late S, and
  G2/M. Every numerator is eligible positive membership in the retained
  configured phase and every denominator is `|A|`.
- `ph3_event_eligibility_qc`: long acquisition-level count rows for imported
  classification rows, eligible 2to4N, identity-invalid, nonfinite DNA,
  below-`b0`, above-`b5`, sub-4N, 4N, each configured phase, and eligible
  configured-phase Unassigned by the retained Slice 3 reason.
- `ph3_4n_boundary_sensitivity_qc`: exactly `primary`, `lower_outward`,
  `lower_inward`, and `upper_inward`, corresponding to closed intervals
  `[b4,b5]`, `[b4-delta,b5]`, `[b4+delta,b5]`, and `[b4,b5-delta]`, where
  `delta = ph3_boundary_sensitivity_fraction * dna_2n_value`. Every predicate
  is additionally restricted to retained eligible `A`, so lower-outward never
  admits an ineligible event. No upper expansion exists.

For all percentages, a positive denominator with a zero numerator gives numeric
zero and `ok`. A zero denominator gives `NA_real_` and
`undefined_zero_denominator`. Counts are stored as integers. Classification
validation rechecks the exact Slice 3 predicates using only fields retained in
that authoritative table: identity validity, DNA finiteness, inclusive
`[b0,b5]` eligibility, `[b0,b4)`/`[b4,b5]` ownership, configured phase
interval ownership, assignment, and the approved exclusion/Unassigned reason
sets. Metrics still consume retained memberships and are never reconstructed
from raw or separately normalized input tables. The validator additionally
requires `dna_finite` to equal `is.finite(dna_norm)` exactly, preserves Slice
3's production invariant that every `identity_valid` value is `TRUE`, and binds
the retained `exact_direct_identity_multiset_containment` method to both active
G1 and pH3-positive validated containment records. Invalid delta,
nonfinite/inverted intervals, a nonfinite percentage outside the approved
zero-denominator case, or failure of the primary-sensitivity invariance check
is fatal.

Every table retains available acquisition, sample/prefix, condition,
biological-replicate, and technical-acquisition keys plus acquisition level,
classification/output schema, positivity, eligibility, interval, 4N, sub-4N,
containment, configuration, export-operation, and manifest provenance. The
established method IDs are:

- `flowjo_owner_approved_positive_population_v1` (the exact manifest-approved
  value is retained in production rows);
- `identity_validated_finite_dna_b0_b5_inclusive_v1`;
- `configured_shared_boundaries_left_closed_v1`;
- `configured_b4_b5_closed_v1`;
- `configured_b0_b4_left_closed_right_open_v1`;
- the exact retained Slice 2 `containment_method_id`;
- `ph3-event-classification-1.0.0` and future output schema `ph3-1.0.0`.

No aggregation, pooled or weighted result, biological-replicate table,
condition statistic, CSV, compatibility alias, plot, report, geometry object,
or later-slice output is attached. Explicit legacy
`legacy_count_only_unverified_v1` behavior and warnings remain on the unchanged
legacy `quantify_ph3()` path.

## Files

Read: the complete Slice 4 handoff; workspace `AGENTS.md`; the complete
scientific-owner approval and implementation plan/test specification; complete
Slice 1, Slice 2, and Slice 3 implementation records; `R/ph3.R`, `R/analysis.R`,
focused manifest/configuration helpers in `R/input.R`, `R/config.R`, and
`R/pseudocolor_helpers.R`; and the focused existing Slice 2/3 and legacy pH3
tests. No experimental project directory or input was inspected.

Modified: `R/ph3.R`, `R/analysis.R`, and
`tests/testthat/test-ph3-input-containment.R`.

Generated as source/test documentation only:
`tests/testthat/test-ph3-acquisition-metrics.R` and this implementation record.
No analysis output, report, image, cache, binary, archive, package build, or
rendered artifact was generated.

Sample mappings, exclusions, gates, positivity thresholds, G1 normalization,
statistical methods, aggregation methods, and biological claims changed:
none. The approved Slice 4 acquisition estimands and QC are newly implemented;
the retained Slice 3 membership, boundaries, eligibility, and reason values are
not mutated.

## SYNTHETIC tests prepared

The focused in-memory test covers all five hand-calculated metrics; distinct
within-region and prevalence denominators; exact `b0`, `b4`, and `b5`
ownership; positive-denominator zero results; zero `A`, `S`, and `F`
denominators; positives outside eligibility; all five configured-phase
prevalences; classification-reason QC and Unassigned reconciliation; exactly
four closed sensitivity rows and no upper expansion; exact-boundary predicates;
primary invariance; invalid/inverted sensitivity; missing, malformed, and
duplicated classification; unchanged retained classification; and production
attachment of only the four approved acquisition tables. Review-driven cases
also cover active configuration/manifest/analysis and acquisition/sample/prefix
binding, each exact Slice 3 predicate, approved reason sets, and explicit
`identity_valid` and `configured_phase_assigned` QC rows/counts. Existing
regression commands below protect Slice 3, Slice 2, legacy pH3, EdU, and POI
behavior.

## Independent scientific-review resolution

The first scientific-integrity review identified two blocking gaps—retained
classification provenance was not rebound to the active analysis, and predicate
validation did not recheck every exact Slice 3 semantic—plus one important QC
traceability gap: explicit identity-valid and configured-phase-assigned rows
were absent. The implementer resolved those findings with the active
configuration/analysis/manifest/acquisition binding, exact retained-field
predicate checks described above, two explicit QC rows, and repeated exact
identity/assignment counts. No metric formula, gate, threshold, normalization,
sample mapping, or biological interpretation changed. The corrected diff is
awaiting renewed independent read-only review and remains unexecuted.

The renewed scientific review then identified coordinated-tampering cases in
which internally consistent altered `dna_finite` or `identity_valid` fields
could evade the first correction, plus an unbound containment method. The
second correction adds the exact finiteness/all-valid invariants and active
two-child containment binding described above. Focused negative tests alter
the dependent eligibility, partition, assignment, and reason fields together
and still require fail-closed rejection; a mismatched active containment method
is also rejected. This second correction likewise changes no metric or
scientific method and remains unexecuted pending re-review and human testing.

The final important review finding was that the two active containment rows'
method/status were bound but their operation, manifest, and prefix linkage was
not rechecked at the Slice 4 boundary. The validator now requires those fields
to exist and requires both child rows' `prefix`, `export_operation_id`, and
`manifest_digest` to equal the retained classification and active manifest
binding. Focused SYNTHETIC tests independently mismatch each value and remove a
required field; every case must fail closed with `active_provenance_mismatch`.
No metric or scientific method changed, and the correction remains unexecuted
pending final re-review and human verification.

## First user-run verification and test-only correction

The user reported the following local results for the reviewed pre-correction
tree:

- focused `ph3-acquisition-metrics`: **PASS**;
- Slice 3 `ph3-event-classification`: **PASS**;
- Slice 2 `ph3-input-containment`: **FAIL 2, WARN 0, PASS 509**.

Both Slice 2 failures were the stale empty-quantitation assertions at the then
current test lines 861 and 867. Production now correctly attaches the four
approved Slice 4 acquisition tables, so the old expectation that quantitation
remain an empty list no longer described the approved contract. Static review
found no production defect.

The test-only correction retitles that high-level case and requires the exact
four approved quantitation names, data-frame types, expected five-metric,
five-phase, four-sensitivity, and sixteen-QC-row per-acquisition shapes across
the exact `SYNTHETIC-reference` and `SYNTHETIC-treatment` acquisition IDs,
acquisition-level markers, and exact total table sizes.
It explicitly rejects legacy `ph3` and all aggregation/condition-statistic/CSV/
plot/report/geometry objects. All prior containment, retained classification,
event-identity, transient verified-event stripping, and rejected manual
`quantify_ph3()` assertions remain. Quantitation is captured before the rejected
manual call and must remain exactly identical afterward.

Agents did not execute the corrected test or any project code. At that point,
the corrected tree remained **NOT RUN — local execution required**; the focused
Slice 4 and Slice 3 passes above were valid user-reported evidence only for the
pre-correction tree, and the Slice 2 result remained the reported failure until
rerun.

## Final user-run verification

The user subsequently ran the required verification locally on the corrected
tree and reported all checks successful:

- corrected Slice 2 `ph3-input-containment`: **PASS**;
- existing legacy pH3 regression: **PASS**;
- R-to-Python boundary regression: **PASS**;
- Python export-contract unit tests: **OK**;
- complete R test suite, including EdU and POI regressions: **PASS**;
- external-temporary-directory package build and `R CMD check --no-manual`:
  **Status: OK**;
- `git diff --check`: no output;
- repository inventory: exactly the five intended Slice 4 paths plus the
  unchanged inherited ignored residue.

The focused Slice 4 acquisition-metrics and Slice 3 event-classification tests
had already passed on the pre-correction tree. The correction changed only the
stale Slice 2 test expectations and this implementation record; no production
calculation changed after those focused passes. Verification was performed by
the user, not by an agent.

## Local verification contract

Agents did not execute R, Python, project code, tests, packages, builds, checks,
renders, examples, analyses, snapshots, or benchmarks. The commands below were
executed locally by the user with the successful results recorded above.

Working directory for every repository command:
`/Users/sansamc/OMRF Dropbox/Chris Sansam/Synchd/Organized Dropbox/Projects and Data/Projects/Cowork_Projects/2026_Analyze_FACS_Data/facs_pseudocolor_workflow`

Prerequisites: R >= 4.2 with the package's declared dependencies plus
`testthat`, `withr`, and `devtools`; supported Python with repository contract
test dependencies for the existing exporter regression. Dependencies must
already be installed. No experimental input is needed or permitted.

1. Slice 4 focused tests:
   `Rscript -e 'devtools::test(filter = "ph3-acquisition-metrics", stop_on_failure = TRUE)'`
2. Slice 3 regression:
   `Rscript -e 'devtools::test(filter = "ph3-event-classification", stop_on_failure = TRUE)'`
3. Slice 2 regression:
   `Rscript -e 'devtools::test(filter = "ph3-input-containment", stop_on_failure = TRUE)'`
4. Existing legacy pH3 regression:
   `Rscript -e 'devtools::test(filter = "^ph3$", stop_on_failure = TRUE)'`
5. Export/R boundary:
   `Rscript -e 'devtools::test(filter = "python-boundary", stop_on_failure = TRUE)'`
6. Repository-only Python contract:
   `python3 -m unittest -v tests/python/test_export_contract.py`
7. Full R suite, including EdU and POI regressions:
   `Rscript -e 'devtools::test(stop_on_failure = TRUE)'`
8. Set
   `ph3_repo="/Users/sansamc/OMRF Dropbox/Chris Sansam/Synchd/Organized Dropbox/Projects and Data/Projects/Cowork_Projects/2026_Analyze_FACS_Data/facs_pseudocolor_workflow"`
   and `ph3_check_dir=$(mktemp -d)`, then run
   `(cd "$ph3_check_dir" && R CMD build "$ph3_repo")`, followed by
   `R CMD check --no-manual "$ph3_check_dir"/facspseudocolor_*.tar.gz --output="$ph3_check_dir"`.
9. Run `git diff --check`, `git status --short --untracked-files=all`, and
   `git status --short --ignored` from the canonical repository.

Expected success: all focused and full tests have zero failures/errors and only
explicitly asserted legacy warnings; Python reports `OK`; package check reports
`Status: OK` with zero errors, warnings, or notes; build/check outputs exist
only below the temporary directory; and repository status contains only the
five exact intended Slice 4 paths—`R/analysis.R`, `R/ph3.R`,
`tests/testthat/test-ph3-acquisition-metrics.R`,
`tests/testthat/test-ph3-input-containment.R`, and this implementation
record—plus the inherited ignored-residue inventory.

Stop on any changed identity, eligibility, positivity, phase, interval,
containment, provenance, legacy, EdU, or POI value; any denominator mismatch;
any fifth or upper-expansion production sensitivity row; any nonzero undefined
percentage; any network attempt; any warning beyond an asserted legacy warning;
any test/check failure or reviewable note; any experimental-input selection;
or any unexpected repository artifact.

Return only the command, `PASS` or `FAIL`, failed test/check names, a short
relevant error excerpt, and unexpected generated files. Do not return large
successful logs.

## Assumptions and unresolved uncertainty

The merged Slice 1–3 direct identity, containment, manifest linkage,
classification schema, configured phase labels, and method IDs remain
authoritative. Manifest mappings are explicit and are not inferred from file
or row order. The current implementation has been locally verified by the user,
and the final independent scientific-integrity and artifact/security reviews
were clean with no blocking, important, or advisory findings.
