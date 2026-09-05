# pH3 Slice 6 output-contract audit and decision memo

**Status:** AUDIT COMPLETE; SCIENTIFIC-OWNER DECISIONS REQUIRED BEFORE IMPLEMENTATION  
**Date:** 2026-08-29

## Scope, source, and audit state

Canonical source:
`/Users/sansamc/OMRF Dropbox/Chris Sansam/Synchd/Organized Dropbox/Projects and Data/Projects/Cowork_Projects/2026_Analyze_FACS_Data/facs_pseudocolor_workflow`

Audit branch: `feature/ph3-output-contract-audit`  
Clean synchronized starting branch: `main`  
Baseline commit: `ead9e13e18c8f6243e9622161e0c9b475a42da0b`

This is documentation and planning only. It does not authorize or implement
analysis, export, plotting, report, Quarto, geometry, GUI, configuration, or
compatibility behavior. The authoritative numerical inputs remain the four
merged Slice 4 acquisition tables and four merged Slice 5 aggregate tables.
Reporting must consume them unchanged and must not reclassify events or
recompute metrics.

Audit checkpoints:

1. Source-of-truth review: complete. The confirmed pH3 owner record, plan,
   final Slice 1–5 records, current pH3 result/report surfaces and focused
   tests, report documentation, and merged EdU output contract were reviewed.
2. Candidate-output matrix: complete; see below.
3. Owner-decision extraction: complete; six bounded decisions remain.
4. Independent review resolution: complete. Both final re-reviews found no
   blocking, important, or advisory findings after the resolutions recorded
   below.

No experimental input or generated report was inspected. Historical reports,
whole-dataset audits, historical QMDs, and EdU method records were not used to
infer pH3 science. The EdU implementation record was used only to compare the
existing explicit-write, overwrite-refusal, and report-input architecture.

Exact bounded file inventory read by the orchestrator or independent reviewers:

- workspace records: `AGENTS.md`,
  `HANDOFF_PH3_SLICE6_OUTPUT_CONTRACT_AUDIT.md`,
  `ph3_scientific_owner_approval_2026-08-24.md`, and
  `ph3_implementation_plan_and_test_spec_5de77a9.md`;
- merged implementation records:
  `docs/implementation/ph3_slice1_export_identity_provenance_2026-08-25.md`,
  `docs/implementation/ph3_slice2_input_validation_containment_2026-08-28.md`,
  `docs/implementation/ph3_slice3_event_classification_2026-08-29.md`,
  `docs/implementation/ph3_slice4_acquisition_metrics_sensitivity_2026-08-29.md`,
  `docs/implementation/ph3_slice5_aggregation_scientific_owner_decision_2026-08-29.md`,
  and
  `docs/implementation/ph3_slice5_replicate_condition_aggregation_2026-08-29.md`;
- focused current APIs/tests: `R/ph3.R`, `R/analysis.R`, `R/results.R`,
  `R/artifacts.R`, `tests/testthat/test-ph3-acquisition-metrics.R`,
  `tests/testthat/test-ph3-replicate-condition-aggregation.R`, and
  `tests/testthat/test-results-api.R`;
- report/template sources and documentation: `inst/quarto/facs_ph3_4n.qmd`,
  `inst/quarto/_report-setup.R`, `docs/PH3_MODE.md`, `docs/REPORTS.md`, and
  `docs/CONFIGURATION.md`;
- architectural comparison only:
  `docs/implementation/edu_output_contract_implementation_2026-08-23.md`.

## Merged table inventory and intended roles

| Table | Level | Fixed content | Proposed Slice 6 role |
|---|---|---|---|
| `ph3_metrics_acquisition` | technical acquisition | Five approved metrics with counts, denominators, value/status, bounds, keys, and provenance | Machine-readable; required report acquisition/QC detail, with placement/prominence decided by S6-D01 |
| `ph3_phase_prevalence` | technical acquisition | Five configured-phase prevalences, each over eligible 2to4N | Machine-readable; secondary acquisition detail, never within-phase positivity |
| `ph3_event_eligibility_qc` | technical acquisition | Imported, eligible, exclusion, partition, phase, and Unassigned reason counts | Machine-readable and required standard QC table/display |
| `ph3_4n_boundary_sensitivity_qc` | technical acquisition | Primary plus three approved in-span perturbations with both denominators | Machine-readable and required standard QC table; plot only after Decision S6-D04 |
| `ph3_metrics_biological_replicate` | biological replicate | Unweighted means of finite technical percentages plus source counts/status | Machine-readable; required standard biological-replicate detail and proposed plot input |
| `ph3_phase_prevalence_biological_replicate` | biological replicate | Unweighted phase-prevalence means plus source counts/status | Machine-readable; secondary biological-replicate detail and possible plot input |
| `ph3_metrics_condition_summary` | condition | Descriptive mean, SD, SEM, biological-replicate counts/status | Machine-readable; required standard condition summary; possible plot input only under Decision S6-D02 |
| `ph3_phase_prevalence_condition_summary` | condition | Descriptive phase-prevalence mean, SD, SEM, biological-replicate counts/status | Machine-readable; secondary condition summary; possible plot input only under Decision S6-D03 |

All eight tables remain attached under `analysis$quantitation` in their merged
order. The two Slice 4 QC tables remain acquisition-level and are not
numerically aggregated. Acquisition counts must never be presented as
biological-replicate counts.

## Candidate standard-output matrix

| Candidate output | Source | Level | Purpose | Classification and rationale |
|---|---|---|---|---|
| Human-readable run summary | analysis metadata, manifest, warnings, all table statuses | run | Identify analysis/profile/schema, conditions, declared acquisitions/replicates, methods, provenance linkage, warnings, and undefined/partial states | Deterministic and required; contains no biological interpretation |
| Five-metric condition table | `ph3_metrics_condition_summary` | condition | Primary descriptive condition results | Approved and deterministic |
| Five-metric biological-replicate table | `ph3_metrics_biological_replicate` | biological replicate | Show independent points and technical-source completeness | Approved and deterministic |
| Five-metric acquisition detail | `ph3_metrics_acquisition` | acquisition | Audit technical values, numerators, denominators, and provenance | Approved; report prominence is lab-standard Decision S6-D01 |
| Phase-prevalence tables | three phase tables | acquisition/biorep/condition | Secondary configured-phase prevalence over eligible 2to4N | Approved numerically; report prominence/plot is Decision S6-D03 |
| Eligibility/reconciliation QC | `ph3_event_eligibility_qc` | acquisition | Exclusions, partition reconciliation, and Unassigned warning | Approved and required standard QC |
| 4N sensitivity table | `ph3_4n_boundary_sensitivity_qc` | acquisition | Show robustness without replacing primary interval | Approved and required standard QC |
| Primary metric plot | biological-replicate and/or condition metric table | biorep/condition | Visualize five approved estimands | Scientific-owner Decision S6-D02; plot/statistic/error-bar convention is not approved |
| Phase-prevalence plot | biological-replicate and/or condition phase table | biorep/condition | Visualize secondary prevalence | Scientific-owner Decision S6-D03 |
| Sensitivity plot | sensitivity QC table | acquisition | Triggered or advanced diagnostic | Approved only as advanced/QC-triggered; Decision S6-D04 must define trigger or advanced-only behavior |
| Unassigned diagnostic plot | eligibility QC plus existing authoritative classifications | acquisition | Diagnose nonzero eligible Unassigned | QC-triggered only; any event-level implementation is deferred unless separately bounded and approved |
| Verified FlowJo geometry overlay | linked geometry plus analysis provenance | acquisition | Nonquantitative gate display | Optional advanced output; must be omitted/unavailable when exact linkage fails; never required for numerical report |
| Derived 4N display region | none approved | display | Proposed display-only rectangle | Advanced method development; out of scope |
| Inferential statistics, p-values, confidence intervals, reference normalization | none | condition | Inference/comparison | Explicitly out of scope; requires new owner approval |
| Legacy all-row pH3 outputs/current `facs_ph3_4n.qmd` semantics | legacy `quantitation$ph3` | legacy | Historical compatibility | Explicitly excluded from the standard production profile; do not redirect to canonical meanings |

## Required standard semantics

Every standard metric label must state its biological question and denominator:

- overall positivity: pH3-positive among eligible 2to4N Single Cells, `|P∩A|/|A|`;
- within-4N positivity: pH3-positive among eligible 4N cells, `|P∩F|/|F|`;
- 4N-positive prevalence: pH3-positive 4N cells among eligible 2to4N cells, `|P∩F|/|A|`;
- within-sub-4N positivity: pH3-positive among eligible sub-4N cells, `|P∩S|/|S|`;
- sub-4N-positive prevalence: pH3-positive sub-4N cells among eligible 2to4N cells, `|P∩S|/|A|`;
- configured-phase prevalence: pH3-positive cells in the named configured phase
  among eligible 2to4N cells, never positivity within that phase.

The report must display exact interval bounds and inclusivity, units as percent,
condition and biological-replicate identities in explicit manifest order, and
technical acquisitions only as technical-source detail. It must not call a
technical acquisition a replicate. A condition with one declared biological
replicate may show its descriptive value and the independent point, must show
SD/SEM as undefined, and must say that replicated inference is unavailable.
No reference condition, reference ratio, comparison direction, or normalized
value may be inferred. If a future approved profile merely annotates a named
reference condition, the annotation must remain descriptive and must not alter
values.

`NA` must render as `Not defined`, accompanied by its exact structured status
and source-unit counts. Numeric zero must remain `0`. `ok_partial_undefined`
must visibly identify omitted undefined source units; it must not look like
`ok`. Failed provenance, missing required rows, mixed schema/method/config,
failed containment, or unreconciled QC must stop before report or export output
is written. Nonzero eligible Unassigned must remain a prominent warning.

Required run-summary provenance includes analysis ID; `ph3-1.0.0`; report
profile/template version; source commit/workflow version; configuration digest;
manifest digest and input-manifest keys; export-operation IDs; workspace and
gate approval/linkage fields; condition/replicate/acquisition mapping counts;
all positivity, containment, eligibility, interval, 4N, sub-4N, classification,
aggregation, and schema method IDs; exact bounds/inclusivity; result-status
counts; exclusions, Unassigned, sensitivity, warnings, and legacy-profile
status. Long or nested provenance should remain losslessly linked to the
analysis RDS or a local structured manifest rather than be flattened
ambiguously into CSV.

## Scientific-owner decisions required

### S6-D01 — standard report table depth

- **A:** Show condition and biological-replicate primary tables in the main
  results, all acquisition/QC tables in standard QC/detail sections, and all
  eight tables in machine-readable export. Consequence: concise main narrative
  with complete traceability.
- **B:** Show all eight tables in the standard report body. Consequence: maximum
  visibility but a substantially longer routine report.

No option may omit the eligibility or sensitivity QC tables from the standard
report, or omit acquisition detail from the durable output contract.

### S6-D02 — primary five-metric plot contract

- **A:** Biological-replicate points grouped by condition, with condition mean
  shown descriptively; no error bars until separately chosen. Consequence:
  independent units are visible and one-replicate conditions are honest.
- **B:** Defer the plot and make the fixed tables the only standard primary
  output. Consequence: no unapproved graphic statistic; later plot approval is
  required.

If a plot is approved, the owner must also confirm layout (five facets versus
separate panels), point positioning, condition ordering from the manifest, and
whether SD or SEM is shown when at least two finite biological replicates
exist. Neither SD nor SEM may be selected by implementation convenience.

### S6-D03 — secondary phase-prevalence display

- **A:** Standard table only; defer a plot. Consequence: the approved secondary
  values remain visible without choosing a composition/line/bar convention.
- **B:** Add a biological-replicate-point display faceted by configured phase,
  subject to explicit confirmation of any mean/error-bar overlay. Consequence:
  easier phase comparison but greater risk of confusing prevalence with
  within-phase positivity unless labels are exact.

Stacked composition is not an available choice because these values are
prevalences over eligible 2to4N and are not approved as a compositional
estimand.

### S6-D04 — sensitivity and Unassigned diagnostics

- **A:** Keep the required sensitivity and eligibility tables standard, with
  plots available only through an explicit advanced option. Consequence: no
  unapproved automatic threshold.
- **B:** Define owner-approved deterministic QC triggers for adding the relevant
  diagnostic plots. Consequence: requires a separate exact trigger definition;
  none currently exists.

### S6-D05 — machine-readable export set

- **A:** Export all eight merged tables as same-named CSVs plus the complete
  analysis RDS; record their common provenance linkage in the run summary.
  Consequence: smallest lossless contract using current authoritative tables.
- **B:** Export the six historically named canonical CSVs only, requiring an
  explicit mapping/renaming decision for the two phase aggregate tables and a
  decision about the two additional Slice 5 condition/replicate tables.
  Consequence: preserves the older six-file list but cannot be implemented
  without resolving name/content ambiguity.

CSV export must be explicit, local, preflight all requested paths before any
write, refuse existing files unless `overwrite = TRUE`, write no index column,
preserve stable types/NA/status, and produce no partial set. RDS is the lossless
source of truth. No automatic destination, timestamped fallback, network path,
or silent overwrite is allowed.

### S6-D06 — standard report profile and geometry boundary

- **A:** Add a canonical pH3 standard profile that requires exactly one of
  explicit config or analysis RDS, consumes the eight tables, and treats a
  verified geometry overlay as a separate optional advanced section.
- **B:** Require verified geometry for the standard report. Consequence:
  numerical reporting would stop when geometry is unavailable even though the
  approved numerical analysis may proceed without it.

The existing `facs_ph3_4n.qmd` is legacy/incompatible with the canonical
contract: it uses the legacy table, all-Single-Cell denominator language,
unqualified G2/M/4N wording, mandatory geometry, and suppressed warnings. It
must not be silently repurposed. Quarto may select a profile and call package
APIs, but must contain no metric, aggregation, mapping, QC, or fallback logic.

## Deterministic implementation sequence after approval

1. Freeze S6-D01–D06 in a dated scientific-owner record, including every plot
   statistic/error-bar choice and exact export names.
2. Add a pH3 output validator that requires all eight tables, exact schemas,
   keys/order/provenance agreement, statuses, and Slice 4/5 reconciliation.
3. Add a pure report-model/bundle builder that consumes validated tables only;
   create deterministic labels, run summary, tables, warnings, and approved
   plots without event reclassification.
4. Extend the explicit result writer with the approved pH3 CSV/RDS set,
   all-path preflight, atomic/no-partial behavior, and overwrite refusal.
5. Add the new canonical Quarto profile as a thin presentation layer with
   exactly-one-input validation and explicit output destinations.
6. Update report/configuration/migration documentation, preserving legacy pH3,
   EdU, and POI APIs and semantics unchanged.

## Focused synthetic/input-free test specification

All future fixtures must be in-memory or unmistakably `SYNTHETIC` test-only
fixtures. No experimental or example-data analysis is permitted.

- Exact eight-table presence, schemas, types, unique keys, deterministic row
  order, `ph3-1.0.0`, and common provenance validation; every missing,
  duplicated, extra, mixed, or malformed case fails before output.
- Report model consumes table values exactly and never calls event
  classification/metric/aggregation functions; mutation of each authoritative
  table value appears unchanged in its corresponding display.
- Every metric/phase label contains the exact question, denominator, units,
  interval, and inclusivity; forbidden ambiguous labels and stacked phase
  composition are absent.
- Acquisition, biological-replicate, and condition levels remain distinct;
  technical acquisitions never become independent points or sample size.
- Zero, zero denominator, `ok_partial_undefined`,
  `undefined_no_finite_values`, and one-biological-replicate cases render with
  exact values/statuses and no implied replicated inference.
- Nonzero Unassigned warning is prominent; required QC tables are present;
  sensitivity never changes a primary metric; no fifth/upper-expansion row.
- Each fatal provenance/QC/schema/input case produces no report, CSV, RDS, or
  partial directory; path collision preflight and overwrite false/true behavior
  are exact.
- Approved CSV filenames and read-back types/NA/status/provenance linkage are
  exact; no deprecated or unexpected files are written.
- Quarto source is thin, requires exactly one of config/RDS, selects only the
  canonical pH3 profile, does not suppress required warnings, and contains no
  scientific calculation, mapping, fallback, or automatic file selection.
- Regression tests prove unchanged legacy pH3, EdU canonical export/report,
  POI, general result saving, and no network/telemetry/credentials/destructive
  I/O or machine-specific production paths.

## Later user-run verification contract

Verification for this audit: **NOT RUN — local execution required**. No command
needs to be run now. After the later implementation, a laboratory member should
run from the canonical source directory above, with supported R and declared
package/test dependencies already installed:

1. `Rscript -e 'devtools::test(filter = "ph3-output-contract", stop_on_failure = TRUE)'`
2. `Rscript -e 'devtools::test(filter = "ph3-reports", stop_on_failure = TRUE)'`
3. Slice 4, Slice 5, legacy pH3, EdU, POI, result-API, and appearance/artifact
   regression filters named by the implementer.
4. `Rscript -e 'devtools::test(stop_on_failure = TRUE)'`
5. Build and `R CMD check --no-manual` only from a newly created external
   temporary directory, followed by `git diff --check` and full repository
   status inspection.

Expected success: zero failures/errors and no unreviewed warnings/notes; exact
approved temporary outputs only; no repository report, CSV, RDS, image, cache,
tarball, check directory, or other artifact. Stop on any value, denominator,
replicate-unit, provenance, QC, legacy, EdU, or POI regression; any partial
output; network attempt; experimental-input access; unexpected warning; or
unexpected generated file. Report only command, `PASS`/`FAIL`, failed test
names, a short error excerpt, and unexpected generated files.

## Independent review and resolution

Scientific-integrity review found no blocking or important finding. Across the
initial review and re-review it raised three advisory wording ambiguities: the
acquisition metrics could be read as optional rather than required retained
detail, and each of the two condition-summary rows used a generic replicate-
count label. The inventory now says the acquisition detail is required and
names both condition-level counts explicitly as biological-replicate counts.
The final scientific re-review found no remaining blocking, important, or
advisory findings and confirmed that no new estimand, replicate unit,
statistic, normalization, biological claim, or interpretation is introduced.

Artifact/security review found no blocking, important, or advisory findings.
It confirmed that the plan accurately treats current sequential writers as an
implementation risk rather than as atomic, requires complete validation and
collision preflight before directory creation or any write, and requires no
partial report/export set. It found no external transmission, machine-specific
production path, hidden fallback, destructive I/O, generated artifact, or
unrelated scope expansion.

## Assumptions, ambiguities, and exclusions

- The merged Slice 1–5 tables and methods are authoritative and unchanged.
- The older approval's six canonical file names predate the confirmed four
  Slice 5 tables; S6-D05 records the resulting ambiguity rather than inventing
  aliases or dropping tables.
- Condition/replicate/acquisition order comes only from explicit manifest
  indices; source-ID JSON lexicographic order is provenance serialization only.
- No QC acceptance threshold, reference condition, reference normalization,
  inferential statistic, error-bar convention, geometry requirement, or plot
  statistic is assumed.
- Derived 4N regions, downstream thresholds, alternative intervals, pooled or
  weighted aggregation, event-level diagnostics, GUI behavior, compatibility
  aliases, and historical artifact rewriting are out of scope.
- No code or test file was changed; no project code, R, Python, package, test,
  build/check, Quarto render, analysis, snapshot, or benchmark was run; no
  experimental input or generated artifact was read, modified, or created.

## Audit change report

- Canonical source directory: the absolute directory recorded above.
- Exact Git state: branch `feature/ph3-output-contract-audit`, baseline and
  current `HEAD` `ead9e13e18c8f6243e9622161e0c9b475a42da0b`; no commit created;
  one untracked audit memo and no other tracked/non-ignored change.
- Experimental inputs read or modified: none.
- Files modified: none (the memo is a new source-documentation file).
- Files generated: this audit memo only,
  `docs/implementation/ph3_slice6_output_contract_audit_2026-08-29.md`.
- Sample mappings, exclusions, gates, thresholds, normalization, statistics,
  or biological claims changed: none.
- Tests/checks/project execution: **NOT RUN — local execution required**;
  `git diff --check` is the only static text check and produced no output.
- Independent reviews: scientific wording advisories resolved and clean final
  re-review; artifact/security review clean.
- Unresolved uncertainty: scientific-owner Decisions S6-D01–S6-D06.
