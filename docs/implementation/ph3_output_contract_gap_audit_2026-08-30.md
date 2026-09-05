# pH3 confirmed output-contract gap audit

**Date:** 2026-08-30

**Status:** AUDIT COMPLETE — IMPLEMENTATION NOT AUTHORIZED

**Audience:** Code Architect and scientific owner

## Contract source and confirmation gate

The sole contract used for this audit is:

`/Users/sansamc/Documents/Codex/2026-08-30/ph3-output-contract-decisions/outputs/2026-08-30_final_ph3_output_qc_report_contract.md`

- Contract status: **OWNER-CONFIRMED** (`line 4`).
- Contract SHA-256:
  `1302f0aa96240849b3a0045372805918f25502586ce0a5bfeac687d9e6d99b0c`.
- The contract explicitly states that confirmation does not authorize code
  implementation (`lines 7 and 295-299`). This record is therefore audit and
  planning only.

No earlier decision draft or unconfirmed output memo was substituted for this
contract.

## Canonical source and baseline

Canonical source directory:

`/Users/sansamc/OMRF Dropbox/Chris Sansam/Synchd/Organized Dropbox/Projects and Data/Projects/Cowork_Projects/2026_Analyze_FACS_Data/facs_pseudocolor_workflow`

Audit branch: `feature/ph3-output-contract-audit`

Baseline/current `HEAD`: `ead9e13e18c8f6243e9622161e0c9b475a42da0b`

Baseline description: `Merge pull request #21 from SansamLab/feature/ph3-replicate-condition-aggregation`

At audit start, `HEAD`, local `main`, and `origin/main` all named this commit.
It is therefore the current merged baseline used for all classifications below.
Relevant merged implementation history is:

- `6fd4bb3` — pH3 biological-replicate and condition aggregation;
- `7accf2f` — pH3 acquisition metrics and sensitivity QC;
- `4a476bf` — pH3 event eligibility classification;
- `c2f7aaa` — pH3 production input validation and containment;
- `ccc2cd7` — pH3 export identity and provenance contract;
- `b963e26` / `3fdc46e` — older focused pH3 geometry/report path; and
- `e2a2fcd` — original pH3 mode.

The starting worktree had no tracked or staged change and one pre-existing,
untracked file:
`docs/implementation/ph3_slice6_output_contract_audit_2026-08-29.md`.
That file is not in the merged baseline, was not treated as authoritative
evidence, and was not modified.

## Scope and protected material

This audit compares the owner-confirmed pH3 contract with current merged source,
tests, configuration, report templates, and relevant Git history. It does not
authorize or perform any change to data handling, sample mapping, event
membership, gating, thresholds, normalization, statistics, biological claims,
production code, tests, configuration, reports, generated artifacts, or Git
state.

Experimental inputs accessed: **none**. No `.fcs`, `.wsp`, `.prism`, experimental
CSV, analysis RDS, generated report, figure, or prior experimental output was
opened. All conclusions come from text source, tests, configuration examples,
documentation, and read-only Git metadata.

## Files inspected

Governance and contract:

- workspace `AGENTS.md`;
- the owner-confirmed contract named above; and
- the pre-existing untracked 2026-08-29 audit memo, only to identify and avoid
  treating it as merged or authoritative.

Current pH3 science and provenance:

- `R/ph3.R`, `R/analysis.R`, `R/input.R`, `R/config.R`,
  `R/pseudocolor_helpers.R`, and `R/phase_quantitation.R`;
- `python/export_contract.py`, `python/export_flowjo_populations.py`,
  `python/export_flowjo_gate_geometry.py`; and
- `tools/flowjo-orchestration.R`.

Current report, plotting, artifact, and configuration surfaces:

- `R/ph3_geometry.R`, `R/results.R`, `R/artifacts.R`, `R/appearance.R`;
- `inst/quarto/facs_ph3_4n.qmd`, `inst/quarto/facs_complete.qmd`, and
  `inst/quarto/_report-setup.R`;
- `examples/config_ph3.yml`, `DESCRIPTION`;
- `docs/PH3_MODE.md`, `docs/REPORTS.md`, `docs/CONFIGURATION.md`, and
  `docs/PYTHON_INTERFACE.md`.

Focused tests and test-fixture guidance:

- `tests/testthat/test-ph3-event-classification.R`;
- `tests/testthat/test-ph3-acquisition-metrics.R`;
- `tests/testthat/test-ph3-replicate-condition-aggregation.R`;
- `tests/testthat/test-ph3-input-containment.R`;
- `tests/testthat/test-ph3.R`, `tests/testthat/test-results-api.R`,
  `tests/testthat/test-manifest.R`, and `tests/testthat/test-python-boundary.R`;
- `tests/python/test_export_contract.py`,
  `tests/python/test_create_synthetic_flowjo_index_fixture.py`, and
  `tests/python/test_flowkit_source_index_semantics.py`; and
- `tests/fixtures/SYNTHETIC_flowjo_index_semantics/README.md`.

## Classification method

- **Implemented as required:** current production path fulfills the confirmed
  clause and has direct code/test evidence.
- **Partially implemented:** a usable prerequisite or subset exists, but the
  complete confirmed behavior or user-visible obligation does not.
- **Absent:** no current production implementation was found.
- **Conflicts with contract:** current behavior would produce a scientifically
  or operationally different result, label, fallback, or output.

The classification applies to the production direct-identity pH3 path. A legacy
plot or report does not satisfy the production contract merely because it has a
similar name.

## Findings matrix

### Definitions, outcomes, and PHD-03 prevalence

| Contract item | Status | Current evidence and gap | Gap type |
|---|---|---|---|
| Primary mitotic fraction is 4N pH3-positive / `eligible_2to4n`, labeled “Analysis singlets (2N–4N)” | **Implemented numerically; report conflicts** | Production uses the 4N-positive numerator and `a_n <- sum(eligible_2to4n)` (`R/ph3.R:638-658`), with a fixed synthetic denominator test (`tests/testthat/test-ph3-acquisition-metrics.R:107-123`). The focused legacy report instead says “all Single Cell events” (`inst/quarto/facs_ph3_4n.qmd:68-74`; `R/ph3_geometry.R:313-318`). | Scientific core implemented; report label/route conflict |
| Four required outcomes A–D | **Partial** | A and B already exist as `ph3_4n_positive_prevalence_within_2to4n_percent` and `ph3_sub_4n_positive_prevalence_within_2to4n_percent` (`R/ph3.R:648-678`). No population-specific pH3 signal outcome C or D exists in the exact production tables (`R/ph3.R:839-844,1437-1445`). | Scientific method/output |
| Preserve approved pH3-positive membership and 4N/below-4N definitions | **Implemented as required** | Direct event identity supplies pH3 membership, while configured boundaries deterministically create eligible, below-4N, and 4N membership (`R/ph3.R:175-371`). Provenance and reconciliation fail closed (`R/ph3.R:393-585`). | Scientific core |
| Display Cells → Singlets → Analysis singlets → pH3-positive → 4N/below-4N, with provenance and `VERIFIED`/`COMPUTED` | **Absent** | Classification retains deterministic membership and geometry status fields (`R/ph3.R:274-347`), but no report displays the hierarchy or status labels. The focused report has neither `VERIFIED` nor `COMPUTED` labels (`inst/quarto/facs_ph3_4n.qmd:1-89`). | Report/geometry provenance |
| Both A and B use the same Analysis-singlets denominator | **Implemented as required** | Both production prevalence specifications use `a_n`; below-4N does not use the below-4N population as its denominator (`R/ph3.R:642-658`). | Scientific core |
| A and B partition total pH3-positive prevalence; percentages are not normalized | **Implemented as required** | Classification makes below-4N and 4N mutually exclusive and exhaustive within eligibility (`R/ph3.R:274-280,353-363`); prevalence is direct `100 * numerator / denominator` (`R/ph3.R:589-599`). | Scientific core |

### PHD-04 background correction, fallback, reference, and basis

| Contract item | Status | Current evidence and gap | Gap type |
|---|---|---|---|
| Per-sample linear regression on eligible pH3-negative events | **Absent; active label conflict** | pH3 normalization copies raw pH3 into `target_norm` and `target_bgsub` and records `g1_dna_only` (`R/ph3.R:46-60`). No pH3 background fit occurs in the production analysis path (`R/analysis.R:175-213`). Raw is therefore presently aliased under a background-subtracted field name. | Blocking scientific method / no-silent-fallback |
| Individual validity: at least 100 negatives; finite inputs; nonzero finite variation/full rank; finite coefficients/predictions; no positive-range extrapolation | **Absent** | No pH3 regression diagnostics or validity evaluator exists. Current `baseline_minimum_negative_events` is generic/EdU-era configuration and is not a pH3 regression implementation (`R/config.R:46-49`; `R/analysis.R:175-213`). | Blocking scientific method |
| C/D sample values are medians of background-corrected pH3 in positive 4N/below-4N populations | **Absent** | Production tables contain percentages, eligibility QC, sensitivity QC, and their percentage aggregates only (`R/ph3.R:630-844,1247-1450`). | Blocking scientific output |
| If any sample fails, fit one pooled negative model from all samples in the biological-replicate set and replace every individual model in that set | **Absent** | No set-wide fit, validity, replacement, or correction-tier field exists. Current `model_group` is an older technical grouping and is not used for pH3 signal correction (`R/pseudocolor_helpers.R:203-213`; `R/analysis.R:175-213`). | Blocking scientific method |
| Pooled model uses the same validity rules, including coverage across the set | **Absent** | No pooled pH3 regression QC exists. | Blocking scientific method |
| If pooling fails, switch the entire set to raw for C/D, displays, exports, reference ratios, and permitted within-basis tests | **Absent; current behavior conflicts** | All pH3 samples are presently raw from the outset and are not identified as a fallback tier (`R/ph3.R:46-60`). There is no set-wide tier decision or propagation to outputs. | Blocking scientific method |
| Never mix raw and corrected within a set; visibly label raw; export failure reason | **Absent; current behavior conflicts** | There is no analytical-basis or failure-reason field. Generic defaults call the display background-subtracted (`R/config.R:61-62`; `R/results.R:684-700,824-825`) even though pH3 is raw. | Blocking safeguard |
| Regression failure leaves A/B valid | **Partial structural support** | A/B are derived only from membership/count fields and currently do not depend on signal (`R/ph3.R:638-678`). No regression failure path or test yet proves that this independence is retained through the fallback engine. | Scientific safeguard |
| Explicit reference sample produces population-specific ratios on the set’s basis; no reference produces direct medians and no invented fold change | **Absent; configuration conflicts** | pH3 configuration rejects replicate references (`R/config.R:181-189`), and the example says pH3 has no reference (`examples/config_ph3.yml:54-55`). Generic ratios are condition-based and bypassed by the production pH3 branch (`R/analysis.R:223-235`; `R/results.R:308-340`). | Blocking scientific interface/method |
| Condition summaries and paired signal statistics never combine individual-corrected, pooled-corrected, and raw bases; mixed bases remain visible and stratified with absence reasons | **Absent** | Current production schemas have no signal basis, C/D value, stratum, or unavailable-summary reason (`R/ph3.R:1404-1445`). | Blocking scientific safeguard |
| Basis isolation does not affect A/B | **Absent as an explicit guard** | Existing A/B are independent, but no future-facing validator/test prevents signal-basis logic from filtering them. | Scientific safeguard |

### PHD-05 display-only offset and positive-domain safeguard

| Contract item | Status | Current evidence and gap | Gap type |
|---|---|---|---|
| Per-sample offset = median(raw eligible negative) − median(corrected eligible negative) | **Conflicts with contract** | Current helper chooses one analysis-wide scalar from anchors, baselines, or positive raw values and rounds it to a power of ten (`R/results.R:654-681`). Configuration also permits an arbitrary numeric offset (`R/config.R:251-259`). | Blocking display/scientific boundary |
| Offset affects plots only, never membership, medians, ratios, thresholds, exports, or tests | **Partially implemented** | `prepare_pseudocolor_signal()` copies the analysis and changes only the copy’s display `target_norm` (`R/results.R:684-700`), a useful isolation boundary. It does not implement the confirmed formula or provide invariance evidence for pH3 C/D, reference, export, and statistics. | Scientific safeguard |
| If any post-offset displayed value is ≤0, shift the entire plotted distribution so its minimum is exactly 1 | **Absent; current filtering conflicts** | No safeguard calculation exists. Log plotting filters nonpositive values from limits and displayed events (`R/pseudocolor_helpers.R:999-1017`) rather than translating the full affected distribution. | Blocking display safeguard |
| Safeguard is visibly flagged, exported with shift, and isolated from analytics | **Absent** | No pH3 QC field or report flag records a trigger, shift, or alignment exception. | Scientific QC/report/export |

### PHD-06 report and machine-readable output package

| Contract item | Status | Current evidence and gap | Gap type |
|---|---|---|---|
| Routine report emphasizes A–D, biological-replicate points, and condition summaries | **Absent** | The focused pH3 QMD shows one legacy 4N gate view (`inst/quarto/facs_ph3_4n.qmd:68-78`). The generic complete QMD prints every available quantitation object (`inst/quarto/facs_complete.qmd:61-97`) and is not a production fallback: its bundle dispatches to legacy pH3 plotting that expects `quantitation$ph3` (`R/artifacts.R:58-77`; `R/results.R:416-430`). | User-visible report |
| Routine report omits detailed counts, sample measurements, and regression diagnostics | **Conflicts with current templates** | The focused report prints a sample summary and full provenance (`inst/quarto/facs_ph3_4n.qmd:80-89`); the complete report dumps quantitation, input report, warnings, provenance, and appearance (`inst/quarto/facs_complete.qmd:70-120`). | User-visible report |
| Concise QC states correction tier, safeguards/fallbacks, geometry status/missingness, source locations, SHA-256 | **Partial inputs; output absent** | Input manifests retain paths and SHA-256 (`R/input.R:1238-1246`), but no curated report state contains all required QC/provenance fields. | Report/export |
| `sample_results.csv` exact role | **Absent** | No pH3 CSV projection/writer exists. The current CSV helper is EdU-only (`R/results.R:872-907,951-959`). | Machine-readable export |
| `condition_summary.csv` exact role, including individual values, basis strata, eligible mean/SD/statistics, and absence reasons | **Absent** | Existing pH3 condition tables contain only prevalence summaries and no C/D, bases, tests, or unavailable reasons (`R/ph3.R:1421-1445`). | Machine-readable export |
| `regression_qc.csv` exact diagnostics, tiers, reasons, offsets, safeguard | **Absent** | No pH3 regression or display QC table exists. | Machine-readable export |
| `provenance_manifest.json` relates report/exports/plots and records definitions, statuses, locations, checksums | **Absent** | Existing `export-manifest.json` is immutable **input** provenance, not the contract output manifest (`R/input.R:724-788,1238-1246`). | Machine-readable export/provenance |
| Canonical HTML is rendered from Quarto | **Partially implemented infrastructure; canonical profile absent** | Quarto HTML templates and exactly-one-input helper exist (`inst/quarto/_report-setup.R:1-8`), but no production contract pH3 QMD exists. The focused QMD is legacy-incompatible and dereferences `analysis$quantitation$ph3` (`inst/quarto/facs_ph3_4n.qmd:87`), which production intentionally omits (`tests/testthat/test-ph3-acquisition-metrics.R:353-403`). | Report |
| Every standard plot exported individually as PDF and PNG | **Absent** | `save_facs_results()` writes at most one aggregate PDF and one aggregate PNG (`R/results.R:934-1005`), not a registered per-plot pair. | Artifact export |
| Export condition panels, per-sample pseudocolor, per-sample three-panel gating, and per-sample DNA distributions | **Partial** | Generic per-sample DNA-versus-pH3 plotting exists (`R/results.R:717-826`). Production A–D panels, FSC/SSC + singlet + analysis gating overviews, and DNA distributions are absent. | Plot/report |
| Output paths cannot overlap or alias immutable inputs/provenance and cannot escape an owned package root | **Absent as an output-package safeguard** | Current writers accept arbitrary resolved paths and compare only requested outputs/existing files (`R/results.R:855-979`; `R/artifacts.R:16-23`). They do not compare destinations with analysis inputs, immutable export-operation artifacts, geometry supplements, configuration/source paths, or symlink-resolved targets. | Artifact/security integrity |

### PHD-07 statistics

| Contract item | Status | Current evidence and gap | Gap type |
|---|---|---|---|
| Only explicitly configured treatment-versus-control comparisons; pair by biological-replicate set; complete pairs only | **Absent** | The config key list has no pH3 comparison structure (`R/config.R:5-29`), and production pH3 stops after descriptive tables (`R/analysis.R:223-228`). | Blocking scientific method/interface |
| Signal pairs must share one analytical basis; mixed bases suppress the overall statistic but retain stratified values | **Absent** | No signal basis or comparison layer exists. | Blocking scientific safeguard |
| Two-sided paired t-test; at least 3 complete pairs; exact untransformed A/B or applicable C/D median/reference ratio; display transforms excluded; no small-n test switching | **Absent** | No production pH3 inferential-statistics implementation or focused test exists. | Blocking scientific method |
| Four independent Holm families, exact adjusted p-values, approved symbols | **Absent** | No pH3 Holm adjustment, symbol mapping, or bracket model exists. | Blocking scientific method/report |
| With fewer than 3 pairs, descriptive only and no p-value/symbol/bracket | **Absent** | No pH3 statistical availability state exists. | Scientific safeguard/report |

### PHD-08 canonical report profile and geometry

| Contract item | Status | Current evidence and gap | Gap type |
|---|---|---|---|
| Ordered A–D condition summary with points/means; SD only at n≥3; eligible stats brackets only | **Absent** | Production condition tables currently cover percentages only (`R/ph3.R:1331-1445`), and no canonical plot/report model exists. Current tables compute descriptive SD at n≥2 (`R/ph3.R:1345-1359`); a future report must independently enforce the display rule without deleting exportable values. | Report plus upstream C/D/stats |
| Per-sample DNA-versus-pH3 pseudocolor | **Partial** | Generic plots exist (`R/results.R:717-826`), but pH3 signal basis, contract offset, safeguard, and visible QC labels are absent. | Plot/scientific-display boundary |
| Per-sample three-panel Cells, Singlets, DNA/pH3 gating overview | **Absent** | Current focused plot is only DNA/pH3 with one pH3 polygon (`R/ph3_geometry.R:230-319`). Production input validation requires exactly `complete`, `g1`, and `ph3_positive` (`R/input.R:584-588`), the example maps only Single Cells/G1/pH3 Positive (`examples/config_ph3.yml:42-51`), and geometry orchestration defaults to `ph3_positive` only (`tools/flowjo-orchestration.R:230-234`). Required Cells/Singlets event backgrounds, channels, and verified gate geometry cannot be inferred. | Upstream provenance/plot |
| Per-sample DNA-content distribution | **Absent** | No production pH3 DNA-distribution plot builder was found. | Plot |
| Concise QC/provenance appendix | **Absent** | Current QMD prints uncurated objects and suppresses warnings globally (`inst/quarto/facs_ph3_4n.qmd:3-4,80-89`). | Report/QC |
| Manifest of accompanying exports | **Absent** | No output artifact registry/manifest exists. | Report/export |
| Verified source gates labeled `VERIFIED`; computed analysis regions shown and labeled `COMPUTED`; sources stated; never conflate them | **Conflicts with contract** | R validates a supplied geometry CSV but does not consume its `geometry-manifest.json`, detached digest, artifact hash, workspace/operation binding, or verified status (`R/ph3_geometry.R:15-111`). The report nevertheless calls it exact FlowJo geometry without status labels (`inst/quarto/facs_ph3_4n.qmd:68-89`). Computed regions are not labeled. | Blocking provenance/report safeguard |
| If verified geometry is missing, issue valid numbers, retain computed displays, mark report incomplete, record status, never reconstruct verified geometry | **Conflicts with contract** | The QMD stops if geometry is absent (`inst/quarto/facs_ph3_4n.qmd:42-48`), and the reader stops if any sample is missing (`R/ph3_geometry.R:45-60`). The focused plot disables computed phase overlays (`R/ph3_geometry.R:256-259`), so missing verified geometry suppresses computed displays and the report. | Blocking geometry/report safeguard |
| Multi-replicate verified geometry can be resolved without losing provenance | **Conflicts with current integration** | Geometry orchestration writes one fixed, immutable `gate_geometry.csv` per export operation/biological replicate (`tools/flowjo-orchestration.R:266-317`), while the QMD accepts one CSV and the reader requires it to cover every analysis prefix (`inst/quarto/facs_ph3_4n.qmd:8,43-48`; `R/ph3_geometry.R:45-60`). Manual concatenation would discard per-operation hash/linkage and is prohibited as a fallback. | Blocking provenance integration |
| Analysis provenance truthfully represents verified geometry | **Conflicts with current integration** | Production input validation requires the base input manifest’s `geometry_overlay_status` to be exactly `not_requested` (`R/input.R:484-556`), and classification copies it (`R/ph3.R:292-335`). A separately verified geometry supplement exists in Python (`python/export_flowjo_gate_geometry.py:229-303`; `python/export_contract.py:487-514`), but R does not consume it. | Blocking provenance integration |

### Standard, conditional, and deferred profile behavior

| Contract item | Status | Current evidence and gap | Gap type |
|---|---|---|---|
| All standard plots, concise appendix, HTML, and four CSV/JSON outputs are routine | **Absent except partial pseudocolor infrastructure** | No canonical production report/output registry exists. | User-visible output |
| Reference ratios only with explicit reference | **Absent; current pH3 reference prohibition conflicts** | See PHD-04 reference row above. | Scientific/interface |
| Brackets only for eligible explicit comparisons | **Absent** | No pH3 comparison/statistics layer exists. | Scientific/report |
| Pooled/raw tier and mixed-basis stratification appear conditionally and visibly | **Absent** | No tier/basis output exists. | Scientific/report/export |
| Missing verified geometry causes incomplete mark but preserves computed display | **Conflicts** | Current focused report stops and disables computed overlays. | Report/provenance |
| Sensitivity and Unassigned presentation remain deferred and are not silently added to the canonical profile | **Partial legacy-surface conflict risk** | Generic pH3 plotting automatically offers assignment and sensitivity plots (`R/results.R:420-430`), and legacy documentation treats them as supported outputs (`docs/PH3_MODE.md:67-80`). They are not in the focused QMD, but a new canonical profile must not reuse the generic bundle wholesale. | Report scope/non-goal |

## Scientific-method gaps versus presentation-only gaps

### Scientific-method and scientific-integrity blockers

1. PHD-04 is not implemented: no individual fit, strict validity, pooled
   set-wide replacement, raw set-wide fallback, C/D medians, basis, or reason.
2. Raw pH3 is presently aliased as `target_bgsub` and can be described as
   background-subtracted. This is a no-silent-fallback and truthful-labeling
   conflict, even though current production numerical tables do not yet report
   signal outcomes.
3. The mapping from current technical acquisitions to the contract’s
   “sample-level median” is not owner-confirmed. Pooling events, averaging
   acquisition medians, or another rule would produce different values and
   must not be inferred.
4. Exact pH3 reference configuration semantics and the behavior for a missing,
   zero, or nonfinite population-specific reference median are unresolved.
5. Comparison configuration and designated-control mapping are absent; PHD-07
   methods cannot be safely attached to inferred pairs.
6. Basis-stratified aggregation and statistical exclusion rules are absent.
7. PHD-05’s formula and positive-domain safeguard are absent. Current silent
   exclusion of nonpositive log-display values violates the confirmed display
   safeguard.
8. Verified versus computed geometry provenance cannot currently be stated
   truthfully by the R production/report path.

### Report, export, artifact, and UI gaps

1. No canonical four-outcome pH3 report profile or ordered report model exists.
2. All four required CSV/JSON output files and their exact projections are
   absent.
3. No named standard-plot registry or individual PDF+PNG export exists.
4. Three-panel gating overviews and DNA-content distributions are absent.
5. The current focused report is legacy-only, requires a geometry CSV that the
   R path does not verify against its immutable supplement, uses the wrong
   denominator text, suppresses warnings, and cannot consume a production
   direct-identity object.
6. No curated QC appendix, incomplete-report banner, export inventory, or
   output checksum manifest exists.
7. The configuration/example/docs/configurator do not expose confirmed
   reference/comparison/output behavior and retain incompatible legacy wording.

These presentation gaps must not be filled by recalculating science inside
Quarto or by renaming legacy fields. Report and export code must consume one
validated scientific output model.

## Dependency-ordered implementation-slice plan

The order below is driven by user-visible contract obligations and the data and
provenance dependencies needed to fulfill them. It is not an authorization to
implement.

### Decision gate D0 — close unresolved interfaces before code

**Scope:** Obtain explicit owner/architect decisions for: (a) technical-
acquisition to sample-level C/D aggregation; (b) a stable software key that
preserves the already-explicit top-level biological-replicate grouping without
using `technical_replicate` or file order as biology; (c) exact reference
declaration and invalid/zero reference-median behavior; (d) empty or nonfinite
positive-population C/D behavior, including stop versus explicit unavailable
status; (e) raw-tier display-offset semantics when no valid corrected-negative
distribution exists; (f) comparison/designated-control configuration and the
unavailable-result behavior for degenerate/nonfinite paired tests; (g) required
event populations and exact FSC/SSC and DNA width/area channel identifiers for
the first two gating panels; (h) multi-operation geometry input declaration;
(i) exact HTML/plot naming and directory layout; and (j) exact long-versus-wide
schemas/column names for the four exports. Decide whether optional within-basis
statistics are in the first implementation; omission is the conservative
default unless explicitly requested.

**Inputs/outputs:** Owner-confirmed addendum or architect interface record;
no analysis output.

**Non-goals:** No new gate, threshold, reference, comparison, filename,
fallback, or aggregation default may be inferred from current file order or
legacy behavior.

### Slice 1 — validated pH3 output/configuration model

**Boundary:** Add explicit schemas and fail-closed validation for the four
outcomes, correction basis/reasons, reference state, comparisons, geometry
state, and required output registry. Preserve current membership and A/B values
unchanged.

**Inputs → outputs:** Current classification, manifest, Slice 4/5 tables, and
approved D0 configuration → one validated in-memory pH3 contract model with no
file writes.

**Likely sources:** `R/config.R`, `R/analysis.R`, new `R/ph3_outputs.R`,
`examples/config_ph3.yml`.

**Focused tests:** New `tests/testthat/test-ph3-output-model.R`; exact keys,
types, order, comparison/reference validation, missing/extra/mixed fields,
stable replicate-set identity, and byte-for-byte preservation of A/B and event
membership. Use only in-memory, explicitly `SYNTHETIC` fixtures.

**Local-only command:**
`Rscript -e 'devtools::test(filter = "ph3-output-model", stop_on_failure = TRUE)'`

**Non-goals:** No regression, statistics, plot, Quarto render, file write, or
legacy-field reinterpretation.

### Slice 2 — regression validity and set-wide fallback engine

**Boundary:** Implement PHD-04 fitting and tier selection from
`eligible_2to4n & !ph3_positive_member`; produce event corrections and complete
individual/pooled diagnostics without changing positivity or DNA membership.

**Inputs → outputs:** Validated event classification, raw DNA/pH3, explicit
replicate-set mapping → per-sample and pooled fit diagnostics, one set-wide
basis decision, corrected or raw event signal, failure reasons.

**Likely sources:** `R/ph3.R`, `R/analysis.R`, possibly new
`R/ph3_background.R`.

**Focused tests:** New `test-ph3-background-regression.R`: 99 versus 100
negatives; nonfinite DNA/pH3; zero variation/rank failure; nonfinite
coefficient/prediction; exact coverage equality versus a coverage gap; any one
sample failure forcing pooled replacement for the entire set; pooled failure
forcing raw for the entire set; exact basis/reason propagation; A/B unchanged.
Also cover empty 4N/below-4N positive populations and nonfinite positive-event
pH3/predictions using the D0-approved stop or explicit-unavailable behavior;
never silently omit or substitute events. All data must be deterministic and
labeled `SYNTHETIC`.

**Local-only command:**
`Rscript -e 'devtools::test(filter = "ph3-background-regression", stop_on_failure = TRUE)'`

**Non-goals:** No reference ratios, condition summaries, statistics, display
offset, report, or file export.

### Slice 3 — C/D sample results, reference handling, and basis-safe summaries

**Boundary:** Compute population-specific sample medians using the D0-approved
acquisition rule; apply optional reference ratios only after basis resolution;
build biological-replicate and condition projections that never combine bases.

**Inputs → outputs:** Slice 2 event/basis/QC plus current A/B and explicit
reference declaration → complete sample results and basis-stratified condition
values, means/SD only when comparable, and explicit absence reasons.

**Likely sources:** `R/ph3.R`, new `R/ph3_outputs.R`.

**Focused tests:** New `test-ph3-signal-basis-comparability.R`: C/D population
membership, direct medians with no reference, population-specific reference
ratios, D0-approved invalid-reference behavior, no mixed-basis mean/SD,
same-basis eligibility, raw values retained visibly, exact absence reasons,
A/B unaffected, empty/nonfinite positive outcomes, and no row-order inference.

**Local-only command:**
`Rscript -e 'devtools::test(filter = "ph3-signal-basis-comparability", stop_on_failure = TRUE)'`

**Non-goals:** No inferential test, plot, display offset, or file write.

### Slice 4 — configured paired statistics and annotation model

**Boundary:** Implement PHD-07 only for explicit configured comparisons and
complete biological-replicate pairs. Produce a data model for report brackets;
do not draw them here.

**Inputs → outputs:** Validated comparison config and Slice 3 analytical values
→ paired-test eligibility, exact raw/adjusted p-values, approved symbols,
pair counts, basis stratum, and unavailability reasons.

**Likely sources:** new `R/ph3_statistics.R`, `R/config.R`,
`R/ph3_outputs.R`.

**Focused tests:** New `test-ph3-statistics.R`: pair by stable replicate key,
not row order; incomplete pairs; n=2 versus n=3; two-sided paired t-test;
untransformed A/B; applicable C/D direct median or ratio; display fields
excluded; four independent Holm families; exact boundary symbols; mixed-basis
suppression; degenerate/zero-variance differences and nonfinite test results
produce the D0-approved unavailable/error state with no alternate test; Holm
family membership remains exact when a configured comparison is unavailable;
no small-sample normality-driven test switch.

**Local-only command:**
`Rscript -e 'devtools::test(filter = "ph3-statistics", stop_on_failure = TRUE)'`

**Non-goals:** No unconfigured comparison, alternate test, confidence interval,
normality-test switching, plot, or file write.

### Slice 5 — display-only offset and positive-domain QC

**Boundary:** Implement the per-sample negative-median offset and deterministic
minimum-one safeguard in a plot-only projection. Record offset, trigger, shift,
and alignment exception for QC.

**Inputs → outputs:** Slice 2 analytical event signal plus the exact eligible
negative mask and affected plot-event set → display values and display QC only.

**Likely sources:** `R/results.R`, `R/pseudocolor_helpers.R`, new
`R/ph3_plots.R` or `R/ph3_outputs.R`.

**Focused tests:** New `test-ph3-display-offset.R`: exact formula per sample;
negative-median alignment when no safeguard triggers; full-distribution shift
and minimum exactly 1 when triggered; visible QC fields; original analysis
object unchanged; membership, analytical medians, ratios, thresholds, A/B,
statistics, and analytical export model byte-for-byte invariant. Exercise the
D0-approved raw-tier offset representation together with the visible
`RAW—NOT BACKGROUND SUBTRACTED` label, failure reason, and safeguard behavior;
do not silently treat raw as corrected or assume an offset of zero.

**Local-only command:**
`Rscript -e 'devtools::test(filter = "ph3-display-offset", stop_on_failure = TRUE)'`

**Non-goals:** No configurable scientific offset, power-of-ten rounding,
silent filtering of nonpositive plot values, or modification of analytical
columns.

### Slice 6 — verified/computed geometry-set provenance

**Boundary:** Resolve all per-operation geometry supplements for the analysis;
verify detached digests, artifact hashes, workspace, operation, sample, gate,
and channel bindings; represent source gates as `VERIFIED`; represent
deterministic analysis regions as `COMPUTED`; treat missing verified geometry
as an incomplete-report state while preserving computed regions.

**Inputs → outputs:** Explicit multi-operation geometry locations, immutable
input manifests/supplements, approved channel/population mapping, computed
membership/boundaries → one per-sample gate/region registry with status, source,
and missingness reasons.

**Likely sources:** `python/export_contract.py`,
`python/export_flowjo_gate_geometry.py`, `tools/flowjo-orchestration.R`,
`R/input.R`, `R/ph3_geometry.R`, `R/ph3_outputs.R`.

**Focused tests:** New `test-ph3-geometry-provenance.R` plus focused Python
contract tests: multi-operation resolution; mutated JSON/digest/artifact;
wrong workspace/operation/sample/gate/channels; partial/missing verified gates;
computed regions remain present and labeled; no reconstruction or manual CSV
concatenation fallback. Fixtures must be explicitly `SYNTHETIC` and test-only.

**Local-only commands:**

- `Rscript -e 'devtools::test(filter = "ph3-geometry-provenance", stop_on_failure = TRUE)'`
- `python3 -m unittest -v tests/python/test_export_contract.py`

**Non-goals:** No gate reconstruction from event extrema, no status promotion
from `COMPUTED` to `VERIFIED`, no geometry-driven change to event membership,
and no reading original experimental FCS/WSP in unit tests.

### Slice 7 — canonical plot registry and report model

**Boundary:** Build the ordered A–D panels, per-sample pseudocolor, three-panel
gating overviews, DNA distributions, concise QC appendix model with required
source locations and SHA-256 checksums, incomplete banner state, and
accompanying-file inventory from validated prior slices.

**Inputs → outputs:** Slices 1–6 model → named editable plot registry and
human-readable report model; no writes.

**Likely sources:** new `R/ph3_plots.R`, `R/ph3_outputs.R`,
`R/ph3_geometry.R`, `R/results.R`, `R/artifacts.R`.

**Focused tests:** New `test-ph3-report-model.R`: exact section/plot order;
biological-replicate points and means; SD whiskers only at n≥3; eligible
brackets only; all basis/raw/incomplete/QC labels; `VERIFIED`/`COMPUTED`
sources; source locations and SHA-256 checksums; computed display survives
missing verified geometry; no detailed counts/regression table in human model;
no sensitivity or Unassigned presentation; plot-ready data exactly match
analytical/display projections.

**Local-only command:**
`Rscript -e 'devtools::test(filter = "ph3-report-model", stop_on_failure = TRUE)'`

**Non-goals:** No scientific calculation inside plotting, no sensitivity or
Unassigned presentation, and no change to legacy pH3/EdU/POI plots.

### Slice 8 — machine-readable projections and complete artifact plan

**Boundary:** Build the exact in-memory projections for `sample_results.csv`,
`condition_summary.csv`, and `regression_qc.csv`; define the final
`provenance_manifest.json` schema—including analysis definitions, gate
provenance, `VERIFIED`/`COMPUTED` status, source locations, SHA-256 checksums,
and report/export/plot relationships—and register every required artifact. Do
not write or finalize any file in this slice.

**Inputs → outputs:** Validated output/report/plot/geometry models plus explicit
artifact naming/layout decisions → exact CSV tables and one complete planned
artifact registry, including the future HTML and every plot PDF/PNG.

**Likely sources:** new `R/ph3_exports.R`, `R/results.R`, `R/artifacts.R`.

**Focused tests:** New `test-ph3-output-contract.R`: exact names, schemas,
ordering, and types; basis and reason retention; all required report/plot/export
relationships present in the plan; exact analysis definitions, gate provenance,
status, source-location, and checksum fields; no index column; no unexpected
artifact; report/export analytical values never contain display offsets or
safeguard shifts except their QC fields.

**Local-only command:**
`Rscript -e 'devtools::test(filter = "ph3-output-contract", stop_on_failure = TRUE)'`

**Non-goals:** No file write, render, checksum finalization, automatic
destination, or export of deferred outputs.

### Slice 9 — canonical Quarto profile and user-facing configuration/docs

**Boundary:** Add a new production pH3 QMD rather than silently repurposing
`facs_ph3_4n.qmd`; render the exact ordered report; expose the approved explicit
configuration; update documentation and configurator/migration guidance while
preserving legacy semantics.

**Inputs → outputs:** Validated analysis RDS or explicit config and Slice 8
report model/registry → renderable canonical HTML source and a controlled
staging-tree HTML render. This slice does not publish the complete package or
finalize the output manifest.

**Likely sources:** new `inst/quarto/facs_ph3.qmd`,
`inst/quarto/_report-setup.R`, `R/config.R`, `R/configurator.R`,
`examples/config_ph3.yml`, `docs/PH3_MODE.md`, `docs/REPORTS.md`,
`docs/CONFIGURATION.md`, `docs/MIGRATION.md`.

**Focused tests:** New `test-ph3-report-contract.R`: exactly one config/RDS;
thin QMD with no science/fallback/mapping logic; exact section order; concise
QC; visible incomplete mark; warnings not globally hidden; exact artifact
inventory and links; local render from an explicitly `SYNTHETIC` test-only
analysis in an external temporary directory.

**Local-only command:**
`Rscript -e 'devtools::test(filter = "ph3-report-contract", stop_on_failure = TRUE)'`

**Non-goals:** No rewriting historical reports/artifacts, no compatibility
alias that changes legacy meaning, no hidden example data, and no sensitivity
or Unassigned sections.

### Slice 10 — transactional complete-package finalizer

**Boundary:** Create a canonical, destination-parent-local sibling staging
directory on the same filesystem; preflight the complete inventory and all
protected paths; write the three CSVs and every plot PDF/PNG; render the HTML;
then compute hashes and write `provenance_manifest.json` last. Validate the
entire staged package before one no-clobber publication operation. A concurrent
claim must fail without replacing either package.

**Inputs → outputs:** Slices 8–9 projections, plot registry, renderable report,
explicit destination, and protected input/provenance/config/source path inventory
→ exactly one complete HTML/CSV/JSON/PDF/PNG package with verified relationships
and SHA-256 checksums. The canonical package finalizer is no-clobber only; an
existing destination requires a different explicit empty destination rather
than replacement.

**Likely sources:** new `R/ph3_exports.R`, `R/results.R`, `R/artifacts.R`, and
`inst/quarto/_report-setup.R`.

**Focused tests:** New `test-ph3-output-artifacts.R`: canonicalized destination
and sibling staging stay under the destination parent and on the same
filesystem; reject symlink/path escapes and any overlap/alias with experimental
inputs, export-operation directories, input manifests, geometry artifacts,
configuration, or source; exact inventory; duplicate paths; checksum/link
mutation; any pre-existing destination is rejected; failure injection before
each finalization stage; concurrent destination claim; no partial publication;
no broad recursive deletion; no unexpected file. The manifest must include the
final HTML, analysis definitions, gate provenance, `VERIFIED`/`COMPUTED` status,
source locations, and every plot/export relationship and checksum before
publication.

**Local-only command:**
`Rscript -e 'devtools::test(filter = "ph3-output-artifacts", stop_on_failure = TRUE)'`

**Non-goals:** No automatic destination, timestamp fallback, overwrite or
replacement of an existing package/directory, direct sequential writes into the
published package, cleanup outside the owned staging directory, or export of
deferred outputs.

## Local-only verification contract

Verification for this audit is **NOT RUN — local execution required**.

For future implementation, run from the canonical source directory recorded
above. Prerequisites are R >= 4.2, the dependencies already declared in
`DESCRIPTION`, an already-installed `devtools`, Python 3 plus the existing
exporter test dependencies when Slice 6 changes Python, and Quarto only for the
canonical report/render slice. Do not install packages or access a network as
part of verification. Use only deterministic, unmistakably `SYNTHETIC`
test-only fixtures; do not use experimental inputs.

Run the focused command recorded under each slice, then the existing pH3
regression filters in dependency order:

1. `Rscript -e 'devtools::test(filter = "ph3-input-containment", stop_on_failure = TRUE)'`
2. `Rscript -e 'devtools::test(filter = "ph3-event-classification", stop_on_failure = TRUE)'`
3. `Rscript -e 'devtools::test(filter = "ph3-acquisition-metrics", stop_on_failure = TRUE)'`
4. `Rscript -e 'devtools::test(filter = "ph3-replicate-condition-aggregation", stop_on_failure = TRUE)'`
5. `Rscript -e 'devtools::test(filter = "ph3", stop_on_failure = TRUE)'`
6. `Rscript -e 'devtools::test(stop_on_failure = TRUE)'`

Expected success for every focused/full test command: zero failures/errors and
only explicitly asserted warnings. Expected repository-generated files: none.
The future report/output tests may create the exact declared HTML,
CSV/JSON/PDF/PNG package only inside a controlled external temporary directory;
they must validate the inventory and checksums before cleanup. A future direct
Quarto command cannot safely be named until Slice 9 establishes the canonical
QMD name and required parameters; the test must own that render in a temporary
directory.

After all tests, a laboratory member may build and check from a new external
temporary directory and then inspect repository status:

```bash
ph3_repo='/Users/sansamc/OMRF Dropbox/Chris Sansam/Synchd/Organized Dropbox/Projects and Data/Projects/Cowork_Projects/2026_Analyze_FACS_Data/facs_pseudocolor_workflow'
ph3_check_dir="$(mktemp -d)"
(cd "$ph3_check_dir" && R CMD build "$ph3_repo")
R CMD check --no-manual "$ph3_check_dir"/facspseudocolor_*.tar.gz --output="$ph3_check_dir"
git -C "$ph3_repo" diff --check
git -C "$ph3_repo" status --short
```

Expected generated files are the package tarball and check directory only
under `ph3_check_dir`; expected repository-generated files are none. Expected
package-check success is `Status: OK`; expected `git diff --check` output is
empty. Stop on any error, warning/note not explicitly reviewed, multiple or
missing tarballs, or repository file not intentionally in scope.

Stop and report `FAIL` on any changed A/B value or event membership; inferred
sample/reference/comparison mapping; mixed signal basis; silent fallback;
display value entering analytics; missing computed geometry when provenance
exists; false `VERIFIED` label; missing/extra/partial artifact; checksum
mismatch; experimental-input access; network attempt; unexpected warning/note;
or unexpected repository file. Report only the command, `PASS` or `FAIL`, failed
test names, a short relevant error excerpt, and unexpected generated files.

## Unresolved decisions and blockers

The following block implementation and must not be guessed:

1. The scientific rule converting one or more current technical acquisitions
   into the contract’s sample-level C/D median.
2. The biological grouping itself is already the explicit top-level replicate
   mapping (`replicate`/`replicate_index`), not a new scientific decision. The
   architect must still define and validate a stable software key that preserves
   that mapping and never promotes `technical_replicate`, `model_group`, or file
   order to biological grouping.
3. The explicit reference declaration and the required behavior for missing,
   zero, or nonfinite 4N or below-4N reference medians.
4. The required stop-versus-explicit-unavailable status and report/export
   representation when a positive population is empty or contains nonfinite
   raw/corrected signal needed for C/D.
5. Raw-tier display-offset semantics when pooled correction is invalid and no
   valid corrected-negative distribution exists.
6. The comparison/designated-control configuration schema and exact validation
   of treatment-versus-control direction.
7. The unavailable/error representation for degenerate or nonfinite paired
   t-test results and their membership in a Holm family; no alternate test may
   be selected silently.
8. Exact required source populations and channel identifiers for FSC/SSC Cells
   and DNA width/area Singlets panels; current inputs are insufficient to infer
   them safely.
9. The multi-operation geometry-set declaration and artifact format. Manual
   concatenation of per-operation CSV files is not acceptable provenance.
10. Exact output schema columns/layout (including long versus wide
   `condition_summary.csv`), canonical HTML/plot names, and output-directory
   layout.
11. Whether optional within-basis statistics are part of the first canonical
   implementation. The contract permits but does not require them; do not add
   them silently.

These are decision blockers, not invitations to use legacy behavior or nearby
files as defaults.

## Independent review and self-review

An independent scientific-integrity reviewer inspected the confirmed contract,
current science/configuration code, and focused tests. It confirmed the PHD-03
production numerator/denominator implementation and identified blocking absence
of PHD-04, PHD-05, PHD-07, basis isolation, and reference handling, plus the raw
signal mislabel and legacy report denominator conflict. It also identified the
technical-acquisition/sample-median and invalid-reference decisions above.

An independent artifact/security reviewer inspected report, geometry, export,
writer, and test sources. It confirmed the absent output package/report profile,
the verified/computed and missing-geometry conflicts, the multi-operation
geometry integration conflict, warning suppression, and the absence of an
atomic complete-package finalizer. It found no network upload, telemetry, or
credential behavior in the pH3 analysis/report path. External `system2()` calls
are confined to explicit FlowJo-preparation orchestration, not the installed
analysis/report path. Existing writers refuse overwrite by default, but future
multi-file publication still needs staging, complete-inventory validation, and
scope-confined cleanup.

The reviewers then inspected the complete draft record. Their blocking
finalization-order finding was resolved by making Slice 8 projection-only,
keeping Slice 9 as an unpublished staged render, and adding Slice 10 to write
the output manifest last and publish only a validated complete package. Their
important findings were resolved by removing an unsupported “verified” label
from the current R plot, preserving the already-explicit biological-replicate
grouping while requiring a stable software key, adding empty/nonfinite positive
outcomes, raw-tier offset semantics, and degenerate statistics to D0/tests,
requiring protected-input and symlink/path preflight, making provenance fields
explicit in Slices 7–10, and adding exact external build/check commands. The
same-filesystem sibling staging and concurrent no-clobber requirements resolve
the transactional advisory.

The orchestrator then rechecked every matrix row against the owner-confirmed
contract and cited surrounding source. The audit does not treat a legacy path,
input manifest, or pre-existing untracked memo as proof of a production output
contract.

## Explicit handback to the Code Architect

The merged baseline has a strong, fail-closed production identity,
classification, A/B prevalence, QC, and percentage-aggregation foundation.
PHD-03 is numerically implemented. The confirmed output contract is not a thin
reporting change: C/D background correction/fallback/basis behavior, PHD-05,
reference semantics, PHD-07, truthful geometry status, missing-geometry
behavior, four exports, plot registry, and canonical report are absent or in
direct conflict.

Recommended next handoff: first obtain the D0 scientific/interface decisions,
then authorize a narrowly scoped implementation of Slices 1–3 with separate
scientific-integrity and artifact/security review. Do not begin with Quarto or
rename legacy fields: that would make user-visible output appear complete before
the analytical basis and provenance safeguards exist.

## Audit change report

- Canonical source directory: recorded above.
- Experimental inputs read: none; modified: none.
- Files modified: none.
- Files generated: this new audit record only,
  `docs/implementation/ph3_output_contract_gap_audit_2026-08-30.md`.
- Pre-existing untracked file preserved unchanged:
  `docs/implementation/ph3_slice6_output_contract_audit_2026-08-29.md`.
- Sample mappings, exclusions, gates, thresholds, normalization, statistics,
  or biological claims changed: none.
- Project code/tests/R/Python/Quarto/analysis/render execution:
  **NOT RUN — local execution required**.
- Read-only/static checks: contract SHA-256 computed; document whitespace scan
  clean; final `HEAD`, branch, staged diff, tracked diff, and worktree status
  rechecked. No project code or test command was used for these checks.
- Commit, push, PR, merge, branch move, or Git-index change: none.
