# EdU output-contract gap audit and implementation plan

**Audit date:** 2026-08-30

**Status:** READ-ONLY AUDIT COMPLETE — IMPLEMENTATION NOT AUTHORIZED BY THIS RECORD

**Verification:** `NOT RUN — local execution required`

## Source contract and authority

The source-of-truth record for this audit is:

`/Users/sansamc/Documents/Codex/2026-08-30/edu-output-first-design/outputs/edu-desired-analyses-and-figures-specification-2026-08-30.md`

Its status is **Scientific-owner confirmed; ready for Code Architect handoff**.
It authorizes documentation and design only, not implementation. It confirms the
every-sample pseudocolor figure, equal-biological-replicate DNA-density figure,
shared provenance contract, and per-sample display-offset convention. It also
incorporates without changing the previously confirmed EdU output definitions.

The incorporated antecedent record inspected for those definitions was:

`../edu_output_contract_approved_2026-08-23.md`

That antecedent confirms the seven canonical EdU tables, their populations and
denominators, stable names, rejected intensity outputs, acquisition-level
retention, and compatibility rules. The 2026-08-30 record supersedes its old
instruction to preserve the former automatic display-offset calculation.

## Canonical baseline and repository state

- Canonical source directory:
  `/Users/sansamc/OMRF Dropbox/Chris Sansam/Synchd/Organized Dropbox/Projects and Data/Projects/Cowork_Projects/2026_Analyze_FACS_Data/facs_pseudocolor_workflow`
- Current branch: `feature/ph3-output-contract-audit`
- Baseline commit: `ead9e13e18c8f6243e9622161e0c9b475a42da0b`
- Baseline subject: `Merge pull request #21 from SansamLab/feature/ph3-replicate-condition-aggregation`
- Baseline date: `2026-08-29T16:14:22-05:00`
- Local `HEAD`, `main`, and locally available `origin/main` resolve to the same
  commit; no remote fetch was performed.
- The approved EdU table contract entered history at
  `c0bd0f997dc5ee1e6aa4206364e87518a520fb64` and was merged by `5de77a9`.
- Two unrelated untracked pH3 audit documents were present before this audit:
  `docs/implementation/ph3_output_contract_gap_audit_2026-08-30.md` and
  `docs/implementation/ph3_slice6_output_contract_audit_2026-08-29.md`. They
  were preserved untouched.
- The requested EdU audit file did not exist at audit start.

## Scope, protected inputs, and audit method

This was a source, documentation, test-source, and Git-history audit only. The
current merged implementation was compared with the confirmed 2026-08-30 record
and its incorporated 2026-08-23 definitions. Statuses mean:

- **Implemented:** current code satisfies the confirmed item.
- **Partially implemented:** a material subset exists, but required behavior,
  labels, warnings, or packaging are missing.
- **Absent:** no current implementation of the confirmed output exists.
- **Conflicting:** current behavior implements a different contract and would
  violate the confirmed standard.

No `.fcs`, `.wsp`, `.prism`, experimental CSV, packaged example CSV, or other
experimental input was opened. No production or example analysis was performed.
No R or Python code, test, package check, Quarto render, report render, snapshot,
benchmark, or package installation/loading was executed. No Git mutation,
commit, push, pull, fetch, branch creation, PR, or merge was performed.

## Relevant files inspected

### Authority, history, package surface, and documentation

- `AGENTS.md`
- the 2026-08-30 source contract named above
- `../edu_output_contract_approved_2026-08-23.md`
- `../edu_current_implementation_audit_254e090.md`
- Git status, commit identity, recent history, and the file list/statistics for
  EdU implementation commit `c0bd0f9`
- `DESCRIPTION`, `NAMESPACE`, `NEWS.md`, `README.md`
- `config.yml`, `inst/config/config_edu.yml`
- `docs/CONFIGURATION.md`, `docs/MIGRATION.md`, `docs/REPORTS.md`
- `docs/implementation/edu_output_contract_implementation_2026-08-23.md`
- `.gitignore`, `.Rbuildignore`

### R source

- `R/analysis.R`
- `R/appearance.R`
- `R/artifacts.R`
- `R/config.R`
- `R/configurator.R`
- `R/input.R`
- `R/normalize.R`
- `R/phase_quantitation.R`
- `R/pseudocolor_helpers.R`
- `R/results.R`

### Reports and export boundary

- `pseudocolor_plots.qmd`
- `inst/quarto/_report-setup.R`
- `inst/quarto/facs_complete.qmd`
- `inst/quarto/facs_pseudocolor.qmd`
- `inst/quarto/facs_quantitation.qmd`
- `inst/quarto/facs_cell_cycle.qmd`
- `inst/quarto/facs_diagnostics.qmd`
- `python/export_contract.py`
- `python/export_flowjo_populations.py`
- `tools/flowjo-orchestration.R`

The Python export contract was inspected only as a local provenance design
comparison. Its production profile is pH3-specific and the EdU R path does not
consume it; it is not evidence that the EdU report contract is implemented.

### Tests

- `tests/testthat.R`
- `tests/testthat/helper-current-workflow.R`
- `tests/testthat/helper-load.R`
- `tests/testthat/test-analysis.R`
- `tests/testthat/test-appearance-and-artifacts.R`
- `tests/testthat/test-config.R`
- `tests/testthat/test-density-and-palettes.R`
- `tests/testthat/test-edu-output-contract.R`
- `tests/testthat/test-input.R`
- `tests/testthat/test-results-api.R`
- `tests/python/test_export_contract.py`
- `tests/testthat/test-python-boundary.R`

## Executive finding

The merged code substantially implements the earlier, table-centered EdU output
contract: all seven canonical aggregate tables and their acquisition-level forms
exist, approved denominators and intensity populations are represented, rejected
all-cell intensity outputs are absent, and quantitative medians exclude the
display offset.

The newly confirmed standard output layer is not complete. The most serious
conflict is the display offset: current code chooses one pooled power-of-ten
constant, with fallback sources, and applies it to every sample. The confirmed
standard requires a separate value for each sample computed from the exact
eligible negative/background events. The equal-replicate DNA-density output is
entirely absent. Pseudocolor panels exist per manifest row but lack complete
sample/acquisition and QC/suppression labels. Current percentage presentation is
limited to the historical five-gate view, while the exact standard percentage-
figure inventory is not enumerated in the supplied authority. Intensity figures
depend on deprecated aliases, and provenance is an incomplete in-memory list
rather than a shared core plus EdU block with a machine manifest and human
appendix.

## Complete contract matrix — scientific methods and estimands

| Confirmed item | Status | Current evidence | Gap or consequence |
|---|---|---|---|
| Eligible Single Cells are validated events with finite metric-required derived values after only explicit, approved, provenance-recorded exclusions; Unassigned remains eligible | **Partially implemented** | Shared EdU classification requires named fields and records metric-specific eligibility (`R/phase_quantitation.R:590-620`). Input reading fails on missing channels/files (`R/input.R:89-160`). | There is no structured EdU exclusion/suppression ledger or provenance-recorded exclusion reason. Eligibility is present at event/metric level but not complete at sample/output level. |
| Historical five-gate assigned composition and denominator remain unchanged | **Implemented** | `edu_assigned_phase_composition` is constructed at `R/results.R:225-241`; denominator logic is at `R/phase_quantitation.R:734-795`; exact tests are at `tests/testthat/test-edu-output-contract.R:89-120`. | No numerical gap found. The current historical plot is driven through the deprecated alias. |
| Whole-Single-Cells composition: G1, Early/Mid/Late S, Negative S, G2/M over eligible Single Cells; Unassigned QC-only | **Implemented numerically** | `R/results.R:243-252`; `R/phase_quantitation.R:743-790`; tests `tests/testthat/test-edu-output-contract.R:89-120`. | The value is table-only and per-sample QC is not surfaced in reports. Whether a distinct figure is standard is an authority gap, not inferred here. |
| Negative S is computed-negative in `[G1 max, G2/M min)` | **Implemented** | `R/phase_quantitation.R:629-655`; exact-boundary test `tests/testthat/test-edu-output-contract.R:75-87`. | No method gap found. |
| Six-gate assigned composition uses its own assigned-count denominator and sums to 100% | **Implemented numerically** | `R/results.R:254-263`; tests `tests/testthat/test-edu-output-contract.R:104-120`. | The value is table-only; whether a distinct figure is standard is unresolved. |
| Early/Mid/Late regional positivity uses all eligible events in each half-open DNA region, with `NA`/`zero_denominator` for empty regions | **Implemented numerically** | `R/phase_quantitation.R:798-835`; tests `tests/testthat/test-edu-output-contract.R:122-190`. | The value is table-only, and non-OK statuses are not promoted to report warnings. Figure status is unresolved. |
| Overall positivity uses eligible events in `[G1 min, G2/M max)`, including phase-Unassigned denominator events | **Implemented numerically** | `R/phase_quantitation.R:837-877`; tests `tests/testthat/test-edu-output-contract.R:122-168`. | Phase-Unassigned and status information is table-only. Figure status is unresolved. |
| Positive-cell regional intensity is median background-subtracted EdU among computed-positive eligible Early/Mid/Late-S cells | **Implemented numerically** | `R/results.R:284-294`; `R/phase_quantitation.R:879-914`; identity tests `tests/testthat/test-edu-output-contract.R:192-230`. | Figure exists only through deprecated `phase_medians`; canonical linkage and status annotations are partial. |
| Whole computed-positive-population intensity is median background-subtracted EdU among all computed-positive eligible cells | **Implemented numerically** | `R/results.R:296-305`; `R/phase_quantitation.R:915-933`; tests `tests/testthat/test-edu-output-contract.R:192-230`. | Figure exists only through deprecated `whole_medians`; canonical linkage and status annotations are partial. |
| Do not calculate standard all-cell regional or whole-Single-Cells intensity | **Implemented** | Absence is asserted at `tests/testthat/test-edu-output-contract.R:50-67`; automatic names are limited at `R/results.R:872-881`. | No conflicting standard output found. |
| Automatically calculate seven canonical tables independent of legacy display switches; retain acquisition-level values | **Implemented** | High-level analysis always calls all quantitation types (`R/analysis.R:216-235`); canonical tables are built at `R/results.R:225-306`; switch invariance is tested at `tests/testthat/test-edu-output-contract.R:274-295`. | No calculation-availability gap found. Figure visibility is assessed separately because the supplied records do not define a complete standard figure mapping. |
| Stable canonical names, explicit metadata, schema version, deprecated aliases only for compatibility | **Partially implemented** | Canonical table names are at `R/results.R:872-881`; output schema and deprecation warning at `R/results.R:345-397`; CSV opt-in rules at `R/results.R:884-907`. | Provenance still points to the 2026-08-23 decision and old offset formula. Figures still consume aliases. Metadata lacks the new density, QC, and per-sample offset contract. |
| One DNA-versus-EdU pseudocolor plot for every biological replicate/sample; no representative sample | **Partially implemented** | All EdU manifest rows display by default (`R/results.R:10-17`); one panel is built per row (`R/results.R:742-749`) and individual plots are keyed by prefix (`R/results.R:804-807`). Per-panel KDE uses that panel's displayed events without downsampling (`R/pseudocolor_helpers.R:1012-1035`). | The assembled grid is compatible with the layout-open contract, but each plot title is condition only (`R/pseudocolor_helpers.R:1098-1110`) and lacks complete sample/replicate and QC/suppression labeling. No separate-file-per-sample requirement is inferred. |
| Every pseudocolor panel identifies condition, biological replicate/sample identity, axes, and applicable QC/suppression | **Partially implemented / conflicting for technical acquisitions** | Condition is the panel title and axes are labeled (`R/pseudocolor_helpers.R:1098-1103`); biological replicate is an external row label (`R/results.R:787-801`). | Standalone plots omit replicate/sample identity and QC. With technical acquisitions, `n_rows` is based on `model_group` while row labels use unique biological-replicate names (`R/results.R:787-790`), which can mislabel or fail (`R/pseudocolor_helpers.R:1165-1169`, `1300-1307`). |
| Failed or unavailable sample is never silently replaced | **Implemented for the current standard input/computed-field path; suppression presentation absent** | Required files/channels fail explicitly in the current input path (`R/input.R:89-160`). The canonical classifier fails when required computed fields are missing (`R/phase_quantitation.R:590-605`). | A lower helper can fall back from missing `edu_computed_positive` to exported `edu_positive` events (`R/phase_quantitation.R:219-229`), but the standard automatic path fails earlier. Treat it as a latent hazard that new figure code must not reuse. No failed/suppressed sample placeholder exists. |
| Per-sample display offset equals `median(raw eligible negative) - median(background-corrected eligible negative)` using the same negative/background events | **Conflicting** | Background fitting constructs per-acquisition `negative_events` (`R/pseudocolor_helpers.R:727-743`) but does not retain the required two medians, membership proof, or resolved offset in the returned model (`R/pseudocolor_helpers.R:820-831`). | Current resolver instead pools other values across samples and returns one rounded scalar (`R/results.R:654-682`). This is a blocking scientific-method conflict. |
| Do not use power-of-ten rounding, a fixed target, or one global add-back | **Conflicting** | `resolve_background_subtracted_offset()` uses `10^round(log10(center))` and fallback sources (`R/results.R:654-682`); config permits a fixed numeric offset (`R/config.R:251-259`); tests explicitly expect both (`tests/testthat/test-results-api.R:77-105`). | Standard behavior directly violates the confirmed 2026-08-30 convention. Any manual/fixed legacy route must not remain the standard path. |
| Display offset affects plotted values only, never positivity, corrected medians, reference normalization, thresholds, statistics, or quantitative exports | **Implemented numerically; operationally coupled** | Classifier separates `quantitative_signal` from `display_signal` (`R/phase_quantitation.R:607-620`, `656-690`); intensities use the quantitative coordinate (`R/phase_quantitation.R:879-933`); invariance is tested (`tests/testthat/test-edu-output-contract.R:234-250`). | Canonical quantitation resolves the display offset before calculating tables (`R/results.R:136-165`), so a display-only offset failure can abort quantitative results. Per-sample replacement must preserve numerical isolation and decouple failure domains. |
| Standard report uses the correct display coordinate and limits | **Conflicting in the bundle/report path** | Figure bundle resolves appearance before applying the display transform (`R/artifacts.R:63-72`). Appearance derives automatic limits from `analysis_y_limits()` (`R/appearance.R:88-90`), which reads pre-display `target_norm` (`R/results.R:20-40`). | With automatic limits, transformed `target_bgsub + offset` panels can be cropped using limits derived from the legacy divided coordinate. |
| Condition DNA density is the pointwise arithmetic mean of independently estimated, unit-area biological-replicate densities on identical deterministic coordinates and bandwidth | **Absent** | No 1D EdU DNA-density table, estimator, plot, config, report section, export, or test exists. The only density helper is a per-panel 2D event-color KDE with sample-specific bandwidth (`R/pseudocolor_helpers.R:85-112`, `1023-1035`). | Blocking missing standard scientific output. The 2D pseudocolor KDE cannot be relabeled or reused as the required 1D estimand. |
| DNA density uses all eligible events, no downsampling, no cross-biological-replicate pooling, and no event-count weighting | **Absent** | No condition-density implementation exists. Standard pseudocolor does not downsample. The only explicit downsampling is the optional gate-assignment diagnostic (`R/phase_quantitation.R:1104-1131`), seeded/restored at `R/results.R:546-569`. | New code must not call the diagnostic sampler or concatenate biological-replicate events. Absence must not be mistaken for compliance. |
| Density grid/bandwidth are deterministic, shared, explicit, reproducible, and invariant to file/input order | **Absent for density** | Current manifest is explicit/config-backed and preserves configured order (`R/input.R:3-16`); replicates must contain the same conditions in the same order (`R/pseudocolor_helpers.R:249-260`). No directory/file-order inference was found in the current EdU path. | No density implementation or permutation-invariance test exists. Exact grid and bandwidth are Code Architect mechanics, not new scientific decisions. |
| Optional faint replicate density curves use the same grid/bandwidth and do not alter the mean | **Absent, optional** | No density plot exists. | Safe initial implementation may keep this option off. Enabling it is presentation-only if it reads the same immutable replicate-density table. |
| No inferential interval, new statistic, threshold, gate, normalization, or QC rule is introduced | **Implemented for current change state** | No new EdU density inference exists; current generic plots offer descriptive SD/SEM/none (`R/phase_quantitation.R:288-306`, `R/config.R:77-84`). | The planned density figure must not reuse those error bars. This audit does not authorize any scientific-method change. |

## Complete contract matrix — figures, reports, exports, and provenance

| Confirmed item | Status | Current evidence | Gap or consequence |
|---|---|---|---|
| Common assay-agnostic provenance core with schema and workflow identities | **Absent** | Base provenance contains package/version, R version, config paths, and input paths only (`R/analysis.R:74-87`). | Missing schema identity/version, workflow identity distinct from package, run ID/time, source integrity, acquisition/instrument metadata, specification identities, status ledger, warnings, and artifact identities. |
| Core identifies samples, conditions, biological replicates, and included/excluded/failed/suppressed outputs with reasons | **Partially implemented in separate objects; absent as provenance** | Sample identities exist in `sample_manifest`; metric statuses exist in canonical tables (`R/phase_quantitation.R:734-939`). | They are not assembled into the common provenance core. Non-OK technical statuses can become `partial` while finite values are averaged (`R/phase_quantitation.R:942-989`) but are not promoted to the warning ledger. |
| Core records source identifiers and checksums where available, processing/gating/correction/normalization/QC spec identities, all interpretation warnings, and generated artifact identities | **Absent** | Current provenance stores absolute paths without hashes (`R/analysis.R:80-83`). Output writer returns paths only (`R/results.R:934-1014`). | Reports can expose machine-specific paths; no portable input identity, artifact checksum/registry, or cross-record reconciliation exists. |
| Clearly separated EdU assay block contains background population, correction, positivity/threshold, medians/reference normalization, density mechanics/rule, per-sample offsets, and EdU QC | **Absent / conflicting** | Current `edu_output_contract` block records old decision IDs, transforms, one scalar offset, and reference-normalization status (`R/results.R:347-388`). | It is not a complete EdU block, points to the old contract, documents the obsolete global offset, and has no density or QC content. |
| Each standard report shows a concise provenance summary and all provenance warnings | **Absent** | Complete report dumps raw input report, warning vector, provenance list, and appearance together (`inst/quarto/facs_complete.qmd:109-120`). Pseudocolor and diagnostics reports similarly dump raw objects. Quantitation and cell-cycle reports have no provenance/warning section. | There is no shared concise summary component, and not every standard report shows warnings. Templates also hide runtime warnings globally with `warning: false`. |
| Complete machine-readable provenance manifest is included in every report/export package | **Absent** | `save_facs_results()` writes PDF, PNG, RDS, and optional CSVs only (`R/results.R:934-1014`). | No manifest serializer, schema validator, checksum, package member registry, or returned manifest path exists. |
| Complete human-readable provenance appendix is included in every report/export package | **Absent** | Reports print raw lists; report helper saves analysis and bundle RDS only (`inst/quarto/_report-setup.R:23-33`). | No complete, consistently rendered appendix exists. |
| Concise summary, full manifest, and appendix are consistent views of one source | **Absent** | No shared provenance schema or renderer exists. | Independent ad hoc printing would risk drift; all views must be generated from one validated object. |
| Generated quantitative exports and report artifacts are identified and reconciled | **Absent** | CSV names are validated before writing (`R/results.R:884-907`, `951-979`), but the writer serially emits files and only returns paths (`R/results.R:988-1014`). | No artifact identity/checksum is inserted into provenance. A mid-write failure can leave a partial package. |
| Existing outputs refuse accidental overwrite and avoid external transmission | **Implemented as a safety control** | Explicit overwrite guards exist (`R/artifacts.R:16-24`; `R/results.R:969-975`). No network, telemetry, upload, credential collection, or remote call was found in the EdU analysis/report/export path. | Publication is not transactional/atomic; this is an artifact-integrity gap, not evidence of external transmission. |
| Full contract has focused deterministic, input-free tests | **Partially implemented** | Existing `SYNTHETIC` tests cover the seven tables, boundaries, offset isolation, technical aggregation, and CSV names (`tests/testthat/test-edu-output-contract.R:1-399`). | Tests lock the obsolete global/fixed offset and do not cover per-sample panels/QC, common density, order invariance, shared provenance, manifest, appendix, warning completeness, or package reconciliation. |

### Requested assessment: percentage and intensity figure sets

This assessment is intentionally outside the confirmed-item status matrix where
the authority is incomplete. The current percentage figure set contains only the
historical five-gate assigned-composition plot (`R/results.R:501-519`), while the
other approved percentage values are table-only in the complete report
(`inst/quarto/facs_complete.qmd:70-96`). The supplied confirmed records approve
the tables and calculation independence but do not identify which additional
percentage figures are standard, their geometry/error display, or whether their
visibility must ignore legacy display switches. Those points are unresolved and
must not be inferred.

The established EdU intensity views are present: regional computed-positive and
whole computed-positive plots correctly state that the display offset is
excluded (`R/results.R:456-499`). They are only partially canonicalized because
they consume deprecated aliases, use generic plot names, and do not display
partial/insufficient status. Relinking them to numerically identical canonical
tables is a presentation/export change; changing their estimand or error display
is not authorized by this audit.

### Scientific-method gaps versus presentation/export gaps

| Gap class | Items |
|---|---|
| Scientific-method/output-estimand gaps | Per-sample negative-median restoration offset; retained proof of the exact negative/background population; equal-weight unit-area biological-replicate DNA density; confirmed technical-acquisition rule for that density; confirmed QC/stop rules governing inclusion and suppression. |
| Report/export/presentation-only gaps | Complete panel identity and QC labels; display-coordinate y-limit bug; canonical linkage/status annotation for existing intensity figures; concise provenance summary; full human appendix; machine manifest; artifact registry/reconciliation; layout and optional faint replicate curves. |
| Authority gap before classifying additional figures | Exact standard percentage/intensity figure inventory, allowed error display, legacy-switch visibility policy, and exact set of repository report templates considered standard. |

## Current behavior that would violate or endanger the confirmed contract

### Blocking conflicts

1. **Global rounded/fixed display offset.** One scalar is selected across all
   samples from positive G1 anchors, then baseline or raw-target fallbacks, and
   rounded to a power of ten (`R/results.R:654-682`). The same value is applied
   to every sample (`R/results.R:690-700`). This directly violates the confirmed
   per-sample raw-negative-median restoration formula.
2. **No equal-replicate DNA-density estimand.** There is no current result to
   validate. Adding a pooled KDE, event-count-weighted mean, downsampled KDE, or
   average on sample-specific coordinates/bandwidth would violate the contract.
3. **Missing per-sample QC/suppression state.** A standard panel cannot display
   required warning/suppression information, and a failed sample cannot be
   represented as a traceable suppressed output.
4. **Technical-acquisition identity/layout conflict.** Current panel row counts
   and labels use different keys. More importantly, the confirmed density
   estimand is biological-replicate-level while the repository supports multiple
   technical acquisitions per biological replicate; silently pooling events or
   averaging acquisition densities would select an unapproved estimand.

### Important latent or presentation risks

- `quant_source_data()` can substitute the exported legacy EdU-positive table
  when computed positivity is missing (`R/phase_quantitation.R:219-229`). The
  standard classifier currently fails first, but new report code must not rely
  on this fallback.
- Current standard pseudocolor plots do not downsample, and current sample
  mapping is explicit/config-backed. The only event downsampling found is in an
  optional gate-assignment diagnostic. That diagnostic must remain visibly
  labeled and must never feed density estimation or canonical quantitation.
- Automatic report y limits may be calculated on the legacy divided coordinate
  and applied after the background-subtracted display transform. This can crop
  the standard plot without changing the underlying quantitative data.
- Metric-level `zero_denominator`, `insufficient_events`, `partial`, and `mixed`
  states are stranded in tables rather than consistently represented as report
  warnings and provenance status/reason records.
- No current numerical EdU table appears to depend on filesystem order. The
  manifest is built from explicit configuration order. The new density result
  still needs permutation tests because deterministic order invariance is part
  of its contract.
- Display-offset isolation is numerically protected today, but the display
  resolver is called during canonical quantitation. Display failure should not
  silently change values or create a substitute; the implementation should
  make the coupling and failure policy explicit.

## Unresolved items and decision boundary

### Blocking scientific or authoritative-record gaps

1. **QC and stop-condition source.** The 2026-08-30 record requires the already
   confirmed QC/stop contract but does not enumerate its sample/output states,
   suppression rules, or reason vocabulary. The current repository has input
   failures and metric statuses but no complete EdU sample/output suppression
   model. Locate the authoritative confirmed record or obtain owner confirmation
   before defining panel suppression, density inclusion, or omission behavior.
2. **Technical acquisitions in biological-replicate density.** For a biological
   sample represented by multiple technical acquisitions, event pooling,
   acquisition-level unit-area density averaging, and another consolidation rule
   are different estimands. The 2026-08-30 record prohibits pooling across
   biological replicates but does not resolve within-biological-replicate
   technical acquisitions. The antecedent explicitly left technical-acquisition
   aggregation outside its scope. This requires a confirmed rule.
3. **Canonical percentage/intensity figure and standard-report inventory.** The incorporated record
   confirms seven standard tables but does not explicitly map every percentage
   table to a required figure or define each figure's allowed error display. The
   user-visible figure set must be recovered from the referenced prior confirmed
   design or confirmed by the owner; do not infer new scientific figures. The
   supplied records also do not identify which repository report templates are
   “standard,” so the per-standard-report provenance requirement cannot be
   bounded to an exact template set until that inventory is confirmed.

The exact retained eligible negative/background population is not a new choice
if implementation can carry forward the same `negative_events` used by the
current event-specific background fit. The implementation must preserve and
prove that identity rather than reconstructing a different computed-negative
population downstream. If that identity cannot be retained unambiguously, stop
and return the population definition to the scientific owner.

### Code Architect mechanics explicitly left open by the confirmed record

- exact deterministic DNA coordinate grid;
- exact deterministic bandwidth rule;
- provenance field names, serialization format, schema versioning, controlled
  vocabularies, and package layout;
- visual layout for all sample panels and warnings; and
- whether optional faint replicate density curves are shown.

These are not scientific blockers when they are explicit, reproducible,
order-invariant, shared across compared replicates, and do not change the
experimental unit, denominator, unit-area normalization, equal-replicate mean,
analytical values, or scientific meaning. Optional replicate curves may remain
off initially.

## Bounded implementation slices in dependency order

All tests proposed below must use only deterministic, unmistakably `SYNTHETIC`
test fixtures stored in test-only locations or constructed in memory. No
experimental inputs are needed or authorized.

Unless a slice states an additional prerequisite, every Slice 1-5 verification
command uses this working directory:

```text
/Users/sansamc/OMRF Dropbox/Chris Sansam/Synchd/Organized Dropbox/Projects and Data/Projects/Cowork_Projects/2026_Analyze_FACS_Data/facs_pseudocolor_workflow
```

Common prerequisites are R 4.2 or newer, `devtools`, `testthat`, and all declared
package dependencies already installed. No production or example experimental
input is permitted. Expected persistent/generated repository files are none;
temporary synthetic files may exist only beneath the test framework's temporary
directory and must be removed at test completion. These working-directory,
prerequisite, and generated-file statements are part of each Slice 1-5 local
verification instruction below.

### Slice 0 — recover the missing confirmed decision inputs

**Purpose:** close only the authoritative gaps that would otherwise force an
unapproved scientific choice.

**Exact scope:** locate or obtain confirmation for the EdU QC/stop-condition
state-and-reason contract, the standard phase-percentage/intensity figure
inventory and standard-report template set, and the treatment of multiple
technical acquisitions in one biological-replicate DNA density. Record decisions
without changing code.

**Likely files affected:** a new owner-confirmed decision record outside or
inside `docs/implementation/` as directed by the scientific owner; no production
source, config, or test file.

**Input/output contract:** input is explicit scientific-owner confirmation;
output is a versioned decision record with no inferred defaults.

**Prerequisite:** the exact authoritative records or explicit scientific-owner
confirmation for every listed blocker.

**Test intent:** none; review the record for internal consistency, owner
authority, and exact scope.

**Local verification command:** none. Documentation review only.

**Expected success and generated file:** all blockers are explicitly resolved
without inferred defaults, and exactly one versioned decision record is created
at the location authorized by the owner.

**Stop conditions:** ambiguity, contradiction, missing owner authority, or any
production source, test, configuration, experimental input, generated report, or
Git-state change.

**Non-goals:** no code, no method selection by the architect, no reinterpretation
of current defaults as approval.

### Slice 1 — every-sample pseudocolor and per-sample display offset

**Purpose:** correct the highest-priority visible EdU output and the current
display-method conflict.

**Exact scope:** retain immutable evidence for the exact negative/background
events used in each sample's background fit; calculate and validate the two
medians and their difference per sample; apply the keyed value only to display
copies; preserve quantitative fields unchanged; produce one keyed plot per
sample with condition, biological replicate, sample/acquisition identity, axes,
and confirmed QC/suppression status; make automatic limits use the actual display
coordinate; remove the global power-of-ten/fallback/fixed-target route from the
standard path. Fail closed on missing identity, population, fields, or status.

**Likely files affected:** `R/pseudocolor_helpers.R`, `R/analysis.R`,
`R/results.R`, `R/appearance.R`, `R/artifacts.R`, `R/config.R`, `config.yml`,
`inst/config/config_edu.yml`, `pseudocolor_plots.qmd`,
`inst/quarto/facs_pseudocolor.qmd`, `inst/quarto/facs_complete.qmd`, a new
`tests/testthat/test-edu-pseudocolor-output-contract.R`, and focused revisions to
`tests/testthat/test-results-api.R` and
`tests/testthat/test-edu-output-contract.R`. Documentation/Rd updates should be
made only through their source.

**Input contract:** quantified EdU analysis with explicit manifest identity,
exact retained eligible-negative membership or immutable membership proof, raw
and corrected EdU values, and the confirmed sample/output QC ledger.

**Output contract:** per-sample offset ledger; one individually identifiable
editable plot per sample; assembled report layout; no changed positivity,
threshold, median, ratio, statistic, or analytical CSV value.

**Synthetic/input-free test intent:** use at least two biological replicates and
two conditions with deliberately different negative medians; hand-check the
formula for every sample; permute event and manifest order; assert complete panel
identity and status labels; assert missing masks/samples/statuses fail; assert no
sample substitution; assert only display values change; assert thresholds,
positivity, canonical tables, reference ratios, and quantitative CSV columns are
identical before/after offset resolution.

**Local-only verification command:** using the common working directory and
prerequisites above, run:

```bash
R -q -e 'devtools::test(filter = "edu-pseudocolor-output-contract", stop_on_failure = TRUE)'
```

Expected success: `FAIL 0`, `WARN 0`, no skipped contract test, and no persistent
file created in the repository. Temporary test files, if any, exist only beneath
the test temporary directory and are removed afterward. Stop on any changed
analytical value, missing panel, population/status fallback, warning, skip,
failure, or unexpected file.

**Non-goals:** no background-fit, positivity, gate, threshold, DNA normalization,
reference normalization, statistic, or new QC-rule change; no density summary.

### Slice 2 — confirmed phase-percentage and EdU-intensity figure set

**Purpose:** after Slice 0 confirms the exact inventory, expose only those
canonical quantitative outputs as standard, traceable figures instead of only
tables or deprecated aliases.

**Exact scope:** after Slice 0 resolves the authoritative inventory, map only
the confirmed canonical tables/value fields to named editable plots; show the
approved source population, denominator, transformation, status, and warning
labels; make figure visibility follow the exact Slice 0 display-policy decision;
keep legacy views explicit and nonstandard. This slice does not assume that
legacy switches must be ignored.

**Likely files affected:** `R/results.R`, preferably a new focused
`R/edu_figures.R`, `R/artifacts.R`, `R/appearance.R`,
`pseudocolor_plots.qmd`, `inst/quarto/facs_complete.qmd`,
`inst/quarto/facs_quantitation.qmd`, `inst/quarto/facs_cell_cycle.qmd`, a new
`tests/testthat/test-edu-standard-figures.R`, and report/API documentation.

**Input contract:** the seven existing canonical aggregate and acquisition
tables plus the resolved figure mapping and status/warning ledger.

**Output contract:** named plots whose plot-ready data exactly match the named
canonical source table and value column. Existing scientific values remain
unchanged.

**Synthetic/input-free test intent:** assert every required and only required
figure exists; trace each layer to its canonical table/value; assert exact
denominator/source labels; assert display-switch behavior exactly matches the
Slice 0 decision and cannot alter plotted values; assert no rejected all-cell
intensity figure; assert offset is absent from quantitative values; assert
partial/insufficient statuses are visible.

**Local-only verification command:** using the common working directory and
prerequisites above, run:

```bash
R -q -e 'devtools::test(filter = "edu-standard-figures", stop_on_failure = TRUE)'
```

Expected success: `FAIL 0`, `WARN 0`, no skipped contract test, and no repository
outputs. Stop on an unsupported plot/statistic, missing label/status, mismatch
with a canonical table, changed scientific value, warning, skip, or unexpected
file.

**Non-goals:** no new figure not confirmed in Slice 0; no gate, threshold,
denominator, median, reference-normalization, error-bar, or statistical change.

### Slice 3 — equal-biological-replicate DNA-density result and figure

**Purpose:** implement the absent condition-level DNA-density standard while
preserving biological replicate as the experimental unit.

**Exact scope:** select and version an exact deterministic common coordinate grid
and shared bandwidth rule; for every biological replicate and condition use all
eligible events, estimate independently on the same coordinates/bandwidth,
normalize numerically to unit area, and then compute the pointwise arithmetic
mean of replicate densities; record event counts, replicate counts, status,
resolved mechanics, and warnings; draw the standard condition comparison.
Optional faint replicate curves remain off initially or read the same immutable
replicate-density table without changing the mean.

**Likely files affected:** preferably a new `R/edu_density.R`, `R/results.R`,
`R/artifacts.R`, `R/appearance.R`, `R/config.R` only if versioned mechanics need
configuration, `inst/quarto/facs_complete.qmd`,
`inst/quarto/facs_quantitation.qmd`, a new
`tests/testthat/test-edu-dna-density-contract.R`, documentation, and
`NAMESPACE`/Rd only if a public API is approved.

**Input contract:** explicit eligible-event membership keyed by condition and
biological replicate/sample, normalized DNA coordinate, confirmed QC ledger,
and the Slice 0 technical-acquisition rule.

**Output contract:** replicate-density and condition-mean tables with common
coordinate, density, unit-area check, event count, replicate count, grid and
bandwidth identifiers/resolved values, status, and warnings; a plot labeled with
condition, DNA axis, unit-area/equal-replicate statement, included replicate
count, and applicable exclusions/QC.

**Synthetic/input-free test intent:** use unequal event counts so the expected
condition mean is exactly `(d1 + d2) / 2` and differs from pooled/event-weighted
KDE; assert each replicate integrates to one within the documented numerical
tolerance; assert identical coordinates/bandwidth; assert all eligible events
are used; prohibit sampling/downsampling; permute events, replicates, and input
order; assert optional overlays do not alter the mean; assert failed/missing
replicates follow the confirmed stop/suppression rule with explicit reasons.

**Local-only verification command:** using the common working directory and
prerequisites above, run:

```bash
R -q -e 'devtools::test(filter = "edu-dna-density-contract", stop_on_failure = TRUE)'
```

Expected success: `FAIL 0`, `WARN 0`, no skipped contract test, exact arithmetic
and permutation assertions, unit-area checks within the documented tolerance,
and no repository output. Stop on pooling, event-count weighting, sampling,
non-unit area, grid/bandwidth mismatch, omitted replicate, order dependence,
warning, skip, failure, or unexpected file.

**Non-goals:** no cross-replicate event pooling, event-count weighting,
downsampling, stochasticity, inferential band/error display, new normalization,
gate, threshold, statistic, QC rule, or advanced output.

### Slice 4 — common provenance core and EdU assay block

**Purpose:** establish one validated provenance source for all reports and
exports.

**Exact scope:** define and version an assay-agnostic core; add run identity/time,
workflow identity/version, explicit source identifiers and integrity identifiers
where available, sample/condition/replicate identities, acquisition/instrument
metadata where available, processing/gating/correction/normalization/QC spec
identities, full included/excluded/failed/suppressed ledger and reasons, all
interpretation/reproducibility warnings, and an artifact registry. Add a clearly
separated EdU block with the retained negative/background definition,
event-specific correction, positivity/resolved thresholds, corrected-median and
reference-normalization specifications, density mechanics and equal-replicate
rule, per-sample offset formula/values, and EdU QC decisions/warnings. Preserve
existing assay-specific pH3 evidence rather than flattening or discarding it.

**Likely files affected:** preferably a new `R/provenance.R`, `R/analysis.R`,
`R/input.R`, `R/results.R`, `R/artifacts.R`, and new
`tests/testthat/test-provenance-contract.R` plus focused assay tests.

**Input contract:** immutable analysis, config, sample/input identities, models,
resolved mechanics, QC decisions, warnings, and planned artifact identities.

**Output contract:** a validated common core and named assay block. Unavailable
metadata is explicitly unavailable; no value is invented or borrowed from a
nearby file. No network or external transmission.

**Synthetic/input-free test intent:** validate required fields/types and
core/assay separation; test permutation-stable sample mapping; calculate hashes
only for explicit synthetic temporary files; ensure unavailable acquisition
metadata is not fabricated; preserve every warning and reason; cross-check EdU
offset/density values against computation records; fail on tampering or missing
required fields; verify POI/pH3 can share the core without semantic regression.

**Local-only verification command:** using the common working directory and
prerequisites above, run:

```bash
R -q -e 'devtools::test(filter = "provenance-contract", stop_on_failure = TRUE)'
```

Expected success: `FAIL 0`, `WARN 0`, no skipped contract test, and no persistent
repository output. Synthetic checksum fixtures may exist only beneath the test
framework's temporary directory during the test and are removed afterward. Stop
on missing or stale fields, invented/defaulted metadata, warning loss, assay
regression, warning, skip, failure, or unexpected file.

**Non-goals:** no analytical, gate, threshold, mapping, exclusion, normalization,
statistical, or biological-claim change; no telemetry, upload, or network call;
no manifest or report rendering yet.

### Slice 5 — standard report/export package, manifest, and appendix

**Purpose:** deliver a complete, consistent user-facing standard package.

**Exact scope:** generate a concise provenance summary and all provenance
warnings from the validated object; generate a complete human-readable appendix;
serialize a complete machine-readable manifest; register and reconcile every
quantitative export and report artifact; assemble every-sample pseudocolor,
confirmed canonical percentage/intensity figures, equal-replicate DNA density,
canonical tables, warnings, manifest, and appendix. Preserve explicit no-clobber
behavior and use staged/finalized publication so a failure does not silently
leave a plausible partial package. Define exact file names, serialization, and
layout as versioned Code Architect mechanics.

**Likely files affected:** `inst/quarto/_report-setup.R`, exactly the report
templates confirmed as standard by Slice 0, `R/results.R`, `R/artifacts.R`, the
new provenance/figure/density source files, a new
`tests/testthat/test-edu-report-export-contract.R`, a test-only renderer under
`tests/fixtures/SYNTHETIC_edu_report_contract/`, and report/config/API
documentation. No unconfirmed report template is silently designated standard.

**Input contract:** complete analyzed object, standard figure bundle, and
validated provenance object with resolved output paths and statuses.

**Output contract:** one requested report/export package containing only
confirmed standard reports, figures, tables, the complete machine manifest, and
the complete human appendix; returned artifact registry exactly reconciles with
the published files and their identities.

**Synthetic/input-free test intent:** static template contract checks and
temporary-directory export tests; exact consistency among summary, manifest,
appendix, and analysis; all warnings visible; every artifact registered;
machine manifest parse/round-trip; no-clobber and interrupted-publication
behavior; no partial package on preflight failure; no network call.

**Local-only verification command:** using the common working directory and
prerequisites above, run the focused export test:

```bash
R -q -e 'devtools::test(filter = "edu-report-export-contract", stop_on_failure = TRUE)'
```

Expected success: `FAIL 0`, `WARN 0`, no skipped contract test, and only declared
files inside test temporary directories. Stop on a missing warning/artifact,
summary/manifest/appendix mismatch, undeclared extra file, partial publication,
warning, skip, failure, or repository artifact.

Then perform an actual input-free render. Additional prerequisite: Quarto and the
declared report-rendering packages are already installed. Slice 5 must provide a
test-only renderer at the path below that creates its analysis only from a
clearly labeled deterministic `SYNTHETIC` fixture and renders exactly the
Slice-0-confirmed standard report set:

```bash
VERIFY_DIR="$(mktemp -d /private/tmp/facspseudocolor-edu-report.XXXXXX)"
Rscript tests/fixtures/SYNTHETIC_edu_report_contract/render.R "$VERIFY_DIR"
```

Expected generated files include
`$VERIFY_DIR/SYNTHETIC-standard-report-package/provenance-manifest.json`,
`$VERIFY_DIR/SYNTHETIC-standard-report-package/provenance-appendix.md`, and the
confirmed standard report, figure, and table files. Because their exact names
are Code Architect packaging mechanics, Slice 5 must freeze and enumerate every
one in its implementation record and in the machine manifest before this command
is considered verification-ready. No undeclared file is allowed. Expected
success is a zero exit status, complete manifest-to-filesystem reconciliation,
all confirmed HTML report files present, and no render warning. Stop on any
warning, render failure, missing/extra/unregistered file, manifest/appendix/report
mismatch, use of non-synthetic input, or repository output.

**Non-goals:** no new scientific export, figure, statistic, QC method, telemetry,
upload, or production-data render; do not edit generated artifacts as source.

## Final local-only integration verification after implementation

Working directory:

```text
/Users/sansamc/OMRF Dropbox/Chris Sansam/Synchd/Organized Dropbox/Projects and Data/Projects/Cowork_Projects/2026_Analyze_FACS_Data/facs_pseudocolor_workflow
```

Prerequisites:

- R 4.2 or newer;
- `devtools`, `testthat`, and declared package dependencies already installed;
- a shell supporting `mktemp`;
- no production or experimental input is required for the synthetic contract
  tests.

Run the full suite:

```bash
R -q -e 'devtools::test(stop_on_failure = TRUE)'
```

Build and check only in a temporary directory:

```bash
PROJECT_DIR='/Users/sansamc/OMRF Dropbox/Chris Sansam/Synchd/Organized Dropbox/Projects and Data/Projects/Cowork_Projects/2026_Analyze_FACS_Data/facs_pseudocolor_workflow'
VERIFY_DIR="$(mktemp -d /private/tmp/facspseudocolor-edu-output.XXXXXX)"
(
  cd "$VERIFY_DIR" &&
  R CMD build "$PROJECT_DIR" &&
  R CMD check --no-manual facspseudocolor_0.1.0.9000.tar.gz
)
```

Expected success conditions:

- full tests report `FAIL 0` and `WARN 0` with no skipped contract test;
- the package tarball and `facspseudocolor.Rcheck/` exist only beneath
  `$VERIFY_DIR`;
- `R CMD check` reports `Status: OK`, zero errors, zero warnings, and no
  unexpected notes;
- no analysis, report, figure, CSV, RDS, manifest, appendix, cache, package, or
  temporary artifact appears in the repository.

Stop rather than weaken validation or change a scientific method if any test
fails or skips, any unexpected warning/note occurs, any confirmed artifact is
missing, any analytical value changes unexpectedly, any population or sample is
substituted/omitted, or any repository file is unexpectedly generated. The local
operator should report only the command, `PASS` or `FAIL`, failed test names, a
short relevant error excerpt, and unexpected generated files—not a large
successful log.

## Independent review findings

### Scientific-integrity review

The independent reviewer inspected the actual source and confirmed these
blocking findings:

- current global rounded/fixed offset conflicts with the per-sample formula;
- negative-event membership used for fitting is not retained for exact offset
  proof;
- the equal-weight unit-area DNA-density estimand is absent;
- technical-acquisition handling for biological-replicate density is unresolved;
- sample-level QC/suppression rules needed for panels and density are missing;
- technical-acquisition panel labels/layout can be ambiguous or fail.

The reviewer also confirmed that current standard pseudocolor is per sample and
does not downsample; sample mappings are explicit rather than inferred from
filesystem order; the seven canonical tables are numerically implemented; and
the display offset is excluded from current quantitative intensity values. The
latent computed-positive-to-exported-positive fallback remains an important
safeguard issue for future figure code.

### Artifact/security review

The independent reviewer confirmed that the shared provenance schema, complete
EdU block, report summary/warning component, machine manifest, human appendix,
artifact registry, and contract tests are absent. Current reports print raw
objects, not consistent concise/full views. Output writes are guarded against
accidental overwrite but are sequential rather than staged and reconciled, so a
failure can leave a partial package.

No network, telemetry, upload, credential collection, broad destructive file
operation, or external experimental-data transmission was found in the EdU
analysis/report/export path. Current generated/source hygiene and explicit
overwrite guards are positive controls.

### Final-document review and resolution

The reviewers then inspected this actual audit file. Their important findings
were resolved as follows:

- removed the unsupported implication that the contract requires one output file
  per sample; the real requirement is one separately identifiable, fully labeled
  plot/panel per sample, and the assembled layout remains an architect choice;
- moved the additional percentage/intensity figure inventory and legacy-switch
  policy out of the confirmed-item matrix and into an explicit authority gap;
- classified the legacy exported-positive fallback as a latent helper hazard,
  not active substitution in the current standard automatic path;
- added the exact standard-report template inventory to Slice 0 rather than
  designating templates by inference;
- removed an unsupported event-free figure-bundle requirement;
- added shared working-directory/prerequisite/generated-file rules to every
  implementation-slice verification, clarified temporary-file containment, and
  added an actual synthetic report-render verification requirement; and
- required Slice 5 to freeze and enumerate its architect-selected package member
  names before local verification is considered ready.

After those resolutions, the scientific reviewer reported no other blocking or
important scientific-method error, and the artifact reviewer reported no other
blocking completeness or security finding.

## Self-review against code and contract

The final document was checked against:

- every required 2026-08-30 contract section: every-sample pseudocolor,
  equal-replicate density, shared provenance, per-sample offset, standard versus
  optional output, architect mechanics, and explicit deferrals;
- every incorporated 2026-08-23 canonical table/population/denominator and
  compatibility rule;
- current source implementations and surrounding code rather than prior
  implementation summaries alone;
- current report and export routes, including the root report and installed
  templates;
- existing focused tests, including tests that presently enforce the obsolete
  offset behavior;
- pooling, downsampling, input/file order, display-versus-quantitative isolation,
  population fallback, warning/suppression, and provenance concerns; and
- both independent reviewers' read-only findings.

No unsupported threshold, gate, normalization, statistic, QC criterion, export,
or biological claim was added. Exact grid, bandwidth, provenance mechanics, and
layout are correctly classified as architect decisions. Missing QC, figure-set,
and technical-acquisition estimand decisions are recorded as blockers rather
than guessed.

## Handback summary for the Code Architect

- **Baseline:** `ead9e13e18c8f6243e9622161e0c9b475a42da0b` in the canonical
  `facs_pseudocolor_workflow/` repository; local `main`/`origin/main` match.
- **Implemented foundation:** seven canonical EdU aggregate tables and their
  acquisition-level forms, approved prior denominators/populations, rejected
  intensity omissions, display/quantitative separation, and canonical CSV
  names.
- **Blocking implementation gaps:** global obsolete display offset; no
  equal-replicate DNA density; no complete sample QC/suppression representation;
  no shared provenance/manifest/appendix; incomplete panel identity.
- **Scientific/authority blockers before relevant slices:** confirmed
  QC/stop-condition record, technical-acquisition rule for biological-replicate
  density, and exact standard canonical percentage/intensity figure inventory.
- **Recommended next task:** complete Slice 0 decision recovery, then authorize
  Slice 1, the every-sample pseudocolor/per-sample display-offset implementation.
- **Experimental inputs accessed:** none.
- **Production code, tests, configs, inputs, reports, artifacts, and Git state
  modified:** none.
- **Audit file created:**
  `docs/implementation/edu_output_contract_gap_audit_2026-08-30.md`.
- **Other files generated:** none.
- **Execution:** `NOT RUN — local execution required`; no R, Python, tests,
  package checks, Quarto, analyses, or renders were executed.
