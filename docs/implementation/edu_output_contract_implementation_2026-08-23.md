# EdU Output-Contract Implementation Report
Status: IMPLEMENTED — LOCALLY VERIFIED (FULL TEST SUITE AND R CMD CHECK)

## Verification status

`PASS — user-executed full test suite and R CMD check`

Codex did not execute project code, R or Python scripts, tests, package checks,
Quarto renders, analysis pipelines, example-data analyses, snapshot generation,
or benchmarks while implementing the contract or preparing this report.

The user supplied the following local verification result. Codex recorded the
result without rerunning or independently reproducing it:

```text
Command: R -q -e 'devtools::test(stop_on_failure = TRUE)'
Result: PASS
Date: 2026-08-23
Summary: PASS 387 | FAIL 0 | WARN 0 | SKIP 0
Duration: 166.8 seconds
Failed tests: None
Relevant error excerpt: None
```

The user subsequently supplied this additional local build and check result.
Codex recorded it without running or independently reproducing either command:

```text
Date: 2026-08-24
Branch: feature/edu-output-contract-v2
Command: R CMD build followed by R CMD check --no-manual
R version: 4.5.0
Platform: aarch64-apple-darwin20
Build result: PASS
Package installation: PASS
Package tests during check: PASS
Check result: Status: OK
Errors: 0
Warnings: 0
Notes: 0
Generated artifacts: temporary directory only
```

Together, these user-executed results verify the full test suite and package
build/check for the tested implementation of the confirmed EdU output contract.
They are not approval of any scientific method beyond that confirmed contract
and do not authorize changes to positivity, background modeling, DNA
normalization, phase boundaries, aggregation, statistics, sample mappings,
exclusions, or biological interpretation.

## Repository baseline and current state

- Canonical source directory: `facs_pseudocolor_workflow/`
- Branch at initial implementation-report creation: `main`
- Current branch and user-reported build/check verification branch:
  `feature/edu-output-contract-v2`
- Starting commit: `254e0900442dae7732744545454f38f6e9df2986`
- Current commit: `254e0900442dae7732744545454f38f6e9df2986`
- Commit state: no implementation commit was created; all changes described
  here are in the working tree.
- Tracked working tree: 12 modified files before this report was added.
- New implementation test file:
  `tests/testthat/test-edu-output-contract.R`.
- New durable report: this file.
- Pre-existing untracked file:
  `PROMPT_WHOLE_DATASET_ORACLE_AUDIT.md`; not read as authority, edited, moved,
  or deleted by this task.
- Ignored build artifact present at report creation:
  `facspseudocolor_0.1.0.9000.tar.gz`. It was generated outside Codex during a
  later local build attempt and was not inspected as scientific evidence,
  modified, or deleted.

## Authoritative records and implementation boundary

The implementation used these records:

1. `../edu_output_contract_approved_2026-08-23.md` — authoritative confirmed
   scientific-owner approval record.
2. `../edu_current_implementation_audit_254e090.md` — implementation audit of
   the starting commit.
3. **Agent Handoff: EdU Output-Contract Implementation Plan and Test
   Specification** — the file-level implementation and deterministic-test
   specification supplied as a task attachment. The handoff instructed the
   planning agent not to save a separate plan unless explicitly requested, so
   no repository plan artifact existed at implementation start.

The approval record was treated as authoritative for scientific definitions.
The handoff constrained file scope, table names, shared classification,
compatibility, provenance, report/export effects, implementation slices, and
the focused test matrix. No historical QMD or publication method was used to
reconstruct or reinterpret an approved decision.

## Files inspected

No `.fcs`, `.wsp`, `.prism`, `.xit`, or exported experimental CSV input was
opened. The following files or supplied records were inspected directly by the
implementer or the required independent reviewers.

### Workspace records

- `AGENTS.md`
- `../edu_output_contract_approved_2026-08-23.md`
- `../edu_current_implementation_audit_254e090.md`
- The supplied EdU implementation-plan/test-specification attachment
- The supplied local-verification transcript attachment

### R source and API files

- `R/analysis.R`
- `R/appearance.R`
- `R/artifacts.R`
- `R/config.R`
- `R/input.R`
- `R/normalize.R`
- `R/phase_quantitation.R`
- `R/pseudocolor_helpers.R`
- `R/results.R`

### Test sources and helpers

- `tests/testthat/helper-current-workflow.R`
- `tests/testthat/helper-load.R`
- `tests/testthat/test-analysis.R`
- `tests/testthat/test-appearance-and-artifacts.R`
- `tests/testthat/test-phase-quantitation.R`
- `tests/testthat/test-results-api.R`
- `tests/testthat/test-edu-output-contract.R` after it was added

### Reports, documentation, and manuals

- `inst/quarto/_report-setup.R`
- `inst/quarto/facs_complete.qmd`
- `inst/quarto/facs_quantitation.qmd`
- Directly relevant portions of `inst/quarto/facs_cell_cycle.qmd`,
  `inst/quarto/facs_pseudocolor.qmd`, `inst/quarto/facs_diagnostics.qmd`, and
  `inst/quarto/facs_configurator.qmd`
- `README.md`
- `NEWS.md`
- `docs/CONFIGURATION.md`
- `docs/MIGRATION.md`
- `docs/REPORTS.md`
- `man/quantify_cell_cycle.Rd`
- `man/plot_facs_quantitation.Rd`
- `man/save_facs_results.Rd`

## Files changed

### Modified production source

- `R/phase_quantitation.R`
  - adds the shared EdU event-classification representation;
  - adds acquisition-level collectors for composition, positivity, and
    positive-cell intensity;
  - adds technical-acquisition aggregation with explicit status and provenance;
  - allows an explicit axis title and caption for the legacy composition plot.
- `R/results.R`
  - integrates the seven canonical EdU tables into automatic quantitation;
  - records output-schema and display provenance;
  - preserves and warns about deprecated aliases;
  - keeps calculation independent of panel/display switches;
  - adds explicit canonical CSV export and opt-in deprecated CSV export;
  - updates scientifically explicit quantitation axes and captions.
- `R/appearance.R`
  - labels background-subtracted event displays as using a display offset;
  - labels the alternative divided display as legacy background-divided.

### Modified reports and documentation

- `inst/quarto/facs_complete.qmd`
- `inst/quarto/facs_quantitation.qmd`
- `README.md`
- `NEWS.md`
- `docs/MIGRATION.md`
- `docs/REPORTS.md`
- `man/quantify_cell_cycle.Rd`
- `man/save_facs_results.Rd`

### Modified and new tests

- `tests/testthat/test-analysis.R` — captures the intentional one-per-analysis
  EdU alias deprecation warning.
- `tests/testthat/test-edu-output-contract.R` — new, test-only, explicitly
  `SYNTHETIC` deterministic fixture and focused contract tests.

### Files generated or created

- Created source/test artifact:
  `tests/testthat/test-edu-output-contract.R`.
- Created documentation artifact:
  `docs/implementation/edu_output_contract_implementation_2026-08-23.md`.
- Generated analysis, report, figure, CSV, RDS, snapshot, or benchmark outputs:
  none by Codex.
- Generated package artifact present but not created by Codex:
  `facspseudocolor_0.1.0.9000.tar.gz`; left untouched.

## Shared event-classification implementation

`classify_edu_events()` in `R/phase_quantitation.R` creates one event-level
representation from the already normalized Single Cells data. It records:

- composition and historical-composition eligibility;
- positivity eligibility;
- regional and whole-positive-population intensity eligibility;
- existing computed positivity;
- Early/Mid/Late DNA-region membership;
- historical five-gate assignment;
- Negative-S membership;
- six-gate assignment;
- Unassigned and positivity-Unassigned status;
- unoffset quantitative background-subtracted signal;
- actual configured display signal and transformation;
- display offset and whether it was applied;
- normalized DNA.

The classifier requires the named derived fields and stops when they are
missing. It does not substitute another column or population. It does not
recompute positivity, background, DNA normalization, or existing boundaries.

For a background-subtracted event display, the existing display-offset resolver
is called and the display coordinate is background-subtracted signal plus that
offset. For a legacy divided display, no unused background-subtracted offset is
resolved; the actual divided coordinate is recorded with offset zero. All
canonical quantitative intensity fields continue to use unoffset
`target_bgsub`.

## Decision-to-implementation mapping

| Decision | Implemented calculation or behavior | Functions and fields | Result tables, reports, and exports |
|---|---|---|---|
| Eligibility | Uses validated Single Cells with finite values required by each metric; no undocumented event substitution or filtering | `classify_edu_events()` eligibility fields | Eligibility/source-population fields and counts in all canonical tables |
| EDU-D05a | Preserves G1, Early S, Mid S, Late S, and G2/M assigned-count denominator | `historical_five_gate_assignment`; `collect_edu_composition()` | `edu_assigned_phase_composition`; value `assigned_phase_composition_pct`; legacy `phase_percentages` meaning retained |
| EDU-D05b | Reports six biological categories as percentages of eligible Single Cells; Unassigned remains QC-only | `six_gate_assignment`, `unassigned`, `unassigned_n`, `unassigned_pct_of_eligible_single_cells` | `edu_single_cells_phase_composition`; value `single_cells_phase_composition_pct`; categories may total below 100% |
| EDU-D05c regional | Early/Mid/Late positivity numerator is computed-positive eligible events in the DNA region; denominator is every eligible event in that region | `collect_edu_regional_positivity()`; half-open interval metadata; `zero_denominator` status | `edu_regional_positivity`; value `regional_edu_positive_pct` |
| EDU-D05c overall | Overall positivity uses `[G1 minimum, G2/M maximum)` and retains phase-Unassigned events in the denominator | `collect_edu_overall_positivity()`; `phase_unassigned_n` | `edu_overall_positivity`; value `overall_edu_positive_pct` |
| EDU-D05d | Negative S is computed-negative in `[G1 maximum, G2/M minimum)`; six assigned categories use their own assigned-count denominator | `negative_s`; `six_gate_assignment`; `collect_edu_composition()` | `edu_six_gate_phase_composition`; value `six_gate_phase_composition_pct` |
| EDU-D06a | Median background-subtracted EdU among computed-positive eligible events in Early, Mid, and Late S; sparse regions retain existing minimum-event behavior and receive structured status | `collect_edu_positive_intensity(regional = TRUE)` | `edu_positive_cell_regional_intensity`; value `positive_cell_regional_edu_bgsub_median`; no all-cell regional intensity table |
| EDU-D06b rejected | Does not calculate whole-Single-Cells intensity | Absence enforced by implementation and focused tests | No automatic table, report, or CSV for this rejected metric |
| EDU-D06c | Preserves the whole computed-positive-population median | `collect_edu_positive_intensity(regional = FALSE)` | `edu_positive_population_intensity`; value `positive_population_edu_bgsub_median`; legacy `whole_medians` remains positive-only |
| EDU-D06d | Preserves event-display offset behavior while keeping counts, classifications, and quantitative medians unoffset | `resolve_background_subtracted_offset()` remains unchanged; `classify_edu_events()` records actual display coordinate; `edu_metric_metadata()` and analysis provenance separate transforms | Explicit event and quantitative axis labels; display fields in tables and provenance |
| EDU-D13a | Calculates all seven canonical tables for EdU independently of legacy plot/display switches and retains acquisition-level forms | EdU path in `quantify_cell_cycle()`; `*_acquisition` collectors; `average_edu_metric_table()` | Seven canonical aggregate names plus seven `*_acquisition` names; complete and quantitation reports display both levels |
| EDU-D13b | Adds canonical fields, retains aliases without semantic redirection, warns once, versions schema, and isolates legacy methods | `edu_output_schema_version()` returns `2L`; one warning stored in `analysis$warnings`; `legacy_background_divided` only when explicitly requested | Canonical CSVs require `output_csv_dir`; deprecated CSV aliases additionally require `include_deprecated_csv = TRUE`; Figure 1 output remains absent |

## Canonical output tables and metadata

The automatic modern EdU output set is:

1. `edu_assigned_phase_composition`
2. `edu_single_cells_phase_composition`
3. `edu_six_gate_phase_composition`
4. `edu_regional_positivity`
5. `edu_overall_positivity`
6. `edu_positive_cell_regional_intensity`
7. `edu_positive_population_intensity`

Each has a matching `*_acquisition` table. Depending on metric applicability,
tables record source population, numerator/denominator or source-population
counts, Unassigned counts, half-open DNA intervals, quantitative and display
transforms, actual display offset, offset application, reference-normalization
status, aggregation level/method, output-schema version, and metric status.

Acquisition statuses include `ok`, `zero_denominator`, and
`insufficient_events`. Aggregated tables preserve mixed acquisition quality
using `partial` or `mixed`, `non_ok_technical_n`, and
`technical_metric_statuses` rather than silently converting mixed status to
`ok`.

Technical-acquisition values continue to use the existing unweighted mean.
Counts are summed as provenance. Therefore an aggregate percentage is
intentionally not reconstructed as aggregate numerator divided by aggregate
denominator; the acquisition tables are the exact count-to-value trace.

## Backward compatibility and schema

- Output schema is recorded as integer version `2` in analysis provenance and
  canonical tables.
- `phase_percentages` continues to mean only historical five-gate assigned
  composition.
- `phase_medians` continues to mean only computed-positive Early/Mid/Late-S
  intensity for the explicitly selected compatibility signal.
- `whole_medians` continues to mean only whole computed-positive-population
  intensity for the explicitly selected compatibility signal.
- Corresponding acquisition aliases remain available.
- Aliases are deprecated for one documented major-release compatibility window.
- One concise warning is emitted and stored per analysis; re-quantifying the
  returned analysis does not emit a second copy.
- Old fields are not redirected to whole-Single-Cells composition, regional
  positivity, overall positivity, or whole-Single-Cells intensity.
- Explicit legacy background-divided quantitation is named
  `legacy_background_divided` and is not calculated unless requested through
  the existing signal compatibility control.
- Legacy Figure 1/FlowJo-positive outputs remain outside the automatic modern
  output set.
- Existing saved analyses and historical reports, figures, CSVs, and RDS files
  are not rewritten.

## Report and export changes

- `inst/quarto/facs_complete.qmd` shows the seven canonical aggregate tables,
  the acquisition-level tables, and a separate deprecated-alias section.
- `inst/quarto/facs_quantitation.qmd` shows canonical and acquisition-level
  tables for EdU rather than presenting only the unqualified legacy table.
- Event displays explicitly identify background-subtracted signal with display
  offset, or the alternative legacy background-divided coordinate.
- Quantitative intensity axes explicitly identify computed-positive source
  populations and state that the display offset is excluded.
- Historical phase composition plots explicitly identify the five-gate
  denominator and legacy meaning.
- `save_facs_results()` accepts an explicit `output_csv_dir` and writes the
  seven canonical aggregate tables plus their acquisition-level forms.
- Deprecated CSV aliases are omitted by default and are exported only when
  `include_deprecated_csv = TRUE` is explicitly requested.
- All requested CSV names are validated before directories or files are
  written, preventing partial export caused by a missing requested alias.

## Focused tests and protected invariants

All newly constructed events are labeled `SYNTHETIC` and live only in
`tests/testthat/test-edu-output-contract.R`. They encode arithmetic and boundary
invariants, not biological conclusions.

| Test | Scientific or operational invariant protected |
|---|---|
| `EdU contract exposes only the seven approved canonical tables` | Exact modern output set, acquisition-level availability, schema version, and absence of rejected intensity/Figure 1 outputs |
| `Negative S uses approved half-open boundaries` | Negative S includes G1 maximum and excludes G2/M minimum without changing historical G1/G2/M assignment |
| `composition denominators and Unassigned QC reconcile exactly` | Historical five-gate denominator, whole-Single-Cells denominator, assignment-specific Unassigned QC, displayed total below 100%, and six-gate total of 100% |
| `regional and overall positivity retain every eligible denominator event` | Regional denominators include all eligible interval events; overall denominator includes phase-Unassigned events |
| `regional and overall intervals are exactly half-open` | Exact inclusion/exclusion at 775, 1265, 1620, and 2125 for the default synthetic scale |
| `zero regional denominators are structured and do not invent values` | Empty regions return `NA` with `zero_denominator`; sparse intensity returns `NA` with `insufficient_events` |
| `canonical intensities preserve positive-only legacy values` | Canonical regional/whole background-subtracted medians exactly equal legacy positive-only values; five-gate values remain equal |
| `display offset is recorded but cannot change quantitative values` | Changing display offset changes only display signal, not canonical medians or quantitative coordinates |
| `legacy divided display records its actual coordinate and zero offset` | Alternative display provenance matches actual divided signal and does not require or claim an unused offset |
| `deprecated aliases warn once and plot switches cannot change calculations` | One-warning-per-analysis behavior and calculation independence from legacy plot switches |
| `technical acquisitions retain unweighted summary aggregation` | Existing unweighted technical-acquisition mean remains unchanged while counts remain available |
| `technical aggregation preserves mixed acquisition status` | A valid acquisition cannot hide another acquisition's zero denominator or insufficient status |
| `modern and explicitly requested legacy signal populations remain separate` | Modern computed-positive outputs remain distinct from legacy divided and Figure 1 populations |
| `CSV export uses canonical names and requires explicit legacy selection` | Canonical export names, aggregate/acquisition availability, and opt-in deprecated CSV behavior |

`tests/testthat/test-analysis.R` was also updated to expect the deliberate EdU
alias deprecation warning from the high-level analysis entry point. No numerical
expectation, fixture, mapping, method, or biological claim in that test was
changed.

The full test suite result recorded above was executed and supplied by the user,
not by Codex. The reported result was `PASS 387`, `FAIL 0`, `WARN 0`, and
`SKIP 0` on 2026-08-23.

## Methods and protected artifacts not changed

- Positivity classification: unchanged. The classifier consumes the existing
  `edu_computed_positive` field and never recomputes the boundary.
- Background modeling: unchanged.
- DNA normalization: unchanged.
- Existing G1, Early S, Mid S, Late S, and G2/M boundaries: unchanged.
- Display-offset calculation: unchanged for background-subtracted displays.
- Quantitative medians: continue to exclude the display offset.
- Technical-acquisition aggregation calculation: unchanged unweighted mean;
  only explicit provenance and mixed-status reporting were added.
- Biological-replicate summary behavior: unchanged.
- Reference normalization behavior: unchanged and reported separately.
- Statistical methods: unchanged; no inferential method was added.
- Sample-to-condition mappings: unchanged.
- Biological and technical replicate membership: unchanged.
- Exclusions: unchanged; no structured exclusion feature was added.
- Channels: unchanged.
- Gates and thresholds: no existing gate or threshold was changed. The approved
  Negative-S output category was added from existing configured G1/G2/M bounds
  and existing computed negativity.
- Historical artifacts: no existing HTML, PDF, PPTX, PNG, CSV, RDS, WSP, FCS,
  Prism, or publication artifact was rewritten.
- Experimental inputs: none read or modified.

## Independent scientific-integrity review

The scientific-integrity reviewer inspected the actual diff and surrounding
code read-only and did not execute project code or open experimental inputs.

Findings and resolutions:

1. **Blocking — historical Unassigned QC undercount.** The first implementation
   used six-gate Unassigned status for every composition table, so historical
   five-gate QC omitted Negative-S events outside that metric's denominator.
   Resolution: Unassigned is now assignment-specific; the synthetic historical
   fixture requires 9 events outside the five-gate assignment.
2. **Important — mixed technical status hidden.** An aggregate was initially
   labeled `ok` whenever any acquisition had a finite value. Resolution:
   `partial`/`mixed`, total acquisition count, non-OK count, and source status
   ledger now preserve the failed/empty acquisition.
3. **Blocking — divided-display provenance mismatch.** The initial shared record
   always described background-subtracted-plus-offset display even when the
   configured display was divided. Resolution: the actual display signal,
   transform, zero offset, and application flag are now recorded for each mode.
4. **Advisory — exact boundary test gaps.** Initial tests did not place events
   exactly at all regional/overall joins. Resolution: synthetic events now test
   exact 775/2125 and 1265/1620 half-open behavior.

Final scientific re-review found no remaining blocking or important findings.
The reviewer noted only that aggregate percentages intentionally use the
approved unweighted mean while aggregate counts are summed; acquisition-level
tables and aggregation metadata make that behavior explicit.

## Independent artifact/security review

The artifact/security reviewer inspected the actual diff and surrounding code
read-only and did not execute project code.

Findings and resolutions:

1. **Blocking — incorrect display provenance for legacy divided mode.** Resolved
   together with the scientific finding by recording the actual divided signal,
   transform, and zero offset.
2. **Important — possible partial deprecated CSV export.** Requested deprecated
   names were initially appended after canonical validation. Resolution: the
   complete requested name set is validated before path creation or any write.
3. **Important — unused offset coupled normalized display to plotting data.**
   Offset resolution initially ran for every EdU quantitation, which could make
   a divided-display analysis fail on an offset it did not use. Resolution: the
   resolver is called only for a configured background-subtracted display;
   divided display uses zero without requiring anchor/baseline/raw fallback
   candidates.

Final artifact/security re-review found no remaining blocking or important
findings. It found no new network calls, telemetry, secrets, shell execution,
destructive I/O, unsafe absolute paths, hidden output fallback, generated/cache
additions to the change, or unrelated tracked edits.

## Assumptions, warnings, and unresolved uncertainty

- The exact schema-version value was a technical implementation choice because
  the approval required versioning without assigning a literal. Integer `2`
  was selected to distinguish the additive contract schema from the prior
  output shape; this does not version or approve the scientific method itself.
- `partial` and `mixed` are operational aggregate statuses used to avoid hiding
  acquisition-level `zero_denominator` or `insufficient_events` states.
- The supplied handoff was used as the implementation/test specification even
  though it was not saved as a repository plan artifact.
- A concise alias deprecation warning is expected once per EdU analysis. A
  repeated warning for the same returned analysis is not expected.
- Local verification is based on the user-supplied full-suite and package
  build/check results recorded above. Codex did not execute or independently
  reproduce those results.
- Passing tests confirm only the tested EdU output-contract implementation;
  they do not approve any scientific method outside the confirmed contract.
- Missing required event fields remain a hard error; no fallback or substitute
  field is used.
- The structured exclusion-mechanism design remains out of scope.
- Positivity-boundary/background-model choices, sparse nonzero-region policy,
  statistical inference, and other decisions explicitly excluded by the
  approval remain unresolved and unchanged.
- Implementation changes remain uncommitted.

## Remaining work

1. Decide whether to commit the implementation and this report; the changes
   remain in the working tree.
2. Retain the exact commands below for reproducible future verification. Any
   later failure, warning, unexpected check note, or unexpected generated file
   must be investigated.
3. Remove or archive the ignored package tarball only through an explicit
   repository-hygiene decision; this report does not authorize deleting it.

## Exact local verification instructions

### Working directory

```text
/Users/sansamc/OMRF Dropbox/Chris Sansam/Synchd/Organized Dropbox/Projects and Data/Projects/Cowork_Projects/2026_Analyze_FACS_Data/facs_pseudocolor_workflow
```

### Prerequisites

- R 4.2 or newer.
- Package development dependencies already installed, including `devtools` and
  `testthat`.
- A shell that supports quoted paths and `mktemp` as shown.
- No production experimental input needs to be added or opened for the focused
  synthetic test.

### Commands

From the working directory, run the focused contract test:

```bash
R -q -e 'devtools::test(filter = "edu-output-contract", stop_on_failure = TRUE)'
```

Then run the full package test suite:

```bash
R -q -e 'devtools::test(stop_on_failure = TRUE)'
```

Build and check in a temporary directory so package artifacts are not written
into the repository:

```bash
PROJECT_DIR='/Users/sansamc/OMRF Dropbox/Chris Sansam/Synchd/Organized Dropbox/Projects and Data/Projects/Cowork_Projects/2026_Analyze_FACS_Data/facs_pseudocolor_workflow'
VERIFY_DIR="$(mktemp -d /private/tmp/facspseudocolor-edu-contract.XXXXXX)"

(
  cd "$VERIFY_DIR" &&
  R CMD build "$PROJECT_DIR" &&
  R CMD check --no-manual facspseudocolor_0.1.0.9000.tar.gz
)
```

### Expected success conditions

- Focused contract tests finish with `FAIL 0` and `WARN 0`.
- Full package tests finish with `FAIL 0` and `WARN 0`.
- `R CMD build` creates
  `$VERIFY_DIR/facspseudocolor_0.1.0.9000.tar.gz`.
- `R CMD check` creates `$VERIFY_DIR/facspseudocolor.Rcheck/`.
- Package check finishes with `Status: OK`, zero errors, and zero warnings.
- No new analysis output, report, figure, CSV, RDS, cache, or package artifact
  appears in the repository working tree.

### Failures or warnings that require stopping

Stop and report rather than continuing or changing scientific methods if any of
the following occurs:

- any failed or skipped focused contract test;
- any failed full-suite test;
- any unexpected warning from tests;
- any package-check error or warning;
- any unexpected package-check note;
- a missing canonical table or metadata field;
- a legacy alias that differs numerically from its approved canonical source;
- any unexpected generated file in the repository;
- any missing input, population, channel, configured boundary, or required
  event field.

Do not weaken validation, alter a boundary, change a denominator, substitute an
input, or suppress a warning to obtain a passing result.

## Local verification result template

```text
EdU output-contract local verification
Date:
Operator:
Working directory:
Starting commit/worktree state:

Focused command:
Result: PASS | FAIL
Failed test names: none | <names>
Relevant error excerpt: none | <short excerpt>

Full-suite command:
Result: PASS | FAIL
Failed test names: none | <names>
Relevant error excerpt: none | <short excerpt>

Build/check command:
Result: PASS | FAIL
Check status: OK | ERROR | WARNING | NOTE
Relevant excerpt: none | <short excerpt>

Expected generated files:
- <temporary tarball path>
- <temporary .Rcheck directory>

Unexpected generated files: none | <paths>
Unexpected warnings or notes: none | <short description>

Overall verification: PASS | FAIL
Scientific-owner follow-up required: no | yes — <reason>
```

Do not paste large successful logs. Report the commands, `PASS` or `FAIL`,
failed test names, a relevant error excerpt, and any unexpected generated files.
