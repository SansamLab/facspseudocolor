# pH3 raw-4N density cutoff — contract amendment and implementation plan

**Date:** 2026-09-01  
**Status:** Owner-confirmed scientific rule and fixed production method settings.  
**Replaces:** the uncommitted FlowJo-label threshold-calibration diagnostic as the proposed positivity-calibration route.

## Purpose

This amendment makes pH3 positivity a deterministic, control-calibrated call
without using the FlowJo `pH3 Positive` gate as a numerical authority or
training label.  It does not change the approved DNA normalization, 2N–4N
eligibility, 4N/below-4N definitions, event-specific pH3 background
subtraction for signal outcomes, or the four pH3 report outcomes.

For every explicit biological-replicate set, derive one cutoff from that
set's matched untreated control and apply the *same raw-signal cutoff* to
every condition in that set, across all DNA regions.  The six-sample Cdk1as
experiment is an expected first use, but no sample names or experimental
files are encoded by this contract.

## Owner-confirmed scientific rule

1. Select the matched untreated control's retained **4N events**.  "4N" means
   the existing owner-approved computed 4N membership; it is not a new
   FlowJo gate and it is not re-estimated here.
2. Use the event's **raw pH3 fluorescence** on its recorded linear scale.  Do
   not log-transform, DNA-background-correct, offset, normalize, rescale, or
   impute it before density estimation.
3. Estimate a one-dimensional density from those values.  Identify the one
   dominant negative-cell peak.  The cutoff is the first right-side local
   minimum after that peak: the point at which the peak returns to baseline.
4. Apply that cutoff unchanged to the control and every matched treatment in
   that replicate set.  Positivity is `raw_pH3 > cutoff`; equality is negative.
   The cutoff applies across all retained DNA regions, including 4N and
   below-4N.
5. If a unique dominant peak or its qualifying right-side local minimum
   cannot be established deterministically, stop the entire replicate set and
   report it as unavailable.  Do not borrow a cutoff, pool controls, use a
   treatment-derived cutoff, use FlowJo membership, choose a nearby valley,
   or fall back to a quantile.
6. DNA-dependent background-corrected pH3 remains the sole basis for C/D
   pH3 signal medians and their reference-relative values.  It must never be
   used for the positivity cutoff under this amendment.  The required
   noncircular background-membership sequence is resolved in the owner-decision
   section below; implementation must not assume that the old FlowJo member
   flag remains a background-negative label.

## Deterministic estimator and selection contract

The following are technical lab-standard defaults, not biological choices.
They must be implemented as immutable, versioned method settings and emitted
in provenance.  They are proposed because the repository has no established
one-dimensional raw-pH3 cutoff estimator; the implementation must not reuse
the two-dimensional pseudocolor bandwidth setting.

| Item | Required deterministic default |
| --- | --- |
| Input rows | Exactly all retained control events with `eligible_2to4n == TRUE`, `four_n_member == TRUE`, and valid direct identity. No FlowJo pH3 membership filter; every selected raw pH3 value must then be finite or the set stops. |
| Raw scale | Recorded numeric linear fluorescence; density calculations do not use `log10()`. |
| Minimum support | At least 100 qualifying 4N raw-signal events. Fewer stops the set. |
| Estimator | `stats::density()` with `kernel = "gaussian"`, `bw = "nrd0"`, `adjust = 1`, `n = 2048`, and `from`/`to` fixed to the observed finite control 4N raw-signal minimum/maximum. |
| Grid extrema | An interior local maximum has `y[i] > y[i-1]` and `y[i] >= y[i+1]`; an interior local minimum has `y[i] < y[i-1]` and `y[i] <= y[i+1]`. This left-to-right tie convention is deterministic. Endpoints never qualify. |
| Dominant peak | The unique local maximum with the greatest density value. Equal greatest values stop the set as ambiguous; a monotonic or boundary-only density also stops. |
| Cutoff | First local minimum at a grid index greater than the dominant-peak index. Its grid `x` value is the raw cutoff. A later, deeper minimum is not substituted. |
| Application universe | Every retained valid Single-Cell event in the matched control or treatment acquisition, across all DNA regions. A raw finite event is positive only when `raw_pH3 > cutoff`; equality is negative. Nonfinite raw pH3 is an explicit unavailable event, never called negative or positive. |

The first local-minimum definition is the owner-confirmed operational meaning
of “returns to baseline.”  No additional density-height fraction is proposed;
adding one would change the approved rule and requires a new owner decision.

### Mandatory QC and stop rules

The caller must fail closed for the replicate set when any condition below is
true:

- explicit replicate-set/control mapping is absent, duplicated, ambiguous, or
  does not contain exactly one control;
- control 4N membership/identity/provenance is not retained and valid;
- fewer than 100 qualifying raw control 4N events;
- any selected control raw pH3 value is nonfinite, or density output/grid is
  nonfinite, non-increasing, wrong-sized, or otherwise invalid;
- the selected control raw-pH3 range has zero width;
- no interior local maximum, no unique density-dominant local maximum, no
  right-side local minimum, or a cutoff outside the observed finite control
  4N raw-signal range;
- any target sample in the set lacks retained direct identity, raw pH3, or
  required DNA memberships needed to apply the frozen cutoff.

No automatic smoothing retry, bandwidth escalation/de-escalation, alternate
peak selection, clipping, or cross-set fallback is permitted.  A technical
diagnostic may show the density and failure reason, but must not create a
callable threshold after a stop.

## Required provenance, QC, and report artifacts

For each replicate set, persist one immutable cutoff record with at least:

- method ID/version and all estimator settings;
- control sample/acquisition identities; matched condition and replicate-set
  identities; raw channel identity; event-selection predicate;
- direct-identity/export-operation/full-manifest/config/analysis bindings;
- qualifying control event count; raw range; density grid size; dominant-peak
  index/raw coordinate/density; cutoff index/raw coordinate/density;
- `cutoff_status` (`available` or a specific stop reason); and
- exact target sample IDs to which the cutoff was applied.

Persist per-acquisition positivity-application QC: total retained events,
finite raw events, nonfinite raw events, cutoff, called-positive count, and
called-negative count.  Reconciliation must show that all finite applicable
events are partitioned once, without changing event order or direct identity.

The report must include, for every replicate set, a raw-pH3 4N control density
plot with the dominant peak and cutoff visibly marked.  A stopped set receives
a clear no-cutoff/failed-QC panel rather than a guessed line.  Outcome panels
must name the cutoff method and show unavailable replicate-set values rather
than silently omitting them.  These are future output additions; this
amendment does not authorize their implementation.

## Effects on the existing Slice 1–5 pipeline

| Current component | Required amendment |
| --- | --- |
| Slice 1 output model (`R/ph3_outputs.R`) | Add a versioned cutoff-method declaration and exact cutoff/provenance schema; do not use `flowjo_owner_approved_positive_population_v1` as positivity authority for the new computed call. Existing direct-identity export binding remains mandatory. |
| Slice 2 background (`R/ph3_background.R`) | The regression estimator itself is unchanged, but its negative membership cannot remain FlowJo-derived. After the raw cutoff call is frozen, fit the DNA-dependent background regression from computed-negative **eligible 2N–4N** events. Raw pH3 remains independently retained for cutoff selection; background output is not an input to the density. |
| Slice 3 signal outcomes (`R/ph3_signal.R`) | Recompute `ph3_positive_member` for every retained valid Single-Cell event from the frozen raw cutoff, retain FlowJo membership only as optional non-authoritative provenance if present, and continue using corrected/raw fallback basis only for C/D signals. |
| Slice 4 report model (`R/ph3_report_model.R`) | Bind every A–D result to the cutoff record/status and make stopped sets explicit. Do not summarize across unavailable sets as though they were zero. |
| Slice 5 plots/Quarto (`R/ph3_report_plots.R`, `inst/quarto/facs_ph3_output_contract.qmd`) | Add required raw-control 4N density/cutoff QC panels and method/provenance text. Existing four outcome panels consume only the validated revised model. |
| Current legacy pH3 route (`R/ph3.R`, `docs/PH3_MODE.md`) | Leave unchanged unless a separately approved migration routes production pH3 through the output-contract model. Legacy FlowJo-gate semantics must remain labeled legacy, never silently changed. |
| Direct identity/export contract | Still required for production use. It establishes event/sample/gate provenance and protects mappings; it no longer needs an exported `pH3 Positive` child population for numerical cutoff calibration. Required populations must be revised deliberately to distinguish upstream DNA/Single-Cell provenance from obsolete pH3-positive membership. |

The last row is an implementation boundary, not authorization to relax the
direct-identity contract.  Any export schema change must keep exact source
event identity, sample binding, ordered channels, workspace coverage,
manifest hashes, and fail-closed containment.

**Transitional compatibility boundary:** current production direct-identity
ingestion still requires `complete`, `g1`, and `ph3_positive` exports. Until a
later validated export-contract migration changes that requirement, retain the
`ph3_positive` child artifact only as optional, non-authoritative provenance;
the new computed call must not claim that it can omit the child before the
validator and manifest schema are changed and verified.

## Deprecation of uncommitted calibration/export-preparation artifacts

The following current untracked work is incompatible with this amendment:

- `R/ph3_threshold_calibration.R`
- `man/diagnose_ph3_threshold_calibration.Rd`
- `man/write_ph3_threshold_calibration.Rd`
- `tests/testthat/test-ph3-threshold-calibration.R`
- `docs/implementation/ph3_threshold_calibration_diagnostic_2026-08-31.md`
- external `test_reports/ph3_cdk1as_verified_export_calibration_2026-08-31/`
  runbook/template artifacts, if present.

They must **not** be staged as production work.  Preserve their history by
leaving them untracked and adding a short dated supersession note only in this
amendment (or in the final implementation record).  Do not delete user files
or derived folders without explicit user instruction.  The synthetic FlowJo
index fixture and direct-identity infrastructure are not deprecated; only the
FlowJo-positive-label calibration/runbook path is.

## Minimal implementation sequence

1. **Decision freeze and contract/schema slice.** Resolve the two named
   owner/method-owner decisions below. Then add the raw-density cutoff method
   ID, frozen estimator settings, cutoff and application-QC schemas, and
   configuration validation. Explicitly distinguish computed pH3 membership
   from legacy FlowJo membership.
2. **Pure cutoff engine slice.** Implement input validation, density fit,
   extrema selection, exact per-set control mapping, frozen application, and
   fail-closed statuses using unmistakably SYNTHETIC in-memory fixtures.
3. **Background/output-model integration slice.** Reconcile computed
   membership to event identities; fit background from computed-negative
   eligible 2N–4N events; regenerate A/B and C/D inputs from it; bind cutoff
   provenance; and prevent mixed or stopped-set summaries.
4. **QC/report slice.** Add density/cutoff panels, per-set status display,
   machine-readable cutoff tables, and thin Quarto consumption.  No numerical
   calculations belong in Quarto.
5. **Migration/real-data readiness slice.** Deliberately revise the verified
   export specification to remove the pH3-positive child as a numerical
   requirement while preserving all upstream provenance.  Prepare a separate
   derived-data runbook; do not alter historical exports, FCS, or WSP files.

Each slice requires focused synthetic tests, local human execution, and
separate read-only scientific-integrity and artifact/security review before
the next slice.  No repository code, test, export, or real-data analysis was
executed for this planning amendment.

## Remaining required decisions and assumptions

The following owner-confirmed science is fixed: raw scale, matched untreated
control, 4N calibration selection, first right-side local minimum, unchanged
application to all conditions in a set, all-DNA-region membership calling, and
fail-closed behavior.

Two decisions are unavoidable before a callable production implementation:

1. **Noncircular background-negative definition.** The existing background
   engine fits negatives as `eligible_2to4n & !ph3_positive_member`. Once
   positivity is recalculated, retaining the FlowJo flag there would violate
   its non-authority; replacing it changes the regression/C-D values. The plan
   proposes this sequence: raw 4N-control cutoff -> computed positivity for
   all retained valid Single-Cell events -> regression on computed-negative
   eligible 2N–4N events -> corrected C/D medians. The owner must confirm this
   consequence before production work.
2. **Versioned density settings.** The numerical settings in the table are
   proposed technical defaults, not historical lab standards. A method owner
   must approve them as `ph3_raw_4n_density_cutoff_v1` (or provide a different
   exact parameter set) before an engine produces a callable real-data cutoff.
   A design-only synthetic prototype could be built before that approval, but
   must not be connected to production analysis or real-data reporting.

Technical-acquisition policy is made fail-closed without a biological choice:
each explicit matched control sample must have exactly one acquisition for
cutoff derivation. Multiple control acquisitions stop the set until a separate
approved policy specifies how to combine them; they must not be event-pooled
or silently selected. The cutoff record therefore includes exactly one control
acquisition ID. A target sample may also have only one acquisition in the
initial implementation; multi-acquisition extension requires its own explicit
application/aggregation contract.

Experimental inputs read: none. Experimental inputs modified: none. Source
or test execution: **NOT RUN — local execution required**. This document does
not authorize code changes, data exports, staging, commits, or deletion.

## 2026-09-01 confirmation addendum

The owner confirmed the two previously pending items: use the noncircular
sequence `raw 4N-control cutoff -> computed positivity -> computed-negative
eligible 2N–4N background regression -> corrected C/D medians`, and freeze
the table's Gaussian/nrd0/adjust-1/2,048-grid/minimum-100 settings as
`ph3_raw_4n_density_cutoff_v1`. The FlowJo-positive-label calibration path is
superseded and must remain untracked rather than being deleted.
