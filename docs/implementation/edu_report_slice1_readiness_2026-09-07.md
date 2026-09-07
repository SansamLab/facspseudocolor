# EdU report Slice 1 readiness and handoff

**Date:** 2026-09-07
**Status:** Owner-authorized, implementation not started

## Purpose

Define the smallest EdU report slice that can be implemented now without
choosing a new scientific method: truthful, individually identifiable
DNA-versus-EdU pseudocolor panels with the owner-confirmed per-sample,
display-only background-restoration offset.

This is deliberately not an EdU report redesign. It does not implement the
later DNA-density, standard quantitative-figure, provenance-manifest, or export
package work.

## Authority

The scientific-owner-confirmed decisions are recorded outside the repository in:

- `/Users/sansamc/Documents/Codex/2026-08-30/edu-output-first-design/outputs/edu-desired-analyses-and-figures-specification-2026-08-30.md`
- `/Users/sansamc/Documents/Codex/2026-08-30/edu-implementation-decision-recovery/outputs/edu-implementation-decision-recovery-2026-08-30.md`
- `/Users/sansamc/Documents/Codex/2026-08-30/edu-implementation-decision-recovery/outputs/edu-code-architect-slice-1-handoff-2026-08-30.md`

The repository audit that maps these decisions to current implementation gaps is
`docs/implementation/edu_output_contract_gap_audit_2026-08-30.md`.

## What the package can report now

The current EdU engine already calculates, without any new mathematical method:

1. Seven canonical aggregate tables and matching acquisition-level tables:
   assigned-phase composition, Single-Cells composition, six-gate composition,
   regional and overall positivity, regional positive-cell intensity, and
   whole-positive-population intensity.
2. Existing per-sample DNA-versus-EdU pseudocolor panels.
3. Existing standard quantitative plots for historical five-gate assigned
   composition, regional computed-positive intensity, and whole
   computed-positive intensity.

The existing automatic display path is not suitable as the standard EdU report
because it uses a global rounded/fixed offset rather than the owner-confirmed
sample-specific restoration rule, and its panel identity/QC presentation is
incomplete.

## Bounded Slice 1 implementation contract

For every plottable sample/acquisition:

- Generate one separately identifiable DNA-versus-EdU pseudocolor panel.
- Identify condition, biological replicate/sample, technical acquisition where
  applicable, axes, and an established QC warning or suppression reason.
- Retain proof of the exact eligible negative/background events used by the
  already-established per-sample background fit.
- Compute only for display:

  `median(raw exact negative events) - median(background-corrected exact negative events)`.

- Apply that keyed offset only to that panel's display coordinate.
- Resolve automatic display limits from the actual displayed coordinate.
- If the required membership proof or offset cannot be resolved, retain the
  analytical result but visibly suppress that panel with its explicit reason;
  do not substitute another population, sample, transform, offset, or panel.

The slice must leave positivity calls, background fits, DNA normalization,
gates, thresholds, corrected medians, reference normalization, statistics,
canonical tables, and quantitative CSV values unchanged.

## Out of scope for Slice 1

- Equal-technical-acquisition then equal-biological-replicate DNA-density
  calculation and plot.
- New or promoted percentage/intensity figures.
- A common provenance schema, manifest, appendix, or atomic report/export
  package.
- New QC thresholds, QC criteria, or reason vocabulary.
- Any experimental input access or real-data render.

## Required implementation and review sequence

1. Inspect current EdU display/background-fit ownership and identity fields;
   stop if the exact fitted negative-event membership cannot be preserved.
2. Implement only the bounded display contract and deterministic `SYNTHETIC`
   tests. Do not access experimental inputs.
3. Obtain independent read-only scientific-integrity review of the complete
   diff.
4. Obtain independent read-only artifact/security review of the complete diff.
5. Under the current temporary execution exception, run only the focused Slice
   1 tests after both reviews, followed by the affected EdU/full suite if the
   focused test passes. All artifacts must remain outside the repository.

## Initial test requirements

Use unmistakably `SYNTHETIC` fixtures to prove:

- different samples receive different offsets exactly matching the formula;
- offset identity cannot be reassigned after manifest/event reordering;
- quantitative values and exports remain byte-for-byte/value-for-value
  unchanged by offset resolution;
- panel labels remain correct with multiple technical acquisitions;
- missing/mutated membership proof, offset inputs, or established status
  suppresses the affected display explicitly and never substitutes data;
- limits use the displayed coordinate; and
- no global, rounded, fixed, pooled, or fallback offset is used in the standard
  path.

Suggested focused command after review:

```sh
Rscript -e 'devtools::test(filter = "edu-pseudocolor-output-contract", stop_on_failure = TRUE)'
```

## Later owner decisions

No additional owner decision is required for Slice 1. The exact density grid
and bandwidth, report package/provenance serialization, and later report layout
are Code Architect mechanics for later slices; no later slice is authorized by
this readiness record.

## Audit facts for this record

- Canonical source: `facs_pseudocolor_workflow/`.
- Experimental inputs read or modified: none.
- Production/test code modified: none.
- Generated data, reports, figures, or package artifacts: none.
- Tests/renders executed: none.
- Scientific methods changed: none.
