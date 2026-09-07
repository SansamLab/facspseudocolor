# EdU report Slice 2: density palette, positivity, and intensity plots

**Date:** 2026-09-07
**Status:** IMPLEMENTED — AWAITING INDEPENDENT REVIEW

## Owner-approved scope

This slice changes the EdU output-contract report only.

- Every DNA-versus-EdU pseudocolor panel uses the project-standard refined
  density palette already used by the preferred pH3 report. This is solely a
  color rendering change.
- A **2N–4N EdU-positive percentage** plot is added. For each technical
  acquisition, its numerator is the computed EdU-positive events whose
  normalized DNA is in the inclusive interval `[2N, 4N]`; its denominator is
  every positivity-eligible event in that same interval.
- A separate **Early/Mid/Late-S EdU-positive percentage** plot is added from
  the canonical `edu_regional_positivity` tables. It retains each established
  regional denominator and is intentionally not called “All S”.
- Two **EdU-positive intensity** figures are added: one for all computed
  EdU-positive events and one for Early/Mid/Late-S computed EdU-positive
  events. Both use only the canonical background-subtracted positive-cell
  intensity tables already retained by the completed analysis.

Both plots retain the established unweighted technical-acquisition mean within
each biological replicate. Condition bars summarize those biological-replicate
values; error bars and individual points follow the existing
`quant_error_bar` and `quant_show_points` configuration settings.

## Validation and traceability

The report refuses to build when either canonical regional positivity table is
missing, incomplete, duplicated, or does not exactly reconcile from the
acquisition rows to the retained biological-replicate table. The newly derived
2N–4N rows require exact manifest order, event-row identity, DNA value,
computed-positivity value, and positivity-eligibility reconciliation between
the retained event classification and its normalized source events.

Every canonical acquisition row is also checked against the current manifest's
exact five-field identity—biological replicate and index, technical
acquisition, condition and index. Regional tables must contain that identity
multiset once per approved Early/Mid/Late-S region; all-positive intensity has
it once per acquisition. Omitted, duplicated, or substituted acquisition rows
therefore stop the report before numerical re-aggregation.

Focused SYNTHETIC regression coverage independently mutates regional positivity,
all-positive intensity, and regional intensity acquisition rows. It confirms
that substituted, omitted, and duplicated manifest identities fail closed;
regional intensity covers all three identity failures explicitly.

No panel recalculates background correction, positivity, regions, or the
canonical regional values. The 2N–4N report table is derived from the existing
completed classification only, with its definition and inclusive endpoints in
the returned provenance. Intensity figures use no event data at all: their
source and biological-replicate tables must exactly reconcile from the existing
canonical acquisition tables or the report stops.

## Files

Modified:

- `R/edu_pseudocolor_output_contract.R`
- `inst/quarto/facs_edu_pseudocolor_output_contract.qmd`
- `tests/testthat/test-edu-pseudocolor-output-contract.R`

Generated:

- this implementation record only.

Experimental inputs read or modified: none.

Sample mappings, exclusions, gates, thresholds, background correction,
normalization, statistics, biological claims, and canonical quantitative
exports changed: none. The only new numerical display is the approved 2N–4N
percentage table/plot, calculated from the retained classification; it does
not alter the analysis object.

## Verification

Executed under the temporary owner-approved execution exception:

```sh
Rscript -e 'devtools::test(filter = "edu-pseudocolor-output-contract", stop_on_failure = TRUE)'
```

Result: PASS — 23 assertions, zero failures, warnings, and skips. This run used
only unmistakably SYNTHETIC in-memory test events. No real analysis, experimental
input, or Quarto render was executed.

The intensity extension was added after this run and has not yet been executed;
the complete corrected-tree verification remains required.

Still required after review: run the existing EdU output-contract regression,
the full package suite, and a temporary-directory report render using an
explicitly configured real experiment. Stop on any changed canonical table,
membership/identity mismatch, unexpected warning, or generated repository
artifact.

## Assumptions and remaining uncertainty

The terms “2N–4N” are implemented as inclusive numeric DNA endpoints. This
does not redefine the existing half-open regional Early/Mid/Late-S intervals.
No unresolved scientific decision remains for this slice.
