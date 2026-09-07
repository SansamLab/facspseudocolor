# EdU report Slice 3: reference-relative intensity

**Date:** 2026-09-07
**Status:** IMPLEMENTED — AWAITING INDEPENDENT REVIEW

## Owner-approved scope

The two EdU-positive intensity report figures now show the completed
background-subtracted biological-replicate medians relative to the explicitly
configured matched `Untreated` reference. The two percentage figures remain
absolute percentages.

This is a report-only fold calculation: for each biological replicate and,
for the regional panel, each established Early/Mid/Late-S category, the
canonical biological-replicate median is divided by the canonical median for
that replicate's configured `Untreated` condition. The matched Untreated value
is therefore 1. Reference-relative calculations occur only after the
established unweighted technical-acquisition mean; neither event data nor a
background fit, positivity call, DNA region, sample mapping, or reference is
recomputed or selected anew.

## Fail-closed validation

The report requires exactly one configured reference in every
biological/technical acquisition pair and one consistent reference condition
within each biological replicate. It stops if the reference is missing,
duplicated, inconsistent across technical acquisitions, or absent from an
intensity/category group. It continues to require the canonical
acquisition-to-biological-replicate reconciliation before any relative
calculation. A finite but nonpositive background-subtracted reference median
does not define a meaningful fold: that group is retained with an explicit
`unavailable_invalid_reference_intensity` status and `NA` fold values, rather
than falling back to raw, display-offset, or direct values. A nonfinite
numerator is similarly explicit and unavailable.

## Files

Modified:

- `R/edu_pseudocolor_output_contract.R`
- `inst/quarto/facs_edu_pseudocolor_output_contract.qmd`
- `tests/testthat/test-edu-pseudocolor-output-contract.R`

Generated:

- this implementation record.

Experimental inputs read or modified: none.

Sample mappings, exclusions, gates, thresholds, background correction,
normalization of canonical results, statistics, biological claims, and
canonical quantitative exports changed: none. The report-only display changes
the two intensity plots from direct completed medians to their approved
within-biological-replicate reference-relative fold values.

## Verification

Focused synthetic verification is required after independent review:

```sh
cd "/Users/sansamc/OMRF Dropbox/Chris Sansam/Synchd/Organized Dropbox/Projects and Data/Projects/Cowork_Projects/2026_Analyze_FACS_Data/facs_pseudocolor_workflow"
Rscript -e 'devtools::test(filter = "edu-pseudocolor-output-contract", stop_on_failure = TRUE)'
```

Expected: zero failures, errors, or unexpected warnings; no project-tree
outputs. Stop on a failed reference identity/coverage/reconciliation check or
any changed canonical table.
