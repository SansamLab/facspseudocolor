# pH3 Legacy CSV Pilot: Shared Clone-Group Y Axis

**Date:** 2026-09-07
**Status:** Implemented; review and focused execution pending.

## Approved display-only change

For the legacy pH3 CSV pilot, every available pH3-versus-DNA pseudocolor panel
within an explicit clone/replicate group uses the same log10 y-axis limits.
The limits are calculated from the union of finite, positive, offset-restored
display values from that group's available members in the existing 1.6N--4.4N
visual DNA window. No values from another clone group are used.

The group-axis method identifier is
`clone_group_union_positive_display_signal_v1`. The full finite minimum and
maximum are used. A 0.02-log10-unit pad is applied only if every contributing
value is identical, solely to make a valid ordered plotting interval; it does
not remove or clip an event.

## Unavailable members

A sample that cannot support the already-required shared display target or
log10 display remains an explicit unavailable placeholder. It contributes no
signal to the group-axis calculation. If all group members are unavailable, the
group axis is recorded unavailable. If only some members are unavailable, valid
members share an axis derived from the available members and carry the explicit
`available_partial_group_members` provenance status.

## Scientific scope

This changes presentation only. It does not change inputs, sample mappings,
gates, raw or corrected values, background regression, cutoff estimation,
positivity calls, eligibility, outcomes, summaries, or statistical methods.
No experimental inputs were read or modified.

## Verification prepared

The SYNTHETIC pilot tests now cover identical within-group y limits, group
isolation, unavailable-group behavior, and non-mutation of event corrections.
After independent reviews, run only:

```sh
cd "/Users/sansamc/OMRF Dropbox/Chris Sansam/Synchd/Organized Dropbox/Projects and Data/Projects/Cowork_Projects/2026_Analyze_FACS_Data/facs_pseudocolor_workflow"
Rscript -e 'devtools::test(filter = "ph3-legacy-csv-pilot", stop_on_failure = TRUE)'
```

Expected result: PASS with no warnings. Do not render or use experimental data
until focused verification is successful.
