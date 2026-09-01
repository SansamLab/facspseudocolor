# pH3 raw-4N density cutoff implementation

**Date:** 2026-09-01  
**Status:** Implemented, independently re-reviewed, pending local verification.

## Background-field compatibility correction (2026-09-01)

The user-reported background-regression verification exposed a malformed-field
test fixture that assigned a scalar character value to a data-frame column.
R recycled that scalar to all retained rows, so the resulting field was valid
and no error was appropriate. The fixture now supplies an event-length factor
column, which is genuinely malformed for the character-only event provenance
contract. `R/ph3_background.R` now rejects every present optional cutoff
status, reason, or provenance field unless it is a character vector with
exactly one value per retained event. Truly absent historical fields still
receive the existing compatibility defaults only for non-computed methods.
Computed raw-4N cutoff inputs still require every field and their exact call
state. No cutoff, background, membership, or fallback rule changed.

Codex ran no R, Python, package, or analysis command for this correction.
Runtime verification remains **NOT RUN — local execution required**.

## Background-regression compatibility correction (2026-09-01)

The owner-reported affected-regression run found that established SYNTHETIC
pre-cutoff models do not contain the event-level fields introduced for the
raw-4N cutoff (`positivity_call_status`, `positivity_call_reason_code`, and
the cutoff provenance fields). Directly passing their absent/zero-length
values into the background event table caused `data.frame()` construction to
stop before the established background tests could run.

`R/ph3_background.R` now supplies deterministic compatibility values only
when a classification does **not** declare the computed raw-4N method:
`called` for status and scalar missing values for the reason/provenance
fields. Present optional fields must still have exactly one value per retained
event. A classification that declares the computed raw-4N method remains
fail-closed: every new status, reason, and provenance field is required at
event length; missing, null, zero-length, mixed-method, or malformed values
stop analysis. The method identifier must be one nonmissing scalar value
repeated across each retained classification; this prevents a malformed table
from entering the legacy path.

The computed raw-4N path intentionally excludes uncalled/nonfinite raw events
from its computed-negative regression membership. The established legacy path
retains its original membership exactly: a nonfinite negative remains a fit
failure and follows the documented per-replicate-set raw fallback rather than
being silently dropped. Focused SYNTHETIC tests now cover compatibility values,
legacy nonfinite fallback, a noncomputed method identifier, malformed
event-field length, missing computed fields, invalid computed status/reason
semantics, and mixed method identifiers.
No cutoff rule, legacy membership, or production fail-closed validation was
weakened.

This correction was inspected statically only. **NOT RUN — local execution
required** after the correction.

## Verification correction (2026-09-01)

The first user-run focused test did not load because
`R/ph3_raw_density_cutoff.R` had an unmatched closing bracket in the two
adjacent scalar density lookups used for `peak_raw_signal` and `peak_density`.
No analysis, cutoff, validation, or test assertion ran. The correction closes
each existing `[[...]]` lookup and changes no method or validation behavior.
The nearby `cutoff` and `cutoff_density` lookups were inspected and already
had balanced delimiters. Verification remains **NOT RUN — local execution
required** after this source correction.

The failed local-test transcript supplied by the scientific owner recorded the
pre-correction expression exactly as:

```r
peak_raw_signal = fit$x[[peak[[1L]]], peak_density = fit$y[[peak[[1L]]],
```

Because this new source file is untracked, Git cannot independently reconstruct
its pre-correction version. The source-history assertion is therefore limited
to that supplied transcript; the current-source review independently verified
only that the two corrected scalar lookups and the adjacent cutoff lookups are
balanced. This limitation does not alter the required local rerun.

## Scope

This amendment replaces the uncommitted FlowJo-positive-label calibration
proposal with the owner-confirmed pH3 positivity method:

1. Within each explicit biological-replicate set, use its one matched untreated
   control's retained 4N events.
2. On raw linear pH3 fluorescence, fit `stats::density()` with Gaussian kernel,
   `bw = "nrd0"`, `adjust = 1`, and a 2,048-point observed-range grid.
3. Require at least 100 qualifying control events, select the unique
   density-dominant interior local maximum, and use the first right-side local
   minimum as the cutoff.
4. Apply `raw_pH3 > cutoff` unchanged to all retained events in that control
   and its matched treatments, across every DNA region.
5. Fit DNA-dependent background regression from the resulting computed-negative
   eligible 2N–4N events. Corrected signal remains the C/D median basis.

The FlowJo pH3-positive child export remains retained direct-identity
provenance during the transition, but is not the numerical positivity
authority when `ph3_positivity_method: ph3_raw_4n_density_cutoff_v1` is set.

## Files changed

- `R/config.R`: adds the explicit pH3 positivity-method selector.
- `R/ph3_raw_density_cutoff.R`: fixed raw-4N density estimator, per-set
  control resolution, frozen application, cutoff provenance, and application
  QC.
- `R/analysis.R`: applies the computed cutoff before production pH3
  quantitation, background correction, signals, and report-model derivation.
- `R/ph3.R`: reconciles computed positivity against raw signal/cutoff while
  retaining and validating the source FlowJo method as non-authoritative
  provenance.
- `R/ph3_outputs.R` and `R/ph3_background.R`: bind cutoff provenance to the
  output model and background regression; regression negatives now follow the
  computed membership when enabled.
- `tests/testthat/test-ph3-raw-4n-density-cutoff.R`: focused SYNTHETIC checks.

The prior uncommitted FlowJo-label calibration diagnostic source, tests,
manual pages, and public API exposure were intentionally removed because the
owner superseded that route with the confirmed density-cutoff method. Its
dated diagnostic record and external export runbook are retained as historical
provenance only; they are not loadable or runnable production functionality.

## Scientific changes

- Positivity threshold: changed only when the new explicit configuration
  method is selected, from the FlowJo child membership to the confirmed
  matched-control raw-4N density cutoff.
- Background-negative membership: computed-negative eligible 2N–4N events,
  after the raw cutoff has been frozen; this is noncircular.
- Sample mappings, exclusions, DNA gates, DNA normalization, reference rules,
  statistics, EdU, POI, and legacy pH3 behavior: **none**.

Experimental inputs read or modified: **none**.

## Local verification — required

Working directory:

```sh
cd "/Users/sansamc/OMRF Dropbox/Chris Sansam/Synchd/Organized Dropbox/Projects and Data/Projects/Cowork_Projects/2026_Analyze_FACS_Data/facs_pseudocolor_workflow"
```

Prerequisites: R >= 4.2, the package's declared dependencies, and locally
installed `testthat`, `withr`, and `devtools`. Run the commands only against
this source checkout; no package installation, export, real-data analysis, or
Quarto render is authorized by this verification step.

Run locally only:

```sh
Rscript -e 'devtools::test(filter = "ph3-raw-4n-density-cutoff", stop_on_failure = TRUE)'
Rscript -e 'devtools::test(filter = "ph3-output-model|ph3-background-regression|ph3-signal-outcomes|ph3-condition-report-model|ph3-condition-report-plots", stop_on_failure = TRUE)'
Rscript -e 'devtools::test(filter = "^ph3$", stop_on_failure = TRUE)'
Rscript -e 'devtools::test(stop_on_failure = TRUE)'
git diff --check
git status --short --untracked-files=all
```

Expected: zero failures/errors and no output from `git diff --check`. Stop on
any nonunique peak, missing right-side minimum, fewer than 100 control 4N
events, nonfinite selected raw signal, missing/ambiguous control mapping,
changed direct identity, changed DNA membership, unexpected warning, or
unexpected generated file. Do not use a report derived from real data until a
separate direct-identity export/migration review confirms its input contract.

Agent execution: **NOT RUN — local execution required.**

## Review findings and resolution

Independent scientific review confirmed the raw KDE estimator and its
matched-control selection, then identified a blocking contract gap: a failed
cutoff must persist as an explicit unavailable replicate-set record and
no-cutoff report panel. The owner resolved this on 2026-09-01: retain all A–D
rows for that set as `NA` / `unavailable_cutoff_failure`; omit it from the
condition mean; retain valid-clone means with explicit partial-cutoff-failure
coverage; and show a density QC panel with the failure reason and no cutoff
line. The implementation now follows that rule without converting a failure
to zero, negative membership, or a fallback threshold.

The review also identified and this change resolves: legacy exact
positivity-method aggregation rejection, nonfinite raw events entering the
computed-negative background fit, incomplete cutoff provenance fields, and
one-to-one cutoff application-QC reconciliation.

Independent artifact/security review found no network, telemetry, secret,
unsafe I/O, or destructive behavior. It identified two calibration-function
exports that referenced deliberately untracked superseded source; those
exports were removed. It also requested explicit local-test prerequisites and
retained density error text; both are addressed above. Neither reviewer ran
project code or accessed experimental files.

After the owner-approved unavailable-set resolution, fresh independent
scientific-integrity re-review found no remaining blocking finding: uncalled
and nonfinite events are excluded from computed-negative regression, failed
sets are retained unavailable without a surrogate cutoff, and report captions
and no-cutoff QC panels expose partial coverage. Fresh artifact/security
re-review found no remaining blocker: unexpected errors are rethrown, while
only enumerated expected QC failures become unavailable-set records; cutoff
and application provenance bind exactly to active retained events. No reviewer
executed code or accessed experimental inputs.

After the parse correction, an independent scientific-integrity review found
no blocking issue in the current source, but recorded the untracked-file
history limitation above as an important provenance observation. An independent
artifact/security review found no blocking or important issue; it confirmed
that no unsafe I/O, network, telemetry, secret, or hidden fallback was added.
Both reviews were read-only and neither executed project code nor accessed
experimental inputs.

## 2026-09-01 signal-outcomes compatibility correction

The user reported that `ph3-background-regression` passed (79 assertions),
the condition report-model and plot contexts passed (19 and 12 assertions),
and the output-model context passed (46 assertions). The subsequent
`ph3-signal-outcomes` context stopped in six nominal SYNTHETIC cases before
their intended assertions because its pre-computed-cutoff fixture lacked the
two retained call-state columns now required at the signal-outcome boundary.
The signal validator correctly rejected the incomplete event provenance; this
was not evidence of a cutoff, membership, or production-calculation failure.

The smallest correction adds explicit legacy-compatible values to the
SYNTHETIC fixture only: every fixture event has
`positivity_call_status = "called"` and
`positivity_call_reason_code = NA_character_`. The fixture continues to
represent an already-called, non-computed-cutoff path. No production source,
cutoff method, direct identity, event membership, background basis, test
assertion, or failure-closed validation was weakened. Nominal cases can now
reach their C/D assertions, while mutation cases continue to target their
specific validator paths.

Codex ran no project code, tests, package operations, analysis, rendering, or
experimental-data access for this correction. `git diff --check` is the only
permitted static check and remains required before staging. Verification after
this correction is **NOT RUN — local execution required**.

Fresh independent scientific-integrity review found no blocking, important, or
advisory finding: the pre-computed fixture's fixed positive membership is now
explicitly called, while the computed raw-4N cutoff path remains fail-closed.
Fresh independent artifact/security review also found no blocking, important,
or advisory finding: the correction is test/documentation-only and introduces
no I/O, network, telemetry, paths, secrets, artifacts, or hidden fallback.

## 2026-09-01 legacy signal-outcome cutoff-status correction

The next local `ph3-signal-outcomes` run reached the intended C/D calculation
boundary but stopped in nominal legacy SYNTHETIC cases with `missing value
where TRUE/FALSE needed`. The legacy events correctly record
`positivity_call_status = "called"` and intentionally retain
`positivity_call_reason_code = NA_character_`; however, the sample-level
cutoff-failure check used `any(reason_code == "unavailable_cutoff_failure")`.
An ordinary missing reason therefore produced `NA` rather than the required
explicit `FALSE` for a non-computed model.

The correction makes the legacy/non-computed path return explicit `FALSE` for
each acquisition only when every retained event is explicitly called and every
call reason is missing. The cutoff-failure check now matches only the named
failure reason and never treats an ordinary missing reason as evidence of a
failed cutoff. Models that declare the computed raw-4N method must instead
retain the exact cutoff object schema, method settings, one cutoff record per
replicate set, one application record per acquisition, matching cutoff status,
and matching event call status; absent or malformed computed-cutoff provenance
stops fail closed. A focused SYNTHETIC assertion covers both the legacy `FALSE`
state and rejection of a computed method without its cutoff provenance.

No cutoff estimator, threshold, event membership, background-regression
selection, C/D formula, sample mapping, exclusion, normalization, statistic,
plot, report, experimental input, or legacy numerical result was changed.
Codex did not run project code, tests, package operations, rendering, or data
access. Verification after this correction is **NOT RUN — local execution
required**; only `git diff --check` may be used as a static check.
