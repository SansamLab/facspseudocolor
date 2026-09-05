# pH3 legacy CSV pilot import — 2026-09-01

## Log-coordinate pseudocolor density correction — 2026-09-02

User visual QC identified a sharp horizontal pseudocolor-density transition
across the restored negative pH3 band. The panel already displayed restored
pH3 on a log10 y-axis, but its two-dimensional KDE had been calculated in
linear displayed-fluorescence units. Applying the log transformation only
after that KDE warps its color field and can create this display artifact.

The panel now calculates KDE density from `dna_norm` and
`log10(displayed_signal)`, after the established finite, positive display
domain filter, while still plotting `displayed_signal` under the same log10
y-axis. This is display-only: event values, correction predictions and
offsets, cutoffs, positivity, A/B prevalence, C/D medians, summaries,
provenance, and sample mapping are unchanged. No fallback or value
substitution is introduced; events at or below zero remain excluded only from
the pre-existing log-display layer and remain recorded in its QC.

The focused unmistakably `SYNTHETIC` test independently recomputes each
panel's density in its log-coordinate display system and requires an exact
match. Codex did not run tests, renders, package checks, or access
experimental/generated inputs. Verification is **NOT RUN — local execution
required**. After independent reviews, run:

```sh
Rscript -e 'devtools::test(filter = "ph3-legacy-csv-pilot", stop_on_failure = TRUE)'
```

Expected: PASS with no warnings. Then rerender only into a new external
temporary directory using the existing pilot configuration.

## Owner-approved log10 corrected-display refinement — 2026-09-02

## Log10 analytic-region shading warning correction — 2026-09-02

User visual QC of the rendered legacy pilot showed `NaNs produced` warnings
while the corrected, offset-restored pseudocolor panels used a log10 y-axis.
The source was display-only: the analytic-region `geom_rect()` layer used
`ymin = -Inf` and `ymax = Inf`, and log10 has no representation for the lower
infinite bound.  The region rectangles now receive the finite, strictly
positive, ordered y limits derived from the same already positive-domain
display events used by the panel.  If those limits cannot be established, the
pilot stops visibly with `invalid_log_display_y_limits`; it does not clamp,
substitute, or change events, background correction, cutoffs, positivity,
summaries, or provenance.

The focused unmistakably `SYNTHETIC` assertion requires the rectangle layer to
contain only finite positive limits with `ymax > ymin`.  No experimental input
or generated report was accessed or modified by Codex.  Verification is **NOT
RUN — local execution required**. After independent reviews, run:

```sh
Rscript -e 'devtools::test(filter = "ph3-legacy-csv-pilot", stop_on_failure = TRUE)'
```

Expected: PASS with no warnings. Then rerender only into a new external
temporary directory using the existing pilot configuration.

### Approved display-only contract

After the existing per-sample display offset is added to the
background-subtracted pH3 signal, the legacy-pilot pseudocolor panels use a
log10 y-axis. This is a visual presentation change only: raw values,
background predictions, background-subtracted values, provisional and final
cutoffs, positivity calls, A/B prevalence, C/D medians, and all summaries are
unchanged.

Because log10 has no representation for values at or below zero, the pilot
does not silently transform, clamp, substitute, or omit such values from the
analysis. Per sample, only finite events in the existing 1.6N--4.4N visual
window are considered for the panel. Restored values greater than zero are
used for plotted density/points; restored values at or below zero are excluded
from that visual layer only. The new `log_display_qc` records the finite visual
window count, positive-domain count, nonpositive excluded count, status, and
reason. The panel subtitle explicitly flags exclusions when they occur.

If no positive-domain event remains, the report produces an explicit
`LOG10 DISPLAY UNAVAILABLE` no-data panel with no point density and no cutoff
overlay. If a restored cutoff is nonpositive, the available event panel is
retained but displays an explicit no-cutoff visual state. Neither case falls
back to raw pH3 or changes any computed analysis result.

The display label is now `Background-subtracted pH3 with display offset; log10
display`, accurately distinguishing this pilot from raw pH3 and from the
unshifted corrected analytical scale.

### Scope and verification

Files modified: `R/ph3_legacy_pilot.R`,
`inst/quarto/facs_ph3_legacy_pilot.qmd`,
`tests/testthat/test-ph3-legacy-csv-pilot.R`, and this record. Experimental
inputs and generated reports were not accessed or modified. New tests use only
unmistakably `SYNTHETIC` in-memory values and cover the log scale, per-sample
nonpositive-display QC, no mutation of analytical outputs, explicit zero
positive-domain placeholders, and no cutoff overlay in that no-data state.

Codex did not run R, Python, package loading, tests, builds, checks, or
renders. Verification is **NOT RUN — local execution required**. After the
required independent scientific-integrity and artifact/security reviews, run
the focused pilot test from the canonical repository:

```sh
Rscript -e 'devtools::test(filter = "ph3-legacy-csv-pilot", stop_on_failure = TRUE)'
```

Expected: PASS with no warnings. Then rerender only in a new external temporary
directory using the existing pilot configuration. The report remains
`PILOT / LIMITED-PROVENANCE — NOT PUBLICATION-GRADE`.

## Current-method notice — 2026-09-02

The raw-only pilot and raw-pseudocolor descriptions retained below are dated
historical chronology. They are superseded for the current pilot by the
background-subtracted final-cutoff method in the next section; they do not
describe the current calculation, report panels, or C/D outcomes.

## Background-subtracted final-cutoff upgrade — 2026-09-02

## Display-only corrected-pseudocolor offset restoration — 2026-09-02

### Owner-approved display contract

For every sample with an available event-level correction and clone-matched
final corrected cutoff, the pseudocolor display now uses the same eligible
final pH3-negative events (`eligible` DNA region and corrected pH3 at or below
that clone's final cutoff) to calculate:

`display_offset = median(raw pH3) - median(background-subtracted pH3)`.

Only the displayed corrected event signal and displayed red final-cutoff line
receive this constant translation. This restores the displayed final-negative
median to its raw-signal median, without power-of-ten rounding or a fixed
target. The display labels explicitly say that pH3 is background-subtracted
with a display offset restored.

The offset cannot alter the background regression, provisional or final cutoff
determination, positivity, eligibility, A/B prevalence, C/D medians, condition
summaries, or any stored analytical value. Per-sample QC records the exact
final-negative count, raw and corrected medians, offset, displayed cutoff,
status, and reason. If correction, final cutoff, or the valid final-negative
subset is unavailable or nonfinite, the report uses a labelled unavailable
placeholder and does not substitute raw signal or invent an offset.

Files changed for this display-only implementation: `R/ph3_legacy_pilot.R`,
`inst/quarto/facs_ph3_legacy_pilot.qmd`,
`tests/testthat/test-ph3-legacy-csv-pilot.R`, and this record. Experimental
inputs were not accessed or modified. Codex did not run R, Python, tests,
builds, checks, or renders; verification is **NOT RUN — local execution
required**.

After independent reviews, run:

```sh
Rscript -e 'devtools::test(filter = "ph3-legacy-csv-pilot", stop_on_failure = TRUE)'
```

Expected: PASS with no warnings. Then rerender only in an external temporary
directory using the existing pilot configuration. No new repository artifact
is expected.

## Clone-shared display-negative target — 2026-09-02

The prior per-sample raw-negative restoration target is superseded for the
legacy pilot display. For each explicit configured clone/biological-replicate
group, the display target is the R `median()` of that group's **per-sample raw
pH3 medians** among eligible final pH3-negative events. Each sample therefore
contributes exactly one value, rather than contributing all of its events; for
an even number of samples, R's deterministic median convention (the arithmetic
mean of the two middle sample medians) applies.

For every sample in that complete configured group:

`display_offset = group_shared_raw_negative_median - sample_corrected_negative_median`.

This makes the displayed corrected final-negative medians agree within a clone
group. The offset and shifted red cutoff remain display-only. Background
regression, provisional/final cutoffs, positivity, A--D outcomes, condition
summaries, mappings, and stored analytical values are unchanged.

There is no partial-group display fallback. Every configured sample in a clone
group must have an available, finite final-negative raw and corrected median.
Otherwise every panel in that group is an unavailable placeholder, with
`shared_target_status = unavailable`,
`shared_target_reason_code = configured_group_member_display_input_unavailable`,
and `display_offset_reason_code = clone_group_shared_raw_negative_median_unavailable`.
The display-offset QC and pseudocolor metadata retain the shared target and its
status/reason for provenance.

This implementation uses only clearly labelled SYNTHETIC tests. Codex did not
run project code, tests, builds, checks, or renders; verification is **NOT RUN
— local execution required**.

### Status

Implemented for local user verification only. This owner-approved update
changes the pilot scientific method. It does not change the explicit Cdk1as
six-sample mapping, historical input locations, `flowjo.rebuild: false`, or
the strict `production_direct_identity_v1` workflow.

### Approved deterministic method implemented

For every pilot sample independently, the existing raw 4N KDE
peak-to-right-valley rule supplies a **provisional** raw-signal split. That
split is an audit intermediate only: its provisional-negative events fit a
linear pH3-versus-normalized-DNA background model for the same sample. The
predicted background is subtracted event-by-event. For each clone, the final
KDE cutoff is then derived only from the **background-subtracted 4N events of
that clone's Untreated control**, and applied only to that clone's Untreated
and matched NMPP1 treatment. Final A/B positivity, C/D positive-event medians,
and pseudocolor panels all use the corrected scale. The report retains the
raw provisional cutoff and regression coefficients/counts as separately named
QC provenance, so they cannot be mistaken for final calls.

The pilot fails closed if a complete-event value is missing/nonfinite, a
per-sample provisional raw cutoff fails, a background fit is invalid, or the
final corrected-control cutoff fails. Any such sample correction failure makes
all outcomes for the affected clone unavailable and excludes those unavailable
values from condition means. No failed clone can contribute a numerical zero.

### Prefix-alignment correction — 2026-09-02

Before the pilot accesses any complete-event table, the names of
`analysis$normalized_data` must exactly match the configured manifest prefixes
in manifest order. A reordered, missing, renamed, or unnamed table list stops
with `normalized_data_prefix_alignment`; list position is never used as an
implicit sample-to-condition or sample-to-clone mapping. A clearly
`SYNTHETIC` reordered-list regression covers this fail-closed boundary.

This correction changes neither the explicit mapping, cutoff or background
method, normalization, report calculations, nor experimental inputs.

### Scope and verification

Code changed: `R/ph3_legacy_pilot.R`,
`inst/quarto/facs_ph3_legacy_pilot.qmd`,
`man/plot_ph3_legacy_pilot_report.Rd`, and
`tests/testthat/test-ph3-legacy-csv-pilot.R`.
Experimental inputs accessed or modified by Codex: none. Synthetic tests only
use unmistakably `SYNTHETIC` in-memory data. Tests are **NOT RUN — local
execution required**. Codex performed no R/Python execution, rendering,
package build/check, or real-data access.

Required local verification from the canonical repository:

```sh
Rscript -e 'devtools::test(filter = "ph3-legacy-csv-pilot", stop_on_failure = TRUE)'
Rscript -e 'devtools::test(filter = "^ph3$", stop_on_failure = TRUE)'
Rscript -e 'devtools::test(stop_on_failure = TRUE)'
git diff --check
git status --short --untracked-files=all
```

Expected: zero failures/errors and no unexpected generated repository files.
After user verification, render only into an external temporary directory. The
report remains `PILOT / LIMITED-PROVENANCE — NOT PUBLICATION-GRADE`: event
identity and FlowJo gate containment are unverified, and no verified geometry
is claimed.

### Corrected-scale SYNTHETIC fixture correction — 2026-09-02

The first local run after the background-subtracted upgrade reported ten
focused-test failures.  The production pilot correctly treated every affected
sample correction as unavailable; the error was the intended-valid
`SYNTHETIC` fixture.  Its treatment arrays had two equal-height point-mass
peaks, and its control arrays were broad uniform ramps.  Those shapes do not
guarantee the approved KDE contract's required unique dominant peak followed
by a right-side local minimum.  Consequently provisional cutoffs, background
fits, corrected final cutoffs, corrected panels, and their downstream
assertions were unavailable.

The fixture now uses deterministic, unmistakably `SYNTHETIC` 160-event
negative and 40-event positive clusters.  Each has a unique dominant negative
peak, a separated positive peak, and repeated finite DNA coordinates spanning
both provisional membership groups.  This supplies a valid sample-level
background regression and a valid corrected-control cutoff without weakening
the KDE, coverage, correction, clone-matching, or fail-closed rules.  The
existing zero-range cutoff-failure and invalid-background-fit cases remain
separate genuine failure cases.

Production code and experimental inputs were not changed or accessed.  The
corrected tree is **NOT RUN — local execution required**.  After independent
scientific-integrity and artifact/security review, run:

```sh
Rscript -e 'devtools::test(filter = "ph3-legacy-csv-pilot", stop_on_failure = TRUE)'
```

Expected: PASS with no warnings.  Stop and report the failed test and a short
error excerpt if it does not pass.

## Status

Implemented for local user verification only. This path is deliberately
separate from `production_direct_identity_v1` and is not publication-grade.

## Raw pH3-versus-DNA pseudocolor QC addition — 2026-09-02

### Compact pseudocolor display layout — 2026-09-02

The pilot report now renders each raw pH3-versus-DNA QC panel at a compact
2.5 x 2.7 inch outer figure size. Its data viewport is forced to a square
aspect ratio, making the plotting area approximately 1 x 1 inch after the
title, subtitle, axes, ticks, and provenance caption. This is a display-only
layout change: the events, raw signal values, density calculation, cutoff,
DNA regions, sample mapping, and all numerical outcomes are unchanged. The
panel continues to show its raw/non-background-subtracted and
limited-provenance labels, matched-control cutoff, and region overlays.

Focused synthetic coverage asserts that every pilot pseudocolor panel retains
the square data-viewport contract. Verification is **NOT RUN — local execution
required**; use the existing pilot focused test followed by the external
temporary-directory render command below. No experimental inputs were read or
modified by Codex.

### Display-bound correction — 2026-09-02

The first local execution of the new `SYNTHETIC` pseudocolor-panel tests
stopped before panel construction because `coord_cartesian()` was mistakenly
given all three display-region boundaries (2N lower bound, 4N lower bound,
and 4N upper bound) as its x limits. This was a display-only implementation
error; no cutoff, positivity call, DNA-region membership, mapping, summary,
or input handling had run or changed. The panel now keeps the three values
for region shading and boundary lines, while supplying only the two outer
2N--4N endpoints to `coord_cartesian()`. A focused, unmistakably `SYNTHETIC`
assertion requires every constructed panel to retain exactly those two limits.

Verification after this correction is **NOT RUN — local execution required**:

```sh
Rscript -e 'devtools::test(filter = "ph3-legacy-csv-pilot", stop_on_failure = TRUE)'
```

Expected: PASS with no warnings. No experimental inputs were accessed or
modified. Independent scientific-integrity and artifact/security reviews of
this correction are required before real-data pilot rendering resumes.

The legacy pilot report now includes one deterministic, display-only raw
pH3-versus-normalized-DNA pseudocolor panel for every explicitly configured
pilot sample (six for the Cdk1as clone 6/7/15 × Untreated/2 h NMPP1 config).
Each panel draws only from that sample's already loaded complete-event table,
uses a 2D KDE of normalized DNA versus log10 raw pH3 for visual density, and
retains raw pH3 values on a log-scaled display. The matching clone's raw 4N
density cutoff is overlaid when available. The 2N--4N analysis bounds, the
below-4N region, and the 4N region are visibly distinguished.

The panels are QC only: they do not change the pilot cutoff, positivity calls,
summaries, sample mapping, or any production direct-identity behavior. They
do not reconstruct or claim verified FlowJo geometry. Every panel retains the
`PILOT / LIMITED-PROVENANCE — NOT PUBLICATION-GRADE` label and explicitly
states that raw pH3 is not background-subtracted. If a cutoff is unavailable,
the panel shows its explicit reason and omits the horizontal cutoff line.

Implementation files: `R/ph3_legacy_pilot.R`,
`inst/quarto/facs_ph3_legacy_pilot.qmd`,
`tests/testthat/test-ph3-legacy-csv-pilot.R`, and
`man/plot_ph3_legacy_pilot_report.Rd`. Experimental inputs read or modified
by Codex: none. New focused coverage is entirely in-memory and unmistakably
`SYNTHETIC`; it verifies a unique panel and matched-control-cutoff provenance
for every configured sample, the raw/limited-provenance labels, and a cutoff
overlay when the cutoff is available.

Verification is **NOT RUN — local execution required**. From the canonical
repository, run:

```sh
Rscript -e 'devtools::test(filter = "ph3-legacy-csv-pilot", stop_on_failure = TRUE)'
Rscript -e 'devtools::test(filter = "^ph3$", stop_on_failure = TRUE)'
Rscript -e 'devtools::test(stop_on_failure = TRUE)'
git diff --check
```

Expected: all tests pass with no unexpected warnings and `git diff --check`
is silent. Only then rerender the pilot in an external temporary directory
using the command in **Local verification** below. Expected render result: one
clearly labelled HTML report containing six raw pH3-versus-DNA QC panels; no
files are generated inside the repository or historical experiment directory.

Independent scientific-integrity review found no blocking or important issue.
Its one advisory—explicitly test the unavailable-cutoff display branch—was
resolved with an input-free `SYNTHETIC` test that requires the reason text and
the absence of a horizontal cutoff overlay. Independent artifact/security
review was clean: no external I/O, network/telemetry, hidden fallback,
machine-specific input selection, generated artifact, or unrelated scope was
introduced. Codex performed only static review and `git diff --check`; tests
and rendering remain **NOT RUN — local execution required**.

## Verification correction

The first user run of `test-ph3-legacy-csv-pilot.R` stopped in its input-list
assertion because its in-memory SYNTHETIC configuration used a relative
`data_dir` without the `config_dir` attribute supplied by file-backed
configuration. This was a test-fixture invocation error, not a pilot-input or
production-validation failure. The assertion now supplies `tempdir()` as the
explicit `data_dir` override solely while inventorying the expected synthetic
population names. The fail-closed production rule remains unchanged: a
relative configured `data_dir` requires a configuration read from a file (or
an explicit caller override).

Corrected-tree verification is **NOT RUN — local execution required**.

## User verification correction — 2026-09-02

The user ran the focused pilot test after the cutoff-failure correction.  The
failed clone itself correctly produced unavailable rows; the remaining five
test failures were in the synthetic summary assertions.  Those assertions
incorrectly compared Untreated and treatment condition means to one pooled
cross-condition mean, and required two finite values even for an intentionally
empty compartment.  The test now compares each condition/outcome summary only
with finite surviving values from that same condition and outcome, and requires
an explicit unavailable summary when no surviving values qualify.  Production
pilot code, cutoff application, mappings, and real inputs are unchanged.

Corrected-tree verification is **NOT RUN — local execution required**.

## Clone-specific cutoff regression coverage

A focused, input-free `SYNTHETIC` regression now gives each configured clone a
distinct Untreated 4N raw-pH3 distribution and confirms that every row for the
clone, including its matched treatment, carries only that clone's cutoff.  It
also uses a clone-1 treatment signal that is positive under clone 1's cutoff
but below clone 2's cutoff, making cross-clone control substitution observable.
This is test coverage only; it does not alter the pilot method, mappings, or
any real input path.

## Cutoff-failure containment correction

A static review of the clone-specific regression found that a failed matched
control cutoff could previously pass a recycled `FALSE` positivity vector to
the A/B percentage helper. That could produce a numerical zero and allow it
into a condition mean. The pilot now sets all four A--D values to `NA` for
every sample in that clone, with `value_status = "unavailable"` and
`reason_code = "unavailable_cutoff_failure"`. The corresponding cutoff field
is `NA`, and only finite values from clones with available cutoffs can enter a
condition summary. A clearly `SYNTHETIC` zero-range-control regression checks
all of these properties while retaining the distinct-cutoff clone-isolation
test.

The pilot now also accepts unavailable results only for the seven explicit
raw-density cutoff failure codes. Any other error (including an internal
programming failure) is rethrown, rather than being converted into a partial
summary. A `SYNTHETIC` mocked internal-error regression verifies that
fail-closed behavior.

This is a fail-closed output correction only: the raw-density method, explicit
six-sample mapping, matched-control rule, and all valid-cutoff calculations
are unchanged. Corrected-tree verification is **NOT RUN — local execution
required**.

## Scope

The Cdk1as pilot config lives outside the repository at
`test_reports/ph3_cdk1as_legacy_csv_pilot_2026-09-01/config.yml`. It reads the
six existing historical `Single Cells` and `G1` CSV exports, with explicit
clone 6/7/15 × Untreated/2 h NMPP1 mappings. It never rebuilds FlowJo exports,
writes into the historical experiment, reads FCS/WSP files, or uses the
historical pH3-positive export as a numerical input.

Each clone's Untreated 4N raw-pH3 density determines its frozen local-minimum
cutoff. That cutoff is applied to its own Untreated and NMPP1 samples. A
cutoff failure produces unavailable rows rather than a guessed value. The
four pilot panels show A/B prevalence and C/D *raw* signal medians; C/D are
explicitly not background-subtracted.

## Limitations

Legacy CSV inputs have neither direct event identity nor verifiable population
containment. Every input warning, result object, plot caption, and Quarto
report labels this as `PILOT / LIMITED-PROVENANCE — NOT PUBLICATION-GRADE`.
The strict direct-identity profile remains unchanged and is required for final
analyses.

## Local verification

Working directory: the canonical repository. Prerequisites: declared R
dependencies, `testthat`, `devtools`, and Quarto for report rendering.

```sh
Rscript -e 'devtools::test(filter = "ph3-legacy-csv-pilot", stop_on_failure = TRUE)'
Rscript -e 'devtools::test(filter = "^ph3$", stop_on_failure = TRUE)'
Rscript -e 'devtools::test(stop_on_failure = TRUE)'
git diff --check
```

Expected: no test failures/errors or unexpected warnings, and no output files
inside the repository.

To create a pilot report only after the tests pass, use an external temporary
directory:

```sh
PH3_REPO="$PWD"
# Set this to the explicit absolute path of the Cdk1as pilot config supplied
# with this project under test_reports/ph3_cdk1as_legacy_csv_pilot_2026-09-01/.
PH3_PILOT_CONFIG="/absolute/path/to/ph3_cdk1as_legacy_csv_pilot_2026-09-01/config.yml"
PH3_PILOT_DIR="$(mktemp -d)"
mkdir -p "$PH3_PILOT_DIR/Rlib"
(
  cd "$PH3_PILOT_DIR" &&
  R CMD build "$PH3_REPO" &&
  R CMD INSTALL -l "$PH3_PILOT_DIR/Rlib" facspseudocolor_*.tar.gz &&
  R_LIBS="$PH3_PILOT_DIR/Rlib" quarto render \
    "$PH3_PILOT_DIR/Rlib/facspseudocolor/quarto/facs_ph3_legacy_pilot.qmd" \
    --output-dir "$PH3_PILOT_DIR/report" -P config="$PH3_PILOT_CONFIG"
)
```

Expected: exactly one pilot-labelled HTML report below the external temporary
directory. Stop on missing files, a channel mismatch, a cutoff failure you do
not want represented as unavailable, or any attempt to rebuild FlowJo exports.

## Provenance

Experimental inputs read or modified by Codex: none. Scientific methods
changed: none in the verified production profile. This pilot is an explicitly
limited real-data import path, not an upgrade of legacy data provenance.

## Display-only pseudocolor refinement (2026-09-02)

The six raw pH3-versus-DNA pilot QC panels use a square data viewport intended
to be approximately 1.25 x 1.25 inches, with a compact relative-density legend.
Their visual DNA window is deterministically 1.6N--4.4N, calculated as
`0.8 * dna_2n_value` through `2.2 * dna_2n_value`. The analytic 2N--4N
boundaries and below-4N/4N shading remain visible but are not changed. This is
strictly display-only: event filtering for metrics, cutoff calibration,
positivity, sample mappings, and summaries are unchanged.

Verification is **NOT RUN — local execution required**. Run the focused legacy
pilot test, then render the report in an external temporary directory using
the commands in the Local verification section above.

## Source-load correction (2026-09-02)

User-local verification reported that `R/ph3_legacy_pilot.R` could not be
parsed at the assignment which stores an available sample's corrected events
under its manifest prefix. The expression had one closing bracket too few in
the nested `corrections[[manifest$prefix[[row]]]]` index. R also requires the
assignment used as an `if` body to be braced. The correction adds only the
missing bracket and braces. The stored key remains exactly the manifest prefix
for the same row; no background fitting, corrected cutoff, sample alignment,
mapping, or report logic changed. Codex did not run R or project code. Verification is
**NOT RUN — local execution required**; run the focused pilot test below.

## Corrected-cutoff provenance assertion correction (2026-09-02)

The user-local focused pilot test then reported one failure in the clearly
`SYNTHETIC` clone-matched cutoff test: it required corrected cutoffs to be
strictly increasing in synthetic clone order. That ordering is not part of the
approved method. Per-sample background regression can validly change the
relative scale and ordering of clone-specific corrected cutoffs.

The test no longer asserts cross-clone monotonicity. It instead verifies the
auditable contract directly: every cutoff record resolves to an existing
manifest prefix, that prefix belongs to the same biological replicate as the
cutoff record, it is the configured `Untreated` control, and every Untreated
and matched-treatment result row for the replicate carries exactly that
replicate's corrected cutoff. No production code, threshold rule, background
subtraction, sample mapping, or failure handling changed.

Experimental inputs accessed or modified by Codex: none. This corrected tree
is **NOT RUN — local execution required**. After independent scientific-
integrity and artifact/security review, run:

```sh
Rscript -e 'devtools::test(filter = "ph3-legacy-csv-pilot", stop_on_failure = TRUE)'
```

Expected: PASS with no warnings. Stop and report the failed test name and a
short error excerpt if it fails.

## Unavailable-correction report containment (2026-09-02)

The user-local Cdk1as pilot render stopped at `clone7_nt` with
`invalid_display_correction`: that sample did not have an available event-level
background correction, and the report previously required every pseudocolor
panel to be drawable from corrected events. This behavior was too broad: an
unavailable correction must make its clone's numerical outcomes unavailable,
but must not prevent visually reviewing other samples.

The report now uses the stored per-sample `provisional_background_qc` row as
the authoritative display-status source. For a sample whose correction status
is `unavailable`, it renders a labelled placeholder containing the explicit
stored status and reason, the limited-provenance warning, and a statement that
no raw events, corrected density, or corrected cutoff are shown. It does not
access that sample's event table to construct a substituted panel. Samples
with an available correction retain the existing corrected-scale pseudocolor
data and clone-matched cutoff behavior. The existing QC table remains in the
report, and condition summaries retain their explicit partial/unavailable
coverage statuses.

An unmistakably `SYNTHETIC` regression invalidates one sample-level background
fit, requires the placeholder and its no-fabrication text, requires no point
or cutoff layer in that panel, confirms a separate valid panel remains a
corrected event display, and confirms unavailable-clone values cannot enter
condition means. No background method, cutoff method, sample mapping, strict
profile, experimental input, or numerical summary formula changed.

Codex did not run R, Python, builds, checks, rendering, or access experimental
data. Verification is **NOT RUN — local execution required**. After required
independent reviews, run the focused pilot test and then render only into an
external temporary directory. Expected render behavior: a report is produced
when some samples have unavailable corrections; those panels are visible
placeholders rather than raw or fabricated corrected plots.

## Display-offset horizontal-line assertion correction (2026-09-02)

The user-local focused pilot test reported six failures in the new
`SYNTHETIC` display-offset assertion: it looked for the final cutoff line's
`yintercept` in `aes_params`, where `ggplot2::geom_hline()` does not retain an
explicitly supplied intercept. Static inspection confirms that the production
panel supplies `metadata_row$displayed_cutoff_signal` directly to
`geom_hline(yintercept = ...)`; `geom_hline()` stores this fixed value in its
one-row layer data.

The correction is test-only. It requires exactly one horizontal-line layer
per available panel and compares that layer's sole stored `yintercept` exactly
with the corresponding `displayed_cutoff_signal` QC value. The display offset
therefore remains explicitly verified while staying isolated from the
background-corrected cutoff, positivity calls, summaries, and all other
analysis values. Production plotting, background subtraction, cutoff
calibration, sample mapping, and experimental inputs are unchanged.

Codex did not run project code, tests, builds, renders, or analyses, and did
not access experimental or generated inputs. Verification is **NOT RUN —
local execution required**. After independent scientific-integrity and
artifact/security review, run:

```sh
Rscript -e 'devtools::test(filter = "ph3-legacy-csv-pilot", stop_on_failure = TRUE)'
```

Expected: PASS with no warnings. Stop and report the failed test name and a
short relevant excerpt if it fails.

## Independent review resolution — 2026-09-05

Two independent, read-only reviews of the full legacy-CSV-pilot change set
were run: a scientific-integrity review and an artifact/security review,
neither authored by the pilot's implementer and independent of each other.

The artifact/security review found no blocking issue across all eight
checked categories: network calls/telemetry, secrets, unsafe execution,
destructive file operations, hardcoded machine-specific paths, hidden
fallbacks, working-tree cleanliness, and render containment to an external
temp directory.

The scientific-integrity review confirmed the Clone 7 fail-closed/exclusion
cascade, every "display-only" claim in this log against the actual code, the
bracket/brace correction, the `geom_hline` test fix, the removed cross-clone
monotonicity assertion, and the unavailable-correction placeholder — all
matching this log's descriptions. It raised one blocking finding: the
pilot's final positivity cutoff is computed from the matched Untreated
control's background-subtracted 4N events rather than the raw-scale cutoff
Rules 2 and 6 of
`docs/implementation/ph3_raw_4n_density_cutoff_contract_amendment_2026-09-01.md`
require for the frozen `ph3_raw_4n_density_cutoff_v1` method, and no
separate, dated, owner-confirmed record for this pilot-specific deviation
existed independent of this log's own prose.

The scientific owner has since directly confirmed that this deviation is
intended. That confirmation is recorded as a dated addendum in
`ph3_raw_4n_density_cutoff_contract_amendment_2026-09-01.md`
("2026-09-05 scope clarification — legacy CSV pilot cutoff exception
(owner-confirmed)"), scoped exclusively to the `legacy_csv_pilot_v1` profile.
It does not change any production code, cutoff calculation, background
method, sample mapping, or test; it closes the review finding with an
owner-confirmed decision record separate from the implementer's own log.

No production code changed as part of this entry, so no test re-run was
required. `git diff --check` remains silent after this addition.
