# Standard EdU report redesign implementation record

**Date:** 2026-09-09
**Verification:** focused report-contract test PASS and corrected real-data
render PASS; broader package verification remains local execution required

## Canonical source and scope

The canonical source is the nested `facs_pseudocolor_workflow` repository,
implemented in the isolated `facs_pseudocolor_workflow-improve-standard-edu-report`
worktree on branch `feature/improve-standard-edu-report` from `origin/main`.
The redesigned standard report is
`inst/quarto/facs_edu_pseudocolor_output_contract.qmd`; its established
scientific report model remains `R/edu_pseudocolor_output_contract.R`.

## Prioritized bench-scientist usability audit

1. **Blocking orientation gap:** the report opened directly on plots without a
   concise experiment, sample, replicate, channel, population, normalization,
   warning, or provenance summary.
2. **Important QC visibility gap:** suppressed panels, non-finite input counts,
   retained warnings, zero denominators, and unavailable reference ratios were
   stored in artifacts but not presented before results.
3. **Important interpretation gap:** plot colors, the display-only EdU offset,
   denominators, technical-acquisition aggregation, and reference-relative
   intensity required prior knowledge.
4. **Important navigation gap:** the short report did not follow the scientist's
   reading sequence or distinguish primary results, QC, methods, and provenance.
5. **Superseded historical export finding:** canonical CSV export existed in
   the package but was not exposed by the original report configuration. An
   interim implementation exposed it; the owner later simplified the standard
   report to an HTML-only contract and removed report-triggered CSV export.
6. **Advisory accessibility gap:** minimal whitespace, generic captions, and
   limited visual hierarchy made rapid scanning difficult. The historical note
   about absent quantitative tables is superseded: interim tables were added,
   then explicitly removed during owner-directed report simplification.

No existing rendered standard EdU HTML or image was present. The gallery images
are POI examples and were inspected only as repository-generated-output
inventory; they were not edited or used as experimental evidence.

## Implementation

The current simplified report has a results-first visual overview, a prominent
pre-results QC status, ordered figures, per-card captions, methods, and detailed
provenance. Historical iterations included full quantitative tables, expanded
QC, comparison matrices, replicate-variation tables, and report-triggered CSV
export; those reader-facing sections are superseded and no longer part of the
standard HTML contract.
Reusable helpers build the overview and QC views from the completed analysis
without altering it. The report no longer accepts `output_csv_dir`, performs
CSV export, or writes quantitative-table files during rendering. Canonical CSV
export remains available through the pre-existing package workflow outside
this report. Full input paths are hidden
by default; explicitly enabling audit tables displays them with a sharing
warning. Structural/report-contract tests
cover section order, configuration controls, identity/provenance visibility,
QC surfacing, and immutability. The focused test explicitly
loads the shared report setup file from `inst/quarto`, matching the report's
helper-loading boundary without moving report-only helpers into the package
namespace.

In a superseded intermediate iteration, the report placed a side-by-side
comparison section before the every-sample audit view. It constructed one
configured-reference-versus-condition matrix per non-reference condition,
with biological replicates as rows. Membership is resolved only from explicit
manifest fields. Every technical acquisition is retained and stacked inside
its condition–replicate cell; missing combinations are rendered as explicit
`MISSING` panels. This is a second presentation of the existing acquisition
plots and performs no aggregation or recalculation.

## Scientific integrity and provenance

Experimental inputs read or modified: **none**. Repository fixture CSVs were
listed as inventory but not opened, analyzed, or modified. Original FCS, WSP,
CSV, Prism, and other source inputs remain untouched.

Sample mappings, exclusions, gates, thresholds, phase boundaries, channel
roles, normalization, background subtraction, denominators, statistics, and
biological interpretations changed: **none**.

The report continues to consume the established FlowJo-derived Single Cells,
G1, and EdU Positive inputs and the validated canonical report tables. It adds
no distribution heuristic or event-count cutoff: event counts are displayed,
and existing warnings/failure states are surfaced. Suspicious distributions
still require scientific visual review of the every-sample panels.

## Files

Modified:

- `R/edu_pseudocolor_output_contract.R` (presentation label only)
- `inst/quarto/facs_edu_pseudocolor_output_contract.qmd`
- `inst/quarto/_report-setup.R`

Added:

- `tests/testthat/test-edu-report-contract.R`
- `docs/implementation/standard_edu_report_redesign_2026-09-09.md`

Generated outputs: **none**.

## Verification required

No R, Python, package, test, analysis, or Quarto command was executed because
the owner's temporary execution authorization expired. From the isolated
worktree, with R package development dependencies and Quarto installed, first
install this modified checkout into an owner-created library outside the
repository, then force every check to use that library:

```sh
cd "<path-to-facspseudocolor-worktree>"
mkdir -p "/absolute/path/outside/repository/edu-verification/library"
R CMD INSTALL --library="/absolute/path/outside/repository/edu-verification/library" .
R_LIBS="/absolute/path/outside/repository/edu-verification/library" Rscript -e 'devtools::test(filter = "edu-report-contract", stop_on_failure = TRUE)'
R_LIBS="/absolute/path/outside/repository/edu-verification/library" Rscript -e 'devtools::test(filter = "edu-pseudocolor-output-contract", stop_on_failure = TRUE)'
R_LIBS="/absolute/path/outside/repository/edu-verification/library" Rscript -e 'devtools::test(filter = "edu-output-contract", stop_on_failure = TRUE)'
R_LIBS="/absolute/path/outside/repository/edu-verification/library" Rscript -e 'devtools::test(stop_on_failure = TRUE)'
```

Expected: package installation and tests complete with zero failures, errors,
or unexpected warnings and no new repository files. The installed package is
created only under the external verification directory. Stop on any changed
canonical table, sample/manifest mismatch, suppressed
panel without an explicit QC row, missing input, unavailable normalization,
unexpected warning, or generated repository artifact.

For a visual check, render only with an owner-selected completed EdU analysis
RDS and put the HTML outside the repository:

```sh
R_LIBS="/absolute/path/outside/repository/edu-verification/library" quarto render inst/quarto/facs_edu_pseudocolor_output_contract.qmd --output-dir "/absolute/path/outside/repository/edu-report-check" -P analysis_rds:"/absolute/path/to/owner-selected-analysis.rds"
```

Expected: one self-contained HTML report in the requested external directory;
no CSV tables or other report exports are expected. Stop if the report reads
experimental CSVs unexpectedly, writes inside the repository, changes any
scientific value, or shows a QC failure that is not understood.

## Assumptions and unresolved issues

- The completed analysis is the scientific source of truth and contains the
  canonical quantitation required by the existing output-contract builder.
- No new automated definition of “low event count” or “suspicious distribution”
  was introduced without owner approval. Counts and existing warnings are made
  prominent; distributions remain a visual QC responsibility.
- The corrected Figure 6 NMPP1 HTML rendered successfully outside the repository;
  detailed visual feedback remains the scientific owner's responsibility.

## Independent review findings and resolutions

Independent scientific-integrity and artifact/security review identified no
change to scientific calculations. Their presentation and operational findings
were resolved as follows:

- Quantitative-table prose now distinguishes retained canonical tables from
  validated report derivations for inclusive 2N–4N positivity and configured-
  reference-relative ratios.
- Baseline provenance states that an independent EdU-negative slope is fitted
  for each acquisition and combined with that acquisition's own G1-derived
  anchor using the configured anchor method.
- Hard-coded “Untreated” presentation labels were replaced with “configured
  matched reference”; no reference identity was selected or changed.
- Full input paths are hidden by default. Explicit audit-path disclosure is
  accompanied by a sharing warning.
- Superseded: an interim Export status listed filenames only. The simplified
  current report removes report-triggered export entirely.
- The dynamic every-sample figure chunk explicitly uses `results: asis`, with a
  structural assertion.
- Verification now installs the modified checkout into an isolated external R
  library and uses that library for tests and rendering.
- Superseded: interim CSV export documentation required an external output
  directory. The current HTML-only report has no CSV export parameter or chunk.

Unresolved: broader package verification has not been run. The rendered layout
awaits scientific-owner feedback. No scientific decision remains unresolved.

### Side-by-side comparison implementation

The comparison-matrix change modified `_report-setup.R`, the standard EdU QMD,
this implementation record, and the focused report-contract test. It read no
experimental inputs and generated no files. Scientific methods changed:
**none**. Tests and render: **NOT RUN — local execution required**. Assumptions:
the validated manifest remains the authoritative source for condition,
reference, biological-replicate, technical-acquisition, and acquisition-prefix
identity; condition display order follows its explicit `condition_index`.

Independent review of the first comparison implementation identified duplicate
replicate-label merging, possible reference/treatment overlap, unequal-stack
sizing, fragile plot-size mutation, and unescaped dynamic headings. These were
resolved by keying rows exclusively with `replicate_index`, excluding every
reference row from treatment cells, weighting each row and the overall figure
by its largest technical-acquisition stack, and escaping condition and prefix
headings at the markup boundary. The report renders each matrix through a
Knitr child chunk with explicit `fig-width` and `fig-height`;
Knitr therefore owns its temporary figure lifecycle. Dynamic headings are
one-line raw HTML with escaped text, so newlines and Markdown control characters
cannot alter document structure. Displayed scientific labels remain unchanged
after browser entity decoding. Focused tests cover duplicate display labels,
unequal stacks, reference exclusion, safe headings, and the child-chunk render
contract. Verification remains **NOT RUN — local execution required**.

The owner's first focused test of the comparison update reported **FAIL 2,
PASS 40, WARN 0**. Both failures were test-only defects: one assertion passed
the integer return value from `anyDuplicated()` to a logical expectation, and
one missing-cell fixture had not constructed its local report object. The
assertion now checks explicitly for integer zero and the fixture constructs the
report before extending it. Production report behavior and scientific
calculations were not implicated. The corrections are **NOT RUN — local
execution required**.

The owner's focused run of the new guide behavior test reported **FAIL 4, WARN
2, PASS 235**. Production legend behavior was not implicated. The test had read
labels from an untrained source manual scale, then compared them with a compact
scale trained during `ggplot_build()`, producing empty/shared-level warnings and
mismatches. The corrected test builds both source and compact plots against
their same existing data before comparing trained labels, limits, and mapped
colors. It continues to verify the blank compact title and unchanged source
scale name. The production report was not changed, and the corrected test has
not been rerun by the implementer.

### Owner-run render finding after implementation

An owner-run render reached report setup but failed because the sourced helper
called the package-internal `validate_analysis_object()` as though it were an
attached export. The helper now resolves that validator explicitly through the
`facspseudocolor` namespace and uses local scalar-string checks instead of
another package-internal helper. A focused test sources `_report-setup.R` into a
clean environment with only base R as its parent, then exercises both the EdU
overview and QC helpers. This correction changes no scientific calculation.
The corrected render was subsequently rerun by the owner and completed cleanly.

### Owner-run focused-test finding after namespace correction

The owner ran the focused `edu-report-contract` test after the namespace fix.
Result: **FAIL 1, PASS 17, no warnings**. The sole failure was in the test's
immutability assertion: it serialized the analysis before adding a synthetic
warning, then removed the warning by assigning `character()`, which did not
necessarily restore the original object representation. The test now sets its
synthetic warning first, serializes that exact pre-helper state, calls the QC
helper, and compares the unchanged post-helper object directly. Production
report code and scientific calculations were not implicated. The corrected
test was rerun by the owner and passed: **PASS 20, FAIL 0, WARN 0, SKIP 0**.

### Owner-run provisional render and network-warning resolution

An owner-run render subsequently completed successfully, but the explicit
Quarto `theme: cosmo` setting attempted to retrieve a font from
`fonts.googleapis.com`. The report must not initiate network access, so the
theme override was removed; Quarto will use its locally available default HTML
styling together with the report's embedded CSS. Structural checks now reject
both the removed theme declaration and a Google Fonts URL. The successful
owner-run report was written outside the repository and is provisional because
it contains the pre-fix network warning. No repository output was generated,
and the corrected no-network render was subsequently rerun by the owner. It
completed successfully with no warnings and produced one self-contained HTML
report in an external temporary output directory.

### Superseded intermediate: simulated focus-group design iteration

Everything in this subsection describes an intermediate design. Its expanded
QC, pairwise comparisons, quantitative and completeness tables, and export
status were subsequently removed by the owner's simplification request. The
current HTML-only contract is documented in “Simplified standard report
structure.”

A structured simulated focus group represented two graduate students, two
postdoctoral fellows, two undergraduate laboratory researchers, two
technicians, and two principal investigators. The coordinator synthesis—not
new experimental evidence—prioritized separation of scientific QC from
software notices, an explicit human visual-review state, a whole-experiment
comparison view, concise readable result tables, less visible implementation
machinery, compact completeness/provenance, and honest export status.

The report now separates data/input checks, normalization/quantitation checks,
required human visual review, and software notices. The retained deprecated-
alias compatibility notice no longer activates the scientific QC alert. A
compact all-condition matrix precedes the pairwise comparisons. Columns follow
explicit manifest condition order, rows use explicit biological-replicate
identity, every technical acquisition is retained, and missing cells remain
visible. Pairwise configured-reference views remain available.

The plot-reading guide states the actual visual contract: DNA limits are shared,
whereas EdU limits can be acquisition-specific unless fixed in configuration;
density color is local, not an absolute cross-panel count scale. Normal reports
hide source code and code controls. Scientist-facing quantitative tables use
selected readable labels while complete tables remain available in an opt-in
audit disclosure. The design-completeness table counts acquisitions without
defining a new expected technical-replicate count. Share-ready provenance uses
only existing metadata and omits directory paths. Export copy explicitly says
when CSV files are neither attached nor generated.

Limitations: this was role-based simulation, not human-subjects research, and
no usability-performance claim is made. No experiment purpose, question,
owner, date, QC threshold, biological interpretation, or new scientific
contrast was inferred. No automated suspicious-distribution rule was added.
Visual and runtime verification of this iteration is **NOT RUN — local
execution required**.

Review follow-up narrowed software-notice routing to an exact match for the
package's known EdU alias-deprecation text. Any unknown retained warning stays
in the data/scientific review lane, including warnings that happen to mention
“software” or “version.” The plot guide now states that vertical position
includes the per-acquisition display-only offset and that this offset is
excluded from positivity and quantitative medians. Readable summary tables now
fail closed when any requested column or display label is absent. The all-
condition overview now rejects missing or empty identities, duplicate prefixes,
panel-order drift, missing ordering identities, and conflicting labels for one
`replicate_index`. Focused regression tests cover each correction. Verification
remains **NOT RUN — local execution required**.

### Direct visual-feedback iteration: all-condition overview

The all-condition overview now suppresses every per-panel density legend and
shows exactly one shared legend per condition block, labeled **Within-panel
relative density**. This wording records that density color is normalized
locally and does not support comparison of absolute density magnitude between
acquisitions. Overview copies of the plots use the compact y-axis title **EdU
signal (display scale)**; the adjacent guide retains the background-subtraction
and per-acquisition display-offset explanation, including that the offset is
excluded from positivity calls and quantitative medians. The every-sample audit
plots are not modified and retain their detailed labels.

Biological-replicate identity is now printed as a compact full-width header
above each replicate row rather than in a width-consuming left column. A new
report parameter, `max_condition_columns` (default `3`), accepts only one
positive integer. Experiments with more conditions are split into consecutive,
explicitly labeled blocks in manifest order. Every replicate row, technical
acquisition, and visible `MISSING` cell is retained in every applicable block.
Focused tests cover the single shared legend contract, local-density label,
compact y title, full-width row-header metadata, preservation of audit labels,
parameter validation, block order, complete condition coverage, and maximum
block width. This iteration changes presentation only; verification is **NOT
RUN — local execution required**.

### Simplified standard report structure

At the owner's explicit request, the standard report now moves directly from
the pre-results QC status into the All-condition plots and Quantitation cards:
the explanatory paragraphs above both galleries were removed while their
headings, group identity, plot captions, accessible labels, and legends remain.
The Experiment overview retains its summary and Samples and replicate structure
table, but no longer displays Design completeness or Input populations and
event counts subsections.

The full Quality control, Side-by-side comparisons, Replicate variation, and
Quantitative tables sections were removed from the visible standard report.
Export was nested under Quantitative tables and was removed with it rather than
left as an orphan subsection; its now-unused parameters and setup export call
were also removed. Artifact review then identified the branch-added comparison,
export, and readable-table report helpers and their focused tests as dormant;
they were removed after repository-wide reference checks confirmed the current
QMD no longer used them. No pre-existing package API was removed. The pre-results QC alert and its retained data/input and
normalization findings remain before all plots, including for non-clear states.
Every-sample plots, detailed quantitative figures, Methods, share-ready
provenance, and optional provenance audit tables remain. Unused comparison and
legacy QC setup objects were removed, while the analysis and report objects and
the `qc_sections` object required by the retained pre-results alert remain.

This is a presentation-only simplification. It does not filter observations or
change sample mappings, exclusions, gates, thresholds, channels,
normalization, denominators, statistics, quantitative objects, or biological
interpretations. Focused structural tests assert the retained heading order,
the absent sections and chunks, the lack of an orphan Export subsection, and
that each gallery heading is followed directly by its rendering chunk.
Verification is **NOT RUN — local execution required**.

### PR 28 macOS CI corrections

The macOS CI log exposed two contract defects. First, intensity provenance had
been generalized so far that it no longer named the actual configured matched
reference. The report model now retains the generic calculation description
and appends the configured reference condition name or names taken from the
already validated overall intensity rows. It verifies that overall and
regional rows expose the same nonmissing reference set before recording that
provenance. No reference is hard-coded, selected, or changed, and intensity
values and normalization calculations are untouched. The focused synthetic
contract now expects its arbitrary configured reference, `Reference`.

Second, the new report-contract test assumed a source-checkout path for files
under `inst/quarto`, which is unavailable in the installed layout used by
`R CMD check`. A test resource resolver now uses `system.file(..., package =
"facspseudocolor")` first and falls back to the relative source-tree resource
only when necessary. All helper sourcing and QMD/helper reads use this resolver;
the prior R-source literal inspection was replaced by inspection of the loaded
plot helper body. No absolute path or platform-specific separator is used.
Verification after these corrections is **NOT RUN — local execution required**.

### Condition-aware compact quantitation card width

Quantitation overview card widths now derive from the explicit manifest
condition count (`C`) and the displayed category count for each independent
row (`P`), never from acquisition count or file order. The presentation-only
initial width is `clamp(140 + 36*C + 18*C*(P-1), 280, 720)` pixels. For the
current compact-verification design, that result is scaled to two thirds and
clamped again to 220–480 pixels. Fractional results are rounded to the nearest
whole CSS pixel so the CSS and raster calculations use exactly the same
deterministic width. Both cards in a row share that result: the overall row
uses `P = 1`, while the regional row uses `P = 3`. Thus the current four-
condition report uses 220-pixel overall cards (the readability floor) and
285-pixel regional cards (428 × 2/3 rounded from 285.333...). The existing
responsive grid wraps cards
rather than shrinking them below the calculated width and permits a card to
use 100% width when its viewport is narrower. The same calculated width also
drives each overview-only ggplot render at the explicit CSS standard of 96
pixels per inch: the current overall and regional plots therefore render at
`220/96` and `285/96` inches, respectively, rather than upscaling a fixed
2.3-inch raster. Condition and category counts must be finite positive
integers; the derived CSS width and inch width must remain finite and positive
or the report stops with a clear error. The independent rows, plot height,
labels, legends, data, scales, mappings, summaries, and scientific methods are
unchanged.

Focused helper tests cover the four-condition scaled widths, rounding, lower and upper clamps,
invalid counts, and finite-input validation. Structural tests require the
report to derive condition count from the explicit sample manifest, declare
the one- and three-category rows, use the validated helper for row width, and
pass its validated 96-pixel-per-inch conversion to the plot render call.
Verification is **NOT RUN — local execution required**.

### Compact quantitative labels and legend keys

The compact Quantitation overview copies now use `% EdU positive` as the y-axis
label for both positivity cards and `EdU Level` for both intensity cards. The
compact fill and colour legend titles are removed because the surrounding card
context already establishes that the unchanged keys encode condition; condition
labels, identities, ordering, and colors remain visible and unchanged.

Both compact positivity cards—the 2N–4N and Regional S-phase percentage
cards—receive the smaller validated legend-key theme: 0.28 cm keys with 0.08 cm
horizontal spacing. This changes legend symbols only, not plotted bars or
biological-replicate points. Both intensity cards retain their existing key
sizes. Full-size plots, axes, legend titles, values, scales,
mappings, colors, statistics, rows, captions, and alt text remain unchanged.
Focused structural tests verify both y-label assignments, title removal, and
positivity-card targeting; helper tests verify key dimensions and reject invalid
sizes. Verification is **NOT RUN — local execution required**.

### Compact embedded-caption removal

The four quantitative source plots carry an internal ggplot caption explaining
what bars or horizontal marks, error bars, and biological-replicate points
represent. In the compact Quantitation overview only, that embedded caption is
now cleared with `labs(caption = NULL)` to reduce repetition. The explicit HTML
metric caption below each card and its figure alt text remain present. Full-size
detailed plots retain their original explanatory captions, and the report's
methods remain unchanged. A structural test requires compact caption removal;
a focused non-mutation test confirms that all four source captions remain
nonempty after presentation copies are created. Scientific calculations,
statistics, error bars, points, bars, scales, and data are unchanged.
Verification is **NOT RUN — local execution required**.

### Robust compact legend-title suppression

A rendered compact plot still displayed `Condition` because the manual scale's
name could survive `labs(fill = NULL, colour = NULL)` through guide resolution.
Compact copies now set `title = NULL` explicitly on both fill and colour
`guide_legend()` objects and apply `legend.title = element_blank()` after the
compact typography theme. This removes the title at both guide and theme levels
without removing condition keys, labels, ordering, or colors. Full-size source
scales retain their original `Condition` title. In addition to structural
checks, a focused behavioral test builds representative fill and colour plots,
confirms the resolved legend-title theme is blank, confirms scale labels remain
identical, and confirms source scale names are not mutated. No scientific value,
mapping, statistic, or plotted mark changed. Verification is **NOT RUN — local
execution required**.

### Overall-intensity x title and accessible mark semantics

The compact all-computed-positive intensity card now hides its redundant
`Population` x-axis title in addition to its already-hidden single-category
tick text and tick mark. The two regional cards retain their category-axis
titles. This remains an overview-only theme change.

The embedded visible ggplot captions remain removed as requested. To preserve
equivalent accessible meaning without adding duplicate visible prose, each
card's figure alt text now states that biological-replicate points are shown—or
explicitly disabled by configuration—that bars or horizontal marks are
condition means as applicable, and that error bars use the configured SD/SEM
display or are disabled when `quant_error_bar: none`. Full-size plot captions
and source plots remain unchanged. No point, bar, error bar, value, mapping,
scale, condition color, statistic, or calculation changed. Focused structural
tests cover both mark types, configuration-aware point/error text, and the two
first-row x-title suppressions. Verification is **NOT RUN — local execution
required**.

Scientific review found the compact intensity label `EdU Level` ambiguous
because both displayed intensity metrics are normalized to their configured
matched reference. Both compact intensity cards now use the concise,
calculation-faithful y-axis label `Relative EdU level`. The detailed captions
continue to state the exact background-subtracted median, technical-acquisition
aggregation, and reference basis. No value, calculation, scale, or full-size
plot changed. Verification is **NOT RUN — local execution required**.

### Final compact quantitation height adjustment

At the owner's request, the compact Quantitation overview base plot allocation
was restored from 4.30 to 2.15 inches while retaining the existing 0.15-inch
allowance per condition-derived legend row. Four conditions therefore return
to a 2.45-inch figure height, and seven conditions use 2.75 inches. The two
independent responsive rows, metric order, bottom legends, x-label suppression,
captions, alt text, source plots, values, scales, and calculations are
unchanged. Focused height expectations were updated. Verification is **NOT RUN
— local execution required**.

### Superseded intermediate: taller two-row quantitative overview

This intermediate iteration used a 4.30-inch base height plus the
condition-count-aware legend allowance. It is superseded by the final
2.15-inch base documented above; it is retained here only as implementation
history. Two independent responsive grids keep
the two overall metrics together in the first row and the two regional metrics
together in the second row. The calculation-faithful first-row titles remain
**2N–4N EdU-positive percentage** and **All computed-positive EdU intensity**;
these metrics have different established populations and denominators. Both
compact overall cards omit redundant single-category tick text and tick marks;
only the 2N–4N percentage card also omits its x-axis title. Full-size plots are
unchanged.

An owner-authorized focused run initially reported **FAIL 1, PASS 198, WARN
0** because a test required integer identity for a correct numeric width of
456. The assertion was corrected without changing production behavior, and the
rerun completed with **PASS 199, FAIL 0, WARN 0, SKIP 0**. A subsequent Quarto
render completed but emitted malformed-fenced-div warnings because the new row
container omitted the class marker before `quantitation-row`; the markup and
its structural assertion were corrected before publication. No data, plot, or
scientific calculation was implicated.

The owner's focused run after this refinement reported **FAIL 1, PASS 198**.
The sole failure was test-only: the two-card maximum width was correctly
calculated as numeric `456`, while the assertion required integer storage type
`456L`. CSS dimensions do not require integer storage, so production code was
left unchanged and the assertion now uses numeric equality. The corrected test
has not been rerun by the implementer.

### Superseded intermediate: double-height quantitation overview

This historical intermediate doubled the compact quantitation data-panel base
allocation from 2.15 to 4.30 inches before adding the existing 0.15-inch
allowance for each condition-derived legend row. Four conditions therefore used
4.60 inches in that iteration. This sizing is superseded: the current base is
2.15 inches and the current four-condition height is 2.45 inches.

The overview now renders two independent responsive grids, each limited to two
columns and using the configured minimum card width. Row one pairs the existing
**2N–4N EdU-positive percentage** plot with the existing **All
computed-positive EdU intensity** plot. These are described as overall metrics,
not as one shared “total S-phase” quantity, because they retain different
established populations and denominators. Row two pairs regional/subphase
positivity with regional/subphase intensity. Because each row has its own grid
container, responsive wrapping cannot intermix overall and regional cards.

On both first-row presentation copies, the redundant single-category x tick
text and tick mark are hidden with supported ggplot2 theme elements. The x-axis
title is hidden only for the 2N–4N percentage card; the overall-intensity x-axis
title is preserved. Explicit captions and alt text retain the complete metric
basis. Regional category axes remain unchanged. The full-size source plots, values,
scales, condition colors, mappings, summaries, and statistics are untouched.
Focused tests cover four- and seven-condition height calculations, two-column
row layout and width, independent-row structure, metric ordering, and the
overview-only axis treatment. Verification is **NOT RUN — local execution
required**.

After artifact review, legend rows and rendered card height were made
condition-count-aware. The guide uses at most two columns and the height grows
deterministically with the number of legend rows; condition labels are neither
truncated nor substituted. Focused tests completed with **PASS 185, FAIL 0,
WARN 0, SKIP 0**. The package installed successfully into the isolated external
library, and the NMPP1 pilot report rendered successfully to the external
temporary output directory. Both independent reviews reported no remaining
findings. Broader package verification was not run.

### Condition-count-aware compact legends

Compact quantitation legend layout now derives from the explicit number of
unique manifest conditions. A validated helper uses at most two legend columns,
computes the required row count without truncating or renaming conditions, and
adds 0.15 inch of figure height per legend row to a 2.15-inch base. Four
conditions therefore retain the reviewed 2.45-inch height; experiments with
more than four conditions receive additional deterministic legend space rather
than a progressively compressed data panel. Only the compact presentation
copies use this layout. Source plots, scales, mappings, labels, condition
identities, values, and statistics remain unchanged. Focused tests cover four
and seven conditions plus invalid counts. Verification is **NOT RUN — local
execution required**.

Follow-up review made the shared legend availability-aware: its source is now
the first manifest-ordered panel whose panel-QC status is `available`, while all
panels remain represented in the overview. If every panel is suppressed, the
overview renders without a density legend and presents an explicit all-panels-
suppressed warning instead of aborting. If an available source unexpectedly
yields no legend, a separate visible legend-unavailable status is shown.
Replicate headers and plot bodies are now separate outer layout rows; every
header receives the same fixed relative height while body height alone responds
to technical-acquisition stack depth. `max_condition_columns` validation now
checks the upper integer bound before conversion, avoiding integer-overflow
coercion. Focused tests cover a suppressed first panel with a later available
legend source, the all-suppressed nonfatal path, retained acquisition coverage,
fixed header allocation under an unequal stack, and an input above
`.Machine$integer.max`. Verification remains **NOT RUN — local execution
required**.

The owner's focused test run after this review reported **FAIL 2, PASS 87,
WARN 0**. Both failures were synthetic-fixture omissions: tests that extended
or replaced the sample manifest and panel list had not extended `panel_qc` in
the same exact prefix order. The fixtures now create matching explicit
available panel-QC rows (while retaining the intended missing condition cell).
Production panel-QC validation was not relaxed. The fixture correction is
**NOT RUN — local execution required**.

The owner's first focused run after the shared-range change reported four test
failures: one stale structural assertion expected a literal
`quantitation-card-` string after ID creation had moved to the validated helper,
and three responsive-overview tests overflowed the node stack because direct
`ggproto` coordinate cloning produced a recursive object. The implementation
now uses the supported `ggplot2::coord_cartesian()` API on the presentation
copy, explicitly preserving the existing x limits, expansion, and clipping and
setting only the approved shared y viewport. Existing y scales,
transformations, and breaks remain attached unchanged, and `default = TRUE`
marks the coordinate replacement as intentional and quiet. The structural
assertion now checks the validated helper's exact `"quantitation-card"` prefix
rather than the obsolete concatenated literal. The corrected tree has not been
rerun by the implementer.

The owner's rerun then reported **FAIL 1, PASS 88, WARN 0** in the five-
condition wrapping fixture. The fixture had used `transform()` inside a loop;
its data-column evaluation retained the template `prefix` rather than reliably
using the loop value. The fixture now assigns character prefixes and status
fields directly, removes irrelevant row names, and asserts exact prefix order
against `as.character(analysis$sample_manifest$prefix)` before exercising the
overview helper. Production validation remains unchanged. This test-only
correction is **NOT RUN — local execution required**.

### Responsive paired-acquisition overview

The fixed all-condition raster overview has been replaced by responsive HTML
plot-card groups. Group membership uses the manifest's explicit `model_group`
field, whose established meaning is one biological-replicate index × one
technical-acquisition identity with a configured reference designation. This
is a reference-first comparison layout only; it defines no shared baseline
fit, reference division, or new scientific contrast. Biological
replicates are not combined. Each full-width header reports the exact replicate,
technical-acquisition, paired-group, and configured-reference identity.
Within each group, the explicit reference acquisition is first; all other
acquisitions follow `condition_index` and manifest order. The global explicit
manifest condition roster supplies visible `MISSING` cards where a normalization
paired acquisition group lacks an expected condition. Every named acquisition remains a separate
card and is never averaged or omitted.

The CSS grid uses `repeat(auto-fit, minmax(...))`, so browser narrowing wraps
cards rather than shrinking them below `min_panel_width` (new validated report
parameter, default 320 px). `max_condition_columns` remains a validated upper
width constraint. Per-card legends are suppressed; each paired acquisition group
uses at most one **Within-panel relative density** legend sourced from its first
available panel. All-suppressed and legend-extraction failure paths retain every
card and show an explicit warning. Compact overview axes and the adjacent
display-offset explanation remain; every-sample audit plots are unchanged.

Focused structural/model tests cover responsive CSS and parameter defaults,
explicit `model_group` separation, reference-first ordering, manifest-order
coverage without acquisition omission, width validation, legend fallback, and
non-mutation of audit labels. Dynamic group/card headings pass through the
existing HTML-escaping boundary. Scientific methods changed: **none**.
Verification is **NOT RUN — local execution required**.

Scientific-language review established the precise UI term **paired acquisition
group**: the explicit biological replicate × technical
acquisition pairing recorded by `model_group`, with a configured reference
designation. The report now says directly that reference-first card ordering is
only a comparison layout and does not define a shared baseline fit, reference
division, or new contrast. Baseline descriptions were corrected throughout:
the standard workflow fits an independent EdU-negative slope for each
acquisition and combines it with that acquisition's G1-derived anchor.

Responsive-group validation now requires nonmissing `is_reference` values, one
shared nonempty `reference_condition` across every row in a group, exactly one
condition match, and exact agreement between that match and `is_reference`.
The CSS maximum-width interpolation now formats validated numeric pixel values
with `%.0f`, avoiding the integer-only `%d` failure after multiplication.
Focused tests cover inconsistent reference declarations/flags, corrected
baseline wording, paired-group language, and the safe width-format contract.
Verification remains **NOT RUN — local execution required**.

### Compact journal-scale overview cards

The responsive all-condition overview now defaults to 220-pixel cards and four
maximum columns, allowing the current four-condition experiment to occupy one
row on a sufficiently wide viewport while retaining responsive wrapping. Each
card renders an actual 2.3 × 2.15 inch figure at 192 DPI; this is not merely a
large raster squeezed into a smaller container. Allowing for compact axes and
card padding, the intended central plot region is approximately 1.5 inches
wide. Exact browser display varies slightly with fonts and label length. Card
images scale to their container without dropping below the configured minimum
before wrapping. The shared legend, compact overview axes, paired acquisition
groups, and all acquisition/MISSING-card coverage are retained. Pairwise and
every-sample audit figures are unchanged. Structural tests cover the new
defaults, physical child-figure dimensions, high-DPI directive, responsive
image sizing, and DPI validation. Verification is **NOT RUN — local execution
required**.

### Caption-below alignment refinement

The compact all-condition overview now removes the embedded ggplot title and
subtitle from its presentation-only panel copies and places the condition,
configured-reference designation, and sample prefix below each image. This
keeps plot tops and data regions vertically aligned when condition text wraps.
Each card is an accessible named group whose `aria-labelledby` target is its
visible condition caption; the sample prefix remains visible immediately below
it. The full every-sample audit plots retain their original title/subtitle
identity labels. Structural tests check caption order, the accessible group
association, removal of overview-only plot titles, and preservation of audit
plot titles. Scientific calculations and source ggplot objects are unchanged.
Verification is **NOT RUN — local execution required**.

### Results-first opening and quantitation overview

The report now opens with the all-condition visual overview, followed
immediately by a compact Quantitation overview, before experiment metadata and
QC. The new quantitation gallery reuses the four existing primary ggplot
objects exactly: regional S-phase positivity, inclusive 2N–4N positivity,
all-computed-positive reference-relative intensity, and regional S-phase
reference-relative intensity. It performs no calculation or resummarization.
Each plot is rendered as a small responsive card using the same configured
minimum width and maximum-column behavior as the acquisition overview. Metric
titles and explicit basis captions appear below the images for aligned plot
tops, and each card is an accessible group named by its visible metric caption.
The later full-size quantitative sections remain unchanged and are explicitly
described as the interpretation/audit views of the same plots. Structural tests
cover section order, all four existing plot-object references, caption order,
and accessible card identifiers. Assumption: “primary quantitative plot
objects” means the four plots already displayed in the detailed EdU report;
no legacy or optional quantitation was added. Scientific methods changed:
**none**. Verification is **NOT RUN — local execution required**.

### Pre-results QC and quantitative-card precision

A prominent pre-results status now precedes the all-condition and quantitation
overviews. It is assembled only from the existing data/input and
normalization/quantitation QC section rows. Any non-`CLEAR` retained state is
shown before the compact quantitative cards with a direct link to detailed QC;
when those retained checks are clear, the status still states that human visual
review is required. It does not add a threshold, filter, or calculation.

The two compact intensity captions now state the exact plotted basis: the
established unweighted mean of technical-acquisition background-subtracted
medians within each biological replicate, relative to that replicate's
configured matched reference; the regional card additionally identifies the
established S-phase region. Each of the four compact scientific plots now has a
concise explicit `fig-alt`. The shared child-plot helper validates alt text and
serializes it safely as a Quarto option. Tests scope the quantitation chunk,
require exactly the four intended plot references and four alt descriptions,
verify its unique indexed card-ID construction and matching `aria-labelledby`
and heading-ID use, and assert that pre-results QC precedes both overview
sections. Scientific methods changed: **none**. Verification is **NOT RUN —
local execution required**.

### Paired-group shared EdU display ranges

With explicit owner approval, the All-condition visual overview now uses one
shared EdU-axis viewport within each manifest-defined paired acquisition
`model_group`. The range is the union of the retained lower and upper display
limits for every available member panel in that exact group. Consequently, no
observation visible under a member's established display range is cropped by
the shared group view. Different groups remain independent.

Only presentation copies receive the shared Cartesian y limits. Existing y
scales, transformations, breaks, event values, display offsets, gates,
quantitation, and full every-sample audit plots are unchanged. Suppressed and
missing panels neither contribute limits nor receive a viewport that could hide
their status text. An available panel with missing, nonfinite, or non-increasing
retained limits stops overview construction instead of accepting a fallback.
Focused tests cover union construction, application to every available group
member, source-coordinate and scale non-mutation, suppressed-panel exclusion,
status preservation, and the all-suppressed path. The scientist-facing report
explains both the within-group shared range and continued independence across
groups.

The prior tautological ARIA assertion was replaced with a validated card-ID
helper. Focused tests require four unique generated quantitation card IDs and
verify that the same selected ID variable is used by both `aria-labelledby` and
the visible heading target. Scientific analytical methods changed: **none**;
this is an approved display-coordinate refinement. Verification is **NOT RUN —
local execution required**.

Artifact review additionally required pairwise validation before constructing a
group union. Every available panel's retained lower and upper limits are now
required to be numeric and finite, with the lower value strictly below the
upper value. Nonnumeric, nonfinite, equal, or reversed pairs stop with a clear
error; equal or reversed errors identify the affected panel prefix. Only after
all member pairs pass is the group union calculated. Focused regressions cover
reversed, equal, and nonnumeric available-panel bounds. This validation neither
changes a valid display range nor introduces a fallback. Verification is
**NOT RUN — local execution required**.

### Owner-authorized verification after paired-group y-axis refinement

The focused report-contract suite initially reported **FAIL 4, PASS 143,
WARN 0**. One stale structural assertion and a recursive coordinate clone were
identified and corrected; no experimental data or scientific calculations were
implicated. The corrected focused suite completed with **PASS 171, FAIL 0,
WARN 0, SKIP 0**. ggplot2 printed informational coordinate-replacement messages
while exercising the overview copies; these were not warnings and are hidden
from the report artifact by the established chunk message setting.

The modified package then installed successfully into the isolated external R
library. The NMPP1 pilot report rendered successfully with Quarto to the
external temporary output directory. The refreshed preview was copied outside
the repository. Broader package verification was not run.

### Compact quantitation legends below plots

Each Quantitation overview card now builds a presentation-only copy with its
legend below the data panel, preventing a right-side legend from compressing
the compact plotting area. Fill and colour legends use supported ggplot2 legend
guides with two columns and row-wise filling so additional condition labels can
wrap. An intermediate iteration reduced compact legend text; that typography
is superseded by the consistency refinement below. The figure height is
increased from 2.15 to 2.45 inches to accommodate the legend without crowding
the metric caption below the image. Full-size detailed plot objects are not
modified. Values, scales, summaries, error bars, replicate points, and
statistics are unchanged. Structural tests scope these requirements to the
compact quantitation chunk. Verification of this latest refinement is **NOT
RUN — local execution required**.

### Superseded intermediate: matched compact typography

Inspection of the All-condition pseudocolor source confirmed that its compact
plots inherit `theme_classic(base_size = 10)`: axis and legend titles resolve
to 10 pt, while tick labels, legend text, and strip text resolve to 8 pt through
the theme's 0.8 relative sizing. An intermediate overview-only theme applied
those same effective sizes to the compact quantitative copies for axis titles,
tick labels, legend titles, legend text, and any facet/strip text. The prior
7.5/7 pt legend overrides were removed. Bottom placement, wrapping, labels,
scales, mappings, bars, points, values, and full-size source plots remain
unchanged. Focused tests verify both the source base theme and every resolved
compact text size. Verification is **NOT RUN — local execution required**.

### Superseded intermediate: 6.5 pt secondary typography

Following visual review, an intermediate Quantitation overview used explicit
8 pt axis and legend titles and 6.5 pt axis tick, legend, and facet/strip text—
approximately 20% below the prior nominal sizes. A validated helper accepts the
title and secondary sizes separately so the requested 6.5 pt secondary value is
not approximated through a relative multiplier. This theme applies only to the
overview plot copies; HTML card captions and full-size plots are unaffected.
Layout, row order, bottom legends, wrapping, labels, scales, mappings, data,
bars, points, and statistics are unchanged. Focused tests verify every explicit
size and reject invalid inputs. Verification is **NOT RUN — local execution
required**.

### Balanced compact typography and accessibility floor

The owner requested smaller compact typography, while accessibility review
identified 6.5 pt secondary text as too small for reliable reading. The final
balance retains 8 pt axis and legend titles and uses 7.5 pt axis ticks and
facet/strip text. For the compact-width verification, legend labels alone are
reduced from 7.5 to 7 pt; legend keys retain their existing configuration. The
helper now enforces an 8 pt minimum for title text, a 7.5 pt minimum for other
secondary text, and a 7 pt minimum for legend labels; nonfinite and missing
sizes also fail validation. These floors apply only to compact Quantitation overview plot
copies. HTML card captions and full-size plots remain unaffected. Focused tests
cover the final sizes and both minimum boundaries. No data, scale, mapping,
label meaning, statistic, or scientific method changed. Verification is **NOT
RUN — local execution required**.
