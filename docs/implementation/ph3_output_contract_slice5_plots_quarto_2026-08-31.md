# pH3 output contract Slice 5 — condition plots and thin Quarto report

**Date:** 2026-08-31  
**Status:** IMPLEMENTED — local verification required  
**Canonical source:** `facs_pseudocolor_workflow`, branch
`feature/ph3-output-contract-signal-outcomes`, based on Slice 4 commit
`3632d07`.

## Purpose and bounded scope

This slice consumes the validated Slice 4 pH3 condition report model and adds
only a presentation layer:

- four editable condition-level panels in fixed A-D order;
- a thin canonical Quarto HTML profile; and
- concise QC/provenance display derived from the same model.

The panels are: 4N pH3-positive prevalence; below-4N pH3-positive
prevalence; pH3 signal in 4N pH3-positive cells; and pH3 signal in
below-4N pH3-positive cells. They show biological-replicate points and the
condition mean. SD is shown only where the validated model supplies it (three
or more finite biological-replicate values). No inferential statistic,
comparison bracket, or p-value is created.

The plot layer never recalculates any event values. It preserves a missing or
partial value as such and states unavailable/mixed-basis statuses in the panel
caption. C/D values use reference-relative medians only when the validated
model specifies a reference; otherwise they use direct background-corrected
medians. Mixed signal bases remain visible in point shape/caption and never
receive a combined condition mean.

Any raw-fallback C/D value is explicitly labeled
`RAW—NOT BACKGROUND SUBTRACTED`, including when it is reference-relative, and
the panel caption identifies affected biological-replicate sets. This is a
presentation safeguard only; basis selection is validated upstream.

## Explicit exclusions

This slice does not add CSV/RDS/JSON writers, image/PDF writers, display
offsets, pseudocolor plots, gate geometry, statistics, comparisons, changes to
background correction, signal values, reference ratios, mappings, gates,
thresholds, normalization, EdU, or POI.

## Files read

- workspace `AGENTS.md`;
- Slice 1-4 pH3 output-contract sources, focused tests, and implementation
  records;
- existing Quarto setup and report conventions.

Experimental `.fcs`, `.wsp`, `.prism`, CSV/RDS inputs, generated reports, and
figures were not read or modified.

## Files modified and generated

- `R/ph3_report_plots.R`;
- `inst/quarto/facs_ph3_output_contract.qmd`;
- `tests/testthat/test-ph3-condition-report-plots.R`; and
- this implementation record.
- `NAMESPACE` and `man/plot_ph3_condition_report.Rd` expose and document the
  plot-builder API; and
- `docs/REPORTS.md` records the canonical template and its boundary.

The three untracked audit documents remain outside this slice and untouched.

## Scientific-change declaration

Sample mappings: none. Exclusions: none. Gates: none. Thresholds: none.
Normalization: none. Background fitting: none. Reference computation: none.
Statistics: none. Biological claims: none. This is a deterministic rendering
layer over validated Slice 4 output tables.

## Tests prepared, not run

Focused in-memory tests use only clearly labelled `SYNTHETIC` values. They
check fixed A-D panel construction; approved prevalence denominator labels;
reference-relative versus direct signal labels; explicit mixed-basis
unavailability; and rejection of tampered report provenance. Codex did not run
R, Python, Quarto, tests, renders, builds, or checks.

From the canonical repository, with R >= 4.2 and declared dependencies, run:

```sh
Rscript -e 'devtools::test(filter = "ph3-condition-report-plots", stop_on_failure = TRUE)'
Rscript -e 'devtools::test(filter = "ph3-condition-report-model|ph3-signal-outcomes|ph3-background-regression|ph3-output-model", stop_on_failure = TRUE)'
Rscript -e 'devtools::test(filter = "^ph3$", stop_on_failure = TRUE)'
Rscript -e 'devtools::test(stop_on_failure = TRUE)'
git diff --check
git status --short --untracked-files=all
```

To render the report only after tests pass, install the current branch and
render the installed template entirely in an external temporary directory.
Replace the analysis-RDS placeholder with one explicit completed pH3 analysis:

```sh
ph3_repo="$PWD"
ph3_render_dir="$(mktemp -d)"
mkdir -p "$ph3_render_dir/Rlib"
(
  cd "$ph3_render_dir" &&
  R CMD build "$ph3_repo" &&
  R CMD INSTALL -l "$ph3_render_dir/Rlib" facspseudocolor_*.tar.gz &&
  R_LIBS="$ph3_render_dir/Rlib" quarto render \
    "$ph3_render_dir/Rlib/facspseudocolor/quarto/facs_ph3_output_contract.qmd" \
    --output-dir "$ph3_render_dir/report" \
    -P analysis_rds="/absolute/path/to/completed_ph3_analysis.rds"
)
```

Expected: zero test failures/errors; no unexpected warnings; a silent
`git diff --check`; one HTML file below `$ph3_render_dir/report` only when the
explicit render command is run; and no generated repository artifacts. Stop on any
changed scientific output, recalculation, status/basis mismatch, missing A-D
panel, unexpected file, or render/test failure. Report command, PASS/FAIL,
failed test names, a short relevant excerpt, and unexpected files only.

## Independent review and resolution

## Package-check hygiene correction (2026-09-01)

User-run external `R CMD check --no-manual` reported one WARNING for
non-ASCII source text in `R/ph3_report_plots.R` and one NOTE for undeclared
plot aesthetics plus unqualified `setNames`. The correction replaces the
three em dashes in R string literals with ASCII source escapes (`\\u2014`), which
preserve the rendered labels; qualifies the two remaining package uses as
`stats::setNames`; and declares the ggplot2 data-mask variables with
`utils::globalVariables`. It changes no plotted data, values, labels at
runtime, methods, APIs, or scientific behavior. Codex did not rerun any
project code, tests, build, check, render, or analysis after this correction;
the external package check must be rerun locally.

The independent scientific-integrity review initially found one blocking
presentation issue: a raw-fallback C/D value could be reference-relative yet
be labeled only as relative signal. The plot layer now gives every such panel
the exact visible label `RAW—NOT BACKGROUND SUBTRACTED` and names affected
replicate-set IDs in the caption. The reviewer also found an important wording
issue in `docs/REPORTS.md`: a report supplied `config` can run an analysis.
The documentation now distinguishes that behavior from the calculation-free
plot layer and the reusable `analysis_rds` path. Final scientific re-review was
clean: no mapping, gate, threshold, normalization, reference, statistic, or
biological-claim change was found.

The independent artifact/security review initially found one blocking issue:
the new Quarto profile suppressed warnings. The warning suppression has been
removed, and retained analysis warnings are also shown in a concise report
section. Final artifact/security re-review was clean: no network, telemetry,
secrets, shell execution, unsafe path fallback, destructive I/O, generated
repository artifact, or unrelated scope expansion was found. Reviewers did not
edit code, run project code, render Quarto, or access experimental inputs.

## Assumptions and unresolved uncertainty

Slice 4’s validated condition report model is authoritative. This slice does
not claim that the existing geometry source is verified, does not create a
geometry display, and does not implement the owner-confirmed display-offset or
positive-domain safeguard; those require a later, separately bounded slice.
