# pH3 threshold-calibration diagnostic — implementation and owner-decision handoff

**Date:** 2026-08-31  
**Status:** Implemented for local execution; no production threshold selected or changed.

## Purpose and scientific boundary

This diagnostic supports calibration of a future deterministic pH3 positivity
rule for the six-sample Cdk1as experiment. It uses the existing FlowJo
`pH3 Positive` membership only as a verified training/QC label. It does not
replace `ph3_positive_member`, calculate a production cutoff, or alter any
existing pH3 output.

Owner-confirmed sample structure:

| Replicate set | Control | Treatment |
| --- | --- | --- |
| Clone 6 | untreated | 2 h NMPP1 |
| Clone 7 | untreated | 2 h NMPP1 |
| Clone 15 | untreated | 2 h NMPP1 |

Owner-confirmed calibration scope: each treatment uses the threshold calibrated
from its own matched biological-replicate untreated control. Pooled controls
must not derive a shared threshold. The diagnostic does not decide the final
candidate false-positive rate or acceptance criterion.

## Preconditions and fail-closed contract

The diagnostic requires a production `ph3_output_model` derived from new,
verified direct-identity exports and their manifest/containment provenance. It
rejects legacy/count-only/ambiguous exports, altered FlowJo labels, raw-signal
fallback bases, guessed sample mappings, incomplete controls, fewer than 100
negative training events in a sample, and nonfinite eligible event signals.

No sample identity, condition, replicate membership, gate, DNA boundary,
existing FlowJo label, background-regression method, or scientific result is
changed. The existing FlowJo gate is not treated as numerical authority.

## Implemented diagnostic outputs

The new `diagnose_ph3_threshold_calibration()` function creates an in-memory,
diagnostic-only object with:

- exact explicit sample/condition manifest and stop-rule results;
- per-event corrected pH3-negative data with stable identities;
- per-sample and pooled descriptive negative distributions (R type-8
  quantiles), including residual signal-versus-DNA slope/correlation QC;
- matched-control-only candidate threshold/concordance curves across explicit
  false-positive-rate candidates, applied unchanged to the corresponding
  control and treatment sample(s);
- overlap with the retained FlowJo training label: counts, observed negative
  false-positive rate, precision, recall, and Jaccard index;
- cross-replicate stability of matched-control candidate thresholds.

`write_ph3_threshold_calibration()` stages all tables, a JSON manifest, and a
concise README, then atomically publishes them only to a new explicit output
directory outside every explicitly named protected root. It refuses overwrite,
output beneath an input/configuration/experiment/repository root, and leaves no
partial final directory if a write fails.

## Required inputs for the Cdk1as run

1. A separate new pH3 configuration whose `ph3_input_profile` is
   `production_direct_identity_v1`, with verified export operations/manifest.
   The historical configuration and CSV/FCS/WSP inputs are not eligible.
2. An explicit six-row sample manifest matching the active configuration:
   `experiment_id`, `replicate_set_id`, `sample_id`, `condition_id`.
3. A user-approved external diagnostic output location that is not inside the
   experiment's input/export directory.
4. Explicit descriptive candidate FPRs. The command below uses a broad
   inspection grid (`0.001, 0.005, 0.01, 0.025, 0.05`) only as candidates; it
   does not approve any of them as the final method.

## Local-only execution

Run from the canonical repository. Codex agents must not run this command.
It reads real data only after verified exports are prepared, and writes only to
the external directory named by `PH3_CALIBRATION_OUTPUT`.

```sh
export PH3_CONFIG="/absolute/path/to/new_verified_cdk1as_ph3_config.yml"
export PH3_EXPERIMENT_ROOT="/absolute/path/to/facspseudocolor_analysis_pH3InCdk1AsCells"
export PH3_EXPORT_ROOT="/absolute/path/to/the_verified_export_operation_directory"
export PH3_CALIBRATION_OUTPUT="/absolute/path/to/new_empty_external_calibration_output"

Rscript - <<'RSCRIPT'
library(facspseudocolor)

config <- read_facs_config(Sys.getenv("PH3_CONFIG"))
analysis <- analyze_facs_experiment(config)
model <- facspseudocolor:::build_ph3_output_model(analysis)
model <- facspseudocolor:::apply_ph3_background_regression(model)

sample_manifest <- data.frame(
  experiment_id = rep("cdk1as_ph3", 6L),
  replicate_set_id = c("clone6", "clone6", "clone7", "clone7", "clone15", "clone15"),
  sample_id = c(
    "clone6_untreated", "clone6_2h_nmpp1",
    "clone7_untreated", "clone7_2h_nmpp1",
    "clone15_untreated", "clone15_2h_nmpp1"
  ),
  condition_id = rep(c("untreated", "nmpp1_2h"), 3L),
  stringsAsFactors = FALSE
)

diagnostic <- diagnose_ph3_threshold_calibration(
  model = model,
  calibration_manifest = sample_manifest,
  candidate_false_positive_rates = c(0.001, 0.005, 0.01, 0.025, 0.05),
  expected_sample_count = 6L
)
write_ph3_threshold_calibration(
  diagnostic,
  output_dir = Sys.getenv("PH3_CALIBRATION_OUTPUT"),
  input_dirs = c(
    Sys.getenv("PH3_EXPERIMENT_ROOT"),
    Sys.getenv("PH3_EXPORT_ROOT"),
    normalizePath(".", winslash = "/", mustWork = TRUE)
  )
)
RSCRIPT
```

Before running, replace the literal IDs in `sample_manifest` with the exact IDs
from the new verified configuration. The function intentionally stops when
they differ. Do not alter historical configuration/data to make them match.

Expected successful output: one new external directory containing
`README.md`, `calibration_manifest.json`, `source_provenance.csv`, `qc_stop_rules.csv`,
`negative_events.csv`, `negative_distribution_by_sample.csv`,
`negative_distribution_pooled.csv`, `candidate_concordance.csv`, and
`control_threshold_stability.csv`. No repository or input-directory file is
expected to be generated.

Stop and report the reason if validation rejects an input, a stop-rule row is
not `pass`, an unexpected file appears beneath the input or repository tree,
or any export/identity/provenance condition is ambiguous.

## Required owner decision after reviewing the diagnostic

The owner must explicitly approve the future production threshold rule before
implementation, including:

1. the fixed candidate FPR (or another exact data-grounded decision rule);
2. numerical stability/concordance QC acceptance limits and what happens when
   they fail;
3. whether the FlowJo-label concordance is an approval gate, descriptive QC,
   or both;
4. the output/report representation of threshold provenance and QC.

No code in this change selects those decisions.

## Change inventory and verification state

Canonical source: `facs_pseudocolor_workflow`.

Modified/generated source and documentation:

- `R/ph3_threshold_calibration.R`
- `NAMESPACE`
- `man/diagnose_ph3_threshold_calibration.Rd`
- `man/write_ph3_threshold_calibration.Rd`
- `tests/testthat/test-ph3-threshold-calibration.R` (SYNTHETIC only)
- this handoff record

Experimental inputs read: none. Experimental inputs modified: none. Existing
pipeline/configuration/results changed: none. New diagnostic outputs generated:
none. Tests/checks: **NOT RUN — local execution required**. The test file is
prepared but must be run locally after implementation review:

```sh
Rscript -e 'devtools::test(filter = "ph3-threshold-calibration", stop_on_failure = TRUE)'
Rscript -e 'devtools::test(filter = "^ph3$", stop_on_failure = TRUE)'
Rscript -e 'devtools::test(stop_on_failure = TRUE)'
git diff --check
git status --short --untracked-files=all
```

Expected verification: no test failures/errors or unexpected warnings, silent
`git diff --check`, and only the files in the inventory above. Do not request
or provide large successful logs; report command, PASS/FAIL, failing test
names/relevant excerpt, and unexpected files.

## Independent review record

Scientific-integrity review initially found that source provenance was validated
transiently then dropped from the diagnostic, that arbitrary nonempty
positivity-method IDs were accepted, and that the negative-event minimum could
be silently lowered. The implementation now retains per-event and
per-acquisition analysis/config/export-operation/manifest/positivity/containment
bindings, requires `flowjo_owner_approved_positive_population_v1`, and enforces
the established 100-event minimum. Final scientific re-review: **clean**, with
no blocking, important, or advisory findings.

Artifact/security review initially found partial-final-output risk and an
incomplete protected-root boundary. The writer now uses same-parent staging and
atomic publish, cleans failed staging, has a synthetic failure-path test, and
requires explicit experiment/export/repository protected roots. Final
artifact/security re-review: **clean**, with no blocking, important, or
advisory findings.
