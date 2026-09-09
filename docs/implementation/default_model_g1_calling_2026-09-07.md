# Default model-only population calling for EdU analysis

**Status:** implemented; verification not run in this implementation turn.

## Architecture decision

G1 model execution remains a pre-analysis repository orchestration step. The
installed R package deliberately does not launch Python, and the immutable
caller needs raw FCS inputs plus external model artifacts that must not be
bundled. R analysis continues to consume deterministic CSVs, but now verifies
the model-gating manifest before normalization whenever the default model mode
is selected.

## Behavior

- EdU `g1_source` and `edu_positive_source` default to `model`.
- Model mode requires an explicit `g1_model_manifest` and fails closed if its
  status, source declaration, prefixes, paths, input/output hashes, artifact
  byte/semantic hashes, schema, or threshold provenance differs.
- The dispatcher invokes only `python/apply_model_only_edu_models.py` in model
  mode. The caller accepts exactly the eight acquisitions in the immutable
  owner-approved Figure 6 mapping and reads no workspace.
- Single Cells uses the approved frozen development rule at 0.205. G1 uses
  `DNA_PROJECTED_G1_CLUSTER001` at 0.42, and EdU Positive uses
  `EDU_POSITIVE_MODEL002` at 0.385. Every artifact, metadata/configuration file,
  schema, feature order, runtime, and threshold is digest-bound.
- G1 features are constructed within predicted Single Cells. EdU features are
  constructed on the complete frame and then scored only within that explicitly
  recorded parent under the branch-only authorization.
- `g1_source: flowjo` explicitly selects the prior contract-aware FlowJo
  Single Cells and G1 export path and rejects a simultaneous model manifest.
- Completed analysis provenance and reader-facing reports state the selected
  G1 source.

The bundled historical EdU fixture and root historical example explicitly set
`g1_source: flowjo`; this preserves their established expected results. New or
migrated configurations that omit the setting enter model mode and must
provide the complete external provenance contract.

## Scientific scope

This changes all three report population memberships when model mode is selected.
It does not fit, tune, recalibrate, download, or bundle a model; alter the
frozen parameters or thresholds; change replicate membership, exclusions,
normalization, statistics, or biological claims; or modify an original input.
Both downstream predictions remain strictly contained in model-derived Single
Cells and are labeled research output rather than expert gates or ground truth.

## Verification plan

From the canonical worktree, with build/check outputs in a new external
temporary directory:

```sh
python3 -m pytest tests/python/test_apply_model_only_edu_models.py tests/python/test_apply_frozen_gate_models.py
Rscript -e 'devtools::test(filter = "config|analysis|python-boundary", stop_on_failure = TRUE)'
Rscript -e 'devtools::test(stop_on_failure = TRUE)'
R CMD build . --output=/new/external/temp
R CMD check --no-manual /new/external/temp/facspseudocolor_*.tar.gz
```

Expected success is zero failures, warnings, or skips in the focused/full
suites and `Status: OK` from package check. Stop on any warning, changed
snapshot, unexpected file, provenance mismatch, or generated repository file.
