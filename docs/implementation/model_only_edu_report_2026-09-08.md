# Experimental model-only EdU report — 2026-09-08

## Scope and status

- Branch: `experiment/model-only-edu-report`
- Base revision: `664087c15794fef3a84319aeaf55aa8b39b67ef8`
- Canonical source: this `facs_pseudocolor_workflow` repository.
- Status: experimental, research-only, non-production.
- The eight explicitly mapped Figure 6 FCS acquisitions were read without modification. No WSP file was read.
- Single Cells, G1, and EdU Positive are all model-derived. No FlowJo population gate is used.

## Frozen decisions

- Single Cells: existing frozen Single Cells classifier and threshold from its immutable metadata/configuration.
- G1: `DNA_PROJECTED_G1_CLUSTER001`, threshold `0.42`, evaluated within predicted Single Cells.
- EdU Positive: `EDU_POSITIVE_MODEL002`, threshold `0.385`; features are computed on the complete source frame, then calls are restricted to the same predicted Single Cells parent.
- Approved sample identities, FCS SHA-256 digests, panel, and channel roles are in `inst/config/model_only_edu_figure6_mapping.json` (SHA-256 `a874815c83816d6d4291054fceeb08fa0df4b0c864499176dafe3a6a400597fb`).
- The G1 and EdU HistGradientBoosting models are evaluated with the audited portable standard-library predictor. The cluster conversion record is `2026-09-08_execute_portable_hgb_export_job6189638.md`; job 6189638 exited 0.

## Portable-model provenance

- Export archive: `cfebfb6c30b3ff847c1c77b14ef3875f7c18d95d1a5a77ad5e48522aaa06ee57`
- Portable G1 JSON: `7fe256e892adfef6415e3957231bf8f68b7ad401dc909ffd2a801fbde8f4b75d`
- Portable EdU JSON: `9e6ed466687efe5f7523dea633e9e9244381c0e613c86f0944f861c06047fdf4`
- Predictor: `2edf1bbd529159e8752d01731e87f19efbc63b514409b1d246d7a4ffc266a2c3`
- Compatibility report: `60695bb0a5d3d728437aa89aa570b9ee9d2c7be7cacfbbb564ab43b4309f3b38`
- Synthetic compatibility: exact calls and repeatability; maximum absolute probability error `2^-53` for G1 (769 rows) and `2^-54` for EdU (850 rows).

## Exact execution

Model inference was run twice with separate output directories:

```sh
python3 python/apply_model_only_edu_models.py \
  /private/tmp/facspseudocolor-model-only-portable.VhGth0/model-only-config.json

python3 python/apply_model_only_edu_models.py \
  /private/tmp/facspseudocolor-model-only-portable.VhGth0/model-only-config-repeat.json
```

The configurations declare output directories `/private/tmp/facspseudocolor-model-only-portable.VhGth0/model-inputs` and `/private/tmp/facspseudocolor-model-only-portable.VhGth0/model-inputs-repeat`, respectively; that is their only material difference. All 24 population CSVs were byte-identical between runs. The primary manifest SHA-256 is `459ebdbcb3b4ac517d7300063a4429af73a25a82d5afff467593b90a8d2b02cf`.

The report was rendered from a clean temporary package installation:

```sh
R_LIBS='/private/tmp/facspseudocolor-model-only-portable.VhGth0/Rlib:/Users/sansamc/Library/R/arm64/4.5/library' \
  '/Users/sansamc/OMRF Dropbox/Chris Sansam/Synchd/Organized Dropbox/Projects and Data/Projects/Cowork_Projects/2026_Analyze_FACS_Data/.quarto-cli/bin/quarto' render \
  inst/quarto/facs_edu_model_gated.qmd \
  --output-dir /private/tmp/facspseudocolor-model-only-portable.VhGth0/rendered-report \
  -P config:/private/tmp/facspseudocolor-model-only-portable.VhGth0/report-config.yml \
  -P model_gate_manifest:/private/tmp/facspseudocolor-model-only-portable.VhGth0/model-inputs/model-gating-manifest.json
```

Rendered HTML: `/private/tmp/facspseudocolor-model-only-portable.VhGth0/rendered-report/facs_edu_model_gated.html`, SHA-256 `20e0d2a97eaed41e6b6adc1763320315619e4cf7ac3bf09de618edfc6a8e5355`.

## Counts and QC

| Acquisition | All | Single Cells | G1 | EdU Positive |
|---|---:|---:|---:|---:|
| jul28_nt | 24591 | 19392 | 8995 | 8972 |
| jul28_auxin | 25075 | 19376 | 9011 | 8811 |
| jul28_nmpp1 | 25178 | 19521 | 9089 | 9455 |
| jul28_auxin_nmpp1 | 24900 | 19455 | 9065 | 9070 |
| jul30_nt | 25257 | 19466 | 9047 | 9584 |
| jul30_auxin | 25332 | 19491 | 9045 | 9283 |
| jul30_nmpp1 | 24641 | 19299 | 9037 | 9874 |
| jul30_auxin_nmpp1 | 24623 | 18995 | 8857 | 9311 |

No model QC warning was emitted for any acquisition.

## Verification record

- Portable synthetic self-test: PASS.
- Six focused synthetic Python model-only tests, invoked directly because `pytest` is unavailable locally: PASS.
- Python compilation check: PASS.
- Full R suite: 1,847 PASS, 0 FAIL, 0 WARN, 0 SKIP (615.7 seconds).
- Post-review focused analysis suite: 26 PASS, 0 FAIL, 0 WARN, 0 SKIP.
- Quarto render: PASS.
- Earlier stopped attempts: two renders correctly rejected order-sensitive JSON list comparisons; validators were changed to compare named values order-insensitively. One package-build invocation emitted warnings for misplaced install-only options; its archive was moved out of the repository and was not used as evidence.

## Scientific interpretation boundary

The composed model pipeline changes the downstream parent distribution relative to the expert-reference parents used during model development. These results do not establish biological ground truth, production performance, or broad generalizability. No sample mapping, exclusion, channel role, normalization, statistical method, or biological claim was changed by this implementation.
