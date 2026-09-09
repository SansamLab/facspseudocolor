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
- G1: `EDU_DOCUMENTED_EDU_G1_RUN001_EDUPANEL001_GENERIC_FL2_FL4`
  (`EDUPANEL001_generic_fl2_fl4`), threshold `0.29`, evaluated within
  model-derived Single Cells. Its 12 frozen features include the documented
  FL2/FL4 panel scatter, DNA-geometry, and EdU terms.
- EdU Positive: `EDU_POSITIVE_MODEL002`, threshold `0.385`; features are computed on the complete source frame, then calls are restricted to the same predicted Single Cells parent.
- Approved sample identities, FCS SHA-256 digests, panel, and channel roles are in `inst/config/model_only_edu_figure6_mapping.json` (SHA-256 `7aa668a606a02ea77a0d9a0a69d90ffdac91a37bda5642351283173b2ad6ab7d`). The only role additions are the documented `FSC-H` and `SSC-H` inputs required by this frozen G1 model; the eight acquisitions and FCS digests are unchanged.
- The G1 and EdU HistGradientBoosting models are evaluated with the audited portable standard-library predictor. The replacement G1 conversion record is `2026-09-08_execute_portable_documented_edu_g1_job6189690.md`; job 6189690 completed with its synthetic portability checks.

## Portable-model provenance

- Documented-EdU G1 export archive: `b989d28c7e773748ea395c359a3febfeb24fe04f5a64b5ea1a1f3d71d3888816`
- Portable G1 JSON: `ea1ccd90328674fb8cc1b6802a66e194f49b3f9f9c1e03c6577ad01fe880c96e`
- Portable EdU JSON: `9e6ed466687efe5f7523dea633e9e9244381c0e613c86f0944f861c06047fdf4`
- Predictor: `2edf1bbd529159e8752d01731e87f19efbc63b514409b1d246d7a4ffc266a2c3`
- G1 compatibility report: `336a7ee3ec3516b6bdd2af1efd11e6300669a68bb335d7b2460096045f4db64f`
- Synthetic compatibility: exact calls and repeatability; maximum absolute probability error `2^-53` for the replacement G1 (860 rows). The retained EdU portability evidence is unchanged.

## Execution and verification after G1 substitution

The replacement was run twice with configs that differ only in their declared
output directory:

```sh
python3 python/apply_model_only_edu_models.py \
  /private/tmp/facspseudocolor-edu-g1-rerun.EWmNlg/model-only-config.json
python3 python/apply_model_only_edu_models.py \
  /private/tmp/facspseudocolor-edu-g1-rerun.EWmNlg/model-only-config-repeat.json
```

All 24 population CSVs were byte-identical between runs. The primary manifest
SHA-256 is `a684cfe104a65cabd6fc4a4cf71716c95c0a8afc565ba423f6b22b38fa61c9c2`.
No model QC warning was emitted.

| Acquisition | All | Single Cells | EdU G1 | EdU Positive |
|---|---:|---:|---:|---:|
| jul28_nt | 24591 | 19392 | 5791 | 8972 |
| jul28_auxin | 25075 | 19376 | 5687 | 8811 |
| jul28_nmpp1 | 25178 | 19521 | 3195 | 9455 |
| jul28_auxin_nmpp1 | 24900 | 19455 | 2742 | 9070 |
| jul30_nt | 25257 | 19466 | 3526 | 9584 |
| jul30_auxin | 25332 | 19491 | 3484 | 9283 |
| jul30_nmpp1 | 24641 | 19299 | 2056 | 9874 |
| jul30_auxin_nmpp1 | 24623 | 18995 | 1824 | 9311 |

Focused verification passed: Python syntax compilation; both portable synthetic
self-tests; eight synthetic model-only tests; and 78 R assertions across the
analysis, configuration, and Python-boundary tests. Direct source installation
and Quarto rendering passed. The package-build attempt was stopped because the
verified long provenance filenames exceed portable tar path-length guidance;
its archive was not used for the report. Generated Python caches were moved out
of the repository.

The complete R regression suite subsequently passed 1,851 assertions with
0 failures, 0 warnings, and 0 skips in 631.7 seconds.

Rendered HTML:
`/private/tmp/facspseudocolor-edu-g1-rerun.EWmNlg/rendered-report/facs_edu_model_gated.html`,
SHA-256 `7addba1ed8e7ad2c5f304323eac9b2c5b09d9f45ccf46e1409d1b54f9069acb8`.

## Scientific interpretation boundary

The composed model pipeline changes the downstream parent distribution relative to the expert-reference parents used during model development. These results do not establish biological ground truth, production performance, or broad generalizability. No sample mapping, exclusion, normalization, statistical method, or biological claim was changed. The channel-role contract was extended only with the documented `FSC-H` and `SSC-H` inputs required by the frozen EdU-G1 model.
