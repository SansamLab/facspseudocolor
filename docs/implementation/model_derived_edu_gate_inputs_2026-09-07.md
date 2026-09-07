# Experimental model-derived Single Cells and G1 inputs for EdU reports

**Status:** EXPERIMENTAL MODEL-DERIVED NON-PRODUCTION

`python/apply_frozen_gate_models.py` reads each explicitly mapped raw FCS and
its explicitly mapped FlowJo workspace into the three population files consumed
by an EdU report. It replaces the report's FlowJo `Single Cells` and `G1`
memberships; the EdU-positive annotation remains the explicitly named FlowJo
population, intersected by direct event identity with the model-derived parent.

The program does not fit, tune, calibrate, or select a model or threshold. It
fails before scoring when an artifact, byte hash, semantic hash, feature order,
channel role, direct identity, source file, or subset relationship is absent or
ambiguous. Outputs go only to a new absent directory and carry an audit
manifest binding source and output hashes, exact role mappings, event counts,
artifacts, and containment results.

The output directory must be an explicit absolute, previously absent path
outside the source repository. Relative paths and the repository itself or any
descendant are rejected. Manifest population paths are path-free basenames and
the R report validator proves their normalized parent is exactly the manifest
directory before hashing them.

FlowKit raw tables may expose the observed two-level
`(PnN detector, PnS stain)` MultiIndex columns.
The implementation explicitly retains the PnN detector component, refuses
blank or duplicate detector names, and validates all configured roles only
after this normalization; any other number of levels stops. When one explicitly mapped acquisition is loaded
from a workspace containing peers, FlowKit may warn that those other samples
were not loaded. Only a `UserWarning` exactly matching
`WSP references <nonempty-peer>, but sample was not loaded.` and not naming the
active FCS is accepted and copied verbatim into provenance. Any added text,
near-miss wording, different warning category, or other warning stops the
operation.

The FlowJo EdU-positive gate table contributes only its canonical direct source
indices. Those indices must be unique, nonnegative integers and an exact subset
of the all-events indices. The emitted EdU-positive table is then selected from
the already normalized all-events measurements by those indices. Thus its
columns, raw values, and event identities are identical to the model-parent
source; FlowKit's population-specific display labels cannot alter the report
channel representation.

## Authorized immutable artifacts

Both rules are the owner-approved frozen `minus_scatter_intensity` development
rules recorded by `manifests/final_development_rule_frozen_v1.json` and
`task_records/2026-08-28_freeze_final_development_rule.md` in the model project.
They use feature schema `fd-feature-v2` in this exact order:

The newer canonical DNA/EdU successor HGB was preferred because it is more
assay-aligned, but its recorded cluster artifact was absent locally and two
authorized read-only authenticated retrieval attempts failed. Nothing was
substituted for, reconstructed from, or fabricated as that HGB. The owner has
therefore made an explicit bounded delegated selection of this older paired
logistic-rule set for this experimental test only. In particular, the older G1
rule was developed against only one
developmental validation experiment and omits the second-marker axes used by
all audited exact expert G1 gates. It is less assay-aligned, experiences a new
parent-distribution shift when evaluated after model-derived rather than
expert-reference Single Cells, and is not recommended for production.

1. `knn25_log_density_2d__fsc_ssc`
2. `asinh150__dna_area`
3. `asinh150__dna_height`
4. `log_ratio__dna_area_height`
5. `linear_tls_signed_distance__dna_height_area`
6. `knn25_log_density_2d__dna_height_area`

| Target | External artifact path | Bytes | Byte SHA-256 | Semantic SHA-256 | Frozen threshold | L2 |
|---|---|---:|---|---|---:|---:|
| Single Cells | `model_code/feature_discovery_runs/COMPACT_STABILITY_SLURM002/local/single_cells_reference/minus_scatter_intensity/frozen_development_rule.json` | 1,069 | `71b459d8fb00930f31a6d289a21f587226fd2d6be4b31ebc782b7bae7839a376` | `f9a53e9fd78d2f39cf5980b8f470c3c44e7a187fa942cbcc0324c563fa68a31a` | 0.205 | 0.01 |
| G1 | `model_code/feature_discovery_runs/COMPACT_STABILITY_SLURM002/local/g1_reference/minus_scatter_intensity/frozen_development_rule.json` | 1,069 | `d653d388d15af3bd22185cb9c8e1addea132e0b5d1d18e60e37c65a7548163b7` | `32723ccc768ed517affa040f2ce62a125dad465399d869f471bc331549f47f6a` | 0.255 | 0.04 |

Paths in the table are relative to the immutable model project
`20260826_ML_Detection_of_G1`; runtime configuration must supply explicit
absolute paths. No artifact or model-development dataset is copied into this
repository.

Artifact selection and limitations were reconciled against the required model
records `task_records/2026-09-07_execute_canonical_dna_edu_successor_model001.md`,
`task_records/2026-09-07_execute_canonical_dna_edu_future_holdout001.md`,
`inventory/CANONICAL_DNA_EDU_GENERATED_METRICS_HETEROGENEITY_REVIEW_001.md`,
and `PROJECT_STATUS.md`, in addition to the frozen-rule manifest and freeze task
record named above.

## Exact bounded Figure 6 mapping

Execution is limited to the existing eight-row configuration mapping; no row
may be selected by file order. Both runs use raw roles `FSC-A` (FSC area),
`SSC-A` (SSC area), `FL2-A` (DNA area), `FL2-H` (DNA height), and `FL4-A`
(EdU), with the exact FlowJo EdU-positive population `S`.

| Run/workspace | Prefix | Exact FCS basename |
|---|---|---|
| 2025-07-28 / `05-Aug-2025.wsp` | `jul28_nt` | `HCT mAC Rif1 CDK1as 7 NT.fcs` |
| 2025-07-28 / `05-Aug-2025.wsp` | `jul28_auxin` | `HCT mAC Rif1 CDK1as 7 5h Auxin.fcs` |
| 2025-07-28 / `05-Aug-2025.wsp` | `jul28_nmpp1` | `HCT mAC Rif1 CDK1as 7 3h NMPP1.fcs` |
| 2025-07-28 / `05-Aug-2025.wsp` | `jul28_auxin_nmpp1` | `HCT mAC Rif1 CDK1as 7 5h Auxin + 3h NMPP1.fcs` |
| 2025-07-30 / `06-Aug-2025.wsp` | `jul30_nt` | `HCT mAC Rif1 CDK1as 7 NT.fcs` |
| 2025-07-30 / `06-Aug-2025.wsp` | `jul30_auxin` | `HCT mAC Rif1 CDK1as 7 5h Auxin.fcs` |
| 2025-07-30 / `06-Aug-2025.wsp` | `jul30_nmpp1` | `HCT mAC Rif1 CDK1as 7 3h NMPP1.fcs` |
| 2025-07-30 / `06-Aug-2025.wsp` | `jul30_auxin_nmpp1` | `HCT mAC Rif1 CDK1as 7 5h Auxin + 3h NMPP1.fcs` |

The report configuration must preserve the two documented biological runs,
their four condition labels, Untreated references, and all established
normalization, positivity, phase, quantitation, and display settings. Only its
data directory changes to the new model-gated output directory; it must set
`flowjo.rebuild: false`, `dna_channel: "DNA content"`, and
`target_channel: "EdU"`.

## Required input contract

Each acquisition must be explicitly mapped to a unique path-free report prefix,
one immutable FCS, one immutable workspace, one unique EdU-positive population,
and five distinct raw channel columns:
FSC area, SSC area, DNA area, DNA height, and EdU. The three emitted CSVs carry
unique `acquisition_id`, `event_identity`, and nonnegative unique integer
`event_index` columns produced from the same acquisition. Every identity equals the
canonical `<acquisition_id>:event_index:<event_index>` spelling, and each input
table must contain exactly one identical acquisition ID. The positive identities
must be a subset of all events.
No mapping is inferred from file order, filename similarity, or nearby files.

The configuration must carry the exact status string
`EXPERIMENTAL MODEL-DERIVED NON-PRODUCTION`, an absent `output_dir`, the two
artifact records (`path`, `byte_sha256`, `semantic_sha256`), and an
`acquisitions` list. A minimal structure is:

```json
{
  "status_label": "EXPERIMENTAL MODEL-DERIVED NON-PRODUCTION",
  "output_dir": "/new/absent/model-gated-inputs",
  "artifacts": {
    "single_cells": {"path": "/exact/single-rule.json", "byte_sha256": "...", "semantic_sha256": "..."},
    "g1": {"path": "/exact/g1-rule.json", "byte_sha256": "...", "semantic_sha256": "..."}
  },
  "acquisitions": [{
    "prefix": "explicit_prefix",
    "acquisition_id": "explicit-stable-id",
    "fcs_path": "/exact/input.fcs",
    "fcs_sha256": "...",
    "workspace_path": "/exact/input.wsp",
    "workspace_sha256": "...",
    "edu_positive_population": "S",
    "channel_roles": {"fsc_area": "FSC-A", "ssc_area": "SSC-A", "dna_area": "FL2-A", "dna_height": "FL2-H", "edu": "FL4-A"}
  }]
}
```

Run `python3 python/apply_frozen_gate_models.py /exact/config.json`, then point
an otherwise unchanged EdU report configuration at the new output directory
with `DNA content` and `EdU` as its report channel names. The dedicated
`facs_edu_model_gated.qmd` template visibly labels the reader-facing HTML
experimental/model-derived and non-production and prints the verified artifact
and input provenance.

## Verification and render history

Verification was executed under the temporary owner-approved local-execution
exception, with builds, installations, generated CSVs, and reports outside the
repository.

- Focused Python model-gating suite: **PASS**, 10 tests.
- Focused R manifest-validation suite: **PASS**, 22 expectations. Two earlier
  focused attempts stopped on synthetic-fixture/validator defects: an OpenSSL
  digest retained its S3 class during JSON serialization, then the synthetic
  rule omitted its threshold. Hashes are now plain scalar strings, the fixture
  is complete, and incomplete rule provenance fails closed.
- Full R suite: **PASS**, 1,834 expectations.
- Package builds and temporary-library installations: **PASS**.
- `R CMD check`: **not completed**. Repository-index warnings triggered the
  required stop; a later attempted invocation also used an incorrect relative
  package-archive path and warned before checking. Neither is recorded as a
  successful package check.
- Real gating v1: **stopped as required**. Expected-but-unhandled FlowKit
  unloaded-peer warnings and the raw MultiIndex identity-column representation
  prevented output acceptance.
- Real gating v2: **gating succeeded; render stopped as required** because the
  population-specific FlowKit EdU-positive columns differed from the normalized
  all-events channel names.
- Real gating v3: **PASS**. FlowJo `S` contributed only validated direct source
  indices; model-derived CSV generation, manifest verification, analysis, and
  HTML rendering completed.

Final local deliverables:

- Manifest:
  `/private/tmp/facspseudocolor-model-gates-verify.NcCX0I/model-gated-inputs-v3/model-gating-manifest.json`
  — SHA-256
  `2cb04f585a773ff0ddfcb16fec28bc0d77bfbebba8360259865f521bba5809e2`.
- Rendered report:
  `/private/tmp/facspseudocolor-model-gates-verify.NcCX0I/rendered-report-v3/facs_edu_model_gated.html`
  — SHA-256
  `57fcdf9a843db3650f4e399d5eaea7417b6b697ba65ca0806c3a9c35741f93a3`.

The completed manifest records this exact event accounting:

| Prefix | All events | Model Single Cells | Model G1 | FlowJo S source | FlowJo S within model parent | G1 contained | S contained |
|---|---:|---:|---:|---:|---:|---|---|
| `jul28_nt` | 24,591 | 19,392 | 7,278 | 8,951 | 8,809 | yes | yes |
| `jul28_auxin` | 25,075 | 19,376 | 7,096 | 8,817 | 8,685 | yes | yes |
| `jul28_nmpp1` | 25,178 | 19,521 | 4,440 | 9,417 | 9,281 | yes | yes |
| `jul28_auxin_nmpp1` | 24,900 | 19,455 | 3,075 | 8,947 | 8,845 | yes | yes |
| `jul30_nt` | 25,257 | 19,466 | 4,599 | 9,519 | 9,362 | yes | yes |
| `jul30_auxin` | 25,332 | 19,491 | 4,317 | 9,083 | 8,954 | yes | yes |
| `jul30_nmpp1` | 24,641 | 19,299 | 2,511 | 9,735 | 9,476 | yes | yes |
| `jul30_auxin_nmpp1` | 24,623 | 18,995 | 1,692 | 9,312 | 9,008 | yes | yes |
| **Total** | **199,597** | **154,995** | **35,008** | **73,781** | **72,420** | **yes** | **yes** |

The v3 preflight and post-run audit confirmed the two workspace hashes and all
eight configured FCS hashes were unchanged. All original FCS, WSP, historical
CSV, model records, and artifacts remained immutable. The v1/v2 output roots
and failed render remain local forensic evidence. No experimental data,
model artifact, generated CSV, manifest, package archive, installed package, or
rendered report is included in the source change or intended for commit.

## Limitations

These logistic rules are development-only, not production validated. Single
Cells developmental validation contained two independent experiments; G1 only
one. The expert gates are provenance-bearing annotations rather than biological
ground truth, and broad biological generalizability is unresolved. G1 is called
only within the model-derived Single Cells parent. Feature extraction is
label-independent but file-local, so each all-events table must represent one
complete acquisition and must not be prefiltered or pooled.

The frozen G1 rule was developed/evaluated within expert-reference Single
Cells. Applying it after the frozen Single Cells model changes the parent event
distribution; this bounded report is the first downstream test of that composed
behavior and must not be described as prior G1 validation performance.
