# Python-to-R interface contract

## Boundary

```text
FlowJo workspace + FCS files
        ↓
python/export_flowjo_populations.py
        ↓
event-level CSV files
        ↓
facspseudocolor R package
```

The installed R package neither invokes Python nor parses FlowJo workspaces.
Python, FlowKit, pandas, and lxml are required only for the optional preprocessing
step in `tools/flowjo-orchestration.R`.

## Production export-operation contract

pH3 production preprocessing uses profile
`production_direct_identity_v1`. One explicitly named operation directory
contains population artifacts, `export-manifest.json`, and the detached
`export-manifest.sha256`. The digest is SHA-256 over the manifest's canonical
UTF-8 JSON bytes (keys sorted, compact separators, one final newline); it is not
embedded in those bytes. Final files use atomic no-clobber publication, so a
concurrent exporter cannot replace a path already claimed by another operation;
a completed operation is never overwritten, and every artifact is re-hashed
before finalization.

Production metadata is an explicit local JSON object. It supplies approval,
transform/compensation context, and an acquisition entry for every workspace
sample, including `acquisition_id`, `sample_id`, configured `prefix`, source FCS
reference, and source FCS SHA-256. Gate owner, approver, date, local approval
record, and positivity method ID/version are mandatory. Missing values stop the
export; usernames and file order are never used as substitutes.

Direct identity is
`<acquisition_id>:event_index:<source-index>`. The CSV also carries the source
`event_index` as character data, identity source/method/version,
`duplicate_occurrence`, operation ID, `export_manifest_digest`, and an explicit
full-manifest reference. `export_manifest_digest` is the SHA-256 of canonical
`manifest_binding_object_v1`: schema, operation/profile/method, immutable
workspace/software/approval contexts, explicit acquisitions, and requested
populations. It excludes creation time, population/artifact ledgers, and the
full-manifest digest, so the same truthful binding digest can be written into
every event row without a hash cycle. The separate `export-manifest.sha256` is
the SHA-256 of the complete finalized canonical manifest bytes and covers the
artifact ledger. The two digests have distinct documented purposes.
The exporter rejects duplicate indices within a population and never generates
an index from output row order. R orchestration explicitly reads all identity
columns as character and has no sequential fallback.

Static inspection establishes that the index returned by FlowKit is used
directly; it does **not** prove that the supported FlowKit `get_gate_events()`
index denotes the same source FCS row across parent and child gates. Production
support therefore remains conditional until the focused SYNTHETIC pinned-
environment test is run successfully. Only then may the operator pass
`--direct-index-semantics-verified`. Composite identity is not implemented.

The pinned test independently compares every gated raw row with the raw source-
acquisition row at the reported index; child-set membership alone is not
sufficient proof. Production also loads only the exact relative FCS references
declared in acquisition metadata, verifies each SHA-256, and requires exact
agreement with workspace sample IDs. It never scans a directory to choose a
production input.

The explicit legacy profile is `legacy_count_only_unverified_v1`. It records a
prominent warning and cannot claim identity, containment, or verified geometry.
Legacy event identity fields are blank, and R orchestration refuses a legacy or
ambiguous profile as production input.

Population artifacts are emitted once per acquisition and requested population
inside the dedicated operation directory. Their ledger entries carry exact
acquisition, sample, population key, gate path, channels, row count, hash, and
intentional-empty status. R orchestration reads and returns these immutable
files directly; it does not create, overwrite, rename, or falsely reference an
unledgered per-sample copy. Immediately before returning them, R invokes the
local standard-library verifier. It revalidates the detached full-manifest
digest and every artifact hash, size, row count, and ordered schema, then checks
every nonempty row's operation ID, manifest-binding digest, acquisition ID,
sample ID, and production profile against the ledger. Mutation or mismatch is
fatal.

Geometry is optional. If not requested, the population manifest records
`geometry_overlay_status: not_requested`. If requested, the geometry exporter
must receive the completed population `--operation-dir` and exact
`--population-key`; it verifies the operation ID, workspace digest, sample,
gate path, and channels, then creates an immutable linked
`geometry-manifest.json`/`.sha256` member. Similar filenames never establish
linkage. Geometry remains nonquantitative.
The geometry artifact and supplement also carry the same structured linkage
scope: one entry per acquisition/sample with exact full gate path and ordered
channels. This scope is reconciled with the geometry rows before finalization;
null scalar linkage cannot stand in for a multi-acquisition artifact.
Verified R geometry orchestration handles one export operation at a time and
writes directly inside that operation directory; it does not rewrite or merge
the ledgered artifact.

## Exact FlowJo gate geometry

`python/export_flowjo_gate_geometry.py` is a separate optional extractor for
two-dimensional FlowJo polygon and rectangle coordinates. It does not replace
or modify `python/export_flowjo_populations.py`.

The repository-level orchestration helper runs the extractor and maps FlowJo
sample IDs to configured prefixes:

```r
library(facspseudocolor)
config <- read_facs_config("config.yml")
source("tools/flowjo-orchestration.R")
prepare_flowjo_gate_geometry_external(
  config,
  output_file = "csv/ph3_gate_geometry.csv",
  population_key = "ph3_positive"
)
```

The sidecar records the workspace, sample ID, configured prefix, unique gate
path, gate type, detector channels, ordered vertices, transformed coordinates,
and inverse-transformed raw coordinates. Extraction stops for ambiguous gates,
unsupported gate types, complement gates, channel mismatches, or transforms
that cannot be inverted exactly.

## Required files

For a configured sample prefix such as `rep1_NT`:

| Mode | Population | Default filename |
|---|---|---|
| EdU | All single cells | `rep1_NT_single_cells.csv` |
| EdU | G1-gated cells | `rep1_NT_g1.csv` |
| EdU | EdU-positive cells | `rep1_NT_edu_positive.csv` |
| POI | All single cells | `rep1_NT_single_cells.csv` |
| PH3 | All single cells | `rep1_NT_single_cells.csv` |
| PH3 | G1-gated cells | `rep1_NT_g1.csv` |
| PH3 | User-gated pH3-positive cells | `rep1_NT_ph3_positive.csv` |

Suffixes may be changed explicitly in YAML. The package never guesses a file
with a different name.

The EdU-positive file is required for every EdU acquisition. It supplies the
positive edge used to identify that acquisition's EdU-negative population;
the package then fits and subtracts the background independently for that
acquisition.

## Required columns and values

Every CSV must contain the exact columns configured as `dna_channel` and
`target_channel`.

- Both columns must be numeric.
- Rows represent individual events.
- The R package treats values as supplied and does not compensate, transform,
  repair, simulate, or substitute them.
- Nonfinite values are counted and reported; scientific functions decide
  whether enough finite events remain for their calculation.
- Additional columns are retained and ignored unless explicitly used.
- Repository orchestration retains `event_index` when available so population
  membership can be traced to the original Single Cell event.

For PH3 mode, an empty pH3-positive population is valid and is exported as a
header-only CSV with the required columns. The resulting percentage is zero.
The pH3-positive gate must be drawn by the user in FlowJo and nested within the
Single Cells population; the R package never estimates a positivity threshold.

The Python exporter writes raw detector values selected by
`dna_source_channel` and `target_source_channel`, renaming them to the configured
R-facing column names.

## Sample mapping

Each YAML sample defines:

- `label`: displayed condition name.
- `prefix`: stem used for output CSV filenames.
- `fcs`: source FCS filename, needed only by optional Python orchestration.

Each replicate defines a `reference` matching exactly one sample label. Prefixes
must be unique across the experiment, and replicate condition order must be
stable and compatible.

## Validation responsibility

The optional orchestration layer checks Python dependencies, workspace paths,
population exports, detector columns, and FCS sample matching. The R package
then independently validates expected filenames, required columns, numeric
types, event counts, references, and configuration structure before analysis.
