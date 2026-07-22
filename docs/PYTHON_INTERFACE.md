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

## Required files

For a configured sample prefix such as `rep1_NT`:

| Mode | Population | Default filename |
|---|---|---|
| EdU | All single cells | `rep1_NT_single_cells.csv` |
| EdU | G1-gated cells | `rep1_NT_g1.csv` |
| EdU | EdU-positive cells | `rep1_NT_edu_positive.csv` |
| POI | All single cells | `rep1_NT_single_cells.csv` |

Suffixes may be changed explicitly in YAML. The package never guesses a file
with a different name.

The EdU-positive file is always required for each replicate reference because
it fits the EdU-negative boundary. It is required for all displayed samples when
EdU apex or EdU-positive quantitation features are requested.

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
