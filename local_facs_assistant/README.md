# Local FACS intake assistant (milestone 1)

This is a deliberately limited local assistant for inspecting an experiment
directory and proposing a structured analysis configuration. It cannot write
files, render a report, execute a shell, access a non-loopback API, or return
FCS event matrices. Its proposal is **not authorization to run an analysis**.

## Architecture and trust boundary

```text
authorized experiment directory
  -> deterministic read-only Python tools
  -> compact metadata JSON (not event measurements)
  -> Ollama on localhost
  -> schema-shaped proposal with analysis_authorized=false
  -> human review
```

The tool host resolves every requested path and rejects paths outside the one
experiment directory named on the command line. The host-side inspection
operations are:

- inventory supported files;
- inspect FCS metadata and channels through FlowKit;
- inspect FlowJo sample references and gate names;
- decode text annotations embedded in FlowJo layouts.

These inspections are mandatory host-side preflight operations, not actions the
model must decide to request. Every inventoried FCS and WSP is inspected;
one failure aborts the proposal. The model then receives the bounded inspection
ledger in a single schema-constrained request with no callable tools. A WSP is
selected only when exactly one candidate exists. With zero or multiple
candidates the input is null and the ambiguity is recorded. QMD files are not
inventoried, read, or used as experimental evidence.

There is no generic file reader or writer. Original FCS and WSP inputs are
never modified. Workspace XML has a size limit and parsing fails closed. The
Ollama URL must use HTTP on `localhost`, `127.0.0.1`, or `::1`.

## Requirements

- Python 3.10 or newer
- Ollama running locally
- FlowKit and `jsonschema` in the Python environment used to launch the assistant

The existing project `.venv` provides FlowKit for the canonical FlowJo export
workflow and `jsonschema` for Draft 2020-12 enforcement. No additional Python
packages are required.

All inspection results and prompts remain local. The selected local Ollama
model can receive FCS/WSP inventory filenames and sizes; FlowJo workspace sample
names, sample references, gate names, version, and decoded layout notes; and FCS
event counts, detector/stain channel names, acquisition
date/time, cytometer/instrument identifiers, project, experiment, and operator
metadata. Event values are never returned to the model. The HTTP client ignores
environment proxy settings and permits only an explicitly validated loopback
endpoint.

FCS metadata inspection rejects files larger than 1 GiB. This conservative cap
limits local memory and parser exposure on a 16 GB Mac and does not alter or
downsample an accepted FCS file.

## Set up the model

From this directory:

```bash
ollama pull qwen3.5:9b
ollama create facs-assistant -f Modelfile
```

The Modelfile uses a 16,384-token context, conservative sampling, a fixed seed,
and scientific-integrity instructions suitable for a 16 GB Apple Silicon Mac.

## Request a proposal

First list the exact observed detector/label pairs without contacting Ollama:

```bash
../.venv/bin/python local_facs_assistant/assistant.py EXPERIMENT_DIR --list-channels
```

Then, from the repository root with Ollama already running, copy exact pairs
from that output into explicit channel-role claims:

```bash
../.venv/bin/python local_facs_assistant/assistant.py \
  ../2026_08_12_CDC45_MTBPdTAG_V1_Tpl_pFOX_TotalCyclinB \
  --channel-role '{"detector":"FL2-A","label":"PE-A","category":"DNA","feature":"DNA content"}' \
  --channel-role '{"detector":"FL4-A","label":"APC-A","category":"POI","feature":"phospho-FOXM1 T600"}'
```

`--channel-role` is repeatable and takes one JSON object with exactly four
non-empty strings: `detector`, `label`, `category`, and `feature`. Category must
be `POI`, `EdU`, `DNA`, or `other`; use `other` for measured features such as
pH3. Detector and label must exactly match a pair observed by FlowKit. Claims
remain unconfirmed in the proposal. If a biological role is unknown, omit the
argument rather than guessing. Channel category/feature is distinct from a
sample's experimental role, which remains null in milestone 1.
At most 32 claims are accepted; each string field is limited to 256 characters.

When exactly one WSP is present, its FCS references must match every inventoried
FCS basename exactly once. Missing, extra, duplicate, or ambiguous references
abort intake before any model request.

The proposal is printed to standard output. Milestone 1 never saves it. Review
all sample and channel mappings manually; filenames and file order are not
accepted as proof of experimental identity. A proposal deliberately contains
`"analysis_authorized": false`.

The Ollama request disables model thinking output and uses Ollama's native JSON
Schema format. The host still rejects Markdown-fenced JSON, surrounding prose,
extra fields, or any proposal that fails schema and inspection-ledger checks.
For each run, a deep copy of the bundled schema is further constrained to the
exact authorized directory name, inventory-ordered FCS list, one unconfirmed
sample row per FCS, and deterministic singleton WSP selection. Top-level
evidence is limited to inspected FCS/WSP paths. Channel rows come only from
repeatable user-supplied claims and are reconciled to exact detector/stain pairs
observed by FlowKit, with evidence set to the FCS paths supporting that pair.
FlowJo layout notes are copied verbatim into schema-pinned `recorded_details`
entries with their workspace and layout names. The reusable schema file itself
is not modified.

Alternative local model or endpoint:

```bash
../.venv/bin/python local_facs_assistant/assistant.py EXPERIMENT_DIR \
  --model qwen3:8b --ollama-url http://127.0.0.1:11434
```

## Tests

Tests use only files named as synthetic fixtures, a mocked FlowKit module, and
mocked Ollama messages. They contain no production observations:

```bash
../.venv/bin/python -m unittest discover -s local_facs_assistant/tests -v
```

## Current limitations

- The bundled Draft 2020-12 schema is enforced by `jsonschema`. Host-side
  reconciliation also requires proposed paths, FCS coverage, channels,
  and evidence to agree exactly with deterministic inspections and explicit
  user channel-role claims.
- Milestone 1 cannot deterministically prove treatment identities, time points,
  sample roles, title, or biological replicate count, so those proposal fields
  are schema-constrained to null and all confirmation flags remain false.
- QMD files are deliberately outside this assistant's inventory, evidence,
  schema, and prompt. They are neither opened nor parsed.
- Workspace inspection reads FlowJo XML metadata but does not evaluate gates or
  load events. Existing FlowKit exporters remain the canonical mechanism for
  reproducing gated populations.
- No configuration generation, writing, rendering, or analysis is included in
  this milestone.
