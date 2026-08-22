# Local FACS intake assistant (milestone 2)

This is a deliberately limited local assistant for inspecting an experiment
directory and assembling a structured analysis configuration deterministically. It cannot write
files, render a report, execute a shell, access a non-loopback API, or return
FCS event matrices. Its proposal is **not authorization to run an analysis**.

## Architecture and trust boundary

```text
authorized experiment directory
  -> deterministic read-only Python tools
  -> compact metadata JSON (not event measurements)
  -> small enum-only advisory request to Ollama on localhost
  -> deterministic host-built proposal with analysis_authorized=false
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
one failure aborts the proposal. The model receives only bounded counts and
returns a tiny schema-constrained advisory status/flag object with no callable
tools. It does not reproduce mappings, paths, recorded details, or authorization
fields. Its role is minimal and advisory: it cannot affect any scientific field,
and the deterministic configuration does not depend on model judgment. The host
computes every mandatory warning flag, permits the model to add only
`human_review_requested`, and rejects omitted warnings or an incoherent status.
The host deterministically constructs the complete proposal. A WSP is
selected only when exactly one candidate exists. With zero or multiple
candidates the input is null and the ambiguity is recorded. QMD files are not
inventoried, read, or used as experimental evidence.

The LLM contribution is optional to the proposal's scientific content, but the
current standard proposal command still requires a successful local advisory
response. Discovery commands do not contact Ollama. A missing or malformed
advisory fails closed rather than changing or partially emitting a proposal.

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
model receives only counts of workspace candidates, FCS files, user channel
claims, user sample claims, user analysis claims, and recorded layout details.
It does not receive filenames, paths, channel identities, sample claims, layout
text, FCS metadata, or event values. The HTTP client ignores
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

The Modelfile uses a 16,384-token context, a 256-token response cap for the tiny
review object, conservative sampling, a fixed seed, and scientific-integrity
instructions suitable for a 16 GB Apple Silicon Mac. Re-run `ollama create`
after changing the Modelfile.

## Request a proposal

First list the exact observed detector/label pairs without contacting Ollama:

```bash
../.venv/bin/python local_facs_assistant/assistant.py EXPERIMENT_DIR --list-channels
../.venv/bin/python local_facs_assistant/assistant.py EXPERIMENT_DIR --list-gates
```

Then, from the repository root with Ollama already running, copy exact pairs
from that output into explicit channel-role claims:

```bash
../.venv/bin/python local_facs_assistant/assistant.py \
  ../2026_08_12_CDC45_MTBPdTAG_V1_Tpl_pFOX_TotalCyclinB \
  --channel-role '{"detector":"FL2-A","label":"PE-A","category":"DNA","feature":"DNA content"}' \
  --channel-role '{"detector":"FL4-A","label":"APC-A","category":"POI","feature":"phospho-FOXM1 T600"}' \
  --sample-map '{"file":"flow_data/1_NoAbCtrl.fcs","condition":"no antibody","time":"0 h","role":"matched_background_control","biological_replicate":1}' \
  --sample-map '{"file":"flow_data/2_NT.fcs","condition":"untreated","time":"0 h","role":"untreated_control","biological_replicate":1}' \
  --analysis '{"name":"pFOXM1 vs DNA","analysis_type":"poi_vs_dna","target_feature":"phospho-FOXM1 T600","dna_feature":"DNA content","population":"Single Cells"}'
```

The sample-map example is intentionally incomplete and will fail: once any
`--sample-map` is supplied, every inventoried FCS must appear exactly once.
This fail-closed behavior prevents partial mappings from silently entering a
proposal. Omit all sample maps for discovery mode.

`--channel-role` is repeatable and takes one JSON object with exactly four
non-empty strings: `detector`, `label`, `category`, and `feature`. Category must
be `POI`, `EdU`, `DNA`, or `other`; use `other` for measured features such as
pH3. Detector and label must exactly match a pair observed by FlowKit. Claims
remain unconfirmed in the proposal. If a biological channel role is unknown,
omit the argument rather than guessing. Channel category/feature is distinct
from the experimental sample role supplied through `--sample-map`.
At most 32 claims are accepted; each string field is limited to 256 characters.

When exactly one WSP is present, its FCS references must match every inventoried
FCS basename exactly once. Missing, extra, duplicate, or ambiguous references
abort intake before any model request.

`--sample-map` is repeatable and requires exactly `file`, `condition`, `time`,
`role`, and `biological_replicate`. File is an exact experiment-relative FCS
path; it is never inferred from order or name. Replicate is a positive integer.
Allowed roles are `matched_background_control`, `experimental_sample`,
`untreated_control`, `vehicle_control`, `no_antibody_control`,
`positive_control`, and `other_control`. These are sample roles, not channel
purposes. Every output row remains unconfirmed and carries
`provenance: user_supplied`.

`--analysis` is repeatable and requires `name`, `analysis_type`,
`target_feature`, `dna_feature`, and `population` (the last two may be null when
the analysis type permits). Types are `poi_vs_dna`, `edu_vs_dna`, `dna_only`,
`feature_vs_dna`, and `other`. `poi_vs_dna`, `edu_vs_dna`, and `dna_only`
strictly require POI, EdU, and DNA target categories respectively;
`feature_vs_dna` is the explicit DNA comparison for category `other`. Feature
names must exactly reference declared channel roles, and every `*_vs_dna`
analysis requires a declared DNA feature. A non-null population must occur in
every inventoried FCS sample in the one reconciled WSP, including background
controls, and it must resolve to the identical full hierarchical path in every
sample. The same terminal name under different parents fails closed. Duplicate
or ambiguous gate names fail closed. Duplicate or
conflicting declarations with the same type, target, DNA feature, and population
fail; distinct population or DNA-feature declarations may coexist under unique
names. No threshold, gate geometry, normalization,
statistic, or background-control relationship is inferred.

The proposal is printed to standard output. Milestone 2 never saves it. Review
all sample and channel mappings manually; filenames and file order are not
accepted as proof of experimental identity. A proposal deliberately contains
`"analysis_authorized": false`.

The Ollama request disables model thinking output and uses a tiny enum-only JSON
Schema. The host rejects Markdown-fenced JSON, truncated output, prose, extra
fields, unknown flags, or overlong arrays. The response is stored only as
`model_review` with `provenance: model_advisory`; it cannot add biological
claims, paths, mappings, methods, or authorization. The host then builds the
full proposal and applies the complete schema and inspection-ledger checks.
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
- Milestone 2 cannot deterministically prove treatment identities, time points,
  title, or biological replicate count without explicit user claims. Sample
  rows are empty in discovery mode; supplied complete mappings are pinned but
  remain unconfirmed. All global authorization flags remain false.
- QMD files are deliberately outside this assistant's inventory, evidence,
  schema, and prompt. They are neither opened nor parsed.
- Workspace inspection reads FlowJo XML metadata but does not evaluate gates or
  load events. Existing FlowKit exporters remain the canonical mechanism for
  reproducing gated populations.
- No configuration generation, writing, rendering, or analysis is included in
  this milestone.
- Discovery flags reject channel, sample, or analysis claim arguments rather
  than silently ignoring them. `--list-gates` reports support separately for
  each experiment-relative FCS path, including each gate's full hierarchical
  path. Analyses still reference terminal population names, which are accepted
  only when that terminal name is unambiguous within every sample and resolves
  to one identical full path across samples.
