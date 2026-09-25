# Minimal FlowJo export contract

This is the common production provenance record for POI, EdU, and pH3 FlowJo
exports. It is deliberately limited to facts needed to identify the source
data and reproduce the export:

- SHA-256 hashes of the source FCS files and workspace;
- an explicit acquisition/sample/prefix mapping and selected FlowJo gate(s);
- the raw DNA and target detector names;
- the FlowKit version that performs the export; and
- a newly generated export operation ID plus SHA-256 hashes of every emitted
  event CSV and count report.

The current minimal orchestration applies the same explicitly listed gate keys
to every explicitly listed source FCS file in an operation. It rejects a
missing or ambiguous sample/gate rather than substituting another population.

The source FCS references are confined relative paths and are hash-checked
before export. The workspace hash and all artifact hashes are stored in the
immutable `export-manifest.json` and its SHA-256 sidecar. The exporter creates
the operation ID; do not provide or reuse one manually.

The minimal metadata file follows
[`examples/contract_minimal_flowjo.json`](../examples/contract_minimal_flowjo.json).
It does not require gate-owner or approver paperwork, compensation or display
transform descriptions, a positivity threshold, or a direct source-event-index
attestation.

## pH3 addition

The pH3 workflow still uses `production_direct_identity_v1` when it makes
claims about direct membership in FlowJo child gates. That extension retains
the synthetic source-index check, direct event identities, child-to-parent
containment checks, and the approved pH3 positivity method. Those requirements
do not apply to continuous-signal POI or ordinary EdU exports.
