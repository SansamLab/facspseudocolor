# SYNTHETIC FlowJo source-index fixture

This is a test-fixture recipe only. It creates no biological model or conclusion.
Every value is deterministic arithmetic and is labeled **SYNTHETIC TEST FIXTURE
— NOT EXPERIMENTAL DATA**. Do not place its files in production data, example,
report, or historical experiment directories, and do not commit the generated
`.fcs` or `.wsp` files.

## 1. Generate the SYNTHETIC acquisitions locally

Work from the canonical `facs_pseudocolor_workflow` directory. Use the same
pinned Python environment intended for the FlowKit verification; it must provide
FlowKit and its existing FlowIO dependency. The generator uses FlowIO's supported
`create_fcs` API to write FCS 3.1 float data and adds no production dependency.

Run exactly:

```bash
python3 tools/create_synthetic_flowjo_index_fixture.py \
  --output-dir /private/tmp/SYNTHETIC_ph3_flowkit_index_fixture
```

Expected output is exactly these two acquisitions (360 events each):

- `SYNTHETIC_INDEX_ACQUISITION_1.fcs`
- `SYNTHETIC_INDEX_ACQUISITION_2.fcs`

The success message begins `Created SYNTHETIC test fixture (not experimental
data):` and lists both paths. Stop if FlowIO is unavailable, the destination
does not contain `SYNTHETIC`, any other file is created, or either expected file
is absent. Existing named files are refused. `--overwrite` replaces only those
two named synthetic FCS files; omit it unless replacement is intentional. The
`/private/tmp` directory is temporary and may not persist.

If overwrite installation and rollback both fail, stop immediately. The error
lists any exact `.backup` path retaining an original file; do not delete it.
Recover it manually before retrying.

The four names FlowJo should display are `SYNTHETIC_FSC_A`,
`SYNTHETIC_FSC_H`, `SYNTHETIC_DNA_A`, and `SYNTHETIC_PH3_A`. Acquisition 2
has a small deterministic offset so acquisition scoping is exercised. Each file
contains 40 separated debris/non-single-cell events and 320 intended single-cell
events: 130 in a G1-like synthetic DNA cluster and 190 in a later-DNA cluster;
70 intended single-cell events have a separated high synthetic pH3 signal.

## 2. Create the SYNTHETIC workspace manually in FlowJo

1. Import both `.fcs` files and confirm all four `SYNTHETIC_*` channels appear.
2. On `SYNTHETIC_FSC_A` versus `SYNTHETIC_FSC_H`, draw one parent polygon named
   exactly `SYNTHETIC Single Cells` around the diagonal group. A reproducible
   test-fixture region is FSC-A 40,000–70,000 and FSC-H approximately 94–106%
   of FSC-A; exclude the low, off-diagonal debris group.
3. Under that parent, draw exactly one `SYNTHETIC G1` gate using
   `SYNTHETIC_DNA_A`, approximately 45,000–65,000.
4. Under the same parent—not under G1—draw exactly one
   `SYNTHETIC pH3 Positive` gate using `SYNTHETIC_PH3_A`, above approximately
   70,000.
5. Apply the identical hierarchy to both acquisitions. Confirm all three
   populations are nonempty in each acquisition.
6. Save as
   `/private/tmp/SYNTHETIC_ph3_flowkit_index_fixture/SYNTHETIC_flowkit_index_semantics.wsp`.

These are test-fixture regions, not biological gates or thresholds. Do not
duplicate any of the three gate names elsewhere: the test requires exactly one
matching path per acquisition.

## 3. Verify the SYNTHETIC source-index behavior locally

From the canonical repository directory, with the same pinned environment:

```bash
export SYNTHETIC_FLOWJO_WORKSPACE="/private/tmp/SYNTHETIC_ph3_flowkit_index_fixture/SYNTHETIC_flowkit_index_semantics.wsp"
export SYNTHETIC_FLOWJO_FCS_DIR="/private/tmp/SYNTHETIC_ph3_flowkit_index_fixture"
export SYNTHETIC_SINGLE_CELLS_GATE="SYNTHETIC Single Cells"
export SYNTHETIC_G1_GATE="SYNTHETIC G1"
export SYNTHETIC_PH3_POSITIVE_GATE="SYNTHETIC pH3 Positive"

python3 -m unittest -v tests/python/test_flowkit_source_index_semantics.py
```

Expected success: one test runs, reports `ok`, and ends with `OK`; no repository
file is generated or modified. The test requires exactly the two generator
filenames, embedded generator/nonexperimental metadata, and direct-child gate
paths under `SYNTHETIC Single Cells`. Stop on missing or ambiguous gates, any
missing or extra acquisition, raw/transformed index mismatch, duplicate parent identities,
source-row measurement mismatch, child containment failure, selection of any
experimental file, or any unexpected repository artifact.

FlowKit adds a `sample_id` bookkeeping column to gated event tables. The proof
requires that this is the only non-detector column, validates its value against
the acquisition, and excludes it only from the source-row measurement
comparison. All four detector columns are still compared at the exact reported
source indices.

Report only the command, `PASS` or `FAIL`, failed test name, short relevant
error excerpt, FlowKit version, and unexpected generated files—not a large
successful log. Until this succeeds, direct-index semantics remain unverified.
