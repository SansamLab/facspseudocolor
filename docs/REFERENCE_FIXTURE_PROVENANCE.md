# Reference fixture provenance

This record covers the 26 experimental event-level CSV fixtures distributed in
`inst/extdata/example` and `inst/extdata/example_poi`. The machine-readable
[fixture manifest](reference-fixture-manifest.csv) gives one row per file,
including its SHA-256 hash and event count. These are genuine derived
experimental data, not synthetic examples.

The manifest uses `unknown` where a value is not known, `not established` where
the available evidence does not establish a relationship, `not applicable`
where a field does not apply, and `restricted—not publicly documented` where a
known identity is intentionally withheld. No blank field should be interpreted
as evidence.

## Shared context and authorization

- Biological material: HCT116, a commonly used human cell line. According to
  the data owner, the fixtures contain no human-subject or donor information and
  cannot identify an individual.
- Data owner/controller: Sansam Lab, Oklahoma Medical Research Foundation
  (OMRF).
- Redistribution: the data owner authorizes redistribution of these exact
  fixtures in the public Git repository and installed R package under the
  [CC BY-NC 4.0 International license](https://creativecommons.org/licenses/by-nc/4.0/).
  The separate [data-license notice](../inst/DATA_LICENSE.md) defines its scope
  and distinguishes the fixture license from the package code's MIT license.
- Related preprint: <https://doi.org/10.64898/2026.05.21.726689>, version posted
  May 22, 2026. Its license was verified as
  [CC BY-NC 4.0 International](https://creativecommons.org/licenses/by-nc/4.0/).
- Privacy review: no formal or informal de-identification review occurred before
  release. No identifiers or metadata were deliberately removed. The data owner
  subsequently reviewed the fixture filenames and column labels and attested
  that they contain no sensitive or identifying information. The fixtures are
  therefore not described as formally de-identified.
- Acquisition: CytoFLEX instrument software; its product name and version are
  unknown. No fluorescence compensation was applied.
- Original event identity: not retained; there is no event-ID column, stable
  source-row identity, or separate crosswalk.
- Git history: all 26 fixtures were introduced in commit `c5e10bf` and were
  unchanged when this record was prepared.

Public documentation intentionally excludes original FCS/workspace hashes,
private absolute paths, and the restricted POI source identity.

## EdU fixtures

Classification: **Genuine derived experimental data; experimental lineage
established by exact hash match.** All 18 package fixtures byte-for-byte match
same-named files in the local publication archive. Fixture hashes and event
counts are public in the manifest; original FCS hashes are intentionally not
published.

- Conditions: untreated, 2-hour auxin, and 4-hour auxin in each of two
  biological replicates. Christopher Sansam, PI of the Sansam Lab, approved the
  sample-to-condition and replicate mapping.
- Populations: `single_cells`, `g1`, and `edu_positive`. Christopher Sansam
  created or approved the gates.
- Columns: `Propidium Iodide PE-A` and `Alexa-647 APC-A`; values are raw
  acquisition measurements selected after gates were reapplied.
- Documented reconstruction environment: archived FlowJo workspaces report
  FlowJo 10.9.0, and local publication records document FlowKit 1.3.2 for the
  later reconstruction that reapplied saved transformations and gates. These
  versions document that reconstruction, not the unknown original acquisition
  software version.
- Exact CSV export date: unknown.

## Protein-of-interest (POI) fixtures

Classification: **Exact provenance unverified.** The data owner identifies the
fixtures as derived from the Total MTBP experiment and confirms that the exact
original FCS files and FlowJo workspace are known and available for private
verification. Their identity is `restricted—not publicly documented`, and, at
the data owner's direction, no controlled source comparison is planned for this
record. No POI fixture hash matched any CSV in the local publication archive;
this hash search does not establish derivation.

- Conditions: replicate-matched untagged HCT116 background control, untreated,
  2-hour auxin, and 4-hour auxin in each of two biological replicates. The files
  are biological replicates only; there are no technical replicates.
- Population and gates: `All Events → Cells → Single Cells`. Christopher Sansam
  created or approved the gates.
- Columns: `DNA content` and `Total MTBP (GFP)`; values are FlowJo-scale values.
- Export: Python with FlowKit. The exact FlowKit version, FlowJo workspace
  version, and CSV export date are unknown.

## Evidence boundaries

Exact hashes, event counts, CSV headers, Git history, EdU publication-archive
matches, and the absence of matching POI CSV hashes in that archive were checked
locally. Biological context, ownership, authorization, privacy assessment,
sample mappings, gate approval, acquisition details, compensation status,
value scale, and retained-event-identity status are data-owner attestations.
Unknown and restricted fields above are intentional; they must not be replaced
with inferred values.
