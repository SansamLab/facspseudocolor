#!/usr/bin/env python3
"""Export events in selected FlowJo populations to CSV files.

The output includes both original FCS measurements and values transformed with
the scales saved in the FlowJo workspace. Event numbers are zero-based FCS row
indices, allowing a row to be traced back to its source event.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import tempfile
from pathlib import Path

import flowkit as fk
import pandas as pd
from lxml import etree

from export_contract import (
    DIRECT_METHOD_ID, LEGACY_PROFILE, MINIMAL_PROFILE, PRODUCTION_PROFILE,
    artifact_record, event_identity_fields, finalize_manifest, load_metadata, new_manifest,
    resolve_production_fcs_files,
)


TRANSFORM_NS = "http://www.isac-net.org/std/Gating-ML/v2.0/transformations"
DEFAULT_POPULATIONS = ("Single Cells", "EDU Positive", "G1", "G2M")
def validate_sample_population_plan(value: object,
                                    population_keys: list[str]) -> dict[str, list[str]]:
    if (not isinstance(value, dict) or not value or
            any(not isinstance(key, str) or not key or
                not isinstance(items, list) or not items or
                len(items) != len(set(items)) or
                any(item not in population_keys for item in items)
                for key, items in value.items())):
        raise ValueError(
            "sample-population plan must map explicit sample IDs to unique, "
            "nonempty requested population-key lists"
        )
    return value


def safe_filename(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("_")


def canonical_flowjo_sample_id(value: str | None) -> str:
    if value is None or re.fullmatch(r"[1-9][0-9]*", value) is None:
        raise RuntimeError("FlowJo sampleID must be a canonical positive integer")
    return value


def flowjo_group_sample_ref_containers(tree: etree._ElementTree) -> list[etree._Element]:
    """Return the exact GroupNode/Group/SampleRefs structures FlowKit consumes."""
    containers = tree.xpath(
        "//*[local-name()='GroupNode']/*[local-name()='Group']"
        "/*[local-name()='SampleRefs']"
    )
    recognized = [
        ref for container in containers
        for ref in container.xpath("./*[local-name()='SampleRef']")
    ]
    all_refs = tree.xpath("//*[local-name()='SampleRef']")
    if len(recognized) != len(all_refs) or any(ref not in recognized for ref in all_refs):
        raise RuntimeError(
            "FlowJo SampleRef occurs outside GroupNode/Group/SampleRefs"
        )
    return containers


def make_flowkit_workspace_copy(source: Path, destination: Path,
                                retained_sample_ids: set[str] | None = None) -> None:
    """Create a FlowKit-only copy, optionally limited to an explicit sample plan.

    FlowJo's minRange is the displayed lower bound, while FlowKit interprets
    the same XML attribute as Gating-ML's A (the magnitude of a negative lower
    bound). Negating a positive FlowJo minimum gives the intended mapping:
    (x - minRange) / (maxRange - minRange).

    Sample removal applies only to ``destination``. The source workspace is
    parsed read-only and is never rewritten.
    """
    tree = etree.parse(str(source))
    post_write_retained_numeric_ids: set[str] | None = None
    if retained_sample_ids is not None:
        samples = tree.xpath(
            "//*[local-name()='SampleList']/*[local-name()='Sample']"
        )
        observed_names: list[str] = []
        observed_numeric_ids: list[str] = []
        for sample in samples:
            sample_nodes = sample.xpath("./*[local-name()='SampleNode']")
            data_sets = sample.xpath("./*[local-name()='DataSet']")
            if len(sample_nodes) != 1 or len(data_sets) != 1:
                raise RuntimeError(
                    "each FlowJo Sample must contain exactly one direct DataSet "
                    "and SampleNode"
                )
            sample_id = sample_nodes[0].get("name", "")
            if not sample_id:
                raise RuntimeError("FlowJo SampleNode is missing its sample name")
            node_numeric_id = canonical_flowjo_sample_id(
                sample_nodes[0].get("sampleID")
            )
            data_numeric_id = canonical_flowjo_sample_id(
                data_sets[0].get("sampleID")
            )
            if node_numeric_id != data_numeric_id:
                raise RuntimeError(
                    "FlowJo DataSet and SampleNode sampleID values do not match"
                )
            observed_names.append(sample_id)
            observed_numeric_ids.append(node_numeric_id)
        if len(observed_names) != len(set(observed_names)):
            raise RuntimeError("FlowJo workspace contains duplicate SampleNode names")
        if len(observed_numeric_ids) != len(set(observed_numeric_ids)):
            raise RuntimeError("FlowJo workspace contains duplicate numeric sampleID values")
        missing = sorted(retained_sample_ids - set(observed_names))
        if missing:
            raise RuntimeError(
                "sample-population plan names samples absent from the FlowJo workspace: "
                + ", ".join(missing)
            )
        all_numeric_ids = set(observed_numeric_ids)
        retained_numeric_ids = {
            numeric_id for sample_name, numeric_id
            in zip(observed_names, observed_numeric_ids)
            if sample_name in retained_sample_ids
        }
        for sample_refs in flowjo_group_sample_ref_containers(tree):
            refs = sample_refs.xpath("./*[local-name()='SampleRef']")
            ref_ids = [canonical_flowjo_sample_id(ref.get("sampleID")) for ref in refs]
            if len(ref_ids) != len(set(ref_ids)):
                raise RuntimeError(
                    "FlowJo group contains duplicate SampleRef sampleID values"
                )
            unresolved = sorted(set(ref_ids) - all_numeric_ids, key=int)
            if unresolved:
                raise RuntimeError(
                    "FlowJo group SampleRef values do not resolve to samples: "
                    + ", ".join(unresolved)
                )
            for ref, numeric_id in zip(refs, ref_ids):
                if numeric_id not in retained_numeric_ids:
                    sample_refs.remove(ref)
        for sample, sample_name in zip(samples, observed_names):
            if sample_name not in retained_sample_ids:
                sample.getparent().remove(sample)
        retained = [
            node.get("name", "")
            for node in tree.xpath(
                "//*[local-name()='SampleList']/*[local-name()='Sample']"
                "/*[local-name()='SampleNode']"
            )
        ]
        if len(retained) != len(retained_sample_ids) or set(retained) != retained_sample_ids:
            raise RuntimeError(
                "temporary FlowJo workspace sample coverage does not exactly match "
                "the sample-population plan"
            )
        retained_refs = [
            canonical_flowjo_sample_id(ref.get("sampleID"))
            for container in flowjo_group_sample_ref_containers(tree)
            for ref in container.xpath("./*[local-name()='SampleRef']")
        ]
        if not set(retained_refs).issubset(retained_numeric_ids):
            raise RuntimeError(
                "temporary FlowJo workspace contains a group reference to a "
                "removed sample"
            )
        post_write_retained_numeric_ids = retained_numeric_ids
    attribute = f"{{{TRANSFORM_NS}}}minRange"
    for element in tree.xpath("//transforms:linear", namespaces={"transforms": TRANSFORM_NS}):
        value = float(element.get(attribute, "0"))
        if value > 0:
            element.set(attribute, str(-value))
    tree.write(str(destination), encoding="UTF-8", xml_declaration=True)
    if post_write_retained_numeric_ids is not None:
        written_tree = etree.parse(str(destination))
        written_refs = [
            canonical_flowjo_sample_id(ref.get("sampleID"))
            for container in flowjo_group_sample_ref_containers(written_tree)
            for ref in container.xpath("./*[local-name()='SampleRef']")
        ]
        if not set(written_refs).issubset(post_write_retained_numeric_ids):
            raise RuntimeError(
                "written temporary FlowJo workspace contains a group reference "
                "to an excluded sample"
            )


def prefix_measurements(frame: pd.DataFrame, prefix: str) -> pd.DataFrame:
    frame = frame.drop(columns=["sample_id"], errors="ignore").copy()
    frame.columns = [f"{prefix}__{column}" for column in frame.columns]
    return frame


def raw_pnn_aliases(frame: pd.DataFrame) -> pd.DataFrame:
    """Copy raw measurements under their exact, unique FCS PnN aliases.

    FlowKit may append a PnS label to a PnN column name.  Retain every
    observed PnN alias instead of imposing an instrument-specific panel: the
    configured analysis channels are separately validated against the export
    contract before analysis.
    """
    raw = frame.drop(columns=["sample_id"], errors="ignore").copy()
    aliases: list[str] = []
    label_pattern = re.compile(
        r"^([^ ]+)(?: ([A-Za-z0-9][A-Za-z0-9+_.()/:~-]*))?$"
    )
    for column in raw.columns:
        label = str(column)
        match = label_pattern.fullmatch(label)
        if match is None:
            raise RuntimeError(
                f"raw FlowKit label cannot be mapped to a PnN alias: {label!r}"
            )
        aliases.append(match.group(1))
    if len(aliases) != len(set(aliases)):
        raise RuntimeError("raw FlowKit labels produce duplicate PnN aliases")
    raw.columns = aliases
    return raw


def write_csv_atomically(frame: pd.DataFrame, path: Path, *, quote_all: bool = False) -> None:
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".partial", dir=str(path.parent)
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            descriptor = -1
            frame.to_csv(handle, index=False,
                         quoting=csv.QUOTE_ALL if quote_all else csv.QUOTE_MINIMAL)
            handle.flush()
            os.fsync(handle.fileno())
        # Publish without clobbering a path concurrently claimed by another run.
        os.link(temporary, path)
        temporary.unlink()
    except BaseException:
        if descriptor >= 0:
            os.close(descriptor)
        temporary.unlink(missing_ok=True)
        raise


def flowjo_saved_counts(workspace_path: Path, populations: list[str]) -> pd.DataFrame:
    tree = etree.parse(str(workspace_path))
    rows: list[dict[str, object]] = []
    for sample_node in tree.xpath("//SampleList/Sample/SampleNode"):
        sample_id = sample_node.get("name", "")
        for population in sample_node.xpath(".//Population"):
            gate_name = population.get("name", "")
            if gate_name in populations:
                rows.append({
                    "sample_id": sample_id,
                    "gate_name": gate_name,
                    "flowjo_saved_count": int(population.get("count", "0")),
                })
    return pd.DataFrame(rows)


def reconcile_planned_counts(analysis: pd.DataFrame, saved: pd.DataFrame,
                             planned_pairs: set[tuple[str, str]]) -> pd.DataFrame:
    analysis_pairs = list(zip(analysis["sample_id"], analysis["gate_name"]))
    saved_pairs = list(zip(saved["sample_id"], saved["gate_name"]))
    if (len(analysis_pairs) != len(planned_pairs) or
            len(set(analysis_pairs)) != len(planned_pairs) or
            set(analysis_pairs) != planned_pairs):
        raise RuntimeError("FlowKit analysis report must contain exactly one row for every planned sample/population pair")
    if (len(saved_pairs) != len(planned_pairs) or
            len(set(saved_pairs)) != len(planned_pairs) or
            set(saved_pairs) != planned_pairs):
        raise RuntimeError("FlowJo saved counts must contain exactly one row for every planned sample/population pair")
    reconciled = analysis.merge(saved, on=["sample_id", "gate_name"], how="left")
    return add_count_qc_columns(reconciled)


def add_count_qc_columns(reconciled: pd.DataFrame) -> pd.DataFrame:
    """Add signed count QC; percentage uses FlowJo saved count as denominator.

    A zero FlowJo saved count has no defined percent denominator and is written
    as an empty/NA percentage, including when both counts are zero.
    """
    reconciled = reconciled.copy()
    for column in ("exported_count", "flowjo_saved_count"):
        values = pd.to_numeric(reconciled[column], errors="coerce")
        if values.isna().any() or (values < 0).any() or (values % 1 != 0).any():
            raise RuntimeError(f"{column} must contain nonnegative integer counts")
        reconciled[column] = values.astype("int64")
    reconciled["count_difference"] = reconciled["exported_count"] - reconciled["flowjo_saved_count"]
    reconciled["absolute_count_difference"] = reconciled["count_difference"].abs()
    denominator = reconciled["flowjo_saved_count"].astype(float)
    reconciled["signed_percent_difference"] = (
        100.0 * reconciled["count_difference"] / denominator
    ).where(denominator != 0, other=float("nan"))
    return reconciled


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workspace", type=Path, help="FlowJo .wsp file")
    parser.add_argument("--fcs-dir", type=Path, help="directory containing the associated FCS files")
    parser.add_argument("--output-dir", type=Path, default=Path("flowjo_population_exports"))
    parser.add_argument("--populations", nargs="+", default=list(DEFAULT_POPULATIONS))
    parser.add_argument("--population-keys", nargs="+",
                        help="manifest keys in the same order as --populations")
    parser.add_argument("--output-suffixes", nargs="+",
                        help="per-acquisition artifact suffixes matching --populations")
    parser.add_argument("--sample-population-plan", type=Path,
                        help="JSON object mapping explicit workspace sample IDs to population keys")
    parser.add_argument("--raw-compatible-columns", action="store_true",
                        help="also emit original raw detector names for direct R consumption")
    parser.add_argument("--contract-metadata", type=Path,
                        help="explicit local JSON provenance and acquisition mapping")
    parser.add_argument("--analysis-channels", nargs=2, metavar=("DNA", "TARGET"),
                        help="raw DNA and target detectors; must match contract metadata")
    parser.add_argument("--export-operation-id")
    parser.add_argument("--profile", choices=(PRODUCTION_PROFILE, MINIMAL_PROFILE, LEGACY_PROFILE),
                        default=LEGACY_PROFILE)
    parser.add_argument("--direct-index-semantics-verified", action="store_true",
                        help="attest the pinned synthetic source-index test passed")
    args = parser.parse_args()

    if args.population_keys and len(args.population_keys) != len(args.populations):
        parser.error("--population-keys must match --populations one-for-one")
    if args.output_suffixes and len(args.output_suffixes) != len(args.populations):
        parser.error("--output-suffixes must match --populations one-for-one")
    population_keys = args.population_keys or [safe_filename(x) for x in args.populations]
    output_suffixes = args.output_suffixes or [f"__{safe_filename(x)}.csv" for x in args.populations]
    if any(Path(value).name != value or not value.endswith(".csv") for value in output_suffixes):
        parser.error("output suffixes must be path-free .csv filename suffixes")
    output_names = [safe_filename(value) for value in args.populations]
    if any(not value for value in output_names) or len(output_names) != len(set(output_names)):
        parser.error("requested population names collide or become empty after filename sanitization")
    if args.profile in (PRODUCTION_PROFILE, MINIMAL_PROFILE) and args.contract_metadata is None:
        parser.error("contract profiles require --contract-metadata")
    if args.profile == MINIMAL_PROFILE and args.export_operation_id is not None:
        parser.error("minimal contract generates --export-operation-id automatically; do not supply one")
    metadata = load_metadata(args.contract_metadata) if args.contract_metadata else {}
    sample_plan = None
    if args.sample_population_plan:
        import json
        sample_plan = json.loads(args.sample_population_plan.read_text(encoding="utf-8"))
        try:
            sample_plan = validate_sample_population_plan(sample_plan, population_keys)
        except ValueError as error:
            parser.error(str(error))

    fcs_dir = args.fcs_dir or args.workspace.parent
    workspace_names = {
        node.get("name", "")
        for node in etree.parse(str(args.workspace)).xpath(
            "//SampleList/Sample/SampleNode")
    }
    operation_sample_ids = set(sample_plan) if sample_plan is not None else workspace_names
    if sample_plan is not None:
        missing_workspace_samples = sorted(operation_sample_ids - workspace_names)
        if missing_workspace_samples:
            parser.error(
                "sample-population plan names samples absent from the FlowJo workspace: " +
                ", ".join(missing_workspace_samples)
            )
    # A contract operation is scoped to the explicit plan, not every sample
    # still referenced by the source workspace. This preserves exact mapping and
    # hash verification for every consumed acquisition while allowing the
    # temporary workspace copy to omit stale, unplanned references.
    if args.profile in (PRODUCTION_PROFILE, MINIMAL_PROFILE) and sample_plan is not None:
        metadata_by_sample = {item["sample_id"] for item in metadata["acquisitions"]}
        missing_metadata_samples = sorted(operation_sample_ids - metadata_by_sample)
        if missing_metadata_samples:
            parser.error(
                "sample-population plan has no contract acquisition mapping for: " +
                ", ".join(missing_metadata_samples)
            )
        metadata = dict(metadata)
        metadata["acquisitions"] = [
            item for item in metadata["acquisitions"]
            if item["sample_id"] in operation_sample_ids
        ]
    if args.profile in (PRODUCTION_PROFILE, MINIMAL_PROFILE):
        if args.analysis_channels is None:
            parser.error("minimal contract requires --analysis-channels DNA TARGET")
        mapping = metadata.get("analysis_mapping", {})
        if [mapping.get("dna_channel"), mapping.get("target_channel")] != args.analysis_channels:
            parser.error("--analysis-channels must exactly match contract analysis_mapping")
    manifest = new_manifest(
        operation_id=args.export_operation_id,
        profile=args.profile, workspace=args.workspace,
        flowkit_version=fk.__version__, metadata=metadata,
        direct_index_semantics_verified=args.direct_index_semantics_verified,
        requested_populations=population_keys,
    )
    if args.profile == MINIMAL_PROFILE:
        # The generated operation ID is also the immutable directory name. A
        # caller supplies only the externally preflighted export root; this
        # prevents a later render from colliding with an earlier operation.
        operation_root = args.output_dir
        operation_root.mkdir(parents=True, exist_ok=True)
        args.output_dir = operation_root / manifest["export_operation_id"]
    acquisition_by_sample = {
        item["sample_id"]: item for item in manifest["acquisitions"]
    }
    if args.profile in (PRODUCTION_PROFILE, MINIMAL_PROFILE):
        fcs_files = [str(path) for path in resolve_production_fcs_files(
            fcs_dir, manifest["acquisitions"], operation_sample_ids
        )]
    else:
        fcs_files = sorted(
            str(path.resolve()) for path in fcs_dir.rglob("*")
            if path.is_file() and path.suffix.casefold() == ".fcs"
        )
    if not fcs_files:
        parser.error(f"no FCS files found under {fcs_dir}")
    if args.profile == LEGACY_PROFILE and workspace_names:
        referenced = [f for f in fcs_files if Path(f).name in workspace_names]
        skipped = [Path(f).name for f in fcs_files if Path(f).name not in workspace_names]
        if referenced:
            if skipped:
                print(f"  ignoring {len(skipped)} FCS file(s) not in the workspace: "
                      f"{', '.join(skipped)}")
            fcs_files = referenced
    if sample_plan is not None:
        by_name = {Path(path).name: path for path in fcs_files}
        if len(by_name) != len(fcs_files):
            parser.error("duplicate FCS basenames prevent explicit sample mapping")
        absent_fcs = sorted(set(sample_plan) - set(by_name))
        if absent_fcs:
            parser.error("sample-population plan FCS files are absent: " +
                         ", ".join(absent_fcs))
        fcs_files = [by_name[sample_id] for sample_id in sample_plan]

    if args.output_dir.exists():
        if args.profile == MINIMAL_PROFILE:
            parser.error(f"refusing to reuse generated operation directory: {args.output_dir}")
        if any(args.output_dir.iterdir()):
            parser.error(f"refusing to write into nonempty operation directory: {args.output_dir}")
    else:
        # For minimal exports this is an atomic no-clobber claim of the UUID
        # named operation directory. Parents already exist above.
        args.output_dir.mkdir(parents=True, exist_ok=False)
    with tempfile.TemporaryDirectory(prefix="flowjo_export_") as temporary_dir:
        adapted_workspace = Path(temporary_dir) / args.workspace.name
        make_flowkit_workspace_copy(
            args.workspace, adapted_workspace,
            set(sample_plan) if sample_plan is not None else None,
        )
        workspace = fk.Workspace(
            str(adapted_workspace),
            fcs_samples=fcs_files,
            filename_as_id=True,
        )
        sample_ids = workspace.get_sample_ids()
        if not sample_ids:
            parser.error("no FCS filenames matched the samples referenced by the workspace")
        if sample_plan is not None:
            if len(sample_ids) != len(sample_plan) or set(sample_ids) != set(sample_plan):
                parser.error(
                    "loaded FlowJo workspace samples do not exactly match the "
                    "sample-population plan"
                )
            sample_ids = list(sample_plan)
        workspace.analyze_samples(use_mp=False)

        if sample_plan is not None:
            population_by_key = dict(zip(population_keys, args.populations))
            gate_errors = []
            for sample_id in sample_ids:
                for population_key in sample_plan[sample_id]:
                    population = population_by_key[population_key]
                    matches = workspace.find_matching_gate_paths(sample_id, population)
                    if not matches:
                        gate_errors.append(
                            f"missing {population!r} ({population_key}) in {sample_id!r}"
                        )
                    elif len(matches) > 1:
                        gate_errors.append(
                            f"ambiguous {population!r} ({population_key}) in {sample_id!r}"
                        )
            if gate_errors:
                raise RuntimeError(
                    "sample-population export plan cannot be satisfied:\n- " +
                    "\n- ".join(gate_errors)
                )

        report = workspace.get_analysis_report()
        selected_report = report[report["gate_name"].isin(args.populations)].copy()
        planned_pairs = None
        if sample_plan is not None:
            population_by_key = dict(zip(population_keys, args.populations))
            planned_pairs = {
                (sample_id, population_by_key[key])
                for sample_id, keys in sample_plan.items() for key in keys
            }
            selected_report = selected_report[
                selected_report.apply(
                    lambda row: (row["sample_id"], row["gate_name"]) in planned_pairs,
                    axis=1,
                )
            ].copy()
        selected_report = selected_report.rename(columns={"count": "exported_count"})
        saved_counts = flowjo_saved_counts(adapted_workspace, args.populations)
        if planned_pairs is not None:
            saved_counts = saved_counts[
                saved_counts.apply(
                    lambda row: (row["sample_id"], row["gate_name"]) in planned_pairs,
                    axis=1,
                )
            ].copy()
            selected_report = reconcile_planned_counts(
                selected_report, saved_counts, planned_pairs
            )
        else:
            selected_report = selected_report.merge(
                saved_counts, on=["sample_id", "gate_name"], how="left"
            )
            selected_report = add_count_qc_columns(selected_report)
        counts_path = args.output_dir / "population_counts.csv"
        write_csv_atomically(selected_report, counts_path)
        manifest["artifacts"].append(artifact_record(
            counts_path, operation_dir=args.output_dir, role="population_count_report",
            row_count=len(selected_report), columns=list(selected_report.columns),
            linkage={"export_operation_id": manifest["export_operation_id"]},
        ))

        emitted_paths: set[Path] = set()
        for population_key, population, output_suffix in zip(
                population_keys, args.populations, output_suffixes):
            for sample_id in sample_ids:
                if sample_plan is not None and population_key not in sample_plan[sample_id]:
                    continue
                matching_paths = workspace.find_matching_gate_paths(sample_id, population)
                if not matching_paths:
                    if args.profile == PRODUCTION_PROFILE or sample_plan is not None:
                        raise RuntimeError(
                            f"requested population {population!r} is missing in acquisition {sample_id!r}"
                        )
                    continue
                if len(matching_paths) > 1:
                    raise RuntimeError(
                        f"population {population!r} is ambiguous in sample {sample_id!r}; "
                        "the script requires a unique gate name"
                    )
                gate_path = matching_paths[0]
                raw = workspace.get_gate_events(sample_id, population, gate_path, source="raw")
                scaled = workspace.get_gate_events(sample_id, population, gate_path, source="xform")
                production = args.profile == PRODUCTION_PROFILE
                if production:
                    source_indices = raw.index.map(str)
                    if any(re.fullmatch(r"0|[1-9][0-9]*", index) is None
                           for index in source_indices):
                        raise RuntimeError(
                            f"source event indices are not canonical nonnegative integers in "
                            f"{sample_id!r} / {population!r}"
                        )
                    if source_indices.has_duplicates:
                        raise RuntimeError(
                            f"duplicate direct source event indices in {sample_id!r} / {population!r}"
                        )
                else:
                    # POI/EdU provenance binds the exported artifact to its
                    # acquisition, but makes no source-row identity claim.
                    source_indices = [""] * len(raw)
                values = prefix_measurements(raw, "raw").join(prefix_measurements(scaled, "scaled"))
                if args.raw_compatible_columns:
                    raw_compatible = raw_pnn_aliases(raw)
                    collisions = set(raw_compatible.columns) & set(values.columns)
                    if collisions:
                        raise RuntimeError("raw-compatible column collision: " +
                                           ", ".join(sorted(map(str, collisions))))
                    values = raw_compatible.join(values)
                acquisition = acquisition_by_sample.get(sample_id, {})
                acquisition_id = acquisition.get("acquisition_id", "unavailable")
                identity = event_identity_fields(args.profile, acquisition_id,
                                                 list(source_indices))
                values.insert(0, "export_manifest_reference",
                              "export-manifest.json + export-manifest.sha256")
                values.insert(0, "export_manifest_digest",
                              manifest["manifest_binding"]["digest"])
                values.insert(0, "export_operation_id", manifest["export_operation_id"])
                values.insert(0, "export_profile", args.profile)
                values.insert(0, "duplicate_occurrence", identity["duplicate_occurrence"])
                values.insert(0, "identity_method_version", identity["identity_method_version"])
                values.insert(0, "identity_method_id", identity["identity_method_id"])
                values.insert(0, "identity_source", identity["identity_source"])
                values.insert(0, "event_identity", identity["event_identity"])
                values.insert(0, "event_index", identity["event_index"])
                values.insert(0, "acquisition_id", identity["acquisition_id"])
                values.insert(0, "sample_id", sample_id)
                values = values.reset_index(drop=True)

                gate = workspace.get_gate(sample_id, population, gate_path)
                parent_path = "/".join(gate_path)
                channels = [str(column).removeprefix("raw__") for column in raw.columns]
                prefix = acquisition.get(
                    "prefix",
                    safe_filename(Path(sample_id).stem) if sample_plan is not None
                    else safe_filename(sample_id),
                )
                population_record = {
                    "population_key": population_key, "gate_name": population,
                    "gate_path": "/".join((*gate_path, population)),
                    "gate_type": gate.gate_type, "parent_population_path": parent_path,
                    "acquisition_id": acquisition_id, "sample_id": sample_id,
                    "prefix": prefix,
                    "channels": channels,
                    "gate_channels": [str(dimension.id) for dimension in gate.dimensions],
                    "row_count": len(values),
                    "identity_field": "event_identity" if production else None,
                    "identity_method_id": DIRECT_METHOD_ID if production else None,
                    "unique_identity_count": int(values["event_identity"].nunique()) if production else None,
                    "duplicate_base_combination_count": 0, "duplicate_row_count": 0,
                    "intentionally_empty": len(values) == 0,
                    "export_operation_id": manifest["export_operation_id"],
                }
                output_path = args.output_dir / f"{prefix}{output_suffix}"
                if output_path.parent.resolve() != args.output_dir.resolve():
                    raise RuntimeError("population artifact path escapes operation directory")
                if output_path in emitted_paths or output_path.exists():
                    raise RuntimeError(f"per-acquisition artifact collision: {output_path.name}")
                emitted_paths.add(output_path)
                write_csv_atomically(values, output_path, quote_all=True)
                artifact = artifact_record(
                    output_path, operation_dir=args.output_dir, role="population_events",
                    row_count=len(values), columns=list(values.columns),
                    identity_columns=(
                        ["acquisition_id", "event_index", "event_identity"]
                        if args.profile == PRODUCTION_PROFILE else []
                    ),
                    intentionally_empty=len(values) == 0,
                    linkage={
                        "export_operation_id": manifest["export_operation_id"],
                        "acquisition_id": acquisition_id,
                        "sample_id": sample_id, "population_key": population_key,
                        "gate_path": population_record["gate_path"],
                        "channels": channels,
                    },
                )
                population_record["artifact_path"] = artifact["path"]
                population_record["artifact_sha256"] = artifact["sha256"]
                manifest["populations"].append(population_record)
                manifest["artifacts"].append(artifact)

    finalize_manifest(manifest, args.output_dir)

    print(f"FLOWJO_OPERATION_DIR={args.output_dir.resolve()}")
    print(f"Exported {len(args.populations)} populations from {len(sample_ids)} samples to {args.output_dir}")
    for population in args.populations:
        count = int(selected_report.loc[selected_report["gate_name"] == population, "exported_count"].sum())
        print(f"  {population}: {count} events")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
