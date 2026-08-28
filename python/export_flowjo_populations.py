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
    DIRECT_METHOD_ID, LEGACY_PROFILE, PRODUCTION_PROFILE,
    artifact_record, event_identity_fields, finalize_manifest, load_metadata, new_manifest,
    resolve_production_fcs_files,
)


TRANSFORM_NS = "http://www.isac-net.org/std/Gating-ML/v2.0/transformations"
DEFAULT_POPULATIONS = ("Single Cells", "EDU Positive", "G1", "G2M")


def safe_filename(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("_")


def make_flowkit_workspace_copy(source: Path, destination: Path) -> None:
    """Adapt FlowJo linear display minima to Gating-ML linear semantics.

    FlowJo's minRange is the displayed lower bound, while FlowKit interprets
    the same XML attribute as Gating-ML's A (the magnitude of a negative lower
    bound). Negating a positive FlowJo minimum gives the intended mapping:
    (x - minRange) / (maxRange - minRange).
    """
    tree = etree.parse(str(source))
    attribute = f"{{{TRANSFORM_NS}}}minRange"
    for element in tree.xpath("//transforms:linear", namespaces={"transforms": TRANSFORM_NS}):
        value = float(element.get(attribute, "0"))
        if value > 0:
            element.set(attribute, str(-value))
    tree.write(str(destination), encoding="UTF-8", xml_declaration=True)


def prefix_measurements(frame: pd.DataFrame, prefix: str) -> pd.DataFrame:
    frame = frame.drop(columns=["sample_id"], errors="ignore").copy()
    frame.columns = [f"{prefix}__{column}" for column in frame.columns]
    return frame


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
    parser.add_argument("--contract-metadata", type=Path,
                        help="explicit local JSON provenance and acquisition mapping")
    parser.add_argument("--export-operation-id")
    parser.add_argument("--profile", choices=(PRODUCTION_PROFILE, LEGACY_PROFILE),
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
    if args.profile == PRODUCTION_PROFILE and (
            args.contract_metadata is None or not args.export_operation_id):
        parser.error("production profile requires --contract-metadata and --export-operation-id")
    metadata = load_metadata(args.contract_metadata) if args.contract_metadata else {}

    fcs_dir = args.fcs_dir or args.workspace.parent
    workspace_names = {
        node.get("name", "")
        for node in etree.parse(str(args.workspace)).xpath(
            "//SampleList/Sample/SampleNode")
    }
    manifest = new_manifest(
        operation_id=args.export_operation_id or "LEGACY-UNVERIFIED",
        profile=args.profile, workspace=args.workspace,
        flowkit_version=fk.__version__, metadata=metadata,
        direct_index_semantics_verified=args.direct_index_semantics_verified,
        requested_populations=population_keys,
    )
    acquisition_by_sample = {
        item["sample_id"]: item for item in manifest["acquisitions"]
    }
    if args.profile == PRODUCTION_PROFILE:
        fcs_files = [str(path) for path in resolve_production_fcs_files(
            fcs_dir, manifest["acquisitions"], workspace_names
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

    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        parser.error(f"refusing to write into nonempty operation directory: {args.output_dir}")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="flowjo_export_") as temporary_dir:
        adapted_workspace = Path(temporary_dir) / args.workspace.name
        make_flowkit_workspace_copy(args.workspace, adapted_workspace)
        workspace = fk.Workspace(
            str(adapted_workspace),
            fcs_samples=fcs_files,
            filename_as_id=True,
        )
        sample_ids = workspace.get_sample_ids()
        if not sample_ids:
            parser.error("no FCS filenames matched the samples referenced by the workspace")
        workspace.analyze_samples(use_mp=False)

        report = workspace.get_analysis_report()
        selected_report = report[report["gate_name"].isin(args.populations)].copy()
        selected_report = selected_report.rename(columns={"count": "exported_count"})
        saved_counts = flowjo_saved_counts(args.workspace, args.populations)
        selected_report = selected_report.merge(saved_counts, on=["sample_id", "gate_name"], how="left")
        selected_report["count_difference"] = (
            selected_report["exported_count"] - selected_report["flowjo_saved_count"]
        )
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
                matching_paths = workspace.find_matching_gate_paths(sample_id, population)
                if not matching_paths:
                    if args.profile == PRODUCTION_PROFILE:
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
                values = prefix_measurements(raw, "raw").join(prefix_measurements(scaled, "scaled"))
                acquisition = acquisition_by_sample.get(sample_id, {})
                acquisition_id = acquisition.get("acquisition_id", "unavailable")
                production = args.profile == PRODUCTION_PROFILE
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
                population_record = {
                    "population_key": population_key, "gate_name": population,
                    "gate_path": "/".join((*gate_path, population)),
                    "gate_type": gate.gate_type, "parent_population_path": parent_path,
                    "acquisition_id": acquisition_id, "sample_id": sample_id,
                    "prefix": acquisition.get("prefix", "unavailable"),
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
                prefix = acquisition.get("prefix", safe_filename(sample_id))
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
                        "acquisition_id": acquisition.get("acquisition_id"),
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

    print(f"Exported {len(args.populations)} populations from {len(sample_ids)} samples to {args.output_dir}")
    for population in args.populations:
        count = int(selected_report.loc[selected_report["gate_name"] == population, "exported_count"].sum())
        print(f"  {population}: {count} events")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
