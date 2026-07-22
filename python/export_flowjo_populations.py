#!/usr/bin/env python3
"""Export events in selected FlowJo populations to CSV files.

The output includes both original FCS measurements and values transformed with
the scales saved in the FlowJo workspace. Event numbers are zero-based FCS row
indices, allowing a row to be traced back to its source event.
"""

from __future__ import annotations

import argparse
import re
import tempfile
from pathlib import Path

import flowkit as fk
import pandas as pd
from lxml import etree


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
    args = parser.parse_args()

    fcs_dir = args.fcs_dir or args.workspace.parent
    fcs_files = sorted(
        str(path.resolve())
        for path in fcs_dir.rglob("*")
        if path.is_file() and path.suffix.casefold() == ".fcs"
    )
    if not fcs_files:
        parser.error(f"no FCS files found under {fcs_dir}")

    # Only load the FCS files the workspace actually references. This avoids
    # loading unrelated files in the folder (e.g. truncated/corrupt exports that
    # are not part of the analysis) which would otherwise abort the whole load.
    workspace_names = {
        node.get("name", "")
        for node in etree.parse(str(args.workspace)).xpath(
            "//SampleList/Sample/SampleNode")
    }
    if workspace_names:
        referenced = [f for f in fcs_files if Path(f).name in workspace_names]
        skipped = [Path(f).name for f in fcs_files if Path(f).name not in workspace_names]
        if referenced:
            if skipped:
                print(f"  ignoring {len(skipped)} FCS file(s) not in the workspace: "
                      f"{', '.join(skipped)}")
            fcs_files = referenced

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
        selected_report.to_csv(args.output_dir / "population_counts.csv", index=False)

        for population in args.populations:
            population_frames: list[pd.DataFrame] = []
            for sample_id in sample_ids:
                matching_paths = workspace.find_matching_gate_paths(sample_id, population)
                if not matching_paths:
                    continue
                if len(matching_paths) > 1:
                    raise RuntimeError(
                        f"population {population!r} is ambiguous in sample {sample_id!r}; "
                        "the script requires a unique gate name"
                    )
                gate_path = matching_paths[0]
                raw = workspace.get_gate_events(sample_id, population, gate_path, source="raw")
                scaled = workspace.get_gate_events(sample_id, population, gate_path, source="xform")
                values = prefix_measurements(raw, "raw").join(prefix_measurements(scaled, "scaled"))
                values.insert(0, "event_index", values.index.astype(int))
                values.insert(0, "sample_id", sample_id)
                population_frames.append(values.reset_index(drop=True))

            output_path = args.output_dir / f"{safe_filename(population)}.csv"
            if population_frames:
                pd.concat(population_frames, ignore_index=True).to_csv(output_path, index=False)
            else:
                pd.DataFrame(columns=["sample_id", "event_index"]).to_csv(output_path, index=False)

    print(f"Exported {len(args.populations)} populations from {len(sample_ids)} samples to {args.output_dir}")
    for population in args.populations:
        count = int(selected_report.loc[selected_report["gate_name"] == population, "exported_count"].sum())
        print(f"  {population}: {count} events")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
