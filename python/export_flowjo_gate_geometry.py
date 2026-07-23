#!/usr/bin/env python3
"""Export exact two-dimensional FlowJo gate geometry to a CSV sidecar.

This script is operationally separate from export_flowjo_populations.py. It
does not export events or modify a workspace. Polygon vertices and complete
rectangle boundaries are exported in both FlowKit's transformed coordinates
and inverse-transformed channel coordinates.
"""

from __future__ import annotations

import argparse
import re
import tempfile
from pathlib import Path

import flowkit as fk
import numpy as np
import pandas as pd
from lxml import etree


TRANSFORM_NS = "http://www.isac-net.org/std/Gating-ML/v2.0/transformations"


def make_flowkit_workspace_copy(source: Path, destination: Path) -> None:
    """Adapt FlowJo linear display minima for FlowKit parsing."""
    tree = etree.parse(str(source))
    attribute = f"{{{TRANSFORM_NS}}}minRange"
    for element in tree.xpath(
        "//transforms:linear",
        namespaces={"transforms": TRANSFORM_NS},
    ):
        value = float(element.get(attribute, "0"))
        if value > 0:
            element.set(attribute, str(-value))
    tree.write(str(destination), encoding="UTF-8", xml_declaration=True)


def inverse_transform(transform, values: np.ndarray) -> np.ndarray:
    """Return exact inverse-transformed values or fail explicitly."""
    values = np.asarray(values, dtype=float)
    if hasattr(transform, "inverse"):
        return np.asarray(transform.inverse(values), dtype=float)
    # FlowKit 1.3 provides apply(), but not inverse(), for this FlowJo scale.
    if transform.__class__.__name__ == "WSPLogTransform":
        return transform.offset * np.power(10.0, values * transform.decades)
    raise TypeError(
        f"transform {transform.__class__.__name__} has no exact inverse"
    )


def gate_vertices(gate) -> list[list[float]]:
    """Represent a supported 2-D gate as an ordered, closed polygon."""
    if gate.gate_type == "PolygonGate":
        if len(gate.dimensions) != 2 or len(gate.vertices) < 3:
            raise ValueError("polygon gate must have two dimensions and >=3 vertices")
        vertices = [[float(x), float(y)] for x, y in gate.vertices]
    elif gate.gate_type == "RectangleGate":
        if len(gate.dimensions) != 2:
            raise ValueError("rectangle gate must have exactly two dimensions")
        x_dim, y_dim = gate.dimensions
        if any(value is None for value in (x_dim.min, x_dim.max, y_dim.min, y_dim.max)):
            raise ValueError("rectangle gate must have finite lower and upper bounds")
        vertices = [
            [float(x_dim.min), float(y_dim.min)],
            [float(x_dim.max), float(y_dim.min)],
            [float(x_dim.max), float(y_dim.max)],
            [float(x_dim.min), float(y_dim.max)],
        ]
    else:
        raise TypeError(
            f"unsupported gate type {gate.gate_type}; expected PolygonGate "
            "or RectangleGate"
        )
    if vertices[0] != vertices[-1]:
        vertices.append(vertices[0])
    return vertices


def safe_filename(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("_")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workspace", type=Path, help="FlowJo .wsp file")
    parser.add_argument(
        "--fcs-dir", type=Path,
        help="directory containing the associated FCS files",
    )
    parser.add_argument(
        "--output", type=Path, default=Path("flowjo_gate_geometry.csv")
    )
    parser.add_argument("--population", required=True, help="unique FlowJo gate name")
    args = parser.parse_args()

    fcs_dir = args.fcs_dir or args.workspace.parent
    fcs_files = sorted(
        str(path.resolve())
        for path in fcs_dir.rglob("*")
        if path.is_file() and path.suffix.casefold() == ".fcs"
    )
    if not fcs_files:
        parser.error(f"no FCS files found under {fcs_dir}")

    workspace_names = {
        node.get("name", "")
        for node in etree.parse(str(args.workspace)).xpath(
            "//SampleList/Sample/SampleNode"
        )
    }
    if workspace_names:
        referenced = [path for path in fcs_files if Path(path).name in workspace_names]
        if referenced:
            fcs_files = referenced

    rows: list[dict[str, object]] = []
    with tempfile.TemporaryDirectory(prefix="flowjo_geometry_") as temporary_dir:
        adapted = Path(temporary_dir) / args.workspace.name
        make_flowkit_workspace_copy(args.workspace, adapted)
        workspace = fk.Workspace(
            str(adapted), fcs_samples=fcs_files, filename_as_id=True
        )
        sample_ids = workspace.get_sample_ids()
        if not sample_ids:
            parser.error("no FCS filenames matched the workspace samples")

        for sample_id in sample_ids:
            paths = workspace.find_matching_gate_paths(sample_id, args.population)
            if not paths:
                raise RuntimeError(
                    f"population {args.population!r} was not found in "
                    f"sample {sample_id!r}"
                )
            if len(paths) > 1:
                raise RuntimeError(
                    f"population {args.population!r} is ambiguous in "
                    f"sample {sample_id!r}"
                )
            gate_path = paths[0]
            gate = workspace.get_gate(sample_id, args.population, gate_path)
            if getattr(gate, "use_complement", False):
                raise ValueError("complement gates cannot be represented as one polygon")
            dimensions = gate.dimensions
            for dimension in dimensions:
                if dimension.compensation_ref not in (None, "uncompensated"):
                    raise ValueError(
                        "compensated gate dimensions are not supported because "
                        "the event exporter writes uncompensated raw channels"
                    )
                if dimension.transformation_ref is None:
                    raise ValueError(
                        f"dimension {dimension.id!r} has no workspace transform"
                    )
            vertices = np.asarray(gate_vertices(gate), dtype=float)
            x_transform = workspace.get_transform(
                sample_id, dimensions[0].transformation_ref
            )
            y_transform = workspace.get_transform(
                sample_id, dimensions[1].transformation_ref
            )
            x_raw = inverse_transform(x_transform, vertices[:, 0])
            y_raw = inverse_transform(y_transform, vertices[:, 1])
            for index in range(vertices.shape[0]):
                rows.append({
                    "sample_id": sample_id,
                    "gate_name": args.population,
                    "gate_path": "/".join(gate_path),
                    "gate_type": gate.gate_type,
                    "vertex_index": index + 1,
                    "x_channel": dimensions[0].id,
                    "y_channel": dimensions[1].id,
                    "x_transformed": vertices[index, 0],
                    "y_transformed": vertices[index, 1],
                    "x_raw": x_raw[index],
                    "y_raw": y_raw[index],
                    "workspace": str(args.workspace.resolve()),
                })

    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.output, index=False)
    print(
        f"Exported {args.population!r} geometry for {len(sample_ids)} samples "
        f"to {args.output}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
