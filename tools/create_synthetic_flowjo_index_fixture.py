#!/usr/bin/env python3
"""Create two deterministic SYNTHETIC FCS acquisitions for a local FlowJo test."""

from __future__ import annotations

import argparse
import os
import tempfile
from pathlib import Path


GENERATOR_ID = "SYNTHETIC_flowjo_index_fixture_v1"
CONSTRUCTION_ID = "SYNTHETIC_DETERMINISTIC_ARITHMETIC_V1_NO_RANDOM_SEED"
CHANNELS = (
    "SYNTHETIC_FSC_A",
    "SYNTHETIC_FSC_H",
    "SYNTHETIC_DNA_A",
    "SYNTHETIC_PH3_A",
)
SAMPLE_IDS = ("SYNTHETIC_INDEX_ACQUISITION_1", "SYNTHETIC_INDEX_ACQUISITION_2")


def replace_file(source: Path, destination: Path) -> None:
    """Narrow wrapper used for atomic fixture installation and restoration."""
    os.replace(source, destination)


def require_synthetic_destination(destination: Path) -> Path:
    """Resolve and validate a user-supplied fixture directory without creating it."""
    resolved = destination.expanduser().resolve(strict=False)
    if "SYNTHETIC" not in str(resolved).upper():
        raise ValueError("destination path must contain SYNTHETIC")
    if resolved == Path(resolved.anchor):
        raise ValueError("destination must not be a filesystem root")
    return resolved


def _cluster(count: int, fsc: float, dna: float, ph3: float, acquisition: int):
    """Return visibly separated deterministic points; values have no biological meaning."""
    points = []
    acquisition_offset = 800.0 * (acquisition - 1)
    for index in range(count):
        fsc_a = fsc + acquisition_offset + float((index % 11) - 5) * 180.0
        fsc_h = fsc_a * (0.985 + float((index % 5) - 2) * 0.002)
        dna_a = dna + acquisition_offset + float((index % 13) - 6) * 220.0
        ph3_a = ph3 + acquisition_offset + float((index % 7) - 3) * 240.0
        points.append((fsc_a, fsc_h, dna_a, ph3_a))
    return points


def build_synthetic_events(acquisition: int):
    """Build 360 events: 40 debris and 320 intended single-cell events."""
    if acquisition not in (1, 2):
        raise ValueError("acquisition must be 1 or 2")
    events = []
    # Debris/non-single-cell group: deliberately off the FSC diagonal.
    for index in range(40):
        fsc_a = 18000.0 + 500.0 * (acquisition - 1) + (index % 9) * 160.0
        events.append((fsc_a, 9000.0 + (index % 7) * 120.0,
                       20000.0 + (index % 8) * 190.0,
                       12000.0 + (index % 6) * 170.0))
    # Intended Single Cells: G1-like/later-DNA crossed with low/high pH3 signal.
    events.extend(_cluster(100, 52000.0, 54000.0, 18000.0, acquisition))
    events.extend(_cluster(30, 54500.0, 56000.0, 92000.0, acquisition))
    events.extend(_cluster(150, 57500.0, 108000.0, 20000.0, acquisition))
    events.extend(_cluster(40, 60000.0, 112000.0, 96000.0, acquisition))
    return tuple(events)


def metadata_for(acquisition: int):
    sample_id = SAMPLE_IDS[acquisition - 1]
    return {
        "fil": f"{sample_id}.fcs",
        "sampleid": sample_id,
        "src": "SYNTHETIC TEST FIXTURE - NOT EXPERIMENTAL DATA",
        "synthetic_fixture": "SYNTHETIC TEST FIXTURE - NOT EXPERIMENTAL DATA",
        "synthetic_generator": GENERATOR_ID,
        "synthetic_construction": CONSTRUCTION_ID,
        "synthetic_acquisition_id": sample_id,
        "synthetic_interpretation": "SYNTHETIC - NO BIOLOGICAL INTERPRETATION",
        "synthetic_provenance": (
            "SYNTHETIC - NO COMPENSATION; NO TRANSFORMATION; GENERATED FLOAT VALUES"
        ),
    }


def output_paths(destination: Path):
    return tuple(destination / f"{sample_id}.fcs" for sample_id in SAMPLE_IDS)


def write_fixture(destination: Path, overwrite: bool = False):
    destination = require_synthetic_destination(destination)
    paths = output_paths(destination)
    existing = [path for path in paths if path.exists()]
    if existing and not overwrite:
        raise FileExistsError("refusing to overwrite existing fixture file(s): " +
                              ", ".join(path.name for path in existing))
    if any(path.is_symlink() for path in paths):
        raise ValueError("refusing to write through a symbolic-link fixture path")
    destination.mkdir(parents=True, exist_ok=True)
    try:
        from flowio import create_fcs
    except ImportError as error:
        raise RuntimeError(
            "FlowIO is required. Use the pinned Python environment that supplies FlowKit/FlowIO."
        ) from error

    staged = {}
    backups = {}
    installed = []
    committed = False
    try:
        # Stage and fsync both acquisitions before changing either final path.
        for acquisition, final_path in enumerate(paths, start=1):
            events = build_synthetic_events(acquisition)
            flattened = [value for event in events for value in event]
            descriptor, temporary_name = tempfile.mkstemp(
                prefix=f".{final_path.name}.", suffix=".tmp", dir=destination
            )
            try:
                with os.fdopen(descriptor, "wb") as handle:
                    create_fcs(handle, flattened, list(CHANNELS),
                               metadata_dict=metadata_for(acquisition))
                    handle.flush()
                    os.fsync(handle.fileno())
            except BaseException:
                Path(temporary_name).unlink(missing_ok=True)
                raise
            staged[final_path] = Path(temporary_name)

        if not overwrite and any(path.exists() for path in paths):
            raise FileExistsError("refusing fixture files created during staging")
        if overwrite:
            for final_path in paths:
                if final_path.exists():
                    descriptor, backup_name = tempfile.mkstemp(
                        prefix=f".{final_path.name}.", suffix=".backup", dir=destination
                    )
                    os.close(descriptor)
                    backup = Path(backup_name)
                    backup.unlink()
                    replace_file(final_path, backup)
                    backups[final_path] = backup
        for final_path in paths:
            replace_file(staged[final_path], final_path)
            staged.pop(final_path)
            installed.append(final_path)
        committed = True
    except BaseException as original_error:
        rollback_errors = []
        for final_path in installed:
            try:
                final_path.unlink(missing_ok=True)
            except OSError as error:
                rollback_errors.append(f"could not remove {final_path}: {error}")
        for final_path, backup in tuple(backups.items()):
            if backup.exists():
                try:
                    replace_file(backup, final_path)
                    backups.pop(final_path)
                except OSError as error:
                    rollback_errors.append(
                        f"retained recoverable backup {backup} for {final_path}: {error}"
                    )
        if rollback_errors:
            retained = [str(path) for path in backups.values() if path.exists()]
            raise RuntimeError(
                "SYNTHETIC fixture rollback was incomplete; do not delete retained backup(s): "
                + ", ".join(retained) + "; " + "; ".join(rollback_errors)
            ) from original_error
        raise
    finally:
        for temporary in staged.values():
            temporary.unlink(missing_ok=True)
        if committed:
            for backup in backups.values():
                backup.unlink(missing_ok=True)
    return paths


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True, type=Path,
                        help="user-supplied test directory whose path contains SYNTHETIC")
    parser.add_argument("--overwrite", action="store_true",
                        help="replace only the two named SYNTHETIC fixture FCS files")
    args = parser.parse_args()
    paths = write_fixture(args.output_dir, overwrite=args.overwrite)
    print("Created SYNTHETIC test fixture (not experimental data):")
    for path in paths:
        print(path)


if __name__ == "__main__":
    main()
