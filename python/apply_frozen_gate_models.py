#!/usr/bin/env python3
"""Apply immutable development-only Single Cells and G1 rules to raw FCS events.

The caller supplies explicit event identities, channel-role mappings, artifact
paths and both artifact hashes.  This program never fits or adjusts a model.
It writes a new, previously absent output directory containing the three CSV
populations expected by an EdU report and an audit manifest.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import stat
import tempfile
import warnings
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree


FEATURE_SCHEMA_VERSION = "fd-feature-v2"
FEATURES = (
    "knn25_log_density_2d__fsc_ssc",
    "asinh150__dna_area",
    "asinh150__dna_height",
    "log_ratio__dna_area_height",
    "linear_tls_signed_distance__dna_height_area",
    "knn25_log_density_2d__dna_height_area",
)
REQUIRED_ROLES = ("fsc_area", "ssc_area", "dna_area", "dna_height", "edu")
IDENTITY_COLUMNS = ("acquisition_id", "event_identity", "event_index")


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_regular_bytes(path: Path, label: str) -> bytes:
    """Read one non-symlink regular file once for hash-and-parse consistency."""
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0)
    if not hasattr(os, "O_NOFOLLOW"):
        raise ValueError("platform cannot enforce no-follow model reads")
    try:
        descriptor = os.open(path, flags | os.O_NOFOLLOW)
    except OSError as exc:
        raise ValueError(f"{label} is missing or unsafe: {path}") from exc
    try:
        if not stat.S_ISREG(os.fstat(descriptor).st_mode):
            raise ValueError(f"{label} is not a regular file: {path}")
        with os.fdopen(descriptor, "rb") as handle:
            descriptor = -1
            return handle.read()
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def semantic_hash(rule: dict) -> str:
    payload = {key: value for key, value in rule.items() if key not in {
        "sha256", "solver", "selected_l2", "solver_iterations"
    }}
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def load_rule(path: Path, expected_byte_sha256: str,
              expected_semantic_sha256: str, expected_target: str) -> dict:
    content = read_regular_bytes(path, "model artifact")
    observed_byte = hashlib.sha256(content).hexdigest()
    if observed_byte != expected_byte_sha256:
        raise ValueError(f"model artifact byte SHA-256 mismatch: {path}")
    rule = json.loads(content.decode("utf-8"))
    required = {"feature_names", "center", "scale", "coefficients", "intercept",
                "threshold", "target", "schema_version", "sha256"}
    missing = sorted(required - set(rule))
    if missing:
        raise ValueError(f"model artifact is missing fields: {missing}")
    if rule["target"] != expected_target:
        raise ValueError(f"expected target {expected_target!r}, found {rule['target']!r}")
    if rule["schema_version"] != FEATURE_SCHEMA_VERSION:
        raise ValueError("unsupported model feature schema")
    if tuple(rule["feature_names"]) != FEATURES:
        raise ValueError("model feature names or order differ from the frozen schema")
    observed_semantic = semantic_hash(rule)
    if rule["sha256"] != observed_semantic or observed_semantic != expected_semantic_sha256:
        raise ValueError("model artifact semantic SHA-256 mismatch")
    sizes = [len(rule[key]) for key in ("center", "scale", "coefficients")]
    if sizes != [len(FEATURES)] * 3:
        raise ValueError("model parameter dimensions differ from the frozen schema")
    numeric = np.asarray(rule["center"] + rule["scale"] + rule["coefficients"] +
                         [rule["intercept"], rule["threshold"]], dtype=float)
    if not np.isfinite(numeric).all() or (np.asarray(rule["scale"]) <= 0).any():
        raise ValueError("model artifact contains invalid numeric parameters")
    if not 0 < float(rule["threshold"]) < 1:
        raise ValueError("model threshold must be strictly between zero and one")
    return rule


def robust_z(values: np.ndarray) -> np.ndarray:
    median = float(np.median(values))
    scale = float(np.median(np.abs(values - median))) * 1.4826
    if not np.isfinite(scale) or scale <= 1e-12:
        scale = 1.0
    return (values - median) / scale


def knn_log_density(points: np.ndarray) -> np.ndarray:
    if len(points) < 3:
        raise ValueError("at least three events are required for density features")
    k = min(25, len(points) - 1)
    distances, _ = cKDTree(points).query(points, k=k + 1, workers=1)
    radius = np.maximum(distances[:, -1], 1e-12)
    return np.log((k / len(points)) / np.maximum(np.pi * radius ** 2, 1e-24))


def extract_frozen_features(events: pd.DataFrame, roles: dict[str, str]) -> pd.DataFrame:
    missing_roles = sorted(set(REQUIRED_ROLES) - set(roles))
    if missing_roles:
        raise ValueError(f"unresolved required channel roles: {missing_roles}")
    if len(set(roles.values())) != len(roles):
        raise ValueError("channel roles must map to distinct columns")
    missing_columns = sorted(set(roles.values()) - set(events.columns))
    if missing_columns:
        raise ValueError(f"mapped channel columns are missing: {missing_columns}")
    raw = {}
    for role in REQUIRED_ROLES:
        values = pd.to_numeric(events[roles[role]], errors="coerce").to_numpy(float)
        if not np.isfinite(values).all():
            raise ValueError(f"channel role {role!r} contains non-finite values")
        raw[role] = values
    fsc_ssc = np.column_stack((robust_z(raw["fsc_area"]), robust_z(raw["ssc_area"])))
    dna_pulse = np.column_stack((robust_z(raw["dna_height"]), robust_z(raw["dna_area"])))
    covariance = np.cov(dna_pulse, rowvar=False, ddof=0)
    values, vectors = np.linalg.eigh(covariance)
    direction = vectors[:, int(np.argmax(values))]
    if direction[0] < 0 or (direction[0] == 0 and direction[1] < 0):
        direction = -direction
    perpendicular = np.array((-direction[1], direction[0]))
    result = pd.DataFrame({
        FEATURES[0]: knn_log_density(fsc_ssc),
        FEATURES[1]: np.arcsinh(raw["dna_area"] / 150.0),
        FEATURES[2]: np.arcsinh(raw["dna_height"] / 150.0),
        FEATURES[3]: np.log((np.abs(raw["dna_area"]) + 1.0) /
                            (np.abs(raw["dna_height"]) + 1.0)),
        FEATURES[4]: dna_pulse @ perpendicular,
        FEATURES[5]: knn_log_density(dna_pulse),
    }, index=events.index)
    if not np.isfinite(result.to_numpy()).all():
        raise ValueError("feature extraction produced non-finite values")
    return result


def predict(features: pd.DataFrame, rule: dict) -> tuple[np.ndarray, np.ndarray]:
    matrix = features.loc[:, rule["feature_names"]].to_numpy(float)
    standardized = (matrix - np.asarray(rule["center"])) / np.asarray(rule["scale"])
    logits = np.clip(float(rule["intercept"]) +
                     standardized @ np.asarray(rule["coefficients"]), -35, 35)
    probability = 1.0 / (1.0 + np.exp(-logits))
    return probability, probability >= float(rule["threshold"])


def validate_identity(frame: pd.DataFrame, label: str) -> str:
    missing = sorted(set(IDENTITY_COLUMNS) - set(frame.columns))
    if missing:
        raise ValueError(f"{label} lacks direct event identity columns: {missing}")
    if frame["event_identity"].isna().any() or (frame["event_identity"].astype(str) == "").any():
        raise ValueError(f"{label} contains blank event identities")
    if frame["event_identity"].duplicated().any():
        raise ValueError(f"{label} contains duplicate event identities")
    event_index = pd.to_numeric(frame["event_index"], errors="coerce")
    if event_index.isna().any() or (event_index < 0).any() or (event_index % 1 != 0).any():
        raise ValueError(f"{label} contains invalid direct event indices")
    if event_index.duplicated().any():
        raise ValueError(f"{label} contains duplicate direct event indices")
    acquisition_ids = frame["acquisition_id"].drop_duplicates().astype(str).tolist()
    if len(acquisition_ids) != 1 or not acquisition_ids[0]:
        raise ValueError(f"{label} must contain exactly one explicit acquisition_id")
    acquisition_id = acquisition_ids[0]
    expected = [f"{acquisition_id}:event_index:{int(value)}" for value in event_index]
    if frame["event_identity"].astype(str).tolist() != expected:
        raise ValueError(f"{label} event_identity does not equal acquisition_id:event_index")
    return acquisition_id


def model_gate_tables(all_events: pd.DataFrame, edu_positive: pd.DataFrame,
                      roles: dict[str, str], single_rule: dict,
                      g1_rule: dict) -> tuple[dict[str, pd.DataFrame], dict]:
    all_acquisition = validate_identity(all_events, "all-events input")
    positive_acquisition = validate_identity(edu_positive, "EdU-positive input")
    if positive_acquisition != all_acquisition:
        raise ValueError("all-events and EdU-positive inputs name different acquisitions")
    if not set(edu_positive["event_identity"]).issubset(set(all_events["event_identity"])):
        raise ValueError("EdU-positive event identities are not a subset of all events")
    features = extract_frozen_features(all_events, roles)
    single_probability, single_mask = predict(features, single_rule)
    g1_probability, raw_g1_mask = predict(features, g1_rule)
    g1_mask = single_mask & raw_g1_mask
    scored = all_events.copy()
    scored["model_single_cells_probability"] = single_probability
    scored["model_g1_probability"] = g1_probability
    single = scored.loc[single_mask].copy()
    g1 = scored.loc[g1_mask].copy()
    positive = edu_positive.loc[
        edu_positive["event_identity"].isin(set(single["event_identity"]))
    ].copy()
    report_columns = {
        roles["dna_area"]: "DNA content",
        roles["edu"]: "EdU",
    }
    tables = {
        "single_cells": single.rename(columns=report_columns),
        "g1": g1.rename(columns=report_columns),
        "edu_positive": positive.rename(columns=report_columns),
    }
    audit = {
        "all_event_count": int(len(all_events)),
        "model_single_cells_count": int(len(single)),
        "model_g1_count": int(len(g1)),
        "source_edu_positive_count": int(len(edu_positive)),
        "model_parent_intersected_edu_positive_count": int(len(positive)),
        "g1_subset_single_cells": bool(set(g1.event_identity).issubset(set(single.event_identity))),
        "edu_positive_subset_single_cells": bool(set(positive.event_identity).issubset(set(single.event_identity))),
    }
    return tables, audit


def write_csv_absent(frame: pd.DataFrame, path: Path) -> None:
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.",
                                                  suffix=".partial", dir=path.parent)
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            descriptor = -1
            frame.to_csv(handle, index=False)
            handle.flush()
            os.fsync(handle.fileno())
        os.link(temporary, path)
        temporary.unlink()
    except BaseException:
        if descriptor >= 0:
            os.close(descriptor)
        temporary.unlink(missing_ok=True)
        raise


def validate_output_dir(output_dir: Path, repository_root: Path) -> Path:
    if not output_dir.is_absolute():
        raise ValueError("output_dir must be an explicit absolute path")
    resolved_output = output_dir.resolve(strict=False)
    resolved_repository = repository_root.resolve(strict=True)
    if resolved_output == resolved_repository or resolved_repository in resolved_output.parents:
        raise ValueError("output_dir must be outside the source repository")
    if resolved_output.exists():
        raise ValueError(f"refusing existing output destination: {resolved_output}")
    return resolved_output


def write_json_absent(value: dict, path: Path) -> None:
    payload = (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".partial", dir=path.parent
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            descriptor = -1
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        os.link(temporary, path)
        temporary.unlink()
    except BaseException:
        if descriptor >= 0:
            os.close(descriptor)
        temporary.unlink(missing_ok=True)
        raise


def normalize_flowkit_raw_columns(frame: pd.DataFrame) -> pd.DataFrame:
    """Use the exact PnN/detector field from FlowKit raw-data columns."""
    normalized = frame.copy()
    if isinstance(normalized.columns, pd.MultiIndex):
        if normalized.columns.nlevels != 2:
            raise ValueError("FlowKit raw columns must have exactly two PnN/PnS levels")
        detector_names = []
        for column in normalized.columns.to_list():
            detector = column[0] if len(column) else None
            if not isinstance(detector, str) or not detector:
                raise ValueError("FlowKit raw column lacks an explicit PnN detector name")
            detector_names.append(detector)
        if len(detector_names) != len(set(detector_names)):
            raise ValueError("FlowKit raw PnN detector names are not unique")
        normalized.columns = detector_names
    elif not all(isinstance(column, str) and column for column in normalized.columns):
        raise ValueError("FlowKit raw columns must be explicit detector-name strings")
    if len(normalized.columns) != len(set(normalized.columns)):
        raise ValueError("FlowKit raw detector-name columns are not unique")
    return normalized


def canonical_source_indices(index: pd.Index, label: str) -> np.ndarray:
    values = pd.to_numeric(index, errors="coerce")
    if np.isnan(values).any() or (values < 0).any() or (values % 1 != 0).any():
        raise ValueError(f"{label} lacks canonical nonnegative integer source indices")
    values = values.astype(np.int64)
    if len(values) != len(set(values.tolist())):
        raise ValueError(f"{label} contains duplicate source indices")
    return values


def positive_from_gate_indices(all_events: pd.DataFrame,
                               gate_index: pd.Index) -> pd.DataFrame:
    """Select gate members from the canonical all-event measurements by index."""
    all_indices = canonical_source_indices(all_events.index, "all-events table")
    gate_indices = canonical_source_indices(gate_index, "EdU-positive gate")
    missing = sorted(set(gate_indices.tolist()) - set(all_indices.tolist()))
    if missing:
        raise ValueError("EdU-positive gate source indices are outside the all-events table")
    indexed = all_events.copy()
    indexed.index = all_indices
    return indexed.loc[gate_indices].copy()


def validate_unloaded_sample_warnings(caught: list[warnings.WarningMessage],
                                      loaded_sample: str) -> list[str]:
    """Accept only FlowKit's expected warnings for deliberately unloaded peers."""
    messages = []
    pattern = re.compile(r"^WSP references (.+), but sample was not loaded\.$")
    for item in caught:
        message = str(item.message)
        match = pattern.fullmatch(message)
        peer = match.group(1).strip() if match else ""
        peer_identity = peer.strip("'\"")
        if (item.category is not UserWarning or not match or not peer_identity or
                peer_identity.casefold() == loaded_sample.casefold()):
            raise RuntimeError(f"unexpected FlowKit warning: {message}")
        messages.append(message)
    return messages


def raw_fcs_inputs(item: dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Read one exact FCS and its explicitly named FlowJo EdU-positive gate."""
    try:
        import flowkit as fk
        from lxml import etree
    except ImportError as exc:
        raise RuntimeError("raw-FCS gating requires flowkit and lxml") from exc
    fcs_path, workspace_path = Path(item["fcs_path"]), Path(item["workspace_path"])
    for path, key in ((fcs_path, "fcs_sha256"), (workspace_path, "workspace_sha256")):
        if not path.is_file():
            raise ValueError(f"required immutable input is missing: {path}")
        if sha256_path(path) != item[key]:
            raise ValueError(f"immutable input SHA-256 mismatch: {path}")
    acquisition_id = item.get("acquisition_id")
    if not isinstance(acquisition_id, str) or not acquisition_id:
        raise ValueError("each raw-FCS acquisition requires an explicit acquisition_id")
    population = item.get("edu_positive_population")
    if not isinstance(population, str) or not population:
        raise ValueError("each acquisition requires an explicit EdU-positive population")
    with tempfile.TemporaryDirectory(prefix="model_gate_wsp_") as directory:
        adapted = Path(directory) / workspace_path.name
        tree = etree.parse(str(workspace_path))
        attribute = "{http://www.isac-net.org/std/Gating-ML/v2.0/transformations}minRange"
        for element in tree.xpath(
                "//transforms:linear",
                namespaces={"transforms": "http://www.isac-net.org/std/Gating-ML/v2.0/transformations"}):
            value = float(element.get(attribute, "0"))
            if value > 0:
                element.set(attribute, str(-value))
        tree.write(str(adapted), encoding="UTF-8", xml_declaration=True)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            workspace = fk.Workspace(str(adapted), fcs_samples=[str(fcs_path)],
                                     filename_as_id=True)
            sample_ids = workspace.get_sample_ids()
            if sample_ids != [fcs_path.name]:
                raise ValueError("workspace did not resolve exactly the explicitly mapped FCS")
            workspace.analyze_samples(use_mp=False)
        unloaded_warnings = validate_unloaded_sample_warnings(caught, fcs_path.name)
        paths = workspace.find_matching_gate_paths(fcs_path.name, population)
        if len(paths) != 1:
            raise ValueError("EdU-positive population is missing or ambiguous")
        gate_path = paths[0]
        sample = workspace.get_sample(fcs_path.name)
        all_events = normalize_flowkit_raw_columns(
            sample.as_dataframe(source="raw")
        )
        gate_events = workspace.get_gate_events(
            fcs_path.name, population, gate_path, source="raw"
        )
        positive = positive_from_gate_indices(all_events, gate_events.index)
    for frame in (all_events, positive):
        event_index = canonical_source_indices(frame.index, "FlowKit raw table")
        frame.insert(0, "event_index", event_index)
        frame.insert(0, "event_identity", [
            f"{acquisition_id}:event_index:{value}" for value in event_index
        ])
        frame.insert(0, "acquisition_id", acquisition_id)
        frame.reset_index(drop=True, inplace=True)
    source = {
        "source_fcs": str(fcs_path.resolve()), "source_fcs_sha256": item["fcs_sha256"],
        "source_workspace": str(workspace_path.resolve()),
        "source_workspace_sha256": item["workspace_sha256"],
        "edu_positive_population": population,
        "edu_positive_gate_path": "/".join((*gate_path, population)),
        "validated_unloaded_workspace_sample_warnings": unloaded_warnings,
    }
    return all_events, positive, source


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path, help="explicit model-gating JSON config")
    args = parser.parse_args()
    config = json.loads(args.config.read_text(encoding="utf-8"))
    if config.get("status_label") != "EXPERIMENTAL MODEL-DERIVED NON-PRODUCTION":
        parser.error("config must carry the exact experimental/non-production status label")
    try:
        output_dir = validate_output_dir(
            Path(config["output_dir"]), Path(__file__).resolve().parents[1]
        )
    except ValueError as exc:
        parser.error(str(exc))
    acquisitions = config.get("acquisitions")
    if not isinstance(acquisitions, list) or not acquisitions:
        parser.error("config must contain a nonempty acquisitions list")
    prefixes = [item.get("prefix") for item in acquisitions]
    if any(not isinstance(x, str) or not x or Path(x).name != x for x in prefixes):
        parser.error("each prefix must be an explicit path-free string")
    if len(prefixes) != len(set(prefixes)):
        parser.error("acquisition prefixes must be unique")
    artifacts = config["artifacts"]
    single_rule = load_rule(Path(artifacts["single_cells"]["path"]),
                            artifacts["single_cells"]["byte_sha256"],
                            artifacts["single_cells"]["semantic_sha256"],
                            "single_cells_reference")
    g1_rule = load_rule(Path(artifacts["g1"]["path"]),
                        artifacts["g1"]["byte_sha256"],
                        artifacts["g1"]["semantic_sha256"], "g1_reference")
    output_dir.mkdir(parents=True)
    manifest = {"status_label": config["status_label"], "feature_schema_version":
                FEATURE_SCHEMA_VERSION, "single_cells_source": "model",
                "g1_source": "model", "artifacts": artifacts, "acquisitions": []}
    try:
        for item in acquisitions:
            all_events, positive, source = raw_fcs_inputs(item)
            tables, audit = model_gate_tables(all_events, positive, item["channel_roles"],
                                              single_rule, g1_rule)
            records = {}
            for key, suffix in (("single_cells", "_single_cells.csv"),
                                ("g1", "_g1.csv"),
                                ("edu_positive", "_edu_positive.csv")):
                path = output_dir / f"{item['prefix']}{suffix}"
                write_csv_absent(tables[key], path)
                records[key] = {"path": path.name, "sha256": sha256_path(path),
                                "rows": int(len(tables[key]))}
            manifest["acquisitions"].append({
                "prefix": item["prefix"], "acquisition_id": item["acquisition_id"],
                **source,
                "channel_roles": item["channel_roles"], "audit": audit,
                "outputs": records,
            })
        manifest_path = output_dir / "model-gating-manifest.json"
        write_json_absent(manifest, manifest_path)
        return 0
    except BaseException:
        # Keep any partial operation visible for forensic inspection; never silently retry.
        raise


if __name__ == "__main__":
    raise SystemExit(main())
