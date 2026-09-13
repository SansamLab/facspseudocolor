#!/usr/bin/env python3
"""Apply the approved immutable experimental Single Cells, G1, and EdU rules.

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
import math
import re
import tempfile
import warnings
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree


FEATURE_SCHEMA_VERSION = "fd-feature-v2"
PROFILE = "single_g1_edu_frozen_v1"
STATUS = "EXPERIMENTAL MODEL-DERIVED NON-PRODUCTION"
FROZEN_THRESHOLDS = {"single_cells": 0.205, "g1": 0.255, "edu_positive": 0.385}
G1_DNA_MINIMUM_FRACTION = 0.35
MAX_ALL_EVENT_ROWS = 5_000_000
MAX_ALL_EVENT_BYTES = 512 * 1024 * 1024
EDU_PORTABLE_SHA256 = "9e6ed466687efe5f7523dea633e9e9244381c0e613c86f0944f861c06047fdf4"
EDU_SOURCE_PICKLE_SHA256 = "6ef5603555a660b8503379cbbcab131b616336a64330ffd42f45dfaa182042cc"
PORTABLE_PREDICTOR_REFERENCE_SHA256 = "2edf1bbd529159e8752d01731e87f19efbc63b514409b1d246d7a4ffc266a2c3"
FROZEN_RULE_HASHES = {
    "single_cells": {
        "byte_sha256": "71b459d8fb00930f31a6d289a21f587226fd2d6be4b31ebc782b7bae7839a376",
        "semantic_sha256": "f9a53e9fd78d2f39cf5980b8f470c3c44e7a187fa942cbcc0324c563fa68a31a",
    },
    "g1": {
        "byte_sha256": "d653d388d15af3bd22185cb9c8e1addea132e0b5d1d18e60e37c65a7548163b7",
        "semantic_sha256": "32723ccc768ed517affa040f2ce62a125dad465399d869f471bc331549f47f6a",
    },
}
EDU_FEATURES = (
    "dna_area_q", "dna_pulse_q", "dna_tls_position", "dna_tls_distance",
    "dna_area_pulse_log_density", "edu_area_q", "dna_edu_log_density",
)
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


def require_regular_file(path: Path, label: str) -> None:
    if path.is_symlink() or not path.is_file():
        raise ValueError(f"{label} must be an existing regular non-symlink file: {path}")


def validate_all_event_export_size(path: Path) -> None:
    if path.stat().st_size > MAX_ALL_EVENT_BYTES:
        raise ValueError("all-events export exceeds the 512 MiB safety limit")


def semantic_hash(rule: dict) -> str:
    payload = {key: value for key, value in rule.items() if key not in {
        "sha256", "solver", "selected_l2", "solver_iterations"
    }}
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def load_rule(path: Path, expected_byte_sha256: str,
              expected_semantic_sha256: str, expected_target: str,
              expected_threshold: float) -> dict:
    require_regular_file(path, "model artifact")
    observed_byte = sha256_path(path)
    if observed_byte != expected_byte_sha256:
        raise ValueError(f"model artifact byte SHA-256 mismatch: {path}")
    rule = json.loads(path.read_text(encoding="utf-8"))
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
    if not np.isclose(float(rule["threshold"]), expected_threshold,
                      rtol=0.0, atol=1e-12):
        raise ValueError(f"{expected_target} threshold differs from frozen profile")
    return rule


def encoded_float(value: Any) -> float:
    if not isinstance(value, str):
        raise ValueError("portable model encoded float is not a string")
    result = float.fromhex(value)
    if not math.isfinite(result):
        raise ValueError("portable model contains a non-finite float")
    return result


def load_edu_model(model_path: Path) -> dict:
    """Load and validate the exact approved PORTABLE_HGB_V1 EdU model."""
    require_regular_file(model_path, "EdU model artifact")
    if sha256_path(model_path) != EDU_PORTABLE_SHA256:
        raise ValueError("EdU model artifact byte SHA-256 mismatch")
    model = json.loads(model_path.read_text(encoding="utf-8"))
    required = {
        "format": "PORTABLE_HGB_V1", "model_id": "EDU_POSITIVE_MODEL002",
        "source_pickle_sha256": EDU_SOURCE_PICKLE_SHA256,
        "feature_order": list(EDU_FEATURES), "winner": "hgb_small",
        "task": "binary_classification", "link": "logit",
        "decision_rule": "positive_probability >= decision_threshold",
        "authoritative_decimal_decision_threshold": "0.385",
    }
    for key, value in required.items():
        if model.get(key) != value:
            raise ValueError(f"frozen portable EdU model field mismatch: {key}")
    if model.get("n_features_in") != len(EDU_FEATURES) or model.get("classes") != [0, 1]:
        raise ValueError("portable EdU model dimensions or classes differ")
    if not np.isclose(encoded_float(model.get("decision_threshold")),
                      FROZEN_THRESHOLDS["edu_positive"], rtol=0.0, atol=1e-15):
        raise ValueError("portable EdU model threshold differs from frozen profile")
    encoded_float(model.get("baseline_raw"))
    iterations = model.get("iterations")
    if not isinstance(iterations, list) or len(iterations) != model.get("n_iterations"):
        raise ValueError("portable EdU model iteration count differs")
    for iteration in iterations:
        if not isinstance(iteration, list) or len(iteration) != 1:
            raise ValueError("portable EdU binary iteration must contain one tree")
        nodes = iteration[0].get("nodes")
        if not isinstance(nodes, list) or not nodes:
            raise ValueError("portable EdU tree has no nodes")
        for index, node in enumerate(nodes):
            if node.get("index") != index:
                raise ValueError("portable EdU tree node indexes are not contiguous")
            encoded_float(node.get("value"))
            encoded_float(node.get("num_threshold"))
    return model


def portable_tree_value(tree: dict, row: list[float], known: list[set[float]]) -> float:
    nodes, index = tree["nodes"], 0
    while True:
        node = nodes[index]
        if node["is_leaf"]:
            return encoded_float(node["value"])
        feature, value = node["feature_idx"], row[node["feature_idx"]]
        if math.isnan(value):
            go_left = node["missing_go_to_left"]
        elif node["is_categorical"]:
            integer = int(value)
            if value < 0:
                go_left = node["missing_go_to_left"]
            elif integer == value and 0 <= integer < 256:
                words = tree["raw_left_cat_bitsets"][node["bitset_idx"]]
                if bool((words[integer >> 5] >> (integer & 31)) & 1):
                    go_left = True
                elif value in known[feature]:
                    go_left = False
                else:
                    go_left = node["missing_go_to_left"]
            elif value in known[feature]:
                go_left = False
            else:
                go_left = node["missing_go_to_left"]
        else:
            go_left = value <= encoded_float(node["num_threshold"])
        index = node["left"] if go_left else node["right"]


def portable_predict_positive(model: dict, rows: np.ndarray) -> np.ndarray:
    if rows.ndim != 2 or rows.shape[1] != len(EDU_FEATURES) or not np.isfinite(rows).all():
        raise ValueError("portable EdU input differs from frozen feature schema")
    known = [set(encoded_float(value) for value in values)
             for values in model["bin_mapper"]["known_categories"]]
    probabilities = []
    for row in rows:
        raw = encoded_float(model["baseline_raw"])
        converted = [float(value) for value in row]
        for iteration in model["iterations"]:
            raw += portable_tree_value(iteration[0], converted, known)
        probabilities.append(1.0 / (1.0 + math.exp(-raw)) if raw >= 0 else
                             math.exp(raw) / (1.0 + math.exp(raw)))
    return np.asarray(probabilities, dtype=float)


def validate_frozen_artifact_config(artifacts: Any) -> None:
    """Reject any artifact provenance outside the approved frozen profile."""
    if not isinstance(artifacts, dict) or set(artifacts) != {
            "single_cells", "g1", "edu_positive"}:
        raise ValueError("profile requires exactly three frozen model artifacts")
    for target in ("single_cells", "g1"):
        if not isinstance(artifacts[target], dict) or set(artifacts[target]) != {
                "path", "byte_sha256", "semantic_sha256"}:
            raise ValueError(f"{target} artifact must contain exact frozen provenance fields")
        if any(artifacts[target][field] != expected
               for field, expected in FROZEN_RULE_HASHES[target].items()):
            raise ValueError(f"{target} artifact differs from the frozen profile")
    edu = artifacts["edu_positive"]
    if not isinstance(edu, dict) or set(edu) != {"path", "byte_sha256"}:
        raise ValueError("edu_positive artifact must contain exact portable provenance fields")
    if edu["byte_sha256"] != EDU_PORTABLE_SHA256:
        raise ValueError("edu_positive artifact hash differs from the frozen profile")


def edu_features(events: pd.DataFrame, roles: dict[str, str]) -> np.ndarray:
    columns = [roles[key] for key in ("dna_area", "dna_height", "edu")]
    raw = events.loc[:, columns].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    if not np.isfinite(raw).all():
        raise ValueError("EdU model channels contain non-finite values")
    ranked = [pd.Series(raw[:, i]).rank(method="average", pct=True).to_numpy(np.float32)
              for i in range(3)]
    dna_area, dna_pulse, edu_area = ranked
    centered = np.c_[dna_area - dna_area.mean(), dna_pulse - dna_pulse.mean()]
    _, _, vt = np.linalg.svd(centered, full_matrices=False)
    axis = vt[0] * (1 if vt[0, 0] >= 0 else -1)
    def density(x: np.ndarray, y: np.ndarray) -> np.ndarray:
        histogram, _, _ = np.histogram2d(x, y, bins=64, range=((0, 1), (0, 1)))
        i = np.minimum((x * 64).astype(int), 63)
        j = np.minimum((y * 64).astype(int), 63)
        return (np.log1p(histogram[i, j]) / np.log1p(len(x))).astype(np.float32)
    result = np.c_[dna_area, dna_pulse, centered @ axis,
                   centered @ np.array([-axis[1], axis[0]]),
                   density(dna_area, dna_pulse), edu_area,
                   density(dna_area, edu_area)].astype(np.float32)
    if result.shape != (len(events), len(EDU_FEATURES)) or not np.isfinite(result).all():
        raise ValueError("EdU feature construction failed")
    return result


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


def model_gate_tables(all_events: pd.DataFrame, roles: dict[str, str],
                      single_rule: dict, g1_rule: dict,
                      edu_model: Any) -> tuple[dict[str, pd.DataFrame], dict]:
    all_acquisition = validate_identity(all_events, "all-events input")
    features = extract_frozen_features(all_events, roles)
    single_probability, single_mask = predict(features, single_rule)
    g1_probability, raw_g1_mask = predict(features, g1_rule)
    scored = all_events.copy()
    scored["model_single_cells_probability"] = single_probability
    scored["model_g1_probability"] = g1_probability
    edu_probability = portable_predict_positive(edu_model, edu_features(all_events, roles))
    if edu_probability.shape != (len(all_events),) or not np.isfinite(edu_probability).all():
        raise ValueError("EdU model returned invalid probabilities")
    scored["model_edu_positive_probability"] = edu_probability
    positive_mask = single_mask & (edu_probability >= FROZEN_THRESHOLDS["edu_positive"])
    # In EdU experiments, G1 is the non-replicating reference population.
    # Exclude EdU-positive events before this table is used for normalization.
    preliminary_g1_mask = single_mask & raw_g1_mask & ~positive_mask
    dna_area = pd.to_numeric(all_events[roles["dna_area"]], errors="coerce").to_numpy()
    positive_candidate_dna = dna_area[
        preliminary_g1_mask & np.isfinite(dna_area) & (dna_area > 0)
    ]
    if not len(positive_candidate_dna):
        raise ValueError("preliminary EdU-negative G1 population has no positive DNA-A values")
    g1_dna_center = float(np.median(positive_candidate_dna))
    g1_dna_minimum = G1_DNA_MINIMUM_FRACTION * g1_dna_center
    if not np.isfinite(g1_dna_center) or g1_dna_center <= 0:
        raise ValueError("preliminary EdU-negative G1 DNA-A center is invalid")
    g1_mask = (preliminary_g1_mask & np.isfinite(dna_area) &
               (dna_area >= g1_dna_minimum))
    single = scored.loc[single_mask].copy()
    g1 = scored.loc[g1_mask].copy()
    positive = scored.loc[positive_mask].copy()
    report_columns = {
        roles["dna_area"]: "DNA content",
        roles["edu"]: "EdU",
    }
    tables = {
        "all_events": scored.rename(columns=report_columns),
        "single_cells": single.rename(columns=report_columns),
        "g1": g1.rename(columns=report_columns),
        "edu_positive": positive.rename(columns=report_columns),
    }
    audit = {
        "all_event_count": int(len(all_events)),
        "model_single_cells_count": int(len(single)),
        "model_g1_count": int(len(g1)),
        "model_edu_positive_count": int(len(positive)),
        "g1_subset_single_cells": bool(set(g1.event_identity).issubset(set(single.event_identity))),
        "edu_positive_subset_single_cells": bool(set(positive.event_identity).issubset(set(single.event_identity))),
        "g1_excludes_edu_positive": bool(set(g1.event_identity).isdisjoint(set(positive.event_identity))),
        "preliminary_g1_count": int(preliminary_g1_mask.sum()),
        "g1_dna_positive_candidate_count": int(len(positive_candidate_dna)),
        "g1_dna_2n_center_raw": g1_dna_center,
        "g1_dna_minimum_fraction": G1_DNA_MINIMUM_FRACTION,
        "g1_dna_minimum_raw": g1_dna_minimum,
        "g1_dna_minimum_excluded_count": int(preliminary_g1_mask.sum() - g1_mask.sum()),
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


def raw_fcs_inputs(item: dict[str, Any]) -> tuple[pd.DataFrame, dict]:
    """Read one exact immutable FCS with no workspace or FlowJo dependency."""
    try:
        import flowkit as fk
    except ImportError as exc:
        raise RuntimeError("raw-FCS model gating requires flowkit") from exc
    fcs_path = Path(item["fcs_path"])
    require_regular_file(fcs_path, "immutable FCS")
    if sha256_path(fcs_path) != item["fcs_sha256"]:
        raise ValueError(f"immutable FCS is missing or has a SHA-256 mismatch: {fcs_path}")
    acquisition_id = item.get("acquisition_id")
    if not isinstance(acquisition_id, str) or not acquisition_id:
        raise ValueError("each raw-FCS acquisition requires an explicit acquisition_id")
    sample = fk.Sample(str(fcs_path))
    all_events = normalize_flowkit_raw_columns(sample.as_dataframe(source="raw"))
    event_index = canonical_source_indices(all_events.index, "FlowKit raw table")
    all_events.insert(0, "event_index", event_index)
    all_events.insert(0, "event_identity", [
        f"{acquisition_id}:event_index:{value}" for value in event_index])
    all_events.insert(0, "acquisition_id", acquisition_id)
    all_events.reset_index(drop=True, inplace=True)
    source = {
        "source_fcs": str(fcs_path.resolve()), "source_fcs_sha256": item["fcs_sha256"],
        "workspace_used": False,
    }
    return all_events, source


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path, help="explicit model-gating JSON config")
    args = parser.parse_args()
    if args.config.suffix.lower() in (".yml", ".yaml"):
        try:
            import yaml
        except ImportError as exc:
            raise RuntimeError("YAML configuration requires PyYAML") from exc
        root_config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
        config = root_config.get("gating", {}) if isinstance(root_config, dict) else {}
    else:
        config = json.loads(args.config.read_text(encoding="utf-8"))
    if config.get("mode") != "model_experimental" or config.get("profile") != PROFILE:
        parser.error(f"gating must select mode=model_experimental and profile={PROFILE}")
    if config.get("status_label") != STATUS:
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
    for item in acquisitions:
        if set(item) != {"prefix", "acquisition_id", "fcs_path", "fcs_sha256",
                         "channel_roles"}:
            parser.error("each acquisition must contain only the exact required fields")
        if not Path(item["fcs_path"]).is_absolute():
            parser.error("each fcs_path must be absolute")
        if re.fullmatch(r"[0-9a-f]{64}", item["fcs_sha256"]) is None:
            parser.error("each fcs_sha256 must be lowercase hexadecimal SHA-256")
        if tuple(item["channel_roles"]) != REQUIRED_ROLES or len(set(
                item["channel_roles"].values())) != len(REQUIRED_ROLES):
            parser.error("channel_roles must contain the five ordered, distinct frozen roles")
    artifacts = config["artifacts"]
    try:
        validate_frozen_artifact_config(artifacts)
    except ValueError as exc:
        parser.error(str(exc))
    single_rule = load_rule(Path(artifacts["single_cells"]["path"]),
                            artifacts["single_cells"]["byte_sha256"],
                            artifacts["single_cells"]["semantic_sha256"],
                            "single_cells_reference", FROZEN_THRESHOLDS["single_cells"])
    g1_rule = load_rule(Path(artifacts["g1"]["path"]),
                        artifacts["g1"]["byte_sha256"],
                        artifacts["g1"]["semantic_sha256"], "g1_reference",
                        FROZEN_THRESHOLDS["g1"])
    edu_model = load_edu_model(Path(artifacts["edu_positive"]["path"]))
    output_dir.mkdir(parents=True)
    manifest = {"status_label": config["status_label"], "gating_mode": config["mode"],
                "profile": config["profile"], "feature_schema_version":
                FEATURE_SCHEMA_VERSION, "frozen_thresholds": FROZEN_THRESHOLDS,
                "g1_dna_minimum": {
                    "method": "median_positive_dna_area_of_preliminary_edu_negative_g1",
                    "fraction": G1_DNA_MINIMUM_FRACTION,
                },
                "edu_predictor": {"format": "PORTABLE_HGB_V1",
                                  "implementation": "embedded dependency-free reference",
                                  "derived_from_reference_sha256":
                                  PORTABLE_PREDICTOR_REFERENCE_SHA256},
                "artifacts": artifacts,
                "acquisitions": []}
    try:
        for item in acquisitions:
            all_events, source = raw_fcs_inputs(item)
            if len(all_events) > MAX_ALL_EVENT_ROWS:
                raise ValueError("all-events export exceeds the five-million-row safety limit")
            tables, audit = model_gate_tables(all_events, item["channel_roles"],
                                              single_rule, g1_rule, edu_model)
            records = {}
            for key, suffix in (("all_events", "_all_events.csv"),
                                ("single_cells", "_single_cells.csv"),
                                ("g1", "_g1.csv"),
                                ("edu_positive", "_edu_positive.csv")):
                path = output_dir / f"{item['prefix']}{suffix}"
                write_csv_absent(tables[key], path)
                if key == "all_events":
                    validate_all_event_export_size(path)
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
