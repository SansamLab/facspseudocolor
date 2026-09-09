#!/usr/bin/env python3
"""Create model-only EdU report inputs from explicitly mapped immutable FCS files.

All three population calls are model-derived research outputs.  No workspace,
expert gate, or annotation is read.  G1 features are constructed inside the
predicted Single Cells parent; EdU features are constructed on the complete
source frame and scored only for that same parent.
"""

from __future__ import annotations

import argparse
import ctypes
import errno
import hashlib
import json
import os
import shutil
import stat
import sys
import types
import tempfile
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from apply_frozen_gate_models import (
    load_rule, extract_frozen_features, predict, normalize_flowkit_raw_columns,
    sha256_path, validate_output_dir, write_csv_absent, write_json_absent,
)
STATUS = "EXPERIMENTAL MODEL-DERIVED NON-PRODUCTION"
PARENT_METHOD = "model.single_cells.frozen_development_rule"
EXPECTED_SINGLE_THRESHOLD = 0.205
EXPECTED_SINGLE_BYTE_SHA256 = "71b459d8fb00930f31a6d289a21f587226fd2d6be4b31ebc782b7bae7839a376"
EXPECTED_SINGLE_SEMANTIC_SHA256 = "f9a53e9fd78d2f39cf5980b8f470c3c44e7a187fa942cbcc0324c563fa68a31a"
EXPECTED_SINGLE_SCHEMA = "fd-feature-v2"
REQUIRED_ROLES = {"fsc_area", "ssc_area", "dna_area", "dna_pulse", "edu_area"}
APPROVAL = "OWNER_AUTHORIZED_PREDICTED_PARENT_EXPERIMENTAL_BRANCH_ONLY"
MAPPING_SHA256 = "a874815c83816d6d4291054fceeb08fa0df4b0c864499176dafe3a6a400597fb"
G1_MODEL_ID = "DNA_PROJECTED_G1_CLUSTER001"
G1_THRESHOLD = 0.42
G1_FEATURES = ["dna_area_q", "dna_pulse_q", "tls_position",
               "tls_signed_distance", "tls_abs_distance",
               "dna_2d_log_density", "pulse_area_log_ratio"]
EDU_MODEL_ID = "EDU_POSITIVE_MODEL002"
EDU_THRESHOLD = 0.385
EDU_FEATURES = ["dna_area_q", "dna_pulse_q", "dna_tls_position",
                "dna_tls_distance", "dna_area_pulse_log_density",
                "edu_area_q", "dna_edu_log_density"]
ORIGINAL_MODEL_SHA256 = {
    G1_MODEL_ID: "9ce7e448455e2d42090a86c80beb10d925a7c0e64e6261e678335385d9714af6",
    EDU_MODEL_ID: "6ef5603555a660b8503379cbbcab131b616336a64330ffd42f45dfaa182042cc",
}
G1_METADATA = {
    "SUMMARY.json": "7a7446d61ca175e56964fc5ab6edbee081f2eb3f16e7925542157bc3719ddb10",
    "PREFLIGHT.json": "bac957a8837774afee0faba2ff2e832d003f1e37586a68e6af4037083448ec80",
    "per_unit_metrics.csv": "650076ecf722f165549731a0c84f1cbdde13d7cd36b53653da4a5b7b34b03792",
    "model_comparison.csv": "54f538ddb00b28a932181d5cd57dda0b3bfe37786a7f4d3357e0a71e9c5c67d2",
    "derived_label_summary.csv": "e60a3c24bd2ef5b80913d1cbb5387b74df0eef2ac4d67b1dfc5945162943380f",
    "SCOPE.csv": "264b0ab4ef806933fc03eb66f644ed4405fe77d196bee516a7b32a98acc6c608",
}
EDU_METADATA = {
    "SUMMARY.json": "dcd16784a80b68070eeb42126d6411db148646af6acc69ded7fb6a529b39d22f",
    "PREFLIGHT.json": "9bee03f4b039d7645140df43ec03cd586ab7f6d620c0f359646caa65f49914db",
}
EDU_CONFIG_SHA256 = "624a839a38a464976b214afcf8964a5f54009ef7694d16d71a91ba20677d44f3"
PORTABLE_PACKAGE_SHA256 = {
    "MANIFEST.json": "c3885b0c49497672ad83657e2cd5a08beb219201039f5b26cffa002a4dd4fd7f",
    "COMPATIBILITY_REPORT.json": "60695bb0a5d3d728437aa89aa570b9ee9d2c7be7cacfbbb564ab43b4309f3b38",
    "DNA_PROJECTED_G1_CLUSTER001.portable.json": "7fe256e892adfef6415e3957231bf8f68b7ad401dc909ffd2a801fbde8f4b75d",
    "EDU_POSITIVE_MODEL002.portable.json": "9e6ed466687efe5f7523dea633e9e9244381c0e613c86f0944f861c06047fdf4",
    "EXPORT_PROVENANCE.json": "e0b80bfc44faf26886f04f9093eac1d8ed6474f89a43a90a5a8184819bbd40f8",
    "NO_PROTECTED_DATA_PROOF.json": "f0b0b3aa0088c70e0f3e29cd75c5bc7e9fc8c809986b86f2ae73dbcb6301111b",
    "SYNTHETIC_EQUIVALENCE_FIXTURES.json": "ffdcdd99b295aa8d7ce946ccd901be1c6d15247f5c3c0726ef1c09eea05cddf4",
    "portable_hgb_predictor.py": "2edf1bbd529159e8752d01731e87f19efbc63b514409b1d246d7a4ffc266a2c3",
    "run_portable_selftest.py": "1d0584164fc36581fe0e2e958258b37ce2d72ec07f01cecbbbd7eabe54420fb4",
}


def canonical_json_sha256(value: object) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()


def read_regular_bytes(path: Path) -> bytes:
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0)
    if not hasattr(os, "O_NOFOLLOW"):
        raise ValueError("platform cannot enforce no-follow artifact reads")
    descriptor = os.open(path, flags | os.O_NOFOLLOW)
    try:
        if not stat.S_ISREG(os.fstat(descriptor).st_mode):
            raise ValueError(f"artifact is not a regular file: {path}")
        with os.fdopen(descriptor, "rb") as handle:
            descriptor = -1
            return handle.read()
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def load_approved_mapping(repository_root: Path) -> dict:
    path = repository_root / "inst/config/model_only_edu_figure6_mapping.json"
    content = read_regular_bytes(path)
    if hashlib.sha256(content).hexdigest() != MAPPING_SHA256:
        raise ValueError("immutable Figure6 run-mapping digest mismatch")
    mapping = json.loads(content)
    if (mapping.get("schema_version") != 1 or
            mapping.get("status") != "OWNER-APPROVED IMMUTABLE FIGURE6 MODEL-ONLY RUN MAPPING" or
            mapping.get("functional_panel") != "PANEL_FL2_DNA_FL4_EDU"):
        raise ValueError("immutable Figure6 run-mapping contract mismatch")
    roles = mapping.get("channel_roles")
    if not isinstance(roles, dict) or set(roles) != REQUIRED_ROLES:
        raise ValueError("immutable Figure6 mapping lacks exact channel roles")
    if any(not isinstance(value, str) or not value for value in roles.values()):
        raise ValueError("channel roles must be explicit nonblank raw parameter names")
    if len(set(roles.values())) != len(roles):
        raise ValueError("channel roles must resolve to distinct raw parameters")
    rows = mapping.get("acquisitions")
    if not isinstance(rows, list) or len(rows) != 8:
        raise ValueError("immutable Figure6 mapping must contain exactly eight acquisitions")
    return mapping


def validate_mapping(item: dict, approved: dict) -> tuple[dict[str, str], dict]:
    forbidden = {"workspace_path", "workspace_sha256", "edu_positive_population",
                 "channel_roles", "functional_panel", "mapping_sha256"}
    if forbidden & set(item):
        raise ValueError("run request must not contain workspace or self-attested mapping fields")
    matches = [x for x in approved["acquisitions"] if x["prefix"] == item.get("prefix")]
    if len(matches) != 1:
        raise ValueError("run acquisition is absent or ambiguous in immutable mapping")
    record = matches[0]
    for key in ("acquisition_id", "fcs_sha256"):
        if item.get(key) != record[key]:
            raise ValueError(f"run {key} differs from immutable Figure6 mapping")
    path = Path(item.get("fcs_path", ""))
    if path.name != record["fcs_basename"]:
        raise ValueError("FCS basename differs from immutable Figure6 mapping")
    return dict(approved["channel_roles"]), record


def snapshot_fcs(item: dict, directory: Path) -> tuple[Path, str]:
    source = Path(item.get("fcs_path", ""))
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0)
    if not hasattr(os, "O_NOFOLLOW"):
        raise ValueError("platform cannot enforce no-follow FCS reads")
    try:
        descriptor = os.open(source, flags | os.O_NOFOLLOW)
    except OSError as exc:
        raise ValueError(f"required immutable FCS is missing or unsafe: {source}") from exc
    snapshot = directory / source.name
    digest = hashlib.sha256()
    try:
        before = os.fstat(descriptor)
        if not stat.S_ISREG(before.st_mode):
            raise ValueError("FCS input is not a regular file")
        with os.fdopen(descriptor, "rb") as source_handle, snapshot.open("xb") as target:
            descriptor = -1
            for block in iter(lambda: source_handle.read(1 << 20), b""):
                digest.update(block)
                target.write(block)
            target.flush()
            os.fsync(target.fileno())
            after = os.fstat(source_handle.fileno())
        if ((before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns) !=
                (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns)):
            raise ValueError("FCS changed while creating immutable consumed snapshot")
    finally:
        if descriptor >= 0:
            os.close(descriptor)
    observed = digest.hexdigest()
    if observed != item.get("fcs_sha256"):
        raise ValueError("immutable FCS SHA-256 mismatch")
    return snapshot, observed


def read_raw_fcs(item: dict, roles: dict[str, str]) -> tuple[pd.DataFrame, str]:
    try:
        import flowkit as fk
    except ImportError as exc:
        raise RuntimeError("raw-FCS model inference requires flowkit") from exc
    with tempfile.TemporaryDirectory(prefix="model_only_fcs_snapshot_") as temp:
        snapshot, consumed_digest = snapshot_fcs(item, Path(temp))
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            sample = fk.Sample(str(snapshot))
            frame = normalize_flowkit_raw_columns(sample.as_dataframe(source="raw"))
    if caught:
        raise RuntimeError("unexpected warning while reading raw FCS: " +
                           " | ".join(str(x.message) for x in caught))
    missing = sorted(set(roles.values()) - set(frame.columns))
    if missing:
        raise ValueError(f"mapped raw channels are missing: {missing}")
    index = pd.to_numeric(frame.index, errors="coerce").to_numpy()
    if (not np.isfinite(index).all() or not np.array_equal(
            index, np.arange(len(frame), dtype=float))):
        raise ValueError("complete FCS frame event indices must be exact contiguous source order")
    index = index.astype(np.int64)
    acquisition = item["acquisition_id"]
    frame = frame.copy(deep=True)
    frame.insert(0, "event_index", index)
    frame.insert(0, "event_identity", [
        f"{acquisition}:event_index:{value}" for value in index
    ])
    frame.insert(0, "acquisition_id", acquisition)
    frame.reset_index(drop=True, inplace=True)
    return frame, consumed_digest


def average_percentile_rank(values: np.ndarray, dtype=float) -> np.ndarray:
    return pd.Series(values).rank(method="average", pct=True).to_numpy(dtype)


def density64(x: np.ndarray, y: np.ndarray, dtype=float) -> np.ndarray:
    histogram, _, _ = np.histogram2d(x, y, bins=64, range=((0, 1), (0, 1)))
    i = np.minimum((x * 64).astype(int), 63)
    j = np.minimum((y * 64).astype(int), 63)
    return (np.log1p(histogram[i, j]) / np.log1p(len(x))).astype(dtype)


def build_g1_features(parent: pd.DataFrame, roles: dict[str, str]) -> np.ndarray:
    area = parent[roles["dna_area"]].to_numpy(float)
    pulse = parent[roles["dna_pulse"]].to_numpy(float)
    if not len(area) or not np.isfinite(area).all() or not np.isfinite(pulse).all():
        raise ValueError("G1 parent DNA roles are empty or non-finite")
    area_q, pulse_q = average_percentile_rank(area), average_percentile_rank(pulse)
    centered = np.c_[area_q, pulse_q] - np.array([area_q.mean(), pulse_q.mean()])
    _, _, axes = np.linalg.svd(centered, full_matrices=False)
    axis = axes[0] * (1 if axes[0, 0] >= 0 else -1)
    signed = centered @ np.array([-axis[1], axis[0]])
    result = np.c_[area_q, pulse_q, centered @ axis, signed, np.abs(signed),
                   density64(area_q, pulse_q),
                   np.log(np.maximum(pulse, 1e-9)) - np.log(np.maximum(area, 1e-9))]
    if result.shape != (len(parent), 7) or not np.isfinite(result).all():
        raise ValueError("G1 frozen feature construction failed")
    return result


def build_edu_features(frame: pd.DataFrame, roles: dict[str, str]) -> np.ndarray:
    raw = frame[[roles["dna_area"], roles["dna_pulse"], roles["edu_area"]]].to_numpy(float)
    if not len(raw) or not np.isfinite(raw).all():
        raise ValueError("EdU full-frame roles are empty or non-finite")
    area, pulse, edu = (average_percentile_rank(raw[:, i], np.float32) for i in range(3))
    centered = np.c_[area - area.mean(), pulse - pulse.mean()]
    _, _, axes = np.linalg.svd(centered, full_matrices=False)
    axis = axes[0] * (1 if axes[0, 0] >= 0 else -1)
    result = np.c_[area, pulse, centered @ axis,
                   centered @ np.array([-axis[1], axis[0]]),
                   density64(area, pulse, np.float32), edu,
                   density64(area, edu, np.float32)].astype(np.float32)
    if result.shape != (len(frame), 7) or not np.isfinite(result).all():
        raise ValueError("EdU frozen feature construction failed")
    return result


def load_portable_package(repository_root: Path):
    root = repository_root / "inst/models/portable_hgb_export001"
    captured = {}
    for name, expected in PORTABLE_PACKAGE_SHA256.items():
        content = read_regular_bytes(root / name)
        if hashlib.sha256(content).hexdigest() != expected:
            raise ValueError(f"portable package digest mismatch: {name}")
        captured[name] = content
    manifest = json.loads(captured["MANIFEST.json"])
    ledger = {item["path"]: item["sha256"] for item in manifest.get("files", [])}
    if (manifest.get("status") != "complete" or
            manifest.get("portable_model_format") != "PORTABLE_HGB_V1" or
            manifest.get("protected_data_present") is not False or
            ledger != {key: value for key, value in PORTABLE_PACKAGE_SHA256.items()
                       if key != "MANIFEST.json"}):
        raise ValueError("portable package manifest contract mismatch")
    predictor = types.ModuleType("facspseudocolor_portable_hgb_predictor")
    exec(compile(captured["portable_hgb_predictor.py"],
                 "portable_hgb_predictor.py", "exec"), predictor.__dict__)
    models = {}
    for model_id, filename, features, threshold in (
        (G1_MODEL_ID, "DNA_PROJECTED_G1_CLUSTER001.portable.json", G1_FEATURES, G1_THRESHOLD),
        (EDU_MODEL_ID, "EDU_POSITIVE_MODEL002.portable.json", EDU_FEATURES, EDU_THRESHOLD)):
        model = json.loads(captured[filename])
        # Validate the exact schema formerly enforced by load_model, but from
        # the same digest-verified bytes that are retained for inference.
        if (model.get("format") != predictor.FORMAT or
                model.get("task") != "binary_classification" or
                model.get("link") != "logit" or
                model.get("model_id") != model_id or
                model.get("feature_order") != features or
                model.get("n_features_in") != len(features) or
                predictor._float(model["decision_threshold"]) != threshold):
            raise ValueError(f"portable frozen model contract mismatch: {model_id}")
        for iteration in model.get("iterations", []):
            if len(iteration) != 1 or not iteration[0].get("nodes"):
                raise ValueError(f"portable tree schema mismatch: {model_id}")
            for index, node in enumerate(iteration[0]["nodes"]):
                if node.get("index") != index:
                    raise ValueError(f"portable node order mismatch: {model_id}")
                predictor._float(node["value"])
                predictor._float(node["num_threshold"])
        models[model_id] = model
    return predictor, models, root


def verify_artifacts(config: dict):
    artifacts = config.get("artifacts")
    if not isinstance(artifacts, dict) or set(artifacts) != {
            "single_cells", "g1", "edu_positive"}:
        raise ValueError("exactly three provenance-distinct model artifacts are required")
    single = artifacts["single_cells"]
    if (single.get("model_id") != "SINGLE_CELLS_FROZEN_DEVELOPMENT_RULE" or
            single.get("byte_sha256") != EXPECTED_SINGLE_BYTE_SHA256 or
            single.get("semantic_sha256") != EXPECTED_SINGLE_SEMANTIC_SHA256 or
            single.get("feature_schema") != EXPECTED_SINGLE_SCHEMA or
            single.get("threshold") != EXPECTED_SINGLE_THRESHOLD):
        raise ValueError("Single Cells artifact contract is not the hard-pinned approved rule")
    single_rule = load_rule(Path(single["path"]), single["byte_sha256"],
                            single["semantic_sha256"], "single_cells_reference")
    if not np.isclose(single_rule["threshold"], EXPECTED_SINGLE_THRESHOLD,
                      rtol=0, atol=1e-15):
        raise ValueError("Single Cells threshold differs from approved frozen contract")
    g1 = artifacts["g1"]
    if (g1.get("model_id") != G1_MODEL_ID or
            g1.get("byte_sha256") != ORIGINAL_MODEL_SHA256[G1_MODEL_ID] or
            g1.get("portable_sha256") != PORTABLE_PACKAGE_SHA256[
                "DNA_PROJECTED_G1_CLUSTER001.portable.json"] or
            g1.get("threshold") != G1_THRESHOLD or
            g1.get("feature_schema") != G1_FEATURES or
            g1.get("metadata_sha256") != G1_METADATA):
        raise ValueError("DNA_PROJECTED_G1_CLUSTER001 artifact contract mismatch")
    g1_path = Path(g1["path"])
    if hashlib.sha256(read_regular_bytes(g1_path)).hexdigest() != g1["byte_sha256"]:
        raise ValueError("DNA_PROJECTED_G1_CLUSTER001 original model digest mismatch")
    for name, expected in G1_METADATA.items():
        if hashlib.sha256(read_regular_bytes(g1_path.parent / name)).hexdigest() != expected:
            raise ValueError(f"DNA_PROJECTED_G1_CLUSTER001 {name} digest mismatch")
    edu = artifacts["edu_positive"]
    if (edu.get("model_id") != EDU_MODEL_ID or
            edu.get("byte_sha256") != ORIGINAL_MODEL_SHA256[EDU_MODEL_ID] or
            edu.get("portable_sha256") != PORTABLE_PACKAGE_SHA256[
                "EDU_POSITIVE_MODEL002.portable.json"] or
            edu.get("threshold") != EDU_THRESHOLD or
            edu.get("feature_schema") != EDU_FEATURES or
            edu.get("config_sha256") != EDU_CONFIG_SHA256 or
            edu.get("metadata_sha256") != EDU_METADATA):
        raise ValueError("EDU_POSITIVE_MODEL002 artifact contract mismatch")
    edu_path = Path(edu["path"])
    if hashlib.sha256(read_regular_bytes(edu_path)).hexdigest() != edu["byte_sha256"]:
        raise ValueError("EDU_POSITIVE_MODEL002 original model digest mismatch")
    if hashlib.sha256(read_regular_bytes(Path(edu["config_path"]))).hexdigest() != EDU_CONFIG_SHA256:
        raise ValueError("EDU_POSITIVE_MODEL002 configuration digest mismatch")
    for name, expected in EDU_METADATA.items():
        if hashlib.sha256(read_regular_bytes(edu_path.parent / name)).hexdigest() != expected:
            raise ValueError(f"EDU_POSITIVE_MODEL002 {name} digest mismatch")
    predictor, models, portable_root = load_portable_package(
        Path(__file__).resolve().parents[1])
    manifest_artifacts = {"single_cells": dict(single), "g1": {
        **g1, "portable_path": str(portable_root /
            "DNA_PROJECTED_G1_CLUSTER001.portable.json")}, "edu_positive": {
        **edu, "portable_path": str(portable_root /
            "EDU_POSITIVE_MODEL002.portable.json")}}
    return manifest_artifacts, single_rule, predictor, models


def positive_probability(predictor, model: dict, features: np.ndarray,
                         feature_order: list[str]) -> np.ndarray:
    result = np.asarray(predictor.predict_positive_proba(
        model, features.tolist(), feature_order), dtype=float)
    if result.shape != (features.shape[0],):
        raise ValueError("portable frozen model returned an invalid probability vector")
    if not np.isfinite(result).all() or np.any((result < 0) | (result > 1)):
        raise ValueError("frozen estimator returned invalid probabilities")
    return result


def identity_sha256(values) -> str:
    identities = [str(value) for value in values]
    if len(identities) != len(set(identities)) or any(not value for value in identities):
        raise ValueError("event identities must be unique nonblank values")
    payload = ("\n".join(identities) + ("\n" if identities else "")).encode()
    return hashlib.sha256(payload).hexdigest()


def model_qc(name: str, probabilities: np.ndarray, calls: np.ndarray) -> list[str]:
    result = []
    if len(probabilities) < 1000:
        result.append(f"{name}:LOW_EVENT_SUPPORT")
    fraction = float(np.mean(calls)) if len(calls) else 0.0
    if fraction < 0.005 or fraction > 0.95:
        result.append(f"{name}:EXTREME_PREDICTED_FRACTION")
    if not np.any(calls):
        result.append(f"{name}:ZERO_PREDICTED_POSITIVES")
    return result


def publish_directory_noreplace(staging: Path, output: Path) -> None:
    """Atomically publish a directory and refuse an existing destination."""
    libc = ctypes.CDLL(None, use_errno=True)
    source = os.fsencode(staging)
    destination = os.fsencode(output)
    if sys.platform == "darwin" and hasattr(libc, "renamex_np"):
        libc.renamex_np.argtypes = [ctypes.c_char_p, ctypes.c_char_p,
                                    ctypes.c_uint]
        libc.renamex_np.restype = ctypes.c_int
        status = libc.renamex_np(source, destination, 0x00000004)  # RENAME_EXCL
    elif sys.platform.startswith("linux") and hasattr(libc, "renameat2"):
        libc.renameat2.argtypes = [ctypes.c_int, ctypes.c_char_p, ctypes.c_int,
                                   ctypes.c_char_p, ctypes.c_uint]
        libc.renameat2.restype = ctypes.c_int
        status = libc.renameat2(-100, source, -100, destination, 1)  # RENAME_NOREPLACE
    else:
        raise RuntimeError("platform lacks an atomic no-replace directory publish primitive")
    if status != 0:
        error = ctypes.get_errno()
        if error == errno.EEXIST:
            raise FileExistsError(output)
        raise OSError(error, os.strerror(error), output)


def score_acquisition(frame: pd.DataFrame, roles: dict[str, str],
                      single_rule, predictor, portable_models):
    single_roles = {
        "fsc_area": roles["fsc_area"], "ssc_area": roles["ssc_area"],
        "dna_area": roles["dna_area"], "dna_height": roles["dna_pulse"],
        "edu": roles["edu_area"],
    }
    single_features = extract_frozen_features(frame, single_roles)
    single_probability, single_mask = predict(single_features, single_rule)
    if not single_mask.any():
        raise ValueError("predicted Single Cells parent is empty")
    parent = frame.loc[single_mask].copy(deep=True)
    parent_ids = parent["event_identity"].to_numpy(copy=True)
    g1_features = build_g1_features(parent, roles)
    g1_probability = positive_probability(
        predictor, portable_models[G1_MODEL_ID], g1_features, G1_FEATURES)
    g1_mask = g1_probability >= G1_THRESHOLD
    edu_features_full = build_edu_features(frame, roles)
    edu_probability = positive_probability(
        predictor, portable_models[EDU_MODEL_ID],
        edu_features_full[single_mask], EDU_FEATURES)
    edu_mask = edu_probability >= EDU_THRESHOLD
    parent["model.single_cells.probability"] = single_probability[single_mask]
    parent["model.dna_projected_g1_cluster001.probability"] = g1_probability
    parent["model.edu_positive_model002.probability"] = edu_probability
    g1 = parent.loc[g1_mask].copy(deep=True)
    positive = parent.loc[edu_mask].copy(deep=True)
    rename = {roles["dna_area"]: "DNA content", roles["edu_area"]: "EdU"}
    tables = {"single_cells": parent.rename(columns=rename),
              "g1": g1.rename(columns=rename),
              "edu_positive": positive.rename(columns=rename)}
    audit = {
        "all_event_count": int(len(frame)),
        "model_single_cells_count": int(len(parent)),
        "model_g1_count": int(len(g1)),
        "model_edu_positive_count": int(len(positive)),
        "g1_subset_predicted_single_cells": bool(set(g1.event_identity) <= set(parent_ids)),
        "edu_positive_subset_predicted_single_cells": bool(set(positive.event_identity) <= set(parent_ids)),
        "parent_event_identity_exact_alignment": bool(
            np.array_equal(parent.event_identity.to_numpy(), parent_ids)),
        "all_event_identity_sha256": identity_sha256(frame.event_identity),
        "predicted_single_cells_identity_sha256": identity_sha256(parent.event_identity),
        "model_g1_identity_sha256": identity_sha256(g1.event_identity),
        "model_edu_positive_identity_sha256": identity_sha256(positive.event_identity),
        "g1_feature_scope": "predicted_single_cells_parent",
        "edu_positive_feature_scope": "complete_source_frame_then_predicted_single_cells_parent",
        "qc_warnings": (
            model_qc("single_cells", single_probability, single_mask) +
            model_qc("g1", g1_probability, g1_mask) +
            model_qc("edu_positive", edu_probability, edu_mask)
        ),
    }
    if not all(audit[key] for key in audit if key.endswith("single_cells") or
               key.endswith("exact_alignment")):
        raise ValueError("model output containment or event alignment failed")
    return tables, audit


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path)
    args = parser.parse_args()
    config_bytes = read_regular_bytes(args.config)
    config = json.loads(config_bytes)
    if config.get("status_label") != STATUS or config.get("predicted_parent_authorization") != APPROVAL:
        parser.error("explicit model-only non-production owner authorization is required")
    if config.get("flowjo_workspace_gate_used") is not False:
        parser.error("flowjo_workspace_gate_used must be false")
    workspace_scan = dict(config)
    workspace_scan.pop("flowjo_workspace_gate_used", None)
    if "workspace" in json.dumps(workspace_scan, sort_keys=True).casefold():
        parser.error("model-only run configuration must contain no workspace fields")
    repository_root = Path(__file__).resolve().parents[1]
    approved = load_approved_mapping(repository_root)
    if config.get("approved_mapping_sha256") != MAPPING_SHA256:
        parser.error("run does not select the immutable approved Figure6 mapping")
    output = validate_output_dir(Path(config["output_dir"]),
                                 repository_root)
    acquisitions = config.get("acquisitions")
    if not isinstance(acquisitions, list) or len(acquisitions) != 8:
        parser.error("the exact eight approved Figure6 acquisitions are required")
    prefixes = [item.get("prefix") for item in acquisitions]
    approved_prefixes = [item["prefix"] for item in approved["acquisitions"]]
    if prefixes != approved_prefixes or len(prefixes) != len(set(prefixes)):
        parser.error("acquisition order/identity differs from immutable Figure6 mapping")
    artifacts, single_rule, portable_predictor, portable_models = verify_artifacts(config)
    output.parent.mkdir(parents=True, exist_ok=True)
    claim = output.parent / f".{output.name}.claim"
    descriptor = None
    claim_created = False
    staging = None
    manifest = {
        "schema_version": 3, "status_label": STATUS,
        "predicted_parent_authorization": APPROVAL,
        "single_cells_source": "model", "g1_source": "model",
        "edu_positive_source": "model", "flowjo_workspace_gate_used": False,
        "approved_mapping": {
            "path": "inst/config/model_only_edu_figure6_mapping.json",
            "sha256": MAPPING_SHA256,
            "status": approved["status"],
        },
        "run_configuration_sha256": hashlib.sha256(config_bytes).hexdigest(),
        "artifacts": artifacts,
        "feature_scopes": {
            "single_cells": "complete_source_frame",
            "g1": "predicted_single_cells_parent",
            "edu_positive": "complete_source_frame_then_predicted_single_cells_parent",
        },
        "acquisitions": [],
    }
    try:
        descriptor = os.open(claim, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
        claim_created = True
        staging = Path(tempfile.mkdtemp(prefix=f".{output.name}.", suffix=".partial",
                                        dir=output.parent))
        for item in acquisitions:
            roles, mapping_record = validate_mapping(item, approved)
            frame, consumed_digest = read_raw_fcs(item, roles)
            with warnings.catch_warnings(record=True) as scoring_warnings:
                warnings.simplefilter("always")
                tables, audit = score_acquisition(
                    frame, roles, single_rule, portable_predictor, portable_models)
            if scoring_warnings:
                raise RuntimeError("unexpected scientific warning during model scoring: " +
                                   " | ".join(str(x.message) for x in scoring_warnings))
            records = {}
            for name, suffix in (("single_cells", "_single_cells.csv"),
                                 ("g1", "_g1.csv"),
                                 ("edu_positive", "_edu_positive.csv")):
                path = staging / f"{item['prefix']}{suffix}"
                write_csv_absent(tables[name], path)
                records[name] = {
                    "path": path.name, "sha256": sha256_path(path),
                    "rows": int(len(tables[name])),
                    "event_identity_sha256": identity_sha256(
                        tables[name].event_identity),
                }
            manifest["acquisitions"].append({
                "prefix": item["prefix"], "acquisition_id": item["acquisition_id"],
                "source_fcs": str(Path(item["fcs_path"]).resolve()),
                "source_fcs_basename": mapping_record["fcs_basename"],
                "source_fcs_sha256": item["fcs_sha256"],
                "consumed_fcs_snapshot_sha256": consumed_digest,
                "functional_panel": approved["functional_panel"],
                "channel_roles": roles, "audit": audit, "outputs": records,
            })
        write_json_absent(manifest, staging / "model-gating-manifest.json")
        publish_directory_noreplace(staging, output)
        staging = None
        return 0
    except BaseException:
        if staging is not None:
            shutil.rmtree(staging, ignore_errors=True)
        raise
    finally:
        if descriptor is not None:
            os.close(descriptor)
        if claim_created:
            claim.unlink(missing_ok=True)


if __name__ == "__main__":
    raise SystemExit(main())
