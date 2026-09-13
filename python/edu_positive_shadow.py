#!/usr/bin/env python3
"""Digest-bound, shadow-only EDU_POSITIVE_MODEL002 inference.

Features are constructed over a complete per-source event frame and predictions
are emitted only for rows selected by a separate expert-reference Single Cells
mask. This module never reads FCS/WSP files or expert EdU annotations and never
constructs or predicts a parent population.
"""
from __future__ import annotations

import argparse
import io
import csv
import hashlib
import json
import os
import pickle
import platform
import re
import sys
import shutil
import stat
import tempfile
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


MODEL_SHA256 = "6ef5603555a660b8503379cbbcab131b616336a64330ffd42f45dfaa182042cc"
CONFIG_SHA256 = "624a839a38a464976b214afcf8964a5f54009ef7694d16d71a91ba20677d44f3"
FEATURES = [
    "dna_area_q", "dna_pulse_q", "dna_tls_position", "dna_tls_distance",
    "dna_area_pulse_log_density", "edu_area_q", "dna_edu_log_density",
]
THRESHOLD = 0.385
MODEL_NAME = "hgb_small"
SCOPE = "expert_reference_single_cells"
TRANSFORMATION_SOURCE_SHA256 = "4233d4ed8adc377d2c56791a323dd9354894ef776b35924484ddccae82d2bbd2"
APPROVED_PANELS = ["PANEL_FL2_DNA_FL4_EDU", "PANEL_PI_DNA_ALEXA647_EDU"]
FUNCTIONAL_MANIFEST_SHA256 = "c42de633eb9ca6d9a7a22373580a9d1f9ea50b81e15b39a3a603f1d9b5089055"


class ShadowStop(RuntimeError):
    def __init__(self, code: str, message: str, **details: Any):
        super().__init__(message)
        self.record = {"status": "hard_stop", "code": code, "message": message,
                       "details": details}


def fail(message: str, code: str = "validation_failed", **details: Any) -> None:
    raise ShadowStop(code, message, **details)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_regular_bytes(path: Path, role: str) -> tuple[bytes, str]:
    if path.is_symlink() or not path.is_file():
        fail(f"{role} must be a regular non-symlink file", "unsafe_input_path")
    flags = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0)
    descriptor = os.open(path, flags)
    if not stat.S_ISREG(os.fstat(descriptor).st_mode):
        os.close(descriptor)
        fail(f"{role} descriptor is not a regular file", "unsafe_input_descriptor")
    with os.fdopen(descriptor, "rb") as handle:
        content = handle.read()
    return content, hashlib.sha256(content).hexdigest()


def _load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        fail(f"expected one JSON object: {path}")
    return value


def validate_config_bytes(content: bytes, digest: str) -> dict[str, Any]:
    if digest != CONFIG_SHA256:
        fail("frozen shadow configuration digest mismatch")
    config = json.loads(content.decode("utf-8"))
    if not isinstance(config, dict):
        fail("configuration must be one JSON object")
    required = {
        "schema_version": 1, "model_id": "EDU_POSITIVE_MODEL002",
        "model_sha256": MODEL_SHA256, "selected_model": MODEL_NAME,
        "threshold": THRESHOLD, "ordered_features": FEATURES,
        "input_scope": SCOPE, "mode": "shadow_only",
        "source_transformation_sha256": TRANSFORMATION_SOURCE_SHA256,
        "approved_functional_panels": APPROVED_PANELS,
        "functional_panel_source_manifest_sha256": FUNCTIONAL_MANIFEST_SHA256,
        "feature_scope": "full_source_frame_then_expert_single_cells_subset",
        "recorded_runtime": {
            "python": "3.10.12", "numpy": "2.2.6", "pandas": "2.3.3",
            "scipy": "1.15.3", "scikit-learn": "1.7.2",
        },
    }
    for key, expected in required.items():
        if config.get(key) != expected:
            fail(f"frozen shadow configuration field mismatch: {key}")
    formula = config.get("feature_construction", {})
    expected_formula = {
        "quantile_rank": "pandas_rank_method_average_pct_true_float32",
        "tls": "centered_D_P_svd_vt0_sign_nonnegative_D_loading",
        "density": "numpy_histogram2d_64_bins_unit_square_log1p_bin_count_over_log1p_n",
    }
    if formula != expected_formula:
        fail("frozen feature-construction declaration mismatch")
    return config


def load_model_checked(path: Path) -> Any:
    if path.is_symlink() or not path.is_file():
        fail("model must be a regular non-symlink file", "unsafe_model_path")
    flags = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0)
    descriptor = os.open(path, flags)
    if not stat.S_ISREG(os.fstat(descriptor).st_mode):
        os.close(descriptor)
        fail("model descriptor is not a regular file", "unsafe_model_descriptor")
    with os.fdopen(descriptor, "rb") as handle:
        model_bytes = handle.read()
    if hashlib.sha256(model_bytes).hexdigest() != MODEL_SHA256:
        fail("EDU_POSITIVE_MODEL002 digest mismatch")
    # The digest allowlist is checked before deserialization. Pickle remains
    # restricted to this exact owner-approved artifact.
    payload = pickle.loads(model_bytes)
    if not isinstance(payload, dict):
        fail("model payload is not the frozen study dictionary")
    if payload.get("features") != FEATURES:
        fail("model feature order mismatch")
    if payload.get("winner") != MODEL_NAME:
        fail("selected model mismatch")
    if not np.isclose(float(payload.get("threshold", np.nan)), THRESHOLD,
                      rtol=0.0, atol=1e-12):
        fail("model threshold mismatch")
    model = payload.get("model")
    if model is None or not callable(getattr(model, "predict_proba", None)):
        fail("model lacks predict_proba")
    return model


def validate_runtime(config: dict[str, Any]) -> dict[str, str]:
    import scipy
    import sklearn
    observed = {
        "python": platform.python_version(), "numpy": np.__version__,
        "pandas": pd.__version__, "scipy": scipy.__version__,
        "scikit-learn": sklearn.__version__,
    }
    expected = config["recorded_runtime"]
    if observed != expected:
        fail("runtime version mismatch", "runtime_mismatch",
             expected=expected, observed=observed)
    return observed


def _rank(values: np.ndarray) -> np.ndarray:
    return pd.Series(values).rank(method="average", pct=True).to_numpy(np.float32)


def _density(x: np.ndarray, y: np.ndarray, bins: int = 64) -> np.ndarray:
    histogram, _, _ = np.histogram2d(
        x, y, bins=bins, range=((0, 1), (0, 1))
    )
    i = np.minimum((x * bins).astype(int), bins - 1)
    j = np.minimum((y * bins).astype(int), bins - 1)
    return (np.log1p(histogram[i, j]) / np.log1p(len(x))).astype(np.float32)


def construct_features(frame: pd.DataFrame, roles: dict[str, str]) -> np.ndarray:
    if not len(frame):
        fail("refusing empty full source event frame")
    required_roles = {"dna_area", "dna_pulse", "edu_area"}
    if set(roles) != required_roles:
        fail("channel roles must be exactly dna_area, dna_pulse, and edu_area")
    columns = [roles[name] for name in ("dna_area", "dna_pulse", "edu_area")]
    if len(set(columns)) != 3 or any(not isinstance(x, str) or not x for x in columns):
        fail("canonical roles must map to three distinct raw columns")
    missing = [name for name in columns if name not in frame.columns]
    if missing:
        fail("mapped raw channel column missing: " + ", ".join(missing))
    raw = frame[columns].to_numpy(float)
    if not np.isfinite(raw).all():
        fail("non-finite canonical role value")
    dna_area, dna_pulse, edu_area = (_rank(raw[:, i]) for i in range(3))
    centered = np.c_[dna_area - dna_area.mean(), dna_pulse - dna_pulse.mean()]
    _, _, vt = np.linalg.svd(centered, full_matrices=False)
    axis = vt[0] * (1 if vt[0, 0] >= 0 else -1)
    result = np.c_[
        dna_area, dna_pulse, centered @ axis,
        centered @ np.array([-axis[1], axis[0]]),
        _density(dna_area, dna_pulse), edu_area,
        _density(dna_area, edu_area),
    ].astype(np.float32)
    if result.shape != (len(frame), len(FEATURES)) or not np.isfinite(result).all():
        fail("feature construction failure")
    return result


def validate_source_mapping(mapping_bytes: bytes, source_id: str,
                            input_digest: str, mask_digest: str,
                            functional_bytes: bytes, functional_digest: str
                            ) -> tuple[dict[str, str], dict[str, Any]]:
    manifest = json.loads(mapping_bytes.decode("utf-8"))
    if not isinstance(manifest, dict):
        fail("mapping manifest must be one JSON object")
    if manifest.get("schema_version") != 1 or manifest.get("input_scope") != SCOPE:
        fail("mapping manifest is not approved for expert-reference Single Cells")
    sources = manifest.get("sources")
    if not isinstance(sources, list) or any(not isinstance(item, dict) for item in sources):
        fail("mapping manifest sources must be a list")
    matches = [item for item in sources if item.get("source_id") == source_id]
    if len(matches) != 1:
        fail("source mapping missing or ambiguous")
    item = matches[0]
    if item.get("functional_panel") not in APPROVED_PANELS:
        fail("source is not assigned to an approved functional panel")
    if item.get("full_event_frame_sha256") != input_digest:
        fail("full event-frame digest mismatch")
    if item.get("single_cells_mask_sha256") != mask_digest:
        fail("Single Cells mask digest mismatch")
    if functional_digest != FUNCTIONAL_MANIFEST_SHA256:
        fail("functional-panel source manifest digest mismatch")
    functional_rows = list(csv.DictReader(io.StringIO(functional_bytes.decode("utf-8"))))
    functional_matches = [row for row in functional_rows
                          if row.get("edu_label_source_id") == source_id]
    if len(functional_matches) != 1:
        fail("source absent or ambiguous in frozen functional-panel manifest")
    frozen_panel = functional_matches[0].get("functional_panel_status")
    if frozen_panel != item.get("functional_panel") or frozen_panel not in APPROVED_PANELS:
        fail("source panel does not match frozen functional-panel manifest")
    roles = item.get("raw_to_canonical_roles")
    if not isinstance(roles, dict):
        fail("raw-to-canonical channel roles missing")
    return roles, item


def read_single_cells_mask(mask_bytes: bytes, frame: pd.DataFrame,
                           identity_column: str) -> np.ndarray:
    mask_frame = pd.read_csv(io.BytesIO(mask_bytes))
    if list(mask_frame.columns) != [identity_column, "expert_single_cells"]:
        fail("mask columns must be event identity then expert_single_cells")
    if len(mask_frame) != len(frame):
        fail("mask/frame row-count mismatch")
    left = frame[identity_column].astype(str).to_numpy()
    right = mask_frame[identity_column].astype(str).to_numpy()
    if not np.array_equal(left, right):
        fail("mask/frame event identities are not exactly aligned")
    values = mask_frame["expert_single_cells"]
    if not values.isin([0, 1, False, True]).all():
        fail("expert Single Cells mask must contain only 0/1")
    mask = values.astype(bool).to_numpy()
    if not mask.any():
        fail("expert Single Cells mask selects no events")
    return mask


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        fail("refusing empty shadow predictions")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def run_shadow(*, input_path: Path, mask_path: Path, source_id: str,
               event_identity_column: str, mapping_path: Path,
               functional_manifest_path: Path, model_path: Path,
               config_path: Path, project_root: Path, shadow_base: Path,
               run_id: str) -> Path:
    for path in (input_path, mask_path, mapping_path, functional_manifest_path,
                 model_path, config_path):
        if not path.is_file():
            fail(f"required input missing: {path}")
    project_root = project_root.resolve()
    shadow_base = shadow_base.resolve()
    if project_root not in shadow_base.parents:
        fail("shadow base must be inside the declared project root")
    if re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", run_id) is None:
        fail("run ID must be one safe path component")
    output_root = shadow_base / "EDU_POSITIVE_MODEL002" / run_id
    if output_root.exists():
        fail("shadow output root must be absent")
    output_root.parent.mkdir(parents=True, exist_ok=True)
    claim = output_root.parent / f".{output_root.name}.claim"
    descriptor = None
    staging = Path(tempfile.mkdtemp(prefix=f".{output_root.name}.",
                                    dir=str(output_root.parent)))
    try:
        descriptor = os.open(claim, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
        config_bytes, config_digest = read_regular_bytes(config_path, "configuration")
        mapping_bytes, mapping_digest = read_regular_bytes(mapping_path, "mapping manifest")
        functional_bytes, functional_digest = read_regular_bytes(
            functional_manifest_path, "functional-panel source manifest")
        frame_bytes, frame_digest = read_regular_bytes(input_path, "full event frame")
        mask_bytes, mask_digest = read_regular_bytes(mask_path, "Single Cells mask")
        config = validate_config_bytes(config_bytes, config_digest)
        observed_runtime = validate_runtime(config)
        roles, source_record = validate_source_mapping(
            mapping_bytes, source_id, frame_digest, mask_digest,
            functional_bytes, functional_digest)
        frame = pd.read_csv(io.BytesIO(frame_bytes))
        if event_identity_column not in frame.columns:
            fail("required direct event-identity column missing")
        identities = frame[event_identity_column]
        if identities.isna().any() or (identities.astype(str).str.len() == 0).any():
            fail("event identity contains missing values")
        if identities.astype(str).duplicated().any():
            fail("event identity is not unique within source")
        mask = read_single_cells_mask(mask_bytes, frame, event_identity_column)
        model = load_model_checked(model_path)
        full_features = construct_features(frame, roles)
        probability = np.asarray(model.predict_proba(full_features[mask]))[:, 1]
        selected_identities = identities[mask]
        if probability.shape != (int(mask.sum()),) or not np.isfinite(probability).all():
            fail("invalid model probability output")
        rows = [{
            "source_id": source_id,
            "event_identity": str(identity),
            "edu_positive_shadow_probability": format(float(value), ".17g"),
            "edu_positive_shadow_prediction": int(value >= THRESHOLD),
            "annotation_status": "RESEARCH_SHADOW_NOT_BIOLOGICAL_GROUND_TRUTH",
        } for identity, value in zip(selected_identities, probability)]
        prediction_path = staging / "edu_positive_shadow_predictions.csv"
        _write_csv(prediction_path, rows)
        provenance = {
            "schema_version": 1, "status": "complete", "mode": "shadow_only",
            "model_id": config["model_id"], "selected_model": MODEL_NAME,
            "model_sha256": MODEL_SHA256, "configuration_sha256": CONFIG_SHA256,
            "mapping_manifest_sha256": mapping_digest,
            "functional_panel_manifest_sha256": FUNCTIONAL_MANIFEST_SHA256,
            "full_event_frame_sha256": frame_digest,
            "single_cells_mask_sha256": mask_digest,
            "input_scope": SCOPE, "source_id": source_id,
            "functional_panel": source_record["functional_panel"],
            "event_identity_column": event_identity_column,
            "full_source_events": len(frame),
            "expert_single_cells_events": int(mask.sum()),
            "predicted_events": len(rows),
            "coverage_complete": len(rows) == int(mask.sum()),
            "threshold": THRESHOLD, "ordered_features": FEATURES,
            "accepted_runtime": observed_runtime,
            "output_interpretation": "research shadow reproduction of expert annotations; not biological ground truth or production gating",
            "expert_annotations_read": False, "expert_annotations_modified": False,
            "parent_predicted": False, "evaluation_metrics_computed": False,
        }
        provenance_path = staging / "provenance.json"
        provenance_path.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n",
                                   encoding="utf-8")
        ledger = [
            {"file": path.name, "bytes": path.stat().st_size, "sha256": sha256_file(path)}
            for path in (prediction_path, provenance_path)
        ]
        _write_csv(staging / "ARTIFACT_SHA256.csv", ledger)
        if output_root.exists():
            fail("shadow output collision")
        os.rename(staging, output_root)
        return output_root
    except BaseException:
        if staging.exists():
            shutil.rmtree(staging)
        raise
    finally:
        if descriptor is not None:
            os.close(descriptor)
            claim.unlink(missing_ok=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True,
                        help="CSV containing the full per-source event frame")
    parser.add_argument("--single-cells-mask", type=Path, required=True)
    parser.add_argument("--source-id", required=True)
    parser.add_argument("--event-identity-column", required=True)
    parser.add_argument("--mapping", type=Path, required=True)
    parser.add_argument("--functional-panel-manifest", type=Path, required=True)
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--shadow-base", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    args = parser.parse_args()
    try:
        run_shadow(input_path=args.input, mask_path=args.single_cells_mask,
                   source_id=args.source_id,
                   event_identity_column=args.event_identity_column,
                   mapping_path=args.mapping,
                   functional_manifest_path=args.functional_panel_manifest,
                   model_path=args.model, config_path=args.config,
                   project_root=args.project_root, shadow_base=args.shadow_base,
                   run_id=args.run_id)
    except ShadowStop as error:
        print(json.dumps(error.record, sort_keys=True), file=sys.stderr)
        raise SystemExit(2) from error


if __name__ == "__main__":
    main()
