"""Fail-closed DNA-projected G1 experimental shadow inference.

This module scores only a caller-supplied, explicitly established Single Cells
event table.  It never discovers gates, channels, or source files.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
import pickle
import platform
import tempfile
import shutil
import re
import stat
from pathlib import Path
from typing import Any

import numpy as np


MODEL_ID = "DNA_PROJECTED_G1_CLUSTER001"
MODEL_SHA256 = "9ce7e448455e2d42090a86c80beb10d925a7c0e64e6261e678335385d9714af6"
METADATA_SHA256 = {
    "SUMMARY.json": "7a7446d61ca175e56964fc5ab6edbee081f2eb3f16e7925542157bc3719ddb10",
    "PREFLIGHT.json": "bac957a8837774afee0faba2ff2e832d003f1e37586a68e6af4037083448ec80",
}
THRESHOLD = 0.42
FEATURES = (
    "dna_area_q", "dna_pulse_q", "tls_position", "tls_signed_distance",
    "tls_abs_distance", "dna_2d_log_density", "pulse_area_log_ratio",
)
NAMESPACE = "shadow.dna_projected_g1_cluster001"
DISPLAY_NAME = "DNA-projected G1 proxy (experimental shadow)"


class ShadowContractError(ValueError):
    """Raised before prediction when the frozen shadow contract is not met."""


def _read_regular_file_once(path: Path) -> bytes:
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0)
    if not hasattr(os, "O_NOFOLLOW"):
        raise ShadowContractError("platform cannot enforce no-follow artifact reads")
    flags |= os.O_NOFOLLOW
    try:
        descriptor = os.open(path, flags)
    except OSError as error:
        raise ShadowContractError(
            f"artifact must be a regular non-symlink file: {path.name}"
        ) from error
    try:
        if not stat.S_ISREG(os.fstat(descriptor).st_mode):
            raise ShadowContractError(
                f"artifact must be a regular non-symlink file: {path.name}"
            )
        with os.fdopen(descriptor, "rb") as handle:
            descriptor = -1
            return handle.read()
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def _average_percentile_rank(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values, kind="mergesort")
    ranked = np.empty(values.size, dtype=float)
    start = 0
    while start < values.size:
        end = start + 1
        while end < values.size and values[order[end]] == values[order[start]]:
            end += 1
        ranked[order[start:end]] = ((start + 1) + end) / 2.0 / values.size
        start = end
    return ranked


def build_features(events: dict[str, Any], channel_roles: dict[str, str]) -> np.ndarray:
    """Build the frozen seven columns from selected Single Cells only."""
    if set(channel_roles) != {"DNA_A", "DNA_PULSE"}:
        raise ShadowContractError("exactly one DNA_A and DNA_PULSE role is required")
    area_name, pulse_name = channel_roles["DNA_A"], channel_roles["DNA_PULSE"]
    if not all(isinstance(x, str) and x for x in (area_name, pulse_name)):
        raise ShadowContractError("channel role mappings must be explicit raw parameter names")
    if area_name == pulse_name:
        raise ShadowContractError("DNA_A and DNA_PULSE must resolve to distinct parameters")
    if area_name not in events or pulse_name not in events:
        raise ShadowContractError("resolved DNA channel is absent from selected events")
    area = np.asarray(events[area_name], dtype=float)
    pulse = np.asarray(events[pulse_name], dtype=float)
    if area.ndim != 1 or pulse.ndim != 1 or area.size != pulse.size or area.size == 0:
        raise ShadowContractError("selected Single Cells channels must be nonempty aligned vectors")
    if not np.isfinite(area).all() or not np.isfinite(pulse).all():
        raise ShadowContractError("selected raw DNA values contain non-finite values")

    area_q = _average_percentile_rank(area)
    pulse_q = _average_percentile_rank(pulse)
    quantiles = np.column_stack((area_q, pulse_q))
    centered = quantiles - quantiles.mean(axis=0)
    _, _, axes = np.linalg.svd(centered, full_matrices=False)
    axis = axes[0].copy()
    if axis[0] < 0:
        axis *= -1
    normal = np.array([-axis[1], axis[0]])
    position = centered @ axis
    signed_distance = centered @ normal
    area_bin = np.minimum((area_q * 64).astype(int), 63)
    pulse_bin = np.minimum((pulse_q * 64).astype(int), 63)
    counts = np.zeros((64, 64), dtype=np.int64)
    np.add.at(counts, (area_bin, pulse_bin), 1)
    density = np.log1p(counts[area_bin, pulse_bin]) / np.log1p(area.size)
    ratio = np.log(np.maximum(pulse, 1e-9)) - np.log(np.maximum(area, 1e-9))
    result = np.column_stack((
        area_q, pulse_q, position, signed_distance, np.abs(signed_distance),
        density, ratio,
    ))
    if result.shape[1] != len(FEATURES) or not np.isfinite(result).all():
        raise ShadowContractError("constructed frozen features are invalid")
    return result


def _validate_parent(parent: dict[str, Any], row_count: int) -> None:
    required = {"population_name", "source", "method", "selection_status"}
    if not isinstance(parent, dict) or not required.issubset(parent):
        raise ShadowContractError("explicit Single Cells parent provenance is required")
    if parent["population_name"] != "Single Cells" or parent["selection_status"] != "established":
        raise ShadowContractError("Single Cells parent is absent or ambiguous")
    if not isinstance(parent["source"], str) or not parent["source"] or not isinstance(parent["method"], str) or not parent["method"]:
        raise ShadowContractError("Single Cells parent source and method are required")
    if row_count == 0:
        raise ShadowContractError("Single Cells parent is empty")


def _load_verified_model(model_path: Path):
    model_bytes = _read_regular_file_once(model_path)
    if hashlib.sha256(model_bytes).hexdigest() != MODEL_SHA256:
        raise ShadowContractError("model artifact is absent or has a SHA-256 mismatch")
    for name, expected in METADATA_SHA256.items():
        path = model_path.parent / name
        content = _read_regular_file_once(path)
        if hashlib.sha256(content).hexdigest() != expected:
            raise ShadowContractError(f"{name} is absent or has a SHA-256 mismatch")
        if name == "SUMMARY.json":
            summary_bytes = content
    summary = json.loads(summary_bytes.decode("utf-8"))
    winner = summary.get("winner", {})
    if (tuple(summary.get("features", ())) != FEATURES or
            winner.get("model") != "hgb_regularized" or
            not np.isclose(winner.get("threshold", np.nan), THRESHOLD, rtol=0, atol=1e-15)):
        raise ShadowContractError("verified model metadata differs from the frozen contract")
    _validate_runtime()
    bundle = pickle.loads(model_bytes)  # trusted bytes verified before execution
    return _validated_estimator(bundle)


def _validated_estimator(bundle: Any):
    if not isinstance(bundle, dict):
        raise ShadowContractError("serialized model bundle schema is invalid")
    if (tuple(bundle.get("features", ())) != FEATURES or
            not np.isclose(bundle.get("threshold", np.nan), THRESHOLD,
                           rtol=0, atol=1e-15) or
            bundle.get("winner") != "hgb_regularized" or
            "model" not in bundle):
        raise ShadowContractError("serialized model bundle differs from the frozen contract")
    return bundle["model"]


def _validate_runtime() -> None:
    import sklearn
    import scipy
    expected = {"python": "3.10.12", "numpy": "2.2.6", "scipy": "1.15.3",
                "scikit_learn": "1.7.2"}
    observed = {"python": platform.python_version(), "numpy": np.__version__,
                "scipy": scipy.__version__,
                "scikit_learn": sklearn.__version__}
    if observed != expected:
        raise ShadowContractError(
            "runtime compatibility is not established: expected frozen "
            f"{expected}, observed {observed}"
        )


def predict_shadow(*, enabled: bool, events: dict[str, Any], event_identity: Any,
                   channel_roles: dict[str, str], parent: dict[str, Any],
                   model_path: Path) -> dict[str, Any]:
    """Return a separate research-only shadow result; disabled means no scoring."""
    if not enabled:
        return {"enabled": False, "namespace": NAMESPACE, "status": "disabled"}
    identities = np.asarray(event_identity, dtype=object)
    if identities.ndim != 1 or identities.size == 0 or any(not isinstance(x, str) or not x for x in identities):
        raise ShadowContractError("stable nonblank event_identity values are required")
    if len(set(identities.tolist())) != identities.size:
        raise ShadowContractError("event_identity values must be unique")
    _validate_parent(parent, identities.size)
    matrix = build_features(events, channel_roles)
    if matrix.shape[0] != identities.size:
        raise ShadowContractError("event identities and selected channels are not aligned")
    model = _load_verified_model(Path(model_path))
    if hasattr(model, "feature_names_in_") and tuple(model.feature_names_in_) != FEATURES:
        raise ShadowContractError("serialized model feature order differs from frozen schema")
    probabilities = np.asarray(model.predict_proba(matrix), dtype=float)
    if probabilities.ndim != 2 or probabilities.shape[0] != identities.size:
        raise ShadowContractError("model returned an invalid probability matrix")
    classes = list(model.classes_)
    if 1 not in classes:
        raise ShadowContractError("model has no positive class 1")
    positive = probabilities[:, classes.index(1)]
    if not np.isfinite(positive).all() or np.any((positive < 0) | (positive > 1)):
        raise ShadowContractError("model returned invalid probabilities")
    calls = positive >= THRESHOLD
    fraction = float(calls.mean())
    warnings = []
    if identities.size < 1000:
        warnings.append("LOW_EVENT_SUPPORT")
    if fraction < 0.005 or fraction > 0.95:
        warnings.append("EXTREME_PREDICTED_FRACTION")
    if not calls.any():
        warnings.append("ZERO_PROXY_POSITIVES")
    return {
        "enabled": True, "status": "complete", "namespace": NAMESPACE,
        "display_name": DISPLAY_NAME, "interpretation": "research-only annotation-reproduction prediction; not an expert gate or biological ground truth",
        "event_identity": identities.tolist(), "probability": positive.tolist(),
        "shadow_positive": calls.tolist(), "warnings": warnings,
        "provenance": {
            "model_id": MODEL_ID, "model_sha256": MODEL_SHA256,
            "threshold": THRESHOLD, "feature_order": list(FEATURES),
            "parent_population": dict(parent), "resolved_channel_roles": dict(channel_roles),
            "selected_event_count": int(identities.size),
            "predicted_positive_count": int(calls.sum()),
            "predicted_positive_fraction": fraction,
            "runtime": runtime_versions(),
        },
    }


def write_shadow_sidecar(result: dict[str, Any], run_root: Path,
                         project_root: Path) -> tuple[Path, Path]:
    """Atomically publish deterministic CSV/JSON beneath a new absent run root."""
    if result.get("status") != "complete" or result.get("namespace") != NAMESPACE:
        raise ShadowContractError("only a complete G1 shadow result can be exported")
    root = Path(run_root)
    project = Path(project_root)
    if project.is_symlink() or not project.is_dir():
        raise ShadowContractError("project_root must be an existing non-symlink directory")
    required_parent = project.resolve() / NAMESPACE
    if (root.parent.resolve() != required_parent or
            re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", root.name) is None):
        raise ShadowContractError(f"run root must be a neutral ID directly under {NAMESPACE}")
    if root.exists() or root.is_symlink():
        raise FileExistsError(root)
    identities = result.get("event_identity")
    probabilities = result.get("probability")
    calls = result.get("shadow_positive")
    if not all(isinstance(x, list) for x in (identities, probabilities, calls)):
        raise ShadowContractError("shadow output arrays must be lists")
    if not len(identities) or len(identities) != len(probabilities) or len(identities) != len(calls):
        raise ShadowContractError("shadow output arrays must have equal nonzero length")
    if any(not isinstance(x, str) or not x for x in identities) or len(set(identities)) != len(identities):
        raise ShadowContractError("shadow event identities must be nonblank unique strings")
    numeric = np.asarray(probabilities, dtype=float)
    if not np.isfinite(numeric).all() or np.any((numeric < 0) | (numeric > 1)):
        raise ShadowContractError("shadow probabilities must be finite values in [0, 1]")
    if any(not isinstance(x, bool) for x in calls):
        raise ShadowContractError("shadow calls must be boolean")
    required_parent.mkdir(exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{root.name}.", suffix=".partial", dir=root.parent))
    csv_path, manifest_path = staging / "predictions.csv", staging / "manifest.json"
    try:
        with csv_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, lineterminator="\n")
            writer.writerow(("event_identity", "dna_projected_g1_probability", "dna_projected_g1_shadow_positive"))
            for identity, probability, call in zip(identities, probabilities, calls):
                writer.writerow((identity, format(probability, ".17g"), "true" if call else "false"))
        prediction_bytes = csv_path.read_bytes()
        manifest = {key: result[key] for key in ("namespace", "display_name", "interpretation", "warnings", "provenance")}
        manifest["artifacts"] = [{
            "role": "event_aligned_shadow_predictions",
            "path": "predictions.csv",
            "sha256": hashlib.sha256(prediction_bytes).hexdigest(),
            "byte_size": len(prediction_bytes),
            "row_count": len(identities),
            "columns": ["event_identity", "dna_projected_g1_probability",
                        "dna_projected_g1_shadow_positive"],
        }]
        content = (json.dumps(manifest, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False) + "\n").encode("utf-8")
        manifest_path.write_bytes(content)
        os.rename(staging, root)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    return root / "predictions.csv", root / "manifest.json"


def runtime_versions() -> dict[str, str]:
    import sklearn
    import scipy
    return {"python": platform.python_version(), "numpy": np.__version__,
            "scipy": scipy.__version__, "scikit_learn": sklearn.__version__}
