"""SYNTHETIC contract tests for the branch-only model orchestration."""

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))
spec = importlib.util.spec_from_file_location(
    "apply_model_only_edu_models", ROOT / "python" / "apply_model_only_edu_models.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class DeterministicPredictor:
    @staticmethod
    def predict_positive_proba(model, rows, feature_order):
        assert model["feature_order"] == feature_order
        return np.linspace(0.1, 0.9, len(rows), dtype=float).tolist()


def synthetic_frame(n=1200):
    index = np.arange(n)
    return pd.DataFrame({
        "acquisition_id": "SYNTHETIC_ACQUISITION",
        "event_identity": [f"SYNTHETIC_ACQUISITION:event_index:{i}" for i in index],
        "event_index": index,
        "FSC-A": 100 + index % 61,
        "SSC-A": 200 + index % 47,
        "FL2-A": 500 + index * 0.5,
        "FL2-H": 250 + index * 0.25 + (index % 7),
        "FL4-A": 50 + index * 0.1 + (index % 13),
    })


def synthetic_roles():
    return {"fsc_area": "FSC-A", "ssc_area": "SSC-A",
            "dna_area": "FL2-A", "dna_pulse": "FL2-H",
            "edu_area": "FL4-A"}


def synthetic_single_rule():
    features = [
        "knn25_log_density_2d__fsc_ssc", "asinh150__dna_area",
        "asinh150__dna_height", "log_ratio__dna_area_height",
        "linear_tls_signed_distance__dna_height_area",
        "knn25_log_density_2d__dna_height_area"]
    return {"feature_names": features, "center": [0] * 6, "scale": [1] * 6,
            "coefficients": [0] * 6, "intercept": 0,
            "threshold": module.EXPECTED_SINGLE_THRESHOLD}


def test_repeated_model_outputs_are_deterministic_and_event_aligned():
    frame = synthetic_frame()
    models = {
        module.G1_MODEL_ID: {"feature_order": module.G1_FEATURES},
        module.EDU_MODEL_ID: {"feature_order": module.EDU_FEATURES},
    }
    first = module.score_acquisition(
        frame, synthetic_roles(), synthetic_single_rule(),
        DeterministicPredictor(), models)
    second = module.score_acquisition(
        frame, synthetic_roles(), synthetic_single_rule(),
        DeterministicPredictor(), models)
    for population in first[0]:
        pd.testing.assert_frame_equal(first[0][population], second[0][population])
        assert first[0][population]["event_identity"].is_unique
    assert first[1] == second[1]
    assert first[1]["g1_subset_predicted_single_cells"] is True
    assert first[1]["edu_positive_subset_predicted_single_cells"] is True
    assert first[1]["parent_event_identity_exact_alignment"] is True


def test_mapping_digest_and_roles_fail_closed():
    approved = {
        "functional_panel": "PANEL_FL2_DNA_FL4_EDU",
        "channel_roles": synthetic_roles(),
        "acquisitions": [{"prefix": "jul28_nt",
            "acquisition_id": "FIG6-EDU-2025-07-28-NT",
            "fcs_basename": "HCT mAC Rif1 CDK1as 7 NT.fcs",
            "fcs_sha256": "a" * 64}],
    }
    item = {"prefix": "jul28_nt",
            "acquisition_id": "FIG6-EDU-2025-07-28-NT",
            "fcs_path": "/outside/HCT mAC Rif1 CDK1as 7 NT.fcs",
            "fcs_sha256": "a" * 64}
    roles, record = module.validate_mapping(item, approved)
    assert roles == synthetic_roles()
    assert record["prefix"] == "jul28_nt"
    item["fcs_sha256"] = "0" * 64
    try:
        module.validate_mapping(item, approved)
    except ValueError as error:
        assert "differs from immutable" in str(error)
    else:
        raise AssertionError("SYNTHETIC mapping mismatch was accepted")


def test_workspace_and_self_attested_mapping_fields_fail_closed():
    approved = {"functional_panel": "PANEL_FL2_DNA_FL4_EDU",
                "channel_roles": synthetic_roles(),
                "acquisitions": [{"prefix": "jul28_nt",
                  "acquisition_id": "A", "fcs_basename": "a.fcs",
                  "fcs_sha256": "b" * 64}]}
    item = {"prefix": "jul28_nt", "acquisition_id": "A",
            "fcs_path": "/outside/a.fcs", "fcs_sha256": "b" * 64,
            "workspace_path": "/forbidden/source.wsp"}
    try:
        module.validate_mapping(item, approved)
    except ValueError as error:
        assert "must not contain workspace" in str(error)
    else:
        raise AssertionError("workspace field was accepted")


def test_identity_digest_is_order_sensitive_and_empty_safe():
    assert module.identity_sha256(["A", "B"]) != module.identity_sha256(["B", "A"])
    assert module.identity_sha256([]) == __import__("hashlib").sha256(b"").hexdigest()


def test_single_cells_artifact_contract_is_hard_pinned_before_loading():
    config = {"artifacts": {
        "single_cells": {"model_id": "SINGLE_CELLS_FROZEN_DEVELOPMENT_RULE",
            "path": "/absent", "byte_sha256": "0" * 64,
            "semantic_sha256": module.EXPECTED_SINGLE_SEMANTIC_SHA256,
            "feature_schema": module.EXPECTED_SINGLE_SCHEMA,
            "threshold": module.EXPECTED_SINGLE_THRESHOLD},
        "g1": {}, "edu_positive": {}}}
    try:
        module.verify_artifacts(config)
    except ValueError as error:
        assert "hard-pinned" in str(error)
    else:
        raise AssertionError("wrong Single Cells digest was accepted")


def test_qc_flags_zero_and_low_support_without_mutating_threshold():
    flags = module.model_qc("g1", np.array([]), np.array([], dtype=bool))
    assert flags == ["g1:LOW_EVENT_SUPPORT", "g1:EXTREME_PREDICTED_FRACTION",
                     "g1:ZERO_PREDICTED_POSITIVES"]
    assert module.G1_THRESHOLD == 0.42
    assert module.EDU_THRESHOLD == 0.385
