"""SYNTHETIC-only tests for immutable model-derived EdU gate inputs."""

from __future__ import annotations

import hashlib
import json
import sys
import tempfile
import unittest
import warnings
from pathlib import Path
from unittest.mock import Mock, patch

import numpy as np
import pandas as pd


PYTHON_DIR = Path(__file__).resolve().parents[2] / "python"
sys.path.insert(0, str(PYTHON_DIR))

from apply_frozen_gate_models import (  # noqa: E402
    EDU_PORTABLE_SHA256, FEATURES, FROZEN_RULE_HASHES, extract_frozen_features,
    load_rule, model_gate_tables, semantic_hash,
    normalize_flowkit_raw_columns, portable_predict_positive,
    positive_from_gate_indices, validate_output_dir,
    validate_frozen_artifact_config, validate_unloaded_sample_warnings,
    validate_all_event_export_size,
)


def synthetic_events(n: int = 40) -> pd.DataFrame:
    x = np.arange(n, dtype=float)
    return pd.DataFrame({
        "acquisition_id": "SYNTHETIC-A",
        "event_identity": [f"SYNTHETIC-A:event_index:{i}" for i in range(n)],
        "event_index": np.arange(n), "FSC-A": 100 + x * 3,
        "SSC-A": 80 + (x % 9) * 7, "DNA-A": 500 + x * 20,
        "DNA-H": 450 + x * 18, "EDU-A": 20 + x,
    })


ROLES = {"fsc_area": "FSC-A", "ssc_area": "SSC-A", "dna_area": "DNA-A",
         "dna_height": "DNA-H", "edu": "EDU-A"}


def synthetic_rule(target: str, intercept: float, threshold: float = 0.5) -> dict:
    rule = {"feature_names": list(FEATURES), "center": [0.0] * 6,
            "scale": [1.0] * 6, "coefficients": [0.0] * 6,
            "intercept": intercept, "threshold": threshold, "target": target,
            "schema_version": "fd-feature-v2"}
    rule["sha256"] = semantic_hash(rule)
    return rule


def synthetic_edu_model() -> dict:
    """SYNTHETIC one-leaf portable model; never used for scientific output."""
    return {
        "baseline_raw": "0x0.0p+0", "n_features_in": 7,
        "bin_mapper": {"known_categories": [[] for _ in range(7)]},
        "iterations": [[{"nodes": [{"index": 0, "is_leaf": 1,
                                      "value": "0x0.0p+0"}]}]],
    }


class FrozenModelGateTests(unittest.TestCase):
    def test_all_event_export_rejects_oversized_csv(self):
        synthetic_path = Mock()
        synthetic_path.stat.return_value.st_size = 5
        with patch("apply_frozen_gate_models.MAX_ALL_EVENT_BYTES", 4):
            with self.assertRaisesRegex(ValueError, "512 MiB safety limit"):
                validate_all_event_export_size(synthetic_path)

    def test_embedded_portable_traversal_probabilities_and_calls(self):
        """SYNTHETIC deterministic contract for the derived predictor logic."""
        model = {
            "baseline_raw": "0x0.0p+0", "n_features_in": 7,
            "bin_mapper": {"known_categories": [[] for _ in range(7)]},
            "iterations": [[{"raw_left_cat_bitsets": [], "nodes": [
                {"index": 0, "is_leaf": 0, "feature_idx": 0,
                 "is_categorical": 0, "missing_go_to_left": 1,
                 "num_threshold": "0x1.0000000000000p-1", "left": 1,
                 "right": 2},
                {"index": 1, "is_leaf": 1, "value": "0x1.0000000000000p+0"},
                {"index": 2, "is_leaf": 1, "value": "-0x1.0000000000000p+0"},
            ]}]],
        }
        rows = np.asarray([[0.25] + [0.0] * 6, [0.75] + [0.0] * 6])
        observed = portable_predict_positive(model, rows)
        expected = np.asarray([1.0 / (1.0 + np.exp(-1.0)),
                               1.0 / (1.0 + np.exp(1.0))])
        np.testing.assert_allclose(observed, expected, rtol=0.0, atol=1e-15)
        self.assertEqual((observed >= 0.385).astype(int).tolist(), [1, 0])

    def test_artifact_config_rejects_substituted_frozen_rule_hashes(self):
        artifacts = {
            target: {"path": f"/SYNTHETIC/{target}.json", **hashes}
            for target, hashes in FROZEN_RULE_HASHES.items()
        }
        artifacts["edu_positive"] = {
            "path": "/SYNTHETIC/edu.portable.json",
            "byte_sha256": EDU_PORTABLE_SHA256,
        }
        self.assertIsNone(validate_frozen_artifact_config(artifacts))
        for target, field in (("single_cells", "byte_sha256"),
                              ("single_cells", "semantic_sha256"),
                              ("g1", "byte_sha256"),
                              ("g1", "semantic_sha256")):
            changed = json.loads(json.dumps(artifacts))
            changed[target][field] = "0" * 64
            with self.subTest(target=target, field=field):
                with self.assertRaisesRegex(ValueError, "differs from the frozen profile"):
                    validate_frozen_artifact_config(changed)

    def test_features_preserve_event_row_identity(self):
        events = synthetic_events()
        events.index = np.arange(100, 140)
        features = extract_frozen_features(events, ROLES)
        self.assertEqual(list(features.index), list(events.index))
        self.assertEqual(list(features.columns), list(FEATURES))
        self.assertTrue(np.isfinite(features.to_numpy()).all())

    def test_g1_excludes_edu_positive_before_export(self):
        events = synthetic_events()
        single = synthetic_rule("single_cells_reference", intercept=1.0)
        g1 = synthetic_rule("g1_reference", intercept=1.0)
        edu_probability = np.r_[np.ones(20), np.zeros(20)]
        with patch("apply_frozen_gate_models.portable_predict_positive",
                   return_value=edu_probability):
            tables, audit = model_gate_tables(
                events, ROLES, single, g1, synthetic_edu_model())
        self.assertEqual(len(tables["all_events"]), len(events))
        self.assertEqual(len(tables["single_cells"]), len(events))
        self.assertEqual(len(tables["g1"]), 20)
        self.assertEqual(len(tables["edu_positive"]), 20)
        self.assertTrue(audit["g1_subset_single_cells"])
        self.assertTrue(audit["edu_positive_subset_single_cells"])
        self.assertTrue(audit["g1_excludes_edu_positive"])
        self.assertTrue(set(tables["g1"].event_identity).isdisjoint(
            set(tables["edu_positive"].event_identity)))
        self.assertEqual(audit["preliminary_g1_count"], 20)
        self.assertEqual(audit["g1_dna_minimum_fraction"], 0.35)
        self.assertEqual(audit["g1_dna_2n_center_raw"], 1090.0)
        self.assertEqual(audit["g1_dna_minimum_raw"], 381.5)
        self.assertIn("DNA content", tables["single_cells"])
        self.assertIn("EdU", tables["single_cells"])

    def test_missing_direct_identity_fails_closed(self):
        events = synthetic_events()
        events.loc[0, "event_identity"] = "SYNTHETIC-OTHER:event_index:0"
        with self.assertRaisesRegex(ValueError, "does not equal"):
            model_gate_tables(events, ROLES,
                              synthetic_rule("single_cells_reference", 1.0),
                              synthetic_rule("g1_reference", 1.0), synthetic_edu_model())

    def test_mixed_acquisitions_and_noncanonical_identity_fail_closed(self):
        events = synthetic_events()
        events.loc[1, "acquisition_id"] = "SYNTHETIC-B"
        with self.assertRaisesRegex(ValueError, "exactly one"):
            model_gate_tables(events, ROLES,
                              synthetic_rule("single_cells_reference", 1.0),
                              synthetic_rule("g1_reference", 1.0), synthetic_edu_model())
        events = synthetic_events()
        events.loc[1, "event_identity"] = "SYNTHETIC-A:event_index:01"
        with self.assertRaisesRegex(ValueError, "does not equal"):
            model_gate_tables(events, ROLES,
                              synthetic_rule("single_cells_reference", 1.0),
                              synthetic_rule("g1_reference", 1.0), synthetic_edu_model())

    def test_children_are_actually_intersected_when_singlet_rejects(self):
        events = synthetic_events()
        single = synthetic_rule("single_cells_reference", 0.0)
        single["coefficients"][0] = 1.0
        single["center"][0] = float(
            extract_frozen_features(events, ROLES)[FEATURES[0]].median()
        )
        single["sha256"] = semantic_hash(single)
        with patch("apply_frozen_gate_models.portable_predict_positive",
                   return_value=np.zeros(len(events))):
            tables, audit = model_gate_tables(
                events, ROLES, single,
                synthetic_rule("g1_reference", 1.0), synthetic_edu_model()
            )
        self.assertEqual(len(tables["all_events"]), len(events))
        self.assertGreater(len(events), len(tables["single_cells"]))
        self.assertGreater(len(tables["single_cells"]), 0)
        self.assertEqual(len(tables["g1"]), len(tables["single_cells"]))
        self.assertLessEqual(len(tables["edu_positive"]), len(tables["single_cells"]))
        self.assertTrue(audit["g1_subset_single_cells"])
        self.assertTrue(audit["edu_positive_subset_single_cells"])
        self.assertTrue(audit["g1_excludes_edu_positive"])

    def test_artifact_requires_byte_and_semantic_hashes(self):
        with tempfile.TemporaryDirectory(prefix="SYNTHETIC_model_gate_") as directory:
            path = Path(directory) / "SYNTHETIC-rule.json"
            rule = synthetic_rule("g1_reference", 0.0)
            path.write_text(json.dumps(rule), encoding="utf-8")
            byte_hash = hashlib.sha256(path.read_bytes()).hexdigest()
            loaded = load_rule(path, byte_hash, rule["sha256"], "g1_reference", 0.5)
            self.assertEqual(loaded["threshold"], 0.5)
            with self.assertRaisesRegex(ValueError, "byte SHA-256 mismatch"):
                load_rule(path, "0" * 64, rule["sha256"], "g1_reference", 0.5)

    def test_artifact_accepts_decimal_roundoff_but_rejects_threshold_change(self):
        with tempfile.TemporaryDirectory(prefix="SYNTHETIC_model_gate_") as directory:
            path = Path(directory) / "SYNTHETIC-rule.json"
            rule = synthetic_rule("single_cells_reference", 0.0,
                                  0.20499999999999996)
            path.write_text(json.dumps(rule), encoding="utf-8")
            byte_hash = hashlib.sha256(path.read_bytes()).hexdigest()
            loaded = load_rule(path, byte_hash, rule["sha256"],
                               "single_cells_reference", 0.205)
            self.assertEqual(loaded["threshold"], 0.20499999999999996)
            with self.assertRaisesRegex(ValueError, "threshold differs"):
                load_rule(path, byte_hash, rule["sha256"],
                          "single_cells_reference", 0.206)

    def test_output_must_be_new_absolute_and_outside_repository(self):
        with tempfile.TemporaryDirectory(prefix="SYNTHETIC_output_boundary_") as directory:
            repository = Path(directory) / "repository"
            repository.mkdir()
            allowed = Path(directory) / "new-external-output"
            self.assertEqual(validate_output_dir(allowed, repository), allowed.resolve())
            with self.assertRaisesRegex(ValueError, "absolute"):
                validate_output_dir(Path("relative-output"), repository)
            with self.assertRaisesRegex(ValueError, "outside"):
                validate_output_dir(repository, repository)
            with self.assertRaisesRegex(ValueError, "outside"):
                validate_output_dir(repository / "generated", repository)
            allowed.mkdir()
            with self.assertRaisesRegex(ValueError, "existing"):
                validate_output_dir(allowed, repository)

    def test_flowkit_multiindex_columns_use_unique_detector_names(self):
        frame = pd.DataFrame(
            [[1.0, 2.0]],
            columns=pd.MultiIndex.from_tuples([
                ("FSC-A", "Forward Scatter"), ("FL2-A", "DNA content")
            ])
        )
        normalized = normalize_flowkit_raw_columns(frame)
        self.assertEqual(list(normalized.columns), ["FSC-A", "FL2-A"])
        duplicate = frame.copy()
        duplicate.columns = pd.MultiIndex.from_tuples([
            ("FL2-A", "DNA content"), ("FL2-A", "other")
        ])
        with self.assertRaisesRegex(ValueError, "not unique"):
            normalize_flowkit_raw_columns(duplicate)
        three_level = frame.copy()
        three_level.columns = pd.MultiIndex.from_tuples([
            ("FSC-A", "Forward Scatter", "extra"),
            ("FL2-A", "DNA content", "extra")
        ])
        with self.assertRaisesRegex(ValueError, "exactly two"):
            normalize_flowkit_raw_columns(three_level)

    def test_positive_gate_uses_only_indices_into_canonical_all_events(self):
        all_events = pd.DataFrame(
            {"FL2-A": [10.0, 20.0, 30.0], "FL4-A": [1.0, 2.0, 3.0]},
            index=pd.Index([0, 1, 2])
        )
        positive = positive_from_gate_indices(all_events, pd.Index([2, 0]))
        self.assertEqual(list(positive.columns), list(all_events.columns))
        self.assertEqual(list(positive.index), [2, 0])
        self.assertEqual(positive["FL2-A"].tolist(), [30.0, 10.0])
        with self.assertRaisesRegex(ValueError, "duplicate"):
            positive_from_gate_indices(all_events, pd.Index([1, 1]))
        with self.assertRaisesRegex(ValueError, "outside"):
            positive_from_gate_indices(all_events, pd.Index([3]))
        duplicate_parent = all_events.copy()
        duplicate_parent.index = pd.Index([0, 0, 2])
        with self.assertRaisesRegex(ValueError, "duplicate"):
            positive_from_gate_indices(duplicate_parent, pd.Index([0]))

    def test_only_unloaded_other_sample_warnings_are_accepted(self):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            warnings.warn(
                "WSP references OTHER.fcs, but sample was not loaded.", UserWarning
            )
        self.assertEqual(
            validate_unloaded_sample_warnings(caught, "LOADED.fcs"),
            ["WSP references OTHER.fcs, but sample was not loaded."]
        )
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            warnings.warn("compensation differs", UserWarning)
        with self.assertRaisesRegex(RuntimeError, "unexpected FlowKit"):
            validate_unloaded_sample_warnings(caught, "LOADED.fcs")
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            warnings.warn(
                "WSP references LOADED.fcs, but sample was not loaded.", UserWarning
            )
        with self.assertRaisesRegex(RuntimeError, "unexpected FlowKit"):
            validate_unloaded_sample_warnings(caught, "LOADED.fcs")
        for message, category in (
            ("WSP references OTHER.fcs, but sample was not loaded. extra", UserWarning),
            ("WSP references OTHER.fcs but sample was not loaded.", UserWarning),
            ("WSP references OTHER.fcs, but sample was not loaded.", RuntimeWarning),
        ):
            with self.subTest(message=message, category=category):
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter("always")
                    warnings.warn(message, category)
                with self.assertRaisesRegex(RuntimeError, "unexpected FlowKit"):
                    validate_unloaded_sample_warnings(caught, "LOADED.fcs")


if __name__ == "__main__":
    unittest.main()
