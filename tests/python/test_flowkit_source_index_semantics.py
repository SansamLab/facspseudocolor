"""Pinned-environment proof test; accepts only an external SYNTHETIC fixture."""

from __future__ import annotations

import os
import unittest
from pathlib import Path

import flowkit as fk


EXPECTED_FILES = (
    "SYNTHETIC_INDEX_ACQUISITION_1.fcs",
    "SYNTHETIC_INDEX_ACQUISITION_2.fcs",
)
EXPECTED_GENERATOR = "SYNTHETIC_flowjo_index_fixture_v1"
EXPECTED_CONSTRUCTION = "SYNTHETIC_DETERMINISTIC_ARITHMETIC_V1_NO_RANDOM_SEED"
EXPECTED_MEASUREMENT_COLUMNS = (
    "SYNTHETIC_FSC_A",
    "SYNTHETIC_FSC_H",
    "SYNTHETIC_DNA_A",
    "SYNTHETIC_PH3_A",
)


class FlowKitSourceIndexSemanticsTests(unittest.TestCase):
    def test_get_gate_events_indices_are_stable_parent_to_children(self):
        workspace_path = Path(os.environ["SYNTHETIC_FLOWJO_WORKSPACE"])
        fcs_dir = Path(os.environ["SYNTHETIC_FLOWJO_FCS_DIR"])
        self.assertIn("SYNTHETIC", str(workspace_path).upper())
        self.assertIn("SYNTHETIC", str(fcs_dir).upper())
        fcs_paths = sorted(fcs_dir.glob("*.fcs"))
        self.assertEqual(tuple(path.name for path in fcs_paths), EXPECTED_FILES,
                         "fixture must contain exactly the two generated acquisitions")
        fcs_files = [str(path) for path in fcs_paths]
        workspace = fk.Workspace(str(workspace_path), fcs_samples=fcs_files,
                                 filename_as_id=True)
        workspace.analyze_samples(use_mp=False)
        gate_names = {
            "single": os.environ["SYNTHETIC_SINGLE_CELLS_GATE"],
            "g1": os.environ["SYNTHETIC_G1_GATE"],
            "positive": os.environ["SYNTHETIC_PH3_POSITIVE_GATE"],
        }
        scoped = set()
        self.assertEqual(set(workspace.get_sample_ids()), set(EXPECTED_FILES))
        for sample_id in workspace.get_sample_ids():
            source_sample = workspace.get_sample(sample_id)
            metadata = source_sample.get_metadata()
            expected_acquisition = Path(sample_id).stem
            self.assertEqual(metadata.get("synthetic_fixture"),
                             "SYNTHETIC TEST FIXTURE - NOT EXPERIMENTAL DATA")
            self.assertEqual(metadata.get("synthetic_generator"), EXPECTED_GENERATOR)
            self.assertEqual(metadata.get("synthetic_construction"), EXPECTED_CONSTRUCTION)
            self.assertEqual(metadata.get("synthetic_acquisition_id"), expected_acquisition)
            self.assertEqual(metadata.get("synthetic_interpretation"),
                             "SYNTHETIC - NO BIOLOGICAL INTERPRETATION")
            self.assertIn("SYNTHETIC", metadata.get("synthetic_provenance", ""))
            source_raw = source_sample.as_dataframe(source="raw", col_multi_index=False)
            self.assertEqual(tuple(source_raw.columns), EXPECTED_MEASUREMENT_COLUMNS)
            indices = {}
            gate_paths = {}
            for key, gate_name in gate_names.items():
                paths = workspace.find_matching_gate_paths(sample_id, gate_name)
                self.assertEqual(len(paths), 1)
                gate_paths[key] = tuple(paths[0])
                raw = workspace.get_gate_events(sample_id, gate_name, paths[0], source="raw")
                transformed = workspace.get_gate_events(
                    sample_id, gate_name, paths[0], source="xform"
                )
                raw_bookkeeping = [column for column in raw.columns
                                   if column not in source_raw.columns]
                transformed_bookkeeping = [column for column in transformed.columns
                                           if column not in source_raw.columns]
                self.assertEqual(raw_bookkeeping, ["sample_id"])
                self.assertEqual(transformed_bookkeeping, ["sample_id"])
                self.assertTrue((raw["sample_id"] == sample_id).all())
                self.assertTrue((transformed["sample_id"] == sample_id).all())
                measurement_columns = [column for column in raw.columns
                                       if column != "sample_id"]
                self.assertEqual(tuple(measurement_columns), EXPECTED_MEASUREMENT_COLUMNS)
                self.assertEqual(measurement_columns,
                                 [column for column in transformed.columns
                                  if column != "sample_id"])
                self.assertEqual(list(map(str, raw.index)),
                                 list(map(str, transformed.index)))
                self.assertEqual(len(raw.index), len(set(map(str, raw.index))))
                source_rows = source_raw.loc[raw.index, measurement_columns]
                gated_measurements = raw.loc[:, measurement_columns]
                self.assertEqual(list(source_rows.index), list(raw.index))
                self.assertTrue(source_rows.equals(gated_measurements),
                                "gated raw rows must equal source-acquisition rows at the index")
                indices[key] = set(map(str, raw.index))
            expected_child_parent = (*gate_paths["single"], gate_names["single"])
            self.assertEqual(gate_paths["g1"], expected_child_parent)
            self.assertEqual(gate_paths["positive"], expected_child_parent)
            self.assertTrue(indices["g1"] < indices["single"])
            self.assertTrue(indices["positive"] < indices["single"])
            self.assertGreater(len(indices["g1"]), 0)
            self.assertGreater(len(indices["positive"]), 0)
            self.assertLess(len(indices["single"]), len(source_raw))
            scoped.update((sample_id, index) for index in indices["single"])
        self.assertEqual(len(scoped), sum(
            len(workspace.get_gate_events(
                sample_id, gate_names["single"],
                workspace.find_matching_gate_paths(sample_id, gate_names["single"])[0],
                source="raw",
            )) for sample_id in workspace.get_sample_ids()
        ))


if __name__ == "__main__":
    unittest.main()
