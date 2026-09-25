"""SYNTHETIC-only tests for FlowJo sample-population export plans."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path


PYTHON_DIR = Path(__file__).resolve().parents[2] / "python"
sys.path.insert(0, str(PYTHON_DIR))

import pandas as pd  # noqa: E402
from export_flowjo_populations import (  # noqa: E402
    DEFAULT_POPULATIONS, make_flowkit_workspace_copy,
    raw_pnn_aliases, reconcile_planned_counts, validate_sample_population_plan,
)
from lxml import etree  # noqa: E402


class FlowjoExportPlanTests(unittest.TestCase):
    def test_exporter_uses_literal_workspace_edu_gate_name(self):
        self.assertEqual(
            DEFAULT_POPULATIONS[:3],
            ("Single Cells", "EDU Positive", "G1"),
        )

    def test_combined_flowkit_labels_preserve_observed_pnn_aliases(self):
        labels = [
            "FSC-H SYNTHETIC-1", "FSC-A SYNTHETIC-2", "SSC-H SYNTHETIC-3",
            "SSC-A SYNTHETIC-4", "FL2-H PE-H", "FL2-A PE-A",
            "FL4-H APC-H", "FL4-A APC-A", "FSC-Width", "Time",
        ]
        raw = pd.DataFrame([range(len(labels))], columns=labels)
        aliased = raw_pnn_aliases(raw)
        self.assertEqual(list(aliased.columns), [
            "FSC-H", "FSC-A", "SSC-H", "SSC-A", "FL2-H", "FL2-A",
            "FL4-H", "FL4-A", "FSC-Width", "Time",
        ])
        self.assertEqual(list(raw.columns), labels)

    def test_raw_pnn_aliases_reject_duplicate_or_malformed_labels(self):
        valid = ["FSC-H SYNTHETIC", "FL2-A PE-A", "FL4-A APC-A"]
        for labels in (
            [valid[0], "FSC-H DUPLICATE", valid[2]],
            ["invalid label with spaces", *valid[1:]],
        ):
            with self.subTest(labels=labels):
                with self.assertRaisesRegex(RuntimeError, "PnN alias"):
                    raw_pnn_aliases(pd.DataFrame([range(len(labels))], columns=labels))

    def test_temporary_workspace_copy_removes_samples_excluded_by_plan(self):
        xml = b"""<?xml version="1.0" encoding="UTF-8"?>
<Workspace xmlns:gating="http://www.isac-net.org/std/Gating-ML/v2.0/gating"
 xmlns:transforms="http://www.isac-net.org/std/Gating-ML/v2.0/transformations">
<Groups><GroupNode name="All Samples"><Group name="All Samples"><SampleRefs>
  <SampleRef sampleID="1"/><SampleRef sampleID="2"/><SampleRef sampleID="10"/>
</SampleRefs></Group><Population name="SYNTHETIC_GROUP_GATE"/></GroupNode></Groups>
<SampleList>
  <Sample><DataSet sampleID="1"/><SampleNode name="SYNTHETIC_CTRL.fcs" sampleID="1"/></Sample>
  <Sample><DataSet sampleID="2"/><SampleNode name="SYNTHETIC_ASYNC.fcs" sampleID="2"/></Sample>
  <Sample><DataSet sampleID="10"/><SampleNode name="SYNTHETIC_NO_AB.fcs" sampleID="10"/></Sample>
</SampleList></Workspace>
"""
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "source.wsp"
            adapted = Path(directory) / "adapted.wsp"
            source.write_bytes(xml)
            make_flowkit_workspace_copy(
                source, adapted,
                {"SYNTHETIC_CTRL.fcs", "SYNTHETIC_ASYNC.fcs"},
            )
            retained = [
                node.get("name")
                for node in etree.parse(str(adapted)).xpath(
                    "//SampleList/Sample/SampleNode"
                )
            ]
            self.assertEqual(
                set(retained),
                {"SYNTHETIC_CTRL.fcs", "SYNTHETIC_ASYNC.fcs"},
            )
            self.assertNotIn("SYNTHETIC_NO_AB.fcs", retained)
            adapted_tree = etree.parse(str(adapted))
            self.assertEqual(
                [ref.get("sampleID") for ref in adapted_tree.xpath("//SampleRef")],
                ["1", "2"],
            )
            self.assertEqual(
                len(adapted_tree.xpath(
                    "//GroupNode/Population[@name='SYNTHETIC_GROUP_GATE']"
                )),
                1,
            )
            self.assertEqual(source.read_bytes(), xml)
            with self.assertRaisesRegex(RuntimeError, "absent from the FlowJo workspace"):
                make_flowkit_workspace_copy(
                    source, Path(directory) / "must-not-exist.wsp",
                    {"SYNTHETIC_MISSING.fcs"},
                )
            self.assertFalse((Path(directory) / "must-not-exist.wsp").exists())

    def test_temporary_workspace_copy_rejects_unresolved_group_sample_ref(self):
        xml = b"""<Workspace><Groups><GroupNode><Group><SampleRefs>
<SampleRef sampleID="99"/></SampleRefs></Group></GroupNode></Groups><SampleList>
<Sample><DataSet sampleID="1"/><SampleNode name="SYNTHETIC.fcs" sampleID="1"/></Sample>
</SampleList></Workspace>"""
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "source.wsp"
            source.write_bytes(xml)
            with self.assertRaisesRegex(RuntimeError, "do not resolve"):
                make_flowkit_workspace_copy(
                    source, Path(directory) / "adapted.wsp", {"SYNTHETIC.fcs"}
                )

    def test_per_sample_g1_plan_is_accepted(self):
        plan = {
            "SYNTHETIC_CTRL.fcs": ["complete", "edu_positive", "g1"],
            "SYNTHETIC_ASYNC.fcs": ["complete", "edu_positive", "g1"],
        }
        self.assertEqual(
            validate_sample_population_plan(
                plan, ["complete", "edu_positive", "g1"]
            ),
            plan,
        )

    def test_plan_rejects_unknown_or_duplicate_population_keys(self):
        for plan in (
            {"SYNTHETIC.fcs": ["complete", "unknown"]},
            {"SYNTHETIC.fcs": ["complete", "complete"]},
            {"SYNTHETIC.fcs": []},
        ):
            with self.subTest(plan=plan):
                with self.assertRaisesRegex(ValueError, "sample-population plan"):
                    validate_sample_population_plan(
                        plan, ["complete", "edu_positive", "g1"]
                    )

    def test_planned_counts_report_disagreement_without_rejecting_rows(self):
        pairs = {("SYNTHETIC.fcs", "Single Cells"),
                 ("SYNTHETIC.fcs", "EDU+")}
        analysis = pd.DataFrame([
            {"sample_id": pair[0], "gate_name": pair[1], "exported_count": 10}
            for pair in sorted(pairs)
        ])
        saved = pd.DataFrame([
            {"sample_id": pair[0], "gate_name": pair[1], "flowjo_saved_count": 10}
            for pair in sorted(pairs)
        ])
        reconciled = reconcile_planned_counts(analysis, saved, pairs)
        self.assertTrue((reconciled["count_difference"] == 0).all())
        with self.assertRaisesRegex(RuntimeError, "exactly one row"):
            reconcile_planned_counts(pd.concat([analysis, analysis.iloc[[0]]]),
                                     saved, pairs)
        saved.loc[0, "flowjo_saved_count"] = 8
        reconciled = reconcile_planned_counts(analysis, saved, pairs)
        self.assertEqual(reconciled.loc[0, "count_difference"], 2)
        self.assertEqual(reconciled.loc[0, "absolute_count_difference"], 2)
        self.assertEqual(reconciled.loc[0, "signed_percent_difference"], 25.0)
        saved.loc[1, "flowjo_saved_count"] = 0
        reconciled = reconcile_planned_counts(analysis, saved, pairs)
        self.assertTrue(pd.isna(reconciled.loc[1, "signed_percent_difference"]))


if __name__ == "__main__":
    unittest.main()
