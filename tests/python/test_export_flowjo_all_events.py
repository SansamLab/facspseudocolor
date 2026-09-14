"""Pure-logic tests for the general-purpose all-events exporter."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

PYTHON_DIR = Path(__file__).resolve().parents[2] / "python"
sys.path.insert(0, str(PYTHON_DIR))

from export_flowjo_all_events import (  # noqa: E402
    find_channel_column, parse_samples, safe_filename,
)


class SafeFilenameTests(unittest.TestCase):
    def test_replaces_unsafe_characters_and_trims_underscores(self):
        self.assertEqual(safe_filename("ATRi 3 hr!"), "ATRi_3_hr")
        self.assertEqual(safe_filename("__ctrl__"), "ctrl")
        self.assertEqual(safe_filename("plain-name.1"), "plain-name.1")


class FindChannelColumnTests(unittest.TestCase):
    def test_finds_unique_pnn_match_regardless_of_pns(self):
        columns = [("FL2-A", "PE-A"), ("FL4-A", "APC-A"), ("FSC-A", "")]
        self.assertEqual(find_channel_column(columns, "FL2-A"), ("FL2-A", "PE-A"))

    def test_missing_channel_raises_with_available_list(self):
        columns = [("FL2-A", "PE-A")]
        with self.assertRaisesRegex(SystemExit, "not found"):
            find_channel_column(columns, "FL9-A")

    def test_ambiguous_channel_raises(self):
        columns = [("FL2-A", "PE-A"), ("FL2-A", "FITC-A")]
        with self.assertRaisesRegex(SystemExit, "ambiguous"):
            find_channel_column(columns, "FL2-A")


class ParseSamplesTests(unittest.TestCase):
    def test_parses_prefix_and_fcs_filename(self):
        self.assertEqual(
            parse_samples(["ctrl:2_NT.fcs", "atri:6_ATRi.fcs"]),
            [("ctrl", "2_NT.fcs"), ("atri", "6_ATRi.fcs")]
        )

    def test_rejects_missing_colon(self):
        with self.assertRaisesRegex(SystemExit, "PREFIX:FCS_FILENAME"):
            parse_samples(["ctrl_2_NT.fcs"])

    def test_rejects_empty_prefix_or_filename(self):
        with self.assertRaisesRegex(SystemExit, "PREFIX:FCS_FILENAME"):
            parse_samples([":2_NT.fcs"])
        with self.assertRaisesRegex(SystemExit, "PREFIX:FCS_FILENAME"):
            parse_samples(["ctrl:"])

    def test_rejects_duplicate_prefixes(self):
        with self.assertRaisesRegex(SystemExit, "unique"):
            parse_samples(["ctrl:2_NT.fcs", "ctrl:3_NT.fcs"])


class MainRefusesOverwriteTests(unittest.TestCase):
    def test_refuses_to_overwrite_an_existing_output_file(self):
        import export_flowjo_all_events as module

        with tempfile.TemporaryDirectory(prefix="all-events-test-") as directory:
            out_dir = Path(directory) / "out"
            out_dir.mkdir()
            (out_dir / "ctrl_all_events.csv").write_text("existing")
            fcs_dir = Path(directory) / "fcs"
            fcs_dir.mkdir()
            (fcs_dir / "2_NT.fcs").write_bytes(b"not a real FCS file")
            argv = [
                "--fcs-dir", str(fcs_dir), "--output-dir", str(out_dir),
                "--sample", "ctrl:2_NT.fcs",
                "--dna-channel", "FL2-A", "--dna-height-channel", "FL2-H",
            ]
            original_argv = sys.argv
            sys.argv = ["export_flowjo_all_events.py", *argv]
            try:
                with self.assertRaisesRegex(SystemExit, "refusing to overwrite"):
                    module.main()
            finally:
                sys.argv = original_argv


if __name__ == "__main__":
    unittest.main()
