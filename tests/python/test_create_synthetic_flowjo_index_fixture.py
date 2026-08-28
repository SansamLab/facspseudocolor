"""Pure-helper and safety tests for the SYNTHETIC FlowJo fixture generator."""

from __future__ import annotations

import importlib.util
import tempfile
import types
import unittest
from unittest import mock
from pathlib import Path


GENERATOR = Path(__file__).parents[2] / "tools" / "create_synthetic_flowjo_index_fixture.py"
SPEC = importlib.util.spec_from_file_location("synthetic_flowjo_fixture", GENERATOR)
fixture = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(fixture)


class SyntheticFlowJoFixtureTests(unittest.TestCase):
    def test_arrays_are_deterministic_labeled_and_distinct(self):
        first = fixture.build_synthetic_events(1)
        second = fixture.build_synthetic_events(2)
        self.assertEqual(first, fixture.build_synthetic_events(1))
        self.assertEqual(len(first), 360)
        self.assertEqual(len(second), 360)
        self.assertNotEqual(first, second)
        self.assertEqual(len(fixture.CHANNELS), 4)
        self.assertTrue(all("SYNTHETIC" in name for name in fixture.CHANNELS))
        self.assertTrue(all("SYNTHETIC" in name for name in fixture.SAMPLE_IDS))

    def test_intended_children_are_nonempty_and_inside_intended_parent(self):
        for acquisition in (1, 2):
            events = fixture.build_synthetic_events(acquisition)
            single = [event for event in events if 40000 < event[0] < 70000 and
                      0.94 < event[1] / event[0] < 1.06]
            g1 = [event for event in single if 45000 < event[2] < 65000]
            positive = [event for event in single if event[3] > 70000]
            self.assertEqual(len(single), 320)
            self.assertEqual(len(g1), 130)
            self.assertEqual(len(positive), 70)

    def test_destination_must_be_explicitly_synthetic(self):
        with self.assertRaisesRegex(ValueError, "must contain SYNTHETIC"):
            fixture.require_synthetic_destination(Path("ordinary_fixture"))
        accepted = fixture.require_synthetic_destination(Path("SYNTHETIC_fixture"))
        self.assertIn("SYNTHETIC", str(accepted).upper())

    def test_existing_outputs_are_refused_before_flowio_is_used(self):
        with tempfile.TemporaryDirectory(prefix="SYNTHETIC_fixture_test_") as directory:
            destination = Path(directory)
            existing = fixture.output_paths(destination)[0]
            existing.write_bytes(b"SYNTHETIC PLACEHOLDER")
            with self.assertRaisesRegex(FileExistsError, "refusing to overwrite"):
                fixture.write_fixture(destination)

    def test_second_writer_failure_leaves_no_partial_fixture(self):
        calls = []

        def failing_writer(handle, event_data, channel_names, metadata_dict=None):
            calls.append(metadata_dict["synthetic_acquisition_id"])
            if len(calls) == 2:
                raise RuntimeError("SYNTHETIC deliberate writer failure")
            handle.write(b"SYNTHETIC STAGED FCS")

        fake_flowio = types.SimpleNamespace(create_fcs=failing_writer)
        with tempfile.TemporaryDirectory(prefix="SYNTHETIC_fixture_test_") as directory:
            destination = Path(directory)
            with mock.patch.dict("sys.modules", {"flowio": fake_flowio}):
                with self.assertRaisesRegex(RuntimeError, "deliberate writer failure"):
                    fixture.write_fixture(destination)
            self.assertFalse(any(path.exists() for path in fixture.output_paths(destination)))
            self.assertEqual(list(destination.iterdir()), [])

    def test_overwrite_staging_failure_preserves_both_existing_files(self):
        calls = []

        def failing_writer(handle, event_data, channel_names, metadata_dict=None):
            calls.append(metadata_dict["synthetic_acquisition_id"])
            if len(calls) == 2:
                raise RuntimeError("SYNTHETIC deliberate overwrite staging failure")
            handle.write(b"SYNTHETIC REPLACEMENT")

        fake_flowio = types.SimpleNamespace(create_fcs=failing_writer)
        with tempfile.TemporaryDirectory(prefix="SYNTHETIC_fixture_test_") as directory:
            destination = Path(directory)
            paths = fixture.output_paths(destination)
            paths[0].write_bytes(b"SYNTHETIC ORIGINAL ONE")
            paths[1].write_bytes(b"SYNTHETIC ORIGINAL TWO")
            with mock.patch.dict("sys.modules", {"flowio": fake_flowio}):
                with self.assertRaisesRegex(RuntimeError, "overwrite staging failure"):
                    fixture.write_fixture(destination, overwrite=True)
            self.assertEqual(paths[0].read_bytes(), b"SYNTHETIC ORIGINAL ONE")
            self.assertEqual(paths[1].read_bytes(), b"SYNTHETIC ORIGINAL TWO")
            self.assertEqual(set(destination.iterdir()), set(paths))

    def test_restore_failure_retains_and_reports_recoverable_backup(self):
        def writer(handle, event_data, channel_names, metadata_dict=None):
            handle.write(b"SYNTHETIC REPLACEMENT")

        fake_flowio = types.SimpleNamespace(create_fcs=writer)
        with tempfile.TemporaryDirectory(prefix="SYNTHETIC_fixture_test_") as directory:
            # macOS may present /var while Path.resolve() canonicalizes it to
            # /private/var; use the generator's canonical destination so the
            # injected target comparisons match the paths write_fixture uses.
            destination = fixture.require_synthetic_destination(Path(directory))
            paths = fixture.output_paths(destination)
            paths[0].write_bytes(b"SYNTHETIC ORIGINAL ONE")
            paths[1].write_bytes(b"SYNTHETIC ORIGINAL TWO")
            real_replace = fixture.replace_file

            def fail_install_and_one_restore(source, target):
                source = Path(source)
                target = Path(target)
                if source.suffix == ".tmp" and target == paths[1]:
                    raise OSError("SYNTHETIC deliberate install failure")
                if source.suffix == ".backup" and target == paths[0]:
                    raise OSError("SYNTHETIC deliberate restore failure")
                return real_replace(source, target)

            with mock.patch.dict("sys.modules", {"flowio": fake_flowio}), \
                    mock.patch.object(fixture, "replace_file", fail_install_and_one_restore):
                with self.assertRaisesRegex(RuntimeError, "do not delete retained backup") as caught:
                    fixture.write_fixture(destination, overwrite=True)
            retained = list(destination.glob(".*.backup"))
            self.assertEqual(len(retained), 1)
            self.assertIn(str(retained[0]), str(caught.exception))
            self.assertEqual(retained[0].read_bytes(), b"SYNTHETIC ORIGINAL ONE")
            self.assertFalse(paths[0].exists())
            self.assertEqual(paths[1].read_bytes(), b"SYNTHETIC ORIGINAL TWO")


if __name__ == "__main__":
    unittest.main()
