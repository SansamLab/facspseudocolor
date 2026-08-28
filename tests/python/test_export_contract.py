"""SYNTHETIC-only tests for the Slice 1 local export contract."""

from __future__ import annotations

import sys
import tempfile
import unittest
import hashlib
from pathlib import Path


PYTHON_DIR = Path(__file__).resolve().parents[2] / "python"
sys.path.insert(0, str(PYTHON_DIR))

from export_contract import (  # noqa: E402
    LEGACY_PROFILE, PRODUCTION_PROFILE, _atomic_write, artifact_record,
    event_identity_fields,
    canonical_manifest_bytes, finalize_geometry_supplement, finalize_manifest,
    manifest_binding_digest, new_manifest, resolve_production_fcs_files,
    scoped_event_identity, validate_artifact_prefix, validate_geometry_coverage,
    validate_manifest_structure,
    validate_geometry_artifact_scope, verify_finalized_manifest,
    verify_population_artifacts,
)


def synthetic_metadata() -> dict:
    return {
        "software": {"supported_flowkit_version": "SYNTHETIC-PINNED"},
        "approval": {
            "gate_owner": "SYNTHETIC OWNER", "approver": "SYNTHETIC APPROVER",
            "approval_date": "2026-08-24",
            "approval_record": "SYNTHETIC/approval.md",
            "positivity_method_id": "SYNTHETIC_METHOD",
            "positivity_method_version": "1",
        },
        "transformation_context": "SYNTHETIC explicit transform context",
        "compensation_context": "SYNTHETIC uncompensated",
        "acquisitions": [{
            "acquisition_id": "SYNTHETIC-A", "sample_id": "SYNTHETIC-A.fcs",
            "prefix": "SYNTHETIC_A", "source_fcs_reference": "SYNTHETIC-A.fcs",
            "source_fcs_sha256": "0" * 64,
        }],
    }


class ExportContractTests(unittest.TestCase):
    def test_atomic_publication_never_clobbers_existing_path(self):
        with tempfile.TemporaryDirectory(prefix="SYNTHETIC_contract_") as directory:
            target = Path(directory) / "SYNTHETIC-immutable.json"
            target.write_bytes(b"SYNTHETIC ORIGINAL")
            with self.assertRaises(FileExistsError):
                _atomic_write(target, b"SYNTHETIC CONCURRENT REPLACEMENT")
            self.assertEqual(target.read_bytes(), b"SYNTHETIC ORIGINAL")
            self.assertEqual(list(Path(directory).iterdir()), [target])

    def test_acquisition_prefix_is_path_free_safe_filename(self):
        self.assertEqual(validate_artifact_prefix("SYNTHETIC_A-1.2"),
                         "SYNTHETIC_A-1.2")
        for value in ("", ".", "..", "../escape", "nested/name", "/absolute"):
            with self.subTest(value=value):
                with self.assertRaisesRegex(ValueError, "path-free safe"):
                    validate_artifact_prefix(value)

    def test_direct_identity_is_stable_across_parent_and_children_and_scoped(self):
        parent = [scoped_event_identity("SYNTHETIC-A", x) for x in ("0", "4", "8")]
        g1 = [scoped_event_identity("SYNTHETIC-A", x) for x in ("0", "4")]
        positive = [scoped_event_identity("SYNTHETIC-A", "8")]
        self.assertTrue(set(g1).issubset(parent))
        self.assertTrue(set(positive).issubset(parent))
        self.assertNotEqual(scoped_event_identity("SYNTHETIC-B", "0"), parent[0])
        self.assertEqual(parent[0], "SYNTHETIC-A:event_index:0")

    def test_missing_source_identity_has_no_sequence_fallback(self):
        with self.assertRaisesRegex(ValueError, "source event index"):
            scoped_event_identity("SYNTHETIC-A", "")

    def test_legacy_identity_fields_are_blank(self):
        fields = event_identity_fields(LEGACY_PROFILE, "SYNTHETIC-A", ["0", "1"])
        self.assertTrue(all(value == "" for values in fields.values() for value in values))

    def test_exact_fcs_hash_binding_and_workspace_coverage(self):
        with tempfile.TemporaryDirectory(prefix="SYNTHETIC_contract_") as directory:
            root = Path(directory)
            source = root / "SYNTHETIC-A.fcs"
            source.write_bytes(b"SYNTHETIC FCS BYTES")
            metadata = synthetic_metadata()["acquisitions"]
            metadata[0]["source_fcs_sha256"] = hashlib.sha256(source.read_bytes()).hexdigest()
            self.assertEqual(resolve_production_fcs_files(
                root, metadata, {"SYNTHETIC-A.fcs"}), [source.resolve()])
            metadata[0]["source_fcs_sha256"] = "0" * 64
            with self.assertRaisesRegex(ValueError, "SHA-256 mismatch"):
                resolve_production_fcs_files(root, metadata, {"SYNTHETIC-A.fcs"})

    def test_production_is_conditional_and_required_metadata_is_fail_closed(self):
        with tempfile.TemporaryDirectory(prefix="SYNTHETIC_contract_") as directory:
            workspace = Path(directory) / "SYNTHETIC.wsp"
            workspace.write_text("SYNTHETIC WORKSPACE", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "conditional"):
                new_manifest(operation_id="SYNTHETIC-OP", profile=PRODUCTION_PROFILE,
                             workspace=workspace, flowkit_version="SYNTHETIC-PINNED",
                             metadata=synthetic_metadata(),
                             direct_index_semantics_verified=False)

    def test_metadata_allowlist_rejects_credentials_and_unknown_fields(self):
        with tempfile.TemporaryDirectory(prefix="SYNTHETIC_contract_") as directory:
            workspace = Path(directory) / "SYNTHETIC.wsp"
            workspace.write_text("SYNTHETIC WORKSPACE", encoding="utf-8")
            metadata = synthetic_metadata()
            metadata["credentials"] = "SYNTHETIC SECRET MUST NOT SERIALIZE"
            with self.assertRaisesRegex(ValueError, "unsupported metadata field"):
                new_manifest(operation_id="SYNTHETIC-OP", profile=PRODUCTION_PROFILE,
                             workspace=workspace, flowkit_version="SYNTHETIC-PINNED",
                             metadata=metadata, direct_index_semantics_verified=True)

    def test_manifest_digest_artifact_mutation_and_no_overwrite(self):
        with tempfile.TemporaryDirectory(prefix="SYNTHETIC_contract_") as directory:
            root = Path(directory)
            workspace = root / "SYNTHETIC.wsp"
            workspace.write_text("SYNTHETIC WORKSPACE", encoding="utf-8")
            artifact = root / "SYNTHETIC.csv"
            artifact.write_text("event_identity\nSYNTHETIC-A:event_index:0\n", encoding="utf-8")
            manifest = new_manifest(
                operation_id="SYNTHETIC-OP", profile=LEGACY_PROFILE,
                workspace=workspace, flowkit_version="SYNTHETIC-PINNED",
                metadata={}, direct_index_semantics_verified=False,
            )
            manifest["artifacts"].append(artifact_record(
                artifact, operation_dir=root, role="population_events", row_count=1,
                columns=["event_identity"], identity_columns=["event_identity"],
                linkage={"export_operation_id": "SYNTHETIC-OP"},
            ))
            finalize_manifest(manifest, root)
            self.assertEqual(verify_finalized_manifest(root)["status"], "complete")
            with self.assertRaises(FileExistsError):
                finalize_manifest(manifest, root)
            artifact.write_text("MUTATED", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "artifact hash mismatch"):
                verify_finalized_manifest(root)

    def test_geometry_linkage_rejects_filename_only_and_operation_mismatch(self):
        with tempfile.TemporaryDirectory(prefix="SYNTHETIC_contract_") as directory:
            root = Path(directory)
            workspace = root / "SYNTHETIC.wsp"
            workspace.write_text("SYNTHETIC WORKSPACE", encoding="utf-8")
            artifact = root / "same-looking-name.csv"
            artifact.write_text("SYNTHETIC", encoding="utf-8")
            manifest = new_manifest(
                operation_id="SYNTHETIC-OP", profile=LEGACY_PROFILE,
                workspace=workspace, flowkit_version="SYNTHETIC-PINNED",
                metadata={}, direct_index_semantics_verified=False,
            )
            finalize_manifest(manifest, root)
            supplement = {
                "export_operation_id": "SYNTHETIC-WRONG",
                "workspace_sha256": manifest["workspace"]["sha256"],
                "artifacts": [],
            }
            with self.assertRaisesRegex(ValueError, "operation ID"):
                finalize_geometry_supplement(supplement, root)

    def test_canonicalization_is_deterministic(self):
        self.assertEqual(canonical_manifest_bytes({"b": 2, "a": 1}),
                         canonical_manifest_bytes({"a": 1, "b": 2}))

    def test_manifest_binding_digest_is_stable_and_distinct_from_full_manifest(self):
        with tempfile.TemporaryDirectory(prefix="SYNTHETIC_contract_") as directory:
            workspace = Path(directory) / "SYNTHETIC.wsp"
            workspace.write_text("SYNTHETIC WORKSPACE", encoding="utf-8")
            manifest = new_manifest(operation_id="SYNTHETIC-OP", profile=LEGACY_PROFILE,
                                    workspace=workspace, flowkit_version="SYNTHETIC-PINNED",
                                    metadata={}, direct_index_semantics_verified=False)
            binding = manifest_binding_digest(manifest)
            manifest["artifacts"].append({"SYNTHETIC": "ledger change"})
            self.assertEqual(manifest_binding_digest(manifest), binding)
            self.assertNotEqual(binding,
                                hashlib.sha256(canonical_manifest_bytes(manifest)).hexdigest())

    def test_production_schema_accepts_gate_present_zero_and_rejects_mixed_missing(self):
        with tempfile.TemporaryDirectory(prefix="SYNTHETIC_contract_") as directory:
            workspace = Path(directory) / "SYNTHETIC.wsp"
            workspace.write_text("SYNTHETIC WORKSPACE", encoding="utf-8")
            metadata = synthetic_metadata()
            metadata["acquisitions"].append({**metadata["acquisitions"][0],
                                              "acquisition_id": "SYNTHETIC-B",
                                              "sample_id": "SYNTHETIC-B.fcs",
                                              "prefix": "SYNTHETIC_B"})
            manifest = new_manifest(
                operation_id="SYNTHETIC-OP", profile=PRODUCTION_PROFILE,
                workspace=workspace, flowkit_version="SYNTHETIC-PINNED", metadata=metadata,
                direct_index_semantics_verified=True, requested_populations=["positive"],
            )
            artifact_fields = {
                "role": "population_events", "path": "SYNTHETIC.csv", "sha256": "0" * 64,
                "byte_size": 1, "row_count": 0, "columns": [], "identity_columns": [],
                "identity_method_id": None, "intentionally_empty": True,
                "export_operation_id": "SYNTHETIC-OP", "acquisition_id": "SYNTHETIC-A",
                "sample_id": "SYNTHETIC-A.fcs", "population_key": "positive",
                "gate_path": "root/positive", "channels": [],
            }
            manifest["artifacts"] = [artifact_fields]
            population = {
                "population_key": "positive", "gate_name": "positive",
                "gate_path": "root/positive", "gate_type": "PolygonGate",
                "parent_population_path": "root", "acquisition_id": "SYNTHETIC-A",
                "sample_id": "SYNTHETIC-A.fcs", "prefix": "SYNTHETIC_A",
                "channels": [], "gate_channels": ["DNA", "pH3"], "row_count": 0,
                "identity_field": "event_identity", "identity_method_id": "direct",
                "unique_identity_count": 0, "duplicate_row_count": 0,
                "intentionally_empty": True, "export_operation_id": "SYNTHETIC-OP",
                "artifact_path": "SYNTHETIC.csv", "artifact_sha256": "0" * 64,
            }
            manifest["populations"] = [population]
            with self.assertRaisesRegex(ValueError, "exactly cover"):
                validate_manifest_structure(manifest)
            second = dict(population, acquisition_id="SYNTHETIC-B",
                          sample_id="SYNTHETIC-B.fcs", prefix="SYNTHETIC_B")
            manifest["populations"].append(second)
            manifest["artifacts"].append(dict(
                artifact_fields, path="SYNTHETIC_B.csv", acquisition_id="SYNTHETIC-B",
                sample_id="SYNTHETIC-B.fcs",
            ))
            validate_manifest_structure(manifest)

    def test_geometry_requires_exact_sample_gate_and_ordered_channels(self):
        populations = [{"population_key": "positive", "sample_id": "SYNTHETIC-A.fcs",
                        "gate_path": "root/positive", "gate_channels": ["DNA", "pH3"]}]
        acquisitions = [{"sample_id": "SYNTHETIC-A.fcs"}]
        validate_geometry_coverage(populations, acquisitions, "positive",
                                   {("SYNTHETIC-A.fcs", "root/positive", ("DNA", "pH3"))})
        with self.assertRaisesRegex(ValueError, "sample/gate/channels"):
            validate_geometry_coverage(populations, acquisitions, "positive",
                                       {("SYNTHETIC-A.fcs", "root/positive", ("pH3", "DNA"))})

    def test_geometry_artifact_scope_reconciles_all_rows(self):
        observed = {
            ("SYNTHETIC-A", "SYNTHETIC-A.fcs", "root/positive", ("DNA", "pH3")),
            ("SYNTHETIC-B", "SYNTHETIC-B.fcs", "root/positive", ("DNA", "pH3")),
        }
        scope = [
            {"acquisition_id": item[0], "sample_id": item[1], "gate_path": item[2],
             "channels": list(item[3])} for item in sorted(observed)
        ]
        validate_geometry_artifact_scope(scope, observed)
        with self.assertRaisesRegex(ValueError, "linkage scope"):
            validate_geometry_artifact_scope(scope[:-1], observed)

    def test_consumption_verifier_reconciles_rows_and_detects_mutation(self):
        with tempfile.TemporaryDirectory(prefix="SYNTHETIC_contract_") as directory:
            root = Path(directory)
            workspace = root / "SYNTHETIC.wsp"
            workspace.write_text("SYNTHETIC WORKSPACE", encoding="utf-8")
            manifest = new_manifest(operation_id="SYNTHETIC-OP", profile=LEGACY_PROFILE,
                                    workspace=workspace, flowkit_version="SYNTHETIC-PINNED",
                                    metadata={}, direct_index_semantics_verified=False)
            artifact = root / "SYNTHETIC.csv"
            columns = ["export_operation_id", "export_manifest_digest", "acquisition_id",
                       "sample_id", "export_profile"]
            artifact.write_text(
                ",".join(columns) + "\n" +
                ",".join(["SYNTHETIC-OP", manifest["manifest_binding"]["digest"], "",
                          "SYNTHETIC-A.fcs", LEGACY_PROFILE]) + "\n",
                encoding="utf-8",
            )
            manifest["artifacts"].append(artifact_record(
                artifact, operation_dir=root, role="population_events", row_count=1,
                columns=columns, linkage={"export_operation_id": "SYNTHETIC-OP",
                                          "acquisition_id": "", "sample_id": "SYNTHETIC-A.fcs",
                                          "population_key": "positive"},
            ))
            finalize_manifest(manifest, root)
            verify_population_artifacts(root, [artifact])
            artifact.write_text("MUTATED", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "hash mismatch"):
                verify_population_artifacts(root, [artifact])

    def test_canonicalization_rejects_nonfinite_json(self):
        with self.assertRaises(ValueError):
            canonical_manifest_bytes({"not_json": float("nan")})

    def test_manifest_verification_rejects_schema_mismatch_and_path_escape(self):
        with tempfile.TemporaryDirectory(prefix="SYNTHETIC_contract_") as directory:
            root = Path(directory)
            workspace = root / "SYNTHETIC.wsp"
            workspace.write_text("SYNTHETIC WORKSPACE", encoding="utf-8")
            manifest = new_manifest(operation_id="SYNTHETIC-OP", profile=LEGACY_PROFILE,
                                    workspace=workspace, flowkit_version="SYNTHETIC-PINNED",
                                    metadata={}, direct_index_semantics_verified=False)
            manifest["manifest_schema"]["version"] = "WRONG"
            with self.assertRaisesRegex(ValueError, "schema mismatch"):
                finalize_manifest(manifest, root)
            manifest["manifest_schema"]["version"] = "1.0.0"
            manifest["artifacts"].append({
                "path": "../escape.csv", "export_operation_id": "SYNTHETIC-OP",
                "sha256": "0" * 64, "role": "population_events", "byte_size": 1,
                "row_count": 1, "columns": [], "identity_columns": [],
                "identity_method_id": None, "intentionally_empty": False,
                "acquisition_id": None, "sample_id": None, "population_key": None,
                "gate_path": None, "channels": [],
            })
            with self.assertRaisesRegex(ValueError, "escapes"):
                finalize_manifest(manifest, root)


if __name__ == "__main__":
    unittest.main()
