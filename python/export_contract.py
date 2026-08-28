"""Local, fail-closed export-operation contract shared by FlowJo exporters.

The manifest digest is stored in ``export-manifest.sha256``.  It is the
SHA-256 of the UTF-8 JSON bytes produced by :func:`canonical_manifest_bytes`;
the digest is deliberately not embedded in the hashed JSON.
"""

from __future__ import annotations

import hashlib
import argparse
import csv
import json
import os
import re
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


SCHEMA_NAME = "facspseudocolor-flowjo-export"
SCHEMA_VERSION = "1.0.0"
EXPORT_PROFILE_ID = "ph3_export_identity_provenance_v1"
EXPORT_PROFILE_VERSION = "1.0.0"
DIRECT_METHOD_ID = "flowkit_source_event_index_scoped_v1"
DIRECT_METHOD_VERSION = "1.0.0"
LEGACY_PROFILE = "legacy_count_only_unverified_v1"
PRODUCTION_PROFILE = "production_direct_identity_v1"

REQUIRED_APPROVAL = (
    "gate_owner", "approver", "approval_date", "approval_record",
    "positivity_method_id", "positivity_method_version",
)
REQUIRED_ACQUISITION = (
    "acquisition_id", "sample_id", "prefix", "source_fcs_reference",
    "source_fcs_sha256",
)
ALLOWED_TOP_LEVEL = {
    "approval", "acquisitions", "software", "transformation_context",
    "compensation_context",
}
ALLOWED_APPROVAL = set(REQUIRED_APPROVAL)
ALLOWED_SOFTWARE = {
    "supported_flowkit_version", "exporter", "exporter_version",
    "source_commit", "flowjo_version",
}
ALLOWED_ACQUISITION = set(REQUIRED_ACQUISITION)


def validate_artifact_prefix(prefix: str) -> str:
    """Require one non-special, path-free filename component from metadata."""
    if (not isinstance(prefix, str) or prefix in ("", ".", "..") or
            re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", prefix) is None or
            Path(prefix).name != prefix):
        raise ValueError(
            "acquisition prefix must be a path-free safe filename component"
        )
    return prefix


def scoped_event_identity(acquisition_id: str, source_index: Any) -> str:
    """Construct direct identity without coercing or renumbering source indices."""
    if not isinstance(acquisition_id, str) or not acquisition_id:
        raise ValueError("acquisition_id is required")
    token = str(source_index)
    if not token:
        raise ValueError("source event index is required")
    return f"{acquisition_id}:event_index:{token}"


def event_identity_fields(profile: str, acquisition_id: str,
                          source_indices: list[str]) -> dict[str, list[str]]:
    production = profile == PRODUCTION_PROFILE
    blank = [""] * len(source_indices)
    return {
        "acquisition_id": [acquisition_id] * len(source_indices) if production else blank,
        "event_index": source_indices if production else blank,
        "event_identity": [scoped_event_identity(acquisition_id, item)
                           for item in source_indices] if production else blank,
        "identity_source": ["flowkit_get_gate_events_index"] * len(source_indices)
                           if production else blank,
        "identity_method_id": [DIRECT_METHOD_ID] * len(source_indices) if production else blank,
        "identity_method_version": [DIRECT_METHOD_VERSION] * len(source_indices)
                                   if production else blank,
        "duplicate_occurrence": ["1"] * len(source_indices) if production else blank,
    }


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def resolve_production_fcs_files(fcs_dir: Path, acquisitions: list[dict[str, Any]],
                                 workspace_sample_ids: set[str]) -> list[Path]:
    by_sample = {item["sample_id"]: item for item in acquisitions}
    if set(by_sample) != workspace_sample_ids:
        raise ValueError("production acquisition mapping does not exactly cover workspace samples")
    root = fcs_dir.resolve()
    paths = []
    for sample_id in sorted(workspace_sample_ids):
        acquisition = by_sample[sample_id]
        reference = Path(acquisition["source_fcs_reference"])
        if reference.is_absolute() or ".." in reference.parts:
            raise ValueError("production source_fcs_reference must be relative and confined")
        path = (root / reference).resolve()
        if root not in path.parents or not path.is_file() or path.name != sample_id:
            raise ValueError(f"explicit source FCS does not match workspace sample {sample_id!r}")
        if sha256_file(path) != acquisition["source_fcs_sha256"]:
            raise ValueError(f"source FCS SHA-256 mismatch for {sample_id!r}")
        paths.append(path)
    return paths


def validate_geometry_coverage(populations: list[dict[str, Any]],
                               acquisitions: list[dict[str, Any]],
                               population_key: str,
                               observed: set[tuple[str, str, tuple[str, ...]]]) -> None:
    expected_populations = [item for item in populations
                            if item["population_key"] == population_key]
    expected_samples = {item["sample_id"] for item in acquisitions}
    if {item[0] for item in observed} != expected_samples:
        raise ValueError("geometry sample coverage does not match the export operation")
    for sample_id, gate_path, channels in observed:
        matches = [item for item in expected_populations
                   if item["sample_id"] == sample_id and item["gate_path"] == gate_path]
        if len(matches) != 1 or tuple(matches[0]["gate_channels"]) != channels:
            raise ValueError("geometry sample/gate/channels do not match population manifest")


def validate_geometry_artifact_scope(scope: list[dict[str, Any]],
                                     observed: set[tuple[str, str, str, tuple[str, ...]]]) -> None:
    declared = {
        (item["acquisition_id"], item["sample_id"], item["gate_path"],
         tuple(item["channels"])) for item in scope
    }
    if declared != observed or len(scope) != len(declared):
        raise ValueError("geometry artifact linkage scope does not match geometry rows")


def canonical_manifest_bytes(manifest: dict[str, Any]) -> bytes:
    return (json.dumps(manifest, sort_keys=True, separators=(",", ":"),
                       ensure_ascii=False, allow_nan=False) + "\n").encode("utf-8")


def manifest_binding_object(manifest: dict[str, Any]) -> dict[str, Any]:
    """Return immutable provenance/method scope, excluding artifacts and time."""
    fields = (
        "manifest_schema", "export_operation_id", "profile", "export_profile",
        "identity_method", "workspace", "software", "transformation_context",
        "compensation_context", "approval", "acquisitions",
        "requested_populations",
    )
    return {field: manifest[field] for field in fields}


def manifest_binding_digest(manifest: dict[str, Any]) -> str:
    return hashlib.sha256(canonical_manifest_bytes(manifest_binding_object(manifest))).hexdigest()


def _atomic_write(path: Path, content: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".partial", dir=str(path.parent)
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        # Same-filesystem hard-link publication is atomic and fails if another
        # exporter has already claimed this immutable final path.
        os.link(temporary, path)
        temporary.unlink()
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def load_metadata(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError("contract metadata must be one JSON object")
    return value


def _require_text(mapping: dict[str, Any], fields: tuple[str, ...], context: str) -> None:
    missing = [field for field in fields
               if not isinstance(mapping.get(field), str) or not mapping[field].strip()]
    if missing:
        raise ValueError(f"missing required {context} field(s): {', '.join(missing)}")


def _reject_extra_keys(mapping: dict[str, Any], allowed: set[str], context: str) -> None:
    extras = sorted(set(mapping) - allowed)
    if extras:
        raise ValueError(f"unsupported {context} field(s): {', '.join(extras)}")


def _artifact_path(operation_dir: Path, relative: str) -> Path:
    if not isinstance(relative, str) or not relative:
        raise ValueError("artifact path must be a nonempty relative string")
    root = operation_dir.resolve()
    path = (root / relative).resolve()
    if root != path and root not in path.parents:
        raise ValueError("artifact path escapes the export operation directory")
    return path


def validate_manifest_structure(manifest: dict[str, Any]) -> None:
    required = {
        "manifest_schema", "export_operation_id", "created_at", "status", "profile",
        "export_profile", "identity_method", "workspace", "software", "approval",
        "acquisitions", "requested_populations", "populations", "artifacts", "geometry_overlay_status",
        "manifest_binding", "manifest_digest",
    }
    missing = sorted(required - set(manifest))
    if missing:
        raise ValueError(f"manifest missing required field(s): {', '.join(missing)}")
    if manifest["manifest_schema"] != {"name": SCHEMA_NAME, "version": SCHEMA_VERSION}:
        raise ValueError("manifest schema mismatch")
    if manifest["profile"] not in (PRODUCTION_PROFILE, LEGACY_PROFILE):
        raise ValueError("manifest profile mismatch")
    binding = manifest["manifest_binding"]
    if binding.get("digest") != manifest_binding_digest(manifest):
        raise ValueError("manifest binding digest mismatch")
    if not isinstance(manifest["artifacts"], list) or not isinstance(manifest["populations"], list):
        raise ValueError("manifest artifacts and populations must be lists")
    artifact_required = {
        "role", "path", "sha256", "byte_size", "row_count", "columns",
        "identity_columns", "identity_method_id", "intentionally_empty",
        "export_operation_id", "acquisition_id", "sample_id", "population_key",
        "gate_path", "channels",
    }
    for artifact in manifest["artifacts"]:
        if not isinstance(artifact, dict) or artifact_required - set(artifact):
            raise ValueError("artifact ledger entry schema mismatch")
    if manifest["profile"] == PRODUCTION_PROFILE:
        if manifest["identity_method"].get("verified") is not True:
            raise ValueError("production identity method is not verified")
        population_required = {
            "population_key", "gate_name", "gate_path", "gate_type",
            "parent_population_path", "acquisition_id", "sample_id", "prefix",
            "channels", "gate_channels", "row_count", "identity_field",
            "identity_method_id", "unique_identity_count", "duplicate_row_count",
            "intentionally_empty", "export_operation_id", "artifact_path",
            "artifact_sha256",
        }
        for population in manifest["populations"]:
            if not isinstance(population, dict) or population_required - set(population):
                raise ValueError("production population ledger entry schema mismatch")
        expected = {(item["sample_id"], key) for item in manifest["acquisitions"]
                    for key in manifest["requested_populations"]}
        observed = {(item["sample_id"], item["population_key"])
                    for item in manifest["populations"]}
        if not manifest["requested_populations"] or observed != expected:
            raise ValueError("production population ledger does not exactly cover requests/acquisitions")
        artifact_scope = {
            (item["sample_id"], item["population_key"])
            for item in manifest["artifacts"] if item["role"] == "population_events"
        }
        if artifact_scope != expected:
            raise ValueError("production population artifacts lack exact acquisition linkage")


def new_manifest(*, operation_id: str, profile: str, workspace: Path,
                 flowkit_version: str, metadata: dict[str, Any],
                 direct_index_semantics_verified: bool,
                 requested_populations: list[str] | None = None) -> dict[str, Any]:
    if profile not in (PRODUCTION_PROFILE, LEGACY_PROFILE):
        raise ValueError("unknown export profile")
    if not operation_id or any(character in operation_id for character in "/\\"):
        raise ValueError("export_operation_id must be nonempty and path-free")
    if not workspace.is_file():
        raise ValueError(f"workspace does not exist: {workspace}")
    approval = metadata.get("approval", {})
    acquisitions = metadata.get("acquisitions", [])
    software = metadata.get("software", {})
    _reject_extra_keys(metadata, ALLOWED_TOP_LEVEL, "metadata")
    if not isinstance(approval, dict) or not isinstance(software, dict):
        raise ValueError("approval and software metadata must be objects")
    _reject_extra_keys(approval, ALLOWED_APPROVAL, "approval")
    _reject_extra_keys(software, ALLOWED_SOFTWARE, "software")
    if not isinstance(acquisitions, list):
        raise ValueError("acquisitions metadata must be a list")
    for acquisition in acquisitions:
        if not isinstance(acquisition, dict):
            raise ValueError("each acquisition must be an object")
        _reject_extra_keys(acquisition, ALLOWED_ACQUISITION, "acquisition")
        if "prefix" in acquisition:
            validate_artifact_prefix(acquisition["prefix"])
    if profile == PRODUCTION_PROFILE:
        _require_text(approval, REQUIRED_APPROVAL, "approval")
        if not isinstance(acquisitions, list) or not acquisitions:
            raise ValueError("production metadata requires explicit acquisitions")
        for acquisition in acquisitions:
            _require_text(acquisition, REQUIRED_ACQUISITION, "acquisition")
            digest = acquisition["source_fcs_sha256"]
            if len(digest) != 64 or any(character not in "0123456789abcdef" for character in digest):
                raise ValueError("source_fcs_sha256 must be lowercase hexadecimal SHA-256")
        _require_text(metadata, ("transformation_context", "compensation_context"),
                      "provenance")
        _require_text(software, ("supported_flowkit_version",), "software")
        if software["supported_flowkit_version"] != flowkit_version:
            raise ValueError(
                "installed FlowKit version does not equal explicit supported_flowkit_version"
            )
        sample_ids = [item["sample_id"] for item in acquisitions]
        acquisition_ids = [item["acquisition_id"] for item in acquisitions]
        if len(sample_ids) != len(set(sample_ids)) or len(acquisition_ids) != len(set(acquisition_ids)):
            raise ValueError("production acquisition and sample IDs must be unique")
        if not direct_index_semantics_verified:
            raise ValueError(
                "production direct identity is conditional: run the pinned FlowKit "
                "source-index verification and explicitly attest success"
            )
    elif direct_index_semantics_verified:
        raise ValueError("legacy profile cannot claim verified direct identity")
    manifest = {
        "manifest_schema": {"name": SCHEMA_NAME, "version": SCHEMA_VERSION},
        "export_operation_id": operation_id,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "status": "draft",
        "profile": profile,
        "legacy_warning": (
            None if profile == PRODUCTION_PROFILE else
            "LEGACY COUNT-ONLY EXPORT: identity, containment, and geometry linkage are unverified"
        ),
        "export_profile": {"id": EXPORT_PROFILE_ID, "version": EXPORT_PROFILE_VERSION},
        "identity_method": {
            "id": DIRECT_METHOD_ID if profile == PRODUCTION_PROFILE else None,
            "version": DIRECT_METHOD_VERSION if profile == PRODUCTION_PROFILE else None,
            "source": "flowkit_get_gate_events_index" if profile == PRODUCTION_PROFILE else None,
            "semantics": "conditional_on_pinned_environment_verification",
            "verified": bool(direct_index_semantics_verified),
        },
        "workspace": {
            "filename": workspace.name,
            "local_path": str(workspace.resolve()),
            "sha256": sha256_file(workspace),
        },
        "software": {
            "exporter": software.get("exporter", "facspseudocolor FlowJo exporters"),
            "exporter_version": software.get("exporter_version", "1.0.0"),
            "source_commit": software.get("source_commit", "unavailable"),
            "flowjo_version": software.get("flowjo_version", "unavailable"),
            "flowkit_version": flowkit_version,
            "supported_flowkit_version": software.get("supported_flowkit_version", "unavailable"),
            "python_version": sys.version.split()[0],
        },
        "transformation_context": metadata.get("transformation_context", "unavailable"),
        "compensation_context": metadata.get("compensation_context", "unavailable"),
        "approval": approval,
        "acquisitions": acquisitions,
        "requested_populations": requested_populations or [],
        "populations": [],
        "artifacts": [],
        "geometry_overlay_status": "not_requested",
        "manifest_binding": {
            "algorithm": "sha256",
            "canonical_object": "manifest_binding_object_v1",
            "excludes": ["created_at", "populations", "artifacts", "full_manifest_digest"],
        },
        "manifest_digest": {
            "algorithm": "sha256",
            "target": "canonical UTF-8 export-manifest.json bytes",
            "storage": "export-manifest.sha256 sidecar",
        },
    }
    manifest["manifest_binding"]["digest"] = manifest_binding_digest(manifest)
    return manifest


def artifact_record(path: Path, *, operation_dir: Path, role: str,
                    row_count: int | None, columns: list[str], linkage: dict[str, Any],
                    identity_columns: list[str] | None = None,
                    intentionally_empty: bool = False) -> dict[str, Any]:
    resolved_root = operation_dir.resolve()
    resolved = path.resolve()
    if resolved_root != resolved and resolved_root not in resolved.parents:
        raise ValueError("artifact path escapes the export operation directory")
    if not path.is_file():
        raise ValueError(f"artifact does not exist: {path}")
    return {
        "role": role, "path": str(resolved.relative_to(resolved_root)),
        "sha256": sha256_file(path), "byte_size": path.stat().st_size,
        "row_count": row_count, "columns": columns,
        "identity_columns": identity_columns or [],
        "identity_method_id": DIRECT_METHOD_ID if identity_columns else None,
        "intentionally_empty": intentionally_empty,
        "export_operation_id": linkage.get("export_operation_id"),
        "acquisition_id": linkage.get("acquisition_id"),
        "sample_id": linkage.get("sample_id"),
        "population_key": linkage.get("population_key"),
        "gate_path": linkage.get("gate_path"),
        "channels": linkage.get("channels", []),
    }


def finalize_manifest(manifest: dict[str, Any], operation_dir: Path) -> tuple[Path, Path]:
    validate_manifest_structure(manifest)
    if manifest.get("status") != "draft":
        raise ValueError("only a draft manifest can be finalized")
    manifest_path = operation_dir / "export-manifest.json"
    digest_path = operation_dir / "export-manifest.sha256"
    if manifest_path.exists() or digest_path.exists():
        raise FileExistsError("refusing to overwrite a completed export operation")
    operation_id = manifest.get("export_operation_id")
    for artifact in manifest.get("artifacts", []):
        if artifact.get("export_operation_id") != operation_id:
            raise ValueError("artifact operation-ID inconsistency")
        path = _artifact_path(operation_dir, artifact["path"])
        if not path.is_file() or sha256_file(path) != artifact.get("sha256"):
            raise ValueError(f"artifact hash mismatch: {artifact.get('path')}")
    manifest = dict(manifest)
    manifest["status"] = "complete"
    content = canonical_manifest_bytes(manifest)
    digest = hashlib.sha256(content).hexdigest()
    _atomic_write(manifest_path, content)
    _atomic_write(digest_path, f"{digest}  export-manifest.json\n".encode("ascii"))
    return manifest_path, digest_path


def verify_finalized_manifest(operation_dir: Path) -> dict[str, Any]:
    manifest_path = operation_dir / "export-manifest.json"
    digest_path = operation_dir / "export-manifest.sha256"
    content = manifest_path.read_bytes()
    expected = digest_path.read_text(encoding="ascii").split()[0]
    if hashlib.sha256(content).hexdigest() != expected:
        raise ValueError("manifest digest mismatch")
    manifest = json.loads(content)
    validate_manifest_structure(manifest)
    if manifest.get("status") != "complete":
        raise ValueError("manifest is not complete")
    for artifact in manifest.get("artifacts", []):
        path = _artifact_path(operation_dir, artifact["path"])
        if sha256_file(path) != artifact["sha256"]:
            raise ValueError(f"artifact hash mismatch: {artifact['path']}")
    return manifest


def verify_population_artifacts(operation_dir: Path,
                                artifact_paths: list[Path]) -> dict[str, Any]:
    """Verify finalized manifest, ledger metadata, and every consumed CSV row."""
    manifest = verify_finalized_manifest(operation_dir)
    binding = manifest["manifest_binding"]["digest"]
    root = operation_dir.resolve()
    ledger = {item["path"]: item for item in manifest["artifacts"]
              if item["role"] == "population_events"}
    for supplied in artifact_paths:
        path = supplied.resolve()
        if root not in path.parents:
            raise ValueError("consumed population artifact escapes operation directory")
        relative = str(path.relative_to(root))
        record = ledger.get(relative)
        if record is None:
            raise ValueError(f"consumed population artifact is absent from ledger: {relative}")
        if path.stat().st_size != record["byte_size"] or sha256_file(path) != record["sha256"]:
            raise ValueError(f"consumed population artifact size/hash mismatch: {relative}")
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames != record["columns"]:
                raise ValueError(f"consumed population artifact schema mismatch: {relative}")
            rows = list(reader)
        if len(rows) != record["row_count"]:
            raise ValueError(f"consumed population artifact row-count mismatch: {relative}")
        expected = {
            "export_operation_id": manifest["export_operation_id"],
            "export_manifest_digest": binding,
            "acquisition_id": record["acquisition_id"],
            "sample_id": record["sample_id"],
            "export_profile": manifest["profile"],
        }
        for row in rows:
            for field, value in expected.items():
                if row.get(field) != value:
                    raise ValueError(f"consumed population artifact row binding mismatch: {relative}")
    return manifest


def finalize_geometry_supplement(supplement: dict[str, Any], operation_dir: Path) -> tuple[Path, Path]:
    """Finalize an immutable geometry member of an existing manifest set."""
    manifest_path = operation_dir / "geometry-manifest.json"
    digest_path = operation_dir / "geometry-manifest.sha256"
    if manifest_path.exists() or digest_path.exists():
        raise FileExistsError("refusing to overwrite completed geometry linkage")
    base = verify_finalized_manifest(operation_dir)
    if supplement.get("export_operation_id") != base.get("export_operation_id"):
        raise ValueError("geometry operation ID does not match population operation")
    if supplement.get("workspace_sha256") != base["workspace"]["sha256"]:
        raise ValueError("geometry workspace hash does not match population operation")
    if not isinstance(supplement.get("geometry_linkage_scope"), list):
        raise ValueError("geometry supplement lacks structured linkage scope")
    for artifact in supplement.get("artifacts", []):
        path = _artifact_path(operation_dir, artifact["path"])
        if artifact.get("export_operation_id") != base["export_operation_id"]:
            raise ValueError("geometry artifact operation-ID inconsistency")
        if sha256_file(path) != artifact.get("sha256"):
            raise ValueError("geometry artifact hash mismatch")
        if artifact.get("linkage_scope") != supplement["geometry_linkage_scope"]:
            raise ValueError("geometry artifact and supplement linkage scopes differ")
    supplement = dict(supplement)
    supplement["status"] = "complete"
    content = canonical_manifest_bytes(supplement)
    digest = hashlib.sha256(content).hexdigest()
    _atomic_write(manifest_path, content)
    _atomic_write(digest_path, f"{digest}  geometry-manifest.json\n".encode("ascii"))
    return manifest_path, digest_path


def _main() -> int:
    parser = argparse.ArgumentParser(description="Verify a finalized local export operation")
    parser.add_argument("--verify-operation", type=Path, required=True)
    parser.add_argument("--artifacts", nargs="+", type=Path, required=True)
    args = parser.parse_args()
    verify_population_artifacts(args.verify_operation, args.artifacts)
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
