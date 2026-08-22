#!/usr/bin/env python3
"""Confirm a validated proposal and generate a non-executable QMD scaffold."""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import os
from pathlib import Path
import re
import stat
import sys
import tempfile
from typing import Any

from jsonschema import Draft202012Validator

try:
    from .assistant import preflight_inspect, schema_for_run, validate_proposal, validate_singleton_workspace_fcs
    from .facs_tools import IntakeError, MAX_FCS_BYTES, MAX_XML_BYTES, ReadOnlyTools
except ImportError:
    from assistant import preflight_inspect, schema_for_run, validate_proposal, validate_singleton_workspace_fcs
    from facs_tools import IntakeError, MAX_FCS_BYTES, MAX_XML_BYTES, ReadOnlyTools


CONFIRMATION_STATEMENT = (
    "I confirm all sample mappings, channel mappings, replicate assignments, "
    "control relationships, and analysis selections listed above."
)
REVIEWED_NAME = "facs_reviewed_config.v1.json"
QMD_NAME = "facs_reviewed_analysis_scaffold.v1.qmd"
MAX_CONFIG_BYTES = 2 * 1024 * 1024
MAX_RELATIONSHIPS = 64
METHOD_AUTHORIZATION = {
    "thresholds_confirmed": False,
    "normalization_confirmed": False,
    "background_quantile_confirmed": False,
    "dna_alignment_confirmed": False,
    "flowjo_export_authorized": False,
    "render_authorized": False,
    "analysis_output_write_authorized": False,
    "overwrite_authorized": False,
    "rebuild_authorized": False,
}


def _schema(name: str) -> dict[str, Any]:
    path = Path(__file__).with_name("schemas") / name
    return json.loads(path.read_text(encoding="utf-8"))


def _root(path: Path) -> Path:
    if path.is_symlink():
        raise IntakeError("experiment root may not be a symlink")
    resolved = path.resolve(strict=True)
    if not resolved.is_dir():
        raise IntakeError("experiment root must be a directory")
    return resolved


def _contained_regular(root: Path, path: Path) -> Path:
    candidate = path if path.is_absolute() else root / path
    if candidate.is_symlink():
        raise IntakeError(f"refusing symlink input: {candidate.name}")
    resolved = candidate.resolve(strict=True)
    if resolved.parent != root and root not in resolved.parents:
        raise IntakeError("input path escapes the experiment root")
    if not resolved.is_file():
        raise IntakeError("input must be a regular file")
    if resolved.stat().st_size > MAX_CONFIG_BYTES:
        raise IntakeError("configuration exceeds safe size limit")
    return resolved


def _read_json(path: Path) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise IntakeError(f"invalid JSON configuration: {exc}") from exc


def _open_nofollow(path: Path) -> int:
    if not hasattr(os, "O_NOFOLLOW"):
        raise IntakeError("this platform cannot safely open no-follow inputs")
    try:
        return os.open(path, os.O_RDONLY | os.O_NOFOLLOW)
    except OSError as exc:
        raise IntakeError(f"cannot safely open input {path.name}: {exc}") from exc


def _fd_identity(value: os.stat_result) -> tuple[int, int, int, int, int]:
    return (value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns, value.st_ctime_ns)


def _read_bounded_once(path: Path, limit: int) -> bytes:
    """Read one regular file buffer from one no-follow descriptor."""
    fd = _open_nofollow(path)
    try:
        before = os.fstat(fd)
        if not stat.S_ISREG(before.st_mode) or before.st_size > limit:
            raise IntakeError(f"input is not a bounded regular file: {path.name}")
        chunks: list[bytes] = []; remaining = before.st_size
        while remaining:
            chunk = os.read(fd, min(1024 * 1024, remaining))
            if not chunk: raise IntakeError(f"input was truncated while reading: {path.name}")
            chunks.append(chunk); remaining -= len(chunk)
        if os.read(fd, 1): raise IntakeError(f"input grew while reading: {path.name}")
        after = os.fstat(fd)
        path_after = os.stat(path, follow_symlinks=False)
        if (_fd_identity(before) != _fd_identity(after)
                or (path_after.st_dev, path_after.st_ino) != (after.st_dev, after.st_ino)
                or stat.S_ISLNK(path_after.st_mode)):
            raise IntakeError(f"input changed or was replaced while reading: {path.name}")
        return b"".join(chunks)
    finally:
        os.close(fd)


def _decode_json_bytes(data: bytes) -> Any:
    try: return json.loads(data.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise IntakeError(f"invalid JSON configuration: {exc}") from exc


def _hash_regular_file(path: Path, limit: int) -> dict[str, Any]:
    fd = _open_nofollow(path)
    try:
        before=os.fstat(fd)
        if not stat.S_ISREG(before.st_mode) or before.st_size > limit:
            raise IntakeError(f"experimental input is not a bounded regular file: {path.name}")
        digest=hashlib.sha256(); total=0
        while True:
            chunk=os.read(fd,1024*1024)
            if not chunk: break
            total += len(chunk)
            if total > limit: raise IntakeError(f"experimental input exceeds safe size limit: {path.name}")
            digest.update(chunk)
        after=os.fstat(fd); path_after=os.stat(path,follow_symlinks=False)
        if (total != before.st_size or _fd_identity(before) != _fd_identity(after)
                or (path_after.st_dev,path_after.st_ino)!=(after.st_dev,after.st_ino)
                or stat.S_ISLNK(path_after.st_mode)):
            raise IntakeError(f"experimental input changed or was replaced while hashing: {path.name}")
        return {"size_bytes":total,"sha256":digest.hexdigest()}
    finally:
        os.close(fd)


def _hash_inputs(root: Path, paths: list[str]) -> dict[str, dict[str, Any]]:
    if not paths or len(paths) > 501 or len(paths) != len(set(paths)):
        raise IntakeError("input provenance requires a bounded unique path set")
    result: dict[str, dict[str, Any]] = {}
    for relative in paths:
        if not isinstance(relative, str) or not relative or Path(relative).is_absolute():
            raise IntakeError("input provenance paths must be nonempty relative strings")
        candidate = root / relative
        if candidate.is_symlink():
            raise IntakeError(f"refusing symlink experimental input: {relative}")
        resolved = candidate.resolve(strict=True)
        if root not in resolved.parents or not resolved.is_file():
            raise IntakeError(f"experimental input escapes authorized root: {relative}")
        file_stat = resolved.stat()
        limit = MAX_FCS_BYTES if resolved.suffix.casefold() == ".fcs" else MAX_XML_BYTES
        if file_stat.st_size > limit:
            raise IntakeError(f"experimental input exceeds safe size limit: {relative}")
        result[relative] = _hash_regular_file(resolved,limit)
    return result


def _validate(schema: dict[str, Any], value: Any, label: str) -> None:
    errors = sorted(Draft202012Validator(schema).iter_errors(value), key=lambda error: list(error.path))
    if errors:
        raise IntakeError(f"{label} fails schema: {errors[0].message}")


def _confirmation_digest(value: dict[str, Any]) -> str:
    """Bind the explicit confirmation to the exact durable reviewed content."""
    bound = copy.deepcopy(value)
    bound.get("confirmation", {}).pop("content_sha256", None)
    payload = json.dumps(bound, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _canonical_digest(value: Any) -> str:
    payload = json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _control_relationships(values: list[str], proposal: dict[str, Any]) -> list[dict[str, Any]]:
    if not values or len(values) > MAX_RELATIONSHIPS:
        raise IntakeError("confirmation requires a bounded explicit control relationship list")
    sample_by_file = {row["file"]: row for row in proposal["sample_mapping"]}
    analysis_by_name = {row["name"]: row for row in proposal["analyses"]}
    poi_names = {row["name"] for row in proposal["analyses"] if row["analysis_type"] == "poi_vs_dna"}
    eligible_samples={path for path,row in sample_by_file.items() if row["role"]!="matched_background_control"}
    required_pairs={(name,path) for name in poi_names for path in eligible_samples}
    rows: list[dict[str, Any]] = []; covered_pairs: set[tuple[str,str]] = set()
    for raw in values:
        try: claim = json.loads(raw)
        except (TypeError, json.JSONDecodeError) as exc: raise IntakeError(f"invalid --control-relationship JSON: {exc}") from exc
        required={"control_file","relationship","applies_to_samples","applies_to_analyses","applies_to_features"}
        if not isinstance(claim,dict) or set(claim)!=required:
            raise IntakeError("--control-relationship has missing or unexpected fields")
        if claim["relationship"] != "matched_background_control":
            raise IntakeError("unsupported control relationship")
        control=claim["control_file"]
        if control not in sample_by_file or sample_by_file[control]["role"] != "matched_background_control":
            raise IntakeError("control_file must be an explicitly mapped matched background control")
        for field in ("applies_to_samples","applies_to_analyses","applies_to_features"):
            value=claim[field]
            if not isinstance(value,list) or not value or len(value)!=len(set(value)) or not all(isinstance(item,str) and item for item in value):
                raise IntakeError(f"{field} must be a nonempty unique string array")
        target_samples=set(claim["applies_to_samples"])
        if not target_samples.issubset(eligible_samples):
            raise IntakeError("control relationship targets an unknown or background-control sample")
        names=set(claim["applies_to_analyses"])
        if not names.issubset(poi_names) or not names:
            raise IntakeError("control relationship cites an unknown or non-POI analysis")
        features={analysis_by_name[name]["target_feature"] for name in names}
        if set(claim["applies_to_features"]) != features:
            raise IntakeError("control relationship feature scope must exactly match its analyses")
        pairs={(name,path) for name in names for path in target_samples}
        overlap=pairs & covered_pairs
        if overlap: raise IntakeError(f"control relationship scope overlaps existing analysis/sample pairs: {sorted(overlap)}")
        covered_pairs.update(pairs)
        rows.append({**copy.deepcopy(claim),"provenance":"operator_supplied","confirmed":True})
    if covered_pairs != required_pairs:
        missing=sorted(required_pairs-covered_pairs); extra=sorted(covered_pairs-required_pairs)
        raise IntakeError(f"control relationships must cover each POI-analysis/eligible-sample pair exactly once; missing={missing}, extra={extra}")
    return rows


def confirm_proposal(root_path: Path, proposal_path: Path, statement: str,
                     proposal_sha256: str, control_relationships: list[str],
                     attestation_context: str = "local_cli_operator") -> dict[str, Any]:
    """Reinspect immutable inputs and convert one complete v2 proposal to reviewed form."""
    if statement != CONFIRMATION_STATEMENT:
        raise IntakeError("confirmation statement must match the required text exactly")
    root = _root(root_path)
    proposal_file = _contained_regular(root, proposal_path)
    proposal_bytes = _read_bounded_once(proposal_file,MAX_CONFIG_BYTES)
    actual_proposal_sha256 = hashlib.sha256(proposal_bytes).hexdigest()
    if proposal_sha256 != actual_proposal_sha256:
        raise IntakeError("proposal SHA-256 challenge does not match the exact proposal bytes")
    proposal = _decode_json_bytes(proposal_bytes)
    input_paths = list(proposal.get("inputs", {}).get("fcs_files", []))
    workspace = proposal.get("inputs", {}).get("workspace")
    if workspace is not None: input_paths.append(workspace)
    before = _hash_inputs(root, input_paths)
    tools = ReadOnlyTools(root)
    ledger, _ = preflight_inspect(tools)
    after = _hash_inputs(root, input_paths)
    if before != after:
        raise IntakeError("experimental input changed during deterministic inspection")
    validate_singleton_workspace_fcs(ledger)
    inventory = next(row["output"] for row in ledger if row["tool"] == "inventory_experiment")
    fcs_files = [row["path"] for row in inventory["files"] if row["suffix"] == ".fcs"]
    wsp_files = [row["path"] for row in inventory["files"] if row["suffix"] == ".wsp"]
    runtime = schema_for_run(
        _schema("experiment-config.schema.json"), root.name, wsp_files, fcs_files,
        proposal.get("channels", []), proposal.get("sample_mapping", []),
        proposal.get("analyses", []), proposal.get("recorded_details", []),
    )
    validate_proposal(proposal, schema=runtime, ledger=ledger, root_name=root.name)
    if not proposal["sample_mapping"] or not proposal["channels"] or not proposal["analyses"]:
        raise IntakeError("confirmation requires complete sample, channel, and analysis claims")
    if proposal["inputs"]["workspace"] is None:
        raise IntakeError("confirmation requires one reconciled FlowJo workspace")
    relationships = _control_relationships(control_relationships, proposal)
    reviewed = {
        "schema_version": "reviewed-1.0",
        "source_proposal_schema": proposal["schema_version"],
        "experiment": copy.deepcopy(proposal["experiment"]),
        "inputs": copy.deepcopy(proposal["inputs"]),
        "input_provenance": before,
        "sample_mapping": [{**copy.deepcopy(row), "confirmed": True} for row in proposal["sample_mapping"]],
        "channels": [{**copy.deepcopy(row), "confirmed": True} for row in proposal["channels"]],
        "analyses": [{**copy.deepcopy(row), "confirmed": True} for row in proposal["analyses"]],
        "control_relationships": relationships,
        "recorded_details": copy.deepcopy(proposal["recorded_details"]),
        "confirmation": {
            "provenance": "operator_attestation_unverified", "attestation_context": attestation_context,
            "statement": statement, "proposal_sha256": proposal_sha256,
            "sample_mapping_confirmed": True, "channel_mapping_confirmed": True,
            "replicate_assignments_confirmed": True, "control_relationships_confirmed": True,
            "analysis_selection_confirmed": True,
        },
        "method_authorization": copy.deepcopy(METHOD_AUTHORIZATION),
        "uncertainties": copy.deepcopy(proposal["uncertainties"]),
        "evidence": copy.deepcopy(proposal["evidence"]),
    }
    reviewed["confirmation"]["content_sha256"] = _confirmation_digest(reviewed)
    _validate(_schema("reviewed-config.schema.json"), reviewed, "reviewed configuration")
    return reviewed


def _json_bytes(value: Any) -> bytes:
    return (json.dumps(value, indent=2, ensure_ascii=False) + "\n").encode("utf-8")


def _exclusive_batch(root: Path, outputs: dict[str, bytes]) -> list[Path]:
    """Create a validated batch without overwriting; roll back this batch on failure."""
    targets: list[Path] = []
    for name in outputs:
        if Path(name).name != name or name in {".", ".."}:
            raise IntakeError("generated target names must be simple basenames")
        target = root / name
        if target.exists() or target.is_symlink():
            raise IntakeError(f"refusing existing generated target: {name}")
        targets.append(target)
    created: list[tuple[Path, int, int]] = []
    try:
        for target in targets:
            fd, temporary = tempfile.mkstemp(prefix=".facs-reviewed-", dir=root)
            try:
                with os.fdopen(fd, "wb") as stream:
                    stream.write(outputs[target.name]); stream.flush(); os.fsync(stream.fileno())
                os.link(temporary, target)
                stat = target.stat(follow_symlinks=False)
                created.append((target, stat.st_dev, stat.st_ino))
            finally:
                try:
                    os.unlink(temporary)
                except FileNotFoundError:
                    pass
    except Exception:
        for target, device, inode in created:
            try:
                stat = target.stat(follow_symlinks=False)
                if stat.st_dev == device and stat.st_ino == inode and not target.is_symlink():
                    target.unlink()
            except FileNotFoundError:
                pass
        raise
    return targets


def save_reviewed(root_path: Path, proposal_path: Path, statement: str, proposal_sha256: str,
                  control_relationships: list[str], attestation_context: str = "local_cli_operator") -> Path:
    root = _root(root_path)
    reviewed = confirm_proposal(root, proposal_path, statement, proposal_sha256,
                                control_relationships, attestation_context)
    return _exclusive_batch(root, {REVIEWED_NAME: _json_bytes(reviewed)})[0]


def _slug(value: str, index: int) -> str:
    slug = re.sub(r"[^a-z0-9]+", "_", value.casefold()).strip("_")[:48]
    return f"facs_analysis_{index:02d}_{slug or 'unnamed'}.identity.v1.json"


def validate_identity(value: Any, reviewed: dict[str, Any], reviewed_name: str,
                      template_sha256: str) -> None:
    """Reject a stale, mixed, or altered generated identity."""
    _validate(_schema("analysis-identity.schema.json"),value,"analysis identity")
    if value["reviewed_config"] != reviewed_name:
        raise IntakeError("analysis identity references a different reviewed configuration")
    if value["reviewed_content_sha256"] != reviewed["confirmation"]["content_sha256"]:
        raise IntakeError("analysis identity has a stale reviewed configuration digest")
    if value["input_provenance_sha256"] != _canonical_digest(reviewed["input_provenance"]):
        raise IntakeError("analysis identity has stale input provenance")
    if value["template_sha256"] != template_sha256:
        raise IntakeError("analysis identity was generated from a different trusted template")
    if value["analysis"] not in reviewed["analyses"]:
        raise IntakeError("analysis identity contains an unreviewed analysis")
    for field in ("channels","sample_mapping","control_relationships"):
        if value[field] != reviewed[field]:
            raise IntakeError(f"analysis identity has mixed or stale {field}")
    if value["method_authorization"] != METHOD_AUTHORIZATION:
        raise IntakeError("analysis identity method authorization must remain false")


def generate_scaffold(root_path: Path, reviewed_path: Path) -> list[Path]:
    """Generate identity-only configs and a fail-closed QMD from a trusted template."""
    root = _root(root_path)
    reviewed_file = _contained_regular(root, reviewed_path)
    reviewed = _read_json(reviewed_file)
    _validate(_schema("reviewed-config.schema.json"), reviewed, "reviewed configuration")
    if reviewed["confirmation"]["content_sha256"] != _confirmation_digest(reviewed):
        raise IntakeError("reviewed configuration changed after operator attestation")
    if reviewed["experiment"]["directory"] != root.name:
        raise IntakeError("reviewed experiment directory does not match authorized root")
    if reviewed["method_authorization"] != METHOD_AUTHORIZATION:
        raise IntakeError("milestone 3 requires every method/render/output authorization false")
    input_paths=list(reviewed["inputs"]["fcs_files"])+[reviewed["inputs"]["workspace"]]
    if _hash_inputs(root,input_paths) != reviewed["input_provenance"]:
        raise IntakeError("experimental input provenance no longer matches reviewed configuration")
    template_path = Path(__file__).with_name("templates") / "reviewed_analysis.qmd.tmpl"
    if template_path.is_symlink() or not template_path.is_file():
        raise IntakeError("trusted QMD template is missing or is a symlink")
    template_bytes = template_path.read_bytes()
    template_sha256 = hashlib.sha256(template_bytes).hexdigest()
    template = template_bytes.decode("utf-8")
    provenance_sha256 = _canonical_digest(reviewed["input_provenance"])
    outputs: dict[str, bytes] = {}
    analysis_names: list[str] = []
    for index, analysis in enumerate(reviewed["analyses"], 1):
        name = _slug(analysis["name"], index)
        if name in outputs:
            raise IntakeError("analysis names produce duplicate generated targets")
        identity = {
            "schema_version": "analysis-identity-1.0",
            "reviewed_config": reviewed_file.name,
            "reviewed_content_sha256": reviewed["confirmation"]["content_sha256"],
            "input_provenance_sha256": provenance_sha256,
            "template_sha256": template_sha256,
            "analysis": analysis,
            "channels": reviewed["channels"],
            "sample_mapping": reviewed["sample_mapping"],
            "control_relationships": reviewed["control_relationships"],
            "method_choices": {
                "threshold_specification": None, "normalization_method": None,
                "background_quantile": None, "dna_alignment_method": None,
            },
            "method_authorization": copy.deepcopy(METHOD_AUTHORIZATION),
        }
        validate_identity(identity,reviewed,reviewed_file.name,template_sha256)
        outputs[name] = _json_bytes(identity); analysis_names.append(name)
    replacements = {
        "__TITLE_JSON__": json.dumps("Reviewed FACS analysis scaffold", ensure_ascii=False),
        "__REVIEWED_CONFIG_JSON__": json.dumps(reviewed_file.name, ensure_ascii=False),
        "__ANALYSIS_CONFIGS_JSON__": json.dumps(analysis_names, ensure_ascii=False),
        "__REVIEWED_DIGEST_JSON__": json.dumps(reviewed["confirmation"]["content_sha256"]),
        "__INPUT_PROVENANCE_DIGEST_JSON__": json.dumps(provenance_sha256),
        "__TEMPLATE_DIGEST_JSON__": json.dumps(template_sha256),
    }
    qmd = template
    for marker, value in replacements.items():
        qmd = qmd.replace(marker, value)
    if "__" in qmd:
        raise IntakeError("trusted template contains an unresolved marker")
    outputs[QMD_NAME] = qmd.encode("utf-8")
    return _exclusive_batch(root, outputs)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="action", required=True)
    confirm = subparsers.add_parser("confirm", help="explicitly confirm a complete schema-v2 proposal")
    confirm.add_argument("experiment_root", type=Path)
    confirm.add_argument("--proposal", type=Path, required=True)
    confirm.add_argument("--proposal-sha256", required=True)
    confirm.add_argument("--statement", required=True)
    confirm.add_argument("--control-relationship", action="append", default=[], required=True, metavar="JSON")
    confirm.add_argument("--attestation-context", choices=["local_cli_operator","conversation_conveyed_operator_attestation"], default="local_cli_operator")
    generate = subparsers.add_parser("generate", help="generate an identity-only QMD scaffold")
    generate.add_argument("experiment_root", type=Path)
    generate.add_argument("--reviewed-config", type=Path, default=Path(REVIEWED_NAME))
    args = parser.parse_args()
    try:
        if args.action == "confirm":
            paths = [save_reviewed(args.experiment_root, args.proposal, args.statement,
                                   args.proposal_sha256, args.control_relationship,
                                   args.attestation_context)]
        else:
            paths = generate_scaffold(args.experiment_root, args.reviewed_config)
    except (IntakeError, OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr); return 2
    for path in paths:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
