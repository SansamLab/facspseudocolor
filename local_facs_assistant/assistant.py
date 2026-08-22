#!/usr/bin/env python3
"""Build a deterministic read-only FACS proposal with optional local advisory classification."""

from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
import sys
from typing import Any
from urllib.error import URLError
from urllib.parse import urlparse
from urllib.request import HTTPRedirectHandler, ProxyHandler, Request, build_opener

from jsonschema import Draft202012Validator

try:
    from .facs_tools import IntakeError, ReadOnlyTools, TOOL_DEFINITIONS, dispatch
except ImportError:  # Direct execution: python local_facs_assistant/assistant.py
    from facs_tools import IntakeError, ReadOnlyTools, TOOL_DEFINITIONS, dispatch

MAX_TOTAL_TOOL_CALLS = 40
MAX_TOOL_OUTPUT_BYTES = 1024 * 1024
MAX_TOTAL_TOOL_OUTPUT_BYTES = 8 * 1024 * 1024
MAX_OLLAMA_RESPONSE_BYTES = 4 * 1024 * 1024
MAX_CHANNEL_ROLES = 32
MAX_ROLE_STRING_CHARS = 256
MAX_ROLE_JSON_CHARS = 2048
MAX_SAMPLE_MAPS = 500
MAX_ANALYSES = 32
MAX_CLAIM_STRING_CHARS = 256
SAMPLE_ROLES = {"matched_background_control", "experimental_sample", "untreated_control",
                "vehicle_control", "no_antibody_control", "positive_control", "other_control"}
ANALYSIS_TYPES = {"poi_vs_dna", "edu_vs_dna", "feature_vs_dna", "dna_only", "other"}
MODEL_REVIEW_SCHEMA = {
    "type": "object",
    "additionalProperties": False,
    "required": ["status", "flags"],
    "properties": {
        "status": {"enum": ["no_additional_uncertainty", "review_required"]},
        "flags": {
            "type": "array", "maxItems": 7, "uniqueItems": True,
            "items": {"enum": ["workspace_uncertain", "no_fcs_inputs", "channels_unspecified",
                                "samples_unspecified", "analyses_unspecified",
                                "recorded_details_missing", "human_review_requested"]},
        },
    },
    "oneOf": [
        {"properties": {"status": {"const": "no_additional_uncertainty"},
                        "flags": {"maxItems": 0}}},
        {"properties": {"status": {"const": "review_required"},
                        "flags": {"minItems": 1}}},
    ],
}


class _NoRedirect(HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        return None


def loopback_api_url(value: str) -> str:
    parsed = urlparse(value)
    if parsed.scheme != "http" or parsed.hostname not in {"localhost", "127.0.0.1", "::1"}:
        raise IntakeError("Ollama API must use plain HTTP on a loopback host")
    if parsed.username or parsed.password or parsed.query or parsed.fragment:
        raise IntakeError("Ollama API URL may not include credentials, query, or fragment")
    base = value.rstrip("/")
    return base if base.endswith("/api/chat") else base + "/api/chat"


def post_json(url: str, payload: dict[str, Any], timeout: float) -> dict[str, Any]:
    request = Request(url, data=json.dumps(payload).encode("utf-8"), headers={"Content-Type": "application/json"}, method="POST")
    try:
        # Explicitly ignore HTTP(S)_PROXY and related environment variables.
        # The only authorized peer is the validated loopback Ollama endpoint.
        with build_opener(ProxyHandler({}), _NoRedirect()).open(request, timeout=timeout) as response:
            if response.status != 200:
                raise IntakeError(f"Ollama returned HTTP {response.status}")
            raw = response.read(MAX_OLLAMA_RESPONSE_BYTES + 1)
            if len(raw) > MAX_OLLAMA_RESPONSE_BYTES:
                raise IntakeError("Ollama response exceeds safe byte limit")
            result = json.loads(raw.decode("utf-8"))
    except (URLError, TimeoutError, json.JSONDecodeError) as exc:
        raise IntakeError(f"local Ollama request failed: {exc}") from exc
    if not isinstance(result, dict) or not isinstance(result.get("message"), dict):
        raise IntakeError("Ollama response lacks a message object")
    return result


def chat_payload(model: str, messages: list[dict[str, Any]],
                 schema: dict[str, Any]) -> dict[str, Any]:
    """Build a non-thinking, schema-constrained Ollama chat request."""
    return {
        "model": model,
        "messages": messages,
        "stream": False,
        "think": False,
        "format": schema,
        "options": {"temperature": 0.2},
    }


def parse_final_json(content: Any) -> dict[str, Any]:
    """Accept exactly one JSON object; never strip prose or Markdown fences."""
    if not isinstance(content, str):
        raise IntakeError("model final response content was not text")
    try:
        proposed = json.loads(content)
    except json.JSONDecodeError as exc:
        raise IntakeError("model final response was not a single JSON object") from exc
    if not isinstance(proposed, dict):
        raise IntakeError("model final response was not a JSON object")
    return proposed


def preflight_inspect(tools: ReadOnlyTools) -> tuple[list[dict[str, Any]], str]:
    """Deterministically inspect every inventoried FCS and WSP or fail."""
    ledger: list[dict[str, Any]] = []
    total_bytes = 0

    def record(name: str, arguments: dict[str, Any]) -> dict[str, Any]:
        nonlocal total_bytes
        if len(ledger) >= MAX_TOTAL_TOOL_CALLS:
            raise IntakeError("deterministic preflight exceeds total inspection-call limit")
        # Preflight errors are intentionally not converted to model-visible
        # error objects: one unreadable expected input aborts the run.
        output = dispatch(tools, name, arguments)
        encoded = json.dumps(output, ensure_ascii=False)
        size = len(encoded.encode("utf-8"))
        if size > MAX_TOOL_OUTPUT_BYTES:
            raise IntakeError("preflight inspection output exceeds per-call byte limit")
        total_bytes += size
        if total_bytes > MAX_TOTAL_TOOL_OUTPUT_BYTES:
            raise IntakeError("preflight inspection outputs exceed total byte limit")
        ledger.append({"tool": name, "arguments": arguments, "output": output})
        return output

    inventory = record("inventory_experiment", {})
    for row in inventory["files"]:
        path, suffix = row["path"], row["suffix"]
        if suffix == ".fcs":
            record("inspect_fcs_metadata", {"path": path})
        elif suffix == ".wsp":
            record("inspect_wsp", {"path": path})
            record("extract_layout_text", {"path": path})
    encoded_ledger = json.dumps(ledger, ensure_ascii=False, separators=(",", ":"))
    if len(encoded_ledger.encode("utf-8")) > MAX_TOTAL_TOOL_OUTPUT_BYTES:
        raise IntakeError("serialized preflight ledger exceeds total byte limit")
    return ledger, encoded_ledger


def selection_uncertainties(wsp_candidates: list[str]) -> list[str]:
    uncertainties: list[str] = []
    if len(wsp_candidates) != 1:
        label = "No workspace candidate" if not wsp_candidates else f"Multiple workspace candidates ({len(wsp_candidates)})"
        uncertainties.append(f"{label}; inputs.workspace is null.")
    return uncertainties


def parse_channel_role(value: str) -> dict[str, str]:
    """Parse one auditable JSON channel-role claim without inference."""
    if not isinstance(value, str) or len(value) > MAX_ROLE_JSON_CHARS:
        raise IntakeError("--channel-role must be a JSON string within the size limit")
    try:
        claim = json.loads(value)
    except json.JSONDecodeError as exc:
        raise IntakeError(f"invalid --channel-role JSON: {exc}") from exc
    required = {"detector", "label", "category", "feature"}
    if not isinstance(claim, dict) or set(claim) != required:
        raise IntakeError("--channel-role requires exactly detector, label, category, and feature")
    if not all(isinstance(claim[field], str) and claim[field].strip() for field in required):
        raise IntakeError("all --channel-role values must be non-empty strings")
    if any(len(claim[field]) > MAX_ROLE_STRING_CHARS for field in required):
        raise IntakeError("--channel-role string exceeds safe length limit")
    if claim["category"] not in {"POI", "EdU", "DNA", "other"}:
        raise IntakeError("channel category must be one of POI, EdU, DNA, or other")
    return claim


def reconcile_channel_roles(values: list[str],
                            support: dict[tuple[str, str], list[str]]) -> list[dict[str, Any]]:
    """Require each user claim to match one exact observed detector/stain pair."""
    if len(values) > MAX_CHANNEL_ROLES:
        raise IntakeError("too many --channel-role claims")
    rows: list[dict[str, Any]] = []
    seen: set[tuple[str, str]] = set()
    seen_features: set[str] = set()
    for value in values:
        claim = parse_channel_role(value)
        key = (claim["detector"], claim["label"])
        if key in seen:
            raise IntakeError(f"duplicate --channel-role mapping for {key[0]!r}/{key[1]!r}")
        if key not in support:
            raise IntakeError(
                f"--channel-role detector/label pair was not observed in FCS metadata: "
                f"{key[0]!r}/{key[1]!r}"
            )
        if claim["feature"] in seen_features:
            raise IntakeError(f"duplicate --channel-role feature: {claim['feature']!r}")
        seen.add(key); seen_features.add(claim["feature"])
        rows.append({
            "category": claim["category"],
            "feature": claim["feature"],
            "detector": claim["detector"],
            "label": claim["label"],
            "provenance": "user_supplied",
            "confirmed": False,
            "evidence": sorted(set(support[key])),
        })
    return rows


def _parse_claim(value: str, option: str, required: set[str]) -> dict[str, Any]:
    if not isinstance(value, str) or len(value) > MAX_ROLE_JSON_CHARS:
        raise IntakeError(f"{option} must be a JSON string within the size limit")
    try:
        claim = json.loads(value)
    except json.JSONDecodeError as exc:
        raise IntakeError(f"invalid {option} JSON: {exc}") from exc
    if not isinstance(claim, dict) or set(claim) != required:
        raise IntakeError(f"{option} requires exactly {', '.join(sorted(required))}")
    return claim


def reconcile_sample_maps(values: list[str], fcs_files: list[str]) -> list[dict[str, Any]]:
    """Validate explicit sample claims; any supplied set must cover all FCS files."""
    if len(values) > MAX_SAMPLE_MAPS:
        raise IntakeError("too many --sample-map claims")
    if not values:
        return []
    required = {"file", "condition", "time", "role", "biological_replicate"}
    rows: list[dict[str, Any]] = []
    seen: set[str] = set()
    fcs_set = set(fcs_files)
    for value in values:
        claim = _parse_claim(value, "--sample-map", required)
        for field in ("file", "condition", "time", "role"):
            if not isinstance(claim[field], str) or not claim[field].strip():
                raise IntakeError(f"--sample-map {field} must be a nonempty string")
            if len(claim[field]) > MAX_CLAIM_STRING_CHARS:
                raise IntakeError(f"--sample-map {field} exceeds safe length limit")
        replicate = claim["biological_replicate"]
        if isinstance(replicate, bool) or not isinstance(replicate, int) or replicate < 1:
            raise IntakeError("--sample-map biological_replicate must be a positive integer")
        if claim["role"] not in SAMPLE_ROLES:
            raise IntakeError(f"unsupported --sample-map role: {claim['role']!r}")
        if claim["file"] not in fcs_set:
            raise IntakeError(f"--sample-map file is not an exact inventoried FCS path: {claim['file']!r}")
        if claim["file"] in seen:
            raise IntakeError(f"duplicate --sample-map for {claim['file']!r}")
        seen.add(claim["file"])
        rows.append({**claim, "provenance": "user_supplied", "confirmed": False,
                     "evidence": [claim["file"]]})
    if seen != fcs_set:
        raise IntakeError(f"--sample-map claims must cover every FCS exactly once; missing={sorted(fcs_set-seen)}")
    return rows


def workspace_gate_support(ledger: list[dict[str, Any]]) -> dict[str, set[str]]:
    """Return exact experiment-relative FCS -> gate names for one reconciled WSP."""
    return {path: set(gates) for path, gates in workspace_gate_path_support(ledger).items()}


def workspace_gate_path_support(ledger: list[dict[str, Any]]) -> dict[str, dict[str, tuple[str, ...]]]:
    """Return exact FCS -> terminal name -> full hierarchy path."""
    inventory = next(entry["output"] for entry in ledger if entry["tool"] == "inventory_experiment")
    wsp_paths = [row["path"] for row in inventory["files"] if row["suffix"] == ".wsp"]
    if len(wsp_paths) != 1:
        raise IntakeError("per-FCS gate support requires exactly one workspace")
    validate_singleton_workspace_fcs(ledger)
    fcs_paths = [row["path"] for row in inventory["files"] if row["suffix"] == ".fcs"]
    by_name = {Path(path).name: path for path in fcs_paths}
    inspection = next(entry["output"] for entry in ledger
                      if entry["tool"] == "inspect_wsp" and entry["output"]["workspace"] == wsp_paths[0])
    support: dict[str, dict[str, tuple[str, ...]]] = {}
    for sample in inspection["samples"]:
        path = by_name[sample["file"]]
        if path in support:
            raise IntakeError(f"duplicate workspace sample reference for {path!r}")
        records = sample["gates"]
        gates = [record["name"] for record in records]
        if len(gates) != len(set(gates)) or gates != sample["gate_names"]:
            raise IntakeError(f"duplicate/ambiguous gate names for {path!r}")
        support[path] = {record["name"]: tuple(record["path"]) for record in records}
    if set(support) != set(fcs_paths):
        raise IntakeError("per-FCS gate support does not cover the exact FCS inventory")
    return support


def workspace_gate_records(ledger: list[dict[str, Any]]) -> dict[str, list[dict[str, Any]]]:
    """Return full hierarchical gate records per experiment-relative FCS."""
    support = workspace_gate_support(ledger)
    inventory = next(entry["output"] for entry in ledger if entry["tool"] == "inventory_experiment")
    wsp_path = next(row["path"] for row in inventory["files"] if row["suffix"] == ".wsp")
    fcs_paths = [row["path"] for row in inventory["files"] if row["suffix"] == ".fcs"]
    by_name = {Path(path).name: path for path in fcs_paths}
    inspection = next(entry["output"] for entry in ledger
                      if entry["tool"] == "inspect_wsp" and entry["output"]["workspace"] == wsp_path)
    records = {by_name[sample["file"]]: sample["gates"] for sample in inspection["samples"]}
    if set(records) != set(support):
        raise IntakeError("hierarchical gate records do not match per-FCS gate support")
    return records


def reconcile_analyses(values: list[str], channel_rows: list[dict[str, Any]],
                       gate_support: dict[str, dict[str, tuple[str, ...]]]) -> list[dict[str, Any]]:
    """Validate explicit analysis declarations without choosing methods or thresholds."""
    if len(values) > MAX_ANALYSES:
        raise IntakeError("too many --analysis claims")
    required = {"name", "analysis_type", "target_feature", "dna_feature", "population"}
    feature_rows = {row["feature"]: row for row in channel_rows}
    rows: list[dict[str, Any]] = []
    names: set[str] = set()
    signatures: set[tuple[Any, ...]] = set()
    for value in values:
        claim = _parse_claim(value, "--analysis", required)
        for field in ("name", "analysis_type", "target_feature"):
            if not isinstance(claim[field], str) or not claim[field].strip():
                raise IntakeError(f"--analysis {field} must be a nonempty string")
            if len(claim[field]) > MAX_CLAIM_STRING_CHARS:
                raise IntakeError(f"--analysis {field} exceeds safe length limit")
        for field in ("dna_feature", "population"):
            if claim[field] is not None and (not isinstance(claim[field], str) or not claim[field].strip()):
                raise IntakeError(f"--analysis {field} must be null or a nonempty string")
            if isinstance(claim[field], str) and len(claim[field]) > MAX_CLAIM_STRING_CHARS:
                raise IntakeError(f"--analysis {field} exceeds safe length limit")
        kind = claim["analysis_type"]
        if kind not in ANALYSIS_TYPES:
            raise IntakeError(f"unsupported analysis_type: {kind!r}")
        if claim["target_feature"] not in feature_rows:
            raise IntakeError("--analysis target_feature must exactly reference a declared channel feature")
        target_category = feature_rows[claim["target_feature"]]["category"]
        if kind == "poi_vs_dna" and target_category != "POI":
            raise IntakeError("poi_vs_dna target_feature must be a POI feature")
        if kind == "edu_vs_dna" and target_category != "EdU":
            raise IntakeError("edu_vs_dna target_feature must be an EdU feature")
        if kind == "feature_vs_dna" and target_category != "other":
            raise IntakeError("feature_vs_dna target_feature must use category other")
        if kind == "dna_only" and target_category != "DNA":
            raise IntakeError("dna_only target_feature must be a DNA feature")
        if kind in {"poi_vs_dna", "edu_vs_dna", "feature_vs_dna"}:
            dna = feature_rows.get(claim["dna_feature"])
            if dna is None or dna["category"] != "DNA":
                raise IntakeError(f"{kind} requires dna_feature referencing a declared DNA channel")
        elif claim["dna_feature"] is not None:
            dna = feature_rows.get(claim["dna_feature"])
            if dna is None or dna["category"] != "DNA":
                raise IntakeError("optional dna_feature must reference a declared DNA channel")
        if claim["population"] is not None:
            if not gate_support:
                raise IntakeError("--analysis population requires exact per-FCS support from one workspace")
            missing_population = sorted(path for path, gates in gate_support.items()
                                        if claim["population"] not in gates)
            if missing_population:
                raise IntakeError(
                    f"--analysis population {claim['population']!r} is absent from mapped FCS files: "
                    f"{missing_population}"
                )
            hierarchy_paths = {gates[claim["population"]] for gates in gate_support.values()}
            if len(hierarchy_paths) != 1:
                raise IntakeError(
                    f"--analysis population {claim['population']!r} resolves to different "
                    f"hierarchy paths across FCS files: "
                    f"{sorted('/'.join(path) for path in hierarchy_paths)}"
                )
        signature = (kind, claim["target_feature"], claim["dna_feature"], claim["population"])
        if claim["name"] in names or signature in signatures:
            raise IntakeError("duplicate or conflicting --analysis declaration")
        names.add(claim["name"]); signatures.add(signature)
        rows.append({**claim, "provenance": "user_supplied", "confirmed": False})
    return rows


def observed_channel_support(ledger: list[dict[str, Any]]) -> dict[tuple[str, str], list[str]]:
    support: dict[tuple[str, str], list[str]] = {}
    for entry in ledger:
        if entry["tool"] != "inspect_fcs_metadata":
            continue
        output = entry["output"]
        for channel in output["channels"]:
            key = (channel["detector"], channel["stain"])
            support.setdefault(key, []).append(output["file"])
    return {key: sorted(set(paths)) for key, paths in support.items()}


def validate_singleton_workspace_fcs(ledger: list[dict[str, Any]]) -> None:
    """Require a singleton WSP to reference every inventoried FCS exactly once."""
    inventory = next(entry["output"] for entry in ledger if entry["tool"] == "inventory_experiment")
    wsp_paths = [row["path"] for row in inventory["files"] if row["suffix"] == ".wsp"]
    if len(wsp_paths) != 1:
        return
    fcs_paths = [row["path"] for row in inventory["files"] if row["suffix"] == ".fcs"]
    fcs_names = [Path(path).name for path in fcs_paths]
    if len(fcs_names) != len(set(fcs_names)):
        raise IntakeError("inventoried FCS basenames are duplicate and cannot be reconciled to FlowJo")
    inspections = [entry["output"] for entry in ledger if entry["tool"] == "inspect_wsp" and entry["output"]["workspace"] == wsp_paths[0]]
    if len(inspections) != 1:
        raise IntakeError("singleton workspace lacks exactly one successful inspection")
    refs = [sample["file"] for sample in inspections[0]["samples"]]
    if any(not isinstance(ref, str) or not ref for ref in refs):
        raise IntakeError("singleton workspace contains an empty FCS reference")
    if len(refs) != len(set(refs)):
        raise IntakeError("singleton workspace contains duplicate FCS references")
    missing = sorted(set(fcs_names) - set(refs))
    extra = sorted(set(refs) - set(fcs_names))
    if missing or extra:
        raise IntakeError(f"singleton workspace/FCS identity mismatch; missing={missing}, extra={extra}")


def schema_for_run(base_schema: dict[str, Any], root_name: str,
                   wsp_candidates: list[str], fcs_files: list[str],
                   channel_rows: list[dict[str, Any]],
                   sample_rows: list[dict[str, Any]], analysis_rows: list[dict[str, Any]],
                   recorded_details: list[dict[str, str]]) -> dict[str, Any]:
    """Deep-copy and pin deterministic per-run values for structured output."""
    schema = copy.deepcopy(base_schema)
    experiment = schema["properties"]["experiment"]["properties"]
    experiment["directory"] = {"const": root_name}
    inputs = schema["properties"]["inputs"]["properties"]
    inputs["fcs_files"] = {"const": fcs_files}
    inputs["workspace"] = {"const": wsp_candidates[0] if len(wsp_candidates) == 1 else None}
    schema["properties"]["sample_mapping"] = {"const": sample_rows}
    schema["properties"]["channels"] = {"const": channel_rows}
    schema["properties"]["analyses"] = {"const": analysis_rows}
    schema["properties"]["recorded_details"] = {"const": recorded_details}
    inspected_paths = fcs_files + wsp_candidates
    top_evidence = schema["properties"]["evidence"]
    top_evidence["items"] = {"enum": inspected_paths} if inspected_paths else {"not": {}}
    if not inspected_paths:
        top_evidence["maxItems"] = 0
    return schema


def mandatory_review_flags(context: dict[str, int]) -> list[str]:
    """Compute all non-discretionary advisory flags from bounded counts."""
    flags: list[str] = []
    if context["workspace_candidate_count"] != 1:
        flags.append("workspace_uncertain")
    if context["fcs_count"] == 0:
        flags.append("no_fcs_inputs")
    if context["channel_claim_count"] == 0:
        flags.append("channels_unspecified")
    if context["sample_claim_count"] == 0:
        flags.append("samples_unspecified")
    if context["analysis_claim_count"] == 0:
        flags.append("analyses_unspecified")
    if context["recorded_detail_count"] == 0:
        flags.append("recorded_details_missing")
    return flags


def validate_model_review(value: Any, mandatory_flags: list[str] | None = None) -> dict[str, Any]:
    """Validate the model's small enum-only advisory response exactly."""
    errors = sorted(Draft202012Validator(MODEL_REVIEW_SCHEMA).iter_errors(value),
                    key=lambda error: list(error.path))
    if errors:
        raise IntakeError(f"model advisory review fails schema: {errors[0].message}")
    mandatory = set(mandatory_flags or [])
    supplied = set(value["flags"])
    if not mandatory.issubset(supplied):
        raise IntakeError(f"model advisory omitted mandatory flags: {sorted(mandatory-supplied)}")
    if not (supplied - mandatory).issubset({"human_review_requested"}):
        raise IntakeError("model advisory added a non-discretionary flag not required by host counts")
    expected_status = "no_additional_uncertainty" if not supplied else "review_required"
    if value["status"] != expected_status:
        raise IntakeError(f"model advisory status must be {expected_status!r} for supplied flags")
    return value


def assemble_proposal(root_name: str, fcs_files: list[str], wsp_candidates: list[str],
                      channel_rows: list[dict[str, Any]], sample_rows: list[dict[str, Any]],
                      analysis_rows: list[dict[str, Any]], recorded_details: list[dict[str, str]],
                      model_review: dict[str, Any]) -> dict[str, Any]:
    """Construct the complete proposal deterministically; model review is advisory only."""
    workspace = wsp_candidates[0] if len(wsp_candidates) == 1 else None
    evidence = list(fcs_files) + ([workspace] if workspace is not None else [])
    return {
        "schema_version": "2.0",
        "experiment": {"directory": root_name, "title": None, "biological_replicates": None},
        "inputs": {"fcs_files": fcs_files, "workspace": workspace},
        "sample_mapping": sample_rows,
        "channels": channel_rows,
        "analyses": analysis_rows,
        "recorded_details": recorded_details,
        "model_review": {
            "provenance": "model_advisory",
            "status": model_review["status"],
            "flags": model_review["flags"],
        },
        "authorization": {
            "sample_mapping_confirmed": False,
            "channel_mapping_confirmed": False,
            "analysis_selection_confirmed": False,
            "analysis_authorized": False,
        },
        "uncertainties": selection_uncertainties(wsp_candidates),
        "evidence": evidence,
    }


def run_agent(root: Path, model: str, api_url: str, timeout: float,
              channel_roles: list[str] | None = None, sample_maps: list[str] | None = None,
              analyses: list[str] | None = None) -> dict[str, Any]:
    schema_path = Path(__file__).with_name("schemas") / "experiment-config.schema.json"
    base_schema = json.loads(schema_path.read_text(encoding="utf-8"))
    tools = ReadOnlyTools(root)
    ledger, _ = preflight_inspect(tools)
    validate_singleton_workspace_fcs(ledger)
    inventory = next(entry["output"] for entry in ledger if entry["tool"] == "inventory_experiment")
    wsp_candidates = [row["path"] for row in inventory["files"] if row["suffix"] == ".wsp"]
    fcs_files = [row["path"] for row in inventory["files"] if row["suffix"] == ".fcs"]
    channel_support = observed_channel_support(ledger)
    channel_rows = reconcile_channel_roles(channel_roles or [], channel_support)
    sample_rows = reconcile_sample_maps(sample_maps or [], fcs_files)
    gate_support = workspace_gate_path_support(ledger) if analyses and len(wsp_candidates) == 1 else {}
    analysis_rows = reconcile_analyses(analyses or [], channel_rows, gate_support)
    recorded_details: list[dict[str, str]] = []
    for entry in ledger:
        if entry["tool"] != "extract_layout_text":
            continue
        output = entry["output"]
        for layout in output["layouts"]:
            for note in layout["text"]:
                recorded_details.append({
                    "workspace": output["workspace"],
                    "layout": layout["layout"],
                    "text": note,
                })
    schema = schema_for_run(
        base_schema, tools.root.name, wsp_candidates, fcs_files, channel_rows,
        sample_rows, analysis_rows,
        recorded_details,
    )
    review_context = {
        "workspace_candidate_count": len(wsp_candidates),
        "fcs_count": len(fcs_files),
        "channel_claim_count": len(channel_rows),
        "sample_claim_count": len(sample_rows),
        "analysis_claim_count": len(analysis_rows),
        "recorded_detail_count": len(recorded_details),
    }
    messages: list[dict[str, Any]] = [{
        "role": "user",
        "content": "Return only the tiny JSON advisory review requested by the supplied schema. The host already validated and will construct all scientific mappings, paths, provenance, and authorization fields. You cannot add or change biological claims, paths, methods, gates, thresholds, normalization, or statistics. Flags are advisory enums only and never authorize analysis. REVIEW_CONTEXT:\n" + json.dumps(review_context, separators=(",", ":")),
    }]
    result = post_json(api_url, chat_payload(model, messages, MODEL_REVIEW_SCHEMA), timeout)
    message = result["message"]
    if message.get("tool_calls"):
        raise IntakeError("model attempted an unexpected tool call after deterministic preflight")
    review = validate_model_review(
        parse_final_json(message.get("content", "")),
        mandatory_review_flags(review_context),
    )
    proposed = assemble_proposal(
        tools.root.name, fcs_files, wsp_candidates, channel_rows, sample_rows,
        analysis_rows, recorded_details, review,
    )
    validate_proposal(proposed, schema=schema, ledger=ledger, root_name=tools.root.name)
    return proposed


def validate_proposal(value: Any, *, schema: dict[str, Any] | None = None,
                      ledger: list[dict[str, Any]] | None = None,
                      root_name: str | None = None) -> None:
    """Minimal fail-closed validator; the bundled JSON Schema remains authoritative."""
    if schema is not None:
        errors = sorted(Draft202012Validator(schema).iter_errors(value), key=lambda error: list(error.path))
        if errors:
            raise IntakeError(f"proposal fails JSON Schema: {errors[0].message}")
    required = {"schema_version", "experiment", "inputs", "sample_mapping", "channels", "analyses", "recorded_details", "model_review", "authorization", "uncertainties", "evidence"}
    if not isinstance(value, dict) or set(value) != required:
        raise IntakeError("proposed configuration has missing or unexpected top-level fields")
    if value["schema_version"] != "2.0":
        raise IntakeError("unsupported proposed configuration schema version")
    expected_objects = {
        "experiment": {"directory", "title", "biological_replicates"},
        "inputs": {"fcs_files", "workspace"},
        "authorization": {"sample_mapping_confirmed", "channel_mapping_confirmed", "analysis_selection_confirmed", "analysis_authorized"},
    }
    for field, keys in expected_objects.items():
        if not isinstance(value[field], dict) or set(value[field]) != keys:
            raise IntakeError(f"{field} has missing or unexpected fields")
    if not isinstance(value["experiment"]["directory"], str) or not value["experiment"]["directory"]:
        raise IntakeError("experiment.directory must be a non-empty string")
    if not isinstance(value["inputs"]["fcs_files"], list) or not all(isinstance(item, str) for item in value["inputs"]["fcs_files"]):
        raise IntakeError("inputs.fcs_files must be an array of strings")
    auth = value.get("authorization")
    if any(auth.get(field) is not False for field in ("analysis_authorized", "analysis_selection_confirmed", "sample_mapping_confirmed", "channel_mapping_confirmed")):
        raise IntakeError("authorization flags must all be false")
    if not all(isinstance(value.get(field), list) for field in ("sample_mapping", "channels", "analyses")):
        raise IntakeError("sample_mapping, channels, and analyses must be arrays")
    if not isinstance(value.get("recorded_details"), list):
        raise IntakeError("recorded_details must be an array")
    review = value.get("model_review")
    if (not isinstance(review, dict)
            or set(review) != {"provenance", "status", "flags"}
            or review.get("provenance") not in {"model_advisory", "host_advisory"}):
        raise IntakeError("model_review must be a structurally complete advisory object")
    for row in value["sample_mapping"]:
        expected = {"file", "condition", "time", "role", "biological_replicate", "provenance", "confirmed", "evidence"}
        if (not isinstance(row, dict) or set(row) != expected
                or row.get("confirmed") is not False
                or not isinstance(row.get("evidence"), list)):
            raise IntakeError("each sample mapping must contain a boolean confirmed field")
    for row in value["channels"]:
        expected = {"category", "feature", "detector", "label", "provenance", "confirmed", "evidence"}
        if (not isinstance(row, dict) or set(row) != expected
                or row.get("confirmed") is not False
                or not isinstance(row.get("evidence"), list)):
            raise IntakeError("each channel mapping must contain a boolean confirmed field")
    for row in value["analyses"]:
        expected = {"name", "analysis_type", "target_feature", "dna_feature", "population", "provenance", "confirmed"}
        if not isinstance(row, dict) or set(row) != expected or row.get("confirmed") is not False:
            raise IntakeError("each analysis declaration must be unconfirmed and structurally complete")
    if not isinstance(value["uncertainties"], list) or not all(isinstance(item, str) for item in value["uncertainties"]):
        raise IntakeError("uncertainties must be an array of strings")
    if not isinstance(value["evidence"], list) or not all(isinstance(item, str) for item in value["evidence"]):
        raise IntakeError("evidence must be an array of strings")
    if ledger is not None:
        reconcile_proposal(value, ledger, root_name)


def reconcile_proposal(value: dict[str, Any], ledger: list[dict[str, Any]], root_name: str | None) -> None:
    """Reject every proposal assertion not supported by deterministic tool output."""
    inventories = [entry["output"] for entry in ledger if entry["tool"] == "inventory_experiment" and "error" not in entry["output"]]
    if not inventories:
        raise IntakeError("proposal requires a successful inventory_experiment call")
    inventory = inventories[-1]
    inventoried = {row["path"] for row in inventory["files"]}
    fcs_set = {row["path"] for row in inventory["files"] if row["suffix"] == ".fcs"}
    if value["experiment"]["directory"] != root_name:
        raise IntakeError(
            f"experiment.directory must exactly equal authorized root name {root_name!r}; "
            f"received {value['experiment']['directory']!r}"
        )
    if value["experiment"]["title"] is not None:
        raise IntakeError("experiment.title is unsupported and must be null")
    if value["experiment"]["biological_replicates"] is not None:
        raise IntakeError("experiment.biological_replicates is unsupported and must be null")
    inputs = value["inputs"]
    if len(inputs["fcs_files"]) != len(set(inputs["fcs_files"])) or set(inputs["fcs_files"]) != fcs_set:
        raise IntakeError("proposed FCS files must exactly match the inventory")
    for field in ("workspace",):
        if inputs[field] is not None and inputs[field] not in inventoried:
            raise IntakeError(f"proposed {field} was not inventoried")
    wsp_candidates = sorted(row["path"] for row in inventory["files"] if row["suffix"] == ".wsp")
    expected_workspace = wsp_candidates[0] if len(wsp_candidates) == 1 else None
    if inputs["workspace"] != expected_workspace:
        raise IntakeError("workspace selection violates the deterministic single-candidate rule")
    required_uncertainties = selection_uncertainties(wsp_candidates)
    if not set(required_uncertainties).issubset(value["uncertainties"]):
        raise IntakeError("proposal omits required workspace candidate uncertainty")
    mapping_files = [row["file"] for row in value["sample_mapping"]]
    if mapping_files and (len(mapping_files) != len(set(mapping_files)) or set(mapping_files) != fcs_set):
        raise IntakeError("supplied sample mappings must uniquely cover every inventoried FCS file")

    fcs_inspections = {entry["output"]["file"]: entry["output"] for entry in ledger if entry["tool"] == "inspect_fcs_metadata" and "error" not in entry["output"]}
    channel_sources: dict[tuple[str, str], set[str]] = {}
    for path, inspection in fcs_inspections.items():
        for channel in inspection["channels"]:
            channel_sources.setdefault((channel["detector"], channel["stain"]), set()).add(path)
    for row in value["channels"]:
        detector, label = row["detector"], row["label"]
        if detector is None and label is None:
            continue
        if detector is None or label is None or (detector, label) not in channel_sources:
            raise IntakeError("proposed channel detector/label was not observed together in inspected FCS metadata")
        if not set(row["evidence"]).issubset(channel_sources[(detector, label)]) or not row["evidence"]:
            raise IntakeError("channel evidence must name inspected FCS files supporting that detector/label")

    wsp_inspections = {entry["output"]["workspace"] for entry in ledger if entry["tool"] == "inspect_wsp" and "error" not in entry["output"]}
    if inputs["workspace"] is not None and inputs["workspace"] not in wsp_inspections:
        raise IntakeError("proposed workspace must have a successful inspect_wsp call")
    expected_details: list[dict[str, str]] = []
    for entry in ledger:
        if entry["tool"] == "extract_layout_text" and "error" not in entry["output"]:
            for layout in entry["output"]["layouts"]:
                for note in layout["text"]:
                    expected_details.append({"workspace": entry["output"]["workspace"], "layout": layout["layout"], "text": note})
    if value["recorded_details"] != expected_details:
        raise IntakeError("recorded_details must exactly equal decoded FlowJo layout text")
    inspected_paths = set(fcs_inspections)
    inspected_paths |= {entry["output"]["workspace"] for entry in ledger if entry["tool"] in {"inspect_wsp", "extract_layout_text"} and "error" not in entry["output"]}
    for row in value["sample_mapping"]:
        if row["file"] not in row["evidence"] or row["file"] not in fcs_inspections:
            raise IntakeError("each sample row must cite its successfully inspected FCS file")
    all_evidence = list(value["evidence"])
    for row in value["sample_mapping"] + value["channels"]:
        all_evidence.extend(row["evidence"])
    if not set(all_evidence).issubset(inspected_paths):
        raise IntakeError("evidence may contain only exact inspected relative file paths")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("experiment_root", type=Path, help="single authorized experiment directory")
    parser.add_argument("--model", default="facs-assistant")
    parser.add_argument("--ollama-url", default="http://127.0.0.1:11434")
    parser.add_argument("--timeout", type=float, default=120.0)
    parser.add_argument(
        "--list-channels", action="store_true",
        help="print observed detector/label pairs and supporting FCS paths, then exit",
    )
    parser.add_argument(
        "--list-gates", action="store_true",
        help="print exact FlowJo gate names observed in the workspace, then exit",
    )
    parser.add_argument(
        "--sample-map", action="append", default=[], metavar="JSON",
        help='repeatable exact claim with file, condition, time, role, biological_replicate',
    )
    parser.add_argument(
        "--analysis", action="append", default=[], metavar="JSON",
        help='repeatable declaration with name, analysis_type, target_feature, dna_feature, population',
    )
    parser.add_argument(
        "--channel-role", action="append", default=[], metavar="JSON",
        help='repeatable exact claim, e.g. {"detector":"FL2-A","label":"PE-A","category":"DNA","feature":"DNA content"}',
    )
    args = parser.parse_args()
    try:
        if (args.list_channels or args.list_gates) and (args.channel_role or args.sample_map or args.analysis):
            raise IntakeError("discovery flags cannot be combined with channel, sample, or analysis claims")
        if args.list_channels and args.list_gates:
            raise IntakeError("choose only one of --list-channels or --list-gates")
        if args.list_channels:
            tools = ReadOnlyTools(args.experiment_root)
            ledger, _ = preflight_inspect(tools)
            rows = [
                {"detector": detector, "label": label, "fcs_files": paths}
                for (detector, label), paths in sorted(observed_channel_support(ledger).items())
            ]
            print(json.dumps(rows, indent=2, ensure_ascii=False))
            return 0
        if args.list_gates:
            tools = ReadOnlyTools(args.experiment_root)
            ledger, _ = preflight_inspect(tools)
            validate_singleton_workspace_fcs(ledger)
            records = workspace_gate_records(ledger)
            rows = [{"fcs_file": path, "gates": gates}
                    for path, gates in sorted(records.items())]
            print(json.dumps(rows, indent=2, ensure_ascii=False))
            return 0
        proposal = run_agent(
            args.experiment_root, args.model, loopback_api_url(args.ollama_url),
            args.timeout, args.channel_role, args.sample_map, args.analysis,
        )
    except (IntakeError, OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(proposal, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
