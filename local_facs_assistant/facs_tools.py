"""Fail-closed, read-only tools for local FACS experiment intake."""

from __future__ import annotations

import html
from html.parser import HTMLParser
import importlib
from pathlib import Path
import re
from typing import Any
from urllib.parse import unquote, unquote_plus, urlparse
import xml.etree.ElementTree as ET

MAX_XML_BYTES = 64 * 1024 * 1024
MAX_FCS_BYTES = 1024 * 1024 * 1024
MAX_INVENTORY_FILES = 2000
MAX_WORKSPACE_SAMPLES = 500
MAX_GATES_PER_SAMPLE = 1000
MAX_CHANNELS = 256
MAX_LAYOUTS = 200
MAX_NOTES_PER_LAYOUT = 100
MAX_NOTE_CHARS = 20000
MAX_TOTAL_NOTE_CHARS = 500000
ALLOWED_INVENTORY_SUFFIXES = {".fcs", ".wsp"}


class IntakeError(RuntimeError):
    """A safe intake operation could not be completed exactly."""


class _TextExtractor(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.parts: list[str] = []

    def handle_data(self, data: str) -> None:
        self.parts.append(data)


class ReadOnlyTools:
    """Read-only inspection bounded to one configured experiment root."""

    def __init__(self, experiment_root: str | Path) -> None:
        root = Path(experiment_root).expanduser().resolve(strict=True)
        if not root.is_dir():
            raise IntakeError("experiment root must be an existing directory")
        self.root = root

    def _resolve(self, value: str, *, suffixes: set[str] | None = None) -> Path:
        candidate = Path(value).expanduser()
        if not candidate.is_absolute():
            candidate = self.root / candidate
        try:
            resolved = candidate.resolve(strict=True)
            resolved.relative_to(self.root)
        except (OSError, ValueError) as exc:
            raise IntakeError(f"path is absent or outside experiment root: {value!r}") from exc
        if not resolved.is_file():
            raise IntakeError(f"path is not a regular file: {value!r}")
        if suffixes is not None and resolved.suffix.casefold() not in suffixes:
            raise IntakeError(f"unsupported file type for this tool: {resolved.suffix}")
        return resolved

    def _relative(self, path: Path) -> str:
        return path.relative_to(self.root).as_posix()

    def inventory_experiment(self) -> dict[str, Any]:
        files: list[dict[str, Any]] = []
        for path in sorted(self.root.rglob("*")):
            if path.is_symlink() or not path.is_file():
                continue
            if path.suffix.casefold() not in ALLOWED_INVENTORY_SUFFIXES:
                continue
            stat = path.stat()
            files.append({"path": self._relative(path), "suffix": path.suffix.casefold(), "bytes": stat.st_size})
            if len(files) > MAX_INVENTORY_FILES:
                raise IntakeError("experiment inventory exceeds safe file-count limit")
        return {"experiment_root": self.root.name, "files": files}

    def inspect_fcs_metadata(self, path: str) -> dict[str, Any]:
        fcs = self._resolve(path, suffixes={".fcs"})
        if fcs.stat().st_size > MAX_FCS_BYTES:
            raise IntakeError("FCS exceeds safe metadata-reading size limit")
        try:
            flowkit = importlib.import_module("flowkit")
        except ImportError as exc:
            raise IntakeError("FlowKit is required to inspect FCS metadata") from exc
        try:
            sample = flowkit.Sample(str(fcs), filename_as_id=True)
            metadata = dict(sample.get_metadata())
            pnn = list(sample.pnn_labels)
            pns = list(sample.pns_labels)
            if len(pnn) != len(pns):
                raise IntakeError("FlowKit returned inconsistent channel labels")
            if len(pnn) > MAX_CHANNELS:
                raise IntakeError("FCS exceeds safe channel-count limit")
            channels = [{"index": i + 1, "detector": str(detector), "stain": str(stain or "")} for i, (detector, stain) in enumerate(zip(pnn, pns))]
            event_count = int(sample.event_count)
        except IntakeError:
            raise
        except Exception as exc:
            raise IntakeError(f"FlowKit could not parse {self._relative(fcs)}: {exc}") from exc
        safe_keys = ("date", "btim", "etim", "cyt", "cytsn", "inst", "proj", "exp", "op")
        selected = {key: str(metadata[key]) for key in safe_keys if key in metadata}
        return {"file": self._relative(fcs), "event_count": event_count, "channels": channels, "metadata": selected}

    def inspect_wsp(self, path: str) -> dict[str, Any]:
        wsp = self._resolve(path, suffixes={".wsp"})
        tree = self._parse_xml(wsp)
        samples: list[dict[str, Any]] = []
        for sample in tree.findall(".//SampleList/Sample"):
            data_set = sample.find("./DataSet")
            node = sample.find("./SampleNode")
            if node is None:
                continue
            uri = data_set.get("uri", "") if data_set is not None else ""
            file_name = Path(unquote(urlparse(uri).path)).name if uri else node.get("name", "")
            gates: list[dict[str, Any]] = []
            terminal_paths: dict[str, tuple[str, ...]] = {}

            def walk_populations(element: ET.Element, parent_path: tuple[str, ...]) -> None:
                for child in list(element):
                    tag = child.tag.rsplit("}", 1)[-1]
                    if tag == "Population":
                        name = child.get("name")
                        if not name:
                            raise IntakeError(f"workspace sample {file_name!r} contains an unnamed Population")
                        gate_path = parent_path + (name,)
                        if name in terminal_paths:
                            raise IntakeError(
                                f"workspace sample {file_name!r} contains duplicate/ambiguous "
                                f"terminal gate name {name!r} at paths "
                                f"{'/'.join(terminal_paths[name])!r} and {'/'.join(gate_path)!r}"
                            )
                        terminal_paths[name] = gate_path
                        gates.append({"name": name, "path": list(gate_path)})
                        if len(gates) > MAX_GATES_PER_SAMPLE:
                            raise IntakeError("workspace sample exceeds safe gate-count limit")
                        walk_populations(child, gate_path)
                    else:
                        walk_populations(child, parent_path)

            walk_populations(node, ())
            workspace_sample_id = sample.get("sampleID") or (data_set.get("sampleID") if data_set is not None else None)
            samples.append({
                "workspace_sample_id": workspace_sample_id,
                "file": file_name,
                "sample_name": node.get("name"),
                "gate_names": [gate["name"] for gate in gates],
                "gates": gates,
            })
            if len(samples) > MAX_WORKSPACE_SAMPLES:
                raise IntakeError("workspace exceeds safe sample-count limit")
        if not samples:
            raise IntakeError("workspace has no readable SampleList/Sample/SampleNode entries")
        return {"workspace": self._relative(wsp), "flowjo_version": tree.getroot().get("flowJoVersion"), "samples": samples}

    def extract_layout_text(self, path: str) -> dict[str, Any]:
        wsp = self._resolve(path, suffixes={".wsp"})
        tree = self._parse_xml(wsp)
        layouts: list[dict[str, Any]] = []
        total_chars = 0
        for layout in tree.findall(".//LayoutEditor//Layout"):
            notes: list[str] = []
            for content in layout.findall(".//Content"):
                decoded = self._decode_layout_content(content.text or "")
                if decoded and not self._is_dynamic_annotation_only(decoded) and decoded not in notes:
                    if len(decoded) > MAX_NOTE_CHARS:
                        raise IntakeError("layout annotation exceeds safe text-length limit")
                    notes.append(decoded)
                    total_chars += len(decoded)
                    if len(notes) > MAX_NOTES_PER_LAYOUT or total_chars > MAX_TOTAL_NOTE_CHARS:
                        raise IntakeError("workspace exceeds safe layout-text limit")
            if notes:
                layouts.append({"layout": layout.get("name", ""), "text": notes})
                if len(layouts) > MAX_LAYOUTS:
                    raise IntakeError("workspace exceeds safe layout-count limit")
        return {"workspace": self._relative(wsp), "layouts": layouts}

    def _parse_xml(self, path: Path) -> ET.ElementTree:
        if path.stat().st_size > MAX_XML_BYTES:
            raise IntakeError("workspace exceeds safe XML-reading size limit")
        xml_bytes = path.read_bytes()
        if b"<!DOCTYPE" in xml_bytes.upper() or b"<!ENTITY" in xml_bytes.upper():
            raise IntakeError("DTD and entity declarations are not accepted")
        try:
            return ET.parse(path)
        except ET.ParseError as exc:
            raise IntakeError(f"invalid XML in {self._relative(path)}: {exc}") from exc

    @staticmethod
    def _decode_layout_content(value: str) -> str:
        decoded = value
        for _ in range(4):
            newer = html.unescape(decoded)
            if newer == decoded:
                break
            decoded = newer
        decoded = re.sub(r"<input\b[^>]*>", " ", decoded, flags=re.IGNORECASE)
        decoded = unquote_plus(decoded)
        parser = _TextExtractor()
        parser.feed(decoded)
        text = html.unescape(" ".join(parser.parts))
        return re.sub(r"\s+", " ", text).strip()

    @staticmethod
    def _is_dynamic_annotation_only(value: str) -> bool:
        return not value or ("fj.anno." in value and len(value.split()) < 8)

TOOL_DEFINITIONS = [
    {"type": "function", "function": {"name": "inventory_experiment", "description": "List only supported scientific files under the authorized experiment directory.", "parameters": {"type": "object", "properties": {}, "additionalProperties": False}}},
    {"type": "function", "function": {"name": "inspect_fcs_metadata", "description": "Read one FCS file's metadata and channel names with FlowKit; never returns event values.", "parameters": {"type": "object", "properties": {"path": {"type": "string"}}, "required": ["path"], "additionalProperties": False}}},
    {"type": "function", "function": {"name": "inspect_wsp", "description": "Read FlowJo workspace sample references and gate names without analyzing events.", "parameters": {"type": "object", "properties": {"path": {"type": "string"}}, "required": ["path"], "additionalProperties": False}}},
    {"type": "function", "function": {"name": "extract_layout_text", "description": "Decode text annotations from FlowJo layouts.", "parameters": {"type": "object", "properties": {"path": {"type": "string"}}, "required": ["path"], "additionalProperties": False}}}
]


def dispatch(tools: ReadOnlyTools, name: str, arguments: dict[str, Any]) -> dict[str, Any]:
    allowed = {definition["function"]["name"] for definition in TOOL_DEFINITIONS}
    if name not in allowed:
        raise IntakeError(f"unknown or unauthorized tool: {name}")
    method = getattr(tools, name)
    return method(**arguments)
