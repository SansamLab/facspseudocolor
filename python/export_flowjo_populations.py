#!/usr/bin/env python3
"""Development-tree compatibility loader for the installed runtime asset."""

from pathlib import Path

_asset = Path(__file__).resolve().parents[1] / "inst" / "python" / "export_flowjo_populations.py"
exec(compile(_asset.read_bytes(), str(_asset), "exec"), globals(), globals())
