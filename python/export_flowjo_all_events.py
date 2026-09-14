#!/usr/bin/env python3
"""Export ALL raw events (no gating) per acquisition to CSV.

This is the general, reusable counterpart to the population exporter
(`export_flowjo_populations.py`): rather than reading gated populations from
a FlowJo workspace, it reads every raw event straight from each FCS file, with
no gating applied. It exists for the EdU Standard v2 report's DNA-A/DNA-H
"Single Cells" backgating display (`flowjo_all_events_dir` /
`flowjo_dna_height_channel` report params; see docs/REPORTS.md).

Combined with `export_flowjo_populations.py --profile
legacy_count_only_unverified_v1` (the default; see
`tools/flowjo-orchestration.R`'s `prepare_flowjo_csvs_external()` for every
`plot_type` other than `ph3`), this covers the full CSV-preparation need for
FlowJo-mode EdU/POI reports without any experiment-specific scripting. Only
`plot_type: "ph3"` needs the separate, stricter production-identity contract
in `export_contract.py`.

`event_index` values are the same raw FCS row indices FlowKit assigns per
acquisition, matching the population exporter's `event_index` column for the
same FCS file -- so the report can identify which "all events" rows fall
inside a gated population.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import flowkit as fk


def safe_filename(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("_")


def find_channel_column(columns, channel: str):
    """Locate the (PnN, PnS) column tuple for an exact PnN channel name.

    FlowKit's `Sample.as_dataframe()` returns MultiIndex columns of
    (short_name, string_name); the PnS half is often blank or
    instrument-specific and should not be guessed at by the caller.
    """
    matches = [column for column in columns if column[0] == channel]
    if not matches:
        available = sorted({column[0] for column in columns})
        raise SystemExit(f"channel {channel!r} not found; available: {available}")
    if len(matches) > 1:
        raise SystemExit(f"channel {channel!r} is ambiguous: {matches}")
    return matches[0]


def parse_samples(raw_samples: list[str]) -> list[tuple[str, str]]:
    samples = []
    for item in raw_samples:
        prefix, sep, fcs_name = item.partition(":")
        if not sep or not prefix or not fcs_name:
            raise SystemExit(
                f"--sample must be PREFIX:FCS_FILENAME, got {item!r}"
            )
        samples.append((prefix, fcs_name))
    prefixes = [prefix for prefix, _ in samples]
    if len(prefixes) != len(set(prefixes)):
        raise SystemExit("--sample prefixes must be unique")
    return samples


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fcs-dir", type=Path, required=True,
                        help="directory containing the FCS files")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--sample", action="append", required=True, dest="samples",
        metavar="PREFIX:FCS_FILENAME",
        help="one per acquisition; may be repeated. PREFIX must match the "
             "config.yml sample prefix; the report reads "
             "<output-dir>/<PREFIX>_all_events.csv."
    )
    parser.add_argument("--dna-channel", required=True,
                        help="FCS PnN short name for DNA area, e.g. FL2-A")
    parser.add_argument("--dna-height-channel", required=True,
                        help="FCS PnN short name for DNA height, e.g. FL2-H")
    parser.add_argument(
        "--dna-column-name", default=None,
        help="output DNA-area column name; must equal the config's "
             "`dna_channel`. Defaults to 'raw__<dna-channel>', matching "
             "export_flowjo_populations.py's native column-naming "
             "convention (no renaming needed if the config uses it too)."
    )
    parser.add_argument(
        "--dna-height-column-name", default=None,
        help="output DNA-height column name; must equal the report's "
             "`flowjo_dna_height_channel`. Defaults to "
             "'raw__<dna-height-channel>'."
    )
    args = parser.parse_args()

    samples = parse_samples(args.samples)
    for prefix, _ in samples:
        if safe_filename(prefix) != prefix:
            raise SystemExit(f"prefix {prefix!r} is not a safe filename component")

    dna_column_name = args.dna_column_name or f"raw__{args.dna_channel}"
    dna_height_column_name = (
        args.dna_height_column_name or f"raw__{args.dna_height_channel}"
    )
    if dna_column_name == dna_height_column_name:
        raise SystemExit("DNA and DNA-height output column names must differ")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    output_paths = {}
    for prefix, fcs_name in samples:
        fcs_path = args.fcs_dir / fcs_name
        if not fcs_path.is_file():
            raise SystemExit(f"missing FCS file: {fcs_path}")
        out_path = args.output_dir / f"{prefix}_all_events.csv"
        if out_path.exists():
            raise SystemExit(f"refusing to overwrite existing file: {out_path}")
        output_paths[prefix] = out_path

    for prefix, fcs_name in samples:
        fcs_path = args.fcs_dir / fcs_name
        sample = fk.Sample(str(fcs_path))
        events = sample.as_dataframe(source="raw")
        dna_col = find_channel_column(events.columns, args.dna_channel)
        height_col = find_channel_column(events.columns, args.dna_height_channel)
        table = events[[dna_col, height_col]].copy()
        table.columns = [dna_column_name, dna_height_column_name]
        table.insert(0, "event_index", events.index.map(str))
        table.insert(0, "sample_id", fcs_name)
        if table["event_index"].duplicated().any():
            raise SystemExit(
                f"duplicate event_index in all-events export for {fcs_name!r}"
            )
        out_path = output_paths[prefix]
        table.to_csv(out_path, index=False)
        print(f"  {prefix:16s} all events  n={len(table):6d}  -> {out_path.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
