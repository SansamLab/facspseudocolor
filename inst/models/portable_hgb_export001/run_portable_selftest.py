#!/usr/bin/env python3
"""Standard-library replay of frozen synthetic reference predictions."""

import argparse
import hashlib
import importlib.util
import json
import math
import pathlib
import sys


def sha256(path):
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--package-dir", required=True)
    args = parser.parse_args()
    root = pathlib.Path(args.package_dir)
    manifest = json.loads((root / "MANIFEST.json").read_text(encoding="utf-8"))
    for item in manifest["files"]:
        if item["path"] != "MANIFEST.json" and sha256(root / item["path"]) != item["sha256"]:
            raise ValueError("package digest mismatch: " + item["path"])
    spec = importlib.util.spec_from_file_location("portable_hgb_predictor", root / "portable_hgb_predictor.py")
    predictor = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(predictor)
    fixtures = json.loads((root / "SYNTHETIC_EQUIVALENCE_FIXTURES.json").read_text(encoding="utf-8"))
    for fixture in fixtures["models"]:
        model_path = root / (fixture["model_id"] + ".portable.json")
        model = predictor.load_model(model_path, sha256(model_path))
        rows = [[float.fromhex(value) for value in row] for row in fixture["rows_float_hex"]]
        observed = predictor.predict_positive_proba(model, rows, fixture["feature_order"])
        expected = [float.fromhex(value) for value in fixture["expected_positive_probability_float_hex"]]
        threshold = float.fromhex(model["decision_threshold"])
        for index, (left, right) in enumerate(zip(observed, expected)):
            if not math.isclose(left, right, rel_tol=1e-12, abs_tol=1e-12):
                raise ValueError("self-test probability mismatch at row {}".format(index))
            if int(left >= threshold) != fixture["expected_calls"][index]:
                raise ValueError("self-test call mismatch at row {}".format(index))
    print("PASS: portable synthetic self-test")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print("STOP: {}".format(exc), file=sys.stderr)
        sys.exit(2)
