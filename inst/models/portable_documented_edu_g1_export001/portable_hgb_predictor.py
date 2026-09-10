#!/usr/bin/env python3
"""Safe, dependency-free predictor for PORTABLE_HGB_V1 JSON artifacts."""

import argparse
import csv
import hashlib
import json
import math
import pathlib
import sys


FORMAT = "PORTABLE_HGB_V1"


def _fail(message):
    raise ValueError(message)


def _float(value):
    if not isinstance(value, str):
        _fail("encoded float is not a string")
    result = float.fromhex(value)
    if not math.isfinite(result):
        _fail("model contains a non-finite float")
    return result


def _sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_model(path, expected_sha256=None):
    path = pathlib.Path(path)
    if expected_sha256 and _sha256(path) != expected_sha256:
        _fail("portable model digest mismatch")
    with path.open("r", encoding="utf-8") as handle:
        model = json.load(handle)
    if model.get("format") != FORMAT:
        _fail("unsupported portable model format")
    if model.get("task") != "binary_classification" or model.get("link") != "logit":
        _fail("only binary logit HGB models are supported")
    features = model.get("feature_order")
    if not isinstance(features, list) or not features or len(set(features)) != len(features):
        _fail("invalid feature order")
    if model.get("n_features_in") != len(features):
        _fail("feature count mismatch")
    if len(model.get("classes", [])) != 2:
        _fail("binary class order is missing")
    _float(model["baseline_raw"])
    _float(model["decision_threshold"])
    for iteration in model.get("iterations", []):
        if len(iteration) != 1:
            _fail("binary model must contain one tree per iteration")
        tree = iteration[0]
        nodes = tree.get("nodes", [])
        if not nodes:
            _fail("tree has no nodes")
        for index, node in enumerate(nodes):
            if node.get("index") != index:
                _fail("non-contiguous node indexes")
            _float(node["value"])
            _float(node["num_threshold"])
    return model


def _bitset_contains(bitsets, bitset_index, value):
    integer = int(value)
    if integer != value or integer < 0 or integer >= 256:
        return False
    words = bitsets[bitset_index]
    return bool((words[integer >> 5] >> (integer & 31)) & 1)


def _tree_value(tree, row, known_categories):
    nodes = tree["nodes"]
    index = 0
    while True:
        node = nodes[index]
        if node["is_leaf"]:
            return _float(node["value"])
        feature = node["feature_idx"]
        value = row[feature]
        if math.isnan(value):
            go_left = node["missing_go_to_left"]
        elif node["is_categorical"]:
            if value < 0:
                go_left = node["missing_go_to_left"]
            elif _bitset_contains(tree["raw_left_cat_bitsets"], node["bitset_idx"], value):
                go_left = True
            elif value in known_categories[feature]:
                go_left = False
            else:
                go_left = node["missing_go_to_left"]
        else:
            go_left = value <= _float(node["num_threshold"])
        index = node["left"] if go_left else node["right"]


def predict_positive_proba(model, rows, feature_order=None):
    if feature_order is not None and list(feature_order) != model["feature_order"]:
        _fail("input feature order does not exactly match the frozen order")
    n_features = model["n_features_in"]
    known = [set(_float(value) for value in values) for values in model["bin_mapper"]["known_categories"]]
    probabilities = []
    for row in rows:
        if len(row) != n_features:
            _fail("input row has the wrong feature count")
        converted = [float(value) for value in row]
        if any(not math.isfinite(value) for value in converted):
            _fail("input contains a non-finite value")
        raw = _float(model["baseline_raw"])
        for iteration in model["iterations"]:
            raw += _tree_value(iteration[0], converted, known)
        if raw >= 0.0:
            probability = 1.0 / (1.0 + math.exp(-raw))
        else:
            exp_raw = math.exp(raw)
            probability = exp_raw / (1.0 + exp_raw)
        probabilities.append(probability)
    return probabilities


def predict_calls(model, rows, feature_order=None):
    threshold = _float(model["decision_threshold"])
    return [int(value >= threshold) for value in predict_positive_proba(model, rows, feature_order)]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", required=True)
    parser.add_argument("--model-sha256", required=True)
    parser.add_argument("--input-csv", required=True)
    parser.add_argument("--output-csv", required=True)
    args = parser.parse_args()
    output = pathlib.Path(args.output_csv)
    if output.exists():
        _fail("output collision")
    model = load_model(args.model, args.model_sha256)
    with pathlib.Path(args.input_csv).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != model["feature_order"]:
            _fail("CSV feature order does not exactly match the frozen order")
        rows = [[row[name] for name in reader.fieldnames] for row in reader]
    probabilities = predict_positive_proba(model, rows, model["feature_order"])
    threshold = _float(model["decision_threshold"])
    with output.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["row_index", "positive_probability", "call"])
        for index, probability in enumerate(probabilities):
            writer.writerow([index, probability.hex(), int(probability >= threshold)])
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print("STOP: {}".format(exc), file=sys.stderr)
        sys.exit(2)
