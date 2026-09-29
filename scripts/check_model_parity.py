#!/usr/bin/env python3
"""Check that protal and the trainer agree: re-profile saved training samples with a model and compare.

Two things must hold for a model to do in protal what its training report says:
- protal (cPMML) gives each taxon the probability the model file gives it in Python
  (model_pmml.PmmlForest, which random_forest_cmdline.py checks against scikit-learn on every
  training row);
- protal computes the features as it did when the training data was collected. Another protal
  version may compute them differently, and the model then sees other inputs than it was trained on.

The samples come from a collect_training_data.py folder: each design point keeps its SAMs, truth
files and training dumps. protal re-profiles the SAMs of a few points (--profile_only, one run) with
--model; the new dumps' probabilities are compared with the model file scored in Python, and their
features with the dumps written during collection.

usage: check_model_parity.py --db DB --model MODEL.xml --training TRAINING_DIR [--points 2] [-o OUT]
Exit code 1 if a probability or a feature differs.
"""

import argparse
import glob
import os
import subprocess
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from model_pmml import PmmlForest  # noqa: E402


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--db", required=True, help="protal database")
    p.add_argument("--model", required=True, help="PMML model to check")
    p.add_argument("--training", required=True, help="output folder of collect_training_data.py")
    p.add_argument("--points", type=int, default=2, help="design points to re-profile (default 2: the first and last)")
    p.add_argument("--protal", default="protal")
    p.add_argument("-o", "--out", help="folder for protal's outputs (default: TRAINING/parity)")
    p.add_argument("-t", "--threads", type=int, default=4)
    return p.parse_args(argv)


def samples_of(point):
    """(sample, SAM path, truth path, collection dump) of each sample of a design point."""
    meta = os.path.join(point, "sim", "protal.meta")
    output_dir, rows = None, []
    with open(meta) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if fields[0] == "#OUTPUT_DIR":
                output_dir = fields[1]
            elif fields[0] == "#SAMPLEID":
                header = [f.lstrip("#") for f in fields]
            elif not line.startswith("#") and line.strip():
                rows.append(dict(zip(header, fields)))
    for row in rows:
        sam = os.path.join(output_dir, "alignments", row["SAM"])
        dumps = glob.glob(os.path.join(output_dir, "**", row["SAMPLEID"] + ".profile.truth_annotated"), recursive=True)
        if os.path.isfile(sam) and dumps:
            yield row["SAMPLEID"], sam, row["PROFILE_TRUTH"], dumps[0]


def read_dump(path):
    return pd.read_csv(path, sep="\t", float_precision="round_trip")


def main(argv=None):
    opts = parse_args(argv)
    points = sorted(glob.glob(os.path.join(opts.training, "points", "*")))
    if not points:
        sys.exit(f"no design points in {opts.training}/points")
    chosen = points if len(points) <= opts.points else \
        [points[round(i * (len(points) - 1) / max(1, opts.points - 1))] for i in range(opts.points)]
    samples = [s for point in chosen for s in samples_of(point)]
    if not samples:
        sys.exit("no samples with a SAM and a training dump in " + ", ".join(chosen))
    out = opts.out or os.path.join(opts.training, "parity")
    os.makedirs(out, exist_ok=True)
    command = [opts.protal, "--db", opts.db, "--model", opts.model, "--profile_only", ",".join(s[1] for s in samples),
               "--profile_truth", ",".join(s[2] for s in samples), "-o", out, "-t", str(opts.threads),
               "--no_strains", "--no_qcmsa"]
    with open(os.path.join(out, "protal.log"), "w") as log:
        rc = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT).returncode
    if rc != 0:
        sys.exit(f"protal failed with exit code {rc}; see {out}/protal.log")

    model = PmmlForest(opts.model)
    problems, rows, feature_diff = [], 0, {}
    for sample, _, _, collected in samples:
        new_dump = os.path.join(out, sample + ".profile.truth_annotated")
        if not os.path.isfile(new_dump):
            problems.append(f"{sample}: protal wrote no {new_dump}")
            continue
        new, old = read_dump(new_dump), read_dump(collected)
        rows += len(new)
        python = model.predict(new)
        diff = np.abs(python - new["probability"].to_numpy())
        if diff.max() > 0:
            problems.append(f"{sample}: protal's probability differs from the model file's by up to {diff.max():.3g} "
                            f"({int((diff > 0).sum())} of {len(new)} taxa)")
        if list(new.columns) != list(old.columns):
            problems.append(f"{sample}: the dump has other columns than during collection (another protal version)")
            continue
        joined = new.merge(old, on="taxon", suffixes=("", "_collected"))
        if len(joined) != len(new) or len(new) != len(old):
            problems.append(f"{sample}: {len(new)} taxa now, {len(old)} during collection")
        for c in model.features:
            a, b = joined[c].to_numpy(dtype=float), joined[c + "_collected"].to_numpy(dtype=float)
            d = np.abs(a - b) / np.maximum(1e-300, np.maximum(np.abs(a), np.abs(b)))
            if (d > 0).any():
                feature_diff[c] = max(feature_diff.get(c, 0.0), float(d.max()))
    if feature_diff:
        problems.append("protal computes features differently than when the training data was collected (largest "
                        "relative difference): " + ", ".join(f"{c} {v:.3g}" for c, v in sorted(feature_diff.items())))
    lines = [f"model {opts.model}: {len(model.trees)} trees, {len(model.features)} features",
             f"re-profiled {len(samples)} samples of {', '.join(os.path.basename(p) for p in chosen)}: {rows} taxa"]
    lines += ["PROBLEM: " + p for p in problems] or [
        "protal's probabilities equal the model file's (and so scikit-learn's) for every taxon; "
        "the features equal those of the training data"]
    with open(os.path.join(out, "parity.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
