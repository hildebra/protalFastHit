#!/usr/bin/env python3
"""Check that protal and the trainer agree: re-profile saved training samples with a model and compare.

Two things must hold for a model to do in protal what its training report says:
- protal (cPMML) gives each taxon the probability the model file gives it in Python
  (model_pmml.load_model: PmmlForest or PmmlBoosted, which machine_learning_cmdline.py checks against scikit-learn on every
  training row);
- protal computes the features as it did when the training data was collected. Another protal
  version may compute them differently, and the model then sees other inputs than it was trained on.

The samples come from a collect_training_data.py folder: each design point keeps its SAMs, truth
files and training dumps. protal re-profiles the SAMs of a few points (--profile_only, one run) with
--model; the new dumps' probabilities are compared with the model file scored in Python, and their
features with the dumps written during collection.

--read_type checks the model of other reads (se, pb, ont: the collection's samples of that read type,
profiled with protal's option for that model, e.g. --model_se). The single-end samples are the paired-end points'
first reads and, since 2026-10-15, the design's Ultima units (se_ultima_r<reads>, drawn like long reads): the points
are chosen among both (choose), so that the se model is checked on both kinds of its reads.

--level C (since 2026-10-15) checks a model of another sparse-first gate (build_gtdb_database.py --model-gates:
model_<type>_c<C>.xml) on the samples the collection profiled at that gate (collect_training_data.py --model_levels:
<sample>__c<C>, in <profile folder>_c<C>/, beside the database gate's), against their dumps there. Any level the
collector names (8, or 8s, 8h, 8g, 1f: with another pick or --flex_indels) is taken the same way. protal needs no
--sparse_first for it: --profile_only aligns nothing, the gate only steers the seeding (whose outcome the SAM holds),
and a --model given scores every sample whatever its gate (Options::ModelGate, SampleModelDbFile).

A feature that differs in its last digits only (relative difference up to ROUNDING) is noted, not a
problem: a sum added up in another order (in chunks on more threads, say) can round so, and that
moves no taxon across a tree's split in practice. protal sums its features so that they come out
the same on any number of threads, so such a note points to a sum that does not yet.

usage: check_model_parity.py --db DB --model MODEL.xml --training TRAINING_DIR [--read_type pe] [--points 2]
                             [--level C]
Exit code 1 if a probability or a feature differs (beyond rounding, for a feature).
"""

import argparse
import glob
import os
import re
import subprocess
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from model_pmml import load_model  # noqa: E402

# protal's option for each read type's model (ReadType.h).
MODEL_OPTIONS = {"pe": "--model", "se": "--model_se", "pb": "--model_pb", "ont": "--model_ont"}
# The largest relative difference of a feature that is rounding (see above): a double's last digits are ~1e-16.
ROUNDING = 1e-12
# collect_training_data.py's names: the scenarios' points (SCENARIO_PREFIX), the drawn units of each read type (units_of,
# ultima_design_units: pb_b<bases>, ont_b<bases>, se_ultima_r<reads>), and a sample at another level (level_dir,
# unit_map_rows: <sample>__c<level> in <profile folder>_c<level>).
SCENARIO_PREFIX = "sc_"
DRAWN_PREFIXES = {"pe": (), "se": ("se_ultima_",), "pb": ("pb_",), "ont": ("ont_",)}
LEVEL_SUFFIX, LEVEL_FOLDER = "__c", "_c"
LEVEL = re.compile(r"\d+[fshg]?")  # collect_training_data.unit_level_ids: a gate, with a pick or --flex_indels or not
LEVEL_WORDS = {"f": ", --flex_indels", "s": ", stride pick", "h": ", hash pick", "g": ", grid pick"}


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--db", required=True, help="protal database")
    p.add_argument("--model", required=True, help="PMML model to check")
    p.add_argument("--training", required=True, help="output folder of collect_training_data.py")
    p.add_argument("--read_type", default="pe", choices=sorted(MODEL_OPTIONS),
                   help="the reads the model is for (default pe)")
    p.add_argument("--points", type=int, default=2, help="design points to re-profile (default 2), shared among the kinds "
                                                        "of the read type's points (choose): the first and last of one "
                                                        "kind, the first of each of two (se: Illumina and Ultima)")
    p.add_argument("--level", help="the samples as the collection profiled them at this sparse-first gate, for a model "
                                   "of that gate (collect_training_data.py --model_levels; default: at the database's)")
    p.add_argument("--protal", default="protal")
    p.add_argument("--weighted_site_shift", action="store_true",
                   help="protal's --weighted_site_shift, as the collection profiled its samples with it "
                        "(collect_training_data.py --weighted_site_shift; since 2026-10-15: without it, protal computed "
                        "the shifted features otherwise than the dumps, and the check failed)")
    p.add_argument("-o", "--out", help="folder for protal's outputs, rewritten on every run (default: TRAINING/parity, "
                                       "or parity_<read type>, and _c<level> with --level)")
    p.add_argument("-t", "--threads", type=int, default=4)
    opts = p.parse_args(argv)
    if opts.level is not None and not LEVEL.fullmatch(opts.level):
        p.error(f"--level {opts.level!r}: a sparse-first gate as the collector names its levels (8, or 8s, 8h, 8g, 1f)")
    return opts


def level_text(level):
    """How a level is said (collect_training_data.level_text): "sparse-first gate 8", "sparse-first gate 1,
    --flex_indels"."""
    suffix = level[-1:] if level[-1:] in LEVEL_WORDS else ""
    return f"sparse-first gate {level[:len(level) - len(suffix)]}{LEVEL_WORDS.get(suffix, '')}"


def default_out(training, read_type, level=None):
    """TRAINING/parity[_<read type>][_c<level>]: where protal's outputs go without -o."""
    return os.path.join(training, ("parity" if read_type == "pe" else "parity_" + read_type) +
                        ("" if level is None else f"{LEVEL_FOLDER}{level}"))


def drawn(point):
    """Whether a design point's samples are drawn reads (long reads, Ultima reads: its sim/samples.tsv), not a
    simulator's map of paired-end samples (sim/protal.meta)."""
    return os.path.isfile(os.path.join(point, "sim", "samples.tsv"))


def sample_entry(folder, sample, sam_name, truth, level=None):
    """[(sample, SAM path, truth path, collection dump)] of a sample whose SAM and training dump are there, else []:
    as collected at the database's gate (`folder`, the SAM `sam_name`), or at `level` (<folder>_c<level>, the sample
    <sample>__c<level> and its SAM with the same extension, as collect_training_data.unit_map_rows names them)."""
    if level is not None:
        ending = sam_name[sam_name.index(sample) + len(sample):] if sample in sam_name else ".sam.gz"
        folder, sample = f"{folder}{LEVEL_FOLDER}{level}", f"{sample}{LEVEL_SUFFIX}{level}"
        sam_name = sample + ending
    sam = os.path.join(folder, "alignments", sam_name)
    dumps = glob.glob(os.path.join(glob.escape(folder), "**", glob.escape(sample) + ".profile.truth_annotated"),
                      recursive=True)
    return [(sample, sam, truth, dumps[0])] if os.path.isfile(sam) and dumps else []


def samples_of(point, read_type="pe", level=None):
    """(sample, SAM path, truth path, collection dump) of each sample of a design point: its paired-end samples,
    or their first reads alone (se, in protal_se), or its drawn samples (pb, ont, and se's Ultima units:
    sim/samples.tsv); with `level`, those samples as profiled at that level (sample_entry)."""
    if drawn(point):
        with open(os.path.join(point, "sim", "samples.tsv")) as fh:
            next(fh)
            for sample, _, truth, _ in (line.rstrip("\n").split("\t") for line in fh if line.strip()):
                yield from sample_entry(os.path.join(point, "protal"), sample, sample + ".sam.gz", truth, level)
        return
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
        if read_type == "se":
            sample = row["SAMPLEID"] + "_se"
            yield from sample_entry(os.path.join(point, "protal_se"), sample, sample + ".sam.gz", row["PROFILE_TRUTH"],
                                    level)
        else:
            yield from sample_entry(output_dir, row["SAMPLEID"], row["SAM"], row["PROFILE_TRUTH"], level)


def points_with_samples(training, read_type="pe", level=None):
    """[(design point, its samples_of)] of the collection's design points that have samples of the read type (at
    `level`, if given), by name. A point need not have them: with single-end setups of their own
    (collect_training_data.py --se_read_setups, since 2026-10-14) a point holds paired-end samples only, single-end
    ones only, or both. The read type's points: the paired-end points (sim/protal.meta) for pe and se, and its drawn
    units (DRAWN_PREFIXES: pb_*, ont_*, and since 2026-10-15 se's Ultima units se_ultima_*). The scenarios' points
    (sc_*, collect_training_data.SCENARIO_PREFIX) hold deep samples, and the design's give the same check."""
    points = []
    for point in sorted(glob.glob(os.path.join(glob.escape(training), "points", "*"))):
        name = os.path.basename(point)
        if name.startswith(SCENARIO_PREFIX):
            continue
        if drawn(point):
            if name.startswith(DRAWN_PREFIXES[read_type]):
                points.append(point)
        elif read_type in ("pe", "se") and os.path.isfile(os.path.join(point, "sim", "protal.meta")):
            points.append(point)
    return [(p, samples) for p, samples in ((p, list(samples_of(p, read_type, level))) for p in points) if samples]


def choose(points, n):
    """n of the points [(point, samples)] (all if they are no more), in name order: the points of each kind (simulated
    by simulate_metagenomes, or drawn: se's Ultima units) a share of n, in turn while a kind has points left, so that
    with two kinds and n 2 each gives one (since 2026-10-15); within a kind spread evenly by name, the first and the
    last of it for 2, its first for 1 (as before 2026-10-15 over all the points)."""
    if len(points) <= n:
        return list(points)
    kinds = [kind for kind in ([p for p in points if not drawn(p[0])], [p for p in points if drawn(p[0])]) if kind]
    shares, left = [0] * len(kinds), n
    while left > 0:
        for i, kind in enumerate(kinds):
            if left > 0 and shares[i] < len(kind):
                shares[i] += 1
                left -= 1
    chosen = []
    for kind, k in zip(kinds, shares):
        chosen += [kind[round(i * (len(kind) - 1) / max(1, k - 1))] for i in range(k)]
    return sorted(chosen, key=lambda p: p[0])


def read_dump(path):
    return pd.read_csv(path, sep="\t", float_precision="round_trip")


def feature_differences(joined, features):
    """Of each feature whose values differ between a sample's new dump and its collected one (`joined`: the dumps
    merged by taxon, the collected columns named <feature>_collected), the largest relative difference."""
    found = {}
    for c in features:
        a, b = joined[c].to_numpy(dtype=float), joined[c + "_collected"].to_numpy(dtype=float)
        d = np.abs(a - b) / np.maximum(1e-300, np.maximum(np.abs(a), np.abs(b)))
        if (d > 0).any():
            found[c] = float(d.max())
    return found


def compare_sample(model, sample, new, old):
    """A sample's problems (protal's probabilities against the model file's; its dump against the one written during
    collection) and the largest relative difference of each feature that differs. `new`, `old`: the sample's dumps
    now and during collection; either may have no taxa (a shallow sample of a database of few genes)."""
    problems = []
    diff = np.abs(model.predict(new) - new["probability"].to_numpy())
    if (diff > 0).any():
        problems.append(f"{sample}: protal's probability differs from the model file's by up to {diff.max():.3g} "
                        f"({int((diff > 0).sum())} of {len(new)} taxa)")
    if list(new.columns) != list(old.columns):
        problems.append(f"{sample}: the dump has other columns than during collection (another protal version)")
        return problems, {}
    joined = new.merge(old, on="taxon", suffixes=("", "_collected"))
    if len(joined) != len(new) or len(new) != len(old):
        problems.append(f"{sample}: {len(new)} taxa now, {len(old)} during collection")
    return problems, feature_differences(joined, model.features)


def split_rounding(differences):
    """The features (name -> largest relative difference) that protal computes differently, and those that differ
    by rounding only (up to ROUNDING)."""
    return ({c: d for c, d in differences.items() if d > ROUNDING},
            {c: d for c, d in differences.items() if d <= ROUNDING})


def main(argv=None):
    opts = parse_args(argv)
    at = "" if opts.level is None else f" at {level_text(opts.level)}"
    points = points_with_samples(opts.training, opts.read_type, opts.level)
    if not points:
        sys.exit(f"no design point in {opts.training}/points with {opts.read_type} samples{at} that have a SAM and a "
                 "training dump")
    chosen = choose(points, opts.points)
    samples = [s for _, of_point in chosen for s in of_point]
    chosen = [p for p, _ in chosen]
    out = opts.out or default_out(opts.training, opts.read_type, opts.level)
    os.makedirs(out, exist_ok=True)
    # --force: a rerun (of the build, say) profiles into the folder its last run filled, which --profile_only
    # otherwise refuses; with --profile_only it aligns nothing, it only writes the profiles again.
    command = [opts.protal, "--db", opts.db, MODEL_OPTIONS[opts.read_type], opts.model,
               "--profile_only", ",".join(s[1] for s in samples),
               "--profile_truth", ",".join(s[2] for s in samples), "-o", out, "-t", str(opts.threads),
               "--no_strains", "--no_qcmsa", "--force"] + (["--weighted_site_shift"] if opts.weighted_site_shift else [])
    with open(os.path.join(out, "protal.log"), "w") as log:
        rc = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT).returncode
    if rc != 0:
        sys.exit(f"protal failed with exit code {rc}; see {out}/protal.log")

    model = load_model(opts.model)
    problems, rows, empty, feature_diff = [], 0, 0, {}
    for sample, _, _, collected in samples:
        new_dump = os.path.join(out, sample + ".profile.truth_annotated")
        if not os.path.isfile(new_dump):
            problems.append(f"{sample}: protal wrote no {new_dump}")
            continue
        new, old = read_dump(new_dump), read_dump(collected)
        rows += len(new)
        empty += len(new) == 0
        found, differences = compare_sample(model, sample, new, old)
        problems += found
        for c, d in differences.items():
            feature_diff[c] = max(feature_diff.get(c, 0.0), d)
    if rows == 0 and not problems:  # nothing compared: no sample of these points has a taxon
        problems.append("no sample has a taxon to compare; check more points (--points)")
    differs, rounded = split_rounding(feature_diff)
    if differs:
        problems.append("protal computes features differently than when the training data was collected (largest "
                        "relative difference): " + ", ".join(f"{c} {v:.3g}" for c, v in sorted(differs.items())))
    lines = [f"model {opts.model} ({opts.read_type} reads): {len(model.trees)} trees, {len(model.features)} features",
             f"re-profiled {len(samples)} samples of {', '.join(os.path.basename(p) for p in chosen)}{at}: {rows} taxa"
             + (f" ({empty} sample{'s' if empty > 1 else ''} without any)" if empty else "")]
    lines += ["PROBLEM: " + p for p in problems] or [
        "protal's probabilities equal the model file's (and so scikit-learn's) for every taxon; "
        "the features equal those of the training data" + (" but for rounding" if rounded else "")]
    if rounded:
        lines.append("Note: features that differ from the training data's in their last digits only, as a sum added "
                     "up in another order does (largest relative difference): " +
                     ", ".join(f"{c} {v:.3g}" for c, v in sorted(rounded.items())))
    with open(os.path.join(out, "parity.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
