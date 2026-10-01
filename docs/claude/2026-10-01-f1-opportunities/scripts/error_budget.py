#!/usr/bin/env python3
"""error_budget.py [BENCH_DIR] - where the F1 of the v0.7.1 benchmark's runs is lost (docs/claude/2026-10-01-v071-benchmark).

For every run of a version (default 0.7.1, $VARIANTS), each species is put in one class:
  TP        present (findable: a database has it), called
  FN, no reads   present, but no taxon with reads in the run's profile.log (nothing to call)
  FN, seen  present, with reads, probability below the knob
  FP        called, not present; by origin: a relative (same genus) of a simulated species no database has or the
            missing-species database leaves out; a congener of a present species; neither
and, per read type and depth (summed over the point's samples):
  F1 at the knob (0.5), the F1-optimal threshold and its F1 (an oracle: chosen on the same samples), and the
  ceiling: F1 if every present species with reads were called and nothing else (recall limited by the reads).
The species seen but missed are broken down by their simulated read pairs and by the genome they were simulated
from (the representative, or another strain, by its distance from the reference). Writes results/error_budget.md.
"""
import collections
import csv
import glob
import math
import os
import re
import sys

B = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/bench071")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "2026-10-01-v071-benchmark", "scripts"))
sys.argv = [sys.argv[0], B]
import score as bench  # noqa: E402  (the benchmark's truth, lineages and held-out sets)

VARIANTS = os.environ.get("VARIANTS", "v071 v070 v060").split()
RUNS = os.path.join(B, "runs")
OUT = os.path.join(HERE, "..", "results")
THRESHOLDS = [i / 100 for i in range(5, 96)]

distance = {}  # accession -> distance of the genome from its species' representative
with open(os.path.join(B, "world", "full", "simulation", "divergence.tsv")) as fh:
    rows = list(csv.DictReader(fh, delimiter="\t"))
div = {r["accession"]: float(r["strain_divergence"]) for r in rows}
rep_of = {}
for r in rows:
    if r["accession"].startswith("GCF_"):
        rep_of[r["species"].removeprefix("s__")] = r["accession"]
for r in rows:
    sp = r["species"].removeprefix("s__")
    distance[r["accession"]] = 0.0 if r["accession"] == rep_of.get(sp) else div[r["accession"]] + div.get(rep_of.get(sp), 0)


def manifest(sample):
    """{species: (read pairs simulated, max distance of its genomes from the reference)} of a sample's community."""
    m = re.match(r"((?:pb|ont)_b\d+)_s_(\d+)$", sample)
    if m:
        with open(os.path.join(bench.point_dir(m.group(1)), "sim", "samples.tsv")) as fh:
            sample = next(r["community"] for r in csv.DictReader(fh, delimiter="\t") if r["sample"] == sample)
    point = sample.rsplit("_s_", 1)[0]
    out = collections.defaultdict(lambda: [0, 0.0])
    with open(os.path.join(bench.point_dir(point), "sim", "manifest.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["sample"] == sample:
                s = out[bench.species_of(r["species"])]
                s[0] += int(r["read_pairs"])
                s[1] = max(s[1], distance.get(r["genome"], 0.0))
    return out


def seen(run):
    """{species: probability} of every taxon with reads in the run (profile.log)."""
    paths = glob.glob(os.path.join(RUNS, run, "**", "*.profile.log"), recursive=True)
    out = {}
    with open(paths[0]) as fh:
        lines = fh.read().splitlines()
    # 0.6.0a writes the same columns without a header line.
    header = lines[0].split("\t")
    if "Lineage" in header:
        lin, prob, lines = header.index("Lineage"), header.index("Probability"), lines[1:]
    else:
        lin, prob = 3, 1
    for line in lines:
        f = line.split("\t")
        if len(f) > max(lin, prob) and f[lin].startswith("d__"):
            out[bench.species_of(f[lin])] = float(f[prob])
    return out


def f1(tp, fp, fn):
    return 2 * tp / max(1, 2 * tp + fp + fn)


def main():
    budget = collections.defaultdict(lambda: collections.Counter())
    sweep = collections.defaultdict(lambda: collections.defaultdict(lambda: [0, 0, 0]))
    missed = collections.defaultdict(lambda: collections.Counter())
    for done in sorted(glob.glob(os.path.join(RUNS, "*.done"))):
        run = os.path.basename(done)[:-5]
        variant, db, reads, sample = run.split(".", 3)
        if variant not in VARIANTS:
            continue
        simulated = bench.truth_of(sample)
        findable = set(simulated) - bench.unknown - (bench.heldout if db == "missing" else set())
        unfindable_genera = {bench.genus.get(s) for s in set(simulated) - findable}
        findable_genera = collections.Counter(bench.genus.get(s) for s in findable)
        probs = seen(run)
        called = {s for s, p in probs.items() if p >= 0.5}
        point = re.sub(r"_s_\d+$", "", sample)
        key = (variant, db, reads, point)
        c = budget[key]
        man = manifest(sample)
        for s in findable:
            if s in called:
                c["TP"] += 1
            elif s not in probs:
                c["FN, no reads"] += 1
            else:
                c["FN, seen"] += 1
                pairs, dist = man.get(s, (0, 0.0))
                pbin = "1" if pairs <= 1 else "2-3" if pairs <= 3 else "4-10" if pairs <= 10 else "11-100" if pairs <= 100 else ">100"
                dbin = "reference" if dist == 0 else "strain <2%" if dist < 0.02 else "strain 2-4%" if dist < 0.04 else "strain >=4%"
                missed[(variant, db, reads)][("pairs", pbin)] += 1
                missed[(variant, db, reads)][("genome", dbin)] += 1
                missed[(variant, db, reads)][("domain", bench.domain.get(s, "?"))] += 1
        for s in called - findable:
            g = bench.genus.get(s)
            origin = "FP, relative of a missing species" if g in unfindable_genera else \
                "FP, congener of a present species" if findable_genera.get(g) else "FP, other"
            c[origin] += 1
        for t in THRESHOLDS:
            calls = {s for s, p in probs.items() if p >= t}
            tp = len(calls & findable)
            sweep[key][t][0] += tp
            sweep[key][t][1] += len(calls - findable)
            sweep[key][t][2] += len(findable) - tp
    lines = ["# Error budget (error_budget.py)", "", "Summed over the samples of each point; FN and FP by class; F1 at the "
             "knob, at the F1-optimal threshold of these samples (an oracle), and the ceiling if every present species "
             "with reads were called and nothing else.", ""]
    classes = ["TP", "FN, no reads", "FN, seen", "FP, relative of a missing species", "FP, congener of a present species", "FP, other"]
    lines.append("| version | db | reads | point | " + " | ".join(classes) + " | F1 at 0.5 | best threshold: F1 | ceiling |")
    lines.append("|" + "---|" * (len(classes) + 7))
    for key in sorted(budget, key=lambda k: (k[0], k[1], ["pe", "se", "pb", "ont"].index(k[2]), k[3])):
        c = budget[key]
        fp = sum(v for k, v in c.items() if k.startswith("FP"))
        fn = c["FN, no reads"] + c["FN, seen"]
        best_t, best = max(((t, f1(*v)) for t, v in sweep[key].items()), key=lambda x: x[1])
        ceiling = f1(c["TP"] + c["FN, seen"], 0, c["FN, no reads"])
        lines.append(f"| {key[0]} | {key[1]} | {key[2]} | {key[3]} | " + " | ".join(str(c[k]) for k in classes) +
                     f" | {f1(c['TP'], fp, fn):.3f} | {best_t:.2f}: {best:.3f} | {ceiling:.3f} |")
    lines += ["", "Present species seen but not called, by what they were simulated from:", ""]
    for key in sorted(missed):
        m = missed[key]
        parts = []
        for kind in ("pairs", "genome", "domain"):
            parts.append(", ".join(f"{label} {n}" for (k, label), n in sorted(m.items()) if k == kind))
        lines.append(f"- {' '.join(key)}: read pairs simulated: {parts[0]}; genome: {parts[1]}; {parts[2]}")
    # Pooled over depths: the threshold sweep per read type, for one threshold per read type.
    lines += ["", "One threshold per read type (pooled over every depth), F1 at 0.5 and at the best:", ""]
    pooled = collections.defaultdict(lambda: collections.defaultdict(lambda: [0, 0, 0]))
    for key, by_t in sweep.items():
        for t, v in by_t.items():
            p = pooled[key[:3]][t]
            for i in range(3):
                p[i] += v[i]
    for key in sorted(pooled):
        best_t, best = max(((t, f1(*v)) for t, v in pooled[key].items()), key=lambda x: x[1])
        lines.append(f"- {' '.join(key)}: {f1(*pooled[key][0.5]):.4f} at 0.5, {best:.4f} at {best_t:.2f}")
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "error_budget.md"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
