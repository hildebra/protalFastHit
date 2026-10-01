#!/usr/bin/env python3
"""compare.py [BENCH_DIR] - 0.7.1 against 0.7.2-dev (the four features and the depth knobs, run.sh) on the v0.7.1
benchmark's samples, and the pipelines' own estimates; and 0.7.3-dev (VERSION=0.7.3dev TAG=v073 run.sh: the record
ratio, without depth knobs) where its runs are there.

Scores every run of v071 (the v0.7.1 benchmark's), v072 and v072k (run.sh) with the benchmark's score.py (the same
truth, findable species and metrics), and writes BENCH_DIR/results072/runs.tsv and summary.md:
- the models: each pipeline's F1 of species held out and on its independent test set, at 0.5 and, for 0.7.2-dev's
  long-read models, at their depth knobs (from trained_model*.metrics.json);
- per read type and database: mean F1, precision, recall, false positives per sample, Bray-Curtis, archaeal
  recall; and the paired difference per sample (0.7.2-dev less 0.7.1; for pb and ont also the features alone,
  v072k, and the depth knobs alone, v072 less v072k), its mean with a 95% bootstrap interval over the samples.
"""
import collections
import glob
import json
import os
import random
import statistics
import sys

B = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/bench071")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "2026-10-01-v071-benchmark", "scripts"))
sys.argv = [sys.argv[0], B]
import score as S  # noqa: E402  (the benchmark's truth and metrics)

OUT = os.path.join(B, "results072")
VARIANTS = {"v071": "0.7.1", "v072": "0.7.2-dev", "v072k": "0.7.2-dev, --knob 0.5", "v073": "0.7.3-dev"}
# Paired comparisons: (variant, against, label, read types or None for all).
COMPARISONS = [("v072", "v071", "0.7.2-dev less 0.7.1", None), ("v072k", "v071", "features alone", ("pb", "ont")),
               ("v072", "v072k", "depth knobs alone", ("pb", "ont")),
               ("v073", "v072", "0.7.3-dev less 0.7.2-dev", ("pe", "se")),
               ("v073", "v072k", "0.7.3-dev less 0.7.2-dev at 0.5", ("pb", "ont"))]
READS = ["pe", "se", "pb", "ont"]
METRICS = ["F1", "precision", "recall", "FP", "bray_curtis"]


def bootstrap(diffs, n=2000, seed=1):
    rng = random.Random(seed)
    means = sorted(statistics.mean(rng.choices(diffs, k=len(diffs))) for _ in range(n))
    return statistics.mean(diffs), means[int(0.025 * n)], means[int(0.975 * n) - 1]


def model_rows():
    rows = []
    for label, folder in (("0.7.1", "V071"), ("0.7.2-dev", "V072dev"), ("0.7.3-dev", "V073dev")):
        for t in READS:
            path = os.path.join(B, folder, "trained_model" + ("" if t == "pe" else "_" + t) + ".metrics.json")
            if not os.path.exists(path):
                continue
            with open(path) as fh:
                m = json.load(fh)
            knobs = m.get("depth_knobs", {}).get("knobs", {})
            rows.append([t, label, S.fmt(m["evaluation"]["species"]["F1"], 4), S.fmt(m["test"]["this one"]["F1"], 4),
                         ",".join(f"{b}:{k:g}" for b, k in sorted(knobs.items())) or "-",
                         S.fmt(m["depth_knobs"]["F1_at_depth_knobs"], 4) if knobs else "-",
                         S.fmt(m["test_depth_knobs"]["F1"], 4) if knobs else "-"])
    return rows


def main():
    runs = sorted(os.path.basename(p)[:-len(".done")] for p in glob.glob(os.path.join(S.RUNS, "*.done")))
    runs = [r for r in runs if r.split(".")[0] in VARIANTS]
    species_rows, rows = [], []
    for r in runs:
        rows.append(S.score(r, species_rows))
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "runs.tsv"), "w") as fh:
        columns = list(rows[0])
        fh.write("\t".join(columns) + "\n")
        fh.writelines("\t".join(str(d[c]) for c in columns) + "\n" for d in rows)

    by = collections.defaultdict(dict)  # (db, reads, variant) -> {sample: row}
    for r in rows:
        by[(r["db"], r["reads"], r["variant"])][r["sample"]] = r
    lines = ["# 0.7.1 and 0.7.2-dev on the v0.7.1 benchmark (compare.py)", ""]
    lines += S.table("The models: F1 of species held out (cross-validation) and on the pipeline's test set",
                     ["reads", "version", "species held out, 0.5", "test set, 0.5", "depth knobs",
                      "species held out, depth knobs", "test set, depth knobs"], model_rows())
    for db in ("full", "missing"):
        body = []
        for reads in READS:
            for v in VARIANTS:
                got = by.get((db, reads, v))
                if not got:
                    continue
                rs = list(got.values())
                archaea = sum(r["archaea"] for r in rs)
                body.append([reads, VARIANTS[v], len(rs)] + [S.fmt(S.mean([r[m] for r in rs]), 2 if m == "FP" else 4)
                                                            for m in METRICS]
                            + [S.fmt(sum(r["archaea_TP"] for r in rs) / archaea, 3) if archaea else "-"])
        lines += S.table(f"Means per sample, {db} database", ["reads", "version", "samples", "F1", "precision", "recall",
                                                               "FP per sample", "Bray-Curtis", "archaea recall"], body)
        body = []
        for reads in READS:
            for a, b, label, types in COMPARISONS:
                if types and reads not in types:
                    continue
                x, y = by.get((db, reads, a)), by.get((db, reads, b))
                if not x or not y:
                    continue
                shared = sorted(set(x) & set(y))
                cells = []
                for m in ("F1", "FP", "bray_curtis"):
                    mean, low, high = bootstrap([x[s][m] - y[s][m] for s in shared])
                    cells.append(f"{mean:+.4f} ({low:+.4f}, {high:+.4f})" if m != "FP" else f"{mean:+.2f} ({low:+.2f}, {high:+.2f})")
                body.append([reads, label, len(shared)] + cells)
        lines += S.table(f"Paired differences per sample, {db} database (mean, 95% bootstrap interval)",
                         ["reads", "comparison", "samples", "F1", "FP per sample", "Bray-Curtis"], body)
    text = "\n".join(lines) + "\n"
    with open(os.path.join(OUT, "summary.md"), "w") as fh:
        fh.write(text)
    print(text)


if __name__ == "__main__":
    main()
