#!/usr/bin/env python3
"""score.py [BENCH_DIR] - the v0.7.1 benchmark's scoring (docs/claude/2026-10-01-v071-benchmark/scripts/score.py:
detection, abundance, archaea, low-abundance taxa, false positives, time and memory) of this benchmark's runs
(profile.sh: BENCH_DIR/runs_v072, the four versions run again together). Writes BENCH_DIR/results_v072/runs.tsv,
species.tsv and summary.md; overall.md: per database, read type and version the means over all samples; and
paired.md: per sample, each version less the one before it and 0.7.2 less 0.7.1,
mean F1, FP per sample, Bray-Curtis and wall time with a 95% bootstrap interval over the samples.
"""
import collections
import os
import random
import statistics
import sys

B = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/bench071")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "2026-10-01-v071-benchmark", "scripts"))
sys.argv = [sys.argv[0], B]
import score as S  # noqa: E402

S.RUNS = os.path.join(B, "runs_v072")
S.OUT = os.path.join(B, "results_v072")
S.VERSIONS = {"v060": "0.6.0a", "v070": "0.7.0", "v071": "0.7.1", "v072": "0.7.2", "v072k": "0.7.2, --knob 0.5"}
PAIRS = [("v070", "v060"), ("v071", "v070"), ("v072", "v071"), ("v072k", "v071"), ("v072", "v072k")]


def bootstrap(diffs, n=2000, seed=1):
    rng = random.Random(seed)
    means = sorted(statistics.mean(rng.choices(diffs, k=len(diffs))) for _ in range(n))
    return statistics.mean(diffs), means[int(0.025 * n)], means[int(0.975 * n) - 1]


def paired():
    by = collections.defaultdict(dict)
    with open(os.path.join(S.OUT, "runs.tsv")) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(header, line.rstrip("\n").split("\t")))
            by[(r["db"], r["reads"], r["variant"])][r["sample"]] = r
    lines = ["# Paired differences per sample (score.py)", "",
             "Mean over the samples both versions ran, 95% bootstrap interval. wall: seconds per run.", ""]
    for db in ("full", "missing"):
        lines += [f"### {db} database", "", "| reads | comparison | samples | F1 | FP per sample | Bray-Curtis | wall (s) |",
                  "|---|---|---|---|---|---|---|"]
        for reads in ("pe", "se", "pb", "ont"):
            for a, b in PAIRS:
                x, y = by.get((db, reads, a)), by.get((db, reads, b))
                if not x or not y:
                    continue
                shared = sorted(set(x) & set(y))
                cells = []
                for key, digits in (("F1", 4), ("FP", 2), ("bray_curtis", 4), ("wall_s", 1)):
                    diffs = [float(x[s][key]) - float(y[s][key]) for s in shared if x[s].get(key) and y[s].get(key)]
                    if not diffs:
                        cells.append("-")
                        continue
                    mean, low, high = bootstrap(diffs)
                    cells.append(f"{mean:+.{digits}f} ({low:+.{digits}f}, {high:+.{digits}f})")
                lines.append(f"| {reads} | {S.VERSIONS[a]} less {S.VERSIONS[b]} | {len(shared)} | " + " | ".join(cells) + " |")
        lines.append("")
    text = "\n".join(lines) + "\n"
    with open(os.path.join(S.OUT, "paired.md"), "w") as fh:
        fh.write(text)
    print(text)


def overall():
    rows = collections.defaultdict(list)
    with open(os.path.join(S.OUT, "runs.tsv")) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(header, line.rstrip("\n").split("\t")))
            rows[(r["db"], r["reads"], r["variant"])].append(r)
    num = lambda r, k: float(r[k]) if r.get(k) not in (None, "", "nan") else float("nan")
    mean = lambda v: statistics.mean([x for x in v if x == x]) if [x for x in v if x == x] else float("nan")
    lines = ["# Means over all samples (score.py)", "",
             "Per run as score.py scores it; archaea recall over all present archaea. Times: /usr/bin/time -v.", ""]
    for db in ("full", "missing"):
        lines += [f"### {db} database", "", "| reads | version | samples | F1 | precision | recall | FP per sample | "
                  "Bray-Curtis | median abs log2 | archaea recall | wall (s) | CPU (s) | peak RSS (GB) |",
                  "|" + "---|" * 13]
        for reads in ("pe", "se", "pb", "ont"):
            for v in S.VERSIONS:
                rs = rows.get((db, reads, v))
                if not rs:
                    continue
                archaea = sum(int(r["archaea"]) for r in rs)
                cells = [mean([num(r, "F1") for r in rs]), mean([num(r, "precision") for r in rs]),
                         mean([num(r, "recall") for r in rs]), mean([num(r, "FP") for r in rs]),
                         mean([num(r, "bray_curtis") for r in rs]), mean([num(r, "median_abs_log2") for r in rs])]
                lines.append(f"| {reads} | {S.VERSIONS[v]} | {len(rs)} | " +
                             " | ".join(f"{c:.4f}" if i != 3 else f"{c:.2f}" for i, c in enumerate(cells)) +
                             f" | {sum(int(r['archaea_TP']) for r in rs) / archaea:.3f}" +
                             f" | {mean([num(r, 'wall_s') for r in rs]):.1f} | {mean([num(r, 'cpu_s') for r in rs]):.1f}"
                             f" | {mean([num(r, 'peak_rss_gb') for r in rs]):.2f} |")
        lines.append("")
    text = "\n".join(lines) + "\n"
    with open(os.path.join(S.OUT, "overall.md"), "w") as fh:
        fh.write(text)
    print(text)


if __name__ == "__main__":
    S.main()
    overall()
    paired()
