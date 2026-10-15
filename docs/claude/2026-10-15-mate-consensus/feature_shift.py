#!/usr/bin/env python3
"""How the profiles' features and the model's calls move between two runs of the same samples (e.g. --mate_guidance
guide against consensus), from their .truth_annotated tables (every taxon with records: truth, the model's call and
probability, and its features).

  OUT.features.tsv  per feature (numeric column): rows in both runs, the share of them where it changed, its mean in
                    each run over present (truth 1) and absent (truth 0) taxa, and the mean change in units of the
                    feature's interquartile range over all rows of run A (its "scale").
  OUT.flips.tsv     the taxa whose call differs between the runs: sample, taxon, truth, the probability in each run
                    and the five features that moved most (by their scale).
  OUT.rows.tsv      taxa with a row in one run only.
  printed: the calls of both runs over the rows (TP/FP/FN), the flips by truth, the probability's mean change, and
  the features that moved most.

usage: feature_shift.py OUT_PREFIX LABEL_A LABEL_B [--features-from MODEL_REPORT] SAMPLE:TABLE_A:TABLE_B [...]
  --features-from: only the model's features, from its trainer report's "N features (...): a, b, ..." line.
Standard library only.
"""
import csv
import re
import statistics
import sys

SKIP = {"truth", "prediction", "probability", "taxon", "taxon_name"}


def model_features(report):
    with open(report) as fh:
        for line in fh:
            found = re.match(r"^\d+ features \([^)]*\): (.*)$", line.strip())
            if found:
                return [f.strip() for f in found.group(1).split(",")]
    sys.exit(f"{report}: no features line")


def table(path):
    with open(path) as fh:
        return {r["taxon"]: r for r in csv.DictReader(fh, delimiter="\t")}


def number(text):
    try:
        return float(text)
    except (TypeError, ValueError):
        return None


def main():
    prefix, label_a, label_b = sys.argv[1:4]
    specs = sys.argv[4:]
    wanted = None
    if specs and specs[0] == "--features-from":
        wanted = model_features(specs[1])
        specs = specs[2:]
    pairs = []  # (sample, row a, row b)
    only = []
    for spec in specs:
        sample, path_a, path_b = spec.split(":")
        a, b = table(path_a), table(path_b)
        for taxon in sorted(set(a) | set(b), key=lambda t: int(t) if t.isdigit() else t):
            if taxon in a and taxon in b:
                pairs.append((sample, a[taxon], b[taxon]))
            else:
                row = a.get(taxon) or b.get(taxon)
                only.append({"sample": sample, "taxon": taxon, "taxon_name": row["taxon_name"], "truth": row["truth"],
                             "in": label_a if taxon in a else label_b, "prediction": row["prediction"],
                             "probability": row["probability"]})
    if not pairs:
        sys.exit("no rows in both runs")
    features = [c for c in pairs[0][1].keys() if c not in SKIP and not c.startswith("meta_")]
    if wanted is not None:
        missing = [f for f in wanted if f not in features]
        if missing:
            print(f"model features without a column: {', '.join(missing)}", file=sys.stderr)
        features = [f for f in features if f in wanted]
    scale = {}
    for f in features:
        values = sorted(v for v in (number(r[f]) for _, r, _ in pairs) if v is not None)
        if len(values) < 4:
            scale[f] = 1.0
            continue
        q = statistics.quantiles(values, n=4)
        scale[f] = (q[2] - q[0]) or (max(values) - min(values)) or 1.0
    rows = []
    for f in features:
        changed = 0
        moves = []
        means = {(lab, t): [] for lab in ("a", "b") for t in ("1", "0")}
        for _, ra, rb in pairs:
            va, vb = number(ra[f]), number(rb[f])
            if va is None or vb is None:
                continue
            changed += va != vb
            moves.append(abs(vb - va) / scale[f])
            means[("a", ra["truth"])].append(va)
            means[("b", rb["truth"])].append(vb)
        mean = lambda xs: statistics.fmean(xs) if xs else float("nan")
        rows.append({"feature": f, "rows": len(moves), "changed_share": changed / len(moves) if moves else 0,
                     f"{label_a}_present": mean(means[("a", "1")]), f"{label_b}_present": mean(means[("b", "1")]),
                     f"{label_a}_absent": mean(means[("a", "0")]), f"{label_b}_absent": mean(means[("b", "0")]),
                     "mean_scaled_change": mean(moves)})
    rows.sort(key=lambda r: -r["mean_scaled_change"])
    columns = ["feature", "rows", "changed_share", f"{label_a}_present", f"{label_b}_present", f"{label_a}_absent",
               f"{label_b}_absent", "mean_scaled_change"]
    with open(prefix + ".features.tsv", "w") as fh:
        fh.write("\t".join(columns) + "\n")
        for r in rows:
            fh.write("\t".join(f"{r[c]:.6g}" if isinstance(r[c], float) else str(r[c]) for c in columns) + "\n")
    flips = []
    for sample, ra, rb in pairs:
        if ra["prediction"] == rb["prediction"]:
            continue
        moved = sorted(((abs(number(rb[f]) - number(ra[f])) / scale[f], f) for f in features
                        if number(ra[f]) is not None and number(rb[f]) is not None), reverse=True)[:5]
        flips.append({"sample": sample, "taxon": ra["taxon"], "taxon_name": ra["taxon_name"], "truth": ra["truth"],
                      f"p_{label_a}": ra["probability"], f"p_{label_b}": rb["probability"],
                      "moved_most": ", ".join(f"{f} {ra[f]}->{rb[f]}" for _, f in moved)})
    columns = ["sample", "taxon", "taxon_name", "truth", f"p_{label_a}", f"p_{label_b}", "moved_most"]
    with open(prefix + ".flips.tsv", "w") as fh:
        fh.write("\t".join(columns) + "\n")
        for r in flips:
            fh.write("\t".join(str(r[c]) for c in columns) + "\n")
    with open(prefix + ".rows.tsv", "w") as fh:
        columns = ["sample", "taxon", "taxon_name", "truth", "in", "prediction", "probability"]
        fh.write("\t".join(columns) + "\n")
        for r in only:
            fh.write("\t".join(str(r[c]) for c in columns) + "\n")
    for label, idx in ((label_a, 1), (label_b, 2)):
        tp = sum(p[idx]["truth"] == "1" and p[idx]["prediction"] == "1" for p in pairs)
        fp = sum(p[idx]["truth"] == "0" and p[idx]["prediction"] == "1" for p in pairs)
        fn = sum(p[idx]["truth"] == "1" and p[idx]["prediction"] == "0" for p in pairs)
        print(f"{label}: rows in both {len(pairs)}, TP {tp} FP {fp} FN {fn}")
    print(f"rows in one run only: {len(only)} ({sum(r['in'] == label_a for r in only)} {label_a}, "
          f"{sum(r['in'] == label_b for r in only)} {label_b}); called there: "
          f"{sum(r['prediction'] == '1' for r in only)}")
    gained = sum(f["truth"] == "1" and number(f[f"p_{label_b}"]) > number(f[f"p_{label_a}"]) for f in flips)
    print(f"calls flipped: {len(flips)} (present taxa {sum(f['truth'] == '1' for f in flips)}, absent "
          f"{sum(f['truth'] == '0' for f in flips)}; present taxa called in {label_b} only: {gained})")
    dp = [number(rb["probability"]) - number(ra["probability"]) for _, ra, rb in pairs]
    print(f"probability change: mean {statistics.fmean(dp):+.5f}, mean |change| {statistics.fmean(abs(x) for x in dp):.5f}, "
          f"rows changed {sum(x != 0 for x in dp)} of {len(dp)}")
    print("features that moved most (mean change in IQR units; share of rows changed):")
    for r in rows[:15]:
        print(f"  {r['feature']}: {r['mean_scaled_change']:.4f} ({r['changed_share']:.3f})")


if __name__ == "__main__":
    main()
