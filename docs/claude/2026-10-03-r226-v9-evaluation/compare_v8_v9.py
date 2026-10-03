#!/usr/bin/env python3
"""v8 (a3e397d, sigma 1.3, 20% held out, nad+depth+divergence) against v9 (26b065c, sigma 1.3/2.0, 30% held out, 30M
point, + unfiltered + priors): the test sets differ (more species missing from v9's database), so compare rates by
class, not F1: sensitivity by fragments and strain status, false positives per absent taxon by class and per 100
simulated missing species, both by depth; the new features' values by class of taxon."""
import csv, gzip, sys
from collections import Counter, defaultdict

root = sys.argv[1]
TYPES = {"pe": "", "se": "_se", "pb": "_pb", "ont": "_ont"}


def rows_of(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def f1(tp, fp, fn):
    return 2 * tp / (2 * tp + fp + fn) if tp else 0.0


def pct(a, b):
    return f"{100 * a / b:.2f}%" if b else "-"


def median(xs):
    xs = sorted(xs)
    return xs[len(xs) // 2] if xs else float("nan")


def frag_bin(f):
    return "1" if f <= 1 else "2" if f <= 2 else "3-10" if f <= 10 else "11-100" if f <= 100 else ">100"


print("# Independent test sets: composition and rates at knob 0.5 (test predictions joined with the dump)")
print("type\tver\tpresent\tabsent\tabsent beside a missing species\tTP\tFP\tFN\tF1\tsensitivity\tFP per absent\tFP per absent beside missing\tFP per other absent\tbest thr")
data = {}
for t, suffix in TYPES.items():
    for v in ("v8", "v9"):
        preds = rows_of(f"{root}/{v}/model_logs/trained_model{suffix}.test_predictions.tsv.gz")
        dump = {(r["meta_sample"], r["taxon"]): r for r in rows_of(f"{root}/{v}/test/training_data{suffix}.tsv")}
        rows = []
        for r in preds:
            d = dump.get((r["meta_sample"], r["taxon"]))
            if d is None:
                continue
            rows.append((r, d))
        data[(t, v)] = rows
        tp = fp = fn = 0
        absent = present = beside = fp_beside = fp_other = 0
        for r, d in rows:
            tr, c = r["truth"] == "1", float(r["p"]) >= 0.5
            miss = (d["meta_novel_level"] or "") != ""  # the closest species in the sample is one the database lacks
            if tr:
                present += 1
                if c: tp += 1
                else: fn += 1
            else:
                absent += 1
                beside += miss
                if c:
                    fp += 1
                    fp_beside += miss
                    fp_other += not miss
        best = max((f1(*[sum(1 for r, _ in rows if (r["truth"] == "1") and (float(r["p"]) >= k)),
                          sum(1 for r, _ in rows if (r["truth"] != "1") and (float(r["p"]) >= k)),
                          sum(1 for r, _ in rows if (r["truth"] == "1") and (float(r["p"]) < k))]), k)
                   for k in [i / 50 for i in range(5, 46)])
        print(f"{t}\t{v}\t{present}\t{absent}\t{beside}\t{tp}\t{fp}\t{fn}\t{f1(tp, fp, fn):.4f}\t{pct(tp, present)}\t{pct(fp, absent)}\t{pct(fp_beside, beside)}\t{pct(fp_other, absent - beside)}\t{best[1]:.2f}")

print("\n# Test set by depth (pe): FN rate of present taxa, FP per 100 absent taxa, v8 -> v9")
for t in ("pe", "se"):
    by = defaultdict(Counter)
    for v in ("v8", "v9"):
        for r, d in data[(t, v)]:
            dep = int(r["meta_read_pairs"])
            tr, c = r["truth"] == "1", float(r["p"]) >= 0.5
            by[dep][f"present_{v}"] += tr
            by[dep][f"absent_{v}"] += not tr
            by[dep][f"FN_{v}"] += tr and not c
            by[dep][f"FP_{v}"] += c and not tr
    print(f"\n{t}\tread pairs\tFN rate v8\tFN rate v9\tFP/100 absent v8\tFP/100 absent v9\t(FN v8/v9, FP v8/v9)")
    for dep in sorted(by):
        c = by[dep]
        print(f"\t{dep}\t{pct(c['FN_v8'], c['present_v8'])}\t{pct(c['FN_v9'], c['present_v9'])}\t{100 * c['FP_v8'] / max(1, c['absent_v8']):.2f}\t{100 * c['FP_v9'] / max(1, c['absent_v9']):.2f}\t({c['FN_v8']}/{c['FN_v9']}, {c['FP_v8']}/{c['FP_v9']})")

print("\n# Test set (pe): misses by the present taxon's fragments and whether it is a strain (another genome than the representative)")
print("fragments\tpresent v8\tFN v8\tFN rate v8\tpresent v9\tFN v9\tFN rate v9\tstrains v9\tFN rate strains v9")
bins = defaultdict(Counter)
for v in ("v8", "v9"):
    for r, d in data[("pe", v)]:
        if r["truth"] != "1":
            continue
        fb = frag_bin(float(d["fragments"]))
        c = float(r["p"]) >= 0.5
        bins[fb][f"present_{v}"] += 1
        bins[fb][f"FN_{v}"] += not c
        if d["meta_rep_genome"] == "0":
            bins[fb][f"strain_{v}"] += 1
            bins[fb][f"FN_strain_{v}"] += not c
for fb in ("1", "2", "3-10", "11-100", ">100"):
    c = bins[fb]
    print(f"{fb}\t{c['present_v8']}\t{c['FN_v8']}\t{pct(c['FN_v8'], c['present_v8'])}\t{c['present_v9']}\t{c['FN_v9']}\t{pct(c['FN_v9'], c['present_v9'])}\t{c['strain_v9']}\t{pct(c['FN_strain_v9'], c['strain_v9'])}")

print("\n# v9 test set (pe): the new features by class of taxon (medians; priors -1 = unknown)")
classes = defaultdict(list)
for r, d in data[("pe", "v9")]:
    tr = r["truth"] == "1"
    miss = (d["meta_novel_level"] or "") != ""
    frag = float(d["fragments"])
    cls = ("present strain" if d["meta_rep_genome"] == "0" else "present representative") if tr else ("absent beside a missing species" if miss else "absent, other")
    classes[cls].append(d)
    if tr and frag <= 2:
        classes["present, <= 2 fragments"].append(d)
    if not tr and frag <= 2:
        classes["absent, <= 2 fragments"].append(d)
feats = ["failed_candidate_rate", "fragments_all", "em_fragments", "fragments", "em_own_share", "cluster_ani_radius", "cluster_mean_ani", "cluster_genomes_log10", "rep_contamination", "rep_duplicate_share"]
print("class\tn\t" + "\t".join(feats))
for cls, rows in classes.items():
    print(f"{cls}\t{len(rows)}\t" + "\t".join(f"{median(float(d[f]) for d in rows):.3g}" for f in feats))
known = sum(1 for r, d in data[("pe", "v9")] if float(d["cluster_ani_radius"]) >= 0)
print(f"\ncluster priors known for {known} of {len(data[('pe', 'v9')])} test rows; rep_duplicate_share > 0 in {sum(1 for r, d in data[('pe', 'v9')] if float(d['rep_duplicate_share']) > 0)} rows")
fcr = [(float(d["failed_candidate_rate"]), r["truth"] == "1") for r, d in data[("pe", "v9")]]
for lo, hi in ((0, 0.0001), (0.0001, 0.2), (0.2, 0.5), (0.5, 0.8), (0.8, 1.01)):
    sel = [t for f, t in fcr if lo <= f < hi]
    print(f"failed_candidate_rate in [{lo}, {hi}): {len(sel)} rows, present {pct(sum(sel), len(sel))}")
