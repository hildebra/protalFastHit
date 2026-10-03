#!/usr/bin/env python3
"""Summarise the v9 feature-group ablation (~/v9_eval/<tag><rt>.test_predictions.tsv.gz, scored on v9's test set):
F1 at knob 0.5 (and at the knob curve where one was fitted), best threshold, FP split by class (absent taxa beside a
species the database lacks / other absent taxa), FN by fragments and strain status."""
import csv, gzip, json, os, sys
from collections import Counter

out, tables = sys.argv[1], sys.argv[2]
TAGS = ["v8set", "no_unfiltered", "no_priors", "default", "no_depth", "noflag", "checkm"]
LABEL = {"v8set": "nad+depth+divergence (v8's set)", "no_unfiltered": "+ priors", "no_priors": "+ unfiltered",
         "default": "+ unfiltered + priors (default)", "no_depth": "default without depth",
         "noflag": "default, priors without the singleton flag", "checkm": "default, priors = CheckM + radius only"}


def rows_of(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def f1(tp, fp, fn):
    return 2 * tp / (2 * tp + fp + fn) if tp else 0.0


for rt in ("", "_se"):
    dump = {(r["meta_sample"], r["taxon"]): r for r in rows_of(f"{tables}/test/training_data{rt}.tsv")}
    print(f"\n# {'pe' if rt == '' else 'se'}: v9 test set ({len(dump)} taxa), knob 0.5")
    print("feature set\tF1\tcurve F1\tbest thr (F1)\tTP\tFP\tFN\tFP beside missing sp.\tFP other\tFN 1 frag\tFN 2 frag\tFN 3-10\tFN strains 1-10\tlog loss\tAP")
    for tag in TAGS:
        path = f"{out}/{tag}{rt}.test_predictions.tsv.gz"
        if not os.path.exists(path):
            print(f"{LABEL[tag]}\t(missing)")
            continue
        preds = rows_of(path)
        c = Counter()
        pts = []
        for r in preds:
            d = dump.get((r["meta_sample"], r["taxon"]))
            if d is None:
                continue
            tr, p = r["truth"] == "1", float(r["p"])
            pts.append((p, tr))
            called = p >= 0.5
            frag = float(d["fragments"])
            miss = (d["meta_novel_level"] or "") != ""
            if tr:
                if called: c["TP"] += 1
                else:
                    c["FN"] += 1
                    c["FN1" if frag <= 1 else "FN2" if frag <= 2 else "FN3" if frag <= 10 else "FNx"] += 1
                    if d["meta_rep_genome"] == "0" and frag <= 10: c["FNs"] += 1
            elif called:
                c["FP"] += 1
                c["FPm" if miss else "FPo"] += 1
        best = max((f1(sum(1 for p, t in pts if t and p >= k), sum(1 for p, t in pts if not t and p >= k), sum(1 for p, t in pts if t and p < k)), k) for k in [i / 50 for i in range(5, 46)])
        metrics = {}
        mpath = f"{out}/{tag}{rt}.metrics.json"
        if os.path.exists(mpath):
            with open(mpath) as fh:
                metrics = json.load(fh)
        test = metrics.get("test", {}).get("this one", {})
        curve = metrics.get("test_depth_knobs", {}).get("F1")
        print(f"{LABEL[tag]}\t{f1(c['TP'], c['FP'], c['FN']):.4f}\t{curve if curve is None else f'{curve:.4f}'}\t{best[1]:.2f} ({best[0]:.4f})\t{c['TP']}\t{c['FP']}\t{c['FN']}\t{c['FPm']}\t{c['FPo']}\t{c['FN1']}\t{c['FN2']}\t{c['FN3']}\t{c['FNs']}\t{test.get('log_loss', '')}\t{test.get('AP', '')}")
