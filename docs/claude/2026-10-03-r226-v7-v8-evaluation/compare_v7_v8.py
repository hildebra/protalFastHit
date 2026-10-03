#!/usr/bin/env python3
"""v7 (nad+divergence, knob curve) against v8 (same binary and samples, + sample depth, knob 0.5) on GTDB r226:
per read type the test set by depth, the held-out calls that flipped, and what the flipped taxa look like."""
import csv, gzip, math, sys
from collections import Counter, defaultdict

root = sys.argv[1]
TYPES = {"pe": "", "se": "_se", "pb": "_pb", "ont": "_ont"}


def rows_of(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def f1(tp, fp, fn):
    return 2 * tp / (2 * tp + fp + fn) if tp else 0.0


def confusion(rows, pcol, knob=0.5):
    tp = fp = fn = tn = 0
    for r in rows:
        t, c = r["truth"] == "1", float(r[pcol]) >= knob
        if t and c: tp += 1
        elif t: fn += 1
        elif c: fp += 1
        else: tn += 1
    return tp, fp, fn, tn


def best_threshold(rows, pcol):
    best = (0, 0.5)
    for k in [i / 100 for i in range(5, 96)]:
        tp, fp, fn, _ = confusion(rows, pcol, k)
        best = max(best, (f1(tp, fp, fn), k))
    return best


def depth_key(r):
    return int(r["meta_read_pairs"])


print("# Independent test set, by read type: v7 vs v8 at knob 0.5 (p of the final model) and the best threshold")
print("type\tversion\tTP\tFP\tFN\tF1@0.5\tbest thr\tF1@best")
test = {}
for t, suffix in TYPES.items():
    for v in ("v7", "v8"):
        rows = rows_of(f"{root}/{v}/model_logs/trained_model{suffix}.test_predictions.tsv.gz")
        test[(t, v)] = rows
        tp, fp, fn, _ = confusion(rows, "p")
        bf, bk = best_threshold(rows, "p")
        print(f"{t}\t{v}\t{tp}\t{fp}\t{fn}\t{f1(tp, fp, fn):.4f}\t{bk:.2f}\t{bf:.4f}")

print("\n# Test set by depth (pe, se): FP and FN at 0.5, v7 -> v8")
for t in ("pe", "se"):
    by = defaultdict(lambda: Counter())
    for v in ("v7", "v8"):
        for r in test[(t, v)]:
            d = depth_key(r)
            tr, c = r["truth"] == "1", float(r["p"]) >= 0.5
            by[d][f"present"] += tr if v == "v7" else 0
            by[d][f"absent"] += (not tr) if v == "v7" else 0
            if tr and not c: by[d][f"FN_{v}"] += 1
            if c and not tr: by[d][f"FP_{v}"] += 1
    print(f"\n{t}\tread pairs\tpresent\tabsent\tFP v7\tFP v8\tFN v7\tFN v8")
    for d in sorted(by):
        c = by[d]
        print(f"\t{d}\t{c['present']}\t{c['absent']}\t{c['FP_v7']}\t{c['FP_v8']}\t{c['FN_v7']}\t{c['FN_v8']}")

print("\n# Test set (pe): the taxa whose call flipped between v7 and v8 at 0.5, joined with the dump's features")
dump = {}
for r in rows_of(f"{root}/v8/test/training_data.tsv"):
    dump[(r["meta_sample"], r["taxon"])] = r
p7 = {(r["meta_sample"], r["taxon"]): float(r["p"]) for r in test[("pe", "v7")]}
flips = Counter()
detail = defaultdict(list)
for r in test[("pe", "v8")]:
    key = (r["meta_sample"], r["taxon"])
    c7, c8 = p7[key] >= 0.5, float(r["p"]) >= 0.5
    if c7 == c8:
        continue
    tr = r["truth"] == "1"
    kind = ("TP lost" if tr else "FP removed") if c7 and not c8 else ("TP gained" if tr else "FP added")
    flips[kind] += 1
    d = dump.get(key)
    if d:
        detail[kind].append((int(d["meta_read_pairs"]), float(d["fragments"]), d["meta_rep_genome"], float(d["identity"]),
                             d["meta_novel_level"], d["meta_relative_rank"], float(d["em_own_share"]) if "em_own_share" in d else -1))
for kind, n in flips.most_common():
    rows = detail[kind]
    frag = sorted(x[1] for x in rows)
    med = frag[len(frag) // 2] if frag else 0
    strains = sum(1 for x in rows if x[2] == "0")
    deep = Counter(x[0] for x in rows)
    ident = sorted(x[3] for x in rows)
    print(f"{kind}: {n}; fragments median {med:g}, <=2 fragments {sum(1 for f in frag if f <= 2)}, strains (other genome than the rep) {strains}; "
          f"identity median {ident[len(ident) // 2] if ident else 0:.3f}; by read pairs {dict(sorted(deep.items()))}")
    if kind in ("TP lost", "FP removed"):
        lv = Counter(x[4] or "-" for x in rows)
        print(f"   novel level of the closest missing species: {dict(lv.most_common())}")

print("\n# Species held out (pe), p_species at 0.5: v7 vs v8 and the flips")
h7 = {(r["meta_sample"], r["taxon"]): r for r in rows_of(f"{root}/v7/model_logs/trained_model.predictions.tsv.gz")}
h8 = rows_of(f"{root}/v8/model_logs/trained_model.predictions.tsv.gz")
for v, rows in (("v7", list(h7.values())), ("v8", h8)):
    tp, fp, fn, _ = confusion(rows, "p_species")
    print(f"{v}: TP {tp} FP {fp} FN {fn} F1 {f1(tp, fp, fn):.4f}")
flips = Counter()
fnbins = defaultdict(Counter)
fpbins = defaultdict(Counter)
for r in h8:
    key = (r["meta_sample"], r["taxon"])
    a = h7[key]
    c7, c8 = float(a["p_species"]) >= 0.5, float(r["p_species"]) >= 0.5
    tr = r["truth"] == "1"
    frag = float(r["fragments"])
    fb = "1" if frag <= 1 else "2" if frag <= 2 else "3-10" if frag <= 10 else "11-100" if frag <= 100 else ">100"
    depth = int(r["meta_read_pairs"])
    if tr:
        fnbins[fb]["present"] += 1
        if not c7: fnbins[fb]["FN v7"] += 1
        if not c8: fnbins[fb]["FN v8"] += 1
        if r["meta_rep_genome"] == "0":
            fnbins[fb]["strains"] += 1
            if not c8: fnbins[fb]["FN v8 strains"] += 1
    else:
        fpbins[depth]["absent"] += 1
        if c7: fpbins[depth]["FP v7"] += 1
        if c8: fpbins[depth]["FP v8"] += 1
        if frag <= 1:
            fpbins[depth]["absent 1 frag"] += 1
            if c7: fpbins[depth]["FP v7 1 frag"] += 1
            if c8: fpbins[depth]["FP v8 1 frag"] += 1
    if c7 != c8:
        kind = ("TP lost" if tr else "FP removed") if c7 else ("TP gained" if tr else "FP added")
        flips[kind] += 1
print("flips:", dict(flips.most_common()))
print("\nFN by the present taxon's fragments:\nfragments\tpresent\tFN v7\tFN v8\tstrains\tFN v8 strains")
for fb in ("1", "2", "3-10", "11-100", ">100"):
    c = fnbins[fb]
    print(f"{fb}\t{c['present']}\t{c['FN v7']}\t{c['FN v8']}\t{c['strains']}\t{c['FN v8 strains']}")
print("\nFP by sample depth (and of absent taxa with one fragment):\nread pairs\tabsent\tFP v7\tFP v8\tabsent 1 frag\tFP v7 1 frag\tFP v8 1 frag")
for d in sorted(fpbins):
    c = fpbins[d]
    print(f"{d}\t{c['absent']}\t{c['FP v7']}\t{c['FP v8']}\t{c['absent 1 frag']}\t{c['FP v7 1 frag']}\t{c['FP v8 1 frag']}")
