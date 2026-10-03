#!/usr/bin/env python3
"""Where does each feature set's gain on v9's test set come from, by the taxon's GTDB cluster size (1 genome: a
singleton species, which the simulation can only present as its own representative) and by fragments? FN rate of
present taxa and FP rate of absent taxa per cell, for the ablation models in ~/v9_eval."""
import csv, gzip, os, sys
from collections import Counter, defaultdict

out, tables = sys.argv[1], sys.argv[2]
TAGS = [("v8set", "v8's set"), ("no_priors", "+ unfiltered"), ("no_unfiltered", "+ priors"), ("default", "default"),
        ("noflag", "default, priors without the singleton flag"), ("checkm", "default, priors = CheckM + radius only")]


def rows_of(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def pct(a, b):
    return f"{100 * a / b:.1f}%" if b else "-"


def size_bin(log10_genomes):
    g = round(10 ** log10_genomes) if log10_genomes >= 0 else 0
    return "1" if g <= 1 else "2-4" if g <= 4 else "5+"


for rt in ("", "_se"):
    dump = {(r["meta_sample"], r["taxon"]): r for r in rows_of(f"{tables}/test/training_data{rt}.tsv")}
    print(f"\n# {'pe' if rt == '' else 'se'}: v9 test set, knob 0.5. Cells: FN rate of present taxa (misses/present); FP rate of absent taxa (FP/absent)")
    header = ["feature set", "FN 1 genome, <=2 frag", "FN 2-4, <=2", "FN 5+, <=2", "FN 1 genome, >2 frag", "FN 2-4, >2", "FN 5+, >2",
              "FP 1 genome, <=2 frag", "FP 2-4, <=2", "FP 5+, <=2", "FP 1 genome, >2", "FP 2-4, >2", "FP 5+, >2"]
    print("\t".join(header))
    for tag, label in TAGS:
        path = f"{out}/{tag}{rt}.test_predictions.tsv.gz"
        if not os.path.exists(path):
            print(f"{label}\t(missing)")
            continue
        cells = defaultdict(Counter)
        for r in rows_of(path):
            d = dump.get((r["meta_sample"], r["taxon"]))
            if d is None:
                continue
            tr, called = r["truth"] == "1", float(r["p"]) >= 0.5
            key = (size_bin(float(d["cluster_genomes_log10"])), "<=2" if float(d["fragments"]) <= 2 else ">2")
            if tr:
                cells[("FN",) + key]["n"] += 1
                cells[("FN",) + key]["err"] += not called
            else:
                cells[("FP",) + key]["n"] += 1
                cells[("FP",) + key]["err"] += called
        row = [label]
        for kind in ("FN", "FP"):
            for fb in ("<=2", ">2"):
                for sb in ("1", "2-4", "5+"):
                    c = cells[(kind, sb, fb)]
                    row.append(f"{pct(c['err'], c['n'])} ({c['err']}/{c['n']})")
        print("\t".join(row))
