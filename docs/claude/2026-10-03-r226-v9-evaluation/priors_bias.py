#!/usr/bin/env python3
"""Do the cluster priors separate present from absent taxa because of biology or because of how the simulated species
were chosen? The simulation pool (download_gtdb.py) favours species with several genomes (strains to simulate from);
GTDB's species are mostly singleton clusters. Shares of singleton clusters among present and absent taxa, the FN rate
of present singleton-cluster species against multi-genome species by fragments, and the FP rate of absent taxa by
cluster size, on v9's test set (final model) and training data (species held out)."""
import csv, gzip, sys
from collections import Counter, defaultdict

root = sys.argv[1]


def rows_of(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def pct(a, b):
    return f"{100 * a / b:.1f}%" if b else "-"


def size_bin(log10_genomes):
    g = round(10 ** log10_genomes) if log10_genomes >= 0 else 0
    return "1 genome" if g <= 1 else "2-4" if g <= 4 else "5-19" if g <= 19 else ">=20"


def frag_bin(f):
    return "1-2" if f <= 2 else "3-10" if f <= 10 else ">10"


def analyse(label, preds, pcol, dump):
    print(f"\n# {label}")
    comp = defaultdict(Counter)
    fn = defaultdict(Counter)
    fp = defaultdict(Counter)
    for r in preds:
        d = dump.get((r["meta_sample"], r["taxon"]))
        if d is None:
            continue
        tr, called = r["truth"] == "1", float(r[pcol]) >= 0.5
        sb = size_bin(float(d["cluster_genomes_log10"]))
        comp["present" if tr else "absent"][sb] += 1
        if tr:
            fb = frag_bin(float(d["fragments"]))
            fn[(sb, fb)]["present"] += 1
            fn[(sb, fb)]["FN"] += not called
        else:
            fp[sb]["absent"] += 1
            fp[sb]["FP"] += called
    print("cluster size\tpresent taxa\tshare\tabsent taxa\tshare")
    for sb in ("1 genome", "2-4", "5-19", ">=20"):
        print(f"{sb}\t{comp['present'][sb]}\t{pct(comp['present'][sb], sum(comp['present'].values()))}\t{comp['absent'][sb]}\t{pct(comp['absent'][sb], sum(comp['absent'].values()))}")
    print("\nFN rate of present taxa by cluster size and fragments:\ncluster size\t1-2 fragments\t3-10\t>10")
    for sb in ("1 genome", "2-4", "5-19", ">=20"):
        cells = []
        for fb in ("1-2", "3-10", ">10"):
            c = fn[(sb, fb)]
            cells.append(f"{pct(c['FN'], c['present'])} ({c['FN']}/{c['present']})")
        print(f"{sb}\t" + "\t".join(cells))
    print("\nFP rate of absent taxa by cluster size:\ncluster size\tabsent\tFP\tFP rate")
    for sb in ("1 genome", "2-4", "5-19", ">=20"):
        c = fp[sb]
        print(f"{sb}\t{c['absent']}\t{c['FP']}\t{pct(c['FP'], c['absent'])}")


test_dump = {(r["meta_sample"], r["taxon"]): r for r in rows_of(f"{root}/v9/test/training_data.tsv")}
analyse("v9 test set (pe, final model)", rows_of(f"{root}/v9/model_logs/trained_model.test_predictions.tsv.gz"), "p", test_dump)
train_dump = {(r["meta_sample"], r["taxon"]): r for r in rows_of(f"{root}/v9/training/training_data.tsv")}
analyse("v9 training data (pe, species held out)", rows_of(f"{root}/v9/model_logs/trained_model.predictions.tsv.gz"), "p_species", train_dump)

# The database's own distribution: every species that appears as a taxon anywhere in the training dump, once.
seen = {}
for d in train_dump.values():
    seen[d["taxon"]] = size_bin(float(d["cluster_genomes_log10"]))
dist = Counter(seen.values())
print(f"\n# species appearing as taxa in the training dump ({len(seen)}), by cluster size: " + ", ".join(f"{k} {pct(v, len(seen))}" for k, v in sorted(dist.items())))
