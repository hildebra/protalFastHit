#!/usr/bin/env python3
"""The mock local test's runs side by side: per sample the calls and scores of each run (score_mocks.py score's
<run>/scores/samples.tsv), every Escherichia species' probability and call in the samples whose truth has E. coli (the
profiles' .profile.log), and the long reads E. coli lost to its congeners (../2026-10-13-ecoli-hifi/lost_reads.py on each
run's SAM: kept, and lost by locus, cycle or mixed).

usage: compare.py MOCKS RUN:DB [RUN:DB ...]     (DB: the run's protal_db folder, for its genome2tiid.tsv)
"""
import csv
import glob
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
LOST = os.path.join(HERE, "..", "2026-10-13-ecoli-hifi", "lost_reads.py")
ECOLI_SAMPLES = ("d6331", "atcc", "strain_mixes")  # sample name prefixes whose truth has E. coli
LONG = ("atcc1003_hifi_SRR11606871", "d6331_hifi_SRR13128012", "d6331_ont_ERR15459768")


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def ecoli_taxid(db):
    with open(os.path.join(db, "genome2tiid.tsv")) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[3].endswith(";s__Escherichia coli"):
                return f[1]
    sys.exit(f"no E. coli in {db}/genome2tiid.tsv")


def main():
    mocks = sys.argv[1]
    dbs = dict(spec.split(":", 1) for spec in sys.argv[2:])
    runs = list(dbs)
    print("# per sample: called tp/fp/fn F1 L1 by run")
    scores = {run: {r["sample"]: r for r in read_tsv(os.path.join(mocks, run, "scores", "samples.tsv"))} for run in runs}
    samples = sorted(set().union(*(s.keys() for s in scores.values())))
    print("sample\t" + "\t".join(runs))
    for sample in samples:
        cells = []
        for run in runs:
            r = scores[run].get(sample)
            cells.append("-" if r is None else
                         f"{r['called']} {r['tp']}/{r['fp']}/{r['fn']} F1 {r['f1'][:5]} L1 {r['l1'][:5]}")
        print(sample + "\t" + "\t".join(cells))
    print("\n# Escherichia: probability (* called) by run")
    for sample in samples:
        if not sample.startswith(ECOLI_SAMPLES):
            continue
        rows = {}
        for run in runs:
            for path in glob.glob(os.path.join(mocks, run, "**", sample + ".profile.log"), recursive=True):
                for r in read_tsv(path):
                    if ";g__Escherichia;" in r["Lineage"]:
                        name = r["Lineage"].rsplit("s__Escherichia ", 1)[-1]
                        rows.setdefault(name, {})[run] = f"{float(r['Probability']):.3f}" + ("*" if r["Predicted"] == "1" else "")
        for name in sorted(rows):
            print(f"{sample}\t{name}\t" + "\t".join(rows[name].get(run, "-") for run in runs))
    print("\n# E. coli's long reads: kept | lost to congeners by locus/cycle/mixed, by run")
    for sample in LONG:
        cells = []
        for run in runs:
            sams = glob.glob(os.path.join(mocks, run, "**", sample + ".sam.zst"), recursive=True)
            if not sams:
                cells.append("-")
                continue
            g2t = os.path.join(dbs[run], "genome2tiid.tsv")
            out = subprocess.run([sys.executable, "-I", LOST, sams[0], g2t, ecoli_taxid(dbs[run])],
                                 capture_output=True, text=True, check=True).stdout.splitlines()
            kept = out[0].rsplit(" ", 1)[-1]
            lost = {line.split("\t")[0]: line.split("\t")[1] for line in out[2:]}
            cells.append(f"{kept} | {lost.get('locus', 0)}/{lost.get('cycle', 0)}/{lost.get('mixed', 0)}")
        print(sample + "\t" + "\t".join(cells))


if __name__ == "__main__":
    main()
