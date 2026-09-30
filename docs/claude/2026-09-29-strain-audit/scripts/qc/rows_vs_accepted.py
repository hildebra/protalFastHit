#!/usr/bin/env python3
"""Per species: samples where the model accepts it (misc statistics) vs MSA rows vs truth (manifest)."""
import sys, os, csv, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from consistency import read_fasta

run = sys.argv[1]
man = {}
for r in csv.DictReader(open(os.path.join(run, "sim", "manifest.tsv")), delimiter="\t"):
    sp = "s__" + r["species"].replace(" ", "_")
    man.setdefault(sp, {}).setdefault(r["sample"], []).append((r["genome"], float(r["vertical_coverage"])))
for st in sorted(glob.glob(os.path.join(run, "protal", "misc", "*.statistics.tsv"))):
    sp = os.path.basename(st)[:-len(".statistics.tsv")]
    rows = list(csv.DictReader(open(st), delimiter="\t"))
    acc = [r["#SAMPLEID"] for r in rows if r["Accepted"] == "1"]
    msa = os.path.join(run, "protal", "strains", sp + ".raw.msa.fna")
    names = read_fasta(msa)[0] if os.path.exists(msa) else []
    qc = os.path.join(run, "protal", "strains", sp + ".msa.fna")
    qnames = read_fasta(qc)[0] if os.path.exists(qc) else []
    truth = man.get(sp, {})
    print(f"{sp}: true {len(truth)} accepted {len(acc)} raw-MSA rows {len(names) - (1 if names else 0)} qcmsa rows {len(qnames) - (1 if qnames else 0)}")
    for s in sorted(set(truth) | set(r['#SAMPLEID'] for r in rows)):
        g = truth.get(s, [])
        rr = [r for r in rows if r["#SAMPLEID"] == s]
        vc = rr[0]["VerticalCoverage"] if rr else "-"
        print(f"   {s:10} truth={[(a, round(c, 1)) for a, c in g]} protal_vcov={vc} accepted={s in acc} raw={s in names} qc={s in qnames}")
