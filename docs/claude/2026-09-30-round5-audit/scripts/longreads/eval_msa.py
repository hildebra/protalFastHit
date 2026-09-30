"""Compares the rows of a strain MSA (protal's <species>.msa.fna and .partition.txt) with the true
gene sequences of each sample's genome.

  eval_msa.py MSA PARTITIONS --gtdb G --db D --taxid T --sample NAME=ACCESSION ...

Only genes whose sequence in the sample's genome has the reference gene's length are compared
(substitutions only, so columns are positions). Per sample: bases equal to the truth, bases that
differ, true SNPs (truth != reference) called / missed, IUPAC codes, N and '-'."""
import argparse
import collections
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lrlib import World, read_fasta_list

IUPAC = set("RYSWKMBDHV")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("msa")
    ap.add_argument("partitions")
    ap.add_argument("--gtdb", required=True)
    ap.add_argument("--db", required=True)
    ap.add_argument("--taxid", type=int, required=True)
    ap.add_argument("--sample", action="append", default=[])
    a = ap.parse_args()
    w = World(a.gtdb, a.db)
    rows = dict(read_fasta_list(a.msa))
    parts = []
    with open(a.partitions) as fh:
        for line in fh:
            m = re.match(r"DNA, gene(\d+) = (\d+)-(\d+)", line.strip())
            if m:
                parts.append((int(m.group(1)), int(m.group(2)) - 1, int(m.group(3))))
    print("rows:", ", ".join(f"{k}({len(v)})" for k, v in rows.items()), "genes:", len(parts))
    for spec in a.sample:
        name, acc = spec.split("=")
        if name not in rows:
            print(f"{name}: no row in the MSA")
            continue
        row = rows[name]
        truth_genes = {}
        for m in w.markers.get(acc, []):
            truth_genes[w.geneid[m[0]]] = w.gene_seq(acc, m)
        c = collections.Counter()
        for gid, s, e in parts:
            ref = w.db_genes.get(f"{a.taxid}_{gid}")
            tru = truth_genes.get(gid)
            seg = row[s:e]
            if ref is None or len(seg) != len(ref):
                c["genes_msa_length_differs_from_ref"] += 1
                continue
            if tru is None or len(tru) != len(ref):
                c["genes_skipped_indel_or_absent"] += 1
                continue
            c["genes_compared"] += 1
            for x, t, r in zip(seg, tru, ref):
                if x in "-":
                    c["gap"] += 1
                elif x == "N":
                    c["N"] += 1
                elif x in IUPAC:
                    c["iupac"] += 1
                    if t != r:
                        c["true_snp_iupac"] += 1
                elif x == t:
                    c["correct"] += 1
                    if t != r:
                        c["true_snp_called"] += 1
                else:
                    c["wrong"] += 1
                    if t != r:
                        c["true_snp_wrong_base"] += 1 if x != r else 0
                        c["true_snp_missed_as_ref"] += 1 if x == r else 0
                    else:
                        c["false_snp"] += 1
                if t != r:
                    c["true_snp_sites"] += 1
        called = c["correct"] + c["wrong"]
        print(f"{name} ({acc}): " + ", ".join(f"{k}={v}" for k, v in sorted(c.items())) +
              (f", error_rate={c['wrong'] / called:.2e}" if called else ""))


if __name__ == "__main__":
    main()
