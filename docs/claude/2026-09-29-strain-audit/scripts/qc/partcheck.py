#!/usr/bin/env python3
"""qcmsa output partition vs its input: each output gene block of every row must be a subsequence
(columns removed, order kept) of the same gene's block of the same row in the raw MSA, with gene
cells that qcmsa masks turned into '-'."""
import sys, os, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from consistency import read_fasta, parse_part

def subseq(small, big):
    it = iter(big)
    return all(any(c == b or c == "-" for b in it) for c in small)

raw_dir, qc_dir = sys.argv[1], sys.argv[2]
for out in sorted(glob.glob(os.path.join(qc_dir, "*.msa.fna"))):
    sp = os.path.basename(out)[:-len(".msa.fna")]
    rn, rs = read_fasta(os.path.join(raw_dir, sp + ".raw.msa.fna"))
    raw = dict(zip(rn, rs))
    rp = {g: (s, e) for g, s, e in parse_part(os.path.join(raw_dir, sp + ".raw.partition.txt"))}
    qn, qs = read_fasta(out)
    qp = parse_part(os.path.join(qc_dir, sp + ".partition.txt"))
    bad = 0
    for n, s in zip(qn, qs):
        for g, a, b in qp:
            ra, rb = rp[g]
            if not subseq(s[a - 1:b], raw[n][ra - 1:rb]):
                bad += 1
    print(f"{sp}: {len(qn)} rows x {len(qp)} genes, blocks not a subsequence of their raw gene: {bad}")
