#!/usr/bin/env python3
"""The mates mate guidance added: primary records in the SAMs written with guidance (B) whose read and mate have none
in those written without it (C), over the paired-end test samples (time_guidance.sh keeps round 1's SAMs). How many,
how many partial (clipped: the mate runs past the gene's end), their identity, and whether they were simulated from
the genome of the taxon they are written on (read names carry the source genome; genome2tiid maps it to a taxid).

    python3 rescued_mates.py TIMING_DIR RUN
"""
import collections
import glob
import os
import re
import subprocess
import sys

timing, run = sys.argv[1], sys.argv[2]
taxid_of_genome = {}
with open(f"{run}/training_db/genome2tiid.tsv") as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) >= 2:
            taxid_of_genome[f[0]] = f[1]


def primaries(path):
    out = {}
    proc = subprocess.Popen(f"zstdcat '{path}' | grep -v '^@'", shell=True, stdout=subprocess.PIPE, text=True)
    for line in proc.stdout:
        f = line.split("\t", 6)
        flag = int(f[1])
        if flag & 4 or flag & 256:
            continue
        out[(f[0], flag & 0xC0)] = (f[2], f[5], int(f[4]))
    proc.wait()
    return out


cigar = re.compile(r"(\d+)([MIDSX=])")
c = collections.Counter()
identities = []
clipped_bases = []
for b in sorted(glob.glob(f"{timing}/sam.B/*.sam.zst")):
    a = primaries(b)
    n = primaries(f"{timing}/sam.C/{os.path.basename(b)}")
    for key, (rname, cig, mapq) in a.items():
        if key in n:
            c["in both"] += 1
            continue
        c["added by guidance"] += 1
        ops = collections.Counter()
        for k, op in cigar.findall(cig):
            ops[op] += int(k)
        aligned = ops["M"] + ops["="] + ops["X"] + ops["I"] + ops["D"]
        identities.append((ops["M"] + ops["="]) / max(1, aligned))
        if ops["S"]:
            c["partial (clipped)"] += 1
            clipped_bases.append(ops["S"])
        genome = re.match(r"(GC[AF]_\d+\.\d+)", key[0])
        if genome:
            c["from the taxon's own genome" if taxid_of_genome.get(genome.group(1)) == rname.split("_")[0]
              else "from another genome"] += 1
    c["only without guidance"] += sum(1 for k in n if k not in a)
print(dict(c))
if identities:
    identities.sort()
    q = lambda p: identities[int(p * (len(identities) - 1))]
    print(f"identity of the added mates: median {q(0.5):.3f}, 10th percentile {q(0.1):.3f}")
if clipped_bases:
    clipped_bases.sort()
    print(f"clipped bases of the partial ones: median {clipped_bases[len(clipped_bases) // 2]}")
