#!/usr/bin/env python3
"""How much of protal's depth comes from overlapping mates counted twice?
For each gene: depth counting reads (as protal does) vs depth counting fragments (a position covered
by both mates of a pair counts once). Only primary records with MAPQ >= 4.

usage: overlap_depth.py <sam.gz> [min_cov]
"""
import sys, re, gzip, collections
sam, min_cov = sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 2
by_gene = collections.defaultdict(lambda: collections.defaultdict(list))  # gene -> qname -> intervals
pairs = overlapping = 0
with gzip.open(sam, "rt") as f:
    for line in f:
        if line.startswith("@"): continue
        t = line.split("\t")
        flag = int(t[1])
        if flag & 0x904 or int(t[4]) < 4: continue
        start = int(t[3]) - 1
        length = sum(int(n) for n, op in re.findall(r"(\d+)([MXD=])", t[5]))
        by_gene[t[2]][t[0]].append((start, start + length))
read_depth_bases = frag_depth_bases = only_by_overlap = covered = 0
for gene, frags in by_gene.items():
    L = max(e for iv in frags.values() for s, e in iv)
    rd = [0] * L; fd = [0] * L
    for q, iv in frags.items():
        if len(iv) == 2:
            pairs += 1
            (s1, e1), (s2, e2) = iv
            if min(e1, e2) > max(s1, s2): overlapping += 1
        seen = set()
        for s, e in iv:
            for p in range(s, e):
                rd[p] += 1
                seen.add(p)
        for p in seen: fd[p] += 1
    for p in range(L):
        read_depth_bases += rd[p]; frag_depth_bases += fd[p]
        if rd[p] >= min_cov: covered += 1
        if rd[p] >= min_cov and fd[p] < min_cov: only_by_overlap += 1
print(f"pairs on one gene: {pairs}, overlapping: {overlapping} ({100 * overlapping / max(pairs, 1):.1f}%)")
print(f"read-depth bases {read_depth_bases}, fragment-depth bases {frag_depth_bases} "
      f"(excess {100 * (read_depth_bases - frag_depth_bases) / frag_depth_bases:.2f}%)")
print(f"positions with depth >= {min_cov}: {covered}; reaching it only through a double-counted overlap: {only_by_overlap}")
