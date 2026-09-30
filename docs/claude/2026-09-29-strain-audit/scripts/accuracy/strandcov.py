#!/usr/bin/env python3
"""strandcov.py - per-position read coverage by strand and by fragment from protal's SAM.

Usage: strandcov.py <protal_out_dir> <sample> [<sample> ...]
For each sample and each gene (all species), computes per reference position: forward reads, reverse
reads, distinct fragments (read names). Prints, by distance from the nearer gene end, the fraction of
positions with >= 4 reads that are covered by one strand only, and the ratio reads/fragments (how
often overlapping mates count one molecule twice). Also the fraction of positions with >= 2 reads
covered by a single fragment (the snp_min_cov 2 + strand filter can then pass on one molecule).
"""
import gzip, os, re, sys
from collections import defaultdict
import numpy as np

CIG = re.compile(r"(\d+)([MIDNSHPX=])")


def main():
    out = sys.argv[1]
    for sample in sys.argv[2:]:
        lens = {}
        fwd, rev, frag = {}, {}, {}
        names = defaultdict(lambda: defaultdict(set))
        with gzip.open(os.path.join(out, "alignments", sample + ".sam.gz"), "rt") as fh:
            for line in fh:
                if line.startswith("@SQ"):
                    f = line.split("\t")
                    n = f[1][3:]
                    L = int(f[2][3:])
                    lens[n] = L
                    fwd[n] = np.zeros(L, np.int32); rev[n] = np.zeros(L, np.int32); frag[n] = np.zeros(L, np.int32)
                    continue
                if line[0] == "@":
                    continue
                f = line.split("\t", 7)
                flag = int(f[1])
                if flag & 4:
                    continue
                g = f[2]
                p0 = int(f[3]) - 1
                rl = sum(int(n) for n, op in CIG.findall(f[5]) if op in "MXD=")
                arr = rev[g] if flag & 16 else fwd[g]
                arr[p0:p0 + rl] += 1
                names[g][f[0]].add((p0, p0 + rl))
        for g, d in names.items():
            for nm, iv in d.items():
                m = np.zeros(lens[g], bool)
                for a, b in iv:
                    m[a:b] = True
                frag[g] += m
        bins = [(0, 50), (50, 100), (100, 200), (200, 300), (300, 450), (450, 10 ** 9)]
        acc = {b: [0, 0, 0, 0, 0] for b in bins}  # positions>=4 reads, one-strand, reads, frags, pos>=2reads single-frag
        tot2, single2 = 0, 0
        for g in lens:
            L = lens[g]
            tot = fwd[g] + rev[g]
            dist = np.minimum(np.arange(L), L - 1 - np.arange(L))
            for b in bins:
                sel = (dist >= b[0]) & (dist < b[1])
                s4 = sel & (tot >= 4)
                acc[b][0] += int(s4.sum())
                acc[b][1] += int((s4 & ((fwd[g] == 0) | (rev[g] == 0))).sum())
                acc[b][2] += int(tot[sel].sum())
                acc[b][3] += int(frag[g][sel].sum())
            s2 = tot >= 2
            tot2 += int(s2.sum())
            single2 += int((s2 & (frag[g] == 1)).sum())
        print(f"== {out} {sample}: positions with >=2 reads covered by one fragment only: {single2}/{tot2} = {single2 / max(tot2, 1):.3f}")
        print("dist_from_end\tpos>=4reads\tone_strand_frac\treads_per_fragment")
        for b in bins:
            a = acc[b]
            print(f"{b[0]}-{b[1]}\t{a[0]}\t{a[1] / max(a[0], 1):.3f}\t{a[2] / max(a[3], 1):.3f}")


if __name__ == "__main__":
    main()
