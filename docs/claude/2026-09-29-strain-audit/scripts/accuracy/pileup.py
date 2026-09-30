#!/usr/bin/env python3
"""pileup.py - explain wrong / ambiguous / blanked MSA cells from protal's own alignments.

Usage: pileup.py <label> <protal_out_dir>
Reads results/<label>.errors.tsv (written by evaluate.py) and, per sample, <out>/alignments/<sample>.sam.gz.
ART read names carry the source genome (<accession>_contig<N>-<id>), so every read at a position can be
attributed to the sample's own genome or to another genome (cross-mapping). Writes
results/<label>.pileup.tsv with, per listed position: reads, distinct fragments, reads from other genomes,
support of the called/true/ref base (reads, fwd, rev, fragments, from own genome), and prints a summary.
"""
import gzip, os, re, sys
from collections import defaultdict, Counter
import pandas as pd

ACC = os.path.expanduser("~/audit5/accuracy")
CIG = re.compile(r"(\d+)([MIDNSHPX=])")


def taxids():
    t = {}
    for line in open(os.path.expanduser("~/audit5/world/protal_db/genome2tiid.tsv")):
        f = line.rstrip("\n").split("\t")
        t[f[3].split(";s__")[-1]] = f[1]
    return t


def main():
    label, out = sys.argv[1], sys.argv[2]
    err = pd.read_csv(os.path.join(ACC, "results", label + ".errors.tsv"), sep="\t")
    tax = taxids()
    res = []
    for sample, g in err.groupby("sample"):
        want = defaultdict(dict)  # refname -> pos -> row index
        for i, r in g.iterrows():
            want[f"{tax[r.species]}_{r.gene}"][int(r.pos)] = i
        pile = defaultdict(list)  # row index -> [(base, fwd, frag, src)]
        with gzip.open(os.path.join(out, "alignments", sample + ".sam.gz"), "rt") as fh:
            for line in fh:
                if line[0] == "@":
                    continue
                f = line.split("\t", 11)
                rn = f[2]
                if rn not in want:
                    continue
                flag = int(f[1])
                if flag & 4:
                    continue
                pos0 = int(f[3]) - 1
                cig = CIG.findall(f[5])
                reflen = sum(int(n) for n, op in cig if op in "MXD=")
                wp = want[rn]
                if not any(pos0 <= p < pos0 + reflen for p in wp):
                    continue
                seq = f[9]
                name = f[0]
                src = name.split("_contig")[0]
                fwd = not (flag & 16)
                q, rp = 0, pos0
                for n, op in cig:
                    n = int(n)
                    if op in "MX=":
                        for k in range(n):
                            if rp + k in wp:
                                pile[wp[rp + k]].append((seq[q + k], fwd, name, src))
                        q += n; rp += n
                    elif op in "IS":
                        q += n
                    elif op == "D":
                        for k in range(n):
                            if rp + k in wp:
                                pile[wp[rp + k]].append(("-", fwd, name, src))
                        rp += n
        for i, r in g.iterrows():
            obs = pile.get(i, [])
            own = r.genome
            def sup(b):
                o = [x for x in obs if x[0] == b]
                return (len(o), sum(1 for x in o if x[1]), sum(1 for x in o if not x[1]),
                        len({x[2] for x in o}), sum(1 for x in o if x[3] == own))
            called_b = r.called if r.called in "ACGT" else None
            rec = dict(r)
            rec.update(reads=len(obs), frags=len({x[2] for x in obs}), other_src=sum(1 for x in obs if x[3] != own),
                       other_genomes=";".join(sorted({x[3] for x in obs if x[3] != own})))
            for tag, b in (("called", called_b), ("true", r.true), ("refb", r.ref)):
                s = sup(b) if b else (0, 0, 0, 0, 0)
                rec.update({f"{tag}_reads": s[0], f"{tag}_fwd": s[1], f"{tag}_rev": s[2], f"{tag}_frags": s[3], f"{tag}_own": s[4]})
            res.append(rec)
    df = pd.DataFrame(res)
    df.to_csv(os.path.join(ACC, "results", label + ".pileup.tsv"), sep="\t", index=False)
    # summary
    df["cross"] = df.other_src > 0
    print(f"=== {label}: pileup of {len(df)} positions")
    for call, g in df.groupby("call"):
        print(f"-- {call}: n={len(g)}  with cross-mapped reads: {g.cross.mean():.2f}  "
              f"called base from other genomes only: {((g.called_own == 0) & (g.called_reads > 0)).mean():.2f}  "
              f"median reads {g.reads.median():.0f}  median frags {g.frags.median():.0f}")
    n = df[df.call == "N"]
    if len(n):
        n = n.assign(one_strand=(n.true_fwd == 0) | (n.true_rev == 0), lt2=n.true_reads < 2)
        print("-- N at true SNPs: true allele reads<2: %.2f, true allele on one strand only: %.2f, "
              "true-allele reads all own genome: %.2f" % (n.lt2.mean(), (n.one_strand & ~n.lt2).mean(),
                                                         (n.true_own == n.true_reads).mean()))
        print(n.groupby("depth").apply(lambda x: pd.Series(dict(n=len(x), lt2=x.lt2.mean(),
                                                                 one_strand=(x.one_strand & ~x.lt2).mean(),
                                                                 cross=x.cross.mean()))).to_string())


if __name__ == "__main__":
    main()
