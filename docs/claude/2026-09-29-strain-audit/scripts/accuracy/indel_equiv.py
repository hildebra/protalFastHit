#!/usr/bin/env python3
"""indel_equiv.py - are protal's indel representations sequence-equivalent to the truth?

Usage: indel_equiv.py <run> <protal_out_dir> <min_depth>
For each true indel event in a pure-sample row at >= min_depth, compare the row's local sequence
(ungapped, reference window [p-8, q+8) including insertion columns) with the true local sequence.
'equivalent' = same bases (placement may differ inside a repeat); 'different' = the row spells another
sequence (misplaced indel, or indel missed); 'incomplete' = window holds N/'-'-only coverage gaps.
Also prints where the missed events sit (distance to the gene end).
"""
import os, sys, pickle
from collections import Counter
import numpy as np
sys.path.insert(0, os.path.expanduser("~/audit5/accuracy"))
import evaluate as ev

FL = 8


def main():
    run, out, mind = sys.argv[1], os.path.abspath(sys.argv[2]), float(sys.argv[3])
    man = ev.load_manifest(run)
    res = Counter()
    edge = Counter()
    examples = []
    for sp in sorted({s for v in man.values() for s in v}):
        name = "s__" + sp.replace(" ", "_")
        p = os.path.join(out, "strains", name + ".raw.msa.fna")
        if not os.path.exists(p):
            continue
        names, seqs = ev.read_fasta(p)
        refrow = seqs[0]
        rows = dict(zip(names[1:], seqs[1:]))
        parts = ev.read_partition(os.path.join(out, "strains", name + ".raw.partition.txt"))
        for gid, s, e in parts:
            seg = refrow[s:e + 1]
            # column index (within gene) of each reference position, and insertion columns before it
            refidx = [i for i, c in enumerate(seg) if c != "-"]
            for sample, gl in man.items():
                if sp not in gl or len(gl[sp]) != 1 or gl[sp][0][1] < mind or sample not in rows:
                    continue
                genome = gl[sp][0][0]
                t = ev.TR.get((genome, gid))
                if t is None or genome.startswith("GCF"):
                    continue
                row = rows[sample][s:e + 1]
                L = len(t["tb"])
                events = []
                q = 0
                while q < L:
                    if t["tb"][q] == "-":
                        r = q
                        while r < L and t["tb"][r] == "-":
                            r += 1
                        events.append(("del", q, r))
                        q = r
                    else:
                        q += 1
                events += [("ins", p0, p0) for p0 in t["ins"]]
                for typ, a, b in events:
                    lo, hi = max(0, a - FL), min(L, b + FL)
                    # truth local sequence
                    tr = []
                    for k in range(lo, hi):
                        if k in t["ins"] and k > lo:
                            tr.append(t["ins"][k])
                        if t["tb"][k] != "-":
                            tr.append(t["tb"][k])
                    if typ == "ins" and a == lo:
                        tr.insert(0, t["ins"][a])
                    tr = "".join(tr)
                    # row local sequence: columns from refidx[lo] (incl. preceding insertion cols if lo==a) to refidx[hi-1]
                    c0 = refidx[lo]
                    if typ == "ins" and a == lo:
                        c0 = refidx[lo - 1] + 1 if lo > 0 else 0
                    c1 = refidx[hi - 1]
                    loc = row[c0:c1 + 1]
                    refloc = seg[c0:c1 + 1]
                    if any(ch == "N" for ch in loc) or any(ch == "-" and rc != "-" for ch, rc in zip(loc, refloc)) and typ == "ins":
                        # '-' at a reference column: deletion or coverage gap
                        pass
                    body = "".join(ch for ch in loc if ch != "-")
                    if "N" in loc or (loc.count("-") - refloc.count("-") > (b - a if typ == "del" else 0) + 3):
                        st = "incomplete"
                    elif body == tr:
                        st = "equivalent"
                    else:
                        st = "different"
                    res[(typ, st)] += 1
                    if st == "different":
                        edge[(typ, min(a, L - b) < 30)] += 1
                        if len(examples) < 6:
                            examples.append((sample, sp, gid, typ, a, b - a, tr, body, loc, refloc))
    print(f"== {run} rows >= {mind}x: indel events by local-sequence agreement")
    for k in sorted(res):
        print(k, res[k])
    print("different events within 30 bp of a gene end:", dict(edge))
    for x in examples:
        print(x)


if __name__ == "__main__":
    main()
