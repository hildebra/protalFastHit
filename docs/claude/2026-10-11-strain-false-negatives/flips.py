"""Paired comparison of two runs' calls of the test strains in the shallow samples: strains called by one run only
(fixed / new misses of the second against the first) and an exact two-sided sign test on them.

usage: flips.py <strain dir> <run A> <run B> [<run A> <run B> ...]
"""
import math
import sys

S, PAIRS = sys.argv[1], sys.argv[2:]
truth = [line.rstrip("\n").split("\t") for line in open(f"{S}/reads/truth.tsv")][1:]


def called(run):
    out = {}
    for sample in sorted({r[0] for r in truth if r[1] == "shallow"}):
        with open(f"{S}/run_{run}/profiles/{sample}.profile.log") as f:
            ci = {c: i for i, c in enumerate(f.readline().rstrip("\n").split("\t"))}
            pred = {r[ci["Name"]]: r[ci["Predicted"]] == "1" for r in (l.rstrip("\n").split("\t") for l in f)
                    if len(r) > ci["Name"]}
        for r in truth:
            if r[0] == sample and r[4] == "strain":
                out[(sample, r[2])] = pred.get(r[3], False)
    return out


print("| A | B | strains | called A | called B | B only (fixed) | A only (new misses) | sign test p |")
print("|---|---|---|---|---|---|---|---|")
for a, b in zip(PAIRS[::2], PAIRS[1::2]):
    ca, cb = called(a), called(b)
    fixed = sum(cb[k] and not ca[k] for k in ca)
    new = sum(ca[k] and not cb[k] for k in ca)
    n, k = fixed + new, min(fixed, new)
    p = min(1.0, 2 * sum(math.comb(n, i) for i in range(k + 1)) / 2 ** n) if n else 1.0
    print(f"| {a} | {b} | {len(ca)} | {sum(ca.values())} | {sum(cb.values())} | {fixed} | {new} | {p:.3g} |")
