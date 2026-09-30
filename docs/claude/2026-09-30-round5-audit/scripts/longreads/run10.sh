#!/bin/bash
# Experiment 10: identity spread of ONT reads against the depth identity margin (0.04): share of a
# species' aligned bases below its own-identity threshold, and profiles against the sampled truth,
# for 1.6% and 5% ONT errors and HiFi.
S=$(dirname "$0")
W=~/audit6/longreads
OUT=$S/out_e10.txt
{
for x in "e1/out/hifi e1/hifi.truth.pkl" "e1/out/ont e1/ont.truth.pkl" "e8/out_xd0/ont5 e8/ont5.truth.pkl"; do
  set -- $x
  echo "== $1"
  head -1 $W/$1.profile.log | tr '\t' '\n' | grep -n -E "low_identity|top_identity|VerticalCoverage|vcov|depth|abundance|identity" | tr '\n' ' '; echo
  cut -f1-3 $W/$1.profile
  python3 - $W/$1.profile.log $W/$2 <<'EOF'
import sys, pickle, collections
log, truth = sys.argv[1], sys.argv[2]
rows = [l.rstrip("\n").split("\t") for l in open(log)]
h = rows[0]
want = [c for c in h if any(k in c.lower() for k in ("identity", "coverage", "abund", "name", "taxid", "species", "hits"))]
for r in rows[1:]:
    d = dict(zip(h, r))
    print("   ", {k: d[k] for k in want[:12]})
T = pickle.load(open(truth, "rb"))
b = collections.Counter()
for t in T["reads"]:
    b[t["acc"]] += t["length"]
lens = {"GCA_999001002.1": 280092, "GCF_999002001.1": 269979, "GCA_999003003.1": 268800}
cov = {a: b[a] / lens[a] for a in b}
tot = sum(cov.values())
print("    truth shares", {a: round(c / tot, 3) for a, c in cov.items()}, "coverage", {a: round(c, 2) for a, c in cov.items()})
EOF
done
} > $OUT 2>&1
echo done
