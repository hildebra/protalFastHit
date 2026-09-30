#!/bin/bash
# The long reads that went through the short-read path (e4): coordinates against the truth; and the
# e1 truth's genome shares (for the profiles).
S=$(dirname "$0")
W=~/audit6/longreads
G=~/strain-build/mini_db/gtdb_r226
D=$W/mini_db
OUT=$S/out_e4b.txt
cd $W/e4
{
for s in l20 l70; do
  echo "######## $s in the se run"
  grep -v '^@' out/mixed_se.sam | awk -v p="^${s}_" '$1 ~ p' > $s.only.sam
  python3 $S/eval_sam.py $s.only.sam $s.fq $s.truth.pkl --gtdb $G --db $D --show 5 | cut -c1-250
done
python3 - <<'EOF'
import pickle, collections
for f in ("/home/falk/audit6/longreads/e1/hifi.truth.pkl", "/home/falk/audit6/longreads/e1/ont.truth.pkl"):
    T = pickle.load(open(f, "rb"))
    b = collections.Counter()
    for t in T["reads"]:
        b[t["acc"]] += t["length"]
    lens = {"GCA_999001002.1": 280092, "GCF_999002001.1": 269979, "GCA_999003003.1": 268800}
    cov = {a: b[a] / lens[a] for a in b}
    tot = sum(cov.values())
    print(f, {a: round(c / tot, 3) for a, c in cov.items()})
EOF
} > $OUT 2>&1
echo done
