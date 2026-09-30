#!/bin/bash
# Reads whose gene at a read end went to another species although their own has the gene.
S=$(dirname "$0")
W=~/audit6/longreads
P=~/strain-build/bin/protal
G=~/strain-build/mini_db/gtdb_r226
D=$W/mini_db
mkdir -p $W/e2c && cd $W/e2c
pick() {  # fq name out
  awk -v n="@$2" 'NR%4==1{keep=($1==n)} keep' $1 > $3.fq
}
pick $W/e2/o130000.fq o130000_81 r81
pick $W/e1/ont.fq ont_124 r124
pick $W/e1/ont.fq ont_152 r152
pick $W/e1/ont.fq ont_300 r300
cat r81.fq r124.fq r152.fq r300.fq > four.fq
$P --db $D -1 four.fq --read_type ont --model_ont $D/model_pe.xml -o out -t 1 --sam_format sam --no_qcmsa --no_strains -m 10 --prefix four > four.log 2>&1
echo "rc=$?"
grep -v '^@' out/four.sam | awk '{h=$6; sub(/[0-9]+[MXID].*[MXID]/,"~",h); print $1,$2,$3,$4,$5,h, length($10)}'
python3 - <<'EOF'
import pickle, sys
sys.path.insert(0, "$(dirname "$0")")
for pk, names in (("/home/falk/audit6/longreads/e2/o130000.truth.pkl", ["o130000_81"]), ("/home/falk/audit6/longreads/e1/ont.truth.pkl", ["ont_124", "ont_152", "ont_300"])):
    T = pickle.load(open(pk, "rb"))
    for t in T["reads"]:
        if t["name"] in names:
            print(t["name"], t["acc"], t["start"], t["length"], t["strand"], t["read_length"])
EOF
