#!/bin/bash
# Re-evaluate e1 with the refined categories; look at o130000_81 (a partial gene at a read end
# given to another species).
set -u
S=$(dirname "$0")
W=~/audit6/longreads
P=~/strain-build/bin/protal
G=~/strain-build/mini_db/gtdb_r226
D=$W/mini_db
cd $W/e1
for s in hifi ont; do
  echo "######## $s"
  python3 $S/eval_sam.py out/$s.sam $s.fq $s.truth.pkl --gtdb $G --db $D --exact 0.95 --show 6
done
cd $W/e2
awk 'NR%4==1{keep=($1=="@o130000_81")} keep' o130000.fq > r81.fq
$P --db $D -1 r81.fq --read_type ont --model_ont $D/model_pe.xml -o out81 -t 1 --sam_format sam --no_qcmsa --no_strains -m 10 --prefix r81 > r81.log 2>&1; echo rc=$?
grep -v '^@' out81/r81.sam | awk '{print $1,$2,$3,$4,$5,substr($6,1,40), length($10)}' | tail -8
grep -c . out81/r81.sam
grep -E '^>(2|3)_35$' -A1 $D/reference.fna | awk '{print substr($0,1,60), length($0)}'
python3 - <<EOF
import pickle
T = pickle.load(open("o130000.truth.pkl","rb"))
t = [r for r in T["reads"] if r["name"]=="o130000_81"][0]
print(t["start"], t["length"], t["strand"], t["read_length"], t["blocks"][:2], t["blocks"][-2:])
EOF
grep -P "^GCF_999002001.1\t" $G/simulation/marker_positions.tsv | awk -v m=$(awk '$2==35{print $1}' $D/gene2geneid.tsv) '$2==m'
