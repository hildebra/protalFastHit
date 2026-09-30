#!/bin/bash
# The one wrong-species record of the 5%-error ONT reads (both x_drop settings): debug build.
S=$(dirname "$0")
W=~/audit6/longreads
PA=$W/build-asan/protal
G=~/strain-build/mini_db/gtdb_r226
D=$W/mini_db
cd $W/e8
{
python3 $S/eval_sam.py out_xd0/ont5.sam ont5.fq ont5.truth.pkl --gtdb $G --db $D --exact 0.95 --show 3 | grep -A1 -E "^wrong_taxon	" | cut -c1-200
r=$(python3 $S/eval_sam.py out_xd0/ont5.sam ont5.fq ont5.truth.pkl --gtdb $G --db $D --exact 0.95 --show 3 | grep -A1 -E "^wrong_taxon	" | tail -1 | sed "s/^    \['\([^']*\)'.*/\1/")
g=$(python3 $S/eval_sam.py out_xd0/ont5.sam ont5.fq ont5.truth.pkl --gtdb $G --db $D --exact 0.95 --show 3 | grep -A1 -E "^wrong_taxon	" | tail -1 | awk -F"'" '{print $6}' | cut -d_ -f2)
echo "read $r gene $g"
awk -v n="@$r" 'NR%4==1{k=($1==n)} k' ont5.fq > one.fq
rm -rf dbg_one
PROTAL_LR_DEBUG=1 ASAN_OPTIONS=detect_leaks=0 $PA --db $D -1 one.fq --read_type ont --model_ont $D/model_pe.xml -o dbg_one --prefix one -t 1 --sam_format sam --no_qcmsa --no_strains --no_profile -m 10 --x_drop 0 > one.log 2>&1
grep -E "^LRDBG" one.log | grep -E "read|_$g |_$g:|align-fail" | cut -c1-300
grep -v '^@' dbg_one/one.sam | awk -v g="_$g\$" '$3 ~ g {print $1,$2,$3,$4,$5,length($10)}'
} > $S/out_e8c.txt 2>&1
