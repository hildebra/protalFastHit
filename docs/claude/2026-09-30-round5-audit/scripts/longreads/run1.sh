#!/bin/bash
# Experiment 1: HiFi and ONT reads of the mini community, aligned with --read_type pb / ont; SAMs
# checked against the truth (eval_sam.py).
set -u
S=$(dirname "$0")
W=~/audit6/longreads
P=~/strain-build/bin/protal
G=~/strain-build/mini_db/gtdb_r226
D=$W/mini_db
cd $W
mkdir -p e1 && cd e1
python3 $S/gen_reads.py community hifi --gtdb $G --db $D --reads 300 --platform hifi --seed 1 --name "m64001_e1/{}/ccs"
python3 $S/gen_reads.py community ont --gtdb $G --db $D --reads 300 --platform ont --seed 2 --median 20000 --sigma 0.6 --name "ont_{}"
/usr/bin/time -v $P --db $D -1 hifi.fq --read_type pb --model_pb $D/model_pe.xml -o out --prefix hifi -t 2 --sam_format sam --no_qcmsa > hifi.log 2> hifi.err; echo "hifi rc=$?"
/usr/bin/time -v $P --db $D -1 ont.fq --read_type ont --model_ont $D/model_pe.xml -o out --prefix ont -t 2 --sam_format sam --no_qcmsa > ont.log 2> ont.err; echo "ont rc=$?"
grep -E "Align the|settled|chunks|Model of|Maximum resident|Elapsed" hifi.log hifi.err ont.log ont.err
ls out
python3 $S/eval_sam.py out/hifi.sam hifi.fq hifi.truth.pkl --gtdb $G --db $D > hifi.eval; cat hifi.eval
python3 $S/eval_sam.py out/ont.sam ont.fq ont.truth.pkl --gtdb $G --db $D > ont.eval; cat ont.eval
