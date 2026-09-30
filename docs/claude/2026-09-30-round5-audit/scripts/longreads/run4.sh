#!/bin/bash
# Experiment 4: long reads after the first 100 short ones, given as single-end reads (no
# --read_type): the length check reads only the first 100. What does the short-read path do with
# reads of 20, 70 and 140 kb? Also: FASTA input of pb and ont reads (qualities), and short reads
# given as pb.
S=$(dirname "$0")
W=~/audit6/longreads
P=${P:-~/strain-build/bin/protal}
G=~/strain-build/mini_db/gtdb_r226
D=$W/mini_db
OUT=$S/out_e4.txt
mkdir -p $W/e4 && cd $W/e4
{
python3 $S/gen_reads.py sliding short --gtdb $G --db $D --acc GCF_999002001.1 --platform hifi --length 150 --offsets 1000:201000:2000 --name "s_{}" --seed 4
python3 $S/gen_reads.py sliding l20 --gtdb $G --db $D --acc GCF_999002001.1 --platform hifi --length 20000 --offsets 0:100000:20000 --name "l20_{}" --seed 5
python3 $S/gen_reads.py sliding l70 --gtdb $G --db $D --acc GCF_999002001.1 --platform hifi --length 70000 --offsets 0:140000:70000 --name "l70_{}" --seed 6
python3 $S/gen_reads.py sliding l140 --gtdb $G --db $D --acc GCF_999002001.1 --platform hifi --length 140000 --offsets 0:1:1 --name "l140_{}" --seed 7
cat short.fq l20.fq l70.fq l140.fq > mixed_se.fq
echo "reads: $(( $(wc -l < mixed_se.fq) / 4 ))"
$P --db $D -1 mixed_se.fq --model_se $D/model_pe.xml -o out -t 2 --sam_format sam --no_qcmsa --no_strains --prefix mixed_se > mixed_se.log 2>&1
echo "se (long reads after 100 short ones) rc=$?"
grep -E "Align the|Error|error|too long" mixed_se.log | cut -c1-200 | head
echo "records per read kind (not secondary):"
grep -v '^@' out/mixed_se.sam | awk '!and($2,256){split($1,a,"_"); n[a[1]]++; if (length($10)>maxl[a[1]]) maxl[a[1]]=length($10)} END{for (k in n) print k, n[k], "longest SEQ", maxl[k]}'
grep -v '^@' out/mixed_se.sam | awk '$1 ~ /^l/ {print $1,$2,$3,$4,$5,substr($6,1,60), length($10)}' | head -20
echo "invalid-alignment messages: $(grep -c 'Invalid' mixed_se.log)"
# the same with only 99 short reads before a long one: the check sees it
head -396 short.fq > s99.fq; cat s99.fq l70.fq > s99_long.fq
$P --db $D -1 s99_long.fq --model_se $D/model_pe.xml -o out -t 2 --sam_format sam --no_qcmsa --no_strains --prefix s99 > s99.log 2>&1
echo "se (a long read within the first 100) rc=$?"; grep -E "Error" s99.log | cut -c1-250
} > $OUT 2>&1
echo done
