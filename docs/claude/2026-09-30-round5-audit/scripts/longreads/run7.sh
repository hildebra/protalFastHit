#!/bin/bash
# Experiment 7: a rerun that finds the SAM of a PacBio run and is told --read_type ont (not
# --profile_only): which model, which -a / --snp_min_af, any warning? And -m 0 on long reads.
S=$(dirname "$0")
W=~/audit6/longreads
P=${P:-~/strain-build/bin/protal}
G=~/strain-build/mini_db/gtdb_r226
D=$W/mini_db
OUT=$S/out_e7.txt
mkdir -p $W/e7 && cd $W/e7
{
cp $W/e1/hifi.fq .
$P --db $D -1 hifi.fq --read_type pb --model_pb $D/model_pe.xml -o out --prefix hifi -t 2 --sam_format sam --no_qcmsa --no_strains > first.log 2>&1; echo "first run (pb) rc=$?"
grep -c . out/hifi.sam
cp out/hifi.profile.log first.profile.log
$P --db $D -1 hifi.fq --read_type ont --model_ont $D/model_pe.xml -o out --prefix hifi -t 2 --sam_format sam --no_qcmsa --no_strains > rerun.log 2>&1; echo "rerun (--read_type ont, SAM exists) rc=$?"
grep -E "Skip|Warning|warning|Model of|holds|profiled as|All alignments" rerun.log | cut -c1-250
head -c 0 /dev/null
grep '^@CO' out/hifi.sam
echo "== --profile_only of the same SAM with --read_type ont"
$P --db $D --profile_only out/hifi.sam --read_type ont --model_ont $D/model_pe.xml -o po --prefix hifi -t 2 --no_qcmsa --no_strains > po.log 2>&1; echo "rc=$?"
grep -E "Warning|holds|Model of" po.log | cut -c1-250
echo "== -m 0 and -m 3 on long reads: records"
for m in 0 1 3; do
  $P --db $D -1 hifi.fq --read_type pb --model_pb $D/model_pe.xml -o out_m$m --prefix hifi -t 2 --sam_format sam --no_qcmsa --no_strains --no_profile -m $m > m$m.log 2>&1
  echo "-m $m rc=$? records=$(grep -vc '^@' out_m$m/hifi.sam) secondary=$(grep -v '^@' out_m$m/hifi.sam | awk 'and($2,256)' | wc -l) max_per_segment=$(grep -v '^@' out_m$m/hifi.sam | awk '{if (!and($2,256)) {if (c>mx) mx=c; c=1} else c++} END{if (c>mx) mx=c; print mx}')"
done
echo "== options print for an ONT sample"
grep -E "max score ani|snp min af" $W/e1/ont.log
} > $OUT 2>&1
echo done
