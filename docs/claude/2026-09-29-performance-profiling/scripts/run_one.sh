#!/bin/bash
# One timed protal run: run_one.sh LABEL DB READS_DIR THREADS [protal args...]
# READS_DIR holds one *_R1.fq[.gz] and one *_R2.fq[.gz]. $BIN is the binary (env.sh).
# Writes $PERF_DIR/runs/LABEL/{stdout.log,stderr.log,time.txt,load.txt} and the outputs.
set -u
source "$(dirname "$0")/env.sh"
label=$1; db=$2; reads=$3; t=$4; shift 4
out=$PERF_DIR/runs/$label
rm -rf $out; mkdir -p $out
r1=$(ls $reads/*_R1.fq* | head -1); r2=$(ls $reads/*_R2.fq* | head -1)
uptime > $out/load.txt
/usr/bin/time -v -o $out/time.txt $BIN --db $db -1 $r1 -2 $r2 -o $out/out -t $t --no_qcmsa "$@" > $out/stdout.log 2> $out/stderr.log
rc=$?
uptime >> $out/load.txt
field() { grep "$1" $out/time.txt | awk '{print $NF}'; }
echo -e "$label\trc=$rc\twall=$(field 'Elapsed (wall')\tuser=$(field 'User time')\tsys=$(field 'System time')\tmaxrss_kb=$(field 'Maximum resident')\tload=$(head -1 $out/load.txt | sed 's/.*average: //')"
grep -E "took|Processed reads|Total alignments|Output alignments" $out/stdout.log | sed 's/^/    /'
