#!/bin/bash
# Alignment stage (--no_profile) of two builds, alternated: OLD, OLD with huge pages from the glibc
# tunable, and NEW. ab_alignment.sh OLD_BIN NEW_BIN DB READS_DIR THREADS REPS
source "$(dirname "$0")/env.sh"
old=$1; new=$2; db=$3; src=$4; t=$5; reps=$6
r1=$(ls $src/*_R1.fq* | head -1); r2=$(ls $src/*_R2.fq* | head -1)
out=$PERF_DIR/runs/ab; mkdir -p $out
for rep in $(seq 1 $reps); do
  for v in old oldthp new; do
    bin=$old; [ $v = new ] && bin=$new
    tun=""; [ $v = oldthp ] && tun=glibc.malloc.hugetlb=1
    rm -rf $out/o
    load=$(cut -d' ' -f1 /proc/loadavg)
    GLIBC_TUNABLES=$tun $bin --db $db -1 $r1 -2 $r2 -o $out/o -t $t --no_qcmsa --no_profile > $out/log 2>&1
    echo -e "$(basename $src)\tt=$t\trep=$rep\t$v\t$(grep 'Aligning reads took' $out/log | sed 's/Aligning reads took //')\tload=$load"
  done
done
