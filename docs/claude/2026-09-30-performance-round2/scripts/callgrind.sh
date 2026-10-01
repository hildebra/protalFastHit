#!/bin/bash
# Instruction profile of a whole run on the first PAIRS pairs, 1 thread: callgrind.sh LABEL DB READS_DIR PAIRS
# Everything is collected: --toggle-collect patterns also match OpenMP's outlined
# RunPairedEnd._omp_fn and static initialisers, so collection ends up off in the parallel loop.
# Stages are separated afterwards by their inclusive costs (cg_stage.sh). ~4 min for 50k pairs.
set -u
source "$(dirname "$0")/env.sh"
label=$1; db=$2; reads=$3; pairs=$4
out=$PERF_DIR/cg/$label
rm -rf $out; mkdir -p $out/reads
zcat $reads/*_R1.fq.gz | head -n $((pairs*4)) | gzip -1 > $out/reads/sub_R1.fq.gz
zcat $reads/*_R2.fq.gz | head -n $((pairs*4)) | gzip -1 > $out/reads/sub_R2.fq.gz
/usr/bin/time -f "$label callgrind %e s %M KB" valgrind --tool=callgrind --callgrind-out-file=$out/callgrind.out \
  $PROF_BIN --db $db -1 $out/reads/sub_R1.fq.gz -2 $out/reads/sub_R2.fq.gz -o $out/out -t 1 --no_qcmsa \
  > $out/stdout.log 2> $out/stderr.log
tail -1 $out/stderr.log
bash $here/cg_stage.sh $label $pairs
