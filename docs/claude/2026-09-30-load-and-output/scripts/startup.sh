#!/bin/bash
# Start-up and shut-down of a run on 1000 pairs, at 1 and 8 threads: protal's own timers, wall
# clock, peak RSS and page faults.  startup.sh READS_DIR DB [DB ...]
# READS_DIR holds *_R1.fq.gz and *_R2.fq.gz; the first 1000 pairs are used.
source "$(dirname "$0")/env.sh"
src=$1; shift
r1=$(ls $src/*_R1.fq.gz | head -1); r2=$(ls $src/*_R2.fq.gz | head -1)
zcat $r1 | head -n 4000 | gzip -1 > $OUT/tiny_R1.fq.gz
zcat $r2 | head -n 4000 | gzip -1 > $OUT/tiny_R2.fq.gz
for db in "$@"; do
  for t in 1 8; do
    rm -rf $OUT/o
    echo "== $(basename $db) t=$t (load $(cut -d' ' -f1 /proc/loadavg))"
    /usr/bin/time -v $BIN --db $db -1 $OUT/tiny_R1.fq.gz -2 $OUT/tiny_R2.fq.gz -o $OUT/o -t $t --no_qcmsa --verbose \
      > $OUT/startup.log 2> $OUT/startup.time
    grep -E "^(Preload genomes took|Load Index took|Aligning reads took|Profiling took|Run protal took)" $OUT/startup.log
    grep -E "Elapsed|Maximum resident|User time|System time|Minor" $OUT/startup.time
  done
done
