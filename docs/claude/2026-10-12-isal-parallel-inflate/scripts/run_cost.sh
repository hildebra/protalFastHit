#!/bin/bash
D=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/0ff4a1d7-3599-4615-a178-e4ca253109af/scratchpad/isal_parallel
cat ~/gzpar/data/R1.member.gz > /dev/null
for s in "$@"; do
  taskset -c 0-3 nice -n 5 $D/bin/isalpar cost ~/gzpar/data/R1.member.gz $D/logs/R1.bounds.tsv $s ~/gzpar/data/R1.fq $D/logs/markers_${s}MiB.tsv 2>&1 | tee $D/logs/cost_${s}MiB.log
done
