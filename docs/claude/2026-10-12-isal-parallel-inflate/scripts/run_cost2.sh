#!/bin/bash
D=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/0ff4a1d7-3599-4615-a178-e4ca253109af/scratchpad/isal_parallel
X="taskset -c 0-3 nice -n 5 $D/bin/isalpar"
G=~/gzpar/data/R1.member.gz
cat $G > /dev/null
for s in 1 16 4; do
  echo "load: $(cat /proc/loadavg)"
  $X seq $G | head -1
  $X cost $G $D/logs/R1.bounds.tsv $s - $D/logs/markers_${s}MiB.tsv
done 2>&1 | tee $D/logs/cost_quiet.log
