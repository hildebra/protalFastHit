#!/bin/bash
D=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/0ff4a1d7-3599-4615-a178-e4ca253109af/scratchpad/isal_parallel
X="taskset -c 0-3 nice -n 5 $D/bin/isalpar"
for cfg in "spec 2 4" "spec 3 4" "lead 2 4" "spec 2 1"; do
  set -- $cfg
  $X par $D/data/small.pigz32.gz $D/logs/small.bounds.tsv 0 $3 $1 $2 0 $D/data/small.fq
done | tee $D/logs/short_pieces_every_block.log | awk -F'\t' '{print $3,$5,"threads",$9,"pieces",$13,"spec",$15,"lead",$17,"OK",$NF}'
