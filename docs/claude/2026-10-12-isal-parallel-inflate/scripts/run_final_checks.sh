#!/bin/bash
D=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/0ff4a1d7-3599-4615-a178-e4ca253109af/scratchpad/isal_parallel
bash $D/build.sh > /dev/null
X="taskset -c 0-3 nice -n 5 $D/bin/isalpar"
G=~/gzpar/data/R1.member.gz
B=$D/logs/R1.bounds.tsv
date
cat $G ~/gzpar/data/R1.fq > /dev/null
# every deflate block its own piece (pieceMiB tiny), with the reference check
for cfg in "spec 2 4" "spec 3 4" "lead 2 4" "spec 2 1"; do
  set -- $cfg
  $X par $G $B 0.001 $3 $1 $2 0 ~/gzpar/data/R1.fq
done | tee $D/logs/every_block_pieces.log | awk -F'\t' '{print $3,$5,"threads",$9,"pieces",$13,"spec",$15,"lead",$17,"GBps",$21,"OK",$NF}'
$X cost $G $B 4 - $D/logs/markers_4MiB_b.tsv | tee $D/logs/cost_4MiB_linetypes.log
date
