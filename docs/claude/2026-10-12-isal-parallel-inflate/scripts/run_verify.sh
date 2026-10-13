#!/bin/bash
D=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/0ff4a1d7-3599-4615-a178-e4ca253109af/scratchpad/isal_parallel
X="taskset -c 0-3 nice -n 5 $D/bin/isalpar"
G=~/gzpar/data/R1.member.gz
cat $G > /dev/null; cat ~/gzpar/data/R1.fq > /dev/null
for s in 1 4 16; do
 for sch in spec lead; do
  for v in 2 3; do
   for t in 1 4; do
     $X par $G $D/logs/R1.bounds.tsv $s $t $sch $v 0 ~/gzpar/data/R1.fq 2>&1
   done
  done
 done
done | tee $D/logs/verify.log | awk -F'\t' '{print $3,$5,$7,$9,"pieces",$13,"spec",$15,"lead",$17,"GBps",$21,"rss",$27,"OK",$NF}'
