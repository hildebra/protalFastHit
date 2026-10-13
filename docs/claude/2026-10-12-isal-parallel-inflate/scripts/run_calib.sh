#!/bin/bash
# N concurrent independent sequential ISA-L decodes (each its own process) on cores 0-3
D=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/0ff4a1d7-3599-4615-a178-e4ca253109af/scratchpad/isal_parallel
G=~/gzpar/data/R1.member.gz
cat $G > /dev/null
for r in 1 2 3; do
 for N in 1 2 3 4; do
  for i in $(seq 1 $N); do
    taskset -c 0-3 nice -n 5 $D/bin/isalpar seq $G > /tmp/calib_$i.txt &
  done
  wait
  agg=$(cat /tmp/calib_*.txt | awk -F'\t' '$2=="gzip_crc"{s+=$12; n++; m=(m==""||$12<m)?$12:m} END{printf "%d\tsum_GBps\t%.3f\tmean_per_proc\t%.3f\tmin\t%.3f", n, s, s/n, m}')
  echo -e "round\t$r\tN\t$agg"
  rm -f /tmp/calib_*.txt
 done
done | tee $D/logs/calib_concurrent_seq.tsv
