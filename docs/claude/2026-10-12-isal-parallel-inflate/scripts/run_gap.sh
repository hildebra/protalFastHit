#!/bin/bash
D=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/0ff4a1d7-3599-4615-a178-e4ca253109af/scratchpad/isal_parallel
X="taskset -c 0-3 nice -n 5 $D/bin/isalpar"
G=~/gzpar/data/R1.member.gz
cat $G > /dev/null
for r in 1 2 3; do
  for t in 2 4; do for gap in 2 4 7 12; do
    $X par $G $D/logs/R1.bounds.tsv 4 $t lead 2 $gap - | sed "s/^/round\t$r\t/"
  done; done
done > $D/logs/gap.tsv
awk -F'\t' '{delete kv; for(i=4;i<NF;i+=2) kv[$i]=$(i+1); k="t" kv["threads"] "_gap" kv["gap"]; g[k]=g[k] " " kv["GBps"]; l[k]=l[k] " " kv["lead_pieces"]; w[k]=w[k] " " kv["main_wait_s"]} END{for(k in g) print k, "GBps", g[k], "lead_pieces", l[k], "wait", w[k]}' $D/logs/gap.tsv | sort
