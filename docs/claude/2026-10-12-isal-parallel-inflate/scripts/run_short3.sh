#!/bin/bash
D=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/0ff4a1d7-3599-4615-a178-e4ca253109af/scratchpad/isal_parallel
X="taskset -c 0-3 nice -n 5 $D/bin/isalpar"
paste -d '\n' <(head -n 200000 ~/gzpar/data/R1.fq) <(base64 -w 150 /dev/urandom | head -n 200000) > $D/data/mixed.txt
taskset -c 0-3 nice -n 5 gzip -6 -c $D/data/mixed.txt > $D/data/mixed.gz
ls -la $D/data/mixed.*
$X bounds $D/data/mixed.gz $D/logs/mixed.bounds.tsv
awk -F'\t' 'NR==1{total=$6; next} {o[++n]=$2} END{o[n+1]=total; for(i=1;i<=n;i++){d=o[i+1]-o[i]; if(d==0) z++; else if(d<32768) s++} print "blocks", n, "empty", z+0, "nonempty_under_32KiB", s+0}' $D/logs/mixed.bounds.tsv
for cfg in "spec 2 4" "spec 3 4" "lead 2 4" "spec 2 1"; do
  set -- $cfg
  $X par $D/data/mixed.gz $D/logs/mixed.bounds.tsv 0 $3 $1 $2 0 $D/data/mixed.txt
done | tee $D/logs/short_pieces_mixed.log | awk -F'\t' '{print $3,$5,"threads",$9,"pieces",$13,"spec",$15,"lead",$17,"markers",$51,"OK",$NF}'
