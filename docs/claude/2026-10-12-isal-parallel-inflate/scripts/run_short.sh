#!/bin/bash
D=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/0ff4a1d7-3599-4615-a178-e4ca253109af/scratchpad/isal_parallel
X="taskset -c 0-3 nice -n 5 $D/bin/isalpar"
mkdir -p $D/data
head -c 40000000 ~/gzpar/data/R1.fq > $D/data/small.fq
taskset -c 0-3 nice -n 5 pigz -b 32 -p 2 -c $D/data/small.fq > $D/data/small.pigz32.gz
ls -la $D/data
$X bounds $D/data/small.pigz32.gz $D/logs/small.bounds.tsv
awk -F'\t' 'NR>1{ if (NR>2) { d=$2-p; if (min==""||d<min) min=d; if (d<32768) s++ } p=$2 } END{print "min_block_out", min, "blocks_under_32KiB", s+0, "of", NR-1}' $D/logs/small.bounds.tsv
for cfg in "spec 2 4" "spec 3 4" "lead 2 4" "lead 2 1" "spec 2 1"; do
  set -- $cfg
  $X par $D/data/small.pigz32.gz $D/logs/small.bounds.tsv 0.0001 $3 $1 $2 0 $D/data/small.fq
done | tee $D/logs/short_pieces.log | awk -F'\t' '{print $3,$5,"threads",$9,"pieces",$13,"spec",$15,"lead",$17,"OK",$NF}'
