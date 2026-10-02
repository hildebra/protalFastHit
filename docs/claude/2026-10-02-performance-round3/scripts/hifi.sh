#!/usr/bin/env bash
# HiFi reads (scripts/hifi_reads.py of HEAD) of ~3 Mb of templates from community rl150_p1000_s_1 (the Nanopore
# sample's), and the reference against the prototype on them, as for the other long reads.
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/1680d2f3-dfad-4cd7-9113-3d5aefceb3b5/scratchpad
W=$HOME/mt-work/perf3; mkdir -p $W/hifi
perl $S/perf3/templates.pl $HOME/bench071/samples/points/rl150_p1000/sim/manifests/rl150_p1000_s_1.tsv 3000000 1 > $W/hifi/templates.fa
$HOME/micromamba/envs/protal-db-build/bin/python $W/ref/scripts/hifi_reads.py --templates $W/hifi/templates.fa --out $W/hifi/hifi_3M.fq.gz --seed 1
echo "reads: $(zcat $W/hifi/hifi_3M.fq.gz | awk 'NR % 4 == 2 { n++; b += length($0) } END { print n, b }')"
