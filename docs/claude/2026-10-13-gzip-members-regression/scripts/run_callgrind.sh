#!/bin/bash
# Where the CPU goes on a file of 64 KiB members: callgrind over 64 MiB of the FASTQ, 1 and 2 gzip threads.
set -uo pipefail
W=~/pgzdiag
cd $W/data
[ -s small64k.gz ] || head -c 64M R1.fq | split -b 64K --filter='gzip -6' > small64k.gz
[ -s small1.gz ] || head -c 64M R1.fq | gzip -6 > small1.gz
ls -la small64k.gz small1.gz
cd $W
for t in 1 2; do
  taskset -c 0-3 nice -n 5 valgrind --tool=callgrind --callgrind-out-file=cg.small64k.t$t.out ./bench_pgz data/small64k.gz $t 2>&1 | grep -E "^threads|Collected"
  callgrind_annotate --inclusive=no cg.small64k.t$t.out 2>/dev/null | head -30
done
