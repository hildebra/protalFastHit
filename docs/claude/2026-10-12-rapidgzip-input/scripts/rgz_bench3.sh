#!/bin/bash
cd ~/rgzbuild
F=~/gzpar/data/R1.member.gz
for round in 1 2; do
  for cfg in "4 1 4" "4 2 4" "4 4 4" "2 1 4" "2 2 4"; do
    taskset -c 0-3 ./bench_rgz_x86-64 $F $cfg 0
  done
done
