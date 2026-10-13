#!/bin/bash
# piece_buffers.sh's binaries V0, V1, V3 at 3 threads, 4 alternated rounds (results/piece_buffers_t3.txt).
cd ~/pgzread/data
for round in 1 2 3 4; do for v in V0 V1 V3; do printf 'r%s\t%s\t' $round $v; taskset -c 0-3 nice -n 5 ~/pgzread/bench_$v one.gz 3 | cut -f1,6-9; done; done
