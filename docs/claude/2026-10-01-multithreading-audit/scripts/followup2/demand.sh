#!/usr/bin/env bash
# Input-path audit: can the input path feed N aligners? readbench (the multithreading audit's harness, built
# against the protal-inexp tree: a212559 + PROTAL_EXPERIMENT_BATCH) with 4 workers that spin a fixed time per
# pair, so that they ask the input path for as many pairs per second as N aligners at 34.7k pairs/s each (the
# full database's rate per thread at 6 threads, ~24 s for 5M pairs): work_us = 28.8 * 4 / N. 4 workers and the
# 2 inflating threads fill the 6 vCPUs. Delivered against demanded pairs/s, for N = 6..192, BGZF / standard gzip /
# uncompressed, batches of 32, 128 and 512 pairs. Rows into ~/mt-work/input/demand.tsv.
set -uo pipefail
SP=$(cd "$(dirname "$0")" && pwd)  # this folder
W=$HOME/mt-work/input; T=$HOME/mt-work/inexp/src
mkdir -p $W/rb2; cp $SP/../readbench.cpp $W/rb2/
g++ -std=c++20 -O3 -march=x86-64-v2 -fopenmp -I$T/src -I$T/src/IO -I$T/src/SequenceUtils -I$T/src/Utilities \
    -I$T/build/zlib-ng -I$T/lib/zlib-ng $W/rb2/readbench.cpp $T/src/IO/FastxReader.cpp \
    $T/build/zlib-ng/libz-ng.a -ldeflate -lpthread -o $W/rb2/readbench || { echo "FAIL readbench build"; exit 1; }
B=$HOME/bench071
declare -A R1=( [bgzf]=$B/samples_deep/points/rl150_p5000000/sim/reads/rl150_p5000000_s_1_R1.fq.gz
                [gzip]=$HOME/mt-audit/plain/R1.std.gz [plain]=$HOME/mt-audit/plain/R1.fq )
declare -A R2=( [bgzf]=$B/samples_deep/points/rl150_p5000000/sim/reads/rl150_p5000000_s_1_R2.fq.gz
                [gzip]=$HOME/mt-audit/plain/R2.std.gz [plain]=$HOME/mt-audit/plain/R2.fq )
quiet() { local i; for i in $(seq 6); do awk -v l=$(cut -d" " -f1 /proc/loadavg) "BEGIN { exit !(l < 6) }" && return; sleep 10; done; }
[ -f $W/demand.tsv ] || printf "input\tbatch\taligners\trep\twork_us\tdemand_pairs_s\tdelivered_pairs_s\tshare\twall_s\tuser_s\tsys_s\tload\n" > $W/demand.tsv
for rep in 1 2; do
  cat ${R1[@]} ${R2[@]} > /dev/null
  for n in 6 12 24 48 96 192; do
    work=$(awk -v n=$n 'BEGIN { printf "%.4f", 28.8 * 4 / n }')
    demand=$(awk -v n=$n 'BEGIN { printf "%.0f", n * 34700 }')
    for batch in 32 128 512; do
      for in in bgzf gzip plain; do
        quiet
        load=$(cut -d' ' -f1 /proc/loadavg)
        line=$(PROTAL_EXPERIMENT_BATCH=$batch nice -n 5 /usr/bin/time -f "%e %U %S" -o $W/rb2.time $W/rb2/readbench ${R1[$in]} ${R2[$in]} 4 $work 2> /dev/null)
        got=$(echo "$line" | awk '{ for (i = 1; i < NF; i++) if ($i == "pairs_per_s") print $(i + 1) }')
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" $in $batch $n $rep $work $demand $got \
          $(awk -v g=$got -v d=$demand 'BEGIN { printf "%.2f", g / d }') "$(tail -1 $W/rb2.time | tr ' ' '\t')" $load | tee -a $W/demand.tsv
      done
    done
  done
done
echo DEMAND DONE
