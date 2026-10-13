#!/bin/bash
# protal_threads.sh at -t 4 only, 4 alternated rounds: 1, 2 and 4 reading threads per file (results/protal_threads_t4.txt).
W=~/pgzfix; cd $W/threads
for round in 1 2 3 4; do for kind in one members; do for i in 1 2 4; do
  out=out_x; rm -rf $out
  PROTAL_INFLATE_THREADS=$i taskset -c 0-3 $W/build/protal --db $W/mini_db/protal_db -t 4 -1 ${kind}_R1.fq.gz -2 ${kind}_R2.fq.gz --prefix s -o $out --no_strains > $out.log 2>&1
  a=$(grep -m1 '^Aligning reads took' $out.log | sed 's/Aligning reads took //'); sh=$(grep -m1 '^Sample s input: ' $out.log | grep -oE '[0-9]+% of')
  printf 'r%s\t%s\t-t 4\tinflate %s\taligning %s\treader %s\n' $round $kind $i "$a" "$sh"
done; done; done
