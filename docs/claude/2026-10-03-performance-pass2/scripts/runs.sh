#!/usr/bin/env bash
# Whole runs of HEAD on the v0.7.3 world: 500k pairs at 1 and 6 threads, with stage timers.
set -uo pipefail
W=$HOME/perf-pass2; B=$W/head/build/protal; DB=$HOME/bench071/V073/protal_db
P=$HOME/bench071/samples/points
R1=$P/rl150_p500000/sim/reads/rl150_p500000_s_1_R1.fq.gz; R2=$P/rl150_p500000/sim/reads/rl150_p500000_s_1_R2.fq.gz
for t in 1 6; do
  rm -rf $W/o.pe500k_t$t
  /usr/bin/time -f "%e s wall, %U s user, %M KB" $B --db $DB -1 $R1 -2 $R2 --read_type pe --prefix s -o $W/o.pe500k_t$t -t $t --no_qcmsa --verbose > $W/o.pe500k_t$t.log 2> $W/o.pe500k_t$t.time
  echo "t=$t: $(tail -1 $W/o.pe500k_t$t.time)"
done
# 100k-pair subset for callgrind
mkdir -p $W/reads
zcat $R1 | head -400000 | gzip -1 > $W/reads/pe100k_R1.fq.gz
zcat $R2 | head -400000 | gzip -1 > $W/reads/pe100k_R2.fq.gz
