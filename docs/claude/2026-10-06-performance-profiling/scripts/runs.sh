#!/usr/bin/env bash
# Whole runs of HEAD on the v0.7.5 world with stage timers: pe 500k, se 500k (R1), PacBio and ONT 90 Mb; 1 and 6 threads.
# Usage: runs.sh <round>   (writes ~/perf6/o.<set>_t<t>.r<round>.{log,time})
set -uo pipefail
R=${1:-1}
W=$HOME/perf6; B=$W/head/build/protal; DB=$HOME/bench071/V075/protal_db
P=$HOME/bench071/samples/points
PE1=$P/rl150_p500000/sim/reads/rl150_p500000_s_1_R1.fq.gz; PE2=$P/rl150_p500000/sim/reads/rl150_p500000_s_1_R2.fq.gz
PB=$P/pb_b90000000/sim/reads/pb_b90000000_s_1.fq.gz
ONT=$P/ont_b90000000/sim/reads/ont_b90000000_s_1.fq.gz
run() { # name threads args...
  local n=$1 t=$2; shift 2
  local o=$W/o.${n}_t$t.r$R
  rm -rf $o
  /usr/bin/time -v $B --db $DB "$@" --prefix s -o $o -t $t --no_qcmsa --verbose > $o.log 2> $o.time
  echo "$n t=$t r=$R: $(grep -E 'Elapsed' $o.time | awk '{print $NF}') wall, $(grep 'User time' $o.time | awk '{print $NF}') user, $(grep 'Maximum resident' $o.time | awk '{print $NF}') KB; load $(cut -d' ' -f1 /proc/loadavg)"
}
for t in ${THREADS:-6 1}; do
  run pe $t -1 $PE1 -2 $PE2 --read_type pe
  run se $t -1 $PE1 --read_type se
  run pb $t -1 $PB --read_type pb
  run ont $t -1 $ONT --read_type ont
done
