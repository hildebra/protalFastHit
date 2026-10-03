#!/usr/bin/env bash
set -uo pipefail
W=$HOME/mt-work/perf4; B=$W/ref/build/protal; DB=$HOME/bench071/V073/protal_db
run() { local n=$1 t=$2; shift 2; rm -rf $W/o.$n
  /usr/bin/time -v $B --db $DB "$@" --prefix s -o $W/o.$n -t $t --no_qcmsa --verbose > $W/o.$n.log 2> $W/o.$n.time
  echo "== $n t=$t: $(grep -E 'Elapsed|User time' $W/o.$n.time | tr -s ' ' | tr '\n' ';')"
  grep -E "^(Thread 0 Profile sample|Profiling|Run protal) took" $W/o.$n.log; }
run prof5M_t8 8 --profile_only $W/o.pe5M_t6/s.sam.zst
run prof5M_t12 12 --profile_only $W/o.pe5M_t6/s.sam.zst
run prof5M_t6b 6 --profile_only $W/o.pe5M_t6/s.sam.zst
