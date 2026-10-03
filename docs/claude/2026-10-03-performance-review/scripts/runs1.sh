#!/usr/bin/env bash
# Timed whole runs of HEAD (perf4/ref) on the V073 database: stage timers with --verbose, /usr/bin/time.
set -uo pipefail
W=$HOME/mt-work/perf4; B=$W/ref/build/protal; DB=$HOME/bench071/V073/protal_db
P=$HOME/bench071/samples/points; PD=$HOME/bench071/samples_deep/points
run() { # name threads args...
  local n=$1 t=$2; shift 2
  rm -rf $W/o.$n
  /usr/bin/time -v $B --db $DB "$@" --prefix s -o $W/o.$n -t $t --no_qcmsa --verbose > $W/o.$n.log 2> $W/o.$n.time
  echo "== $n t=$t: $(grep -E 'Elapsed|Maximum resident|User time|System time' $W/o.$n.time | tr -s ' ' | tr '\n' ';')"
}
uptime
run pe500k_t1 1 -1 $P/rl150_p500000/sim/reads/rl150_p500000_s_1_R1.fq.gz -2 $P/rl150_p500000/sim/reads/rl150_p500000_s_1_R2.fq.gz --read_type pe
run pe500k_t6 6 -1 $P/rl150_p500000/sim/reads/rl150_p500000_s_1_R1.fq.gz -2 $P/rl150_p500000/sim/reads/rl150_p500000_s_1_R2.fq.gz --read_type pe
run pe5M_t6 6 -1 $PD/rl150_p5000000/sim/reads/rl150_p5000000_s_1_R1.fq.gz -2 $PD/rl150_p5000000/sim/reads/rl150_p5000000_s_1_R2.fq.gz --read_type pe
run ont90M_t6 6 -1 $P/ont_b90000000/sim/reads/ont_b90000000_s_1.fq.gz --read_type ont
run pb90M_t6 6 -1 $P/pb_b90000000/sim/reads/pb_b90000000_s_1.fq.gz --read_type pb
uptime
