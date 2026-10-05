#!/usr/bin/env bash
# Profiles the synthetic SAMs (from their .sam.zst, no alignment) at the given thread counts; prints the profiling timers.
set -uo pipefail
W=$HOME/perf-gtdb; X=$W/samread; N=${PROTAL:-$W/work/build/protal}; DB=$W/e2e/db/database.protal; TAG=${TAG:-base}
P=$HOME/bench071/samples/points; R1=$P/rl150_p500000/sim/reads/rl150_p500000_s_1_R1.fq.gz; R2=$P/rl150_p500000/sim/reads/rl150_p500000_s_1_R2.fq.gz
for f in ${FILES:-mapped mixed}; do for t in ${THREADS:-6 1}; do
  o=$X/out.$TAG.$f.t$t; rm -rf $o; mkdir -p $o; map=$o.map
  { printf '#OUTPUT_DIR\t%s\n' "$o"; printf '#SAMPLEID\tFIRST\tSECOND\tSAM\tPREFIX\tREAD_TYPE\n'; printf 's\t%s\t%s\t%s\ts\tpe\n' "$R1" "$R2" "$X/$f.sam.zst"; } > $map
  /usr/bin/time -f "%e s wall %M KB" $N --db $DB --map $map -t $t --no_strains > $o.log 2> $o.time
  echo "$TAG $f t$t: $(grep -h '^Profiling sample' $o.log | cut -c1-160) | $(tail -1 $o.time)"
done; done
