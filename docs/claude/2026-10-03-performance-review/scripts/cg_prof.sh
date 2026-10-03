#!/usr/bin/env bash
# callgrind of the profiling stage alone (--profile_only on the 500k-pair SAM), one thread; plus a timed plain run.
set -uo pipefail
W=$HOME/mt-work/perf4; B=$W/ref/build/protal; DB=$HOME/bench071/V073/protal_db
SAM=$W/o.pe500k_t1/s.sam.zst
rm -rf $W/p.pe500k $W/cgp.pe500k
/usr/bin/time -v $B --db $DB --profile_only $SAM --prefix s -o $W/p.pe500k -t 1 --no_qcmsa --verbose > $W/p.pe500k.log 2> $W/p.pe500k.time
echo "== profile_only t1: $(grep -E 'Elapsed|User time' $W/p.pe500k.time | tr -s ' ' | tr '\n' ';')"
grep -E "took" $W/p.pe500k.log | tail -12
valgrind --tool=callgrind --callgrind-out-file=$W/cg.prof.out --compress-strings=no --compress-pos=no \
  $B --db $DB --profile_only $SAM --prefix s -o $W/cgp.pe500k -t 1 --no_qcmsa > $W/cg.prof.log 2>&1
callgrind_annotate --inclusive=yes $W/cg.prof.out 2>/dev/null | head -120 > $W/cg.prof.incl.txt
callgrind_annotate $W/cg.prof.out 2>/dev/null | head -120 > $W/cg.prof.self.txt
echo "== callgrind profiling done"
head -40 $W/cg.prof.self.txt
