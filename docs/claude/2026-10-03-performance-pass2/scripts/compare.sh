#!/usr/bin/env bash
# head (0e4c5a4) against the working tree: SAM text of 500k pairs (1 thread), the profile outputs for the same SAM,
# stage timers, then callgrind of aligning 100k pairs and profiling the SAM.
set -uo pipefail
W=$HOME/perf-pass2; H=$W/head/build/protal; N=$W/work/build/protal; DB=$HOME/bench071/V073/protal_db
P=$HOME/bench071/samples/points
R1=$P/rl150_p500000/sim/reads/rl150_p500000_s_1_R1.fq.gz; R2=$P/rl150_p500000/sim/reads/rl150_p500000_s_1_R2.fq.gz
for round in 1 2; do
 for v in head work; do
  B=$H; [ $v = work ] && B=$N
  rm -rf $W/a.$v.$round
  /usr/bin/time -f "%e s wall, %U s user" $B --db $DB -1 $R1 -2 $R2 --read_type pe --no_profile --prefix s -o $W/a.$v.$round -t 1 --no_qcmsa --verbose > $W/a.$v.$round.log 2> $W/a.$v.$round.time
  echo "$v round $round: $(tail -1 $W/a.$v.$round.time); $(grep -E '^Aligning reads took' $W/a.$v.$round.log); $(grep -E '^\s*Raw alignment took' $W/a.$v.$round.log)"
 done
done
if cmp -s <(zstd -dc $W/a.head.1/s.sam.zst) <(zstd -dc $W/a.work.1/s.sam.zst); then echo "== SAM text identical ($(zstd -dc $W/a.work.1/s.sam.zst | grep -vc '^@') records)"; else echo "== SAM DIFFERS"; diff <(zstd -dc $W/a.head.1/s.sam.zst) <(zstd -dc $W/a.work.1/s.sam.zst) | head -10; fi
SAM=$W/a.head.1/s.sam.zst
for round in 1 2; do
 for v in head work; do
  B=$H; [ $v = work ] && B=$N
  rm -rf $W/c.$v.$round
  /usr/bin/time -f "%e s wall, %U s user" $B --db $DB --profile_only $SAM --prefix s -o $W/c.$v.$round -t 1 --no_qcmsa --verbose > $W/c.$v.$round.log 2> $W/c.$v.$round.time
  echo "profile $v round $round: $(tail -1 $W/c.$v.$round.time); $(grep -E '^Profiling took' $W/c.$v.$round.log)"
 done
done
if diff -r $W/c.head.1 $W/c.work.1 > $W/c.diff; then echo "== profile outputs identical"; else echo "== PROFILE OUTPUTS DIFFER"; head -10 $W/c.diff; fi
