#!/usr/bin/env bash
# Items 4a and 4d: rebuild the working tree, all unit tests, then long reads (PacBio 90 Mb, ONT 90 Mb at 6 threads): HEAD,
# the working tree (4a on), the working tree with --long_read_budget 20 and 50; short reads 500k at 1 thread for identity.
set -uo pipefail
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/37e6327e-2743-4622-91c7-90dfd0315ee2/scratchpad
bash $S/build_work.sh '*' 14 || exit 1
echo "== long-read chain test"; $HOME/perf-gtdb/work/build/tests/protal_tests --gtest_filter='AnchoredAlignment.LongReads*:AlignmentScreen.*' 2>&1 | grep -E 'long reads|candidates|PASSED|FAILED'
W=$HOME/perf-gtdb; H=$W/head/build/protal; N=$W/work/build/protal; DB=$HOME/bench071/V073/protal_db
PB=$HOME/bench071/samples_lr073/points/pb_b90000000/sim/reads/pb_b90000000_s_1.fq.gz
ONT=$HOME/bench071/samples_lr073/points/ont_b90000000/sim/reads/ont_b90000000_s_1.fq.gz
P=$HOME/bench071/samples/points; R1=$P/rl150_p500000/sim/reads/rl150_p500000_s_1_R1.fq.gz; R2=$P/rl150_p500000/sim/reads/rl150_p500000_s_1_R2.fq.gz
run() { local name=$1 bin=$2 t=$3 o=$4; shift 4; rm -rf $o
  /usr/bin/time -f "%e s wall, %U s user" $bin --db $DB "$@" --prefix s -o $o -t $t --no_qcmsa > $o.log 2> $o.time
  echo "$name: $(tail -1 $o.time); $(grep -E '^Aligning reads took' $o.log); $(grep -E '^Profiling took' $o.log)"
  grep -E '^Sample s: .*candidate|gene hits fit' $o.log | sed 's/^/    /'; }
cmpset() { if cmp -s <(zstd -dc $1/s.sam.zst | sort) <(zstd -dc $2/s.sam.zst | sort); then echo "    $3: SAM identical as sets"; else echo "    $3: SAM differs: $(diff <(zstd -dc $1/s.sam.zst | sort) <(zstd -dc $2/s.sam.zst | sort) | grep -c '^<') records of $(zstd -dc $1/s.sam.zst | grep -vc '^@') differ"; fi
  diff -q $1/s.profile $2/s.profile > /dev/null && echo "    $3: profile identical" || echo "    $3: profile differs ($(diff $1/s.profile $2/s.profile | grep -c '^<') lines)"; }
for rt in pb ont; do
  R=$PB; [ $rt = ont ] && R=$ONT
  echo "== $rt 90 Mb, 6 threads"
  for round in 1 2; do
    run "$rt head $round" $H 6 $W/l.$rt.head.$round -1 $R --read_type $rt
    run "$rt work $round" $N 6 $W/l.$rt.work.$round -1 $R --read_type $rt
  done
  run "$rt work budget20" $N 6 $W/l.$rt.b20 -1 $R --read_type $rt --long_read_budget 20
  run "$rt work budget50" $N 6 $W/l.$rt.b50 -1 $R --read_type $rt --long_read_budget 50
  cmpset $W/l.$rt.head.1 $W/l.$rt.work.1 "head vs work (4a)"
  cmpset $W/l.$rt.work.1 $W/l.$rt.b20 "work vs budget 20"
  cmpset $W/l.$rt.work.1 $W/l.$rt.b50 "work vs budget 50"
done
echo "== 500k pairs, 1 thread (short reads unchanged by 4a)"
run "pe head" $H 1 $W/s2.head -1 $R1 -2 $R2 --read_type pe
run "pe work" $N 1 $W/s2.work -1 $R1 -2 $R2 --read_type pe
cmp -s <(zstd -dc $W/s2.head/s.sam.zst) <(zstd -dc $W/s2.work/s.sam.zst) && echo "    SAM identical" || echo "    SAM DIFFERS"
diff -r -x '*_runtime.tsv' -x '*.sam.zst' $W/s2.head $W/s2.work > /dev/null && echo "    other outputs identical" || echo "    OUTPUTS DIFFER"
