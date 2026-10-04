#!/usr/bin/env bash
# After the failed-candidate fold: rebuild + tests (ProfileThreads, ReadEvidence, SampleContext, FalsePositiveFeatures), then
# the deep SAM profiled by HEAD and the new binary at 6 threads (outputs identical?), big2 at 6 and 1 threads (time).
set -uo pipefail
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/37e6327e-2743-4622-91c7-90dfd0315ee2/scratchpad
bash $S/build_work.sh 'ProfileThreads.*:ReadEvidence.*:SampleContext.*:FalsePositiveFeatures.*:AlignmentScreen.*' 6 || exit 1
W=$HOME/perf-gtdb; H=$W/head/build/protal; N=$W/work/build/protal; DB=$HOME/bench071/V073/protal_db; X=$W/unmapped
prof() { local name=$1 bin=$2 t=$3 f=$4 o=$5; rm -rf $o; /usr/bin/time -f "%e s wall, %U s user" $bin --db $DB --profile_only $f --read_type pe --prefix s -o $o -t $t --no_qcmsa > $o.log 2> $o.time
  echo "$name (t=$t): $(tail -1 $o.time); $(grep -E '^Profiling sample' $o.log | cut -c1-110)"; }
prof "deep .sam.zst head" $H 6 $W/d.work.1/s.sam.zst $X/ph
prof "deep .sam.zst work" $N 6 $W/d.work.1/s.sam.zst $X/pw
diff -r -x '*_runtime.tsv' $X/ph $X/pw > /dev/null && echo "    deep profile outputs identical (head vs work, 6 threads)" || { echo "    DEEP OUTPUTS DIFFER"; diff -r -x '*_runtime.tsv' $X/ph $X/pw | head -5; }
prof "deep .sam.zst work" $N 1 $W/d.work.1/s.sam.zst $X/pw1
diff -r -x '*_runtime.tsv' $X/pw $X/pw1 > /dev/null && echo "    deep profile outputs identical (6 vs 1 thread)" || echo "    6 vs 1 THREAD OUTPUTS DIFFER"
prof "big2 head" $H 6 $X/big2.sam $X/bh
prof "big2 work" $N 6 $X/big2.sam $X/bw
diff -r -x '*_runtime.tsv' $X/bh $X/bw > /dev/null && echo "    big2 profile outputs identical (head vs work)" || { echo "    BIG2 OUTPUTS DIFFER"; diff -r -x '*_runtime.tsv' $X/bh $X/bw | head -5; }
prof "big2 work" $N 1 $X/big2.sam $X/bw1
prof "big work" $N 6 $X/big.sam $X/bb
