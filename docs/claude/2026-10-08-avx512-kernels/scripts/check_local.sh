#!/usr/bin/env bash
# check.sh: the work build (AVX-512 kernels, the reworked syncmer fill) against the reference build (HEAD): unit tests,
# whole runs (SAM records sorted, every other output file byte for byte; also PROTAL_SIMD=scalar), callgrind of
# aligning 100k pairs and PacBio 3 Mb. 4 cores, niced.
set -uo pipefail
W=$HOME/avx512; R=$W/ref/build; K=$W/work/build
DB=$HOME/bench071/V075/protal_db; P=$HOME/bench071/samples/points; RD=$HOME/perf6/reads
T="taskset -c 0-3 nice -n 5"
mkdir -p $W/cmp
ls -d $DB $P $RD > /dev/null || exit 1
echo "== unit tests, work build ($(date +%T))"
$T $K/tests/protal_tests --gtest_brief=1 > $W/unit.work.log 2>&1; echo "exit $?"
grep -E "^\[  (PASSED|FAILED|SKIPPED)|tests ran|FAILED  \]" $W/unit.work.log | head -30
$T $K/tests/protal_tests --gtest_filter='FlexScan.*:Syncmers.*:PackedIndex.*:PackedSequence.*:AlignmentScreen.*' > $W/unit.new.log 2>&1
echo "kernel suites exit $?"; grep -E "^\[       OK|FAILED|this CPU|no AVX|blocks, " $W/unit.new.log
echo "-- the same suites with PROTAL_SIMD=scalar (default levels scalar)"
PROTAL_SIMD=scalar $T $K/tests/protal_tests --gtest_filter='FlexScan.*:Syncmers.*:PackedIndex.*:PackedSequence.*' > $W/unit.scalar.log 2>&1
echo "exit $?"; grep -E "FAILED  \]|^\[  PASSED" $W/unit.scalar.log

echo "== whole runs, 4 threads ($(date +%T))"
run() { local b=$1 o=$2; shift 2; rm -rf $o; $T $b --db $DB "$@" --prefix s -o $o -t 4 --no_qcmsa --verbose > $o.log 2>&1; echo "    $(basename $o): exit $?"; }
for set in pe pb ont; do
  case $set in
    pe) args="-1 $P/rl150_p500000/sim/reads/rl150_p500000_s_1_R1.fq.gz -2 $P/rl150_p500000/sim/reads/rl150_p500000_s_1_R2.fq.gz --read_type pe" ;;
    pb) args="-1 $P/pb_b90000000/sim/reads/pb_b90000000_s_1.fq.gz --read_type pb" ;;
    ont) args="-1 $P/ont_b90000000/sim/reads/ont_b90000000_s_1.fq.gz --read_type ont" ;;
  esac
  echo "-- $set"
  run $R/protal $W/cmp/$set.ref $args
  run $K/protal $W/cmp/$set.work $args
  PROTAL_SIMD=scalar run $K/protal $W/cmp/$set.scalar $args
  a=$W/cmp/$set.ref
  sa=$(zstd -dc $a/s.sam.zst | grep -v '^@PG\|^@CO' | LC_ALL=C sort | md5sum | cut -c1-12)
  n=$(zstd -dc $a/s.sam.zst | grep -vc '^@')
  for v in work scalar; do
    b=$W/cmp/$set.$v
    sb=$(zstd -dc $b/s.sam.zst | grep -v '^@PG\|^@CO' | LC_ALL=C sort | md5sum | cut -c1-12)
    diffs=$(cd $a && find . -type f ! -name '*.sam.zst' ! -name '*.err' ! -name '*_runtime.tsv' | while read -r f; do cmp -s "$f" "$b/$f" || echo "$f"; done | tr '\n' ' ')
    echo "    ref vs $v: SAM records $([ "$sa" = "$sb" ] && echo "identical ($n)" || echo "DIFFER ($sa vs $sb)"); other files differing: ${diffs:-none}"
  done
  for v in ref work scalar; do
    grep -E "Kmer extraction took|K-mer extraction took|Seeding took|Taking the k-mers" $W/cmp/$set.$v.log | cut -c1-110 | sed "s/^/    $v: /"
  done
done

echo "== callgrind, one thread ($(date +%T))"
cg() { local b=$1 tag=$2; shift 2; rm -rf $W/cmp/cg.$tag
  $T valgrind --tool=callgrind --callgrind-out-file=$W/cmp/cg.$tag.out $b --db $DB "$@" --prefix s -o $W/cmp/cg.$tag -t 1 --no_qcmsa --no_profile > $W/cmp/cg.$tag.log 2>&1
  callgrind_annotate --inclusive=yes $W/cmp/cg.$tag.out 2>/dev/null > $W/cmp/cg.$tag.incl.txt
  echo "    $tag: $(grep -m1 'PROGRAM TOTALS' $W/cmp/cg.$tag.incl.txt)"; }
for build in ref work; do
  b=$R/protal; [ $build = work ] && b=$K/protal
  cg $b pe.$build -1 $RD/pe100k_R1.fq.gz -2 $RD/pe100k_R2.fq.gz --read_type pe
  cg $b pb.$build -1 $P/pb_b3000000/sim/reads/pb_b3000000_s_1.fq.gz --read_type pb
done
for set in pe pb; do
  for f in 'RunPairedEnd' 'RunLongReads' 'ScanClosedSyncmers' 'ScanWindowsAvx2' 'FillAvx2' 'GetFromLookup'; do
    h=$(grep -F "$f" $W/cmp/cg.$set.ref.incl.txt | head -1 | awk '{print $1}'); w=$(grep -F "$f" $W/cmp/cg.$set.work.incl.txt | head -1 | awk '{print $1}')
    [ -n "$h$w" ] && echo "    $set $f: ref ${h:--}, work ${w:--}"
  done
done
echo "== done ($(date +%T))"
