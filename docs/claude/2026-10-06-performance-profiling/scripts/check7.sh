#!/usr/bin/env bash
# check7.sh: the work build (this session's changes) against the reference build (HEAD): unit tests, whole runs
# (SAM records sorted, every other output file byte for byte, also with PROTAL_FLEX_SCAN=scalar), callgrind of
# aligning 100k pairs and PacBio 3 Mb, measure_performance.sh's new columns, the end-to-end tests on the mini database.
set -uo pipefail
W=$HOME/perf7; R=$W/ref/build; K=$W/work/build
DB=$HOME/bench071/V075/protal_db; P=$HOME/bench071/samples/points; RD=$HOME/perf6/reads
mkdir -p $W/cmp
echo "== unit tests, work build ($(date +%T))"
$K/tests/protal_tests --gtest_brief=1 > $W/unit.work.log 2>&1; echo "exit $?"
grep -E "^\[  (PASSED|FAILED|SKIPPED)|tests ran|FAILED  \]" $W/unit.work.log | head -30
$K/tests/protal_tests --gtest_filter='FlexScan.*:PackedIndex.Lookups*:AlignmentScreen.*:GenomeLoaderPacked.TheFlatTables*' > $W/unit.new.log 2>&1
grep -E "^\[       OK|FAILED|blocks,|screens:|candidates," $W/unit.new.log

echo "== whole runs, 6 threads ($(date +%T))"
run() { local b=$1 o=$2; shift 2; rm -rf $o; $b --db $DB "$@" --prefix s -o $o -t 6 --no_qcmsa --verbose > $o.log 2>&1; echo "    $(basename $o): exit $?"; }
for set in pe se pb ont; do
  case $set in
    pe) args="-1 $P/rl150_p500000/sim/reads/rl150_p500000_s_1_R1.fq.gz -2 $P/rl150_p500000/sim/reads/rl150_p500000_s_1_R2.fq.gz --read_type pe" ;;
    se) args="-1 $P/rl150_p500000/sim/reads/rl150_p500000_s_1_R1.fq.gz --read_type se" ;;
    pb) args="-1 $P/pb_b90000000/sim/reads/pb_b90000000_s_1.fq.gz --read_type pb" ;;
    ont) args="-1 $P/ont_b90000000/sim/reads/ont_b90000000_s_1.fq.gz --read_type ont" ;;
  esac
  echo "-- $set"
  run $R/protal $W/cmp/$set.ref $args
  run $K/protal $W/cmp/$set.work $args
  PROTAL_FLEX_SCAN=scalar run $K/protal $W/cmp/$set.scalar $args
  for v in work scalar; do
    a=$W/cmp/$set.ref; b=$W/cmp/$set.$v
    sa=$(zstd -dc $a/s.sam.zst | grep -v '^@PG\|^@CO' | LC_ALL=C sort | md5sum | cut -c1-12)
    sb=$(zstd -dc $b/s.sam.zst | grep -v '^@PG\|^@CO' | LC_ALL=C sort | md5sum | cut -c1-12)
    n=$(zstd -dc $a/s.sam.zst | grep -vc '^@')
    diffs=$(cd $a && find . -type f ! -name '*.sam.zst' ! -name '*.err' ! -name '*_runtime.tsv' | while read -r f; do cmp -s "$f" "$b/$f" || echo "$f"; done | tr '\n' ' ')
    echo "    ref vs $v: SAM records $([ "$sa" = "$sb" ] && echo "identical ($n)" || echo "DIFFER ($sa vs $sb)"); other files differing: ${diffs:-none}"
  done
  grep -E "seeding:|K-mer screen took|Extending Anchors took|Seeding took" $W/cmp/$set.work.log | cut -c1-200 | sed 's/^/    /'
  grep -E "Extending Anchors took|Seeding took|Alignment handler took" $W/cmp/$set.ref.log | cut -c1-120 | sed 's/^/    ref: /'
done

echo "== callgrind, one thread ($(date +%T))"
cg() { local b=$1 tag=$2; shift 2; rm -rf $W/cmp/cg.$tag
  valgrind --tool=callgrind --callgrind-out-file=$W/cmp/cg.$tag.out $b --db $DB "$@" --prefix s -o $W/cmp/cg.$tag -t 1 --no_qcmsa --no_profile > $W/cmp/cg.$tag.log 2>&1
  callgrind_annotate --inclusive=yes $W/cmp/cg.$tag.out 2>/dev/null > $W/cmp/cg.$tag.incl.txt
  echo "    $tag: $(grep -m1 'PROGRAM TOTALS' $W/cmp/cg.$tag.incl.txt)"; }
for build in ref work; do
  b=$R/protal; [ $build = work ] && b=$K/protal
  cg $b pe.$build -1 $RD/pe100k_R1.fq.gz -2 $RD/pe100k_R2.fq.gz --read_type pe
  cg $b pb.$build -1 $P/pb_b3000000/sim/reads/pb_b3000000_s_1.fq.gz --read_type pb
done
for set in pe pb; do
  for f in 'RunPairedEnd' 'RunLongReads' 'ChainAnchorFinder<protal::KmerLookupSM>::operator()' 'KmerLookupSM::GetFromLookup' 'SimpleAlignmentHandler::AlignAnchor' 'AlignmentScreen::MayAlign'; do
    h=$(grep -m1 -F "$f" $W/cmp/cg.$set.ref.incl.txt | awk '{print $1}'); w=$(grep -m1 -F "$f" $W/cmp/cg.$set.work.incl.txt | awk '{print $1}')
    [ -n "$h$w" ] && echo "    $set $f: ref $h, work $w"
  done
done

echo "== measure_performance.sh ($(date +%T))"
rm -rf $W/mp
REPEATS=1 PERF=0 COHORT=0 THREADS=6 PROTAL=$K/protal bash $W/work/scripts/measure_performance.sh $W/mp $DB \
  pe:$RD/pe100k_R1.fq.gz:$RD/pe100k_R2.fq.gz > $W/mp.log 2>&1; echo "exit $?"
head -1 $W/mp/runs.tsv | tr '\t' '\n' | tail -3 | tr '\n' ' '; echo; tail -1 $W/mp/runs.tsv | tr '\t' '\n' | tail -3 | tr '\n' ' '; echo
grep -E "K-mer screen" $W/mp/stages.tsv | head -3

echo "== end-to-end tests on the mini database ($(date +%T))"
cd $W/work && PROTAL_TEST_DB=$HOME/buildmem/new/data/mini_db/protal_db PROTAL=$K/protal SIMULATE=$K/simulate_metagenomes \
  timeout 3000 python3 -m unittest tests/e2e/test_protal_e2e.py > $W/e2e.log 2>&1; echo "exit $?"
tail -5 $W/e2e.log
echo "CHECKDONE $(date +%T)"
