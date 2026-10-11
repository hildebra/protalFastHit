#!/usr/bin/env bash
# stk_check.sh unit|identity|timing|callgrind|e2e: the Singletrack checks in ~/stk, on cores 0-3 (niced).
#  unit:      the work build's unit tests (and the new ones' counts of other operations)
#  identity:  pe, pb and ont runs (4 threads) of ref (WFA2-lib v2.3.6), work (fork, classic: the database has no
#             build_settings.tsv) and work --wfa_traceback singletrack; SAM records and every other file compared
#  timing:    alternated runs of work classic and singletrack, 4 threads, 3 rounds: stage timers and peak RSS
#  callgrind: one thread on subsets (pe 50k pairs, pb and ont ~3 Mb): instructions in total and in WFA2
set -uo pipefail
MODE=$1
W=$HOME/stk; R=$W/ref/build; K=$W/work/build; D=$W/data; C=$W/cmp
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/6a0fc5a0-51ae-4ff5-a1a7-da994d3168e8/scratchpad
DB=$HOME/realworld_s1/protal_db
mkdir -p $C
pin() { taskset -c 0-3 nice -n 5 "$@"; }
args_of() {
  case $1 in
    pe) echo "-1 $D/pe/reads/pe_s_1_R1.fq.gz -2 $D/pe/reads/pe_s_1_R2.fq.gz --read_type pe" ;;
    pb) echo "-1 $D/long/pb.fq.gz --read_type pb" ;;
    ont) echo "-1 $D/long/ont.fq.gz --read_type ont --no_profile" ;;
    pe50k) echo "-1 $D/sub/pe50k_R1.fq.gz -2 $D/sub/pe50k_R2.fq.gz --read_type pe --no_profile" ;;
    pb3) echo "-1 $D/sub/pb3.fq.gz --read_type pb --no_profile" ;;
    ont3) echo "-1 $D/sub/ont3.fq.gz --read_type ont --no_profile" ;;
    pe3m) echo "-1 $D/pe3m/reads/pe3m_s_1_R1.fq.gz -2 $D/pe3m/reads/pe3m_s_1_R2.fq.gz --read_type pe" ;;
    pbbig) echo "-1 $D/big/pb.fq.gz --read_type pb" ;;
    ontbig) echo "-1 $D/big/ont.fq.gz --read_type ont --no_profile" ;;
  esac
}
run() { local b=$1 o=$2 set=$3; shift 3; rm -rf $o
  pin /usr/bin/time -f "TIME %e s, max RSS %M KB" $b --db $DB $(args_of $set) "$@" --prefix s -o $o -t 4 --no_qcmsa --verbose > $o.log 2>&1
  local rc=$?
  echo "    $(basename $o): exit $rc, $(grep -m1 '^TIME' $o.log)"; }
SETS=${SETS:-pe pb ont}

case $MODE in
unit)
  echo "== unit tests, work build ($(date +%T))"
  pin $K/tests/protal_tests --gtest_brief=1 > $C/unit.work.log 2>&1; echo "exit $?"
  grep -E "^\[  (PASSED|FAILED|SKIPPED)|tests ran|FAILED  \]" $C/unit.work.log | head -30
  pin $K/tests/protal_tests --gtest_filter='WfaTraceback.*:BuildSettings.*:WFA2Wrapper.*' --gtest_output=xml:$C/unit.new.xml > $C/unit.new.log 2>&1
  grep -E "^\[       OK|FAILED" $C/unit.new.log
  grep -o '<property name="[a-z0-9_]*other_operations" value="[^"]*"' $C/unit.new.xml
  ;;
identity)
  for set in $SETS; do
    echo "-- $set ($(date +%T))"
    run $R/protal $C/$set.ref $set
    run $K/protal $C/$set.classic $set
    run $K/protal $C/$set.single $set --wfa_traceback singletrack
    grep -h "WFA2 traceback" $C/$set.classic.log $C/$set.single.log | sed 's/^/    /'
    grep -c "Warning: this run aligns with" $C/$set.classic.log $C/$set.single.log | sed 's/^/    warnings: /'
    for pair in "ref classic" "classic single"; do
      set -- $pair; a=$C/$set.$1; b=$C/$set.$2
      diffs=$(cd $a && find . -type f ! -name '*.sam.zst' ! -name '*.err' ! -name '*_runtime.tsv' ! -name '*.log' | while read -r f; do cmp -s "$f" "$b/$f" || echo "$f"; done | tr '\n' ' ')
      echo "    $1 vs $2: other files differing: ${diffs:-none}"
      python3 $S/compare_sams.py $a/s.sam.zst $b/s.sam.zst --examples 2 | sed 's/^/      /'
    done
  done
  ;;
timing)
  for set in $SETS; do
    echo "-- $set ($(date +%T))"
    for round in 1 2 3 4 5; do
      for v in classic single; do
        extra=""; [ $v = single ] && extra="--wfa_traceback singletrack"
        run $K/protal $C/t.$set.$v.$round $set --no_profile $extra
        grep -E "Processing all samples took|Alignment handler took|Long-read alignment took|Aligning took" $C/t.$set.$v.$round.log | head -3 | cut -c1-110 | sed 's/^/      /'
      done
    done
  done
  ;;
callgrind)
  mkdir -p $D/sub
  [ -f $D/sub/pe50k_R1.fq.gz ] || { zcat $D/pe/reads/pe_s_1_R1.fq.gz | head -200000 | gzip -1 > $D/sub/pe50k_R1.fq.gz
                                    zcat $D/pe/reads/pe_s_1_R2.fq.gz | head -200000 | gzip -1 > $D/sub/pe50k_R2.fq.gz
                                    zcat $D/long/pb.fq.gz | head -1000 | gzip -1 > $D/sub/pb3.fq.gz
                                    zcat $D/long/ont.fq.gz | head -1600 | gzip -1 > $D/sub/ont3.fq.gz; }
  for set in pe50k pb3 ont3; do
    for v in classic single; do
      extra=""; [ $v = single ] && extra="--wfa_traceback singletrack"
      o=$C/cg.$set.$v; rm -rf $o
      pin valgrind --tool=callgrind --cache-sim=yes --callgrind-out-file=$o.out $K/protal --db $DB $(args_of $set) $extra --prefix s -o $o -t 1 --no_qcmsa > $o.log 2>&1
      callgrind_annotate --inclusive=yes $o.out 2>/dev/null > $o.incl.txt
      # events: Ir Dr Dw I1mr D1mr D1mw ILmr DLmr DLmw
      echo "    $set $v:"
      grep -m1 'PROGRAM TOTALS' $o.incl.txt | sed 's/^/      /'
      grep -m1 -E ':wavefront_align ' $o.incl.txt | cut -c1-230 | sed 's/^/      /'
      grep -m1 -E ':wavefront_unialign ' $o.incl.txt | cut -c1-230 | sed 's/^/      /'
      grep -m1 -E ':wavefront_backtrace_affine(_m_only)? ' $o.incl.txt | cut -c1-230 | sed 's/^/      /'
    done
  done
  ;;
esac
echo "CHECKDONE $MODE $(date +%T)"
