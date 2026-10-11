#!/usr/bin/env bash
# identity.sh: whole protal runs of the new build (4 threads) against the base's: classic against the base's dumped
# runs (1 thread, setup.sh), Singletrack against base Singletrack runs (4 threads); SAM records and every other output
# file compared. Then the e2e tests of the new build on the mini database. On cores 0-3.
set -uo pipefail
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/6a0fc5a0-51ae-4ff5-a1a7-da994d3168e8/scratchpad/wfaspeed
W=$HOME/wfaspeed; B=$W/base/build; N=$W/new/build; D=$W/data; O=$W/out
DB=$HOME/realworld_s1/protal_db
pin() { taskset -c 0-3 nice -n 5 "$@"; }
args_of() {
  case $1 in
    pe) echo "-1 $D/pe/reads/pe_s_1_R1.fq.gz -2 $D/pe/reads/pe_s_1_R2.fq.gz --read_type pe" ;;
    pb) echo "-1 $D/pb.fq.gz --read_type pb" ;;
    ont) echo "-1 $D/ont.fq.gz --read_type ont" ;;
  esac
}
run() { local b=$1 o=$2 set=$3; shift 3; rm -rf $o
  pin /usr/bin/time -f "TIME %e s, max RSS %M KB" $b --db $DB $(args_of $set) "$@" --prefix s -o $o -t 4 --no_qcmsa > $o.log 2>&1
  echo "    $(basename $o): exit $?, $(grep -m1 '^TIME' $o.log)"; }
compare() { local a=$1 b=$2
  local diffs=$(cd $a && find . -type f ! -name '*.sam.zst' ! -name '*.err' ! -name '*_runtime.tsv' ! -name '*.log' | while read -r f; do cmp -s "$f" "$b/$f" || echo "$f"; done | tr '\n' ' ')
  echo "    $(basename $a) vs $(basename $b): $(cd $a && find . -type f | wc -l) files, other files differing: ${diffs:-none}"
  python3 $S/compare_sams.py $a/s.sam.zst $b/s.sam.zst --examples 2 | sed 's/^/      /'; }
for set in pe pb ont; do
  echo "-- $set ($(date +%T))"
  run $N/protal $O/$set.new $set
  compare $O/$set.base $O/$set.new
  run $B/protal $O/$set.base_single $set --wfa_traceback singletrack
  run $N/protal $O/$set.new_single $set --wfa_traceback singletrack
  compare $O/$set.base_single $O/$set.new_single
done
echo "== e2e ($(date +%T))"
cd $W/new
PROTAL=build/protal pin bash scripts/mini_db/build_mini_db.sh data/mini_db > $W/mini_db.log 2>&1; echo "mini db exit $?"
PROTAL_TEST_DB=data/mini_db/protal_db PROTAL=build/protal SIMULATE=build/simulate_metagenomes PROTAL_TESTS_REQUIRED=1 \
  pin python3 -m unittest tests/e2e/test_protal_e2e.py > $W/e2e.log 2>&1; echo "e2e exit $?"
tail -4 $W/e2e.log
echo "IDENTITYDONE $(date +%T)"
