#!/usr/bin/env bash
# stk_e2e_fresh.sh: a fresh mini database built by the work tree (scripts/mini_db/build_mini_db.sh, as just mini-db), then
# the whole end-to-end suite of the work build on it, and the same suite of the reference build on the same database.
set -uo pipefail
C=$HOME/stk/cmp; M=$HOME/stk/minidb
cd $HOME/stk/work
rm -rf $M
PROTAL=$HOME/stk/work/build/protal taskset -c 0-3 nice -n 5 bash scripts/mini_db/build_mini_db.sh $M > $C/minidb.log 2>&1; echo "mini db exit $?"
tail -3 $C/minidb.log; ls $M/protal_db | head -30
for t in work ref; do
  B=$HOME/stk/$t/build
  cd $HOME/stk/$t
  PROTAL_TEST_DB=$M/protal_db PROTAL=$B/protal SIMULATE=$B/simulate_metagenomes taskset -c 0-3 nice -n 5 timeout 3000 \
    python3 -m unittest tests/e2e/test_protal_e2e.py > $C/e2e.fresh.$t.log 2>&1; echo "$t suite exit $?"
  tail -3 $C/e2e.fresh.$t.log
  grep -E '^(ERROR|FAIL):' $C/e2e.fresh.$t.log | sort > $C/e2e.fresh.$t.fails
done
echo "failing in work only:"; comm -13 $C/e2e.fresh.ref.fails $C/e2e.fresh.work.fails
echo "failing in ref only:"; comm -23 $C/e2e.fresh.ref.fails $C/e2e.fresh.work.fails
echo "E2EFRESHDONE $(date +%T)"
