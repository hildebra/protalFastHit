#!/usr/bin/env bash
# All unit tests of the packed-index build, then the end-to-end tests on a fresh mini database
# (build_packed.sh first). Output in WORK/verify.log.
SRC=${SRC:-/mnt/c/Users/hildebra/Documents/locDev/protal}
WORK=${WORK:-$HOME/protal-pack}
cd "$WORK" || exit 1
{
  echo "== unit tests =="
  ctest --test-dir build --output-on-failure 2>&1 | tail -8
  echo "== mini database =="
  rm -rf mini && PROTAL="$WORK/build/protal" bash src/scripts/mini_db/build_mini_db.sh mini > mini.log 2>&1 && echo "mini db built" || { echo "mini db failed"; tail -20 mini.log; }
  grep -E "Index in memory|Load index" mini/build.log 2>/dev/null | head -3
  echo "== e2e =="
  cd src && PROTAL_TEST_DB="$WORK/mini/protal_db" PROTAL="$WORK/build/protal" SIMULATE="$WORK/build/simulate_metagenomes" \
    python3 -m unittest tests/e2e/test_protal_e2e.py 2>&1 | tail -6
} > verify.log 2>&1
cat "$WORK/verify.log"
