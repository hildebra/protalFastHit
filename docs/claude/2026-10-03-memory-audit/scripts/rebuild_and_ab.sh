#!/usr/bin/env bash
# Rebuilds the packed-index tree (build_packed.sh), reruns the end-to-end database tests, then the A/B
# against the baseline (ab_packed.sh).
HERE=$(cd "$(dirname "$0")" && pwd)
WORK=${WORK:-$HOME/protal-pack}
bash "$HERE/build_packed.sh" | tail -3 || exit 1
cd "$WORK/src" || exit 1
echo "== e2e: database tests =="
PROTAL_TEST_DB="$WORK/mini/protal_db" PROTAL="$WORK/build/protal" SIMULATE="$WORK/build/simulate_metagenomes" \
  python3 -m unittest tests.e2e.test_protal_e2e.CompressedDatabaseTest 2>&1 | tail -4
echo "== A/B =="
bash "$HERE/ab_packed.sh" 6 1
