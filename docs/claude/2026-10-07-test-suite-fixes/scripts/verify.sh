#!/bin/bash
# Every suite as CI runs it, on the content of one commit (git archive), on CPUs 0-3, PROTAL_TESTS_REQUIRED=1.
# Usage: verify.sh COMMIT
set -o pipefail
COMMIT=$1
REPO=/mnt/c/Users/hildebra/Documents/locDev/protal-testaudit
V=~/ti/verify
rm -rf "$V" && mkdir -p "$V/archive"
tar -x -C "$V/archive" -f /mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/44100026-0bec-4c69-b174-7a4fb9e06969/scratchpad/verify.tar || exit 1  # git archive COMMIT, made on Windows
# Into the existing build tree (an incremental build of exactly this content).
rsync -a --checksum --no-times --delete --exclude=/build*/ --exclude=/data "$V/archive/" ~/testaudit/src/
SRC=~/testaudit/src
cd "$SRC" || exit 1
export PROTAL_TESTS_REQUIRED=1
export PATH=~/ta/venv14/bin:$PATH   # python3 with numpy, pandas, scikit-learn 1.4.1 and joblib, as Ubuntu 24.04's
T="taskset -c 0-3 nice -n 5"
step() {
  local name=$1; shift
  local start=$(date +%s)
  "$@" > "$V/$name.log" 2>&1
  local rc=$?
  echo "$name: exit $rc, $(( $(date +%s) - start )) s"
  [ $rc -ne 0 ] && tail -40 "$V/$name.log"
  return $rc
}
echo "commit $COMMIT, python3 $(python3 -c 'import sklearn; print("sklearn", sklearn.__version__)')"
step build $T ninja -C build -j4 protal simulate_metagenomes protal_tests || exit 1
step unit $T ctest --test-dir build --output-on-failure -j4
grep -E 'tests passed|tests failed' "$V/unit.log"
step mini_db env PROTAL=$SRC/build/protal $T bash scripts/mini_db/build_mini_db.sh "$V/data/mini_db"
step e2e env PROTAL_TEST_DB="$V/data/mini_db/protal_db" PROTAL=$SRC/build/protal SIMULATE=$SRC/build/simulate_metagenomes \
  $T python3 -m unittest -v tests/e2e/test_protal_e2e.py
tail -4 "$V/e2e.log"
step example env PROTAL=$SRC/build/protal THREADS=4 $T bash examples/mini_db/run.sh "$V/example"
grep -E 'PASS|FAIL' "$V/example.log" | tail -2
step mini_db_tests env PROTAL=$SRC/build/protal SIMULATE=$SRC/build/simulate_metagenomes \
  $T bash -c 'python3 -m unittest -v scripts/mini_db/test_*.py'
tail -4 "$V/mini_db_tests.log"
step script_tests $T python3 -m unittest -v scripts/test_insilico_strains.py scripts/test_trace_relatives.py \
  scripts/test_error_reads.py scripts/test_profile_scripts.py scripts/test_strain_scripts.py scripts/test_model_pmml.py
tail -4 "$V/script_tests.log"
grep -h -E '^Ran |skipped' "$V"/e2e.log "$V"/mini_db_tests.log "$V"/script_tests.log
echo done
