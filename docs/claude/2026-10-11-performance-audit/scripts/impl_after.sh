#!/usr/bin/env bash
# After the timed runs: the mini-database and GTDB-build Python tests on the changed protal, the unit tests under
# ASan + UBSan (as CI's sanitizer job), and measure_performance.sh's new columns on one pe run. 4 pinned cores.
set -uo pipefail
A=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/8e87f9e0-3682-453c-a7a2-636d4ec2ef61/scratchpad
REPO=/mnt/c/Users/hildebra/Documents/locDev/protal
W=$HOME/impl; S=$W/work/src; B=$W/work/build
until grep -q "ALL DONE" "$A/impl_perf_queue.out" 2> /dev/null; do sleep 15; done
export PATH=$HOME/micromamba/envs/protal-db-build/bin:$PATH
pin=(taskset -c 0-3 nice -n 5)

echo "[$(date +%T)] rebuild (the deferred tables' timer) and the unit tests"
bash "$A/impl_build.sh" 2>&1 | tail -3
bash "$A/impl_unit.sh" 2>&1 | tail -6
echo "[$(date +%T)] the timed runs again (runs2)"
OUT=$HOME/impl/perf/runs2 bash "$A/impl_perf_runs.sh" 2>&1 | tail -40

echo "[$(date +%T)] e2e tests"
start=$(date +%s)
cd "$S" || exit 1
rsync -a --checksum --no-times "$REPO/tests/e2e/test_protal_e2e.py" "$S/tests/e2e/test_protal_e2e.py"
PROTAL_TEST_DB=$W/mini_db/protal_db PROTAL=$B/protal SIMULATE=$B/simulate_metagenomes "${pin[@]}" python3 -m unittest tests/e2e/test_protal_e2e.py > "$W/e2e.log" 2>&1
echo "e2e exit $? in $(( $(date +%s) - start )) s"
tail -3 "$W/e2e.log"
grep -E "^(FAIL|ERROR):" "$W/e2e.log" | head -20

echo "[$(date +%T)] mini_db tests"
cd "$S" || exit 1
start=$(date +%s)
PROTAL=$B/protal SIMULATE=$B/simulate_metagenomes "${pin[@]}" python3 -m unittest -v scripts/mini_db/test_*.py > "$W/minidb_tests.log" 2>&1
echo "mini_db tests exit $? in $(( $(date +%s) - start )) s"
tail -4 "$W/minidb_tests.log"
grep -E "^(FAIL|ERROR):" "$W/minidb_tests.log" | head -20

echo "[$(date +%T)] ASan unit tests"
FLAGS="-fsanitize=address,undefined -fno-omit-frame-pointer"
start=$(date +%s)
if "${pin[@]}" cmake -S "$S" -B "$W/asan" -G Ninja -DCMAKE_BUILD_TYPE=Debug -DPROTAL_BUILD_TESTS=ON \
       -DCMAKE_CXX_FLAGS="$FLAGS -DPROTAL_NO_CLONES" -DCMAKE_C_FLAGS="$FLAGS" -DCMAKE_EXE_LINKER_FLAGS="$FLAGS" > "$W/asan_configure.log" 2>&1 &&
   "${pin[@]}" cmake --build "$W/asan" -j 4 --target protal_tests > "$W/asan_build.log" 2>&1; then
  echo "asan build in $(( $(date +%s) - start )) s"
  start=$(date +%s)
  ASAN_OPTIONS=detect_leaks=1 UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=1 "${pin[@]}" ctest --test-dir "$W/asan" --output-on-failure -j4 > "$W/asan_tests.log" 2>&1
  echo "asan ctest exit $? in $(( $(date +%s) - start )) s"
  tail -6 "$W/asan_tests.log"
  grep -E "ERROR: (Address|Leak)Sanitizer|runtime error|FAILED" "$W/asan_tests.log" | head -20
else
  echo "asan build FAIL"; grep -E "error" -A4 "$W/asan_build.log" | head -40
fi

echo "[$(date +%T)] measure_performance.sh"
P=$W/perf; rm -rf "$P/measure"
PROTAL=$B/protal THREADS=4 REPEATS=1 COHORT=0 PERF=0 "${pin[@]}" bash "$REPO/scripts/measure_performance.sh" "$P/measure" "$P/db_bin/database.protal" \
    "pe:$P/samples/pe1m/reads/pe1m_1_R1.fq.gz:$P/samples/pe1m/reads/pe1m_1_R2.fq.gz" "pb:$P/samples/long/pb.fq.gz" > "$W/measure.log" 2>&1
echo "measure exit $?"
tail -8 "$W/measure.log"
awk -F'\t' 'NR == 1 { n = NF } { print NF ": " $1 " " $(NF - 2) " " $(NF - 1) " " $NF }' "$P/measure/runs.tsv"
echo "[$(date +%T)] AFTER DONE"
