#!/usr/bin/env bash
# The commit as committed (git archive 0a1beaa): Release build, unit tests, the mini database and the e2e tests. 4 pinned cores.
set -uo pipefail
A=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/8e87f9e0-3682-453c-a7a2-636d4ec2ef61/scratchpad
V=$HOME/vc; rm -rf "$V"; mkdir -p "$V/src"; tar -xf "$A/archives/commit.tar" -C "$V/src"
pin=(taskset -c 0-3 nice -n 5)
start=$(date +%s)
"${pin[@]}" cmake -S "$V/src" -B "$V/build" -G Ninja -DCMAKE_BUILD_TYPE=Release -DPROTAL_BUILD_TESTS=ON > "$V/configure.log" 2>&1 &&
"${pin[@]}" cmake --build "$V/build" -j 4 --target protal protal_tests simulate_metagenomes > "$V/build.log" 2>&1 || { echo "BUILD FAIL"; grep -E "error" -A3 "$V/build.log" | head -30; exit 1; }
echo "build OK in $(( $(date +%s) - start )) s"
bin=$(find "$V/build" -name protal_tests -type f | head -1)
"${pin[@]}" "$bin" --gtest_brief=1 2>&1 | tail -3
export PATH=$HOME/micromamba/envs/protal-db-build/bin:$PATH
cd "$V/src" || exit 1
PROTAL=$V/build/protal "${pin[@]}" bash scripts/mini_db/build_mini_db.sh "$V/mini_db" > "$V/mini_db.log" 2>&1 || { echo "MINI DB FAIL"; tail "$V/mini_db.log"; }
PROTAL_TEST_DB=$V/mini_db/protal_db PROTAL=$V/build/protal SIMULATE=$V/build/simulate_metagenomes "${pin[@]}" python3 -m unittest tests/e2e/test_protal_e2e.py > "$V/e2e.log" 2>&1
echo "e2e exit $?"; tail -3 "$V/e2e.log"; grep -a -E "^(FAIL|ERROR):" "$V/e2e.log" | head
echo "VERIFY DONE"
