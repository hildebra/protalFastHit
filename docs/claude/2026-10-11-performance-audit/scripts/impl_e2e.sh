#!/usr/bin/env bash
# The mini database built by the changed protal, then the end-to-end tests on it; 4 pinned cores.
set -uo pipefail
W=$HOME/impl; S=$W/work/src; B=$W/work/build
export PATH=$HOME/micromamba/envs/protal-db-build/bin:$PATH
pin=(taskset -c 0-3 nice -n 5)
"${pin[@]}" cmake --build "$B" -j 4 --target simulate_metagenomes > "$W/sim_build.log" 2>&1 || { echo "simulator build FAIL"; exit 1; }
rm -rf "$W/mini_db"
cd "$S" || exit 1
start=$(date +%s)
PROTAL=$B/protal "${pin[@]}" bash scripts/mini_db/build_mini_db.sh "$W/mini_db" > "$W/mini_db.log" 2>&1
echo "mini db exit $? in $(( $(date +%s) - start )) s"
grep -E "column_weights|strain_alleles|Wrote .*\.bin" "$W/mini_db.log" | head -8
ls -la "$W/mini_db/protal_db" | head -20
start=$(date +%s)
PROTAL_TEST_DB=$W/mini_db/protal_db PROTAL=$B/protal SIMULATE=$B/simulate_metagenomes "${pin[@]}" python3 -m unittest tests/e2e/test_protal_e2e.py > "$W/e2e.log" 2>&1
echo "e2e exit $? in $(( $(date +%s) - start )) s"
tail -5 "$W/e2e.log"
grep -E "^(FAIL|ERROR):" "$W/e2e.log" | head -20
