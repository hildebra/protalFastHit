#!/usr/bin/env bash
# HEAD moved to 0dea1ff (Singletrack merged, touching RunProtal.h, Build.h, Options.h, Database.h): after impl_runs3.sh,
# both trees rebuilt from the new HEAD (ref) and the new HEAD with the changes (work); the unit tests, the mini database
# and e2e tests, the mini-database and GTDB-build Python tests, and the identity runs again (db_bin remade by the new
# work build's --compress_db). 4 pinned cores.
set -uo pipefail
A=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/8e87f9e0-3682-453c-a7a2-636d4ec2ef61/scratchpad
until grep -q "RUNS3 DONE" "$A/impl_runs3.out" 2> /dev/null; do sleep 15; done
W=$HOME/impl; S=$W/work/src; B=$W/work/build
pin=(taskset -c 0-3 nice -n 5)
echo "[$(date +%T)] build ref and work on $(cat "$A/archives/head.commit")"
REF=1 TARGETS="protal protal_tests simulate_metagenomes" bash "$A/impl_build.sh" 2>&1 | tail -4
bash "$A/impl_verify_src.sh" 2>&1 | head -4
echo "[$(date +%T)] unit tests"
bash "$A/impl_unit.sh" 2>&1 | tail -3
echo "[$(date +%T)] mini database and e2e"
bash "$A/impl_e2e.sh" 2>&1 | tail -8
echo "[$(date +%T)] mini_db tests"
cd "$S" || exit 1
start=$(date +%s)
( export PATH=$HOME/micromamba/envs/protal-db-build/bin:$PATH
  PROTAL=$B/protal SIMULATE=$B/simulate_metagenomes "${pin[@]}" python3 -m unittest -v scripts/mini_db/test_*.py > "$W/minidb_tests_merged.log" 2>&1 )
echo "mini_db tests exit $? in $(( $(date +%s) - start )) s"
tail -3 "$W/minidb_tests_merged.log"
grep -a -E "^(FAIL|ERROR):" "$W/minidb_tests_merged.log" | head -20
echo "[$(date +%T)] identity"
rm -rf "$W/id/db_bin" "$W/id/runs"
bash "$A/impl_identity.sh" 2>&1 | tail -24
echo "[$(date +%T)] MERGED DONE"
