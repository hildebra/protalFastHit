#!/usr/bin/env bash
# The final build (all changes, the binary headers read without their frames included), its unit and e2e tests, then
# impl_post.sh's checks and impl_asan.sh's sanitizer run, one after another. 4 pinned cores.
set -uo pipefail
A=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/8e87f9e0-3682-453c-a7a2-636d4ec2ef61/scratchpad
REPO=/mnt/c/Users/hildebra/Documents/locDev/protal
W=$HOME/impl; S=$W/work/src; B=$W/work/build
pin=(taskset -c 0-3 nice -n 5)
echo "[$(date +%T)] rebuild and the unit tests"
bash "$A/impl_build.sh" 2>&1 | tail -3
bash "$A/impl_verify_src.sh" 2>&1 | head -5
bash "$A/impl_unit.sh" 2>&1 | tail -3
echo "[$(date +%T)] e2e"
cd "$S" || exit 1
start=$(date +%s)
( export PATH=$HOME/micromamba/envs/protal-db-build/bin:$PATH
  PROTAL_TEST_DB=$W/mini_db/protal_db PROTAL=$B/protal SIMULATE=$B/simulate_metagenomes "${pin[@]}" python3 -m unittest tests/e2e/test_protal_e2e.py > "$W/e2e_round2.log" 2>&1 )
echo "e2e exit $? in $(( $(date +%s) - start )) s"
tail -3 "$W/e2e_round2.log"
grep -a -E "^(FAIL|ERROR):" "$W/e2e_round2.log" | head -20
echo "[$(date +%T)] FINAL DONE (round 2)"
