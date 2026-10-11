#!/usr/bin/env bash
# After impl_cg_align.sh: the build with the --ncbi_taxids check before the start-up's threads, the unit tests, the e2e
# suite, and the NCBI taxids test five times more (its exit 2 had raced the ancestry sites' loader once). 4 pinned cores.
set -uo pipefail
A=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/8e87f9e0-3682-453c-a7a2-636d4ec2ef61/scratchpad
REPO=/mnt/c/Users/hildebra/Documents/locDev/protal
until grep -q "CG ALIGN DONE" "$A/impl_cg_align.out" 2> /dev/null; do sleep 15; done
W=$HOME/impl; S=$W/work/src; B=$W/work/build
export PATH=$HOME/micromamba/envs/protal-db-build/bin:$PATH
pin=(taskset -c 0-3 nice -n 5)
echo "[$(date +%T)] rebuild and the unit tests"
bash "$A/impl_build.sh" 2>&1 | tail -3
bash "$A/impl_unit.sh" 2>&1 | tail -3
"${pin[@]}" cmake --build "$B" -j 4 --target simulate_metagenomes > "$W/sim_build.log" 2>&1 || echo "simulator build FAIL"
cd "$S" || exit 1
rsync -a --checksum --no-times "$REPO/tests/e2e/test_protal_e2e.py" "$S/tests/e2e/test_protal_e2e.py"
echo "[$(date +%T)] e2e"
start=$(date +%s)
PROTAL_TEST_DB=$W/mini_db/protal_db PROTAL=$B/protal SIMULATE=$B/simulate_metagenomes "${pin[@]}" python3 -m unittest tests/e2e/test_protal_e2e.py > "$W/e2e_final.log" 2>&1
echo "e2e exit $? in $(( $(date +%s) - start )) s"
tail -3 "$W/e2e_final.log"
grep -a -E "^(FAIL|ERROR):" "$W/e2e_final.log" | head -20
for i in 1 2 3 4 5; do
  PROTAL_TEST_DB=$W/mini_db/protal_db PROTAL=$B/protal SIMULATE=$B/simulate_metagenomes "${pin[@]}" \
      python3 -m unittest tests.e2e.test_protal_e2e.ReadTypeModelTest.test_ncbi_taxids > "$W/ncbi_$i.log" 2>&1
  echo "ncbi test $i: exit $?"
done
echo "[$(date +%T)] FINAL DONE"
