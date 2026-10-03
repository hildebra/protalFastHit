#!/usr/bin/env bash
# Mini database with the fp-anatomy binary, the e2e tests, the mini-db tests (GtdbBuildTest: the whole pipeline with
# its parity check), then a build-and-train pipeline on the benchmark world (~/bench071/world) with the new defaults.
set -uo pipefail
B=$HOME/protal-fp
cd $B/src
export PATH=$HOME/micromamba/envs/protal-db-build/bin:$PATH
PROTAL=$B/build/protal bash scripts/mini_db/build_mini_db.sh $B/mini > $B/mini.log 2>&1 || { echo "mini db failed"; tail -30 $B/mini.log; exit 1; }
echo "mini db done $(date +%T)"
grep -E "Suspect copies|Gene conservation:" $B/mini.log | head -3
PROTAL_TEST_DB=$B/mini/protal_db PROTAL=$B/build/protal SIMULATE=$B/build/simulate_metagenomes \
  timeout 5400 python3 -m unittest tests/e2e/test_protal_e2e.py > $B/e2e.log 2>&1
echo "e2e exit $? $(date +%T)"; tail -5 $B/e2e.log
PROTAL=$B/build/protal SIMULATE=$B/build/simulate_metagenomes PROTAL_TRAIN_PYTHON=$HOME/micromamba/envs/protal-db-build/bin/python \
  timeout 5400 python3 -m unittest scripts.mini_db.test_mini_db > $B/minitest.log 2>&1
echo "mini-db tests exit $? $(date +%T)"; tail -5 $B/minitest.log

# The benchmark world pipeline (as docs/claude/2026-10-03-denoising-implementation/bench.sh, pipeline "new").
W=$HOME/bench071/world
OUT=$HOME/fp_bench
mkdir -p $OUT/logs
common=(--gtdb $W/release_p --protal $B/build/protal --simulator $B/build/simulate_metagenomes -t 6 --seed 1
        --extra-genomes $W/full/simulation/genomes_nonreps --no-binary-check --read-types pe,se --samples 8
        --read-pairs 1000,20000,200000,1000000:4 --test-samples 4 --test-read-pairs 500,10000,100000,1000000:2
        --evaluation basic)
if [ ! -s $OUT/new/model_logs/summary.txt ]; then
  echo "$(date +%T) pipeline new"
  /usr/bin/time -v -o $OUT/logs/new.time python3 $B/src/scripts/build_gtdb_database.py "${common[@]}" --outdir $OUT/new \
    > $OUT/logs/new.log 2>&1 || { echo "pipeline new failed"; tail -30 $OUT/logs/new.log; exit 1; }
  grep -E "Elapsed|Maximum resident" $OUT/logs/new.time
fi
cat $OUT/new/model_logs/summary.txt
grep -E "Suspect copies" $OUT/logs/new.log $OUT/new/*.log 2>/dev/null | head -3
echo "BENCHDONE $(date +%T)"
