#!/usr/bin/env bash
# Run the container-focused e2e tests against the mini DB (2 threads via internal calls).
set -u
cd ~/strain-build/src
export PROTAL=~/strain-build/bin/protal
export SIMULATE=~/strain-build/bin/simulate_metagenomes
export PROTAL_TEST_DB=~/strain-build/mini_db/protal_db
export OMP_NUM_THREADS=2
python3 -m unittest -v tests.e2e.test_protal_e2e.CompressedDatabaseTest tests.e2e.test_protal_e2e.ReadTypeModelTest 2>&1 | tail -70
