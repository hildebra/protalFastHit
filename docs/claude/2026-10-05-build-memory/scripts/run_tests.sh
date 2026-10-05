#!/usr/bin/env bash
# The test suites on the change's tree (~/buildmem/new, build_change.sh): all unit tests, the mini database
# built with this protal, the end-to-end tests against it, and the GTDB build pipeline test (GtdbBuildTest).
cd $HOME/buildmem/new || exit 1
nice cmake --build build --target protal protal_tests simulate_metagenomes -j 6 > tests-build.log 2>&1 || { tail -30 tests-build.log; exit 1; }
echo "== unit tests"; ctest --test-dir build --output-on-failure -j 4 2>&1 | tail -4
echo "== mini database"; rm -rf data/mini_db
PROTAL=build/protal bash scripts/mini_db/build_mini_db.sh data/mini_db > mini_db.log 2>&1 && echo built || { tail -20 mini_db.log; exit 1; }
echo "== e2e"; PROTAL_TEST_DB=data/mini_db/protal_db PROTAL=build/protal SIMULATE=build/simulate_metagenomes \
    python3 -m unittest tests/e2e/test_protal_e2e.py > e2e.log 2>&1; tail -4 e2e.log
# The trainer needs scikit-learn: here the micromamba environment protal-db-build.
echo "== GtdbBuildTest"; PROTAL=build/protal SIMULATE=build/simulate_metagenomes \
    PROTAL_TRAIN_PYTHON=${PROTAL_TRAIN_PYTHON:-$HOME/micromamba/envs/protal-db-build/bin/python} \
    python3 -m unittest scripts.mini_db.test_mini_db.GtdbBuildTest > gtdb_build_test.log 2>&1; tail -4 gtdb_build_test.log
