#!/bin/bash
# Every test suite of the working tree, on the copy built by build_pair.sh new ($PERF_DIR/pair/new):
# unit tests (Release), mini-db generator tests, end-to-end tests on a fresh mini database, the full
# GTDB build test, the model export test, and the unit tests under ASan+UBSan (Debug, as CI runs them).
# Results go to $PERF_DIR/pair/tests/*.log; the last line of each is summarised in summary.txt. TREE picks another
# copy under pair/ (default new), JOBS the parallel jobs (default 4, the ASan build all cores).
source "$(dirname "$0")/env.sh"
T=$PERF_DIR/pair/tests; rm -rf $T; mkdir -p $T
N=$PERF_DIR/pair/${TREE:-new}
cd $N/src
run() { name=$1; shift; echo "$(date +%T) $name" >> $T/progress.txt; "$@" > $T/$name.log 2>&1; echo "$name rc=$? $(tail -3 $T/$name.log | tr '\n' ' ' | cut -c1-200)" >> $T/summary.txt; }
run unit_release ctest --test-dir $N/build -j${JOBS:-4}
run mini_db_test python3 -m unittest scripts/mini_db/test_mini_db.py
rm -rf data/mini_db
run mini_db_build env PROTAL=$N/build/protal bash scripts/mini_db/build_mini_db.sh data/mini_db
run e2e env PROTAL_TEST_DB=$N/src/data/mini_db/protal_db PROTAL=$N/build/protal SIMULATE=$N/build/simulate_metagenomes \
    python3 -m unittest -v tests/e2e/test_protal_e2e.py
# A Python with scikit-learn, joblib, numpy and pandas (TRAIN_PYTHON; the first of these that has them).
for py in ${TRAIN_PYTHON:-} $HOME/protal-train/bin/python $HOME/micromamba/envs/protal-db-build/bin/python; do
  [ -x "$py" ] && "$py" -c "import sklearn, joblib, numpy, pandas" 2> /dev/null && { TRAIN_PYTHON=$py; break; }
done
if [ -x "${TRAIN_PYTHON:-}" ]; then
  run gtdb_build_test env PROTAL=$N/build/protal SIMULATE=$N/build/simulate_metagenomes PROTAL_TRAIN_PYTHON=$TRAIN_PYTHON \
      python3 -m unittest -v scripts.mini_db.test_mini_db.GtdbBuildTest
  run model_test $TRAIN_PYTHON -m unittest -v scripts/test_model_pmml.py
else
  echo "gtdb_build_test, model_test SKIPPED: no Python with scikit-learn (set TRAIN_PYTHON)" >> $T/summary.txt
fi
FLAGS="-fsanitize=address,undefined -fno-omit-frame-pointer"
cmake -S . -B $N/build-asan -G Ninja -DCMAKE_BUILD_TYPE=Debug -DPROTAL_BUILD_TESTS=ON -DCMAKE_CXX_FLAGS="$FLAGS" -DCMAKE_C_FLAGS="$FLAGS" -DCMAKE_EXE_LINKER_FLAGS="$FLAGS" > $T/asan_cmake.log 2>&1
run asan_build cmake --build $N/build-asan --target protal_tests -j "${JOBS:-$(nproc)}"
run unit_asan env ASAN_OPTIONS=detect_leaks=1 UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=1 ctest --test-dir $N/build-asan -j${JOBS:-4}
echo "ALL TESTS RUN" >> $T/summary.txt
