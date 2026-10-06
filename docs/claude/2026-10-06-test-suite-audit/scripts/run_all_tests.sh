#!/usr/bin/env bash
# Build protal at 3bcf535 (git archive, no foreign uncommitted edits) and run every test suite with timings.
set -u
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/57f6bd73-f81f-460d-bd11-1b9ef875fd6f/scratchpad
W=$HOME/testaudit
OUT=$S/results
mkdir -p "$OUT"
rm -rf "$W"; mkdir -p "$W/src"
tar -xf "$S/head.tar" -C "$W/src"
cd "$W/src"
echo "== configure $(date)"
nice -n 10 cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DPROTAL_BUILD_TESTS=ON > "$OUT/configure.log" 2>&1 || { echo configure failed; tail -30 "$OUT/configure.log"; exit 1; }
echo "== build $(date)"
nice -n 10 cmake --build build --target protal simulate_metagenomes protal_tests -- -j4 > "$OUT/build.log" 2>&1 || { echo build failed; tail -40 "$OUT/build.log"; exit 1; }
echo "== gtest $(date)"
( time nice -n 10 build/tests/protal_tests --gtest_output=json:"$OUT/gtest.json" ) > "$OUT/gtest.log" 2>&1
echo "gtest exit $?"
tail -5 "$OUT/gtest.log"
echo "== ctest $(date)"
( time nice -n 10 ctest --test-dir build -j4 ) > "$OUT/ctest.log" 2>&1
echo "ctest exit $?"; tail -6 "$OUT/ctest.log"
echo "== python versions"
python3 --version; python3 -c 'import numpy; print("numpy", numpy.__version__)'
python3 -c 'import sklearn; print("sklearn", sklearn.__version__)' 2>&1 | tail -1
ls ~/protal-train/bin/python 2>&1; command -v art_illumina || echo "no art_illumina"
echo "== mini_db unit tests $(date)"
( time PROTAL=build/protal SIMULATE=build/simulate_metagenomes nice -n 10 python3 -m unittest -v --durations 0 scripts/mini_db/test_mini_db.py ) > "$OUT/mini_db_tests.log" 2>&1
echo "mini_db exit $?"; tail -5 "$OUT/mini_db_tests.log"
echo "== insilico + trace $(date)"
( time nice -n 10 python3 -m unittest -v --durations 0 scripts/test_insilico_strains.py scripts/test_trace_relatives.py ) > "$OUT/script_tests.log" 2>&1
echo "scripts exit $?"; tail -5 "$OUT/script_tests.log"
echo "== model pmml $(date)"
PYT=python3; [ -x ~/protal-train/bin/python ] && PYT=~/protal-train/bin/python
( time nice -n 10 $PYT -m unittest -v --durations 0 scripts/test_model_pmml.py ) > "$OUT/pmml_tests.log" 2>&1
echo "pmml exit $?"; tail -5 "$OUT/pmml_tests.log"
echo "== mini db build $(date)"
( time PROTAL=build/protal nice -n 10 bash scripts/mini_db/build_mini_db.sh data/mini_db ) > "$OUT/mini_db_build.log" 2>&1
echo "mini-db build exit $?"; tail -3 "$OUT/mini_db_build.log"
echo "== e2e $(date)"
( time PROTAL_TEST_DB=data/mini_db/protal_db PROTAL=build/protal SIMULATE=build/simulate_metagenomes nice -n 10 python3 -m unittest -v --durations 0 tests/e2e/test_protal_e2e.py ) > "$OUT/e2e.log" 2>&1
echo "e2e exit $?"; tail -5 "$OUT/e2e.log"
echo "== done $(date)"
