#!/bin/bash
# Run every suite as CI does, timed, on CPUs 0-3, niced. Usage: run_all.sh SRC BUILD RESULTS [suites...]
# suites: unit mini e2e example minidb scripts (default: all)
SP=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/65c3e396-1c71-480e-8173-2f2320bc22de/scratchpad
SRC=${1:-~/testspeed/src}
B=${2:-~/testspeed/build}
R=${3:-~/testspeed/results}
shift 3
suites=${*:-unit mini e2e example minidb scripts}
mkdir -p $R
export PATH=~/micromamba/envs/protal-db-build/bin:$PATH
export PROTAL_TESTS_REQUIRED=1
export OMP_NUM_THREADS=4
TS="taskset -c 0-3 nice -n 5"
cd $SRC
python3 -c "import sys, sklearn, numpy, pandas; print(sys.version.split()[0], 'sklearn', sklearn.__version__, 'numpy', numpy.__version__, 'pandas', pandas.__version__)" > $R/python.txt 2>&1

stamp() { echo "$1 $(( $(date +%s%N) / 1000000 ))" >> $R/times.txt; }

uptime >> $R/times.txt
for s in $suites; do
  case $s in
  unit)
    stamp unit_start
    $TS ctest --test-dir $B -j4 --output-on-failure > $R/ctest.out 2>&1; echo "ctest exit $?" >> $R/times.txt
    stamp unit_end
    $TS $B/tests/protal_tests --gtest_brief=0 > $R/gtest_serial.out 2>&1; echo "gtest serial exit $?" >> $R/times.txt
    stamp unit_serial_end
    ;;
  mini)
    stamp mini_start
    rm -rf $SRC/data/mini_db
    PROTAL=$B/protal $TS bash scripts/mini_db/build_mini_db.sh data/mini_db > $R/mini.out 2>&1; echo "mini exit $?" >> $R/times.txt
    stamp mini_end
    ;;
  ciexample)
    stamp ciexample_start
    rm -rf $SRC/data/mini_db
    PROTAL=$B/protal THREADS=4 $TS bash examples/mini_db/run.sh $SRC/data/mini_db > $R/ciexample.out 2>&1; echo "ciexample exit $?" >> $R/times.txt
    stamp ciexample_end
    ;;
  e2e)
    stamp e2e_start
    PROTAL_TEST_DB=$SRC/data/mini_db/protal_db PROTAL=$B/protal SIMULATE=$B/simulate_metagenomes \
      $TS python3 $SP/timed_unittest.py $R/e2e.tsv tests/e2e/test_protal_e2e.py > $R/e2e.out 2>&1; echo "e2e exit $?" >> $R/times.txt
    stamp e2e_end
    ;;
  example)
    stamp example_start
    rm -rf $R/example_work
    PROTAL=$B/protal THREADS=4 $TS bash examples/mini_db/run.sh $R/example_work > $R/example.out 2>&1; echo "example exit $?" >> $R/times.txt
    stamp example_end
    ;;
  minidb)
    stamp minidb_start
    PROTAL=$B/protal SIMULATE=$B/simulate_metagenomes \
      $TS python3 $SP/timed_unittest.py $R/minidb.tsv scripts/mini_db/test_*.py > $R/minidb.out 2>&1; echo "minidb exit $?" >> $R/times.txt
    stamp minidb_end
    ;;
  scripts)
    stamp scripts_start
    PROTAL=$B/protal SIMULATE=$B/simulate_metagenomes \
      $TS python3 $SP/timed_unittest.py $R/scripts.tsv scripts/test_insilico_strains.py scripts/test_trace_relatives.py \
      scripts/test_error_reads.py scripts/test_composition_accuracy.py scripts/test_profile_scripts.py \
      scripts/test_strain_scripts.py scripts/test_model_pmml.py scripts/test_ancestry_sites.py scripts/test_foreign_rates.py \
      > $R/scripts.out 2>&1; echo "scripts exit $?" >> $R/times.txt
    stamp scripts_end
    ;;
  esac
done
uptime >> $R/times.txt
echo done >> $R/times.txt
