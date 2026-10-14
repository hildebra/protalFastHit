#!/bin/bash
# The baseline's end-to-end tests a second time (tmp_base2), to compare the baseline with itself.
C=~/split-cmp
mkdir -p $C/tmp_base2
cd ~/split-build/src
TMPDIR=$C/tmp_base2 PROTAL_TEST_KEEP=1 PROTAL_TEST_DB=$C/db_base PROTAL=$HOME/split-base/build/protal SIMULATE=$HOME/split-base/build/simulate_metagenomes PROTAL_TESTS_REQUIRED=1 \
    taskset -c 0-3 nice -n 5 python3 -m unittest tests/e2e/test_protal_e2e.py > $C/e2e_base2.log 2>&1
echo "e2e base again: exit $?; $(tail -1 $C/e2e_base2.log)"
bash $(cd "$(dirname "$0")" && pwd)/compare_outputs.sh $C tmp_base tmp_base2
