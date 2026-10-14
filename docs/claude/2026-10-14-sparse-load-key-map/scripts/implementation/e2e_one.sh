#!/bin/bash
# The compact keys e2e test alone, with its duration.
T=~/keysfix/new
cd $T/src
PROTAL_TEST_DB=$T/mini_db/protal_db PROTAL=$T/build/protal SIMULATE=$T/build/simulate_metagenomes PROTAL_TESTS_REQUIRED=1 \
  taskset -c 4,5 python3 -m unittest --durations 5 -v tests.e2e.test_protal_e2e.SparseFirstGatesTest.test_compact_keys_give_the_same_records 2>&1 | tail -12
