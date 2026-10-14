#!/bin/bash
# The rebased split (split-units on 5369c99) against 5369c99: both built (the baseline anew in ~/split-base, the split
# incrementally in ~/split-build), the split's unit tests, then the database and end-to-end comparison.
S=$(cd "$(dirname "$0")" && pwd)
rm -rf ~/split-base
cp $S/base2.tar $S/base.tar
echo "== baseline 5369c99"
bash $S/build_base.sh protal simulate_metagenomes
echo "== split"
bash $S/build_split.sh protal simulate_metagenomes protal_tests
bash $S/unit_split.sh
echo "== comparison"
bash $S/compare_runs.sh
bash $S/compare_outputs.sh
