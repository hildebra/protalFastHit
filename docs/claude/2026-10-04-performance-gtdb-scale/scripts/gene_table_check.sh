#!/usr/bin/env bash
set -uo pipefail
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/37e6327e-2743-4622-91c7-90dfd0315ee2/scratchpad
bash $S/build_work.sh 'GeneTables.*:UniqueKmers.*:ReferenceMap.*:ReferenceFingerprint.*' 12 || exit 1
echo "== bench on the working tree (parallel adds), 30000 taxa x 168 genes"
PROTAL_GENE_TABLE_TAXA=30000 $HOME/perf-gtdb/work/build/tests/protal_tests --gtest_filter='GeneTables.BenchLoadOfLargeTables' 2>&1 | grep -E 'genes in|thread\(s\)|PASSED|FAILED'
