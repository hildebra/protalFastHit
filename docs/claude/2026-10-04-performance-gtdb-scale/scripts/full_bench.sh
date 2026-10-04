#!/usr/bin/env bash
# The gene-table bench at GTDB r226 size (143,614 genomes x 168 genes = 24.1M rows per table), after a build with tests.
set -uo pipefail
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/37e6327e-2743-4622-91c7-90dfd0315ee2/scratchpad
bash $S/build_work.sh 'GeneTables.*:UniqueKmers.*:ReferenceMap.*' 6 || exit 1
echo "== r226-sized bench (working tree)"
PROTAL_GENE_TABLE_TAXA=143614 $HOME/perf-gtdb/work/build/tests/protal_tests --gtest_filter='GeneTables.BenchLoadOfLargeTables' 2>&1 | grep -E 'genes in|thread\(s\)|Gene tables|PASSED|FAILED'
