#!/usr/bin/env bash
# compare_memory_runs.sh BEFORE AFTER... - the outputs of memory.sh's runs of tag BEFORE against each AFTER tag, for
# cohorts of 4 and 16, with and without strain MSAs (docs/claude/2026-10-07-sam-combine).
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
M=${M:-$HOME/samcombine/mem}
BEFORE=$1; shift
for after in "$@"; do
  for n in 4 16; do
    for s in yes no; do
      bash "$HERE/compare_runs.sh" "$M/run_${BEFORE}_${n}_$s" "$M/run_${after}_${n}_$s" | tail -4
    done
  done
done
