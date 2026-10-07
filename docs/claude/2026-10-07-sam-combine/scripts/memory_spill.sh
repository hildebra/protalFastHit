#!/usr/bin/env bash
# memory_spill.sh PROTAL - memory.sh with --strain_spill into $M/spill (tag "spill"), then the outputs against the
# binary of 12b059b's ("before") and what the spill folder holds after the runs (docs/claude/2026-10-07-sam-combine).
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
M=${M:-$HOME/samcombine/mem}
rm -rf "$M/spill"
EXTRA="--strain_spill $M/spill" bash "$HERE/memory.sh" "$1" spill 1 2 4 8 16
echo "files left in the spill folder: $(find "$M/spill" -type f | wc -l)"
bash "$HERE/compare_memory_runs.sh" before spill
