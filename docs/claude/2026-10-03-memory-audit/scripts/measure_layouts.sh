#!/usr/bin/env bash
# Builds index_layout (index_layout_db: measure_db_tool.sh, it needs a Linux-side copy of src/) and runs them on the local worlds (WSL paths; edit the top).
# Raw copies of the single-file databases are made once with protal --decompress_db.
SRC=${SRC:-/mnt/c/Users/hildebra/Documents/locDev/protal}
HERE=$SRC/docs/claude/2026-10-03-memory-audit/scripts
WORK=${WORK:-$HOME/protal-mem2}
PROTAL=${PROTAL:-$HOME/r226-build/build/protal}
mkdir -p "$WORK" && cd "$WORK" || exit 1
g++ -O2 -std=c++17 -o index_layout "$HERE/index_layout.cpp" || exit 1
echo "== db060_full (raw index, 0.6 world) =="
./index_layout ~/bench071/db060_full/index.prx
for d in db900n dbdense; do
  if [ ! -f "raw_$d/index.prx" ]; then
    mkdir -p "raw_$d" && cp ~/protal-perf/$d/database.protal "raw_$d/" && "$PROTAL" --decompress_db --db "raw_$d" -t 6 > "raw_$d.log" 2>&1
  fi
  echo "== $d (raw) =="
  ./index_layout "raw_$d/index.prx"
done
