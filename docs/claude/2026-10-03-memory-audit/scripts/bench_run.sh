#!/usr/bin/env bash
# Builds bench_lookup_layouts (first with AddressSanitizer on the small world) and runs it on the local raw indexes.
SRC=${SRC:-/mnt/c/Users/hildebra/Documents/locDev/protal}
HERE=$SRC/docs/claude/2026-10-03-memory-audit/scripts
WORK=${WORK:-$HOME/protal-mem2}
cd "$WORK" || exit 1
cp "$HERE/bench_lookup_layouts.cpp" . || exit 1
g++ -O1 -g -std=c++17 -mpopcnt -fsanitize=address,undefined -o bench_asan bench_lookup_layouts.cpp || exit 1
echo "== ASan, dbdense, 200k lookups =="
./bench_asan raw_dbdense/index.prx 200000 2>&1 | head -40
g++ -O2 -std=c++17 -mpopcnt -o bench_lookup_layouts bench_lookup_layouts.cpp || exit 1
for db in ~/bench071/db060_full/index.prx raw_dbdense/index.prx; do
  echo "== $db, 4M lookups =="
  ./bench_lookup_layouts "$db" 4000000
done
