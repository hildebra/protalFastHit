#!/usr/bin/env bash
# Builds index_layout_db from a Linux-side copy of src/ (compiling from /mnt/c resolves <zstd.h> to
# protal's Zstd.h: the mount is case-insensitive) and checks it against index_layout on the raw copies;
# builds and runs the lookup benchmark on the 0.6 world's raw index.
SRC=${SRC:-/mnt/c/Users/hildebra/Documents/locDev/protal}
HERE=$SRC/docs/claude/2026-10-03-memory-audit/scripts
WORK=${WORK:-$HOME/protal-mem2}
mkdir -p "$WORK/src" && cd "$WORK" || exit 1
rsync -a --delete "$SRC/src/" "$WORK/src/" || exit 1
cp "$HERE"/*.cpp "$HERE"/*.h "$WORK/" || exit 1
g++ -O2 -std=c++20 -I src -I src/Utilities -I src/Hash index_layout_db.cpp -o index_layout_db -lzstd -pthread || exit 1
for d in db900n dbdense; do
  echo "== $d (database.protal, chunk by chunk) =="
  /usr/bin/time -f "%e s, %M KB" ./index_layout_db ~/protal-perf/$d/database.protal 6
done
g++ -O2 -std=c++17 -mpopcnt -o bench_lookup_layouts bench_lookup_layouts.cpp || exit 1
echo "== lookup benchmark, db060_full =="
./bench_lookup_layouts ~/bench071/db060_full/index.prx 4000000
echo "== lookup benchmark, dbdense =="
./bench_lookup_layouts raw_dbdense/index.prx 4000000
