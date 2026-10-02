#!/bin/bash
# The zstd level of a finished database: size and time of protal --compress_db at levels 3-19 on one
# database's raw files (WSL, 6 threads, nice 10). Usage: compress_levels.sh DATABASE.protal|RAW_FOLDER PROTAL [THREADS]
set -e
DB=${1:-~/tune/V3/protal_db/database.protal}
PROTAL=${2:-~/progress-test/bin/protal}
T=${3:-6}
W=~/compress-levels
rm -rf "$W" && mkdir -p "$W/raw"
if [ -d "$DB" ]; then  # a folder of raw files
    cp "$DB"/* "$W/raw/"
else
    cp "$DB" "$W/raw/database.protal"
    "$PROTAL" --decompress_db --db "$W/raw" -t "$T" > "$W/decompress.log" 2>&1
fi
ls -la "$W/raw"
raw=$(du -sb "$W/raw" | cut -f1)
index=$(stat -c %s "$W/raw/index.prx")
echo "raw files: $raw bytes, index.prx $index bytes"
for level in 3 6 9 12 15 19; do
    rm -rf "$W/l$level" && cp -r "$W/raw" "$W/l$level"
    start=$(date +%s.%N)
    nice -n 10 "$PROTAL" --compress_db --db "$W/l$level" --compress_level "$level" -t "$T" > "$W/l$level.log" 2>&1
    end=$(date +%s.%N)
    size=$(stat -c %s "$W/l$level/database.protal")
    awk -v l=$level -v s=$size -v r=$raw -v a=$start -v b=$end \
        'BEGIN{printf "level %2d: database.protal %d bytes, %.3f of the raw files, %.1f s\n", l, s, s/r, b-a}'
    rm -rf "$W/l$level"
done
