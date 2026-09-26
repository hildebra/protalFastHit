#!/usr/bin/env bash
# build_mini_db.sh - simulate a sparse GTDB release and build a protal DB from it.
#
#   scripts/mini_db/build_mini_db.sh [OUTDIR] [extra simulate_gtdb_release.py args]
#
# OUTDIR (default: data/mini_db, git-ignored) receives
#   gtdb_r226/   the synthetic GTDB release (see simulate_gtdb_release.py)
#   protal_db/   the protal database; use it with --db or $PROTAL_DB_PATH
#   build.log    protal --build output
# The protal binary is $PROTAL (default: build/protal); $PROTAL_BUILD_ARGS adds options
# to protal --build (e.g. --no_compress). By default the build writes index.prx.zst
# (~10 MB here) and replaces reference.fna by reference.fna.zst; uncompressed, index.prx
# is ~3 GB whatever the reference size (the k-mer key map is a fixed 2^30 slots).
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
out=${1:-data/mini_db}
shift || true
protal=${PROTAL:-build/protal}
python=${PYTHON:-python3}

[ -x "$protal" ] || { echo "protal binary not found at '$protal' (set PROTAL=...)" >&2; exit 1; }

gtdb="$out/gtdb_r226"
db="$out/protal_db"
rm -rf "$gtdb" "$db"

"$python" "$here/simulate_gtdb_release.py" --outdir "$gtdb" --release 226 "$@"
"$python" "$here/gtdb_to_protal_db.py" --gtdb "$gtdb" --outdir "$db"

# --no_profile: without it, build mode falls through to profiling an empty sample list.
# shellcheck disable=SC2086
"$protal" --build --no_profile -t 1 ${PROTAL_BUILD_ARGS:-} \
    --db "$db" \
    --reference "$db/reference.fna" \
    --full_reference "$db/full_reference.fna" \
    > "$out/build.log" 2>&1 || { tail -30 "$out/build.log" >&2; exit 1; }

for f in index.prx reference.fna unique_kmers.tsv; do
    [ -s "$db/$f" ] || [ -s "$db/$f.zst" ] || {
        echo "protal --build did not write $db/$f (see $out/build.log)" >&2; exit 1; }
done
echo "Mini protal DB ready: $db"
ls -la "$db"
