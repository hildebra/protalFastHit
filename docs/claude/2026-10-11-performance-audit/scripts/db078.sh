#!/usr/bin/env bash
# 0.7.8's own database of the archaea world's genes: 0.8.1's finished database unpacked (reference.fna, reference.map,
# taxonomy, tables), its index and bundle removed, and 0.7.8 --build on it (as build_gtdb_database.py calls --build,
# without --full_reference: 0.7.8 has no strain alleles). 4 pinned cores.
set -uo pipefail
W=$HOME/perf9; D=$W/db078
rm -rf "$D"; mkdir -p "$D"
"$W/v081/build/protal" --unpack_db --db "$W/arc/protal_db" --unpack_dir "$D" > "$W/db078.unpack.log" 2>&1 || { echo "unpack FAIL"; tail -5 "$W/db078.unpack.log"; exit 1; }
rm -f "$D/index.prx.zst" "$D/index.prx" "$D/database.protal" "$D/gene_table.bin"
ls -la "$D"
start=$(date +%s)
/usr/bin/time -f "%e %U %S %M" -o "$W/db078.time" taskset -c 0-3 nice -n 5 "$W/v078/build/protal" --build --no_profile -t 4 \
    --db "$D" --reference "$D/reference.fna" > "$W/db078.build.log" 2>&1
echo "build exit $? in $(( $(date +%s) - start )) s: $(cat "$W/db078.time")"
tail -5 "$W/db078.build.log"
ls -la "$D"
