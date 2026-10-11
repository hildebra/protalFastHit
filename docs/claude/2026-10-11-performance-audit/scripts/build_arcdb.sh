#!/usr/bin/env bash
# The archaea world's finished database, with the arguments build_gtdb_database.py gave it (alleles from every genome,
# 16 per species and gene, none indexed), at normal priority on 4 pinned cores.
set -uo pipefail
D=$HOME/perf9/arc/protal_db
start=$(date +%s)
/usr/bin/time -f "%e %U %S %M" -o "$HOME/perf9/arcdb.time" taskset -c 0-3 "$HOME/perf9/v081/build/protal" --build --no_profile -t 4 \
    --db "$D" --reference "$D/reference.fna" --compress_level 9 --full_reference "$D/full_reference.fna.zst" \
    --strain_alleles 16 --allele_genome_share 1 --index_alleles 1 > "$HOME/perf9/arcdb.build.log" 2>&1
echo "exit $? in $(( $(date +%s) - start )) s: $(cat "$HOME/perf9/arcdb.time")"
grep -E "took|Memory after" "$HOME/perf9/arcdb.build.log" | tail -30
ls -la "$D"
ls ~/realworld_s1/model_logs | grep -i xml
