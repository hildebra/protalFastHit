#!/bin/bash
# gene_neighbours.py's CPU per representative genome: the tuning world converted, its first N representatives'
# real-sized genomes (pad_genomes.sh), one thread.
# usage: gene_neighbours_cost.sh OUT [N]
set -euo pipefail
out=$1 n=${2:-100}
scripts=${SCRIPTS:-$HOME/build-prof/src/scripts}
mkdir -p "$out"
[ -f "$out/db/reference.map" ] || python3 "$scripts/mini_db/gtdb_to_protal_db.py" --gtdb ~/tune/release_p --outdir "$out/db" \
  --release 226 -t 4 > "$out/convert.log" 2>&1
awk -F'\t' -v n="$n" '$1 ~ /^GCF_/ && k++ < n' ~/bprof/world/genomes.tsv > "$out/reps.tsv"
/usr/bin/time -f "wall %e s, user %U s, sys %S s" -o "$out/time.txt" python3 "$scripts/mini_db/gene_neighbours.py" \
  --db "$out/db" --genome_table "$out/reps.tsv" -t 1 --output "$out/gene_neighbours.tsv" > "$out/log" 2>&1
echo "$n representatives: $(cat "$out/time.txt")"
tail -3 "$out/log"
