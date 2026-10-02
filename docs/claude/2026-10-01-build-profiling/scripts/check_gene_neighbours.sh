#!/bin/bash
# Fix 5: gene_neighbours.py with the k-mer index against the version before it (HEAD's, by find), on the tuning
# world converted (gene_neighbours_cost.sh made ~/bprof/gn/db) and its 765 representatives at real size: the same
# gene_neighbours.tsv and positions, and the time of each, on THREADS threads.
# usage: check_gene_neighbours.sh NEW_SCRIPT OUT [THREADS]
set -euo pipefail
new=$1 out=$2 threads=${3:-4}
repo=/mnt/c/Users/hildebra/Documents/locDev/protal
rm -rf "$out" && mkdir -p "$out/old"
git -C "$repo" show HEAD:scripts/mini_db/gene_neighbours.py > "$out/old/gene_neighbours.py"
cp "$repo/scripts/mini_db/gtdb_to_protal_db.py" "$out/old/"
awk -F'\t' '$1 ~ /^GCF_/' ~/bprof/world/genomes.tsv > "$out/reps.tsv"
for v in old new; do
  script=$([ $v = old ] && echo "$out/old/gene_neighbours.py" || echo "$new")
  /usr/bin/time -o "$out/$v.time" -f "$v: %e s wall, %U s user, %S s system" python3 "$script" --db ~/bprof/gn/db \
    --genome_table "$out/reps.tsv" -t "$threads" --output "$out/$v.tsv" --positions "$out/$v.positions.tsv" > "$out/$v.log" 2>&1
  cat "$out/$v.time"
done
cmp -s "$out/old.tsv" "$out/new.tsv" && echo "gene_neighbours.tsv identical ($(wc -l < "$out/new.tsv") lines)" || echo "gene_neighbours.tsv DIFFERS"
cmp -s "$out/old.positions.tsv" "$out/new.positions.tsv" && echo "positions identical ($(wc -l < "$out/new.positions.tsv") lines)" || echo "positions DIFFER"
diff <(grep -v "^Wrote\|took" "$out/old.log") <(grep -v "^Wrote\|took" "$out/new.log") > /dev/null && echo "logs identical" || echo "logs differ"
