#!/usr/bin/env bash
# Foreign genes where there are some: foreign_world.py's samples (a database species T whose gene A another genome
# also carries, among other neighbours; hgt<i>, and ctl<i> without that genome) profiled against the operon world's
# full database (~/opw/b_gn2/protal_db, its gene neighbours derived again with species lines), with
# --keep_foreign_genes (keep) and without (drop, the default), with strain MSAs (no qcmsa), for pe, pb and ont reads.
# usage: foreign_world.sh [OUT (default ~/fgw)] [THREADS]
set -euo pipefail
out=${1:-$HOME/fgw} threads=${2:-6}
build=${BUILD:-$HOME/opw/b_gn2}
P=${PROTAL:-$HOME/protal-hap/build/protal}
S=${SCRIPTS:-$HOME/protal-hap/src/scripts}
py=$HOME/micromamba/envs/protal-db-build/bin/python3
here=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$out"
db=$out/full_species
if [ ! -f "$db/gene_neighbours.tsv.plain" ]; then
  rm -rf "$db" && mkdir -p "$db"
  "$P" --unpack_db --db "$build/protal_db" --unpack_dir "$db" > "$out/unpack.log" 2>&1
  mv "$db/gene_neighbours.tsv" "$db/gene_neighbours.tsv.plain"
  "$py" "$S/mini_db/gene_neighbours.py" --db "$db" --from_positions "$db/gene_positions.tsv" \
    --output "$db/gene_neighbours.tsv" > "$out/derive.log" 2>&1
fi
if [ ! -f "$out/world/samples.tsv" ]; then
  "$py" "$here/foreign_world.py" "$S" "$db" "$out/world" 8 > "$out/world.log" 2>&1
fi
samples=$(cut -f1 "$out/world/samples.tsv" | tail -n +2 | sort -uV | paste -sd,)
list() { echo $samples | tr , '\n' | sed "s|.*|$out/world/reads/&$1|" | paste -sd,; }
for type in pe pb ont; do
  for arm in keep drop; do
    dir=$out/runs/$type.$arm
    [ -f "$dir.done" ] && continue
    rm -rf "$dir" && mkdir -p "$dir"
    extra=()
    [ $arm = keep ] && extra=(--keep_foreign_genes)
    if [ $type = pe ]; then reads=(-1 "$(list _R1.fq.gz)" -2 "$(list _R2.fq.gz)"); else reads=(-1 "$(list _$type.fq.gz)" --read_type $type); fi
    nice "$P" --db "$db" "${reads[@]}" --prefix "$samples" -o "$dir" -t "$threads" --no_qcmsa "${extra[@]}" > "$dir.log" 2>&1
    rm -f "$dir"/*.sam.zst
    touch "$dir.done"
  done
done
"$py" "$here/foreign_world_summary.py" "$out"
