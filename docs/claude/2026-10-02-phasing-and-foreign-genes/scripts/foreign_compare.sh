#!/usr/bin/env bash
# Foreign genes on the operon world's test samples (build ~/opw/b_gn2 of 2026-10-02, against its training database,
# which lacks the held-out species, with the models it trained), as the collector profiles them:
#   plain_keep    the database as built (gene_neighbours.tsv without species lines), --keep_foreign_genes
#   species_keep  its gene neighbours derived again with species lines (gene_neighbours.py --from_positions on its
#                 own gene_positions.tsv), --keep_foreign_genes
#   species_drop  the same, foreign genes left out (the default)
# Samples: the 6 paired-end samples of 100,000 pairs, the 6 PacBio and 6 Nanopore samples of 30 Mb.
# usage: foreign_compare.sh [OUT (default ~/fg)] [THREADS]
set -euo pipefail
out=${1:-$HOME/fg} threads=${2:-6}
build=${BUILD:-$HOME/opw/b_gn2}
P=${PROTAL:-$HOME/protal-hap/build/protal}
S=${SCRIPTS:-$HOME/protal-hap/src/scripts}
py=$HOME/micromamba/envs/protal-db-build/bin/python3
here=$(cd "$(dirname "$0")" && pwd)
pattern='^(rl100_p100000|pb_b30000000|ont_b30000000)_s_[0-9]+	'
mkdir -p "$out"
db=$out/training_species
if [ ! -f "$db/gene_neighbours.tsv.plain" ]; then
  rm -rf "$db" && mkdir -p "$db"
  "$P" --unpack_db --db "$build/training_db" --unpack_dir "$db" > "$out/unpack.log" 2>&1
  mv "$db/gene_neighbours.tsv" "$db/gene_neighbours.tsv.plain"
  "$py" "$S/mini_db/gene_neighbours.py" --db "$db" --from_positions "$db/gene_positions.tsv" \
    --output "$db/gene_neighbours.tsv" > "$out/derive.log" 2>&1
fi
for arm in plain_keep species_keep species_drop; do
  dir=$out/$arm
  [ -f "$dir/done" ] && continue
  rm -rf "$dir" && mkdir -p "$dir"
  { echo -e "#OUTPUT_DIR\t$dir"; grep '^#SAMPLEID' "$build/test/profile_all/samples.map"
    grep -v '^#' "$build/test/profile_all/samples.map" | grep -P "$pattern" | \
      awk -F'\t' -v d="$dir" 'BEGIN { OFS = "\t" } { $4 = d "/" $1 ".sam.zst"; $5 = d "/" $1; $6 = d "/" $1 ".profile"; print }'
  } > "$dir/samples.map"
  case $arm in
    plain_keep) dbarg=$build/training_db; extra=--keep_foreign_genes ;;
    species_keep) dbarg=$db; extra=--keep_foreign_genes ;;
    species_drop) dbarg=$db; extra= ;;
  esac
  /usr/bin/time -f "%e s, %M kB" -o "$dir/time" nice "$P" --db "$dbarg" --map "$dir/samples.map" -t "$threads" \
    --no_strains --no_qcmsa --model "$build/trained_model.xml" --model_se "$build/trained_model_se.xml" \
    --model_pb "$build/trained_model_pb.xml" --model_ont "$build/trained_model_ont.xml" $extra > "$dir/protal.log" 2>&1
  rm -f "$dir"/*.sam.zst
  touch "$dir/done"
done
"$py" "$here/foreign_summary.py" "$out" plain_keep species_keep species_drop
