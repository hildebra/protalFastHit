#!/bin/bash
# The strain world's databases: the converter once (gtdb_to_protal_db.py on ~/strainfn/world), then one protal --build
# per variant on the same representatives, each with its own full reference (the test half never; the allele half's
# first N genomes per focal species), into ~/strainfn/db_<variant>. Build time and peak memory by /usr/bin/time.
# usage: make_variants.sh <protal binary> <variant>... (variant = name:N:protal flags, flags comma-separated)
set -eu
HERE=$(cd "$(dirname "$0")" && pwd)
PROTAL=$1; shift
PY=~/micromamba/envs/protal-db-build/bin/python
S=~/strainfn
BASE=$S/converted
if [ ! -s $BASE/reference.fna ]; then
  rm -rf $BASE
  taskset -c 0-3 nice -n 5 $PY -I $HERE/../../../scripts/mini_db/gtdb_to_protal_db.py --gtdb $S/world --outdir $BASE -t 4 \
      > $S/convert.log 2>&1
fi
FULL=$BASE/full_reference.fna.zst
for spec in "$@"; do
  IFS=: read -r name n flags <<< "$spec"
  DB=$S/db_$name
  [ -s $DB/database.protal ] && { echo "== $name exists"; continue; }
  rm -rf $DB; mkdir -p $DB
  for f in reference.fna reference.map internal_taxonomy.dmp gene2geneid.tsv genome2tiid.tsv species_priors.tsv \
           species_tree.nwk column_msa; do
    [ -e $BASE/$f ] && ln -s $BASE/$f $DB/$f
  done
  # Keep every genome but the test half and the allele half's genomes of rank above N (N=all: every genome, the test
  # half too: the oracle whose strains are their own alleles).
  awk -F'\t' -v n=$n 'NR > 1 && n != "all" && ($2 == "test" || $5 > n) { print substr($1, 4) }' $S/world/halves.tsv \
      > $DB/left_out.txt
  zstdcat $FULL | awk 'FILENAME == ARGV[1] { out[$1] = 1; next } /^>/ { keep = !($2 in out) } keep' $DB/left_out.txt - \
      > $DB/full_reference.fna
  echo "== $name: N=$n flags '${flags//,/ }': $(grep -c '^>' $DB/full_reference.fna) full-reference records, $(wc -l < $DB/left_out.txt) genomes left out"
  /usr/bin/time -v -o $DB/build.time taskset -c 0-3 nice -n 5 $PROTAL --build --no_profile -t 4 --db $DB \
      --reference $DB/reference.fna --full_reference $DB/full_reference.fna ${flags//,/ } > $DB/build.log 2>&1
  grep -E "^(Strain alleles|Index alleles)" $DB/build.log | cut -c1-400
  grep -E "Maximum resident|Elapsed" $DB/build.time
  ls -la $DB/database.protal
done
echo "== VARIANTS DONE $(date +%T)"
