#!/bin/bash
# The reduced database's inputs folder for build_gtdb_database.py --inputs, as download_gtdb.py writes it: release/
# from the local GTDB r226 files (~/GTDB/r226), cut to the genomes of select_mock_taxa.py (every species of the mock
# benchmark's truth tables, up to 3 congeners per genus, up to 5 strains per species) by extract_subset.py (their marker
# genes, GTDB-Tk's alignment rows of the representatives, the trees, the taxonomy and metadata rows), then the genomes
# to simulate from, fetched from NCBI by download_gtdb.py's own genome step (fetch_genomes.py).
# usage: make_inputs.sh <protal clone> <out folder>     (WSL)
set -eu
SRC=$1; OUT=$2
HERE=$(cd "$(dirname "$0")" && pwd)  # this folder (select_mock_taxa.py); the clone has the repository's scripts
PY=$HOME/micromamba/envs/protal-db-build/bin/python
R=$OUT/release
mkdir -p $R/auxillary_files $OUT/genomes $OUT/taxa
taskset -c 0-3 $PY $HERE/select_mock_taxa.py ~/GTDB/r226 $SRC/scripts/mock_benchmark/truth $OUT/taxa 3 5 1
for set in bac120 ar53; do
  if [ ! -e $OUT/subset_$set.done ]; then
    taskset -c 0-3 nice -n 5 $PY $SRC/docs/claude/2026-10-10-real-ancestry/extract_subset.py ~/GTDB/r226 \
        $OUT/taxa/${set}_taxa.tsv $OUT/subset_$set $set
    touch $OUT/subset_$set.done
  fi
  cp -r $OUT/subset_$set/. $R/
done
cp ~/GTDB/r226/auxillary_files/sp_clusters_r226.tsv $R/auxillary_files/
echo "226.0" > $R/VERSION.txt
export PATH=$HOME/micromamba/envs/protal-db-build/bin:$PATH  # NCBI's datasets CLI for what the FTP server lacks
taskset -c 0-3 nice -n 5 $PY $SRC/docs/claude/2026-10-10-real-small-world/fetch_genomes.py -o $OUT --species 150 \
    --per_species 2 --rep_only_species 800 --no_tech_lookup --seed 1 -t 2 --connections 6
