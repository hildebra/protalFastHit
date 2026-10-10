#!/bin/bash
# The small real world's inputs folder for build_gtdb_database.py --inputs, as download_gtdb.py writes it: release/
# from the two subsets of docs/claude/2026-10-10-real-ancestry (extract_subset.py: five bacterial families and one
# archaeal, their genomes' marker genes extracted, their representatives' rows of GTDB-Tk's alignments, the trees, the
# taxonomy and metadata rows), then the genomes to simulate from, fetched from NCBI by download_gtdb.py's own genome
# step (fetch_genomes.py), which reads only these genomes' metadata: every species with strains (all of its strains,
# at most 20 each in the selection) and every other species from its representative. Archaeal genomes an earlier
# download left in ~/archaea_inputs/genomes are reused. usage: make_inputs.sh <out folder>
set -eu
OUT=$1
HERE=$(cd "$(dirname "$0")" && pwd)
R=$OUT/release
mkdir -p $R/auxillary_files $OUT/genomes
for set in bac120 ar53; do
  cp -r ~/real_ancestry/$set/. $R/
done
[ -e $R/auxillary_files/sp_clusters_r226.tsv ] || cp ~/GTDB/r226/auxillary_files/sp_clusters_r226.tsv $R/auxillary_files/
echo "226.0" > $R/VERSION.txt
# The archaeal genomes already fetched.
reused=0
for acc in $(zcat $R/ar53_metadata_r226.tsv.gz | tail -n +2 | cut -f1 | sed 's/^[A-Z][A-Z]_//'); do
  if [ -s ~/archaea_inputs/genomes/$acc.fna.gz ] && [ ! -e $OUT/genomes/$acc.fna.gz ]; then
    cp ~/archaea_inputs/genomes/$acc.fna.gz $OUT/genomes/
    reused=$((reused + 1))
  fi
done
echo "archaeal genomes reused: $reused"
~/micromamba/envs/protal-db-build/bin/python $HERE/fetch_genomes.py -o $OUT --species 200 --per_species 20 \
    --rep_only_species 200 --no_tech_lookup --seed 1 -t 4 --connections 6
