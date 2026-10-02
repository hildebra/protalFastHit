#!/bin/bash
# build_gtdb_database.py on the operon world (~/opw: 900 species, 765 in the release, 3 genomes each, 8% archaea,
# marker clusters broken up per family; make_operon_world.sh) with gene neighbours from all its genomes, in a reduced
# design: 4 samples per design point, 3 read-pair counts, 3 long-read depths, 2 test samples per design point.
# usage: opw_build.sh OUT   (scripts of the working tree in ~/bprof/gnb/src, binaries in ~/bprof/gnb/build)
set -euo pipefail
out=$1
R=$HOME/bprof/gnb/src B=$HOME/bprof/gnb/build env=$HOME/micromamba/envs/protal-db-build
export PATH=$env/bin:$PATH
rm -rf "$out"
mkdir -p "$(dirname "$out")"
/usr/bin/time -v -o "$out.time" "$env/bin/python3" "$R/scripts/build_gtdb_database.py" --gtdb ~/opw/release_p \
  --extra-genomes ~/opw/world/simulation/genomes_nonreps --outdir "$out" --protal "$B/protal" \
  --simulator "$B/simulate_metagenomes" -t 5 --seed 1 --holdout 0.1 --samples 4 --read-pairs 1000,10000,100000 \
  --species-per-sample 10-40 --archaea 2 --long-read-bases 300000,3000000,30000000 --test-samples 2 \
  --test-read-pairs 1000,10000,100000 --test-species-per-sample 10-60 --test-long-read-bases 300000,3000000,30000000 \
  --pbsim "$env/bin/pbsim" --pbsim-models "$env/data" "${@:2}" > "$out.log" 2>&1
echo "rc $?" >> "$out.log"
