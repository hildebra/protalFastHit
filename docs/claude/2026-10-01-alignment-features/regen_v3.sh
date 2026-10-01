#!/bin/bash
# Regenerates the GTDB-like tuning world of the V2 build with its seeds (../2026-09-29-model-training-tuning/setup.sh)
# and runs the V2 build command (../2026-09-30-read-type-models.md) with another protal build and the working tree's
# scripts, into OUTDIR (log: OUTDIR.log).
#     bash regen_v3.sh BUILD_DIR OUTDIR
set -eu
export PATH=$HOME/micromamba/envs/protal-db-build/bin:$PATH
B=$1
OUT=$2
T=$HOME/tune
R=/mnt/c/Users/hildebra/Documents/locDev/protal
TD=$R/docs/claude/2026-09-29-model-training-tuning
mkdir -p $T
if [ ! -d $T/release_p ]; then
  python3 $R/scripts/mini_db/gtdb_like_lineages.py --species 900 --archaea 0.08 --seed 1 > $T/lineages.txt
  /usr/bin/time -f "release %e s" python3 $R/scripts/mini_db/simulate_gtdb_release.py --outdir $T/world \
      --lineages $T/lineages.txt --seed 21 --genomes_per_species 3 --strain_divergence 0.002-0.02 \
      --species_divergence 0.015-0.06 2> $T/world.log
  tail -2 $T/world.log
  python3 $TD/make_release_p.py $T/world $T/release_p $T/unknown.txt 0.15 7
  du -sh $T/world $T/release_p
fi
/usr/bin/time -v python3 $R/scripts/build_gtdb_database.py --gtdb $T/release_p --outdir $OUT \
    --protal $B/protal --simulator $B/simulate_metagenomes -t 6 --seed 1 \
    --extra-genomes $T/world/simulation/genomes_nonreps --samples 4 \
    --read-pairs 1000,20000,200000 --test-samples 2 --test-read-pairs 500,10000,500000 \
    --long-read-bases 300000,6000000,60000000 --test-long-read-bases 150000,3000000,90000000 \
    --pbsim $HOME/micromamba/envs/protal-db-build/bin/pbsim > $OUT.log 2>&1
echo "build done: $(grep -c . $OUT.log) log lines"
du -sh $OUT
