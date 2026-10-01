#!/bin/bash
# Regenerates the benchmark data of round 2 after the work folder was removed: a build of HEAD
# ($PERF_DIR/pair/base), the 900-species world from its lineages (simulate_gtdb_release.py, as in
# scripts/mini_db/gtdb_like_lineages.py's example) and its database db900n, the read sets w900 and mix,
# and the dense world (dense_world.sh) with dbdense, dense_w and dense_mix; plain copies of all four.
#   regen_data.sh LINEAGES   (default ~/tune/lineages.txt)
set -e
source "$(dirname "$0")/env.sh"
R1=$here/../../2026-09-29-performance-profiling/scripts
lineages=${1:-$HOME/tune/lineages.txt}
mkdir -p $PERF_DIR; cd $PERF_DIR
log=$PERF_DIR/regen.log; : > $log
step() { echo "$(date +%T) $*" >> $log; }
step build base
rm -rf pair/base; mkdir -p pair/base/src
(cd $PROTAL_SRC && git archive HEAD) | tar -x -C pair/base/src
(cd pair/base && cmake -S src -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DPROTAL_BUILD_TESTS=ON > cmake.log 2>&1 &&
 cmake --build build --target protal protal_avx2 simulate_metagenomes protal_tests -j "$(nproc)" > build.log 2>&1)
mkdir -p build-rel; for b in protal protal_avx2 simulate_metagenomes; do ln -sf $PERF_DIR/pair/base/build/$b build-rel/$b; done
step world900
rm -rf world900
python3 $PROTAL_SRC/scripts/mini_db/simulate_gtdb_release.py --outdir world900 --lineages $lineages \
  --strain_divergence 0.002-0.02 --species_divergence 0.015-0.06 > world900.log 2>&1
step db900n
PROTAL_BIN=$PERF_DIR/pair/base/build/protal bash $here/prep_db.sh db900n $PERF_DIR/world900 >> $log 2>&1
step reads w900 mix
export SIM=$PERF_DIR/pair/base/build/simulate_metagenomes
bash $R1/prep_reads.sh w900 $PERF_DIR/world900/simulation/genomes.tsv 60 >> $log 2>&1
bash $R1/prep_reads.sh --mix mix w900 >> $log 2>&1
step dense
bash $here/dense_world.sh 4 60 >> $log 2>&1
PROTAL_BIN=$PERF_DIR/pair/base/build/protal bash $here/prep_db.sh dbdense $PERF_DIR/dense/world >> $log 2>&1
bash $R1/prep_reads.sh dense_w $PERF_DIR/dense/world/simulation/genomes.tsv 20 >> $log 2>&1
bash $R1/prep_reads.sh --mix dense_mix dense_w >> $log 2>&1
step plain copies
for ds in mix w900 dense_w dense_mix; do
  g=reads/$ds; [ -d $g/reads ] && g=$g/reads
  mkdir -p reads/${ds}_plain
  for m in R1 R2; do zcat $g/*_$m.fq.gz > reads/${ds}_plain/${ds}_$m.fq; done
  mkdir -p reads/${ds}_300k
  for m in R1 R2; do head -n 1200000 reads/${ds}_plain/${ds}_$m.fq > reads/${ds}_300k/${ds}_300k_$m.fq; done
done
step DONE
