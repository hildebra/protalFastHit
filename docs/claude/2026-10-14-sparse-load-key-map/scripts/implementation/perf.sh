#!/bin/bash
# protal before (HEAD) and after (this session's changes), one thread on core 4, alternating, 3 rounds:
#  load: the strain world's grid database, sparse level 2 (C 4 grid gate) direct and compact, and whole direct
#        (--sequential_load, so that "Load Index took" is the index alone; a few thousand read pairs);
#  seed: the reduced r226 database, ungated, direct and compact: D6331 pe (100k pairs) and HiFi (1,000 reads);
#        ns per lookup = the Seeding row of misc/<prefix>_runtime.tsv over the k-mers looked up.
K=~/keysfix
OUT=$K/perf
mkdir -p $OUT $K/data
cd $OUT
PE1=$K/data/pe_1.fq; PE2=$K/data/pe_2.fq; HIFI=$K/data/hifi.fq
[ -f $PE1 ] || zcat ~/Mocks/reads/d6331_pe_labs/d6331_hku_SRR23926922/SRR23926922_1.fq.gz | head -400000 > $PE1
[ -f $PE2 ] || zcat ~/Mocks/reads/d6331_pe_labs/d6331_hku_SRR23926922/SRR23926922_2.fq.gz | head -400000 > $PE2
[ -f $HIFI ] || zcat ~/Mocks/reads/d6331_hifi/d6331_hifi_SRR13128012/SRR13128012.fq.gz | head -4000 > $HIFI
GRID=~/strainfn/db_grid/database.protal
R226=~/Mocks/reduced_db_r226/database.protal
S1=~/strainfn/reads/deep_01_1.fq.gz; S2=~/strainfn/reads/deep_01_2.fq.gz
run() {  # name binary args...
  name=$1; bin=$2; shift 2
  rm -rf $OUT/$name
  taskset -c 4 $K/$bin/build/protal "$@" -o $OUT/$name -t 1 --no_strains --no_qcmsa > $OUT/$name.log 2>&1
  echo "$name exit=$?" >> $OUT/runs.txt
}
for round in 1 2 3; do
  for bin in base new; do
    run load_${bin}_whole_$round $bin --db $GRID -1 $S1 -2 $S2 --prefix s --sequential_load --no_profile
    run load_${bin}_sparse_$round $bin --db $GRID -1 $S1 -2 $S2 --prefix s --sequential_load --no_profile --sparse_first 4 --sparse_first_pick grid
    run load_${bin}_sparsec_$round $bin --db $GRID -1 $S1 -2 $S2 --prefix s --sequential_load --no_profile --sparse_first 4 --sparse_first_pick grid --compact_keys
  done
  for bin in base new; do
    for keys in direct compact; do
      extra=(); [ $keys = compact ] && extra=(--compact_keys)
      run seed_${bin}_pe_${keys}_$round $bin --db $R226 -1 $PE1 -2 $PE2 --prefix pe --read_type pe --sparse_first 0 --no_profile "${extra[@]}"
      run seed_${bin}_hifi_${keys}_$round $bin --db $R226 -1 $HIFI --prefix hifi --read_type pb --sparse_first 0 --no_profile "${extra[@]}"
    done
  done
done
echo done >> $OUT/runs.txt
