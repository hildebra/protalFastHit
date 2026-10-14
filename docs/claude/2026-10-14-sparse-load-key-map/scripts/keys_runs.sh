#!/bin/bash
# The key map step (core -> slots) in four layouts, one thread on core 4: real reads on the reduced r226 database
# (paired-end read ahead with a 2 us gap, HiFi streamed), the strain world's grid database whole and at level 4, and
# synthetic key maps at r226's density (whole) and a C 16 sparse one.
K=~/sparseprof/build/keybench
OUT=~/sparseprof/runs/keys
mkdir -p $OUT
R226=~/Mocks/reduced_db_r226/database.protal
GRID=~/strainfn/db_grid/database.protal
PE=~/Mocks/reads/d6331_pe_labs/d6331_hku_SRR23926922/SRR23926922_1.fq.gz
HIFI=~/Mocks/reads/d6331_hifi/d6331_hifi_SRR13128012/SRR13128012.fq.gz
STRAIN=$OUT/strain_deep_1.fq.gz
[ -f $STRAIN ] || cat ~/strainfn/reads/deep_*_1.fq.gz > $STRAIN
run() { name=$1; shift; echo "== $name: $*" >> $OUT/log.txt; taskset -c 4 $K "$@" > $OUT/$name.tsv 2>> $OUT/log.txt; cat $OUT/$name.tsv >> $OUT/log.txt; }
run r226red_pe keys $R226 $PE 300000 1 3 0 2000
run r226red_hifi keys $R226 $HIFI 3000 0 3 0
run grid_pe keys $GRID $STRAIN 100000 1 5 0 2000
run grid4_pe keys $GRID $STRAIN 100000 1 5 4 2000
run synth_r226_pe synth 0.092 44 8000000 3 27 1 2000
run synth_r226_long synth 0.092 44 8000000 3 3000 0
run synth_c16_pe synth 0.055 3 8000000 3 27 1 2000
run synth_c16_long synth 0.055 3 8000000 3 3000 0
echo done >> $OUT/log.txt
