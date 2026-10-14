#!/bin/bash
# The key map step in every layout, one thread on core 4, the direct map on huge pages; each case twice: the compact
# map's huge-page advice as protal has it (KEYBENCH_OLD_ADVICE=1: refused, 4 KB pages) and fixed.
K=~/sparseprof/build/keybench
OUT=~/sparseprof/runs/keys2
mkdir -p $OUT
[ -d ~/sparseprof/runs/keys ] && [ ! -d ~/sparseprof/runs/keys_v1 ] && mv ~/sparseprof/runs/keys ~/sparseprof/runs/keys_v1
R226=~/Mocks/reduced_db_r226/database.protal
GRID=~/strainfn/db_grid/database.protal
PE=~/Mocks/reads/d6331_pe_labs/d6331_hku_SRR23926922/SRR23926922_1.fq.gz
HIFI=~/Mocks/reads/d6331_hifi/d6331_hifi_SRR13128012/SRR13128012.fq.gz
STRAIN=~/sparseprof/runs/keys_v1/strain_deep_1.fq.gz
run() {
  name=$1; shift
  for advice in old new; do
    echo "== $name $advice: $*" >> $OUT/log.txt
    if [ $advice = old ]; then e="KEYBENCH_OLD_ADVICE=1"; else e="KEYBENCH_NEW_ADVICE=1"; fi
    env $e taskset -c 4 $K "$@" > $OUT/${name}_$advice.tsv 2>> $OUT/log.txt
    cat $OUT/${name}_$advice.tsv >> $OUT/log.txt
  done
}
run synth_r226_pe synth 0.092 44 8000000 3 27 1 2000
run synth_r226_long synth 0.092 44 8000000 3 3000 0
run r226red_pe keys $R226 $PE 300000 1 3 0 2000
run r226red_hifi keys $R226 $HIFI 3000 0 3 0
run grid_pe keys $GRID $STRAIN 100000 1 5 0 2000
run synth_c16_pe synth 0.055 3 8000000 3 27 1 2000
run synth_c16_long synth 0.055 3 8000000 3 3000 0
echo done >> $OUT/log.txt
