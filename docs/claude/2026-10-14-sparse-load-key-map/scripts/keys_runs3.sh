#!/bin/bash
# After the second load A/B: the key map step, fixed huge-page advice, the baseline build (popcount a libgcc call, as
# protal's PrepareLookups path) against the x86-64-v3 build (popcnt inline, as a target_clones path), core 4.
bash /mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/cf462e00-4c32-4801-b980-f96668bbdc3a/scratchpad/load_ab2.sh
OUT=~/sparseprof/runs/keys3
mkdir -p $OUT
R226=~/Mocks/reduced_db_r226/database.protal
PE=~/Mocks/reads/d6331_pe_labs/d6331_hku_SRR23926922/SRR23926922_1.fq.gz
HIFI=~/Mocks/reads/d6331_hifi/d6331_hifi_SRR13128012/SRR13128012.fq.gz
run() {
  name=$1; shift
  for build in keybench keybench_v3; do
    echo "== $name $build: $*" >> $OUT/log.txt
    taskset -c 4 ~/sparseprof/build/$build "$@" > $OUT/${name}_$build.tsv 2>> $OUT/log.txt
    cat $OUT/${name}_$build.tsv >> $OUT/log.txt
  done
}
run synth_r226_pe synth 0.092 44 8000000 3 27 1 2000
run synth_r226_long synth 0.092 44 8000000 3 3000 0
run r226red_pe keys $R226 $PE 300000 1 3 0 2000
run r226red_hifi keys $R226 $HIFI 3000 0 3 0
echo done >> $OUT/log.txt
