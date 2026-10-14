#!/usr/bin/env bash
# Unique cores counted from the alignment (UniqueCores.h, 2026-10-14), part 2: the small real world built twice again by
# the new protal (~/lfv2: mix, 15 samples a point at 75/100/150/250 bp; b150, 60 at 150 bp; as lfv_build.sh), each
# build's ungated model (new) trained again with the previous set (prev, the eight unique k-mer rates) and with both (the
# new set plus the eight rates), and the six models profiling the strain world's run_u_* alignments (uc_strain.sh), scored
# by evaluate.py and at equal false positives. 4 cores, niced.
set -euo pipefail
S=$(dirname "$(readlink -f "$0")")
D=/mnt/c/Users/hildebra/Documents/locDev/protal/docs/claude/2026-10-14-length-free-features
export PATH=$HOME/micromamba/envs/protal-db-build/bin:$PATH
V=~/lfv2
W=~/strainfn
PROTAL=$V/bin/protal
PREVIOUS="normalized+adjacency+distance+depth+divergence+unfiltered+ref+complexity+consistency+shape+neighbourhood+ancestry+gaps+untried+alleles+polymorphic+weights+strain_identity"
RATES=lu_per_kb,lsu_per_kb,lu_gene_rate,lsu_gene_rate,lu_gene_rate2,lsu_gene_rate2,lu_gene_rate3,lsu_gene_rate3
build() {
  local name=$1 samples=$2 setups=$3 out=$V/$1
  if [ ! -s $out/protal_db/database.protal ]; then
    rm -rf $out
    taskset -c 0-3 nice -n 5 python3 $V/src/scripts/build_gtdb_database.py --inputs ~/realworld --outdir $out \
        --protal $V/bin/protal --simulator $V/bin/simulate_metagenomes -t 4 --seed 1 \
        --samples $samples --read-pairs 2000,10000,50000,200000,1000000 --read-setups $setups \
        --species-per-sample 15-60 --archaea 2 --congeners 0.5:2-5 --read-types pe \
        --test-samples 15 --test-read-pairs 10000,100000,1000000 --test-species-per-sample 15-60 \
        --scenarios none --holdout-clades genus:3 --holdout-max-share 0.1 \
        --error-reads none --evaluation basic --no-binary-check --release 226 --model forest > $out.log 2>&1
  fi
  echo "== $name built $(date +%T); $(grep -h -m1 '^Unique cores:' $out/logs/*.log | cut -c1-120)"
}
build mix 15 75:HS25:250:40,100:HS20:300:40,150:HSXt:350:50,250:MSv3:550:50
build b150 60 150:HSXt:350:50
mkdir -p $V/models
train() {
  local b=$1 name=$2; shift 2
  [ -s $V/models/${b}_$name.xml ] && return
  local command
  command=$(grep -h -m1 -o "machine_learning_cmdline.py --truth-file .*" $V/$b/logs/classifier_training_c1.log)
  if [ "$1" = --features ]; then
    command=$(SET="$2" OUT="$V/models/${b}_${name}" perl -pe 's/--features \S+/--features $ENV{SET}/; s/--output-prefix \S+/--output-prefix $ENV{OUT}/' <<< "$command")
  else
    perl $D/lfv_inputs.pl $V/models/${b}_new.xml $V/models/${b}_${name}_inputs.xml "$@"
    command=$(IN="$V/models/${b}_${name}_inputs.xml" OUT="$V/models/${b}_${name}" perl -pe 's/--features \S+/--reference-pmml $ENV{IN}/; s/--output-prefix \S+/--output-prefix $ENV{OUT}/' <<< "$command")
  fi
  (cd $V/$b && taskset -c 0-3 nice -n 5 python3 $V/src/scripts/$command > $V/models/${b}_${name}.log 2>&1)
}
for b in mix b150; do
  if [ ! -s $V/models/${b}_new.xml ]; then
    rm -rf $V/models/${b}_unpack
    $PROTAL --unpack_db --db $V/$b/protal_db --unpack_dir $V/models/${b}_unpack > $V/models/${b}_unpack.log 2>&1
    cp $V/models/${b}_unpack/model_pe_c1.xml $V/models/${b}_new.xml
  fi
  train $b prev --features "$PREVIOUS"
  train $b both --add $RATES
  echo "== models of $b: $(ls $V/models/${b}_*.xml | xargs -n1 basename | tr '\n' ' ')"
done
MODELS="mix_new mix_prev mix_both b150_new b150_prev b150_both"
samples=$(cut -f1 $W/reads/samples.tsv)
for run in l150 l75 g4; do
  sams=$(for s in $samples; do printf '%s,' $W/run_u_$run/alignments/$s.sam.zst; done); sams=${sams%,}
  truths=$(for s in $samples; do printf '%s,' $W/cf/truth/$s.txt; done); truths=${truths%,}
  for model in $MODELS; do
    R=$W/run_x_${run}_$model
    [ -s $R/done ] && continue
    rm -rf $R; mkdir -p $R/profiles
    ln -s $W/run_u_$run/alignments $R/alignments
    taskset -c 0-3 nice -n 5 $PROTAL --db $W/db_cores --profile_only "$sams" --profile_truth "$truths" -o $R/profiles -t 4 \
        --no_strains --no_qcmsa --model $V/models/$model.xml > $R/protal.log 2>&1
    date > $R/done
  done
done
runs=$(for run in l150 l75 g4; do for model in $MODELS; do printf 'x_%s_%s ' $run $model; done; done)
taskset -c 0-3 nice -n 5 python3 -I $W/cf/evaluate.py $W $W/eval_x $runs > $W/eval_x.log 2>&1
awk 'index($0, "## calls") == 1 {p=1; next} /^## /{p=0} p' $W/eval_x.log | grep -E '^\| x_[a-z0-9_]+ \| all'
awk 'index($0, "## knobs") == 1 {p=1; next} /^## /{p=0} p' $W/eval_x.log | grep -E '^\| x_'
perl $D/matched_fp.pl $W $runs
F=/mnt/c/Users/hildebra/Documents/locDev/protal/docs/claude/2026-10-14-sparse-index-calls/feature_shift.pl
for B in l75 g4; do
  echo "== feature shift l150 -> $B (mix_new dumps)"
  perl $F $W/run_x_l150_mix_new/profiles $W/run_x_${B}_mix_new/profiles 0 | grep -E '^(feature|rows|unique2?_seed_[a-z]+|lu_per_kb|lsu_per_kb|lu_gene_rate3|lsu_gene_rate3) '
done
echo "== done $(date +%T)"
