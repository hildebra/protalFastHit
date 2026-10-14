#!/usr/bin/env bash
# Validation of the length-free features (2026-10-14), part 2: the models of lfv_build.sh's builds mix and b150 (each
# trained with the new default set, with the previous one and as the variants below) on the strain world's samples
# aligned by the new protal:
# ungated at 150, 100 and 75 bp (the reads cut as len_runs.sh cut them) and the grid at C 4 on its sparse index at 150
# bp. One alignment per setting, profiled with each model (--profile_only, training dumps), scored by the strain study's
# evaluate.py. 4 cores, niced.
set -euo pipefail
export PATH=$HOME/micromamba/envs/protal-db-build/bin:$PATH
V=~/lfv
W=~/strainfn
PROTAL=$V/bin/protal
PREVIOUS="normalized+adjacency+distance+depth+divergence+unfiltered+ref+complexity+consistency+shape+neighbourhood+ancestry+gaps+untried+alleles+polymorphic+weights+strain_identity"
mkdir -p $V/models
# The models: each build's ungated one (model_pe_c1.xml, --model-gates; the new default set) and its table trained again
# with the previous set.
for b in mix b150; do
  if [ ! -s $V/models/${b}_new.xml ]; then
    rm -rf $V/models/${b}_unpack
    $PROTAL --unpack_db --db $V/$b/protal_db --unpack_dir $V/models/${b}_unpack > $V/models/${b}_unpack.log 2>&1
    cp $V/models/${b}_unpack/model_pe_c1.xml $V/models/${b}_new.xml  # trained ungated, as the strain study's model
  fi
  if [ ! -s $V/models/${b}_prev.xml ]; then
    command=$(grep -h -m1 -o "machine_learning_cmdline.py --truth-file .*" $V/$b/logs/classifier_training_c1.log || true)
    [ -n "$command" ] || { echo "no trainer command for $b"; exit 1; }
    # The same command with the previous default set and its own output prefix.
    command=$(PREV="$PREVIOUS" OUT="$V/models/${b}_prev" perl -pe 's/--features \S+/--features $ENV{PREV}/; s/--output-prefix \S+/--output-prefix $ENV{OUT}/' <<< "$command")
    (cd $V/$b && taskset -c 0-3 nice -n 5 python3 $V/src/scripts/$command > $V/models/${b}_prev.log 2>&1)
  fi
  echo "== models of $b: $(ls $V/models/${b}_*.xml | tr '\n' ' ')"
done
# Variants of a build's new model, trained on its table: by the inputs of a stub (lfv_inputs.pl, --reference-pmml) or,
# with --features, by a set. mix nomapq: without the MAPQ features (mean_mapq, uniqueness, low_mapq_share); mix shares:
# the plain unique seed shares (ZU / ZS, ZT / ZS) in place of the enrichments (no gene unique share, no seed expectation).
# Then, after the first results (the new set's extra false positives: taxa of 1-5 fragments without a unique seed), for
# both builds: both, the new set with the eight unique k-mer rates it replaces added back; lf, the previous set with
# length_free (its replacements, without the read_length features).
S=$(dirname "$(readlink -f "$0")")
RATES=lu_per_kb,lsu_per_kb,lu_gene_rate,lsu_gene_rate,lu_gene_rate2,lsu_gene_rate2,lu_gene_rate3,lsu_gene_rate3
variant() {
  local b=$1 name=$2; shift 2
  [ -s $V/models/${b}_$name.xml ] && return
  local command
  command=$(grep -h -m1 -o "machine_learning_cmdline.py --truth-file .*" $V/$b/logs/classifier_training_c1.log)
  if [ "$1" = --features ]; then
    command=$(SET="$2" OUT="$V/models/${b}_${name}" perl -pe 's/--features \S+/--features $ENV{SET}/; s/--output-prefix \S+/--output-prefix $ENV{OUT}/' <<< "$command")
  else
    perl $S/lfv_inputs.pl $V/models/${b}_new.xml $V/models/${b}_${name}_inputs.xml "$@"
    command=$(IN="$V/models/${b}_${name}_inputs.xml" OUT="$V/models/${b}_${name}" perl -pe 's/--features \S+/--reference-pmml $ENV{IN}/; s/--output-prefix \S+/--output-prefix $ENV{OUT}/' <<< "$command")
  fi
  (cd $V/$b && taskset -c 0-3 nice -n 5 python3 $V/src/scripts/$command > $V/models/${b}_${name}.log 2>&1)
  echo "== model ${b}_$name: $(grep -m1 -o '^[0-9]* features ([^:]*' $V/models/${b}_${name}.report.txt || true)"
}
variant mix nomapq --drop mean_mapq,uniqueness,low_mapq_share
variant mix shares --swap unique_seed_enrichment:unique_seed_share,unique2_seed_enrichment:unique2_seed_share
for b in mix b150; do
  variant $b both --add $RATES
  variant $b lf --features "$PREVIOUS+length_free"
done
MODELS="mix_new mix_prev mix_nomapq mix_shares mix_both mix_lf b150_new b150_prev b150_both b150_lf"
# The alignments: one run per setting, no profile.
samples=$(cut -f1 $W/reads/samples.tsv)
align() {
  local name=$1 reads=$2; shift 2
  local R=$W/run_v_$name
  [ -s $R/done ] && return
  rm -rf $R; mkdir -p $R/alignments
  {
    printf '#OUTPUT_DIR\t%s\n#SAMPLEID\tFIRST\tSECOND\tSAM\tPREFIX\tPROFILE\tREAD_TYPE\n' $R
    for s in $samples; do
      printf '%s\t%s\t%s\t%s\t%s\t%s\tpe\n' $s $reads/${s}_1.fq.gz $reads/${s}_2.fq.gz $R/alignments/$s.sam.zst $R/$s $R/$s.profile
    done
  } > $R/samples.map
  taskset -c 0-3 nice -n 5 $PROTAL --db $W/db_grid --map $R/samples.map -t 4 --no_strains --no_qcmsa --no_profile "$@" > $R/protal.log 2>&1
  date > $R/done
  echo "== aligned $name"
}
align l150 $W/reads --sparse_first 0
align l100 $W/run_len100/reads --sparse_first 0
align l75 $W/run_len75/reads --sparse_first 0
align g4 $W/reads --sparse_first 4 --sparse_first_pick grid
# Each alignment profiled with each model.
for run in l150 l100 l75 g4; do
  sams=$(for s in $samples; do printf '%s,' $W/run_v_$run/alignments/$s.sam.zst; done); sams=${sams%,}
  truths=$(for s in $samples; do printf '%s,' $W/cf/truth/$s.txt; done); truths=${truths%,}
  for model in $MODELS; do
    R=$W/run_v_${run}_$model
    [ -s $R/done ] && continue
    rm -rf $R; mkdir -p $R/profiles
    ln -s $W/run_v_$run/alignments $R/alignments
    taskset -c 0-3 nice -n 5 $PROTAL --db $W/db_grid --profile_only "$sams" --profile_truth "$truths" -o $R/profiles -t 4 \
        --no_strains --no_qcmsa --model $V/models/$model.xml > $R/protal.log 2>&1
    date > $R/done
  done
  echo "== profiled $run"
done
runs=$(for run in l150 l100 l75 g4; do for model in $MODELS; do printf 'v_%s_%s ' $run $model; done; done)
taskset -c 0-3 nice -n 5 python3 -I $W/cf/evaluate.py $W $W/eval_v $runs > $W/eval_v.log 2>&1
awk 'index($0, "## calls") == 1 {p=1; next} /^## /{p=0} p' $W/eval_v.log | grep -E '^\| v_[a-z0-9_]+ \| all'
awk 'index($0, "## knobs") == 1 {p=1; next} /^## /{p=0} p' $W/eval_v.log | grep -E '^\| v_'
echo "== strain world done $(date +%T)"
