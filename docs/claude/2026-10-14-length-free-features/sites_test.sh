#!/usr/bin/env bash
# Can per-nucleotide sites stand in for the unique k-mers? mix's and b150's new models trained again without the two
# unique seed enrichments (no unique k-mer feature left; the ancestry site features stay), on the strain world's
# l150, l75 and g4 alignments of lfv_strain.sh, scored by evaluate.py and matched_fp.pl beside new and prev.
# 4 cores, niced.
set -euo pipefail
export PATH=$HOME/micromamba/envs/protal-db-build/bin:$PATH
D=/mnt/c/Users/hildebra/Documents/locDev/protal/docs/claude/2026-10-14-length-free-features
V=~/lfv
W=~/strainfn
PROTAL=$V/bin/protal
for b in mix b150; do
  name=${b}_nouniq
  if [ ! -s $V/models/$name.xml ]; then
    perl $D/lfv_inputs.pl $V/models/${b}_new.xml $V/models/${name}_inputs.xml --drop unique_seed_enrichment,unique2_seed_enrichment
    command=$(grep -h -m1 -o "machine_learning_cmdline.py --truth-file .*" $V/$b/logs/classifier_training_c1.log)
    command=$(IN="$V/models/${name}_inputs.xml" OUT="$V/models/$name" perl -pe 's/--features \S+/--reference-pmml $ENV{IN}/; s/--output-prefix \S+/--output-prefix $ENV{OUT}/' <<< "$command")
    (cd $V/$b && taskset -c 0-3 nice -n 5 python3 $V/src/scripts/$command > $V/models/$name.log 2>&1)
  fi
  echo "== $name: $(grep -m1 -o '^[0-9]* features' $V/models/$name.report.txt)"
done
samples=$(cut -f1 $W/reads/samples.tsv)
for run in l150 l75 g4; do
  sams=$(for s in $samples; do printf '%s,' $W/run_v_$run/alignments/$s.sam.zst; done); sams=${sams%,}
  truths=$(for s in $samples; do printf '%s,' $W/cf/truth/$s.txt; done); truths=${truths%,}
  for model in mix_nouniq b150_nouniq; do
    R=$W/run_v_${run}_$model
    [ -s $R/done ] && continue
    rm -rf $R; mkdir -p $R/profiles
    ln -s $W/run_v_$run/alignments $R/alignments
    taskset -c 0-3 nice -n 5 $PROTAL --db $W/db_grid --profile_only "$sams" --profile_truth "$truths" -o $R/profiles -t 4 \
        --no_strains --no_qcmsa --model $V/models/$model.xml > $R/protal.log 2>&1
    date > $R/done
  done
done
runs=$(for run in l150 l75 g4; do for m in mix_new mix_prev mix_nouniq b150_new b150_prev b150_nouniq; do printf 'v_%s_%s ' $run $m; done; done)
taskset -c 0-3 nice -n 5 python3 -I $W/cf/evaluate.py $W $W/eval_s $runs > $W/eval_s.log 2>&1
awk 'index($0, "## calls") == 1 {p=1; next} /^## /{p=0} p' $W/eval_s.log | grep -E '^\| v_[a-z0-9_]+ \| all'
awk 'index($0, "## knobs") == 1 {p=1; next} /^## /{p=0} p' $W/eval_s.log | grep -E '^\| v_'
perl $D/matched_fp.pl $W $runs
grep -n -E "^(unique|ancestry|lu_|lsu_)" $V/models/mix_nouniq.varimp.tsv | head -12
echo "== done $(date +%T)"
