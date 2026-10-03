#!/usr/bin/env bash
# The evaluation of this report, after bench.sh: the branch's trainer with the feature groups of
# ablation_features.patch (na+em, na+distance, na+rank) copied to ~/denoise_bench/rf2, then every model of
# eval2_cross.sh (two training designs x two feature sets x two test sets x pe, se) and eval2_ablation.sh (the groups,
# and models trained on samples of up to 200,000 read pairs) into ~/denoise_bench/eval2, and the summaries.
set -euo pipefail
R=/mnt/c/Users/hildebra/Documents/locDev/protal-denoise/docs/claude/2026-10-03-denoising-implementation
REPO=/mnt/c/Users/hildebra/Documents/locDev/protal-denoise
PY=$HOME/micromamba/envs/protal-db-build/bin/python
RF=$HOME/denoise_bench/rf2
mkdir -p $RF
for f in random_forest_cmdline.py lineages.py model_features.py model_pmml.py; do cp $REPO/scripts/$f $RF/; done
patch -d $RF -p1 < $R/ablation_features.patch
rsync -a --checksum --no-times --exclude __pycache__ $REPO/scripts/ $HOME/protal-denoise/src/scripts/
bash $R/eval2_cross.sh
bash $R/eval2_ablation.sh
$PY $R/eval2_summary.py $HOME/denoise_bench/eval2 $HOME/denoise_bench $HOME/protal-denoise/src/scripts > $R/eval_summary.txt
$PY $R/singleton.py $HOME/denoise_bench $HOME/denoise_bench/eval2 > $R/singleton_count_only.txt
$PY $R/singleton_refined.py $HOME/denoise_bench $HOME/denoise_bench/eval2 > $R/singleton_refined.txt
$PY $R/singleton_features.py $HOME/denoise_bench > $R/singleton_features.txt
echo all-done
