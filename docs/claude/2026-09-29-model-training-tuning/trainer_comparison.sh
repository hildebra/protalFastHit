#!/bin/bash
# Old trainer (git HEAD, sklearn2pmml + Java) vs the new one, both worlds; every model scored by protal
# (--profile_only on the independent test samples' SAMs).
set -u
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/fbe611b7-59a8-4a3b-b159-58842b2426b7/scratchpad/audit4
R=/mnt/c/Users/hildebra/Documents/locDev/protal/scripts
PY=~/protal-train/bin/python; PROTAL=~/audit4/bin/protal; W=~/audit4/eval
mkdir -p $W
for world in toy hard; do
  case $world in
    toy)  train=~/audit3/model_train/training_data.tsv; testdir=~/audit3/model_test; db=~/audit3/db20; tax=~/audit3/db20/internal_taxonomy.dmp;;
    hard) train=~/audit4/data_train/training_data.tsv; testdir=~/audit4/data_test; db=~/audit4/db48; tax=~/audit4/db48_files/internal_taxonomy.dmp;;
  esac
  M=$W/$world; rm -rf $M; mkdir -p $M; cd $S
  /usr/bin/time -f "%e" -o $M/old_gtdb.time $PY old_random_forest_cmdline.py --truth-file $train --output-prefix $M/old_gtdb \
      --features normalized --ntree 512 --maxnodes 128 --seed 1 --threads 4 > $M/old_gtdb.log 2>&1 || echo "old_gtdb failed"
  /usr/bin/time -f "%e" -o $M/old_default.time $PY old_random_forest_cmdline.py --truth-file $train --output-prefix $M/old_default \
      --seed 1 --threads 4 > $M/old_default.log 2>&1 || echo "old_default failed"
  /usr/bin/time -f "%e" -o $M/new_fit.time $PY $R/random_forest_cmdline.py --truth-file $train --output-prefix $M/new_fitonly \
      --taxonomy $tax --threads 4 --evaluation none > $M/new_fitonly.log 2>&1 || echo "new fit failed"
  /usr/bin/time -f "%e" -o $M/new.time $PY $R/random_forest_cmdline.py --truth-file $train --output-prefix $M/new \
      --taxonomy $tax --threads 4 > $M/new.log 2>&1 || echo "new failed"
  /usr/bin/time -f "%e" -o $M/new_all.time $PY $R/random_forest_cmdline.py --truth-file $train --output-prefix $M/new_all \
      --features all --taxonomy $tax --threads 4 --evaluation basic > $M/new_all.log 2>&1 || echo "new_all failed"
  cp /mnt/c/Users/hildebra/Documents/locDev/protal/scripts/random_forest.xml $M/shipped.xml
  for model in shipped old_gtdb old_default new new_all; do
    [ -s $M/$model.xml ] || continue
    for point in $testdir/points/*; do
      meta=$point/sim/protal.meta
      sams=$(awk -F'\t' -v d=$point/protal/alignments '!/^#/{printf "%s%s/%s", s, d, $4; s=","}' $meta)
      truths=$(awk -F'\t' '!/^#/{printf "%s%s", s, $7; s=","}' $meta)
      $PROTAL --db $db --profile_only $sams --profile_truth $truths -o $M/score_$model/$(basename $point) --model $M/$model.xml \
          -t 4 --no_strains --no_qcmsa > $M/score_$model.$(basename $point).log 2>&1 || echo "scoring $world $model $(basename $point) failed"
    done
  done
  # protal's run time with each model on one small point (best of 3), the shipped one included
  point=$(ls -d $testdir/points/rl150_p1000); meta=$point/sim/protal.meta
  sams=$(awk -F'\t' -v d=$point/protal/alignments '!/^#/{printf "%s%s/%s", s, d, $4; s=","}' $meta)
  for model in shipped old_gtdb old_default new; do
    xml=$M/$model.xml; [ $model = shipped ] && xml=/mnt/c/Users/hildebra/Documents/locDev/protal/scripts/random_forest.xml
    [ -s $xml ] || continue
    for i in 1 2 3; do
      /usr/bin/time -f "%e %M" -a -o $M/runtime_$model.txt $PROTAL --db $db --profile_only $sams -o /tmp/rt_$world --model $xml -t 4 --no_strains --no_qcmsa > /dev/null 2>&1
    done
    echo "$model $(stat -c %s $xml) $(sort -n $M/runtime_$model.txt | head -1)" >> $M/runtime.txt
  done
done
$PY $S/eval_models.py $W > $W/summary.txt 2>&1; cat $W/summary.txt
