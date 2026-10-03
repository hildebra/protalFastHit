#!/bin/bash
# From the repository root, on the runs extracted into local/ (git-ignored): local/v5 and local/v6 from
# local/protal_r226_v{5,6}_{logs,training}.tgz (tar xzf each into its folder), and local/r226_v2_new from the earlier
# report. Retrains the four models on v5's tables with two more feature sets (about 20 minutes on 6 threads), moves the
# paired-end outputs to trained_model.* as the other runs have them, and writes compare_models_output.txt.
# Usage: PY=python3 bash docs/claude/2026-10-03-r226-v5-v6-training/run_all.sh
set -e
PY=${PY:-python3}
HERE=$(cd "$(dirname "$0")" && pwd)
L=local
for run in "relatives normalized+adjacency+relatives" "adjacency normalized+adjacency"; do
    set -- $run
    mkdir -p $L/v5_$1
    for spec in "training_data:512:" "training_data_se:512:_se" "training_data_pb:128:_pb" "training_data_ont:128:_ont"; do
        IFS=: read f leaves suffix <<< "$spec"
        $PY scripts/random_forest_cmdline.py --truth-file $L/v5/training/$f.tsv --output-prefix $L/v5_$1/trained_model$suffix \
            --features $2 --ntree 64 --maxnodes $leaves --seed 1 --threads ${THREADS:-6} \
            --taxonomy $L/v5/internal_taxonomy.dmp --evaluation full --depth-knobs --fdr-calls \
            --test-file $L/v5/test/$f.tsv > $L/v5_$1/classifier_training$suffix.log 2>&1
    done
done
$PY "$HERE/compare_models.py" v2-new=$L/r226_v2_new v6=$L/v6 v5-adj-local=$L/v5_adjacency v5-distance=$L/v5 \
    v5-relatives=$L/v5_relatives > "$HERE/compare_models_output.txt"
for suffix in "" _se _pb _ont; do
    for r in v6 v5 v5_relatives; do $PY "$HERE/error_budget.py" $L/$r $L/v5/test $suffix; done
done > "$HERE/error_budget_output.txt"
