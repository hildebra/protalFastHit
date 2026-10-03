#!/usr/bin/env bash
# Scores every set (score.py, strain_score.py, the v0.7.3 benchmark's speed.py) and copies the result tables, the
# pipeline's evaluation and the benchmark log into this report's results/ folder.
HERE=$(cd "$(dirname "$0")" && pwd)
B=${BENCH:-$HOME/bench071}
PY=${PY:-$HOME/micromamba/envs/protal-db-build/bin/python}
D=$HERE/../results
cd $B || exit 1
[ -d runs_lr074 ] && [ ! -d runs_lr075 ] && mv runs_lr074 runs_lr075   # an earlier name of the folder
rm -rf '$D'  # a mistyped copy once made this folder
$PY $HERE/score.py $B > logs/v075_score.log 2>&1 || { echo "scoring failed"; tail -8 logs/v075_score.log; exit 1; }
$PY $HERE/strain_score.py $B > logs/v075_strain_score.log 2>&1 || { echo "strain scoring failed"; tail -5 logs/v075_strain_score.log; exit 1; }
$PY $HERE/../../2026-10-02-v073-benchmark/scripts/speed.py $B results_v075 results_lr075 > logs/v075_speed.md 2>&1 || { echo "speed failed"; tail -3 logs/v075_speed.md; }
mkdir -p $D/short $D/long_reads $D/strains $D/pipeline_0.7.5
cp results_v075/*.md $D/short/
cp results_lr075/*.md $D/long_reads/
cp results_strains_v075/summary.md $D/strains/
cp logs/v075_speed.md $D/speed.md
cp V075/model_logs/summary.txt V075/model_logs/holdout.txt $D/pipeline_0.7.5/
grep -E 'Elapsed|Maximum' logs/pipeline_0.7.5.time > $D/pipeline_0.7.5/time.txt
cp logs/v075_benchmark.log $D/benchmark.log
ls $D $D/long_reads
cat results_lr075/overall.md results_lr075/paired.md
cat logs/v075_speed.md | tail -8
