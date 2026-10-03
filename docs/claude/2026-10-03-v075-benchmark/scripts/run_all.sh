#!/usr/bin/env bash
# The whole benchmark, in order, logging to $B/logs/v075_benchmark.log: 0.7.5 built, its pipeline, the paired-end
# and single-end runs, the long reads, the strains, then the scores (the env's python for the scorers).
HERE=$(cd "$(dirname "$0")" && pwd)
B=${BENCH:-$HOME/bench071}
PY=${PY:-$HOME/micromamba/envs/protal-db-build/bin/python}
mkdir -p $B/logs
{
  echo "== $(date) start =="
  bash $HERE/profile.sh && echo "== profile done $(date +%T)" || echo "== profile FAILED"
  bash $HERE/long_reads.sh && echo "== long reads done $(date +%T)" || echo "== long reads FAILED"
  bash $HERE/strain_runs.sh && echo "== strains done $(date +%T)" || echo "== strains FAILED"
  $PY $HERE/score.py $B > $B/logs/v075_score.log 2>&1 && echo "== scored" || { echo "== scoring FAILED"; tail -20 $B/logs/v075_score.log; }
  $PY $HERE/strain_score.py $B > $B/logs/v075_strain_score.log 2>&1 && echo "== strains scored" || { echo "== strain scoring FAILED"; tail -20 $B/logs/v075_strain_score.log; }
  $PY $HERE/../../2026-10-02-v073-benchmark/scripts/speed.py $B results_v075 results_lr075 > $B/logs/v075_speed.md 2>&1
  echo "== $(date) end =="
} > $B/logs/v075_benchmark.log 2>&1
