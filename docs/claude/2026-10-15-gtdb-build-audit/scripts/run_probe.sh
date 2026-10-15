#!/bin/bash
# usage: run_probe.sh NAME [extra build options]
# A mini build (probe_build.py) into ~/bgaudit/probe/NAME, its console in ~/bgaudit/probe_NAME.out. Needs protal and
# simulate_metagenomes in ~/bgaudit/build, the scripts in ~/bgaudit/src (or PROBE_SRC=<src>/scripts/mini_db), and the
# protal-db-build environment (envs/protal-db-build.yaml) in ~/micromamba/envs.
HERE=$(cd "$(dirname "$0")" && pwd)
export PATH="$HOME/micromamba/envs/protal-db-build/bin:$PATH"
export PROTAL=$HOME/bgaudit/build/protal SIMULATE=$HOME/bgaudit/build/simulate_metagenomes
export PROTAL_TRAIN_PYTHON=$HOME/micromamba/envs/protal-db-build/bin/python3
cd ~/bgaudit || exit 1
name=$1; shift
taskset -c 0-3 nice -n 5 python3 "$HERE/probe_build.py" ~/bgaudit/probe "$name" "$@" > ~/bgaudit/probe_"$name".out 2>&1
rc=$?
echo "probe $name exit $rc"
tail -c 4000 ~/bgaudit/probe_"$name".out
