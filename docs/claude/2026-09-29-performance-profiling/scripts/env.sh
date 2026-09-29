#!/bin/bash
# Settings shared by the profiling scripts. Override on the command line, e.g.
#   PERF_DIR=/scratch/me/perf PROTAL_SRC=~/protal bash scripts/build.sh
# PERF_DIR    work folder on a local (Linux) file system: builds, databases, reads, runs
# PROTAL_SRC  the protal checkout to build (default: the one this report is in)
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PERF_DIR=${PERF_DIR:-$HOME/protal-perf}
PROTAL_SRC=${PROTAL_SRC:-$(cd "$here/../../../.." && pwd)}
BIN=${BIN:-$PERF_DIR/build-rel/protal_avx2}
PROF_BIN=${PROF_BIN:-$PERF_DIR/build-prof/protal_avx2}
SIM=${SIM:-$PERF_DIR/build-rel/simulate_metagenomes}
