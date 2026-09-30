#!/bin/bash
# Settings shared by the scripts. Override on the command line, e.g.
#   PERF_DIR=/scratch/me/perf BIN=~/protal/build/protal_avx2 bash scripts/startup.sh DB
# PERF_DIR    work folder on a local (Linux) file system: databases, reads, runs
# BIN         the protal binary (Release build)
# PROTAL_SRC  the protal checkout, for the benchmarks that include its headers
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PERF_DIR=${PERF_DIR:-$HOME/protal-perf}
BIN=${BIN:-$PERF_DIR/build-rel/protal_avx2}
PROTAL_SRC=${PROTAL_SRC:-$(cd "$here/../../../.." && pwd)}
OUT=${OUT:-$PERF_DIR/io}
mkdir -p "$OUT"
# The function part of a callgrind_annotate line: count, then the name without file and arguments.
cg_lines() { grep -v "'2" | awk '{ c = $1; $1 = ""; $2 = ""; s = $0; sub(/^ +/, "", s); sub(/^[^:]*:/, "", s); sub(/\(.*/, "", s); sub(/ \[.*/, "", s); printf "%16s  %s\n", c, substr(s, 1, 120) }'; }
