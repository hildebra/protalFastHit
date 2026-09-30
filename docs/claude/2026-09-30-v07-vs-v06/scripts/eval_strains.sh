#!/usr/bin/env bash
# Score run A's strain MSAs of both versions against the true genotypes, with the strain audit's
# tools: evaluate.py (copied here, with each version's qcmsa re-run) and summarize.py.
set -euo pipefail
B=${BENCH:-$HOME/bench07}
SRC6=${SRC6:-$HOME/protal-0.6.0a/src} SRC7=${SRC7:-$HOME/fix-build/src}
HERE=$(cd $(dirname $0) && pwd)
AUDIT=$HERE/../../2026-09-29-strain-audit/scripts/accuracy
for v in 06 07; do
  src=$([ $v = 06 ] && echo $SRC6 || echo $SRC7)/scripts/qcmsa.py
  python3 $HERE/make_colmap.py $src $B/qcmsa_colmap_v$v.py
  # 0.6.0a writes 0-based partitions
  QCMSA_COLMAP=$B/qcmsa_colmap_v$v.py PARTITION_OFFSET=$([ $v = 06 ] && echo 1 || echo 0) python3 $HERE/evaluate.py A $B/runs/strainA.v$v bench_v$v $B/runs/strainA.v$v.log
done
python3 $AUDIT/summarize.py bench_v06 bench_v07 | tee $B/strains_summary.txt
