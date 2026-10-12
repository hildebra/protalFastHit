#!/usr/bin/env bash
# --sparse_first's kept path changes nothing (docs/claude/2026-10-12-sparse-first): with every core sampled (C = 1) and
# any score counting (score 0, 1 hit), a pair is given up only if none of its k-mers' cores is in the index, and such a
# pair has no seeds anyway; so the SAM records (as a sorted set) and the profile must be those of a run without it.
#   bash sf_exact.sh PROTAL DATABASE R1 R2 WORK
set -uo pipefail
protal=$1 db=$2 r1=$3 r2=$4 work=$5
mkdir -p "$work"
pin=(taskset -c 0-3 nice -n 5)
run() {  # name, extra arguments
  local name=$1; shift
  rm -rf "$work/$name"
  "${pin[@]}" "$protal" --db "$db" -1 "$r1" -2 "$r2" --prefix s -o "$work/$name" -t 4 --no_strains "$@" > "$work/$name.log" 2>&1
  echo "$name exit $?"
}
records() { zstdcat "$1" | grep -v '^@' | sort | md5sum | cut -c1-16; }
run plain
run all_kept --sparse_first 1 --sparse_first_score 0 --sparse_first_hits 1
a=$(records "$work/plain/s.sam.zst"); b=$(records "$work/all_kept/s.sam.zst")
echo "SAM records: $a vs $b: $([ "$a" = "$b" ] && echo same || echo DIFFER)"
cmp -s "$work/plain/s.profile" "$work/all_kept/s.profile" && echo "profile: same" || echo "profile: DIFFERS"
grep -h "sparse first" "$work/all_kept.log"
