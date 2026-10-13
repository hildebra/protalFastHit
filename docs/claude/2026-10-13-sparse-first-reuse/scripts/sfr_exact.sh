#!/usr/bin/env bash
# KeepSparseLookups changes no output: base (5334f2c) and reuse (+ the patch) on the archaea world's 1M pairs, paired-end
# and single-end (R1 alone), at --sparse_first 1 (the path that changes) and 2 (unchanged); the SAM records as a sorted
# set, the profile and the run's seeding and sparse-first lines must be the same.
#   bash sfr_exact.sh BASE_PROTAL REUSE_PROTAL DATABASE R1 R2 WORK
set -uo pipefail
base=$1 reuse=$2 db=$3 r1=$4 r2=$5 work=$6
mkdir -p "$work"
pin=(taskset -c 0-3 nice -n 5)
run() {  # binary, name, reads and options...
  local bin=$1 name=$2; shift 2
  rm -rf "$work/$name"
  "${pin[@]}" "$bin" --db "$db" --prefix s -o "$work/$name" -t 4 --no_strains "$@" > "$work/$name.log" 2>&1
  echo "$name exit $?"
}
records() { zstdcat "$1" | grep -v '^@' | sort | md5sum | cut -c1-16; }
lines() { grep -E " seeding: |sparse first|records written" "$1"; }
for c in 1 2; do
  for mode in pe se; do
    if [ "$mode" = pe ]; then reads=(-1 "$r1" -2 "$r2"); else reads=(-1 "$r1"); fi
    run "$base" "base_${mode}_c$c" "${reads[@]}" --sparse_first "$c"
    run "$reuse" "reuse_${mode}_c$c" "${reads[@]}" --sparse_first "$c"
    a=$(records "$work/base_${mode}_c$c/s.sam.zst"); b=$(records "$work/reuse_${mode}_c$c/s.sam.zst")
    echo "$mode C$c SAM records: $a vs $b: $([ "$a" = "$b" ] && echo same || echo DIFFER)"
    cmp -s "$work/base_${mode}_c$c/s.profile" "$work/reuse_${mode}_c$c/s.profile" && echo "$mode C$c profile: same" || echo "$mode C$c profile: DIFFERS"
    if diff <(lines "$work/base_${mode}_c$c.log") <(lines "$work/reuse_${mode}_c$c.log") > /dev/null; then
      echo "$mode C$c seeding and sparse-first lines: same"
    else
      echo "$mode C$c seeding and sparse-first lines: DIFFER"; diff <(lines "$work/base_${mode}_c$c.log") <(lines "$work/reuse_${mode}_c$c.log")
    fi
    lines "$work/reuse_${mode}_c$c.log"
  done
done
