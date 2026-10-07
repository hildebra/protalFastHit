#!/usr/bin/env bash
# no_outdir.sh - protal without -o: --profile_only and a plain -1/-2 run, one and two samples, in empty folders
# (docs/claude/2026-10-07-sam-combine). Needs run.sh's work folder (its reads and SAMs).
set -u
PROTAL=${PROTAL:-$HOME/samcombine/bin/protal}
DB=${DB:-$HOME/samcombine/db/database.protal}
W=${W:-$HOME/samcombine/work}
RUN="taskset -c 0-3 nice -n 5"
R=$W/reads/run1
cd "$W" || exit 1
rm -rf no_o; mkdir -p no_o
try() {
  local name=$1; shift
  mkdir -p "no_o/$name"
  (cd "no_o/$name" && $RUN "$PROTAL" --db "$DB" -t 4 --no_qcmsa "$@") > "no_o/$name.log" 2>&1
  local rc=$?
  echo "$name: exit $rc; files: $(cd "no_o/$name" && find . -type f | sort | tr '\n' ' ')"
  grep -E "what\(\)|Error" "no_o/$name.log" | head -2 | sed 's/^/    /'
}
cp "$W/run1/alignments/sa.sam.zst" no_o/sa_copy.sam.zst
try profile_only_one --profile_only "$W/no_o/sa_copy.sam.zst"
try reads_one -1 "$R/sa_1.fq" -2 "$R/sa_2.fq" --prefix sa
try reads_two -1 "$R/sa_1.fq,$R/sb_1.fq" -2 "$R/sa_2.fq,$R/sb_2.fq" --prefix sa,sb
try reads_two_no_strains -1 "$R/sa_1.fq,$R/sb_1.fq" -2 "$R/sa_2.fq,$R/sb_2.fq" --prefix sa,sb --no_strains
try reads_two_dot -1 "$R/sa_1.fq,$R/sb_1.fq" -2 "$R/sa_2.fq,$R/sb_2.fq" --prefix sa,sb -o .
