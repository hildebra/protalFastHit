#!/usr/bin/env bash
# The whole comparison (docs/claude/2026-10-13-sparse-first-reuse): base and reuse binaries on the archaea world of
# docs/claude/2026-10-12-seeding-profile (WSL ~/seed): the exactness checks, callgrind of the paired-end loop on 100,000
# pairs (the sparse-first report's sf_callgrind.sh) and 7 alternated rounds of timing on 1M pairs.
#   bash run_all.sh BASE_PROTAL REUSE_PROTAL WORK
set -uo pipefail
here=$(cd "$(dirname "$0")" && pwd)
cg=$here/../../2026-10-12-sparse-first/scripts/sf_callgrind.sh
base=$1 reuse=$2 work=$3
db=~/seed/db/database.protal
r1=~/seed/samples/pe1m/reads/pe1m_1_R1.fq.gz r2=~/seed/samples/pe1m/reads/pe1m_1_R2.fq.gz
s1=~/seed/samples/pe100k_R1.fq.gz s2=~/seed/samples/pe100k_R2.fq.gz
echo "== exact $(date +%T)"
bash "$here/sfr_exact.sh" "$base" "$reuse" "$db" "$r1" "$r2" "$work/exact"
echo "== callgrind $(date +%T)"
bash "$cg" "$base" "$db" "$s1" "$s2" "$work/cg" base_c1 --sparse_first 1
bash "$cg" "$reuse" "$db" "$s1" "$s2" "$work/cg" reuse_c1 --sparse_first 1
bash "$cg" "$base" "$db" "$s1" "$s2" "$work/cg" plain
echo "== speed $(date +%T)"
bash "$here/sfr_speed.sh" "$base" "$reuse" "$db" "$work/speed" 7 -- -1 "$r1" -2 "$r2"
echo "== done $(date +%T)"
