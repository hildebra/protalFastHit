#!/usr/bin/env bash
# Trial B2: callgrind of a profile-only run of the local pe SAM on one thread (the SAM parse's instructions).
set -uo pipefail
W=~/perf085
P=${PROTAL:-$W/bin/protal_base}
TAG=${TAG:-base}
DB=~/Mocks/reduced_db_r226/database.protal
SAM=$W/runs/pe_hku/SRR23926922.sam.zst
cd $W/runs
rm -rf cg_$TAG
taskset -c 0-3 nice -n 5 valgrind --tool=callgrind --callgrind-out-file=$W/runs/cg_$TAG.out $P --db $DB --profile_only $SAM --prefix hku -t 1 --no_strains --force -o $W/runs/cg_$TAG > $W/runs/cg_$TAG.log 2>&1
echo "callgrind exit $?"
callgrind_annotate --inclusive=yes $W/runs/cg_$TAG.out 2>/dev/null | head -70 | cut -c1-220
echo "=== self"
callgrind_annotate $W/runs/cg_$TAG.out 2>/dev/null | head -45 | cut -c1-220
