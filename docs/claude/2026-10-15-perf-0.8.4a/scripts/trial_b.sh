#!/usr/bin/env bash
# Trial B: the profiling's SAM parse. One pe alignment of the D6331 HKU mock on the reduced r226 database (its SAM
# kept), then profile-only runs of that SAM on 4 threads, their "SAM read" lines.
set -uo pipefail
W=~/perf085
P=${PROTAL:-$W/bin/protal_base}
DB=~/Mocks/reduced_db_r226/database.protal
R=~/Mocks/reads/d6331_pe_labs/d6331_hku_SRR23926922
OUT=$W/runs/pe_hku
mkdir -p $W/runs
if [ ! -d $OUT ]; then
  taskset -c 0-3 nice -n 5 $P --db $DB -1 $R/SRR23926922_1.fq.gz -2 $R/SRR23926922_2.fq.gz -t 4 --no_strains -o $OUT > $W/runs/pe_hku.align.log 2>&1
  echo "align exit $?"
fi
SAM=$(find $OUT -name "*.sam*" | head -1)
ls -la $SAM
grep -hE "records written|Aligning reads took|Profiling sample|SAM read:" $W/runs/pe_hku.align.log | cut -c1-600
for rep in 1 2 3; do
  taskset -c 0-3 nice -n 5 $P --db $DB --profile_only $SAM --prefix hku -t 4 --no_strains --force -o $W/runs/po_$rep > $W/runs/po_$rep.log 2>&1
  echo "profile_only $rep exit $?"
  grep -hE "Profiling sample|SAM read:|Profiling took|Run protal took" $W/runs/po_$rep.log | cut -c1-700
done
