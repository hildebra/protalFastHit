#!/usr/bin/env bash
# A four-sample map run (500k pairs each) at 6 threads, strains on: times the cohort stages (profiling, strain MSAs).
set -uo pipefail
W=$HOME/mt-work/perf4; B=$W/ref/build/protal; DB=$HOME/bench071/V073/protal_db
P=$HOME/bench071/samples/points/rl150_p500000/sim/reads
while ! grep -q "Run protal took" $W/o.prof5M_t6b.log 2>/dev/null; do sleep 5; done
rm -rf $W/o.map4; mkdir -p $W/o.map4
{ printf '#OUTPUT_DIR\t%s\n#INPUT_DIR\t%s\n#SAMPLEID\tFIRST\tSECOND\tPREFIX\n' $W/o.map4 $P
  for i in 1 2 3 4; do printf 's%d\trl150_p500000_s_%d_R1.fq.gz\trl150_p500000_s_%d_R2.fq.gz\ts%d\n' $i $i $i $i; done; } > $W/map4.tsv
/usr/bin/time -v $B --db $DB --map $W/map4.tsv --read_type pe -t 6 --no_qcmsa --verbose > $W/o.map4.log 2> $W/o.map4.time
echo "== map4 t=6: $(grep -E 'Elapsed|User time|Maximum resident' $W/o.map4.time | tr -s ' ' | tr '\n' ';')"
grep -E "^(Load Index|Aligning reads|Processing all samples|Profiling|Strain-level MSAs|Run protal) took|Profile sample took" $W/o.map4.log
ls $W/o.map4/strains | head; echo "strain species: $(($(wc -l < $W/o.map4/strains/species.tsv)-1))"
