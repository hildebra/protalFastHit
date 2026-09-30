#!/bin/bash
# Profile 48 test_congeners samples of the tuning world (6 points) with the current strain-fixes build,
# strain MSAs included, into ~/genes_study/cong. Usage: genes_run.sh [threads]
set -u
T=${1:-6}
BIN=$HOME/audit5/bin/protal_P2
DB=$HOME/tune/H2/protal_db
OUT=$HOME/genes_study/cong
mkdir -p $OUT
MAP=$OUT/samples.map
{
  printf '#OUTPUT_DIR\t%s\n' "$OUT/out"
  printf '#SAMPLEID\tFIRST\tSECOND\tSAM\tPREFIX\tPROFILE\n'
  for rl in 100 150 250; do for p in 30000 100000; do
    pt=$HOME/tune/test_congeners/points/rl${rl}_p$p/sim
    awk -F'\t' -v dir="$pt/reads" '!/^#/ && NF >= 6 {print $1 "\t" dir "/" $2 "\t" dir "/" $3 "\t" $4 "\t" $5 "\t" $6}' $pt/protal.meta
  done; done
} > $MAP
echo "$(grep -vc '^#' $MAP) samples"
cat $HOME/tune/test_congeners/points/rl{100,150,250}_p{30000,100000}/sim/manifest.tsv | awk 'NR==1 || $1 != "sample"' > $OUT/manifest.tsv
/usr/bin/time -v $BIN --db $DB --map $MAP -t $T > $OUT/protal.log 2> $OUT/protal.err
echo "protal exit $?"
grep -E "Elapsed|Maximum resident" $OUT/protal.err
ls $OUT/out/strains | wc -l
