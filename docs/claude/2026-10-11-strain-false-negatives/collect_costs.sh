#!/bin/bash
# Each variant's cost: alleles and edits kept, indexed allele seeds, the reference's index k-mers, database.protal size,
# build and run peak memory (/usr/bin/time), into costs.tsv next to this script.
HERE=$(cd "$(dirname "$0")" && pwd)
S=~/strainfn
{
  printf 'variant\talleles\tedits\tallele_seeds\treference_kmers\tdb_bytes\tbuild_peak_kb\tbuild_wall\n'
  for d in $S/db_*; do
    v=${d##*/db_}
    al=$(grep -m1 "^Strain alleles: [0-9]" $d/build.log | sed -E 's/^Strain alleles: ([0-9]+) alleles \(([0-9]+) edits\).*/\1\t\2/')
    [ -n "$al" ] || al=$'0\t0'
    seeds=$(grep -m1 "^Index alleles: [0-9]" $d/build.log | sed -E 's/^Index alleles: ([0-9]+) seeds.*/\1/')
    km=$(grep -m1 "^Index cores: " $d/build.log | sed -E 's/^Index cores: ([0-9]+) k-mers.*/\1/')
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' $v "$al" ${seeds:-0} $km $(stat -c %s $d/database.protal) \
        $(grep "Maximum resident" $d/build.time | awk '{print $NF}') $(grep "Elapsed" $d/build.time | awk '{print $NF}')
  done
} > $HERE/costs.tsv
{
  printf 'run\trun_peak_kb\trun_wall\n'
  for r in $S/run_*; do
    printf '%s\t%s\t%s\n' ${r##*/run_} $(grep "Maximum resident" $r/run.time | awk '{print $NF}') $(grep "Elapsed" $r/run.time | awk '{print $NF}')
  done
} > $HERE/run_costs.tsv
cat $HERE/costs.tsv $HERE/run_costs.tsv
