#!/usr/bin/env bash
# samdiff.pl over all samples of the strain world, g0 against g4, summed; and the primary mapped records of each run.
S=$(cd "$(dirname "$0")" && pwd)
cd ~/strainfn
for s in $(cut -f1 reads/samples.tsv); do perl $S/samdiff.pl run_g0/alignments/$s.sam.zst run_g4/alignments/$s.sam.zst; done |
perl -ne 'while (/(\w+)=(\d+)(?:\/(\d+))?/g) { if (defined $3) { $a{$1} += $2; $b{$1} += $3 } else { $n{$1} += $2 } }
  END { print join(" ", map { "$_=$n{$_}" } sort keys %n), "\n"; printf "%s: %d -> %d (%.3f)\n", $_, $a{$_}, $b{$_}, $b{$_}/$a{$_} for sort keys %a }'
for r in g0 g4; do
  printf "%s primary mapped records: " $r
  for s in $(cut -f1 reads/samples.tsv); do zstdcat run_$r/alignments/$s.sam.zst | grep -v '^@' | awk -F'\t' 'and($2,4)==0 && and($2,256)==0' | wc -l; done | paste -sd+ | bc
done
