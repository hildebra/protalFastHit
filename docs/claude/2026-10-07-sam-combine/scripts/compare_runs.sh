#!/usr/bin/env bash
# compare_runs.sh A B - compares two runs' outputs file by file (docs/claude/2026-10-07-sam-combine): the strain
# folder's MSAs, partitions, meta and SNP tables, haplotypes and species list, misc's per-species SNP and coverage
# tables, and the profiles. Prints the files that differ and a count.
set -u
A=$1; B=$2
n=0; d=0
while IFS= read -r f; do
  n=$((n + 1))
  if ! cmp -s "$A/$f" "$B/$f"; then d=$((d + 1)); echo "  differs: $f"; fi
done < <(cd "$A" && find . -type f \( -path './strains/*' -o -name '*.snps_*.tsv' -o -name '*.hcov.tsv' -o -name '*.profile*' \) | sort)
missing=$(cd "$B" && find . -type f \( -path './strains/*' -o -name '*.snps_*.tsv' -o -name '*.hcov.tsv' -o -name '*.profile*' \) | sort |
          comm -13 <(cd "$A" && find . -type f \( -path './strains/*' -o -name '*.snps_*.tsv' -o -name '*.hcov.tsv' -o -name '*.profile*' \) | sort) - | wc -l)
echo "$A vs $B: $n files, $d differ, $missing only in the second"
