#!/usr/bin/env bash
# evidence_size.sh RUN SAMPLE - what one sample keeps for the strain MSAs, from a run's outputs
# (docs/claude/2026-10-07-sam-combine): the species it enters MSAs of, their genes' bases (one coverage value each), and
# the positions with variants (bins), from misc/<species>.snps_total.tsv (positions with a valid allele) and
# .snps_filtered.tsv (multi-allelic plus positions without one) less .snps_multiallelic.tsv.
set -u
RUN=$1; S=${2:-s1}
awk -F'\t' -v s="$S" '$1 == s { n++; bases += $14 } END { printf "species rows %d, gene bases %d (coverage at 4 B: %.1f MB)\n", n, bases, bases * 4 / 1e6 }' \
  <(cat "$RUN"/strains/*.meta.tsv)
total=$(cat "$RUN"/misc/*.snps_total.tsv | awk -F'\t' -v s="$S" '$1 == s { t += $2 } END { print t + 0 }')
noisy=$(cat "$RUN"/misc/*.snps_filtered.tsv | awk -F'\t' -v s="$S" '$1 == s { t += $2 } END { print t + 0 }')
multi=$(cat "$RUN"/misc/*.snps_multiallelic.tsv | awk -F'\t' -v s="$S" '$1 == s { t += $2 } END { print t + 0 }')
echo "positions with variants (all reads): $((total + noisy - multi)); with a valid allele $total, multi-allelic $multi"
echo "species with MSAs: $(ls "$RUN"/strains/*.raw.msa.fna 2>/dev/null | wc -l)"
