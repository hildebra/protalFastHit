#!/usr/bin/env bash
# Re-profile an existing run's alignments with relaxed / changed filters (protal reuses existing SAMs).
# Usage: run_variants.sh <run> <variant> [<variant> ...]
#   variants: nostrand mincov1 mincov1_nostrand noqual af10 maxall1 qc_nopars qc_nogenegate qc_nomrate2 qc_all
set -euo pipefail
ACC=$HOME/audit5/accuracy
DB=$HOME/audit5/world/protal_db
export PROTAL_QCMSA_SCRIPT=$HOME/audit5/src/scripts/qcmsa.py
run=$1; shift
for v in "$@"; do
  case $v in
    nostrand)         args=(--snp_no_strand);;
    mincov1)          args=(--snp_min_cov 1);;
    mincov1_nostrand) args=(--snp_min_cov 1 --snp_no_strand);;
    noqual)           args=(--snp_min_phred_sum 0 --snp_min_mean_qual 0);;
    af10)             args=(--snp_min_af 0.1);;
    maxall1)          args=(--snp_max_alleles 1);;
    qc_nopars)        args=(--qcmsa_args "--min-parsimony-samples 0");;
    qc_nogenegate)    args=(--qcmsa_args "--gene-min-samples 0 --gene-min-hcov 0 --gene-min-mean-depth 0");;
    qc_nomrate2)      args=(--qcmsa_args "--min-bad 1000000 --no-mask-cell-outliers");;
    combo)            args=(--snp_min_cov 1 --snp_no_strand --snp_min_af 0.15);;
    qc_all)           args=(--qcmsa_args "--min-parsimony-samples 0 --gene-min-samples 0 --gene-min-hcov 0 --gene-min-mean-depth 0 --min-bad 1000000 --no-mask-cell-outliers");;
    *) echo "unknown variant $v"; exit 1;;
  esac
  out=$ACC/prot_${run}_$v
  if [ -d "$out/strains" ]; then continue; fi
  mkdir -p "$out/alignments"
  for f in $ACC/prot_$run/alignments/*.sam.gz; do ln -sf "$f" "$out/alignments/"; done
  /usr/bin/time -f "protal ${run}_$v: %e s" $HOME/audit5/bin/protal --db $DB --map $ACC/sim_$run/protal.meta -o "$out" -t 2 "${args[@]}" \
      > $ACC/logs/protal_${run}_$v.log 2>&1 || echo "protal ${run}_$v exit $?"
done
