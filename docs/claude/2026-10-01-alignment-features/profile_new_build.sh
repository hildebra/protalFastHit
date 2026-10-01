#!/bin/bash
# Aligns and profiles the training and test samples of a build_gtdb_database.py run (default ~/tune/V2) with a new
# protal build against that run's training database, as collect_training_data.py profiles them (--no_strains
# --no_qcmsa, each sample with its READ_TYPE and truth), into $OUT/<set>/: SAMs, profiles and training dumps
# (<sample>.profile.truth_annotated).
#     bash profile_new_build.sh PROTAL OUT [RUN]
set -eu
PROTAL=$1
OUT=$2
RUN=${3:-$HOME/tune/V2}
for set in test training; do
  out=$OUT/$set
  mkdir -p $out/sam $out/prof
  awk -F'\t' -v OFS='\t' -v O="$out" '
    /^#OUTPUT_DIR/ {print "#OUTPUT_DIR", O "/prof"; next}
    /^#/ {print; next}
    {n=split($4, p, "/"); $4 = O "/sam/" p[n]; $5 = O "/prof/" $1; $6 = O "/prof/" $1 ".profile"; print}' \
    $RUN/$set/profile_all/samples.map > $out/samples.map
  echo "$set: $(grep -vc '^#' $out/samples.map) samples"
  /usr/bin/time -v $PROTAL --db $RUN/training_db --map $out/samples.map -t 6 --no_strains --no_qcmsa \
      > $out/protal.log 2> $out/protal.err || { echo "protal failed for $set"; tail -20 $out/protal.err; exit 1; }
  grep -E "Elapsed|Maximum resident" $out/protal.err
  echo "dumps: $(ls $out/prof/*.truth_annotated 2>/dev/null | wc -l)"
done
