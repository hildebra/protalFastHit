#!/bin/bash
# Realigns the V2 training and test samples with secondary alignments (-m 3: the best and the two next
# candidates that protal aligns by default, align_top 3), alignment only, into ~/fpexp/<set>/.
set -eu
PROTAL=~/protal-head/build/protal
for set in test training; do
  out=~/fpexp/$set
  mkdir -p $out/sam $out/prof
  awk -F'\t' -v OFS='\t' -v O="$out" '
    /^#OUTPUT_DIR/ {print "#OUTPUT_DIR", O "/prof"; next}
    /^#/ {print; next}
    {n=split($4, p, "/"); $4 = O "/sam/" p[n]; $5 = O "/prof/" $1; $6 = O "/prof/" $1 ".profile"; print}' \
    ~/tune/V2/$set/profile_all/samples.map > $out/samples.map
  echo "$set: $(grep -vc '^#' $out/samples.map) samples"
  /usr/bin/time -v $PROTAL --db ~/tune/V2/training_db --map $out/samples.map -t 6 -m 3 --no_profile \
      > $out/protal.log 2> $out/protal.err || { echo "protal failed for $set"; tail -5 $out/protal.err; exit 1; }
  grep -E "Elapsed|Maximum resident" $out/protal.err
done
echo done
