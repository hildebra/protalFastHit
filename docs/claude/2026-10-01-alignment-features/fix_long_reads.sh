#!/bin/bash
# After the SettleByRead fix (a gene whose best hit is clearly another taxon's keeps it, with MAPQ 0): a copy of the
# consensus build's outputs (~/fpexp/cons) with the PacBio and ONT samples aligned and profiled again by the fixed
# binary, the tables made again, and the long-read comparison run again (paired-end and single-end are unchanged).
#     bash fix_long_reads.sh [PROTAL]
set -eu
export PATH=$HOME/micromamba/envs/protal-db-build/bin:$PATH
D=$(cd "$(dirname "$0")" && pwd)
PROTAL=${1:-$HOME/protal-cons/build/protal}
RUN=$HOME/tune/V3
F=$HOME/fpexp/fix
rm -rf $F
cp -a $HOME/fpexp/cons $F
for set in test training; do
  out=$F/$set
  sed -i "s#/fpexp/cons/#/fpexp/fix/#g" $out/samples.map
  awk -F'\t' '/^#/ || $NF == "pb" || $NF == "ont"' $out/samples.map > $out/samples_long.map
  for s in $(awk -F'\t' '!/^#/ {print $1}' $out/samples_long.map); do
    rm -rf $out/sam/$s.* $out/prof/$s.* $out/prof/*/$s.* $out/prof/*/${s}_*
  done
  echo "$set: $(grep -vc '^#' $out/samples_long.map) long-read samples"
  /usr/bin/time -v $PROTAL --db $RUN/training_db --map $out/samples_long.map -t 6 --no_strains --no_qcmsa \
      > $out/protal_long.log 2> $out/protal_long.err || { echo "protal failed for $set"; tail -20 $out/protal_long.err; exit 1; }
  grep -E "Elapsed|Maximum resident" $out/protal_long.err
  grep -h "gene hits fit several taxa" $out/protal_long.log
done
cd $D
V2=$RUN python3 tables_new_build.py $F
READ_TYPES="pb ont" python3 run_new_build.py $F 5
