#!/bin/bash
# Align and profile the simulated samples (simulate_strain_reads.py's samples.tsv) against database variants, one
# protal run per variant (the database loaded once), with the small real world's paired-end model (seed 1, 371ac84) so
# that every variant is judged by the same model. Peak memory by /usr/bin/time.
# usage: run_samples.sh <protal binary> <run>... (run = name:database variant:extra protal flags, comma-separated)
set -eu
PROTAL=$1; shift
S=~/strainfn
MODEL=$S/model_s1/model_pe.xml
if [ ! -s $MODEL ]; then
  mkdir -p $S/model_s1
  $PROTAL --unpack_db --db ~/realworld_sp_s1/protal_db --unpack_dir $S/model_s1 > $S/model_s1/unpack.log 2>&1
fi
for spec in "$@"; do
  IFS=: read -r name db flags <<< "$spec"
  R=$S/run_$name
  [ -s $R/done ] && { echo "== $name exists"; continue; }
  rm -rf $R; mkdir -p $R/alignments $R/profiles
  {
    printf '#OUTPUT_DIR\t%s\n#SAMPLEID\tFIRST\tSECOND\tSAM\tPREFIX\tPROFILE\tREAD_TYPE\n' $R
    while IFS=$'\t' read -r sample f1 f2; do
      printf '%s\t%s\t%s\t%s\t%s\t%s\tpe\n' $sample $f1 $f2 $R/alignments/$sample.sam.zst $R/$sample $R/profiles/$sample.profile
    done < $S/reads/samples.tsv
  } > $R/samples.map
  /usr/bin/time -v -o $R/run.time taskset -c 0-3 nice -n 5 $PROTAL --db $S/db_$db --map $R/samples.map -t 4 \
      --no_strains --no_qcmsa --model $MODEL ${flags//,/ } > $R/protal.log 2>&1
  echo "== $name ($db ${flags//,/ }): $(grep -E 'Maximum resident' $R/run.time) $(grep -E 'Elapsed' $R/run.time)"
  date > $R/done
done
echo "== RUNS DONE $(date +%T)"
