#!/usr/bin/env bash
# The strain world's reads at other lengths: each mate cut to its first L bases (the same fragments, positions and
# errors), aligned and profiled ungated with the same binary and model as cf_runs.sh, with training dumps.
# usage: len_runs.sh L...   (150: the reads as they are)
set -eu
W=~/strainfn
CF=$W/cf
PROTAL=$CF/protal
MODEL=$W/model_s1/model_pe.xml
samples=$(cut -f1 $W/reads/samples.tsv)
for L in "$@"; do
  R=$W/run_len$L
  [ -s $R/done ] && { echo "== len $L exists"; continue; }
  rm -rf $R; mkdir -p $R/alignments $R/profiles $R/reads
  {
    printf '#OUTPUT_DIR\t%s\n#SAMPLEID\tFIRST\tSECOND\tSAM\tPREFIX\tPROFILE\tREAD_TYPE\tPROFILE_TRUTH\n' $R
    while IFS=$'\t' read -r sample f1 f2; do
      if [ $L = 150 ]; then r1=$f1; r2=$f2; else
        r1=$R/reads/${sample}_1.fq.gz; r2=$R/reads/${sample}_2.fq.gz
        for p in "$f1 $r1" "$f2 $r2"; do set -- $p
          gzip -dc $1 | perl -pe 'BEGIN { $L = shift @ARGV } if ($. % 2 == 0) { chomp; $_ = substr($_, 0, $L) . "\n" }' $L | gzip -1 > $2 &
        done
        wait
      fi
      printf '%s\t%s\t%s\t%s\t%s\t%s\tpe\t%s\n' $sample $r1 $r2 $R/alignments/$sample.sam.zst $R/$sample $R/profiles/$sample.profile $CF/truth/$sample.txt
    done < $W/reads/samples.tsv
  } > $R/samples.map
  /usr/bin/time -v -o $R/run.time taskset -c 0-3 nice -n 5 $PROTAL --db $W/db_grid --map $R/samples.map -t 4 \
      --no_strains --no_qcmsa --model $MODEL --sparse_first 0 > $R/protal.log 2>&1
  echo "== len $L: $(grep -c 'Write truth to' $R/protal.log) profiles, $(grep -E 'Elapsed' $R/run.time)"
  date > $R/done
done
