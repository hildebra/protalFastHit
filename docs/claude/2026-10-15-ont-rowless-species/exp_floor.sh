#!/bin/bash
# The small build's ont_b150000000 communities at 6, 30 and 150 Mb (same seeds), simulated with each read's origin and
# errors in its name; profiled on the training database as the collector does, at protal's ONT identity floor (-a)
# 0.85 (the default), 0.80 and 0.75; their marker reads traced. Cores 0-3, niced.
set -e
export PATH=$HOME/micromamba/envs/protal-db-build/bin:$PATH
P=$(dirname "$(readlink -f "$0")")  # this folder: the trace scripts and the simulator patch
R=~/ontu/rw_scratch
# ~/ontu/srcpos: the source at 55f850d; the study's simulator patch applied once (skipped if it is)
patch -d ~/ontu/srcpos -p1 --forward --reject-file=- < $P/read_origin_names.patch || true
touch ~/ontu/srcpos/src/RandomForest/LongReadSimulator.cpp
taskset -c 0-3 nice -n 5 ninja -C ~/ontu/buildpos -j4 simulate_metagenomes > ~/ontu/ninjapos2.log 2>&1 || { tail -20 ~/ontu/ninjapos2.log; exit 1; }
for D in 6000000 30000000 150000000; do
  T=~/ontu/trace/d$D
  rm -rf $T && mkdir -p $T/reads
  awk -F'\t' -v T=$T -v D=$D 'BEGIN{OFS="\t"} NR==1{print; next} {sub(/ont_b150000000/, "ont_b" D, $1); $2=T"/reads/"$1".fq.zst"; $3=D; print}' \
    ~/ontu/keep/ont_b150000000/samples.tsv > $T/samples.tsv
  awk -F'\t' -v D=$D 'BEGIN{OFS="\t"} NR==1{print; next} {sub(/ont_b150000000/, "ont_b" D, $1); print}' \
    ~/ontu/keep/ont_b150000000/genomes.tsv > $T/genomes.tsv
  PROTAL_SIM_READ_ORIGIN=1 taskset -c 0-3 nice -n 5 ~/ontu/buildpos/simulate_metagenomes \
    --long_samples $T/samples.tsv --long_genomes $T/genomes.tsv \
    --long_setup qshmm:QSHMM-ONT-HQ:8000:6000:0.97:39/24/36 --long_stats $T/stats.tsv -t 4 \
    --long_model $HOME/micromamba/envs/protal-db-build/data/QSHMM-ONT-HQ.model --genome_store $R/genome_store \
    --gene_positions $R/gene_positions.tsv --long_marker_reads $T/marker_reads.tsv > $T/simulate.log 2>&1
  for A in 0.85 0.80 0.75; do
    O=$T/out_a$A
    mkdir -p $O
    {
      echo -e "#OUTPUT_DIR\t$O"
      echo -e "#SAMPLEID\tFIRST\tSECOND\tSAM\tPREFIX\tPROFILE\tPROFILE_TRUTH\tREAD_TYPE\tUNMAPPED_READS\tSPARSE_FIRST"
      awk -F'\t' -v T=$T -v O=$O -v D=$D 'BEGIN{OFS="\t"} $1 ~ /^ont_b150000000_s_/ {s=$1; sub(/ont_b150000000/, "ont_b" D, s);
        print s, T"/reads/"s".fq.zst", "-", O"/"s".sam.gz", O"/"s, O"/"s".profile", $7, $8, $9, $10}' $R/training/profile_all/run1/samples.map
    } > $O/samples.map
    /usr/bin/time -f "protal -a $A d$D: %e s, %M KB" taskset -c 0-3 nice -n 5 ~/ontu/build/protal --db $R/training_db \
      --map $O/samples.map -t 4 --no_strains --no_qcmsa --sparse_index emulate -a $A > $O/protal.log 2>$O/time.txt
    cat $O/time.txt | tail -1
    for f in $O/*.profile.truth_annotated; do
      n=$(basename $f .profile.truth_annotated)
      awk -F'\t' 'NR>1{print $4}' $f > $O/$n.rows
      taskset -c 0-3 nice -n 5 python3 $P/trace_marker_reads.py --reads $T/reads/$n.fq.zst --sam $O/$n.sam.gz \
        --genomes $T/genomes.tsv --sample $n --gene-positions $R/gene_positions.tsv --rows $O/$n.rows \
        --heldout ~/ontu/rw/model_logs/heldout_species.txt --taxonomy ~/ontu/rw/work/internal_taxonomy.dmp \
        --out $O/$n > /dev/null
    done
    # rows: present and absent taxa (truth column 1)
    awk -F'\t' 'FNR>1{n++; if($1==1)p++; else a++} END{printf "rows d%s a%s: %d (present %d, absent %d)\n", D, A, n, p, a}' \
      D=$D A=$A $O/*.profile.truth_annotated
    python3 $P/aggregate_marker_reads.py $O/*.reads.tsv.gz > $O/tables.txt
    grep -E "^### species with|species with one|species with two" $O/tables.txt
  done
done
echo done
