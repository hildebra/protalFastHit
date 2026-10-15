#!/bin/bash
# Re-simulate the small build's ont_b150000000 training samples with each read's origin in its name (the study's
# patched simulator; same seeds, same reads), profile them as the collector does, and trace their marker reads.
set -e
export PATH=$HOME/micromamba/envs/protal-db-build/bin:$PATH
P=$(dirname "$(readlink -f "$0")")  # this folder: the trace scripts and the simulator patch
T=~/ontu/trace/b150
R=~/ontu/rw_scratch
rm -rf $T && mkdir -p $T/reads $T/out
awk -F'\t' -v T=$T 'BEGIN{OFS="\t"} NR==1{print; next} {$2=T"/reads/"$1".fq.zst"; print}' ~/ontu/keep/ont_b150000000/samples.tsv > $T/samples.tsv
cp ~/ontu/keep/ont_b150000000/genomes.tsv $T/genomes.tsv
time PROTAL_SIM_READ_ORIGIN=1 taskset -c 0-3 nice -n 5 ~/ontu/buildpos/simulate_metagenomes \
  --long_samples $T/samples.tsv --long_genomes $T/genomes.tsv \
  --long_setup qshmm:QSHMM-ONT-HQ:8000:6000:0.97:39/24/36 --long_stats $T/stats.tsv -t 4 \
  --long_model $HOME/micromamba/envs/protal-db-build/data/QSHMM-ONT-HQ.model --genome_store $R/genome_store \
  --gene_positions $R/gene_positions.tsv --long_marker_reads $T/marker_reads.tsv > $T/simulate.log 2>&1
if cmp -s <(sort $T/marker_reads.tsv) <(sort $R/training/points/ont_b150000000/sim/marker_reads.tsv); then
  echo "marker_reads.tsv identical to the build's"
else
  echo "marker_reads.tsv DIFFERS from the build's"
fi
zstd -dc $T/reads/ont_b150000000_s_1.fq.zst | head -1
# the map: the collector's rows of these samples, with this folder's paths
{
  echo -e "#OUTPUT_DIR\t$T/out"
  echo -e "#SAMPLEID\tFIRST\tSECOND\tSAM\tPREFIX\tPROFILE\tPROFILE_TRUTH\tREAD_TYPE\tUNMAPPED_READS\tSPARSE_FIRST"
  awk -F'\t' -v T=$T 'BEGIN{OFS="\t"} $1 ~ /^ont_b150000000_s_/ {s=$1; print s, T"/reads/"s".fq.zst", "-", T"/out/"s".sam.gz", T"/out/"s, T"/out/"s".profile", $7, $8, $9, $10}' \
    $R/training/profile_all/run1/samples.map
} > $T/samples.map
grep -c ont_b150000000 $T/samples.map
time taskset -c 0-3 nice -n 5 ~/ontu/build/protal --db $R/training_db --map $T/samples.map -t 4 --no_strains --no_qcmsa \
  --sparse_index emulate > $T/protal.log 2>&1
# rows: the same as the build's?
for s in $(seq 1 12); do
  n=ont_b150000000_s_$s
  a=$(cut -f5 $T/out/$n.profile.truth_annotated | sort | md5sum | cut -c1-8)
  b=$(cut -f5 $R/training/points/ont_b150000000/protal/profiles/$n.profile.truth_annotated | sort | md5sum | cut -c1-8)
  [ "$a" = "$b" ] || echo "$n rows differ from the build's"
  awk -F'\t' 'NR>1{print $4}' $T/out/$n.profile.truth_annotated > $T/out/$n.rows
  python3 $P/trace_marker_reads.py --reads $T/reads/$n.fq.zst --sam $T/out/$n.sam.gz --genomes $T/genomes.tsv \
    --sample $n --gene-positions $R/gene_positions.tsv --rows $T/out/$n.rows --out $T/out/$n
done
echo done
