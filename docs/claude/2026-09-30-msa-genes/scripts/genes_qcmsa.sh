#!/bin/bash
# qcmsa on every raw MSA of the congener run, as protal runs it; then the gene analysis.
S=$(cd "$(dirname "$0")" && pwd)
D=$HOME/genes_study/cong/out/strains
n=0
for msa in $D/*.raw.msa.fna; do
  sp=$(basename $msa .raw.msa.fna)
  python3 $HOME/audit5/bin/qcmsa_P2.py $msa $D/$sp.raw.partition.txt $D/$sp.meta.tsv --prefix $D/$sp --reapply-hcov 1000 \
    > $D/$sp.qcmsa.log 2>&1 || echo "qcmsa failed: $sp"
  n=$((n + 1))
done
echo "$n species"
python3 $S/msa_genes.py $HOME/genes_study/cong/out $HOME/genes_study/tuneH2
