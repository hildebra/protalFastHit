# --discard-constant output given to IQ-TREE as the justfile does (GTR+G) and with +ASC.
cd ~/audit5/qc
S=~/audit5/smoke/protal/strains; sp=s__Mockella_gamma
mkdir -p asc
python3 ~/audit5/src/scripts/qcmsa.py $S/$sp.raw.msa.fna $S/$sp.raw.partition.txt $S/$sp.meta.tsv --prefix asc/$sp.dc \
   --reapply-hcov 1000 --discard-constant --min-parsimony-samples 0 2>&1 | grep "Site cleanup"
for m in GTR+G GTR+G+ASC; do
  iqtree2 -s asc/$sp.dc.msa.fna -m $m -B 1000 -T 2 --seqtype DNA --prefix asc/$sp.$m -redo > asc/$sp.$m.log 2>&1
  echo "$m rc=$?"; grep -a "ERROR\|Total tree length\|constant sites" asc/$sp.$m.log asc/$sp.$m.iqtree 2>/dev/null | head -4
done
