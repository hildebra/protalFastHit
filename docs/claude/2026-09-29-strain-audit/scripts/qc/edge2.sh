set -u
Q=~/audit5/qc; SRC=$Q/runs/lowcov
export PROTAL_QCMSA_SCRIPT=$HOME/audit5/src/scripts/qcmsa.py
rm -rf $Q/runs/e7; mkdir -p $Q/runs/e7; cp -r $SRC/protal $Q/runs/e7/protal
sed "s#^\#OUTPUT_DIR.*#\#OUTPUT_DIR\t$Q/runs/e7/protal#" $SRC/sim/protal.meta > $Q/runs/e7/map.tsv
touch -d '2020-01-01' $Q/runs/e7/protal/strains/*
~/audit5/bin/protal --db ~/audit5/world/protal_db --map $Q/runs/e7/map.tsv -t 2 --map_range 4-6 > $Q/runs/e7/r.log 2>&1
echo "rc=$?"
tr '\r' '\n' < $Q/runs/e7/r.log | grep -a "across samples"
echo "files still dated 2020 (not rewritten by this run):"
ls -l --time-style=+%Y $Q/runs/e7/protal/strains/ | awk '$6==2020{print $7}' | sed 's/\..*//' | sort | uniq -c
grep -c '>' $Q/runs/e7/protal/strains/s__Mockella_alpha.msa.fna
