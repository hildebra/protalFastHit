# usage: simrun.sh <name> <simulate args...>
name=$1; shift
D=~/audit5/qc/runs/$name
mkdir -p $D
cd $D
~/audit5/bin/simulate_metagenomes --genome_table ~/audit5/world/gtdb_r226/simulation/genomes.tsv \
  -o $D/sim --protal_metafile $D/protal -t 2 "$@" > sim.log 2>&1 || { echo "sim failed"; tail sim.log; exit 1; }
PROTAL_QCMSA_SCRIPT=$HOME/audit5/src/scripts/qcmsa.py ~/audit5/bin/protal --db ~/audit5/world/protal_db \
  --map $D/sim/protal.meta -t 2 > protal.log 2>&1
echo "protal rc=$?"
tr '\r' '\n' < protal.log | grep -a "qcmsa\|WARNING\|ERROR\|Error" | grep -v "^\[qcmsa\] python3\|^\[qcmsa\] /" | head -40
ls $D/protal/strains | head -80
