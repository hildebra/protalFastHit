#!/bin/bash
# Kill the collector's protal run mid-way, then rerun the collector: does it resume to the same table?
exec > ~/audit6/gtdb/collect_kill.out 2>&1
G=~/audit6/gtdb
O=$G/b2
PY=~/protal-train/bin/python
S=$G/src/scripts
cp $O/training/training_data.tsv $G/kc_ref.tsv
for p in rl100_p20000 rl150_p20000; do
  find $O/training/points/$p -name '*.truth_annotated' -delete
  find $O/training/points/$p -name '*.sam.zst' -delete
done
rm -f $O/training/training_data.tsv
ARGS="--db $O/training_db --genome_table $O/genomes.tsv -o $O/training --protal $HOME/strain-build/bin/protal --simulator $HOME/strain-build/bin/simulate_metagenomes --samples 2 --read_pairs 1000,20000 --read_setups 100:HS20:300:40,150:HS25:350:50 --archaea 1 --species_per_sample 8-12 --seed 1 -t 2 --taxonomy $O/internal_taxonomy.dmp --congeners 0 --novel_species $O/heldout_species.txt --novel_clades 1"
$PY $S/collect_training_data.py $ARGS > $G/kc1.log 2>&1 &
C=$!
until ls $O/training/points/*/protal/alignments/*.partial > /dev/null 2>&1; do sleep 0.2; done
sleep 1
echo "$(date +%T) killing protal while it writes:"; ls $O/training/points/*/protal/alignments/
pkill -9 -f -- "--map $O/training/profile_all/samples.map"
wait $C; rc=$?; echo "collector exit $rc"; tail -2 $G/kc1.log
echo "left behind:"; ls $O/training/points/*/protal/alignments/ | grep -v '^$'
find $O/training/points -name '*.truth_annotated' | wc -l
echo "--- rerun"
$PY $S/collect_training_data.py $ARGS > $G/kc2.log 2>&1; echo "collector exit $?"; tail -3 $G/kc2.log
grep -c Skip $O/training/profile_all/protal.log
cmp $G/kc_ref.tsv $O/training/training_data.tsv && echo "resumed table identical to the uninterrupted run's"
