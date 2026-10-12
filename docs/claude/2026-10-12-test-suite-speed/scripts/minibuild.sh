#!/bin/bash
# The mini DB's build with other threads and levels, alternated twice; database.protal's size.
export PATH=~/micromamba/envs/protal-db-build/bin:$PATH
W=~/testspeed/minibuild; rm -rf $W; mkdir -p $W
cd ~/testspeed/src
python3 scripts/mini_db/simulate_gtdb_release.py --outdir $W/gtdb --release 226 > /dev/null 2>&1
python3 scripts/mini_db/gtdb_to_protal_db.py --gtdb $W/gtdb --outdir $W/conv > /dev/null 2>&1
python3 scripts/mini_db/gene_neighbours.py --db $W/conv --genome_table $W/gtdb/simulation/genomes.tsv > /dev/null 2>&1
P=~/testspeed/build/protal
for round in 1 2; do
for v in "t1_l19:-t 1" "t4_l19:-t 4" "t4_l3:-t 4 --compress_level 3" "t4_l1:-t 4 --compress_level 1"; do
  name=${v%%:*}; opts=${v#*:}
  rm -rf $W/$name; cp -r $W/conv $W/$name
  s=$(date +%s.%N)
  taskset -c 0-3 nice -n 5 $P --build --no_profile $opts --db $W/$name --reference $W/$name/reference.fna --full_reference $W/$name/full_reference.fna > $W/$name.log 2>&1
  echo "$name round $round: $(echo "$(date +%s.%N) - $s" | bc | cut -c1-5) s, database.protal $(stat -c %s $W/$name/database.protal) bytes"
done
done
rm -rf $W
