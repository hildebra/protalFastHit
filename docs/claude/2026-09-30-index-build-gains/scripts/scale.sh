#!/bin/bash
# How the phases of protal --build grow with the database: the tuning world at 25, 50 and 100% of
# its species, 56217e1 with timers (base) and with the performance branch's syncmer scan and huge
# pages (perf); base with huge pages through the glibc tunable (base_thp) at 100%.
set -e
S=$(cd "$(dirname "$0")" && pwd)  # this folder; head56.tar, syncmer.patch, huge.patch: see README
B=$HOME/protal-scale
if [ ! -x $B/perf/build/protal ]; then
  rm -rf $B; mkdir -p $B
  for v in base perf; do
    mkdir -p $B/$v/src; tar -xf $S/head56.tar -C $B/$v/src
    python3 $S/instrument.py $B/$v/src
  done
  (cd $B/perf/src && patch -p1 -s < $S/syncmer.patch && patch -p1 -s < $S/huge.patch)
  for v in base perf; do
    cmake -S $B/$v/src -B $B/$v/build -DCMAKE_BUILD_TYPE=Release > $B/$v/configure.log 2>&1
    cmake --build $B/$v/build -j 8 --target protal > $B/$v/make.log 2>&1 || { tail -30 $B/$v/make.log; exit 1; }
  done
  echo "built"
fi
W=/tmp/scale; mkdir -p $W
R=$B/base/src
if [ ! -f $W/in100/reference.fna ]; then
  cut -f2 ~/tune/release_p/*_taxonomy_r226.tsv | sed 's/.*;s__/s__/' | sort -u > $W/species.txt
  python3 - $W <<'EOF'
import random, sys
w = sys.argv[1]
species = open(f"{w}/species.txt").read().split("\n")[:-1]
random.Random(1).shuffle(species)
for keep in (50, 25):
    open(f"{w}/exclude{keep}.txt", "w").write("\n".join(species[len(species) * keep // 100:]) + "\n")
print(len(species), "species")
EOF
  for p in 100 50 25; do
    ex=""; [ $p != 100 ] && ex="--exclude_species $W/exclude$p.txt"
    python3 $R/scripts/mini_db/gtdb_to_protal_db.py --gtdb ~/tune/release_p --outdir $W/in$p -t 8 $ex 2> $W/convert$p.log
    echo "in$p: $(grep -c '>' $W/in$p/reference.fna) reference genes, $(du -m $W/in$p/reference.fna | cut -f1) MB;" \
         "$(grep -c '>' $W/in$p/full_reference.fna) full genes, $(du -m $W/in$p/full_reference.fna | cut -f1) MB"
  done
fi
run() {  # variant size threads rep
  local v=$1 p=$2 t=$3 bin=$B/${1%_thp}/build/protal env=""
  [ "${v%_thp}" != "$v" ] && env="GLIBC_TUNABLES=glibc.malloc.hugetlb=1"
  local D=$W/db; rm -rf $D; cp -r $W/in$p $D
  env $env /usr/bin/time -f "%M" -o $W/rss $bin --build --no_profile --no_compress -t $t --db $D \
      --reference $D/reference.fna --full_reference $D/full_reference.fna 2>&1 | python3 $S/stamp.py > $W/log_${v}_${p}_${t}_$4
  python3 $S/phases.py "$v p$p t$t" $W/log_${v}_${p}_${t}_$4
  echo "    rss=$(( $(cat $W/rss) / 1024 )) MB output=$(cat $D/unique_kmers.tsv $D/index.prx | md5sum | cut -c1-8)"
  rm -rf $D
}
for p in 25 50 100; do
  for t in 1 8; do
    for v in base perf; do run $v $p $t 1; done
  done
done
for t in 1 8; do
  for v in base_thp perf base base_thp; do run $v 100 $t 2; done
done
