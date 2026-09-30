#!/bin/bash
# Strain audit setup: build protal at the committed HEAD in ~/audit5/src and a synthetic strain
# world (8 species x 6 genomes; the DB holds each species' representative) in ~/audit5/world.
set -euo pipefail
A=$HOME/audit5
SRC=/mnt/c/Users/hildebra/Documents/locDev/protal
mkdir -p $A/bin
rm -rf $A/src; mkdir -p $A/src
git -C $SRC archive HEAD | tar -x -C $A/src
git -C $SRC rev-parse --short HEAD > $A/COMMIT
cd $A/src
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DPROTAL_BUILD_TESTS=ON > $A/cmake.log 2>&1
cmake --build build --target protal simulate_metagenomes protal_tests -j 8 > $A/build.log 2>&1
cp build/protal build/simulate_metagenomes $A/bin/
ctest --test-dir build 2>&1 | tail -2
cat > $A/lineages.txt <<'L'
d__Bacteria;p__Simulatota;c__Simulatia;o__Simulales;f__Simulaceae;g__Mockella;s__Mockella alpha
d__Bacteria;p__Simulatota;c__Simulatia;o__Simulales;f__Simulaceae;g__Mockella;s__Mockella beta
d__Bacteria;p__Simulatota;c__Simulatia;o__Simulales;f__Simulaceae;g__Mockella;s__Mockella gamma
d__Bacteria;p__Fictota;c__Fictia;o__Fictales;f__Fictaceae;g__Fakibacter;s__Fakibacter gamma
d__Bacteria;p__Fictota;c__Fictia;o__Fictales;f__Fictaceae;g__Fakibacter;s__Fakibacter delta
d__Bacteria;p__Mockota;c__Mockia;o__Mockales;f__Mockaceae;g__Testella;s__Testella one
d__Bacteria;p__Solota;c__Solia;o__Solales;f__Solaceae;g__Dummya;s__Dummya solo
d__Archaea;p__Thermosimota;c__Thermosimia;o__Thermosimales;f__Thermosimaceae;g__Calidella;s__Calidella fervens
L
W=$A/world; rm -rf $W; mkdir -p $W
python3 scripts/mini_db/simulate_gtdb_release.py --outdir $W/gtdb_r226 --release 226 --seed 21 \
  --lineages $A/lineages.txt --genomes_per_species 6 --strain_divergence 0.001-0.01 > $W/release.log 2>&1
python3 scripts/mini_db/gtdb_to_protal_db.py --gtdb $W/gtdb_r226 --outdir $W/protal_db > $W/convert.log 2>&1
$A/bin/protal --build --no_profile -t 4 --db $W/protal_db --reference $W/protal_db/reference.fna \
  --full_reference $W/protal_db/full_reference.fna > $W/build.log 2>&1
ls -la $W/protal_db
echo "SETUP DONE $(cat $A/COMMIT)"
