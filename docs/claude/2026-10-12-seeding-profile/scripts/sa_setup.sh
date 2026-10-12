#!/usr/bin/env bash
# The seeding profile's local world: HEAD (git archive in archives/seed_head.tar) built Release with -g (the same code,
# debug info for callgrind's line counts); every r226 archaeal species converted by HEAD's build_gtdb_database.py
# (stopped once the conversion is done), the finished database built with the pipeline's --build arguments (alleles
# from every genome, 16 per species and gene, none indexed); 1M pairs of 150 bp from 120 species with strains and
# congeners, and their first 100,000 pairs. 4 pinned cores.
set -uo pipefail
A=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/8e87f9e0-3682-453c-a7a2-636d4ec2ef61/scratchpad
W=$HOME/seed; H=$W/head
pin=(taskset -c 0-3 nice -n 5)
mkdir -p "$W"
if [ ! -x "$H/build/protal" ]; then
  rm -rf "$H"; mkdir -p "$H/src"; tar -xf "$A/archives/seed_head.tar" -C "$H/src"
  start=$(date +%s)
  "${pin[@]}" cmake -S "$H/src" -B "$H/build" -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-g" -DCMAKE_C_FLAGS="-g" > "$H/configure.log" 2>&1 &&
  "${pin[@]}" cmake --build "$H/build" -j 4 --target protal simulate_metagenomes > "$H/build.log" 2>&1 || { echo "BUILD FAIL"; grep -E "error" -A3 "$H/build.log" | head -20; exit 1; }
  echo "build OK in $(( $(date +%s) - start )) s"
fi
export PATH=$HOME/micromamba/envs/protal-db-build/bin:$PATH
if [ ! -s "$W/db/database.protal" ]; then
  rm -rf "$W/arc"
  "${pin[@]}" python3 "$H/src/scripts/build_gtdb_database.py" --gtdb "$HOME/archaea_inputs/release" --extra-genomes "$HOME/archaea_inputs/genomes" \
      --simulate-species "$HOME/archaea_inputs/simulation_species.txt" --outdir "$W/arc" --protal "$H/build/protal" \
      --simulator "$H/build/simulate_metagenomes" -t 4 --seed 1 --read-types pe --scenarios none --error-reads none \
      --evaluation basic --no-binary-check --release 226 > "$W/arc_pipeline.log" 2>&1 &
  until grep -q "5/11 training database" "$W/arc/console.log" 2> /dev/null; do sleep 5; done
  for p in $(pgrep -u "$USER" -f "build_gtdb_database.py|collect_training_data.py|seed/arc"); do
    [ "$p" = "$$" ] && continue
    grep -q "seed/" "/proc/$p/cmdline" 2> /dev/null && kill -TERM "$p"
  done
  sleep 3
  for p in $(pgrep -u "$USER" -f "seed/arc"); do [ "$p" = "$$" ] || kill -KILL "$p" 2> /dev/null; done
  sleep 1
  D=$W/arc/protal_db; rm -f "$D/database.protal"
  start=$(date +%s)
  "${pin[@]}" "$H/build/protal" --build --no_profile -t 4 --db "$D" --reference "$D/reference.fna" --compress_level 9 \
      --full_reference "$D/full_reference.fna.zst" --strain_alleles 16 --allele_genome_share 1 --index_alleles 1 > "$W/db.build.log" 2>&1 \
      || { echo "DB BUILD FAIL"; tail "$W/db.build.log"; exit 1; }
  echo "database built in $(( $(date +%s) - start )) s"
  mkdir -p "$W/db"; mv "$D/database.protal" "$W/db/"
fi
S=$W/samples
if [ ! -s "$S/pe100k_R2.fq.gz" ]; then
  rm -rf "$S"; mkdir -p "$S"
  "${pin[@]}" "$H/build/simulate_metagenomes" --genome_table "$W/arc/work/genomes.tsv" -o "$S/pe1m" -n 1 --sample_prefix pe1m --total_read_pairs 1000000 \
      --species_per_sample 120 --strains_per_species 0.3 --congener_groups 0.25:2-5 --pln_sigma 1.5 --read_length 150 --sequencer HSXt --seed 11 -t 4 > "$S/pe1m.log" 2>&1
  for e in 1 2; do zcat "$S/pe1m/reads/pe1m_1_R${e}.fq.gz" | head -400000 | gzip -1 > "$S/pe100k_R$e.fq.gz"; done
fi
ls -la "$W/db" "$S"
echo "SETUP DONE"
