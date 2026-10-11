#!/usr/bin/env bash
# The performance comparison's world, after the identity runs: every r226 archaeal species converted by HEAD's
# build_gtdb_database.py (stopped once the conversion is done), the finished database built by HEAD's protal with the
# pipeline's --build arguments (db_ref: text tables only), a copy given the binary forms by the changed protal
# (--compress_db: db_bin); samples from the 3,259 genomes. 4 pinned cores.
set -uo pipefail
A=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/8e87f9e0-3682-453c-a7a2-636d4ec2ef61/scratchpad
W=$HOME/impl/perf; REF=$HOME/impl/ref; WORK=$HOME/impl/work
until grep -q "vs cohort_work_bin_ahead" "$A/impl_identity.out" 2> /dev/null; do sleep 15; done
export PATH=$HOME/micromamba/envs/protal-db-build/bin:$PATH
pin=(taskset -c 0-3 nice -n 5)
mkdir -p "$W"
if [ ! -s "$W/db_ref/database.protal" ]; then
  rm -rf "$W/arc"
  "${pin[@]}" python3 "$REF/src/scripts/build_gtdb_database.py" --gtdb "$HOME/archaea_inputs/release" --extra-genomes "$HOME/archaea_inputs/genomes" \
      --simulate-species "$HOME/archaea_inputs/simulation_species.txt" --outdir "$W/arc" --protal "$REF/build/protal" \
      --simulator "$REF/build/simulate_metagenomes" -t 4 --seed 1 --read-types pe,pb --scenarios none --error-reads none \
      --evaluation basic --no-binary-check --release 226 > "$W/arc_pipeline.log" 2>&1 &
  until grep -q "5/11 training database" "$W/arc/console.log" 2> /dev/null; do sleep 5; done
  me=$$
  for p in $(pgrep -u "$USER" -f "impl/perf/arc|build_gtdb_database.py|collect_training_data.py"); do
    [ "$p" = "$me" ] && continue
    grep -q "impl/perf" "/proc/$p/cmdline" 2> /dev/null && kill -TERM "$p"
  done
  sleep 3
  for p in $(pgrep -u "$USER" -f "impl/perf/arc"); do [ "$p" = "$me" ] || kill -KILL "$p" 2> /dev/null; done
  sleep 1
  D=$W/arc/protal_db
  rm -f "$D/database.protal"
  /usr/bin/time -f "%e %M" -o "$W/db_ref.time" "${pin[@]}" "$REF/build/protal" --build --no_profile -t 4 --db "$D" --reference "$D/reference.fna" \
      --compress_level 9 --full_reference "$D/full_reference.fna.zst" --strain_alleles 16 --allele_genome_share 1 --index_alleles 1 > "$W/db_ref.build.log" 2>&1 \
      || { echo "build FAIL"; tail "$W/db_ref.build.log"; exit 1; }
  mkdir -p "$W/db_ref" "$W/db_bin"
  mv "$D/database.protal" "$W/db_ref/"
  cp "$W/db_ref/database.protal" "$W/db_bin/"
  M=$HOME/realworld_s1/model_logs
  for t in pe:trained_model.xml se:trained_model.xml pb:trained_model_pb.xml ont:trained_model_pb.xml; do
    for d in db_ref db_bin; do "$REF/build/protal" --add_model "$M/${t#*:}" --read_type "${t%%:*}" --db "$W/$d" > /dev/null 2>&1; done
  done
  /usr/bin/time -f "%e %M" -o "$W/db_bin.time" "${pin[@]}" "$WORK/build/protal" --compress_db --db "$W/db_bin/database.protal" -t 4 > "$W/db_bin.log" 2>&1 \
      || { echo "compress FAIL"; tail "$W/db_bin.log"; exit 1; }
  echo "db_ref built in $(cat "$W/db_ref.time"); db_bin in $(cat "$W/db_bin.time")"; grep -E "^Wrote|^Rewrite" "$W/db_bin.log"
  ls -la "$W/db_ref" "$W/db_bin"
fi
S=$W/samples
if [ ! -s "$S/long/pb.fq.gz" ]; then
  SIM=$REF/build/simulate_metagenomes T=$W/arc/work/genomes.tsv
  rm -rf "$S"; mkdir -p "$S/long"
  "${pin[@]}" "$SIM" --genome_table "$T" -o "$S/pe1m" -n 1 --sample_prefix pe1m --total_read_pairs 1000000 --species_per_sample 120 \
      --strains_per_species 0.3 --congener_groups 0.25:2-5 --pln_sigma 1.5 --read_length 150 --sequencer HSXt --seed 11 -t 4 > "$S/pe1m.log" 2>&1
  "${pin[@]}" "$SIM" --genome_table "$T" -o "$S/cohort" -n 24 --sample_prefix c --total_read_pairs 50000 --species_per_sample 30-150 \
      --strains_per_species 0.3 --congener_groups 0.25:2-5 --pln_sigma 1.5 --read_length 150 --sequencer HSXt --seed 21 -t 4 > "$S/cohort.log" 2>&1
  M=$S/pe1m/manifests/pe1m_1.tsv
  printf 'sample\tout\tbases\tseed\npb\t%s\t100000000\t31\n' "$S/long/pb.fq.gz" > "$S/long/pb_samples.tsv"
  { printf 'sample\tgenome\tfasta\tweight\thost\n'; awk -F'\t' 'NR > 1 { print "pb\t" $2 "\t" $11 "\t" $6 "\t0" }' "$M"; } > "$S/long/pb_genomes.tsv"
  "${pin[@]}" "$SIM" --long_samples "$S/long/pb_samples.tsv" --long_genomes "$S/long/pb_genomes.tsv" --long_setup hifi:15000:3000:3 -t 4 > "$S/pb.log" 2>&1
fi
ls -la "$S" "$S/long"
echo SETUP DONE
