#!/usr/bin/env bash
# A deep sample: 5M read pairs (2x150, ART HS25) of 150 species of the tuning world, with second and
# third strains (0.3, 0.1), simulated from the genome table of ~/tune/V2 into $BENCH/points/deep5m, then
# profiled by 0.6.0a, 0.7, and 0.7 with --depth_identity_margin 0.08 (v07dm008), as profile.sh's runs.
set -uo pipefail
B=${BENCH:-$HOME/bench07}
P6=${P6:-$HOME/protal-0.6.0a/src/build/protal} P7=${P7:-$HOME/fix-build/bin/protal}
SIM=${SIM:-$HOME/fix-build/bin/simulate_metagenomes}
T=${T:-6}
P=$B/points/deep5m R=$B/runs
if [ ! -f $P/sim/protal.meta ]; then
  rm -rf $P; mkdir -p $P
  /usr/bin/time -v -o $B/logs/deep5m.simulate.time $SIM --genome_table $HOME/tune/V2/genomes.tsv -o $P/sim -n 1 \
    --sample_prefix deep5m_s_ --total_read_pairs 5000000 --species_per_sample 150 --read_length 150 --sequencer HS25 \
    --fragment_mean 350 --fragment_stdev 50 --strains_per_species 0.3,0.1 --seed 7 -t $T --protal_metafile $P/protal \
    > $B/logs/deep5m.simulate.log 2>&1 || { echo "simulation failed"; exit 1; }
fi
run() {
  local name=$1; shift
  [ -f $R/$name.done ] && return
  rm -rf $R/$name; mkdir -p $R/$name
  echo "$(date +%T) $name"
  /usr/bin/time -v -o $R/$name.time "$@" -o $R/$name -t $T > $R/$name.log 2>&1 && touch $R/$name.done || echo "  $name failed"
}
id=deep5m_s__1  # the simulator appends _<n> to --sample_prefix
r=$P/sim/reads/$id
pair=(-1 ${r}_R1.fq.gz -2 ${r}_R2.fq.gz --prefix $id --no_qcmsa)
run v06.$id $P6 --db $B/tune/db06 "${pair[@]}"
run v07.$id $P7 --db $B/tune/db07 "${pair[@]}"
run v07dm008.$id $P7 --db $B/tune/db07 --depth_identity_margin 0.08 "${pair[@]}"
