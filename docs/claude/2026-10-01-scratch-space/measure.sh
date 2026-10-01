#!/bin/bash
# Bytes per read pair (ART, the three default read setups) and per base (pbsim3, the default PacBio and
# Nanopore setups) of what collect_training_data.py writes. Run in WSL on 2026-10-01; see README.md.
W=~/progress-test/measure
mkdir -p "$W"
SIM=~/progress-test/bin/simulate_metagenomes         # copy of the build of c575aa5's tree
GT=~/protal-perf/world900/simulation/genomes.tsv     # the synthetic 900-species world
ENV=/home/falk/micromamba/envs/protal-db-build/bin   # pbsim3 (envs/protal-db-build.yaml)

for setup in 100:HS20:300:40 150:HSXt:350:50 250:MSv3:550:50; do
    IFS=: read L P FM FS <<< "$setup"
    out="$W/pe_$L"
    rm -rf "$out"
    nice -n 10 "$SIM" --genome_table "$GT" -o "$out" -n 1 --sample_prefix s --total_read_pairs 50000 --species_per_sample 10 \
        --read_length "$L" --sequencer "$P" --fragment_mean "$FM" --fragment_stdev "$FS" --seed 3 -t 2 > "$out.log" 2>&1
    gz=$(cat "$out"/reads/*_R?.fq.gz | wc -c)
    plain=$(zcat "$out"/reads/*_R?.fq.gz | wc -c)
    pairs=$(( $(zcat "$out"/reads/*_R1.fq.gz | wc -l) / 4 ))
    echo "pe $setup: $pairs pairs, gz $gz bytes ($(( gz / pairs )) per pair), plain $plain ($(( plain / pairs )) per pair)"
done

MODELS=$(dirname "$(find /home/falk/micromamba/envs/protal-db-build -name "QSHMM-ONT-HQ.model" | head -1)")
G=$(cut -f3 "$GT" | sed -n 2p)   # the first genome (line 1 is the header)
zcat -f "$G" > "$W/g.fna"
for spec in "errhmm ERRHMM-SEQUEL 15000 3000 0.999 -" "qshmm QSHMM-ONT-HQ 8000 6000 0.97 39:24:36"; do
    read M MOD LM LS ACC R <<< "$spec"
    d="$W/lr_$M"; rm -rf "$d"; mkdir -p "$d"
    extra=(); [ "$R" != "-" ] && extra=(--difference-ratio "$R")
    ( cd "$d" && nice -n 10 "$ENV/pbsim" --strategy wgs --method "$M" --"$M" "$MODELS/$MOD.model" --genome "$W/g.fna" \
        --depth 20 --length-mean "$LM" --length-sd "$LS" --accuracy-mean "$ACC" --seed 1 --prefix "$d/r" "${extra[@]}" > "$d/log" 2>&1 )
    bases=$(zcat -f "$d"/r_*.fq.gz | awk 'NR%4==2{n+=length($0)} END{print n}')
    gz1=$(zcat -f "$d"/r_*.fq.gz | gzip -1 | wc -c)   # the collector merges a sample's reads at gzip level 1
    echo "$M: $bases bases; pbsim fq.gz $(cat "$d"/*.fq.gz | wc -c), maf.gz $(cat "$d"/*.maf.gz | wc -c), ref $(cat "$d"/*.ref | wc -c) bytes (genome $(wc -c < "$W/g.fna")); merged gzip -1 $gz1 bytes"
done
