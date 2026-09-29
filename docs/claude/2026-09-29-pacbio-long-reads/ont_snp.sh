#!/usr/bin/env bash
# SNP noise of long reads: two samples each of HiFi and ONT reads of the mock community (single
# strains), strain MSAs: HiFi with the SNP filters' defaults, ONT with --snp_min_af 0, 0.2 (its default) and 0.3.
set -euo pipefail
E=$(cd "$(dirname "$0")" && pwd)
B=${B:-$HOME/protal-lr-build}  # holds src/ (a checkout), build/protal, mini_db/, mini_db_close/
P=$B/build/protal
MINI=$B/mini_db
G=$MINI/gtdb_r226/simulation
C=$B/src/examples/mini_db/community.tsv
W=$B/snp_noise
mkdir -p "$W"
cd "$W"
db=$W/db
if [ ! -s "$db/model_ONT.xml" ]; then
    "$P" --unpack_db --db "$MINI/protal_db/database.protal" --unpack_dir "$db" -t 8 > unpack.log 2>&1
    for m in model_PB.xml model_ONT.xml; do cp "$db/model.xml" "$db/$m"; done
fi
for platform in hifi ont; do
    for seed in 5 6; do
        [ -s ${platform}_$seed.fq ] || python3 "$E/simulate_long_reads.py" --platform $platform --genomes "$G/genomes.tsv" \
            --community "$C" --out ${platform}_$seed --bases 10e6 --seed $seed > /dev/null
    done
done
run() {  # platform read_type tag extra...
    local platform=$1 type=$2 tag=$3; shift 3
    rm -rf "out_$tag"
    "$P" --db "$db" -1 ${platform}_5.fq,${platform}_6.fq --prefix ${platform}_a,${platform}_b --read_type $type \
        -o "out_$tag" -t 8 --no_qcmsa "$@" > "$tag.log" 2>&1
    echo "== $tag"
    python3 "$E/msa_noise.py" --strains "out_$tag/strains"
}
run hifi pb hifi
run ont ont ont_af0 --snp_min_af 0
run ont ont ont_af0.2
run ont ont ont_af0.3 --snp_min_af 0.3
