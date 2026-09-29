#!/usr/bin/env bash
# Native long reads vs reads split into short pieces: simulated HiFi reads from the mini database's
# mock community (examples/mini_db/community.tsv), profiled with protal --read_type pacbio and, cut
# into 150, 250 and 1000 bp pieces, as single-end reads. The database's model.xml stands in for
# model_pacbio.xml and model_se.xml. Usage: run.sh [bases] [seed]
set -euo pipefail
B=${B:-$HOME/protal-lr-build}  # holds src/ (a checkout), build/protal, mini_db/
P=$B/build/protal
E=$(cd "$(dirname "$0")" && pwd)
BASES=${1:-20e6}
SEED=${2:-5}
MINI=${MINI:-$B/mini_db}
W=$B/lr_eval_$(basename "$MINI")_$SEED
G=$MINI/gtdb_r226/simulation
COMMUNITY=$B/src/examples/mini_db/community.tsv
rm -rf "$W"
mkdir -p "$W"
cd "$W"

db=$W/db
"$P" --unpack_db --db "$MINI/protal_db/database.protal" --unpack_dir "$db" -t 8 > unpack.log 2>&1
cp "$db/model.xml" "$db/model_se.xml"
cp "$db/model.xml" "$db/model_pacbio.xml"

python3 "$E/simulate_hifi.py" --genomes "$G/genomes.tsv" --community "$COMMUNITY" --out hifi --bases "$BASES" --seed "$SEED"

evaluate() {
    python3 "$E/evaluate.py" --sam "out/$1.sam" --truth hifi.truth.tsv --taxonomy "$db/internal_taxonomy.dmp" \
        --genomes "$G/genomes.tsv" --markers "$G/marker_positions.tsv" --profile_log "out/$1.profile.log" \
        --community "$COMMUNITY" --label "$1" --log "$1.log" ${2:-}
}

"$P" --db "$db" -1 hifi.fq --read_type pacbio --prefix native -o out -t 8 --no_qcmsa --no_strains > native.log 2>&1
evaluate native --header > results.tsv
for length in 150 250 1000; do
    python3 "$E/split_reads.py" --fastq hifi.fq --out "split$length.fq" --length "$length" > /dev/null
    "$P" --db "$db" -1 "split$length.fq" --prefix "split$length" -o out -t 8 --no_qcmsa --no_strains > "split$length.log" 2>&1
    evaluate "split$length" >> results.tsv
done
column -t -s $'\t' results.tsv
ls -la out/*.sam | awk '{print $5, $9}'
