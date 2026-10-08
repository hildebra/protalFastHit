#!/bin/bash
# The composition check (docs/claude/2026-10-08-sample-composition): a synthetic world of 40 species whose background DNA
# is drawn per species from 0.2-2 Mb, two databases of it (every species; every 4th species held out, as species a
# database lacks), six samples of 200,000 read pairs, protal on both; evaluate.py compares each sample's
# <profile>.composition with the simulator's truth. Steps whose outputs exist are skipped.
#
#   SRC=~/protal-composition/src BIN=~/protal-composition/build bash run_validation.sh OUTDIR
set -euo pipefail
SRC=${SRC:-$HOME/protal-composition/src}
BIN=${BIN:-$HOME/protal-composition/build}
OUT=${1:-$HOME/protal-composition/val}
T=${THREADS:-4}
SPECIES=${SPECIES:-40}
SAMPLES=${SAMPLES:-6}
PAIRS=${PAIRS:-200000}
PER_SAMPLE=${PER_SAMPLE:-12}
LENGTHS=${LENGTHS:-200000-2000000}
mkdir -p "$OUT"
cd "$OUT"

if [ ! -d gtdb ]; then
    python3 "$SRC/scripts/mini_db/gtdb_like_lineages.py" --species "$SPECIES" --archaea 0.1 --seed 5 > lineages.txt
    python3 "$SRC/scripts/mini_db/simulate_gtdb_release.py" --outdir gtdb --lineages lineages.txt --genomes_per_species 2 \
        --genome_length "$LENGTHS" --seed 5 > release.log 2>&1
fi
if [ ! -s db_all/database.protal ]; then
    python3 "$SRC/scripts/mini_db/gtdb_to_protal_db.py" --gtdb gtdb --outdir db_all -t "$T" > convert.log 2>&1
    awk -F';' 'NR % 4 == 0 { print $7 }' lineages.txt > held_out.txt
    python3 "$SRC/scripts/mini_db/gtdb_to_protal_db.py" --from_db db_all --exclude_species held_out.txt --outdir db_held \
        > derive.log 2>&1
    for db in db_all db_held; do
        "$BIN/protal" --build --no_profile -t "$T" --db "$db" --reference "$db/reference.fna" \
            --full_reference "$db/full_reference.fna" > "$db.build.log" 2>&1
    done
fi
if [ ! -s sims/protal.meta ]; then
    "$BIN/simulate_metagenomes" --genome_table gtdb/simulation/genomes.tsv --output_dir sims --samples "$SAMPLES" \
        --sample_prefix s --total_read_pairs "$PAIRS" --species_per_sample "$PER_SAMPLE" --distribution lognormal \
        --strains_per_species 0.3 --seed 7 -t "$T" --protal_metafile sims/protal > simulate.log 2>&1
fi

# Absolute output folders: with --map, a relative -o came out doubled (run_db_all/run_db_all/...) on 2026-10-08.
for db in db_all db_held; do
    rm -rf "$OUT/run_$db"
    "$BIN/protal" --db "$db" --map sims/protal.meta -t "$T" -o "$OUT/run_$db" --no_strains --no_qcmsa > "run_$db.log" 2>&1
done
# The profiles again from the SAMs alone (--profile_only reads the scanned reads from the SAM header).
rm -rf "$OUT/rerun_db_held"
"$BIN/protal" --db db_held --profile_only "$OUT/run_db_held/alignments/*.sam.zst" -t "$T" -o "$OUT/rerun_db_held" --no_strains \
    --no_qcmsa > rerun_db_held.log 2>&1
for f in "$OUT"/run_db_held/profiles/*.profile; do
    cmp -s "$f" "$OUT/rerun_db_held/$(basename "$f")" || echo "rerun differs: $(basename "$f")" >&2
done

python3 "$SRC/docs/claude/2026-10-08-sample-composition/scripts/evaluate.py" --manifest sims/manifest.tsv \
    --held_out held_out.txt --run all=run_db_all --run held=run_db_held | tee evaluation.tsv
