#!/usr/bin/env bash
# ONT-like reads (simulate_long_reads.py --platform ont) of the mock community, profiled natively
# (--read_type ont) and, cut into 1 kb pieces, as single-end reads. The database's model.xml stands
# in for model_ONT.xml and model_se.xml. EXTRA: more protal options for the native run; TAG names
# the run. Usage: run_ont.sh [bases] [seed]
set -euo pipefail
B=${B:-$HOME/protal-lr-build}  # holds src/ (a checkout), build/protal, mini_db/
P=$B/build/protal
E=$(cd "$(dirname "$0")" && pwd)
BASES=${1:-20e6}
SEED=${2:-5}
MINI=${MINI:-$B/mini_db}
TAG=${TAG:-default}
ERR=${ERR:-}
W=$B/ont_eval_$(basename "$MINI")_$SEED${ERR:+_err$ERR}
G=$MINI/gtdb_r226/simulation
COMMUNITY=$B/src/examples/mini_db/community.tsv
mkdir -p "$W"
cd "$W"

db=$W/db
if [ ! -s "$db/model_ONT.xml" ]; then
    "$P" --unpack_db --db "$MINI/protal_db/database.protal" --unpack_dir "$db" -t 8 > unpack.log 2>&1
    cp "$db/model.xml" "$db/model_se.xml"
    cp "$db/model.xml" "$db/model_ONT.xml"
fi
[ -s ont.fq ] || python3 "$E/simulate_long_reads.py" --platform ont --genomes "$G/genomes.tsv" --community "$COMMUNITY" --out ont --bases "$BASES" --seed "$SEED" ${ERR:+--error $ERR}

evaluate() {
    python3 "$E/evaluate.py" --sam "out/$1.sam" --truth ont.truth.tsv --taxonomy "$db/internal_taxonomy.dmp" \
        --genomes "$G/genomes.tsv" --markers "$G/marker_positions.tsv" --profile_log "out/$1.profile.log" \
        --community "$COMMUNITY" --label "$1" --log "$1.log" ${2:-}
}

rm -f "out/native_$TAG".*
"$P" --db "$db" -1 ont.fq --read_type ont --prefix "native_$TAG" -o out -t 8 --no_qcmsa --no_strains ${EXTRA:-} > "native_$TAG.log" 2>&1
evaluate "native_$TAG" --header
grep -h 'settled by their read\|seeded in chunks' "native_$TAG.log" || true
if [ -n "${SPLIT:-}" ]; then
    [ -s split1000.fq ] || python3 "$E/split_reads.py" --fastq ont.fq --out split1000.fq --length 1000 > /dev/null
    rm -f out/split1000.*
    "$P" --db "$db" -1 split1000.fq --prefix split1000 -o out -t 8 --no_qcmsa --no_strains > split1000.log 2>&1
    evaluate split1000
fi
