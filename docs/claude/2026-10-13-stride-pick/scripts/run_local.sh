#!/bin/bash
# Local --stage1_probe runs for the stride pick (docs/claude/2026-10-13-stride-pick): 4 cores, niced.
# P: a protal with the stride rules; the worlds of earlier reports in WSL's home (~/strainfn, ~/seed).
set -uo pipefail
P=${P:-$HOME/stride-build/build/protal}
R=${R:-$HOME/stride-run}
here=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$R"
# 1. Strain world (docs/claude/2026-10-11-strain-false-negatives): the 40 samples' reads, of test strains the
#    representatives' database lacks, as one paired-end sample.
if [ ! -s "$R/strains_1.fq.gz" ]; then
    cat $(ls ~/strainfn/reads/*_1.fq.gz | sort) > "$R/strains_1.fq.gz"
    cat $(ls ~/strainfn/reads/*_2.fq.gz | sort) > "$R/strains_2.fq.gz"
fi
echo "strain pairs: $(( $(zcat "$R/strains_1.fq.gz" | wc -l) / 4 ))"
run() {  # name db out extra...
    local name=$1 db=$2 out=$3; shift 3
    rm -rf "$out"
    /usr/bin/time -v taskset -c 0-3 nice -n 5 "$P" --db "$db" "$@" --prefix "$name" -o "$out" -t 4 --no_strains --no_qcmsa > "$out.log" 2>&1
    echo "$name rc=$? $(grep -m1 'Elapsed (wall' "$out.log" | awk '{print $NF}')"
    grep -m1 "^Sample $name stage-1 probe:" "$out.log"
}
STRAINS=(-1 "$R/strains_1.fq.gz" -2 "$R/strains_2.fq.gz" --model "$HOME/strainfn/model_s1/model_pe.xml")
run strains ~/strainfn/db_reps "$R/strains_probe" "${STRAINS[@]}" --stage1_probe
run strains ~/strainfn/db_reps "$R/strains_plain" "${STRAINS[@]}"
cmp "$R/strains_probe/strains.profile" "$R/strains_plain/strains.profile" && echo "strains: profiles identical"
# 2. Archaea world (docs/claude/2026-10-12-seeding-profile): 100k simulated pairs, most outside the marker genes.
run arc ~/seed/db "$R/arc_probe" -1 ~/seed/samples/pe100k_R1.fq.gz -2 ~/seed/samples/pe100k_R2.fq.gz --stage1_probe
# 3. The long-read path: the archaea world's first mates typed as PacBio reads (one chunk each).
run arcpb ~/seed/db "$R/arcpb_probe" -1 ~/seed/samples/pe100k_R1.fq.gz --read_type pb --stage1_probe
python3 "$here/stride_compare.py" --dir "$R/strains_probe/misc" --sample strains --profile "$R/strains_probe/strains.profile" > "$R/strains_compare.md"
python3 "$here/stride_compare.py" --dir "$R/arc_probe/misc" --sample arc --profile "$R/arc_probe/arc.profile" > "$R/arc_compare.md"
python3 "$here/stride_compare.py" --dir "$R/arcpb_probe/misc" --sample arcpb > "$R/arcpb_compare.md"
