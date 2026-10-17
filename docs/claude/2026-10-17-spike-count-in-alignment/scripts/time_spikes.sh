#!/usr/bin/env bash
# Before/after timing of the spike check: protal at HEAD (ref) and the working tree (work) on two local real mock
# samples with the reduced r226 database, 4 threads on cores 0-3, rounds alternating the binaries' order. Then
# --profile_only of work's last SAM (the SAM pass). Writes ~/spikecount/timing.tsv and compares the outputs.
set -uo pipefail
ROOT=$HOME/spikecount
DB=$HOME/Mocks/reduced_db_r226/database.protal
ROUNDS=${ROUNDS:-3}
declare -A R1 R2
R1[dohlman]=$HOME/Mocks/reads/d6331_dna/d6331_dna_SRR38158090/SRR38158090_1.fq.gz
R2[dohlman]=$HOME/Mocks/reads/d6331_dna/d6331_dna_SRR38158090/SRR38158090_2.fq.gz
R1[hku]=$HOME/Mocks/reads/d6331_pe_labs/d6331_hku_SRR23926922/SRR23926922_1.fq.gz
R2[hku]=$HOME/Mocks/reads/d6331_pe_labs/d6331_hku_SRR23926922/SRR23926922_2.fq.gz
runs=$ROOT/runs; rm -rf "$runs"; mkdir -p "$runs"
out=$ROOT/timing.tsv
printf 'round\tbin\tsample\twall_s\tuser_s\tmax_rss_kb\taligning_s\tprofiling_s\tspike_s\tprofile_sample_s\n' > "$out"
took() { grep -m1 -oE "^$1 took .*" "$2" | sed "s/^$1 took //" | awk '{ s = 0; for (i = 1; i <= NF; i++) { n = $i + 0;
    if ($i ~ /ms$/) s += n / 1000; else if ($i ~ /m$/) s += n * 60; else if ($i ~ /s$/) s += n } printf "%.3f", s }'; }
stage() { awk -F'\t' -v s="$1" '$1 == s { printf "%.4f", $2 }' "$2"; }
for r in $(seq 1 "$ROUNDS"); do
    if [ $((r % 2)) = 1 ]; then order="ref work"; else order="work ref"; fi
    for bin in $order; do
        for s in dohlman hku; do
            o=$runs/$bin.$s.$r
            taskset -c 0-3 /usr/bin/time -f '%e %U %M' -o "$o.time" "$ROOT/$bin/build/protal" --db "$DB" -1 "${R1[$s]}" -2 "${R2[$s]}" \
                --prefix "$s" -o "$o" -t 4 --no_qcmsa --no_strains > "$o.log" 2>&1 || { echo "$bin $s $r failed"; tail -5 "$o.log"; }
            read -r wall user rss < "$o.time"
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$r" "$bin" "$s" "$wall" "$user" "$rss" "$(took 'Aligning reads' "$o.log")" \
                "$(took Profiling "$o.log")" "$(stage 'Profiling: spike genes' "$o/misc/${s}_runtime.tsv")" \
                "$(stage 'Profiling: sample' "$o/misc/${s}_runtime.tsv")" >> "$out"
            [ "$r" = "$ROUNDS" ] || rm -f "$o"/*.sam*
        done
    done
done
column -t "$out"
echo
for s in dohlman hku; do
    a=$runs/ref.$s.$ROUNDS; b=$runs/work.$s.$ROUNDS
    echo "== $s: spike line ref / work"
    grep -E "^Sample $s: (no )?spike mode" "$a.log" | sed 's/; in .*//'
    grep -E "^Sample $s: (no )?spike mode" "$b.log" | sed 's/; in .*//'
    for f in "$s.profile" "$s.profile.spike_genes.tsv" "$s.profile.composition"; do
        if [ -e "$a/$f" ] || [ -e "$b/$f" ]; then cmp -s "$a/$f" "$b/$f" && echo "$f identical" || echo "$f DIFFERS"; fi
    done
    sam=$(ls "$b/$s".sam* | head -1)
    taskset -c 0-3 "$ROOT/work/build/protal" --db "$DB" --profile_only "$sam" --prefix "$s" -o "$runs/only.$s" -t 4 --no_qcmsa \
        --no_strains > "$runs/only.$s.log" 2>&1 || { echo "profile_only $s failed"; tail -5 "$runs/only.$s.log"; }
    echo "profile_only (the SAM pass): spike genes $(stage 'Profiling: spike genes' "$runs/only.$s/misc/${s}_runtime.tsv") s"
    grep -E "^Sample $s: (no )?spike mode" "$runs/only.$s.log" | sed 's/; in .*//'
    for f in "$s.profile" "$s.profile.spike_genes.tsv"; do
        if [ -e "$b/$f" ] || [ -e "$runs/only.$s/$f" ]; then cmp -s "$b/$f" "$runs/only.$s/$f" && echo "profile_only $f identical" || echo "profile_only $f DIFFERS"; fi
    done
    zcmp() { if [ "${1%.zst}" != "$1" ]; then zstd -dcq "$1"; else cat "$1"; fi; }
    asam=$(ls "$a/$s".sam* | head -1)
    cmp -s <(zcmp "$asam" | grep -v '^@') <(zcmp "$sam" | grep -v '^@') && echo "SAM records identical" || echo "SAM records DIFFER"
done
