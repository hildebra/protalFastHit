#!/bin/bash
# Local runs of the sparse-first gate with each pick (--sparse_first_pick hash|stride; docs/claude/2026-10-13-stride-pick):
# the archaea world's 100k pairs, and their first mates typed as PacBio reads, ungated and at C 2, 8 and 16. Per run:
# the pairs (reads) given up, the sampled lookups, the records written and the aligning time. 4 cores, niced.
set -uo pipefail
P=${P:-$HOME/stride-build/build/protal}
R=${R:-$HOME/stride-run/gate}
mkdir -p "$R"
printf 'sample\tgate\tpick\tgiven_up\tof\tsampled_lookups\trecords\taligning\n' > "$R/gate.tsv"
run() {  # sample gate pick reads...
    local sample=$1 gate=$2 pick=$3; shift 3
    local name="${sample}_c${gate}_${pick}" out="$R/${sample}_c${gate}_${pick}"
    rm -rf "$out"
    local extra=()
    [ "$gate" != 0 ] && extra=(--sparse_first "$gate" --sparse_first_pick "$pick")
    taskset -c 0-3 nice -n 5 "$P" --db ~/seed/db "$@" --prefix "$name" -o "$out" -t 4 --no_strains --no_qcmsa "${extra[@]}" > "$out.log" 2>&1 \
        || { echo "$name failed, see $out.log"; return; }
    local line given of sampled records aligning
    line=$(grep -m1 "^Sample $name sparse first (" "$out.log")
    given=$(echo "$line" | grep -oE 'hit\(s\)\): [0-9]+' | grep -oE '[0-9]+$')
    of=$(echo "$line" | grep -oE ' of [0-9]+ (pairs|reads)' | grep -oE '[0-9]+')
    sampled=$(echo "$line" | grep -oE '[0-9]+ sampled lookups' | grep -oE '^[0-9]+')
    records=$(grep -m1 "^Sample $name: .* candidate alignments tried" "$out.log" | grep -oE '[0-9]+ records written' | grep -oE '^[0-9]+')
    aligning=$(grep -m1 '^Aligning reads took' "$out.log" | sed -E 's/.*took //')
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$sample" "$gate" "$pick" "${given:-NA}" "${of:-NA}" "${sampled:-NA}" "${records:-NA}" "${aligning:-NA}" >> "$R/gate.tsv"
}
PE=(-1 ~/seed/samples/pe100k_R1.fq.gz -2 ~/seed/samples/pe100k_R2.fq.gz)
PB=(-1 ~/seed/samples/pe100k_R1.fq.gz --read_type pb)
run arc 0 none "${PE[@]}"
run arcpb 0 none "${PB[@]}"
for gate in 2 8 16; do
    for pick in hash stride; do
        run arc "$gate" "$pick" "${PE[@]}"
        run arcpb "$gate" "$pick" "${PB[@]}"
    done
done
cat "$R/gate.tsv"
