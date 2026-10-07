#!/bin/bash
# The FP and FN records of a build's error reads (model_logs/error_reads/<type>/<set>/<point>/<sample>.sam.zst, as
# error_reads.py wrote them), without the reads of "unseen" species and without QUAL, into one archive to copy back.
#
#   bash extract_error_records.sh [ERROR_READS_DIR] [OUT.tar.gz]
#   TYPES="pe se" JOBS=8 bash extract_error_records.sh ...
#
# A record is kept if its xe tag names an FP or FN taxon (FP:<taxid>, FN:<taxid>), or seeded:/source: of a taxid that
# the sample's taxa.tsv lists as FN; seeded:/source: of unseen species are left out. The header keeps only the @SQ
# lines of the genes the kept records name. The taxa tables keep their FP and FN rows; summary.tsv is copied.
set -uo pipefail

SRC=${1:-/hpc-home/hildebra/DB/protal/protal0.7.8_r226_v15/model_logs/error_reads}
OUT=${2:-$PWD/v15_error_records.tar.gz}
TYPES=${TYPES:-pe se pb ont}
JOBS=${JOBS:-4}
WORK=$(mktemp -d "${TMPDIR:-/tmp}/error_records.XXXXXX")
trap 'rm -rf "$WORK"' EXIT
DEST=$WORK/v15_error_records
export SRC DEST

one() {
    local sam=$1 rel base taxa out
    rel=${sam#"$SRC"/}
    base=${sam%.sam.zst}
    base=${base%.FP}
    base=${base%.FN}
    taxa=$base.taxa.tsv
    out=$DEST/${rel%.sam.zst}.sam.zst
    mkdir -p "$(dirname "$out")"
    if ! zstd -dcq "$sam" | awk -F'\t' -v OFS='\t' -v taxa="$taxa" '
        BEGIN { while ((getline line < taxa) > 0) { split(line, f, "\t"); if (f[3] == "FN") fn[f[4]] = 1 } }
        /^@/ { head[++nh] = $0; next }
        {
            xe = ""
            for (i = 12; i <= NF; i++) if (substr($i, 1, 5) == "xe:Z:") { xe = substr($i, 6); break }
            n = split(xe, why, ","); keep = 0
            for (j = 1; j <= n; j++) {
                split(why[j], p, ":")
                if (p[1] == "FP" || p[1] == "FN" || ((p[1] == "seeded" || p[1] == "source") && (p[2] in fn))) keep = 1
            }
            if (!keep) next
            if (NF >= 11) $11 = "*"
            rec[++nr] = $0
            if ($3 != "*") gene[$3] = 1
            if ($7 != "*" && $7 != "=") gene[$7] = 1
        }
        END {
            for (i = 1; i <= nh; i++) {
                if (substr(head[i], 1, 4) == "@SQ\t") { split(head[i], s, "\t"); if (!(substr(s[2], 4) in gene)) continue }
                print head[i]
            }
            print "@CO\textract_error_records.sh: FP and FN records only (seeded/source of FN taxa kept, of unseen species not), QUAL *"
            for (i = 1; i <= nr; i++) print rec[i]
        }' | zstd -q -19 -o "$out"; then
        echo "incomplete (kept what decompressed): $rel" >&2
    fi
    if [ -f "$taxa" ]; then
        awk -F'\t' 'NR == 1 || $3 == "FP" || $3 == "FN"' "$taxa" > "$DEST/${taxa#"$SRC"/}"
    fi
}
export -f one

for t in $TYPES; do
    if [ ! -d "$SRC/$t" ]; then
        echo "$t: no $SRC/$t"
        continue
    fi
    n=$(find "$SRC/$t" -name '*.sam.zst' ! -name '*.unseen.sam.zst' | wc -l)
    echo "$t: $n SAMs, $(du -sh "$SRC/$t" | cut -f1)"
    find "$SRC/$t" -name '*.sam.zst' ! -name '*.unseen.sam.zst' -print0 | xargs -0 -r -P "$JOBS" -I{} bash -c 'one "$1"' _ {}
    mkdir -p "$DEST/$t"
    [ -f "$SRC/$t/summary.tsv" ] && cp "$SRC/$t/summary.tsv" "$DEST/$t/"
    echo "$t: $(find "$DEST/$t" -name '*.sam.zst' | wc -l) SAMs, $(du -sh "$DEST/$t" | cut -f1) after filtering"
done

tar -C "$WORK" -czf "$OUT" v15_error_records
ls -lh "$OUT"
