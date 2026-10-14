#!/bin/bash
# The end-to-end tests' kept folders of the two binaries (TMPDIR tmp_base, tmp_split), paired by test class (mkdtemp's
# 8-character suffix cut off): every file compared byte for byte (a .zst or .gz decompressed first); a SAM's records
# also as a set (the alignment threads write them in the order they finish); otherwise the differing ones listed, with
# whether they still differ once folder paths and times are masked.
C=${1:-~/split-cmp}
A=${2:-tmp_base}
B=${3:-tmp_split}
cd $C
text() { if [[ $1 == *.zst ]]; then zstd -dcq "$1"; elif [[ $1 == *.gz ]]; then zcat "$1"; else cat "$1"; fi 2>/dev/null; }
mask() {
    text "$1" | sed -E "s#$C/(tmp|db)_(base|split)[A-Za-z0-9_/.]*#PATH#g; s#/tmp[A-Za-z0-9_./]*#PATH#g; s#$HOME/split-(base|build)/build/#BIN/#g" |
        sed -E 's/[0-9]+(\.[0-9]+)?( ?(s|ms|GB|MB|kB|%))\b/N\2/g; s/[0-9]+\.[0-9]+/F/g'
}
same=0; sam_order=0; masked_same=0; differ=0; missing=0; dirs=0
for bdir in $A/*/; do
    name=$(basename $bdir)
    prefix=${name:0:${#name}-8}
    sdir=$(ls -d $B/${prefix}*/ 2>/dev/null | head -1)
    [ -z "$sdir" ] && { echo "no split folder for $prefix"; missing=$((missing + 1)); continue; }
    dirs=$((dirs + 1))
    while IFS= read -r f; do
        if [ ! -e "$sdir/$f" ]; then echo "missing in split: $prefix/$f"; missing=$((missing + 1)); continue; fi
        if cmp -s "$bdir/$f" "$sdir/$f" || cmp -s <(text "$bdir/$f") <(text "$sdir/$f"); then same=$((same + 1)); continue; fi
        if [[ $f == *.sam* ]] && cmp -s <(text "$bdir/$f" | grep '^@') <(text "$sdir/$f" | grep '^@') &&
           cmp -s <(text "$bdir/$f" | grep -v '^@' | sort) <(text "$sdir/$f" | grep -v '^@' | sort); then
            sam_order=$((sam_order + 1)); continue
        fi
        if cmp -s <(mask "$bdir/$f") <(mask "$sdir/$f"); then masked_same=$((masked_same + 1)); continue; fi
        differ=$((differ + 1)); echo "DIFFERS: $prefix/$f"
    done < <(cd $bdir && find . -type f | sort)
done
echo "$dirs folder pairs: identical $same, SAMs with the same header and records in another order $sam_order, other files equal but for paths and times $masked_same, differing $differ, missing $missing"
