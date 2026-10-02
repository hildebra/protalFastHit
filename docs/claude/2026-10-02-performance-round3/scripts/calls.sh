#!/usr/bin/env bash
# calls.sh CGFILE: calls to wavefront_align (from callgrind's call records) and WFA2's self cost by function.
f=$1
awk '/^fn=/ { fn = $0 } /^cfn=/ { cfn = $0 } /^calls=/ { split($1, a, "="); if (cfn ~ /wavefront_align$/ || cfn ~ /\) wavefront_align$/) c += a[2] } END { print "calls to wavefront_align:", c }' $f
callgrind_annotate --threshold=100 $f 2>/dev/null | grep -E '^ *[0-9,]+ .*(wavefront|wfa|slab|cigar|mm_allocator|memset|memcpy)' | head -24 | sed -E 's/\[\/[^]]*\]//' | awk '{ $1 = $1; print }' | cut -c1-110
