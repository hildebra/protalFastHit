#!/bin/bash
# Rebuilds the one- and four-thread s2.fq.gz of the failing run from the gtest log and compares their content.
set -u
SP=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/99a2c9fe-d722-41e1-bc32-3fc6baabdeb2/scratchpad
W=~/lrdet/fail39
mkdir -p $W
cd $W
perl $SP/unescape.pl $SP/fail278.txt 4 one_s2.fq.gz
perl $SP/unescape.pl $SP/fail278.txt 6 four_s2.fq.gz
for f in one_s2 four_s2; do
    gzip -t $f.fq.gz && echo "$f: gzip OK"
    zcat $f.fq.gz > $f.fq
    echo "$f: $(wc -c < $f.fq) bytes, $(($(wc -l < $f.fq) / 4)) reads, md5 $(md5sum < $f.fq | cut -c1-12)"
done
cmp one_s2.fq four_s2.fq && echo "same FASTQ content"
cmp one_s2.fq.gz four_s2.fq.gz | head -2
# the BGZF blocks of each file: offset, compressed size, content size
for f in one_s2 four_s2; do
    perl -e '
        local $/; open(my $h, "<:raw", $ARGV[0]) or die; my $d = <$h>; my $o = 0; my $k = 0; my $c = 0;
        while ($o < length $d) {
            my $bs = unpack("v", substr($d, $o + 16, 2)) + 1;
            my $isize = unpack("V", substr($d, $o + $bs - 4, 4));
            printf "%s block %d at %d: %d bytes, content %d at %d\n", $ARGV[0], $k++, $o, $bs, $isize, $c;
            $o += $bs; $c += $isize;
        }' $f.fq.gz > $f.blocks
    wc -l < $f.blocks
done
diff <(cut -d" " -f3- one_s2.blocks) <(cut -d" " -f3- four_s2.blocks)
