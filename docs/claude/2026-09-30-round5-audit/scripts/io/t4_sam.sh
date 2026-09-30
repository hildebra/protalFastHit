#!/bin/bash
# T4: reading SAMs back: --profile_only on damaged, empty, zero-gene, CRLF and foreign SAMs in each format.
set -u
source $(dirname "$0")/lib.sh
T=$W/t4; rm -rf $T; mkdir -p $T/in; cd $T
I=$T/in
T1=$W/t1
BASEPROF=$T1/sam_t1/sa.profile
# random reads: nothing aligns
python3 - <<'EOF'
import random
r = random.Random(11)
with open("/home/falk/audit6/io/t4/in/rnd_R1.fq", "w") as a, open("/home/falk/audit6/io/t4/in/rnd_R2.fq", "w") as b:
    for i in range(200):
        for fh, m in ((a, 1), (b, 2)):
            s = "".join(r.choice("ACGT") for _ in range(100))
            fh.write(f"@r{i}/{m}\n{s}\n+\n{'I'*100}\n")
EOF
for f in sam gz zst; do
  prun rnd_$f.log --db $DB -1 $I/rnd_R1.fq -2 $I/rnd_R2.fq --prefix rnd -o rnd_$f -t 2 --no_qcmsa --sam_format $f > rnd_$f.rc
  s=$(ls rnd_$f/rnd.sam* | grep -v err)
  echo "zero-gene align $f: $(cat rnd_$f.rc) file=$s size=$(stat -c %s $s) header: $(samtext $s | tr '\n' '|' | cut -c1-120)"
  grep -a -i -E 'error|warn|no usable' rnd_$f.log | head -3
done
po() { # po NAME SAMFILE : --profile_only, report rc, errors, whether the profile equals the plain run's
  local name=$1 sam=$2
  prun po_$name.log --db $DB --profile_only $sam --prefix $name -o po -t 2 --no_qcmsa --no_strains > po_$name.rc
  local p=po/$name.profile same=-
  [ -e $p ] && same=$(cmp -s $p $BASEPROF && echo SAME || echo "diff($(wc -l < $p) lines)")
  printf '%-16s %s profile=%s | %s\n' $name "$(cat po_$name.rc)" "$same" \
    "$(grep -a -i -E 'error|warn|truncat|corrupt|missing|no usable|empty' po_$name.log | grep -a -v -E '^(snp|msa)' | head -3 | tr '\n' '|' | cut -c1-330)"
}
for f in sam gz zst; do po zero_$f $(ls rnd_$f/rnd.sam* | grep -v err); done
SAM=$T1/sam_t1/sa.sam; GZ=$T1/gz_t2/sa.sam.gz; ZST=$T1/zst_t2/sa.sam.zst
po ok_sam $SAM; po ok_gz $GZ; po ok_zst $ZST
# damaged copies
zs=$(stat -c %s $ZST); gs=$(stat -c %s $GZ); ss=$(stat -c %s $SAM)
head -c $((zs/2)) $ZST > $I/zmid.sam.zst
# frame boundaries of the zst file (seek table): cut after the last data frame = drop the seek table
python3 - "$ZST" "$I" <<'EOF'
import struct, sys
p, out = sys.argv[1], sys.argv[2]
d = open(p, "rb").read()
# seek table footer: last 9 bytes = frames(4) desc(1) magic(4)
nframes = struct.unpack("<I", d[-9:-5])[0]
table = 8 + nframes * 8 + 9
body = d[:-table]
open(out + "/znotable.sam.zst", "wb").write(body)
# drop the last data frame too (cut at an earlier frame boundary), keep no table
sizes = [struct.unpack("<I", d[len(d)-table+8+8*i: len(d)-table+12+8*i])[0] for i in range(nframes)]
open(out + "/zcutframe.sam.zst", "wb").write(body[:len(body) - sizes[-1]])
# drop the last data frame but keep a (now wrong) table
open(out + "/zdropframe.sam.zst", "wb").write(body[:len(body) - sizes[-1]] + d[-table:])
print("frames", nframes, "sizes tail", sizes[-3:])
EOF
cp $ZST $I/zflip.sam.zst; printf '\x55' | dd of=$I/zflip.sam.zst bs=1 seek=$((zs/2)) conv=notrunc 2>/dev/null
head -c $((gs/2)) $GZ > $I/gmid.sam.gz
head -c $((gs-28)) $GZ > $I/gnoeof.sam.gz
cp $GZ $I/gflip.sam.gz; printf '\x55' | dd of=$I/gflip.sam.gz bs=1 seek=$((gs/2)) conv=notrunc 2>/dev/null
head -c $((ss/2)) $SAM > $I/smid.sam                     # cut mid-line
head -c $((ss/2)) $SAM | sed '$d' > $I/sline.sam          # cut at a line boundary (undetectable)
: > $I/empty.sam; : > $I/empty.sam.gz; : > $I/empty.sam.zst
cp $SAM $I/plain_named.sam.zst; cp $SAM $I/plain_named.sam.gz
gzip -c $SAM > $I/plaingzip.sam.gz                        # single-member gzip (older protal, pigz)
zstd -qc $SAM > $I/plainzstd.sam.zst                      # zstd CLI, no marker
head -c $(( $(stat -c %s $I/plainzstd.sam.zst) - 100 )) $I/plainzstd.sam.zst > $I/plainzstd_cut.sam.zst
sed 's/$/\r/' $SAM > $I/crlf.sam
grep -v '^@SQ' $SAM > $I/nosq.sam                          # records name genes the header lacks
awk 'BEGIN{OFS="\t"} /^@SQ/ && !d {sub(/LN:[0-9]+/, "LN:1"); d=1} {print}' $SAM > $I/badlen.sam
awk '/^@SQ/ && !d {print "@SQ\tSN:999_1\tLN:100"; d=1} {print}' $SAM > $I/foreign.sam
head -c -1 $SAM > $I/nonl.sam
for v in zmid.sam.zst znotable.sam.zst zcutframe.sam.zst zdropframe.sam.zst zflip.sam.zst gmid.sam.gz gnoeof.sam.gz gflip.sam.gz smid.sam sline.sam empty.sam empty.sam.gz empty.sam.zst plain_named.sam.zst plain_named.sam.gz plaingzip.sam.gz plainzstd.sam.zst plainzstd_cut.sam.zst crlf.sam nosq.sam badlen.sam foreign.sam nonl.sam; do
  po ${v//./_} $I/$v
done
