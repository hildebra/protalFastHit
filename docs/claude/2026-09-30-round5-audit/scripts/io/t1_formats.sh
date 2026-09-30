#!/bin/bash
# T1: determinism across thread counts and SAM formats; rerun reuse; --profile_only per format.
set -u
source $(dirname "$0")/lib.sh
T=$W/t1; rm -rf $T; mkdir -p $T; cd $T
for cfg in "sam 1" "sam 2" "gz 1" "gz 2" "zst 1" "zst 2"; do
  set -- $cfg
  d=${1}_t$2
  echo "--- $d"; prun $d.log --db $DB -1 $R/sa_R1.fq -2 $R/sa_R2.fq --prefix sa -o $d -t $2 --no_qcmsa --sam_format $1
  ls $d/
done
for d in sam_t1 sam_t2 gz_t1 gz_t2 zst_t1 zst_t2; do
  f=$(ls $d/sa.sam* | head -1)
  records $f > $d.rec
  samtext $f | grep '^@' > $d.hdr
  samtext $f | grep -v '^@' > $d.rec_unsorted
  profiles $d > $d.prof
  echo "$d $(wc -l < $d.rec) records, $(wc -l < $d.hdr) header lines, md5 rec $(md5sum < $d.rec | cut -c1-8) hdr $(md5sum < $d.hdr | cut -c1-8) unsorted $(md5sum < $d.rec_unsorted | cut -c1-8) prof $(md5sum < $d.prof | cut -c1-8)"
done
# strain outputs across configs
for d in sam_t1 sam_t2 gz_t1 gz_t2 zst_t1 zst_t2; do
  echo "$d strains: $(find $d/strains -type f 2>/dev/null | sort | xargs cat 2>/dev/null | md5sum | cut -c1-8) files $(find $d/strains -type f 2>/dev/null | wc -l)"
done
echo "--- validity checks"
gzip -t gz_t2/sa.sam.gz && echo "gzip -t ok"
zstd -tq zst_t2/sa.sam.zst && echo "zstd -t ok"
ls -la */
echo "--- header sample"; head -5 zst_t2.hdr; tail -2 zst_t2.hdr
echo "--- @SQ genes vs genes named in records"
for d in sam_t1 zst_t2; do
  cut -f3 $d.rec | grep -v '^\*$' | sort -u > $d.rname
  cut -f7 $d.rec | grep -v -e '^\*$' -e '^=$' | sort -u > $d.rnext
  sort -u $d.rname $d.rnext > $d.named
  grep '^@SQ' $d.hdr | sed 's/.*SN:\([^\t]*\).*/\1/' | sort -u > $d.sq
  echo "$d named $(wc -l < $d.named) sq $(wc -l < $d.sq) diff: $(comm -3 $d.named $d.sq | wc -l)"
done
echo "--- rerun (skip) each format"
for d in sam_t1 gz_t2 zst_t2; do
  fmt=${d%_t*}
  m1=$(stat -c %Y.%s $d/sa.sam*)
  cp $d.prof $d.prof.first
  prun $d.rerun.log --db $DB -1 $R/sa_R1.fq -2 $R/sa_R2.fq --prefix sa -o $d -t 2 --no_qcmsa --sam_format $fmt
  grep -i -E 'skip|error|warn' $d.rerun.log | head -5
  m2=$(stat -c %Y.%s $d/sa.sam*)
  profiles $d > $d.prof2
  echo "$d rerun: sam unchanged=$([ "$m1" = "$m2" ] && echo yes || echo no) profile same=$(cmp -s $d.prof.first $d.prof2 && echo yes || echo no)"
done
echo "--- profile_only each format"
for d in sam_t1 gz_t2 zst_t2; do
  f=$(ls $d/sa.sam* | head -1)
  prun po_$d.log --db $DB --profile_only $f --prefix sa -o po_$d -t 2 --no_qcmsa
  profiles po_$d > po_$d.prof
  # compare the profile files by content (names identical: sa)
  echo "po_$d prof md5 $(md5sum < po_$d.prof | cut -c1-8) vs run $(md5sum < $d.prof | cut -c1-8)"
  grep -i -E 'error|warn|fail' po_$d.log | head -5
done
diff <(sed 's/^== .*//' sam_t1.prof) <(sed 's/^== .*//' po_zst_t2.prof) | head -20
