#!/bin/bash
# T7: gene loading: the arena from each kind of reference (raw, bundle, zstd CLI single frame,
# seekable with tiny frames so genes span frames), --preload_genomes_off, --full_sam_header; 1 vs 2 threads.
set -u
source $(dirname "$0")/lib.sh
T=$W/t7; rm -rf $T; mkdir -p $T; cd $T
BASE=$W/t1/sam_t1.rec; BASEPROF=$W/t1/sam_t1/sa.profile
mkdb() { # mkdb DIR: a copy of the unpacked db without reference.fna
  mkdir -p $1; for f in $DB/*; do [ "$(basename $f)" = reference.fna ] || ln -sf $f $1/; done
}
mkdb db_cli; zstd -q -c $DB/reference.fna > db_cli/reference.fna.zst
mkdb db_seek; python3 $S/mkreads.py seekable $DB/reference.fna db_seek/reference.fna.zst 997
zstd -lv db_seek/reference.fna.zst 2>&1 | grep -E 'Frames|Decompressed' | head -2
check() { # check NAME DB THREADS [extra]
  local name=$1 db=$2 t=$3; shift 3
  prun $name.log --db $db -1 $R/sa_R1.fq -2 $R/sa_R2.fq --prefix sa -o $name -t $t --no_qcmsa --no_strains --sam_format zst "$@" > $name.rc
  local s=$name/sa.sam.zst
  if [ -e $s ]; then
    records $s > $name.rec
    printf '%-14s %s records=%s vs_base=%s profile=%s | %s\n' $name "$(cat $name.rc)" $(wc -l < $name.rec) \
      "$(cmp -s $name.rec $BASE && echo SAME || echo diff)" "$(cmp -s $name/sa.profile $BASEPROF && echo SAME || echo diff)" \
      "$(grep -a -E 'Preload|\[ERROR\]|Cannot|Invalid' $name.log | head -2 | tr '\n' '|')"
  else
    printf '%-14s %s no SAM | %s\n' $name "$(cat $name.rc)" "$(grep -a -E '\[ERROR\]|Cannot|Invalid|rror' $name.log | head -3 | tr '\n' '|')"
  fi
}
check raw1 $DB 1; check raw2 $DB 2
check bundle1 $W/bundle.protal 1; check bundle2 $W/bundle.protal 2
check cli1 $T/db_cli 1; check cli2 $T/db_cli 2
check seek1 $T/db_seek 1; check seek2 $T/db_seek 2
check nopreload $DB 2 --preload_genomes_off
check fullhdr $DB 2 --full_sam_header
samtext fullhdr/sa.sam.zst | grep -c '^@SQ'; echo "genes in map: $(wc -l < $DB/reference.map)"
prun po_full.log --db $DB --profile_only fullhdr/sa.sam.zst --prefix fh -o po_full -t 2 --no_qcmsa --no_strains
echo "profile_only of the full-header SAM: $(cmp -s po_full/fh.profile $BASEPROF && echo SAME || echo diff)"
