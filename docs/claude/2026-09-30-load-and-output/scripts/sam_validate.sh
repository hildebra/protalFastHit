#!/bin/bash
# The SAM output of NEW against OLD (before the formats and the trimmed header), 1 and 8 threads:
#   sam_validate.sh OLD_BIN NEW_BIN DB READS_DIR NAME [READS_DIR NAME ...]
# READS_DIR holds NAME_R1.fq.gz and NAME_R2.fq.gz. Checks: --full_sam_header reproduces every output
# of OLD; the default changes only the SAM header, which lists exactly the genes the records name;
# .sam.gz and .sam.zst decompress to the plain SAM, pass gzip -t / zstd -t and profile alike, also
# at 8 threads; a rerun finds them and profiles them again.
source "$(dirname "$0")/env.sh"
old=$1; new=$2; db=$3; shift 3
O=$OUT/sam_validate; rm -rf $O; mkdir -p $O
same_outputs() {  # every output of two folders but runtime.tsv and the files matching $3
  local a=$1 b=$2 skip=${3:-^$} d=0
  for f in $(cd $a && find . -type f ! -name '*_runtime.tsv' | grep -vE "$skip" | sort); do
    cmp -s $a/$f $b/$f || { d=1; echo "    differs: $f"; }
  done
  [ $(cd $a && find . -type f | wc -l) = $(cd $b && find . -type f | wc -l) ] || { d=1; echo "    different file lists"; }
  [ $d = 0 ] && echo "    identical" || true
}
map() {  # map NAME DIR DS SAMNAME: a map for one sample, named as the -1/-2 runs name it
  printf '#OUTPUT_DIR\t%s\n#SAMPLEID\tPREFIX\tFIRST\tSECOND\tSAM\n%s\t%s\t%s\t%s\t%s\n' "$O/$1" "${3}_R" "${3}_R" \
    "$2/${3}_R1.fq.gz" "$2/${3}_R2.fq.gz" "$4" > $O/$1.map
  echo $O/$1.map
}
while [ $# -ge 2 ]; do
  src=$1; ds=$2; shift 2
  r1=$src/${ds}_R1.fq.gz; r2=$src/${ds}_R2.fq.gz; S=${ds}_R.sam
  echo "===== $ds"
  $old --db $db -1 $r1 -2 $r2 -o $O/old_$ds -t 1 --no_qcmsa > $O/old_$ds.log 2>&1
  $new --db $db -1 $r1 -2 $r2 -o $O/full_$ds -t 1 --no_qcmsa --full_sam_header > $O/full_$ds.log 2>&1
  $new --db $db -1 $r1 -2 $r2 -o $O/new_$ds -t 1 --no_qcmsa > $O/new_$ds.log 2>&1
  echo "  --full_sam_header vs OLD, 1 thread, every output:"; same_outputs $O/old_$ds $O/full_$ds
  echo "  default vs OLD, 1 thread, every output but the SAM:"; same_outputs $O/old_$ds $O/new_$ds "\.sam$"
  cmp -s <(grep -v '^@' $O/old_$ds/$S) <(grep -v '^@' $O/new_$ds/$S) && echo "  records identical" || echo "  RECORDS DIFFER"
  grep '^@SQ' $O/new_$ds/$S | sort > $O/sq_new; grep '^@SQ' $O/old_$ds/$S | sort > $O/sq_old
  grep -v '^@' $O/new_$ds/$S | awk -F'\t' '{ print $3; if ($7 != "=" && $7 != "*") print $7 }' | sort -u > $O/used
  cut -f2 $O/sq_new | sed 's/^SN://' | sort > $O/listed
  echo "  @SQ lines: $(wc -l < $O/sq_new) of $(wc -l < $O/sq_old); not in the old header: $(comm -23 $O/sq_new $O/sq_old | wc -l); genes named by records: $(wc -l < $O/used); named but not listed: $(comm -23 $O/used $O/listed | wc -l); listed but not named: $(comm -13 $O/used $O/listed | wc -l)"
  echo "  header bytes: $(grep '^@' $O/old_$ds/$S | wc -c) -> $(grep '^@' $O/new_$ds/$S | wc -c)"
  for ext in gz zst; do
    m=$(map ${ext}1_$ds $src $ds $S.$ext)
    $new --db $db --map $m -t 1 --no_qcmsa > $O/${ext}1_$ds.log 2>&1; echo "  .sam.$ext, 1 thread: rc=$?"
    f=$O/${ext}1_$ds/alignments/$S.$ext
    if [ $ext = gz ]; then gzip -t $f && echo "    gzip -t ok"; zcat $f > $O/dec; else zstd -tq $f && echo "    zstd -t ok"; zstd -dcq $f > $O/dec; fi
    cmp -s $O/dec $O/new_$ds/$S && echo "    decompressed = the plain SAM" || echo "    DECOMPRESSED SAM DIFFERS"
    for p in "" .log .genes.log .gene.log; do cmp -s $O/new_$ds/${ds}_R.profile$p $O/${ext}1_$ds/${ds}_R.profile$p || echo "    profile$p DIFFERS"; done
    echo "    $(stat -c %s $f) bytes (plain $(stat -c %s $O/new_$ds/$S)); temporary files left: $(find $O/${ext}1_$ds -name '*partial*' | wc -l)"
    cp $O/${ext}1_$ds/${ds}_R.profile $O/${ext}1_${ds}.first.profile
    $new --db $db --map $m -t 1 --no_qcmsa > $O/${ext}1_${ds}_rerun.log 2>&1
    grep -q "All alignments are present" $O/${ext}1_${ds}_rerun.log && echo "    rerun: alignment skipped" || echo "    RERUN ALIGNED AGAIN"
    cmp -s $O/${ext}1_${ds}.first.profile $O/${ext}1_$ds/${ds}_R.profile && echo "    rerun: profile identical" || echo "    RERUN PROFILE DIFFERS"
  done
  for ext in gz zst; do
    m=$(map ${ext}8_$ds $src $ds $S.$ext)
    $new --db $db --map $m -t 8 --no_qcmsa > $O/${ext}8_$ds.log 2>&1
    f=$O/${ext}8_$ds/alignments/$S.$ext
    if [ $ext = gz ]; then zcat $f > $O/dec; else zstd -dcq $f > $O/dec; fi
    cmp -s <(grep -v '^@' $O/dec | sort) <(grep -v '^@' $O/new_$ds/$S | sort) && r="records identical" || r="RECORDS DIFFER"
    cmp -s <(grep '^@' $O/dec) <(grep '^@' $O/new_$ds/$S) && h="header identical" || h="HEADER DIFFERS"
    p=identical; cmp -s $O/new_$ds/${ds}_R.profile $O/${ext}8_$ds/${ds}_R.profile || p=DIFFERS
    echo "  .sam.$ext, 8 threads: $r (sorted), $h, profile $p"
  done
done
