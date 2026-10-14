#!/usr/bin/env bash
# Unique cores counted from the alignment (UniqueCores.h, 2026-10-14), part 1: the strain world's database built again by
# the new protal (~/strainfn/db_cores: unique_cores.bin, ancestry sites version 2), its reads aligned at 150 and 75 bp
# and by the grid at C 4 on its sparse index, and each record's ZS, ZU and ZT compared between the whole index and the
# sparse one (uc_tags.pl). The binaries and scripts are copied to ~/lfv2. 4 cores, niced.
set -euo pipefail
S=$(dirname "$(readlink -f "$0")")
V=~/lfv2
W=~/strainfn
mkdir -p $V/bin $V/src
cp ~/lone-build/build/protal ~/lone-build/build/simulate_metagenomes $V/bin/
rsync -a --delete ~/lone-build/src/scripts/ $V/src/scripts/
PROTAL=$V/bin/protal
if [ ! -s $W/db_cores/database.protal ]; then
  rm -rf $W/db_cores
  # db_grid's members (reference.fna, reference.map, the taxonomy, the tables a build writes again) unpacked, the files
  # its build kept beside the database copied.
  $PROTAL --unpack_db --db $W/db_grid --unpack_dir $W/db_cores > $W/db_cores_unpack.log 2>&1
  for f in $W/db_grid/*; do
    case $(basename $f) in database.protal|build.log|build.time) ;; *) cp -r $f $W/db_cores/ ;; esac
  done
  /usr/bin/time -v -o $W/db_cores/build.time taskset -c 0-3 nice -n 5 $PROTAL --build --no_profile -t 4 --db $W/db_cores \
      --reference $W/db_cores/reference.fna --full_reference $W/db_cores/full_reference.fna > $W/db_cores/build.log 2>&1
fi
echo "== db_cores: $(grep -m1 Elapsed $W/db_cores/build.time | sed 's/.*: //')"
grep -E "^Unique cores:|^Ancestry sites:" $W/db_cores/build.log | sed 's/^/  new: /'
grep -E "^Ancestry sites:" $W/db_grid/build.log | sed 's/^/  db_grid: /' || true
samples=$(cut -f1 $W/reads/samples.tsv)
align() {
  local name=$1 reads=$2; shift 2
  local R=$W/run_u_$name
  [ -s $R/done ] && return
  rm -rf $R; mkdir -p $R/alignments
  {
    printf '#OUTPUT_DIR\t%s\n#SAMPLEID\tFIRST\tSECOND\tSAM\tPREFIX\tPROFILE\tREAD_TYPE\n' $R
    for s in $samples; do
      printf '%s\t%s\t%s\t%s\t%s\t%s\tpe\n' $s $reads/${s}_1.fq.gz $reads/${s}_2.fq.gz $R/alignments/$s.sam.zst $R/$s $R/$s.profile
    done
  } > $R/samples.map
  taskset -c 0-3 nice -n 5 $PROTAL --db $W/db_cores --map $R/samples.map -t 4 --no_strains --no_qcmsa --no_profile "$@" > $R/protal.log 2>&1
  date > $R/done
}
align l150 $W/reads --sparse_first 0
align l75 $W/run_len75/reads --sparse_first 0
align g4 $W/reads --sparse_first 4 --sparse_first_pick grid
grep -m1 "^Unique cores:" $W/run_u_g4/protal.log
grep -c "seed expectation" $W/run_u_g4/protal.log || true
echo "== whole index against the grid's sparse index (counts from the alignment)"
perl $S/uc_tags.pl $W/run_u_l150 $W/run_u_g4
echo "== the same before (counts from the anchors, lfv_strain.sh's runs)"
perl $S/uc_tags.pl $W/run_v_l150 $W/run_v_g4
echo "== 150 against 75 bp (counts from the alignment)"
perl $S/uc_tags.pl $W/run_u_l150 $W/run_u_l75 | head -2
echo "== done $(date +%T)"
