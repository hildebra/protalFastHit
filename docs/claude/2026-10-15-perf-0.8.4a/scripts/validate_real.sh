#!/usr/bin/env bash
# The hashed key map on real reads: the D6331 HKU pe mock (1M pairs) and the D6331 HiFi mock against the reduced r226
# database, compact (default), --hashed_keys and --direct_keys: the same SAM records, each run's key map line and
# seeding time; then the whole unit suite and the e2e suite. Cores 0-3.
set -uo pipefail
W=~/hkeys
P=$W/build/protal
DB=~/Mocks/reduced_db_r226/database.protal
PE=~/Mocks/reads/d6331_pe_labs/d6331_hku_SRR23926922
HIFI=~/Mocks/reads/d6331_hifi/d6331_hifi_SRR13128012/SRR13128012.fq.gz
OUT=$W/real
mkdir -p $OUT
for kind in compact hashed direct; do
  case $kind in compact) extra=();; hashed) extra=(--hashed_keys);; direct) extra=(--direct_keys);; esac
  rm -rf $OUT/pe_$kind $OUT/hifi_$kind
  taskset -c 0-3 $P --db $DB -1 $PE/SRR23926922_1.fq.gz -2 $PE/SRR23926922_2.fq.gz -t 4 --no_strains --no_profile \
    --prefix pe_$kind -o $OUT/pe_$kind "${extra[@]}" > $OUT/pe_$kind.log 2>&1
  echo "pe $kind exit $?"
  taskset -c 0-3 $P --db $DB -1 $HIFI --read_type pb -t 4 --no_strains --no_profile \
    --prefix hifi_$kind -o $OUT/hifi_$kind "${extra[@]}" > $OUT/hifi_$kind.log 2>&1
  echo "hifi $kind exit $?"
done
for s in pe hifi; do
  for kind in compact hashed direct; do
    f=$(find $OUT/${s}_$kind -name "*.sam*" | head -1)
    case $f in *.zst) zstd -dc $f;; *) cat $f;; esac | grep -v '^@' | cut -f1-9 | sort | md5sum | sed "s#-#$s $kind records#"
    grep -h "key map\|hashed key map, made\|Aligning reads took" $OUT/${s}_$kind.log | sed 's/.*file: [0-9.]* GB; //' | cut -c1-200
  done
done
cd $W/build/tests && taskset -c 0-3 ./protal_tests 2>&1 | tail -3
cd $W/src && PROTAL_TEST_DB=data/mini_db/protal_db PROTAL=$W/build/protal SIMULATE=$W/build/simulate_metagenomes PROTAL_TESTS_REQUIRED=1 \
  taskset -c 0-3 python3 -m unittest tests/e2e/test_protal_e2e.py 2>&1 | tail -4
echo VALIDATE_DONE
