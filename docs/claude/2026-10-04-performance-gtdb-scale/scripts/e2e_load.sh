#!/usr/bin/env bash
# The binary gene table and the concurrent start-up end to end, on the v0.7.3 world: a copy of its database.protal upgraded
# with --compress_db (adds gene_table.bin); runs of the working tree (concurrent and --sequential_load) against the reference
# build on the original file; an unpacked folder; outputs compared; start-up timers over alternated rounds.
set -uo pipefail
W=$HOME/perf-gtdb; H=$W/head/build/protal; N=$W/work/build/protal; X=$W/e2e; rm -rf $X; mkdir -p $X/db $X/folder
cp $HOME/bench071/V073/protal_db/database.protal $X/db/
P=$HOME/bench071/samples/points; R1=$P/rl150_p500000/sim/reads/rl150_p500000_s_1_R1.fq.gz; R2=$P/rl150_p500000/sim/reads/rl150_p500000_s_1_R2.fq.gz
PB=$HOME/bench071/samples_lr073/points/pb_b90000000/sim/reads/pb_b90000000_s_1.fq.gz
echo "== --compress_db on the single file"
$N --compress_db --db $X/db/database.protal -t 6 --compress_level 9 2>&1 | grep -vE '^\s*$' | tail -6
echo "== again (current: kept)"
$N --compress_db --db $X/db/database.protal -t 6 2>&1 | tail -1
echo "== unpack"
$N --unpack_db --db $X/db/database.protal --unpack_dir $X/folder -t 6 2>&1 | grep -E 'Skip|Unpacked' | cut -c1-160
ls $X/folder | tr '\n' ' '; echo
run() { local name=$1 bin=$2 db=$3 o=$4; shift 4; rm -rf $o
  /usr/bin/time -f "%e s wall, %U s user" $bin --db $db "$@" --prefix s -o $o -t 6 --no_qcmsa > $o.log 2> $o.time
  echo "$name: $(tail -1 $o.time); $(grep -hE '^Run protal took' $o.log)"; }
same() { cmp -s <(zstd -dc $1/s.sam.zst | sort) <(zstd -dc $2/s.sam.zst | sort) && echo "    SAM identical" || echo "    SAM DIFFERS"
  diff -r -x '*_runtime.tsv' -x '*.sam.zst' -x '*.statistics.tsv' $1 $2 > $X/d.txt && echo "    other outputs identical" || { echo "    OUTPUTS DIFFER"; head -6 $X/d.txt; }; }
echo "== outputs: pe 500k"
run "head, original file" $H $HOME/bench071/V073/protal_db/database.protal $X/h -1 $R1 -2 $R2 --read_type pe
run "work, gene_table.bin, concurrent" $N $X/db/database.protal $X/c -1 $R1 -2 $R2 --read_type pe
run "work, gene_table.bin, --sequential_load" $N $X/db/database.protal $X/s -1 $R1 -2 $R2 --read_type pe --sequential_load
run "work, unpacked folder" $N $X/folder $X/f -1 $R1 -2 $R2 --read_type pe
for v in c s f; do echo "  head vs $v:"; same $X/h $X/$v; done
grep -hE '^Gene tables|took|^Note' $X/c.log | grep -vE 'Profiling sample' | sed 's/^/    c: /'
grep -hE '^Gene tables|took' $X/s.log | grep -vE 'Profiling sample' | sed 's/^/    s: /'
grep -hE '^Gene tables' $X/f.log | sed 's/^/    f: /'
diff <(grep -vE 'took|Gene tables|Output dir|sam file|profile file|output prefix|Freeing' $X/s.log) <(grep -vE 'took|Gene tables|Output dir|sam file|profile file|output prefix|Freeing' $X/c.log) > /dev/null && echo "  logs of concurrent and sequential identical but for timers" || { echo "  LOGS DIFFER"; diff <(grep -vE 'took|Gene tables|Output dir|sam file|profile file|output prefix|Freeing' $X/s.log) <(grep -vE 'took|Gene tables|Output dir|sam file|profile file|output prefix|Freeing' $X/c.log) | head -10; }
echo "== outputs: PacBio"
run "head pb" $H $HOME/bench071/V073/protal_db/database.protal $X/hp -1 $PB --read_type pb
run "work pb concurrent" $N $X/db/database.protal $X/cp -1 $PB --read_type pb
same $X/hp $X/cp
echo "== start-up, 3 alternated rounds (500k pairs, 6 threads, whole runs)"
for round in 1 2 3; do
  for v in concurrent sequential; do
    extra=""; [ $v = sequential ] && extra="--sequential_load"
    o=$X/t.$v.$round; rm -rf $o
    $N --db $X/db/database.protal -1 $R1 -2 $R2 --read_type pe --prefix s -o $o -t 6 --no_qcmsa $extra > $o.log 2>&1
    echo "$v $round: $(grep -hE '^(Loading the gene tables|Preload genomes|Loading the taxonomy|Load Index|Profiling|Run protal) took' $o.log | sed 's/ took /=/' | tr '\n' ';')"
  done
done
