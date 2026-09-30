# usage: tree.sh <msa> <outprefix> [extra iqtree args]  -- the justfile's strain-trees command, -T 2
msa=$1; out=$2; shift 2
mkdir -p $(dirname $out)
iqtree2 -s "$msa" -m GTR+G -B 1000 -T 2 --seqtype DNA --prefix "$out" -redo "$@" > "$out.stdout.log" 2>&1
echo "rc=$? $out"
grep -i "warn\|error\|identical\|gap\|ambig\|constant\|parsimony-informative\|singleton" "$out.stdout.log" | head -20
