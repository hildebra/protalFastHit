#!/usr/bin/env bash
# What a sparse first stage would keep, at GTDB r226 scale on a cluster node (experiment only;
# docs/claude/2026-10-12-stage1-probe).
#
#   bash probe_cluster.sh PROTAL_CHECKOUT DATABASE.protal PE_R1.fq.gz PE_R2.fq.gz HIFI.fq.gz OUT_DIR
#
# Builds protal at the checkout's HEAD (Release) in OUT_DIR/src and OUT_DIR/build, then runs it with --stage1_probe on
# the whole paired-end sample (OUT_DIR/pe, prefix pe1) and on the HiFi reads (OUT_DIR/pb, prefix pb2), 32 threads,
# without strain MSAs, and writes probe_summary.py's summary of each (OUT_DIR/pe1_probe.md, OUT_DIR/pb2_probe.md) next
# to the logs. The probe changes no output; it adds misc/<sample>_stage1_probe*.tsv and _first_lookups_probe.tsv.
# Run it on a whole node, e.g.
#   sbatch --exclusive -c 32 --mem 100G --wrap "bash probe_cluster.sh ~/dev/protalFastHit DB.protal R1.fq.gz R2.fq.gz hifi.fq.gz runs/stage1_probe"
set -uo pipefail
here=$(cd "$(dirname "$0")" && pwd)
checkout=$1 db=$2 r1=$3 r2=$4 hifi=$5 out=$6
jobs=${SLURM_CPUS_PER_TASK:-$(nproc)}
threads=$(( jobs < 32 ? jobs : 32 ))
mkdir -p "$out"
out=$(cd "$out" && pwd)
log=$out/probe_cluster.log
exec > >(tee -a "$log") 2>&1
echo "date: $(date -Iseconds); host: $(hostname); job: ${SLURM_JOB_ID:-none}; checkout: $(git -C "$checkout" rev-parse --short HEAD)"
lscpu | grep -E "Model name|^CPU\(s\)|L3 cache"

if [ ! -x "$out/build/protal" ]; then
  rm -rf "$out/src"; mkdir -p "$out/src"
  git -C "$checkout" archive HEAD | tar -x -C "$out/src"
  { cmake -S "$out/src" -B "$out/build" -DCMAKE_BUILD_TYPE=Release && cmake --build "$out/build" -j "$jobs" --target protal; } \
      > "$out/build.log" 2>&1 || { echo "protal build failed:"; tail -30 "$out/build.log"; exit 1; }
fi
"$out/build/protal" --version 2>&1 | head -1
"$out/build/protal" --help 2>&1 | grep -q stage1_probe || { echo "this protal has no --stage1_probe"; exit 1; }

run() {  # name, then the read arguments
  local name=$1; shift
  rm -rf "$out/$name"
  echo "== $name"
  /usr/bin/time -v "$out/build/protal" --db "$db" "$@" -o "$out/$name" -t "$threads" --no_strains --stage1_probe \
      > "$out/$name.log" 2>&1
  echo "exit $?"
  grep -E "stage-1 probe|Aligning reads took|Maximum resident" "$out/$name.log" | cut -c1-800
}
run pe --prefix pe1 -1 "$r1" -2 "$r2"
run pb --prefix pb2 -1 "$hifi" --read_type pb

for s in pe:pe1 pb:pb2; do
  dir=${s%%:*} name=${s##*:}
  python3 "$here/probe_summary.py" --misc "$out/$dir/misc" --sample "$name" --profile "$out/$dir/$name.profile" \
      > "$out/${name}_probe.md" && echo "== $out/${name}_probe.md" && cat "$out/${name}_probe.md"
done
