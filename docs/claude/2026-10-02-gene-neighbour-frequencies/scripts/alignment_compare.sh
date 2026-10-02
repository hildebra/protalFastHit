#!/bin/bash
# Some test samples of a build_gtdb_database.py run profiled again against its training database, as the collector
# did, with and without --no_gene_neighbours: what the gene neighbours did on the reads (protal's log, summed over
# the samples) and the calls (truth-annotated profiles, with the models the build trained: the training database
# holds placeholders).
# usage: [BIN=build dir] [ARMS="with without"] [NAME=alignment_comparison] alignment_compare.sh OUT [THREADS]
#        [SAMPLE_REGEX]   -> OUT/$NAME/{with,without}/
set -euo pipefail
out=$1 threads=${2:-5} pattern=${3:-^(rl150_p100000|pb_b30000000|ont_b30000000)_}
B=${BIN:-$HOME/bprof/gnb/build}
cmp=$out/${NAME:-alignment_comparison}
rm -rf "$cmp" && mkdir -p "$cmp"
for arm in ${ARMS:-with without}; do
  dir=$cmp/$arm
  mkdir -p "$dir"
  { echo -e "#OUTPUT_DIR\t$dir"; grep '^#SAMPLEID' "$out/test/profile_all/samples.map"
    grep -v '^#' "$out/test/profile_all/samples.map" | grep -E "$pattern" | \
      awk -F'\t' -v d="$dir" 'BEGIN { OFS = "\t" } { $4 = d "/" $1 ".sam.zst"; $5 = d "/" $1; $6 = d "/" $1 ".profile"; print }'
  } > "$dir/samples.map"
  extra=$([ $arm = without ] && echo --no_gene_neighbours || true)
  /usr/bin/time -f "%e s, %M kB" -o "$dir/time" "$B/protal" --db "$out/training_db" --map "$dir/samples.map" -t "$threads" \
    --no_strains --no_qcmsa --model "$out/trained_model.xml" --model_se "$out/trained_model_se.xml" \
    --model_pb "$out/trained_model_pb.xml" --model_ont "$out/trained_model_ont.xml" $extra > "$dir/protal.log" 2>&1
done
python3 - "$cmp" <<'PY'
import csv, glob, os, re, sys
cmp = sys.argv[1]
rows = {}
for arm in [a for a in ("with", "without") if os.path.isdir(os.path.join(cmp, a))]:
    log = open(os.path.join(cmp, arm, "protal.log"), errors="replace").read()
    paired = sum(int(n) for n in re.findall(r"Gene neighbours: (\d+) fragments paired across", log))
    guided = sum(int(n) for n in re.findall(r"; (\d+) guided mates found on the gene next to", log))
    found = sum(int(n) for n in re.findall(r"Gene neighbours: (\d+) genes found on reads", log))
    print(f"## {arm} ({open(os.path.join(cmp, arm, 'time')).read().strip()}): {paired} fragments paired across two genes, "
          f"{guided} guided mates found on the next gene, {found} long-read genes found where the neighbours put them")
    for path in sorted(glob.glob(os.path.join(cmp, arm, "*.profile.truth_annotated"))):
        sample = os.path.basename(path).split(".profile")[0]
        kind = sample.split("_")[0]
        with open(path) as fh:
            table = list(csv.DictReader(fh, delimiter="\t"))
        tp = sum(r["truth"] == "1" and r["prediction"] == "1" for r in table)
        fp = sum(r["truth"] == "0" and r["prediction"] == "1" for r in table)
        fn = sum(r["truth"] == "1" and r["prediction"] == "0" for r in table)
        acc = rows.setdefault((kind, arm), [0, 0, 0, 0])
        acc[0] += 1; acc[1] += tp; acc[2] += fp; acc[3] += fn
print("| reads | arm | samples | TP | FP | FN | F1 |")
print("|---|---|---|---|---|---|---|")
for (kind, arm), (n, tp, fp, fn) in sorted(rows.items()):
    f1 = 2 * tp / (2 * tp + fp + fn) if tp else 0
    print(f"| {kind} | {arm} | {n} | {tp} | {fp} | {fn} | {f1:.4f} |")
PY
