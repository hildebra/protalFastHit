#!/bin/bash
# PacBio reads of the same templates (make_templates.py) by pbsim3 as the collector made them until 2026-10-02
# (--strategy templ, the ERRHMM-SEQUEL error model, accuracy 0.999) and by hifi_reads.py (Q 30 +- 3): their time,
# base qualities and errors, and what protal makes of them: each sample profiled against a build's training database
# with the PacBio model that build trained on pbsim3's reads, the taxa's features and calls (truth-annotated).
# usage: compare_reads.sh BUILD_OUT PE_POINT BASES OUT [SAMPLES [SETUP]]   (SETUP of the templates' lengths: hifi:15000:3000:3)
#   (scripts of ~/bprof/gnb/src, protal of ~/bprof/gnb/build; e.g. ~/opw/b_gn2 rl100_p100000 30000000 ~/hifi 2)
set -euo pipefail
build=$1 point=$2 bases=$3 out=$4 samples=${5:-2} setup=${6:-hifi:15000:3000:3}
R=$HOME/bprof/gnb/src env=$HOME/micromamba/envs/protal-db-build P=$HOME/bprof/gnb/build/protal
here=$(cd "$(dirname "$0")" && pwd)
py=$env/bin/python3
mkdir -p "$out"
"$py" "$here/make_templates.py" "$R/scripts" "$build" "$point" "$bases" "$out/templates" "$samples" "$setup"
for t in "$out"/templates/*.fa; do
  s=$(basename "$t" .fa)
  mkdir -p "$out/pbsim/$s" "$out/hifi"
  /usr/bin/time -f "%e %U %S" -o "$out/pbsim/$s.time" "$env/bin/pbsim" --strategy templ --method errhmm \
    --errhmm "$env/data/ERRHMM-SEQUEL.model" --template "$t" --accuracy-mean 0.999 --seed 7 \
    --prefix "$out/pbsim/$s/r" --id-prefix r > "$out/pbsim/$s.log" 2>&1
  mv "$out/pbsim/$s"/r*.fq.gz "$out/pbsim/$s.fq.gz" 2>/dev/null || gzip -c "$out/pbsim/$s"/r*.fq > "$out/pbsim/$s.fq.gz"
  /usr/bin/time -f "%e %U %S" -o "$out/hifi/$s.time" "$py" "$R/scripts/hifi_reads.py" --templates "$t" \
    --out "$out/hifi/$s.fq.gz" --seed 7 > /dev/null
  truth=$build/test/points/$point/sim/protal_goldstd/$s.profile_truth
  for arm in pbsim hifi; do
    "$P" --db "$build/training_db" --model_pb "$build/trained_model_pb.xml" --read_type pb -1 "$out/$arm/$s.fq.gz" \
      --prefix "$s" -o "$out/$arm/protal" --profile_truth "$truth" --no_strains --no_qcmsa -t 5 \
      > "$out/$arm/$s.protal.log" 2>&1
  done
done
"$py" - "$out" <<'PY'
import csv, glob, gzip, os, statistics, sys
out = sys.argv[1]
def quals(path, limit=2000):
    q, n = [], 0
    with gzip.open(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 3:
                q.extend(ord(c) - 33 for c in line.rstrip("\n"))
                n += 1
                if n >= limit:
                    break
    return q
print("| sample | reads by | time s (wall, CPU) | Mb/s per CPU s | base Q median | share >= Q30 | taxa present: TP / FN | FP | present taxa's excess_median (median) |")
print("|---|---|---|---|---|---|---|---|---|")
for t in sorted(glob.glob(os.path.join(out, "templates", "*.fa"))):
    s = os.path.basename(t)[:-3]
    bases = sum(len(l.strip()) for l in open(t) if not l.startswith(">"))
    for arm in ("pbsim", "hifi"):
        wall, user, system = map(float, open(os.path.join(out, arm, s + ".time")).read().split()[-3:])
        q = sorted(quals(os.path.join(out, arm, s + ".fq.gz")))
        rows = list(csv.DictReader(open(os.path.join(out, arm, "protal", s + ".profile.truth_annotated")), delimiter="\t"))
        tp = sum(r["truth"] == "1" and r["prediction"] == "1" for r in rows)
        fn = sum(r["truth"] == "1" and r["prediction"] == "0" for r in rows)
        fp = sum(r["truth"] == "0" and r["prediction"] == "1" for r in rows)
        excess = [float(r["excess_median"]) for r in rows if r["truth"] == "1"]
        print(f"| {s} | {arm} | {wall:.1f}, {user + system:.1f} | {bases / 1e6 / (user + system):.1f} | {q[len(q) // 2]} | "
              f"{sum(v >= 30 for v in q) / len(q):.3f} | {tp} / {fn} | {fp} | {statistics.median(excess):.4f} |")
# The HiFi reads' quality (Phred of their bases' mean error probability) by read length.
import math
bins = ((0, 5000), (5000, 10000), (10000, 15000), (15000, 20000), (20000, 25000), (25000, 10 ** 9))
by_bin = {b: [] for b in bins}
for path in glob.glob(os.path.join(out, "hifi", "*.fq.gz")):
    with gzip.open(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 3:
                quality = line.rstrip("\n")
                q = -10 * math.log10(sum(10 ** (-(ord(c) - 33) / 10) for c in quality) / len(quality))
                by_bin[next(b for b in bins if b[0] <= len(quality) < b[1])].append(q)
print()
print("| read length | reads | read Q: median | 10-90% |")
print("|---|---|---|---|")
for (low, high), qs in by_bin.items():
    if qs:
        qs.sort()
        print(f"| {low // 1000}-{high // 1000 if high < 10 ** 9 else ''} kb | {len(qs)} | {qs[len(qs) // 2]:.1f} | "
              f"{qs[len(qs) // 10]:.1f}-{qs[9 * len(qs) // 10]:.1f} |")
PY
