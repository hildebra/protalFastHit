#!/bin/bash
# The trainer's previous-procedure study before and after it was made cheaper (2026-10-02): the paired-end training
# table of a build, as it is and 5 times as large (enlarge_table.sh), trained with --evaluation full by the trainer
# before (BEFORE: a scripts folder with the earlier random_forest_cmdline.py, which ran the study with --evaluation
# full) and after (--previous-procedure), 4 threads; the study's seconds and results from each metrics.json, and each
# training's CPU seconds (user + system, /usr/bin/time), which a busy machine changes less than its wall time.
# usage: previous_procedure_cost.sh TABLE BEFORE_DIR AFTER_SCRIPTS_DIR OUT
set -euo pipefail
table=$1 before=$2 after=$3 out=$4
here=$(cd "$(dirname "$0")" && pwd)
py=$HOME/micromamba/envs/protal-db-build/bin/python3
mkdir -p "$out"
cp "$table" "$out/x1.tsv"
bash "$here/enlarge_table.sh" "$table" 5 "$out/x5.tsv" > /dev/null
for size in x1 x5; do
  # The trainer before ran the study only with --evaluation full: there it is timed on its own (timing_seconds).
  [ -f "$out/${size}_before.metrics.json" ] || /usr/bin/time -f "%e %U %S" -o "$out/${size}_before.time" \
    "$py" "$before/random_forest_cmdline.py" --truth-file "$out/$size.tsv" \
    --output-prefix "$out/${size}_before" --evaluation full --threads 4 --seed 1 > "$out/${size}_before.log" 2>&1
  [ -f "$out/${size}_after.metrics.json" ] || /usr/bin/time -f "%e %U %S" -o "$out/${size}_after.time" \
    "$py" "$after/random_forest_cmdline.py" --truth-file "$out/$size.tsv" \
    --output-prefix "$out/${size}_after" --evaluation full --previous-procedure --threads 4 --seed 1 \
    > "$out/${size}_after.log" 2>&1
  # Without the study: what the rest of the training costs.
  [ -f "$out/${size}_none.metrics.json" ] || /usr/bin/time -f "%e %U %S" -o "$out/${size}_none.time" \
    "$py" "$after/random_forest_cmdline.py" --truth-file "$out/$size.tsv" \
    --output-prefix "$out/${size}_none" --evaluation full --threads 4 --seed 1 > "$out/${size}_none.log" 2>&1
done
"$py" - "$out" <<'PY'
import json, os, sys
out = sys.argv[1]
print("| table | trainer | rows | grid | study s | training s | training CPU s | the study's CPU s | previous: mtry "
      "| previous F1, species held out | this one F1 |")
print("|---|---|---|---|---|---|---|---|---|---|---|")
for size in ("x1", "x5"):
    for arm in ("before", "after"):
        m = json.load(open(os.path.join(out, f"{size}_{arm}.metrics.json")))
        pp, t = m["previous_procedure"], m["timing_seconds"]
        res = {(r["procedure"], r["judged on"]): r for r in pp["results"]}
        grid = (f"{len(pp['grid_values'])} values x {pp['grid_folds']} folds on {pp['grid_rows']} rows"
                if "grid_values" in pp else "20 values x 5 folds on all rows")
        prev = res.get(("previous", "species held out"), {}).get("F1")
        new = res.get(("this one", "species held out"), {}).get("F1")
        cpu = lambda run: sum(map(float, open(os.path.join(out, f"{size}_{run}.time")).read().split()[-2:]))
        print(f"| {size} | {arm} | {m['data']['taxa']} | {grid} | {t['previous_procedure']:.0f} | {t['total']:.0f} | "
              f"{cpu(arm):.0f} | {cpu(arm) - cpu('none'):.0f} | {pp['mtry']} | {prev:.4f} | {new:.4f} |")
PY
