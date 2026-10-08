#!/usr/bin/env python3
"""F1 of the r226 builds' presence models side by side, from each build's model_logs/summary.txt: per read type the
F1 with species held out (knob 0.5), on the independent test set and on the scenarios' hold-out samples (at the model's
knob), with the knob and the FP and FN behind each. The summaries changed format (v10-v11 without rates or knob, v12
without knob), so each line is parsed from its numbers' count.

    python3 f1_history.py --builds local/v10,local/v11,...,local/v18 > f1_history.txt
"""
import argparse
import os
import re
import sys

import pandas as pd

ROW = re.compile(r"^(pe|se|pb|ont)\s+(.+?)\s{2,}([\d.\s]+)$")
WANTED = ["species held out", "independent test set", "soil: hold-out", "soil_shallow: hold-out", "moderate: hold-out",
          "gut: hold-out", "host: hold-out"]


def parse(path):
    """{(read type, evaluated on): {knob, TP, FP, FN, F1}} of one summary.txt."""
    rows = {}
    with open(path) as fh:
        for line in fh:
            m = ROW.match(line.rstrip("\n"))
            if not m:
                continue
            rt, label, rest = m.group(1), m.group(2).strip(), m.group(3).split()
            if label not in WANTED:
                continue
            values = [float(x) for x in rest]
            if len(values) == 13:  # knob, taxa, TP, FP, TN, FN, sens, spec, precision, F1, FP rate, FN rate, FP/sample
                knob, counts, f1 = values[0], values[1:6], values[9]
            elif len(values) == 12:  # without the knob (v12)
                knob, counts, f1 = 0.5, values[0:5], values[8]
            elif len(values) == 10:  # without the knob and the rates (v10, v11)
                knob, counts, f1 = 0.5, values[0:5], values[8]
            else:
                continue
            _, tp, fp, _, fn = counts
            rows[(rt, label)] = {"knob": knob, "TP": int(tp), "FP": int(fp), "FN": int(fn), "F1": f1}
    return rows


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--builds", required=True, help="build folders, comma-separated, oldest first")
    opts = p.parse_args(argv)
    table = []
    for folder in opts.builds.split(","):
        name = os.path.basename(folder.rstrip("/"))
        for (rt, label), r in parse(os.path.join(folder, "model_logs", "summary.txt")).items():
            table.append({"build": name, "read type": rt, "evaluated on": label, **r})
    t = pd.DataFrame(table)
    pd.set_option("display.width", 250)
    order = [os.path.basename(f.rstrip("/")) for f in opts.builds.split(",")]
    for rt in ("pe", "se", "pb", "ont"):
        sub = t[t["read type"] == rt]
        if sub.empty:
            continue
        f1 = sub.pivot(index="evaluated on", columns="build", values="F1").reindex(index=WANTED, columns=order)
        print(f"\n# {rt}: F1 (species held out at 0.5; test and hold-outs at the model's knob)\n")
        print(f1.dropna(how="all").to_string(float_format=lambda x: f"{x:.4f}"))
        knob = sub[sub["evaluated on"] == "independent test set"].set_index("build")["knob"].reindex(order)
        print("\nknob: " + ", ".join(f"{b} {k:g}" for b, k in knob.items() if pd.notna(k)))
        errs = sub[sub["evaluated on"].isin(["independent test set", "soil: hold-out", "soil_shallow: hold-out"])]
        errs = errs.assign(fpfn=errs["FP"].astype(str) + "/" + errs["FN"].astype(str))
        print("\nFP/FN at the knob:")
        print(errs.pivot(index="evaluated on", columns="build", values="fpfn").reindex(columns=order).dropna(how="all")
              .to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
