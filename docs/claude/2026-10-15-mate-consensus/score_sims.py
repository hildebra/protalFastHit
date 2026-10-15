#!/usr/bin/env python3
"""The calls of run_sims.sh's runs against the simulated truth, per mode, point and sample.

  calls: a profile's species against its PROFILE_TRUTH (the species simulated); FN only of species the database has
         (genome2tiid.tsv). TP/FP/FN, F1.
  rows:  the .truth_annotated table (every taxon with records, truth 1/0 and the model's call at its knob, as the
         build's "independent test set" counts): TP/FP/FN over the rows; a truth species without a row is "unseen".

usage: score_sims.py RUNS_FOLDER DATABASE_FOLDER OUT_PREFIX MODE [MODE ...]
  writes OUT_PREFIX.samples.tsv (a row per mode, point and sample) and OUT_PREFIX.summary.tsv (per mode and point,
  and all points), and prints the summary.
Standard library only.
"""
import csv
import glob
import os
import sys


def species_of(taxonomy):
    return taxonomy.strip().split(";")[-1]


def db_species(folder):
    names = set()
    with open(os.path.join(folder, "genome2tiid.tsv")) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 4:
                names.add(species_of(f[3]))
    return names


def read_map(path):
    rows = []
    with open(path) as fh:
        header = None
        for line in fh:
            if line.startswith("#SAMPLEID"):
                header = line[1:].rstrip("\n").split("\t")
            elif not line.startswith("#") and line.strip():
                rows.append(dict(zip(header, line.rstrip("\n").split("\t"))))
    return rows


def profile_calls(path):
    calls = set()
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2 and f[0] != "?":
                calls.add(species_of(f[1]))
    return calls


def truth_species(path):
    with open(path) as fh:
        return {species_of(line) for line in fh if line.strip() and not line.startswith("#")}


def annotated(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def f1(tp, fp, fn):
    return 2 * tp / (2 * tp + fp + fn) if tp else 0.0


def main():
    runs, database, prefix = sys.argv[1:4]
    modes = sys.argv[4:]
    in_db = db_species(database)
    samples = []
    for mode in modes:
        for map_path in sorted(glob.glob(os.path.join(runs, mode, "*", "map.tsv"))):
            folder = os.path.dirname(map_path)
            point = os.path.basename(folder)
            for row in read_map(map_path):
                profile = os.path.join(folder, "profiles", row["PROFILE"])
                if not os.path.exists(profile):
                    profile = os.path.join(folder, row["PROFILE"])
                if not os.path.exists(profile):
                    print(f"missing {profile}", file=sys.stderr)
                    continue
                truth = truth_species(row["PROFILE_TRUTH"])
                truth_in_db = truth & in_db
                calls = profile_calls(profile)
                tp, fp, fn = len(calls & truth), len(calls - truth), len(truth_in_db - calls)
                table = annotated(profile + ".truth_annotated") if os.path.exists(profile + ".truth_annotated") else []
                rtp = sum(r["truth"] == "1" and r["prediction"] == "1" for r in table)
                rfp = sum(r["truth"] == "0" and r["prediction"] == "1" for r in table)
                rfn = sum(r["truth"] == "1" and r["prediction"] == "0" for r in table)
                with_row = {r["taxon_name"] for r in table}
                unseen = len(truth_in_db - with_row)
                samples.append({"mode": mode, "point": point, "sample": row["SAMPLEID"], "truth": len(truth),
                                "truth_in_db": len(truth_in_db), "TP": tp, "FP": fp, "FN": fn,
                                "F1": f"{f1(tp, fp, fn):.4f}", "rows": len(table), "row_TP": rtp, "row_FP": rfp,
                                "row_FN": rfn, "unseen": unseen})
    columns = ["mode", "point", "sample", "truth", "truth_in_db", "TP", "FP", "FN", "F1", "rows", "row_TP", "row_FP",
               "row_FN", "unseen"]
    with open(prefix + ".samples.tsv", "w") as fh:
        fh.write("\t".join(columns) + "\n")
        for s in samples:
            fh.write("\t".join(str(s[c]) for c in columns) + "\n")
    summary = []
    points = sorted({s["point"] for s in samples}) + ["all"]
    for mode in modes:
        for point in points:
            rows = [s for s in samples if s["mode"] == mode and (point == "all" or s["point"] == point)]
            if not rows:
                continue
            t = {k: sum(r[k] for r in rows) for k in ("TP", "FP", "FN", "row_TP", "row_FP", "row_FN", "unseen")}
            summary.append({"mode": mode, "point": point, "samples": len(rows), **t,
                            "F1": f"{f1(t['TP'], t['FP'], t['FN']):.4f}",
                            "row_F1": f"{f1(t['row_TP'], t['row_FP'], t['row_FN']):.4f}"})
    columns = ["mode", "point", "samples", "TP", "FP", "FN", "F1", "row_TP", "row_FP", "row_FN", "row_F1", "unseen"]
    with open(prefix + ".summary.tsv", "w") as fh:
        fh.write("\t".join(columns) + "\n")
        for s in summary:
            fh.write("\t".join(str(s[c]) for c in columns) + "\n")
    print("\t".join(columns))
    for s in summary:
        print("\t".join(str(s[c]) for c in columns))


if __name__ == "__main__":
    main()
