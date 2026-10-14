#!/usr/bin/env python3
"""The mock samples' measure_performance.sh runs (measure.sh) side by side, by level: dbgates (the database's gates, the
timed run), C0, C4s and C4g (the --test_sparse runs).
  resources.tsv  per sample and level: wall and CPU seconds, peak memory, the index load and the alignment, records
                 written, the share of the reads given up, and the species against the dbgates run's profile
  each level's profiles copied into a run folder of the mock benchmark, MOCKS/<label>_<level>/<set>/profiles/ (with a
  versions.tsv), and scored against the truth tables by score_mocks.py score; its by-run tables copied here, and
  scores.tsv: per sample and level the calls (TP/FP/FN), F1, L1 and Bray-Curtis against the truth.
Standard library only.

usage: collect.py MEASURE_DIR MOCKS LABEL DATABASE OUT_DIR
"""
import csv
import os
import platform
import re
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.join(HERE, "..", "..", "..")
LEVELS = (("dbgates", "plain"), ("C0", "c0"), ("C4s", "c4s"), ("C4g", "c4g"))  # (run folder suffix, file tag)


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path, columns, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def took(log, stage):
    """Seconds of a stage's "<stage> took 1m 2s 345ms" line in a protal log, or NA."""
    if not os.path.exists(log):
        return "NA"
    with open(log, errors="replace") as fh:
        for line in fh:
            found = re.search(re.escape(stage) + r" took (.*)", line)
            if found:
                seconds = 0.0
                for value, unit in re.findall(r"(\d+)(ms|h|m|s)\b", found.group(1).split(" (")[0].split(" mean")[0]):
                    seconds += int(value) * {"ms": 0.001, "s": 1, "m": 60, "h": 3600}[unit]
                return f"{seconds:.2f}"
    return "NA"


def cpu(row):
    try:
        return f"{float(row['user_s']) + float(row['sys_s']):.2f}"
    except (KeyError, ValueError):
        return "NA"


def main():
    measure, mocks, label, database, out = sys.argv[1:6]
    os.makedirs(out, exist_ok=True)
    samples = []
    with open(os.path.join(measure, "samples.tsv")) as fh:
        for line in fh:
            s, name, read_type, _, _ = line.rstrip("\n").split("\t")
            samples.append({"set": s, "sample": name, "read_type": read_type})
    resources = []
    folders = {suffix: os.path.join(mocks, f"{label}_{suffix}") for suffix, _ in LEVELS}
    for suffix, folder in folders.items():
        if os.path.exists(folder):
            shutil.rmtree(folder)
    for s in samples:
        base = os.path.join(measure, s["sample"])
        runs = read_tsv(os.path.join(base, "runs.tsv")) if os.path.exists(os.path.join(base, "runs.tsv")) else []
        if not runs:
            print(f"{s['sample']}: no timed run", file=sys.stderr)
            continue
        plain = runs[0]
        name = plain["sample"]  # measure_performance.sh's name (pe1, pb1, ont1)
        rows = {"dbgates": {"wall_s": plain["wall_s"], "cpu_s": cpu(plain), "max_rss_gb": plain["max_rss_gb"],
                            "load_index_s": plain["load_index_s"], "aligning_s": plain["aligning_s"],
                            "written": plain["written"], "given_up_share": "NA", "species": "NA",
                            "species_lost": "NA", "species_gained": "NA", "bray_curtis_to_dbgates": "0"}}
        for r in read_tsv(os.path.join(base, "sparse.tsv")) if os.path.exists(os.path.join(base, "sparse.tsv")) else []:
            log = os.path.join(base, "logs", f"{name}.sparse{r['sparse_first']}.log")
            rows["C" + r["sparse_first"]] = {"wall_s": r["wall_s"], "cpu_s": cpu(r), "max_rss_gb": r["max_rss_gb"],
                                             "load_index_s": took(log, "Load Index"), "aligning_s": r["aligning_s"],
                                             "written": r["written"], "given_up_share": r["given_up_share"],
                                             "species": r["species"], "species_lost": r["species_lost"],
                                             "species_gained": r["species_gained"],
                                             "bray_curtis_to_dbgates": r["bray_curtis"]}
        for suffix, tag in LEVELS:
            if suffix in rows:
                resources.append({"set": s["set"], "sample": s["sample"], "read_type": s["read_type"], "level": suffix,
                                  **rows[suffix]})
            profile = os.path.join(base, "sparse", f"{name}.{tag}.profile")
            if not os.path.exists(profile):
                continue
            target = os.path.join(folders[suffix], s["set"], "profiles")
            os.makedirs(target, exist_ok=True)
            shutil.copy(profile, os.path.join(target, s["sample"] + ".profile"))
            new = not os.path.exists(os.path.join(folders[suffix], "versions.tsv"))
            with open(os.path.join(folders[suffix], "versions.tsv"), "a") as fh:
                if new:
                    fh.write("sample\tset\ttool\tversion\tdatabase\tdate\thost\n")
                fh.write(f"{s['sample']}\t{s['set']}\tprotal\t{label} {suffix}\t{os.path.abspath(database)}\t"
                         f"{time.strftime('%Y-%m-%d %H:%M:%S')}\t{platform.node()}\n")
    write_tsv(os.path.join(out, "resources.tsv"),
              ["set", "sample", "read_type", "level", "wall_s", "cpu_s", "max_rss_gb", "load_index_s", "aligning_s",
               "written", "given_up_share", "species", "species_lost", "species_gained", "bray_curtis_to_dbgates"],
              resources)
    # The truth: score_mocks.py score of the four run folders (the samples' table with the two extra D6331 pe samples).
    table = os.path.join(mocks, "samples_plus.tsv")
    with open(table, "w") as fh:
        for path in (os.path.join(REPO, "scripts", "mock_benchmark", "samples.tsv"),
                     os.path.join(HERE, "..", "2026-10-13-mock-local-test", "samples_extra.tsv")):
            with open(path) as src:
                fh.write(src.read())
    names = [os.path.basename(f) for f in folders.values() if os.path.isdir(f)]
    result = subprocess.run([sys.executable, os.path.join(REPO, "scripts", "mock_benchmark", "score_mocks.py"), "score",
                             "--reads", mocks, "--samples", table, "--runs", ",".join(names), "--no-sylph"],
                            capture_output=True, text=True)
    os.remove(table)
    if result.returncode:
        sys.exit(result.stderr)
    for f in ("summary_by_run.tsv", "samples_by_run.tsv"):
        shutil.copy(os.path.join(mocks, "scores", f), os.path.join(out, f))
    scored = read_tsv(os.path.join(out, "samples_by_run.tsv"))
    rows = []
    for r in scored:
        rows.append({"set": r["set"], "sample": r["sample"], "read_type": r["read_type"], "kind": r["kind"],
                     "level": r["run"].rsplit("_", 1)[-1], "called": r["called"], "tp": r["tp"], "fp": r["fp"],
                     "fn": r["fn"], "f1": r["f1"], "l1": r["l1"], "bray_curtis": r["bray_curtis"],
                     "false_species": r["false_species"]})
    order = {suffix: i for i, (suffix, _) in enumerate(LEVELS)}
    rows.sort(key=lambda r: (r["set"], r["sample"], order.get(r["level"], 9)))
    write_tsv(os.path.join(out, "scores.tsv"), ["set", "sample", "read_type", "kind", "level", "called", "tp", "fp", "fn",
                                                "f1", "l1", "bray_curtis", "false_species"], rows)
    for r in rows:
        print("\t".join(str(r[c]) for c in ("sample", "level", "called", "tp", "fp", "fn", "f1", "l1")))


if __name__ == "__main__":
    main()
