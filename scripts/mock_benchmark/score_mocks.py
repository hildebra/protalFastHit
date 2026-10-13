#!/usr/bin/env python3
"""Run protal on the mock-community samples that download_mocks.py fetched, and score the profiles against the truth
tables (truth/, prepare_tables.py) and against each other (samples.tsv's groups).

  run    writes a protal map for each set (or each chunk of a set: --chunks N, the chunk from --chunk or SLURM's
         $SLURM_ARRAY_TASK_ID) and runs protal on it: OUT/protal/<set>[/chunk_<i>]/ (map.tsv, protal.log, profiles/,
         misc/). A sample whose profile exists is not run again (--rerun: all). protal runs with --no_strains unless
         --strains; --protal-args adds options (quoted).
  score  reads every sample's <sample>.profile under OUT/protal and writes OUT/scores/:
           samples.tsv     per sample: the species called, protal's unknown share ("?") and the one the truth expects,
                           and by the truth's kind
                             abundance  true and false positives and negatives, precision, recall, F1, recall of the
                                        species at >= 0.1 % and >= 1 % of the placed cells, the share of the profile on
                                        false positives, and the L1 distance and Bray-Curtis dissimilarity of the
                                        profile and the truth (both over the placed species, adding up to 1)
                             closed     the same, without the distances (members' shares unknown)
                             dropout    the strains left out that were called, the ones kept that were missed
                             negative   the species called (all false)
           species.tsv     per sample and species of the truth or the profile: truth share, profile share, status
           replicates.tsv  per pair of samples of a group: Bray-Curtis of their profiles, Jaccard of their calls
           groups.tsv      per group: samples, mean and largest Bray-Curtis, mean and smallest Jaccard
           summary.tsv     per set, read type and truth kind: medians of the above (printed too)
         Shares compare cells: protal's abundances are shares of genomes, the truth's genome copies (prepare_tables.py).

Usage on a SLURM cluster (r226 needs ~40 GB of memory):
  sbatch -c 32 --mem 64G --wrap "python3 scripts/mock_benchmark/score_mocks.py run --reads /scratch/mocks \\
      --db /path/protal_db -t 32 --sets d6331_dna,meslier,d6331_hifi"
  sbatch -c 16 --mem 48G --array 0-9 --wrap "python3 scripts/mock_benchmark/score_mocks.py run \\
      --reads /scratch/mocks --db /path/protal_db -t 16 --sets hcom2_invitro --chunks 10"
  python3 scripts/mock_benchmark/score_mocks.py score --reads /scratch/mocks
Standard library only.
"""
import argparse
import collections
import csv
import glob
import itertools
import os
import shlex
import statistics
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
NAN = float("nan")
THRESHOLDS = (1e-3, 1e-2)
SAMPLE_COLUMNS = ["set", "sample", "read_type", "truth", "kind", "group", "reads", "bases", "cut", "called", "unknown",
                  "unknown_expected", "tp", "fp", "fn", "precision", "recall", "f1", "recall_0.1pct", "recall_1pct",
                  "fp_share", "l1", "bray_curtis", "absent_called", "present_missed", "false_species"]


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path, columns, rows):
    with open(path + ".part", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: fmt(v) for k, v in r.items()})
    os.replace(path + ".part", path)


def fmt(value):
    if isinstance(value, float):
        return "NA" if value != value else f"{value:.6g}"
    return value


def chosen(rows, sets):
    if not sets:
        return rows
    known = set(r["set"] for r in rows)
    unknown = [s for s in sets if s not in known]
    if unknown:
        sys.exit(f"unknown set(s): {', '.join(unknown)}")
    return [r for r in rows if r["set"] in sets]


# ---------------------------------------------------------------------------------------------------------------
# run
# ---------------------------------------------------------------------------------------------------------------

def profiles_of(out):
    """{sample: path of its .profile} under OUT/protal."""
    found = {}
    for path in glob.glob(os.path.join(out, "protal", "**", "profiles", "*.profile"), recursive=True):
        found[os.path.basename(path)[:-len(".profile")]] = path
    return found


def run(args):
    downloaded_path = os.path.join(args.reads, "reads", "downloaded.tsv")
    if not os.path.exists(downloaded_path):
        sys.exit(f"{downloaded_path} not found: run download_mocks.py first")
    samples = chosen(read_tsv(downloaded_path), [s for s in args.sets.split(",") if s])
    out = os.path.abspath(args.out or args.reads)
    done = profiles_of(out)
    chunk = args.chunk if args.chunk is not None else int(os.environ.get("SLURM_ARRAY_TASK_ID", "0"))
    if args.chunks > 1 and not 0 <= chunk < args.chunks:
        sys.exit(f"--chunk {chunk} is not in 0..{args.chunks - 1}")
    failed = []
    for set_name in dict.fromkeys(s["set"] for s in samples):
        members = [s for s in samples if s["set"] == set_name]
        if args.chunks > 1:
            members = members[chunk::args.chunks]
        todo = [s for s in members if args.rerun or s["sample"] not in done]
        if not todo:
            print(f"{set_name}: all {len(members)} samples have profiles", file=sys.stderr)
            continue
        folder = os.path.join(out, "protal", set_name)
        if args.chunks > 1:
            folder = os.path.join(folder, f"chunk_{chunk}")
        os.makedirs(folder, exist_ok=True)
        map_path = os.path.join(folder, "map.tsv")
        with open(map_path, "w") as fh:
            fh.write(f"#OUTPUT_DIR\t{folder}\n#SAMPLEID\tPREFIX\tFIRST\tSECOND\tREAD_TYPE\n")
            for s in todo:
                fh.write(f"{s['sample']}\t{s['sample']}\t{s['first']}\t{s['second']}\t{s['read_type']}\n")
        command = [args.protal, "--db", args.db, "--map", map_path, "-t", str(args.threads)]
        if not args.strains:
            command.append("--no_strains")
        command.append("--force")  # a SAM without its profile may be a killed run's
        command += shlex.split(args.protal_args)
        print(f"{set_name}: {len(todo)} samples: {' '.join(shlex.quote(c) for c in command)}", file=sys.stderr)
        if args.map_only:
            continue
        start = time.time()
        with open(os.path.join(folder, "protal.log"), "a") as log:
            log.write(f"# {time.strftime('%Y-%m-%d %H:%M:%S')} {' '.join(shlex.quote(c) for c in command)}\n")
            log.flush()
            code = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT).returncode
        missing = [s["sample"] for s in todo if not os.path.exists(os.path.join(folder, "profiles",
                                                                                  s["sample"] + ".profile"))]
        print(f"{set_name}: protal exit {code} after {time.time() - start:.0f} s; "
              f"{len(todo) - len(missing)} of {len(todo)} profiles written (log: {folder}/protal.log)", file=sys.stderr)
        if code or missing:
            failed.append(set_name)
    if failed:
        sys.exit(f"failed or incomplete: {', '.join(failed)}")


# ---------------------------------------------------------------------------------------------------------------
# score
# ---------------------------------------------------------------------------------------------------------------

def read_profile(path):
    """({GTDB species: share}, unknown share or NaN) of a protal .profile."""
    species = {}
    unknown = NAN
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3 or line.startswith("#"):
                continue
            try:
                share = float(f[2])
            except ValueError:
                continue
            if f[0] == "?" and f[1] == "?":
                unknown = share
                continue
            name = f[1].split(";")[-1]
            species[name] = species.get(name, 0.0) + share
    return species, unknown


def read_truth(path):
    """(kind, {species: share, NaN for present without a share, 0 for must be absent}, share of '?')."""
    kind = None
    rows = {}
    unknown = 0.0
    with open(path) as fh:
        lines = [l.rstrip("\n").split("\t") for l in fh if l.strip()]
    for f in lines:
        if f[0] == "#kind":
            kind = f[1]
    body = [f for f in lines if not f[0].startswith("#")]
    for f in body[1:]:
        share = NAN if f[1] == "NA" else float(f[1])
        if f[0] == "?":
            unknown = share
        else:
            rows[f[0]] = share
    if kind not in ("abundance", "closed", "dropout", "negative"):
        sys.exit(f"{path}: no #kind line of abundance, closed, dropout or negative")
    return kind, rows, unknown


def normalised(shares):
    total = sum(v for v in shares.values() if v == v and v > 0)
    return {k: v / total for k, v in shares.items() if v == v and v > 0} if total > 0 else {}


def bray_curtis(a, b):
    """Bray-Curtis dissimilarity of two compositions adding up to 1 (half their L1 distance)."""
    if not a and not b:
        return 0.0
    if not a or not b:
        return 1.0
    return sum(abs(a.get(k, 0.0) - b.get(k, 0.0)) for k in set(a) | set(b)) / 2


def ratio(n, d):
    return n / d if d else NAN


def score_sample(profile, unknown, kind, truth, truth_unknown):
    """The scores of one profile against its truth (columns of SAMPLE_COLUMNS), and its species rows."""
    called = set(profile)
    p = normalised(profile)
    s = {"called": len(called), "unknown": unknown, "unknown_expected": truth_unknown if kind else NAN}
    species_rows = []
    if kind in ("abundance", "closed"):
        members = {k for k, v in truth.items() if v != v or v > 0}
        tp, fp, fn = called & members, called - members, members - called
        s.update(tp=len(tp), fp=len(fp), fn=len(fn), precision=ratio(len(tp), len(called)),
                 recall=ratio(len(tp), len(members)), fp_share=sum(p.get(k, 0.0) for k in fp),
                 false_species=",".join(sorted(fp, key=lambda k: -p.get(k, 0.0))[:10]))
        s["f1"] = ratio(2 * s["tp"], 2 * s["tp"] + s["fp"] + s["fn"])
        if kind == "abundance":
            t = normalised(truth)
            for threshold, column in zip(THRESHOLDS, ("recall_0.1pct", "recall_1pct")):
                big = {k for k, v in t.items() if v >= threshold}
                s[column] = ratio(len(big & called), len(big))
            s["l1"] = sum(abs(p.get(k, 0.0) - t.get(k, 0.0)) for k in set(p) | set(t))
            s["bray_curtis"] = bray_curtis(p, t)
        for k in sorted(members | called):
            status = "TP" if k in tp else ("FP" if k in fp else "FN")
            species_rows.append({"species": k, "truth_share": truth.get(k, 0.0) if kind == "abundance"
                                 else ("member" if k in members else "-"), "profile_share": p.get(k, 0.0),
                                 "status": status})
    elif kind == "dropout":
        absent = {k for k, v in truth.items() if v == 0}
        present = {k for k, v in truth.items() if v != v or v > 0}
        s["absent_called"] = ",".join(sorted(absent & called)) or "-"
        s["present_missed"] = ",".join(sorted(present - called)) or "-"
        for k in sorted(absent | present):
            status = ("absent, called" if k in called else "absent") if k in absent else (
                "present, called" if k in called else "present, missed")
            species_rows.append({"species": k, "truth_share": 0.0 if k in absent else "kept",
                                 "profile_share": p.get(k, 0.0), "status": status})
    elif kind == "negative":
        s.update(fp=len(called), false_species=",".join(sorted(called, key=lambda k: -p.get(k, 0.0))[:10]))
        for k in sorted(called):
            species_rows.append({"species": k, "truth_share": 0.0, "profile_share": p.get(k, 0.0), "status": "FP"})
    return s, species_rows


def median(values):
    values = [v for v in values if isinstance(v, (int, float)) and v == v]
    return statistics.median(values) if values else NAN


def score(args):
    out = os.path.abspath(args.out or args.reads)
    samples = chosen(read_tsv(args.samples), [s for s in args.sets.split(",") if s])
    downloaded = {}
    path = os.path.join(args.reads, "reads", "downloaded.tsv")
    if os.path.exists(path):
        downloaded = {r["sample"]: r for r in read_tsv(path)}
    profiles = profiles_of(out)
    truths = {}
    sample_rows, species_rows, compositions = [], [], {}
    for s in samples:
        if s["sample"] not in profiles:
            continue
        profile, unknown = read_profile(profiles[s["sample"]])
        kind, truth, truth_unknown = None, {}, NAN
        if s["truth"] != "-":
            if s["truth"] not in truths:
                truths[s["truth"]] = read_truth(os.path.join(args.truth, s["truth"] + ".tsv"))
            kind, truth, truth_unknown = truths[s["truth"]]
        scores, rows = score_sample(profile, unknown, kind, truth, truth_unknown)
        d = downloaded.get(s["sample"], {})
        row = {k: NAN for k in SAMPLE_COLUMNS}
        row.update({"set": s["set"], "sample": s["sample"], "read_type": s["read_type"], "truth": s["truth"],
                    "kind": kind or "-", "group": s["group"], "reads": d.get("reads", "NA"),
                    "bases": d.get("bases", "NA"), "cut": d.get("cut", "NA"), "absent_called": "-",
                    "present_missed": "-", "false_species": "-"})
        row.update(scores)
        sample_rows.append(row)
        for r in rows:
            species_rows.append({"set": s["set"], "sample": s["sample"], **r})
        compositions[s["sample"]] = (s, normalised(profile))
    if not sample_rows:
        sys.exit(f"no profiles of the chosen samples under {out}/protal")

    pairs = []
    by_group = collections.defaultdict(list)
    for name, (s, _) in compositions.items():
        if s["group"] != "-":
            by_group[(s["set"], s["group"])].append(name)
    group_rows = []
    for (set_name, group), names in sorted(by_group.items()):
        if len(names) < 2:
            continue
        bcs, jaccards = [], []
        for a, b in itertools.combinations(sorted(names), 2):
            pa, pb = compositions[a][1], compositions[b][1]
            bc = bray_curtis(pa, pb)
            union = set(pa) | set(pb)
            jaccard = ratio(len(set(pa) & set(pb)), len(union)) if union else 1.0
            bcs.append(bc)
            jaccards.append(jaccard)
            pairs.append({"set": set_name, "group": group, "sample_a": a, "sample_b": b,
                          "read_type_a": compositions[a][0]["read_type"], "read_type_b": compositions[b][0]["read_type"],
                          "bray_curtis": bc, "jaccard": jaccard})
        group_rows.append({"set": set_name, "group": group, "samples": len(names), "pairs": len(bcs),
                           "mean_bray_curtis": statistics.fmean(bcs), "max_bray_curtis": max(bcs),
                           "mean_jaccard": statistics.fmean(jaccards), "min_jaccard": min(jaccards)})

    summary = []
    keys = collections.OrderedDict()
    for r in sample_rows:
        keys.setdefault((r["set"], r["read_type"], r["kind"]), []).append(r)
    for (set_name, read_type, kind), rows in keys.items():
        summary.append({"set": set_name, "read_type": read_type, "kind": kind, "samples": len(rows),
                        "called": median([r["called"] for r in rows]), "tp": median([r["tp"] for r in rows]),
                        "fp": median([r["fp"] for r in rows]), "fn": median([r["fn"] for r in rows]),
                        "precision": median([r["precision"] for r in rows]),
                        "recall": median([r["recall"] for r in rows]), "f1": median([r["f1"] for r in rows]),
                        "recall_1pct": median([r["recall_1pct"] for r in rows]),
                        "fp_share": median([r["fp_share"] for r in rows]), "l1": median([r["l1"] for r in rows]),
                        "unknown": median([r["unknown"] for r in rows]),
                        "dropouts_called": sum(r["absent_called"] not in ("-", "NA") for r in rows)})

    folder = os.path.join(out, "scores")
    os.makedirs(folder, exist_ok=True)
    write_tsv(os.path.join(folder, "samples.tsv"), SAMPLE_COLUMNS, sample_rows)
    write_tsv(os.path.join(folder, "species.tsv"), ["set", "sample", "species", "truth_share", "profile_share",
                                                    "status"], species_rows)
    write_tsv(os.path.join(folder, "replicates.tsv"), ["set", "group", "sample_a", "sample_b", "read_type_a",
                                                       "read_type_b", "bray_curtis", "jaccard"], pairs)
    write_tsv(os.path.join(folder, "groups.tsv"), ["set", "group", "samples", "pairs", "mean_bray_curtis",
                                                   "max_bray_curtis", "mean_jaccard", "min_jaccard"], group_rows)
    summary_columns = ["set", "read_type", "kind", "samples", "called", "tp", "fp", "fn", "precision", "recall", "f1",
                       "recall_1pct", "fp_share", "l1", "unknown", "dropouts_called"]
    write_tsv(os.path.join(folder, "summary.tsv"), summary_columns, summary)
    print("\t".join(summary_columns))
    for r in summary:
        print("\t".join(str(fmt(r[c])) for c in summary_columns))
    missing = len(samples) - len(sample_rows)
    print(f"{len(sample_rows)} samples scored{f' ({missing} without a profile)' if missing else ''}, "
          f"{len(pairs)} replicate pairs; tables in {folder}", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="command", required=True)
    for name in ("run", "score"):
        p = sub.add_parser(name)
        p.add_argument("--reads", required=True, help="download_mocks.py's output folder (-o)")
        p.add_argument("--out", default="", help="where protal/ and scores/ go (default: --reads)")
        p.add_argument("--sets", default="", help="comma-separated sets (default: all)")
        p.add_argument("--samples", default=os.path.join(HERE, "samples.tsv"))
        if name == "run":
            p.add_argument("--db", required=True, help="the protal database")
            p.add_argument("--protal", default="protal", help="the protal binary (default: protal on PATH)")
            p.add_argument("-t", "--threads", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", "8")))
            p.add_argument("--chunks", type=int, default=1, help="split each set into N chunks (SLURM arrays)")
            p.add_argument("--chunk", type=int, default=None, help="the chunk to run (default $SLURM_ARRAY_TASK_ID)")
            p.add_argument("--strains", action="store_true", help="write strain MSAs (default --no_strains)")
            p.add_argument("--rerun", action="store_true", help="run samples that have profiles again")
            p.add_argument("--protal-args", default="", help="more protal options, quoted")
            p.add_argument("--map-only", action="store_true", help="write the maps, do not run protal")
        else:
            p.add_argument("--truth", default=os.path.join(HERE, "truth"), help="folder of the truth tables")
    args = ap.parse_args()
    if args.command == "run":
        run(args)
    else:
        score(args)


if __name__ == "__main__":
    main()
