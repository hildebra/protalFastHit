#!/usr/bin/env python3
"""score.py - detection, abundance, time and memory of every profiling run of profile.sh.

Usage: score.py [RUNS_DIR] [TEST_POINTS_DIR] > runs.tsv   (then summary tables on stderr)

Truth: a sample's species are those of its simulator truth file (protal_goldstd/<sample>.profile_truth,
the labels the models are trained on); their abundance is the simulator's relative_abundance summed
over the species' genomes (manifest.tsv), renormalised over the species. A long-read or single-end
sample has the truth of the paired-end community it replays.

Per run: TP, FP, FN, precision, recall, F1 of the species reported; Bray-Curtis dissimilarity between
the true and the reported relative abundances (each renormalised to 1; a species missing on one side
counts as 0); the median |log2(reported / true)| over the true positives (both renormalised over them);
wall time, user+system CPU and peak RSS of the whole run (/usr/bin/time -v), and protal's own index
load and alignment times where its log has them.
"""
import collections
import csv
import glob
import math
import os
import re
import statistics
import sys

RUNS = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/bench07/runs")
# Design points: the tuning world's test set, then those of $EXTRA_POINTS (colon-separated folders of
# points laid out alike, e.g. the deep sample of deep.sh).
POINTS = [os.path.expanduser(sys.argv[2] if len(sys.argv) > 2 else "~/tune/V2/test/points")] + \
    [os.path.expanduser(p) for p in os.environ.get("EXTRA_POINTS", "~/bench07/points").split(":") if p]


def point_dir(point):
    return next((os.path.join(p, point) for p in POINTS if os.path.isdir(os.path.join(p, point))),
                os.path.join(POINTS[0], point))


def species_of(lineage):
    return lineage.strip().split(";")[-1].removeprefix("s__")


def truth_of(sample):
    """(present species, {species: true relative abundance}) of a sample (see the docstring)."""
    m = re.match(r"((?:pb|ont)_b\d+)_s_(\d+)$", sample)
    if m:
        with open(os.path.join(point_dir(m.group(1)), "sim", "samples.tsv")) as fh:
            community = next(r["community"] for r in csv.DictReader(fh, delimiter="\t") if r["sample"] == sample)
        return truth_of(community)
    point = sample.rsplit("_s_", 1)[0]
    with open(os.path.join(point_dir(point), "sim", "protal_goldstd", sample + ".profile_truth")) as fh:
        present = {species_of(line) for line in fh if line.strip()}
    abundance = collections.Counter()
    with open(os.path.join(point_dir(point), "sim", "manifest.tsv")) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["sample"] == sample:
                abundance[row["species"]] += float(row["relative_abundance"])
    total = sum(abundance[s] for s in present) or 1.0
    return present, {s: abundance[s] / total for s in present}


def reported(run_dir):
    """{species: abundance} of the run's profile (lineage in a column starting with d__, then the
    abundance)."""
    paths = [p for p in glob.glob(os.path.join(run_dir, "**", "*.profile"), recursive=True)]
    if len(paths) != 1:
        raise SystemExit(f"{run_dir}: expected one .profile, found {paths}")
    result = {}
    with open(paths[0]) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            at = next((i for i, f in enumerate(fields) if f.startswith("d__")), None)
            if at is None or at + 1 >= len(fields):
                continue
            try:
                result[species_of(fields[at])] = float(fields[at + 1])
            except ValueError:
                continue
    return result


def timing(run):
    t = {}
    with open(os.path.join(RUNS, run + ".time")) as fh:
        for line in fh:
            key, _, value = line.strip().rpartition(": ")
            if key.startswith("Elapsed (wall clock)"):
                parts = [float(x) for x in value.split(":")]
                t["wall_s"] = sum(p * 60 ** i for i, p in enumerate(reversed(parts)))
            elif key == "Maximum resident set size (kbytes)":
                t["peak_rss_gb"] = int(value) / 1024 ** 2
            elif key in ("User time (seconds)", "System time (seconds)"):
                t["cpu_s"] = t.get("cpu_s", 0.0) + float(value)
    log = open(os.path.join(RUNS, run + ".log"), errors="replace").read()

    def took(label):
        m = re.search(label + r" took (?:(\d+)m )?(?:(\d+)s )?(?:(\d+)ms)?", log)
        if not m or not any(m.groups()):
            return ""
        mins, secs, ms = (int(g) if g else 0 for g in m.groups())
        return round(mins * 60 + secs + ms / 1000, 2)

    t["load_index_s"] = took("Load Index")
    t["align_s"] = took("Aligning reads")
    return t


def heldout():
    """Species the training database leaves out ($HELDOUT, heldout_species.txt): runs against it (variants
    with "ho") are scored on the other species, renormalised; reads of the others may land on relatives."""
    path = os.path.expanduser(os.environ.get("HELDOUT", "~/tune/V2/heldout_species.txt"))
    with open(path) as fh:
        return {line.split("\t")[0].strip().removeprefix("s__") for line in fh if line.strip()}


def score(run):
    variant, sample = run.split(".", 1)
    kind = "pe"
    if sample.split(".")[0] in ("se", "pb", "ont"):
        kind, sample = sample.split(".", 1)
    present, true_ab = truth_of(sample)
    if "ho" in variant:
        present = present - heldout()
        total = sum(true_ab[s] for s in present) or 1.0
        true_ab = {s: true_ab[s] / total for s in present}
    got = reported(os.path.join(RUNS, run))
    tp, fp, fn = len(present & set(got)), len(set(got) - present), len(present - set(got))
    precision = tp / (tp + fp) if tp + fp else float("nan")
    recall = tp / (tp + fn) if tp + fn else float("nan")
    f1 = 2 * precision * recall / (precision + recall) if tp else 0.0
    total = sum(got.values()) or 1.0
    est = {s: v / total for s, v in got.items()}
    bray = 0.5 * sum(abs(true_ab.get(s, 0.0) - est.get(s, 0.0)) for s in set(true_ab) | set(est))
    both = present & set(got)
    t_tp, e_tp = sum(true_ab[s] for s in both), sum(est[s] for s in both)
    log2 = [abs(math.log2((est[s] / e_tp) / (true_ab[s] / t_tp))) for s in both if est[s] > 0 and true_ab[s] > 0]
    point = re.sub(r"_s_\d+$", "", sample)
    return dict(run=run, variant=variant, reads=kind, scenario=point, sample=sample, species=len(present), TP=tp, FP=fp,
                FN=fn, precision=round(precision, 4), recall=round(recall, 4), F1=round(f1, 4),
                bray_curtis=round(bray, 4), median_abs_log2=round(statistics.median(log2), 3) if log2 else "",
                **{k: (round(v, 2) if isinstance(v, float) else v) for k, v in timing(run).items()})


def main():
    runs = sorted(os.path.basename(p)[:-len(".done")] for p in glob.glob(os.path.join(RUNS, "*.done"))
                  if not os.path.basename(p).startswith("strain"))
    rows = [score(r) for r in runs]
    columns = list(rows[0])
    print("\t".join(columns))
    for row in rows:
        print("\t".join(str(row[c]) for c in columns))
    # Means over a scenario's samples, per variant.
    groups = collections.defaultdict(list)
    for row in rows:
        groups[(row["scenario"], row["reads"], row["variant"])].append(row)
    keys = ["species", "precision", "recall", "F1", "bray_curtis", "median_abs_log2", "wall_s", "cpu_s", "peak_rss_gb",
            "load_index_s", "align_s"]
    sys.stderr.write("| scenario | reads | version | " + " | ".join(keys) + " |\n|" + "---|" * (len(keys) + 3) + "\n")
    for (scenario, reads, variant), group in sorted(groups.items()):
        vals = []
        for k in keys:
            xs = [g[k] for g in group if g[k] != ""]
            vals.append(f"{statistics.mean(xs):.3g}" if xs else "-")
        sys.stderr.write(f"| {scenario} | {reads} | {variant} | " + " | ".join(vals) + " |\n")


if __name__ == "__main__":
    main()
