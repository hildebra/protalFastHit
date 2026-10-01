#!/usr/bin/env python3
"""score.py [BENCH_DIR] - detection, abundance, archaea, low-abundance taxa, false positives, time and memory of
every run of profile.sh. Writes BENCH_DIR/results/runs.tsv (one row per run), species.tsv (one row per run and
species present or reported) and summary.md (the tables of the report).

Truth: a sample's species are those of its simulator truth file (protal_goldstd/<sample>.profile_truth); their
abundance is the simulator's relative_abundance summed over the species' genomes (manifest.tsv). A long-read or
single-end sample has the truth of the paired-end community it replays. Species no database has (the world's
135 unknown ones, world/unknown.txt) and, against a missing-species database, the species it leaves out
(V071/heldout_species.txt) cannot be found: they are not counted as present, the others' true abundances are
renormalised over them, and reads of theirs that land on relatives make false positives.

Per run, at the model's knob (0.5): TP, FP, FN, precision, recall, F1; Bray-Curtis between the true and the
reported relative abundances (each renormalised to 1; a species missing on one side counts as 0); the median
|log2(reported / true)| over the true positives (both renormalised over them); the archaeal species' recall and
false positives; false positives per sample, the share of the reported abundance on them, and how many are
congeners of a species in the sample (of any species simulated, unknown and held-out ones included); wall time,
CPU and peak RSS (/usr/bin/time -v).
"""
import collections
import csv
import glob
import math
import os
import re
import statistics
import sys

B = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/bench071")
RUNS = os.path.join(B, "runs")
POINTS = [os.path.join(B, "samples", "points"), os.path.join(B, "samples_deep", "points")]
OUT = os.path.join(B, "results")
VERSIONS = {"v060": "0.6.0a", "v070": "0.7.0", "v070m008": "0.7.0 with margin 0.08", "v071": "0.7.1",
            "v071m070": "0.7.1 with 0.7.0's model"}
BINS = [(0, 1e-4, "< 0.01%"), (1e-4, 1e-3, "0.01-0.1%"), (1e-3, 1e-2, "0.1-1%"), (1e-2, 2, ">= 1%")]


def species_of(lineage):
    return lineage.strip().split(";")[-1].removeprefix("s__")


lineage = {}
with open(os.path.join(B, "world", "full", "simulation", "genomes.tsv")) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        lineage[species_of(row["gtdb_taxonomy"])] = row["gtdb_taxonomy"].split(";")
domain = {s: lin[0].removeprefix("d__") for s, lin in lineage.items()}
genus = {s: lin[5] for s, lin in lineage.items()}
with open(os.path.join(B, "world", "unknown.txt")) as fh:
    unknown = {species_of(line) for line in fh if line.strip()}
with open(os.path.join(B, "V071", "heldout_species.txt")) as fh:
    heldout = {species_of(line.split("\t")[0]) for line in fh if line.strip()}
for v in ("V070",):
    path = os.path.join(B, v, "heldout_species.txt")
    if os.path.exists(path):
        with open(path) as fh:
            other = {species_of(line.split("\t")[0]) for line in fh if line.strip()}
        if other != heldout:
            sys.exit(f"{v} held out other species than V071: the missing-species databases differ")


def point_dir(point):
    return next((os.path.join(p, point) for p in POINTS if os.path.isdir(os.path.join(p, point))), None)


def truth_of(sample):
    """({species simulated: true relative abundance over all simulated species})."""
    m = re.match(r"((?:pb|ont)_b\d+)_s_(\d+)$", sample)
    if m:
        with open(os.path.join(point_dir(m.group(1)), "sim", "samples.tsv")) as fh:
            community = next(r["community"] for r in csv.DictReader(fh, delimiter="\t") if r["sample"] == sample)
        return truth_of(community)
    point = sample.rsplit("_s_", 1)[0]
    with open(os.path.join(point_dir(point), "sim", "protal_goldstd", sample + ".profile_truth")) as fh:
        simulated = {species_of(line) for line in fh if line.strip()}
    abundance = collections.Counter()
    with open(os.path.join(point_dir(point), "sim", "manifest.tsv")) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["sample"] == sample:
                abundance[species_of(row["species"])] += float(row["relative_abundance"])
    total = sum(abundance[s] for s in simulated) or 1.0
    return {s: abundance[s] / total for s in simulated}


def reported(run_dir):
    """{species: abundance} of the run's profile (the lineage in a column starting with d__, then the abundance)."""
    paths = glob.glob(os.path.join(run_dir, "**", "*.profile"), recursive=True)
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
    with open(os.path.join(RUNS, run + ".load")) as fh:
        t["load_before"] = float(fh.read().split()[0])
    return t


def score(run, species_rows):
    variant, db, reads, sample = run.split(".", 3)
    simulated = truth_of(sample)
    findable = set(simulated) - unknown - (heldout if db == "missing" else set())
    total = sum(simulated[s] for s in findable) or 1.0
    true_ab = {s: simulated[s] / total for s in findable}
    got = reported(os.path.join(RUNS, run))
    tp, fp, fn = findable & set(got), set(got) - findable, findable - set(got)
    precision = len(tp) / (len(tp) + len(fp)) if tp or fp else float("nan")
    recall = len(tp) / len(findable) if findable else float("nan")
    f1 = 2 * precision * recall / (precision + recall) if tp else 0.0
    shown = sum(got.values()) or 1.0
    est = {s: v / shown for s, v in got.items()}
    bray = 0.5 * sum(abs(true_ab.get(s, 0.0) - est.get(s, 0.0)) for s in set(true_ab) | set(est))
    t_tp, e_tp = sum(true_ab[s] for s in tp), sum(est[s] for s in tp)
    log2 = [abs(math.log2((est[s] / e_tp) / (true_ab[s] / t_tp))) for s in tp if est[s] > 0 and true_ab[s] > 0]
    sample_genera = {genus.get(s) for s in simulated}
    fp_congeners = sum(1 for s in fp if genus.get(s) in sample_genera)
    archaea = {s for s in findable if domain.get(s) == "Archaea"}
    for s in findable | set(got):
        species_rows.append(dict(run=run, variant=variant, db=db, reads=reads, sample=sample, species=s,
                                 domain=domain.get(s, "?"), true=round(true_ab.get(s, 0.0), 8),
                                 reported=round(est.get(s, 0.0), 8), present=int(s in findable), called=int(s in got)))
    point = re.sub(r"_s_\d+$", "", sample)
    return dict(run=run, variant=variant, db=db, reads=reads, scenario=point, sample=sample, species=len(findable),
                TP=len(tp), FP=len(fp), FN=len(fn), precision=round(precision, 4), recall=round(recall, 4),
                F1=round(f1, 4), bray_curtis=round(bray, 4),
                median_abs_log2=round(statistics.median(log2), 3) if log2 else "",
                archaea=len(archaea), archaea_TP=len(archaea & tp),
                archaea_FP=sum(1 for s in fp if domain.get(s) == "Archaea"),
                FP_abundance=round(sum(est[s] for s in fp), 4), FP_congeners=fp_congeners,
                **{k: (round(v, 2) if isinstance(v, float) else v) for k, v in timing(run).items()})


def mean(xs):
    xs = [x for x in xs if x != "" and not (isinstance(x, float) and math.isnan(x))]
    return statistics.mean(xs) if xs else float("nan")


def table(title, header, rows):
    lines = [f"\n### {title}\n", "| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    lines += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return lines


def fmt(x, digits=3):
    return "-" if isinstance(x, float) and math.isnan(x) else f"{x:.{digits}f}"


def main():
    runs = sorted(os.path.basename(p)[:-len(".done")] for p in glob.glob(os.path.join(RUNS, "*.done")))
    species_rows, rows = [], []
    for r in runs:
        rows.append(score(r, species_rows))
    os.makedirs(OUT, exist_ok=True)
    for name, data in (("runs.tsv", rows), ("species.tsv", species_rows)):
        with open(os.path.join(OUT, name), "w") as fh:
            columns = list(data[0])
            fh.write("\t".join(columns) + "\n")
            fh.writelines("\t".join(str(d[c]) for c in columns) + "\n" for d in data)

    order = [v for v in VERSIONS if any(r["variant"] == v for r in rows)]
    groups = collections.defaultdict(list)
    for r in rows:
        groups[(r["reads"], r["scenario"], r["db"], r["variant"])].append(r)
    scenarios = sorted({(r["reads"], r["scenario"]) for r in rows},
                       key=lambda x: (["pe", "se", "pb", "ont"].index(x[0]), x[1]))
    out = ["# Summary (score.py)"]
    for db in ("full", "missing"):
        body = []
        for reads, scenario in scenarios:
            cells = []
            for v in order:
                g = groups.get((reads, scenario, db, v))
                cells.append(f"{mean(x['F1'] for x in g):.3f} / {mean(x['bray_curtis'] for x in g):.3f}" if g else "-")
            body.append([reads, scenario] + cells)
        out += table(f"F1 / Bray-Curtis, {db} database (mean over samples)", ["reads", "scenario"] +
                     [VERSIONS[v] for v in order], body)
    for db in ("full", "missing"):
        body = []
        for reads, scenario in scenarios:
            cells = []
            for v in order:
                g = groups.get((reads, scenario, db, v))
                if not g:
                    cells.append("-")
                    continue
                cells.append(f"{mean(x['FP'] for x in g):.1f} / {mean(x['precision'] for x in g):.3f} / "
                             f"{mean(x['FP_abundance'] for x in g):.4f} / {sum(x['FP_congeners'] for x in g)}"
                             f" of {sum(x['FP'] for x in g)}")
            body.append([reads, scenario] + cells)
        out += table(f"False positives per sample / precision / share of the reported abundance / congeners of a "
                     f"species in the sample, {db} database", ["reads", "scenario"] + [VERSIONS[v] for v in order], body)
    for db in ("full", "missing"):
        body = []
        for reads, scenario in scenarios:
            cells = []
            for v in order:
                g = groups.get((reads, scenario, db, v))
                if not g:
                    cells.append("-")
                    continue
                n, tp, fpa = sum(x["archaea"] for x in g), sum(x["archaea_TP"] for x in g), sum(x["archaea_FP"] for x in g)
                cells.append(f"{tp}/{n} = {tp / n:.2f}, FP {fpa}" if n else f"0 present, FP {fpa}")
            body.append([reads, scenario] + cells)
        out += table(f"Archaea: found of present (recall), archaeal false positives, {db} database",
                     ["reads", "scenario"] + [VERSIONS[v] for v in order], body)
    # Low-abundance taxa: recall and abundance error by true relative abundance, per read type and depth.
    by = collections.defaultdict(list)
    for s in species_rows:
        if s["present"]:
            by[(s["reads"], re.sub(r"_s_\d+$", "", s["sample"]), s["db"], s["variant"])].append(s)
    for db in ("full", "missing"):
        body = []
        for reads, scenario in scenarios:
            for low, high, label in BINS:
                cells = []
                for v in order:
                    xs = [s for s in by.get((reads, scenario, db, v), []) if low <= s["true"] < high]
                    if not xs:
                        cells.append("-")
                        continue
                    called = [s for s in xs if s["called"] and s["reported"] > 0]
                    err = statistics.median(abs(math.log2(s["reported"] / s["true"])) for s in called) if called else float("nan")
                    cells.append(f"{len(called) / len(xs):.2f} / {fmt(err, 2)} (n={len(xs)})")
                if any(c != "-" for c in cells):
                    body.append([reads, scenario, label] + cells)
        out += table(f"By true relative abundance: recall / median abs log2(reported / true) of those found, {db} database",
                     ["reads", "scenario", "abundance"] + [VERSIONS[v] for v in order], body)
    body = []
    for reads, scenario in scenarios:
        cells = []
        for v in order:
            g = groups.get((reads, scenario, "full", v))
            cells.append(f"{mean(x['wall_s'] for x in g):.1f} s / {mean(x['cpu_s'] for x in g):.0f} s / "
                         f"{mean(x['peak_rss_gb'] for x in g):.2f} GB" if g else "-")
        body.append([reads, scenario] + cells)
    out += table("Wall time / CPU time / peak RSS, full database (mean over samples)",
                 ["reads", "scenario"] + [VERSIONS[v] for v in order], body)
    loads = [r["load_before"] for r in rows]
    out.append(f"\nLoad average before the runs: median {statistics.median(loads):.2f}, max {max(loads):.2f}.")
    text = "\n".join(out) + "\n"
    with open(os.path.join(OUT, "summary.md"), "w") as fh:
        fh.write(text)
    print(text)


if __name__ == "__main__":
    main()
