#!/usr/bin/env python3
"""design.py - the samples and the tables of the gene neighbours experiment (run.sh, README.md).

  design.py heldout GENOMES                  species left out of the database (8, of genera with congeners)
  design.py community GENOMES HELDOUT        accession, relative abundance of the 30 genomes sampled
  design.py truth COMMUNITY GENOMES TAXONOMY the species present, as GTDB lineages (protal --profile_truth)
  design.py report WORK                      the tables of README.md from the runs in WORK
"""

import collections
import csv
import glob
import math
import os
import random
import re
import statistics
import subprocess
import sys


def genomes(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def species_of(row):
    return row["gtdb_taxonomy"].split(";")[-1]


def heldout(path):
    rows = genomes(path)
    by_genus = collections.defaultdict(set)
    for r in rows:
        by_genus[r["gtdb_taxonomy"].split(";")[-2]].add(species_of(r))
    pool = sorted(s for genus in by_genus.values() if len(genus) >= 2 for s in genus)
    for s in sorted(random.Random(4).sample(pool, 8)):
        print(s)


def community(path, heldout_path):
    rows = genomes(path)
    with open(heldout_path) as fh:
        held = {line.strip() for line in fh if line.strip()}
    rng = random.Random(6)
    in_db = sorted({species_of(r) for r in rows} - held)
    chosen = rng.sample(in_db, 22)
    print("accession\trelative_abundance")
    for i, sp in enumerate(sorted(chosen) + sorted(held)):
        rep = i % 2 == 0 or sp in held  # half of the database's species from their other genome (a strain)
        acc = next(r["accession"] for r in rows if species_of(r) == sp and (r["gtdb_representative"] == "t") == rep)
        print(f"{acc}\t{math.exp(rng.gauss(0, 1)):.4f}")


def truth(community_path, genomes_path, _taxonomy):
    lineage = {r["accession"]: r["gtdb_taxonomy"] for r in genomes(genomes_path)}
    with open(community_path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            print(lineage[row["accession"]])


def sam_lines(path):
    text = subprocess.run(["zstd", "-dcq", path], check=True, stdout=subprocess.PIPE, text=True).stdout \
        if path.endswith(".zst") else open(path).read()
    return [l.split("\t") for l in text.splitlines() if l and not l.startswith("@")]


def auc(positives, negatives):
    if not positives or not negatives:
        return float("nan")
    wins = sum((p > n) + 0.5 * (p == n) for p in positives for n in negatives)
    return wins / (len(positives) * len(negatives))


def report(work):
    os.chdir(work)
    with open("heldout.txt") as fh:
        held = {line.strip() for line in fh if line.strip()}
    names = sorted({os.path.basename(p)[4:-4] for p in glob.glob("log_*.txt")})
    times = collections.defaultdict(list)
    for name in names:
        wall, user, rss = open(f"time_{name}.txt").read().split()[-3:]
        times[name.rsplit("_", 1)[0]].append((float(wall), float(user), int(rss)))
    print("| read type | gene neighbours | wall s (median of 3) | user s | max RSS MB | runs |")
    print("|---|---|---|---|---|---|")
    for key in sorted(times):
        t = times[key]
        print(f"| {key.split('_')[0]} | {key.split('_')[1]} | {statistics.median(x[0] for x in t):.1f} | "
              f"{statistics.median(x[1] for x in t):.1f} | {max(x[2] for x in t) / 1024:.0f} | "
              f"{', '.join(f'{x[0]:.1f}' for x in t)} |")
    print()
    for read_type in ("pe", "pb", "ont"):
        print(f"### {read_type}\n")
        print("| | with | without |\n|---|---|---|")
        cells = collections.defaultdict(dict)
        for arm in ("with", "without"):
            name = f"{read_type}_{arm}_1"
            log = open(f"log_{name}.txt").read()
            for pattern, label in ((r"Gene neighbours: (\d+) fragments paired across", "fragments paired across two genes"),
                                   (r"; (\d+) guided mates found on the gene next", "guided mates found on the next gene"),
                                   (r"(\d+) fragments had one mate sure", "fragments with a guiding mate"),
                                   (r"Gene neighbours: (\d+) genes found on reads", "genes found where the neighbours put them"),
                                   (r"(\d+) gene hits fit several taxa", "gene hits that fit several taxa"),
                                   (r"and (\d+) had no hit of it", "gene hits without one of the read's taxon")):
                m = re.search(pattern, log)
                if m:
                    cells[label][arm] = m.group(1)
            sam = (glob.glob(f"out/{name}/{read_type}*.sam.zst") + glob.glob(f"out/{name}/{read_type}*.sam"))[0]
            records = sam_lines(sam)
            primary = [r for r in records if not int(r[1]) & 0x100]
            cells["primary records"][arm] = len(primary)
            cells["aligned bases (primary, M/X/=)"][arm] = sum(
                int(n) for r in primary for n, op in re.findall(r"(\d+)([MX=])", r[5]))
            if read_type == "pe":
                both = sum(1 for r in primary if int(r[1]) & 0x40 and int(r[1]) & 0x2)
                cells["fragments written as a proper pair"][arm] = both
                cells["of them on two genes"][arm] = sum(1 for r in primary if int(r[1]) & 0x40 and int(r[1]) & 0x2 and r[6] not in ("=", "*"))
            else:
                per_read = collections.Counter(r[0] for r in primary)
                cells["genes per read with a record"][arm] = f"{statistics.mean(per_read.values()):.2f}"
            with open(glob.glob(f"out/{name}/{read_type}.profile.truth_annotated")[0]) as fh:
                header = fh.readline().rstrip("\n").split("\t")
                rows = [dict(zip(header, l.rstrip("\n").split("\t"))) for l in fh]
            called = [r for r in rows if r.get("prediction") in ("1", "1.0")]
            tp = [r for r in called if r["truth"] == "1"]
            fp = [r for r in called if r["truth"] != "1"]
            cells["species called: true / false"][arm] = f"{len(tp)} / {len(fp)}"
            present = [r for r in rows if r["truth"] == "1"]
            absent = [r for r in rows if r["truth"] != "1"]
            for feature in ("adjacent_expected_share", "adjacent_unlikely_share"):
                p = [float(r[feature]) for r in present]
                a = [float(r[feature]) for r in absent]
                cells[f"{feature}: present / absent taxa (mean)"][arm] = \
                    f"{statistics.mean(p) if p else 0:.3f} / {statistics.mean(a) if a else 0:.3f}"
                cells[f"{feature}: AUC present vs absent"][arm] = f"{auc(p, a):.3f}"
            cells["taxa with a record: present / absent"][arm] = f"{len(present)} / {len(absent)}"
        for label, values in cells.items():
            print(f"| {label} | {values.get('with', '')} | {values.get('without', '')} |")
        print()


if __name__ == "__main__":
    command, args = sys.argv[1], sys.argv[2:]
    {"heldout": heldout, "community": community, "truth": truth, "report": report}[command](*args)
