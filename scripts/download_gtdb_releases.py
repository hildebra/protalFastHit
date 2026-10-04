#!/usr/bin/env python3
"""Download the inputs of protal databases for several GTDB releases, on a node with internet.

Usage: python3 scripts/download_gtdb_releases.py -o INPUTS --releases 220,226 [-t 8] [download_gtdb.py options]

The download phase of a release-by-release build (the build phase, which needs no internet, is
scripts/build_gtdb_releases.py --inputs INPUTS). For each release it runs scripts/download_gtdb.py
into INPUTS/gtdb_r<release> (GTDB's files, the genomes to simulate from, the species pool), with
the options it does not know itself passed on (--species, --per_species, --mirror, ...). A rerun
downloads only what is missing. At the end it writes INPUTS/download_summary.tsv and prints it:
per release, the version, what the taxonomy holds (genomes, species, genera, families, orders,
classes, phyla, per domain), the genomes delivered and missing, the strains and species to simulate
from, the marker files, the bytes on disk and the time this run took on it. A release that fails
is reported and the others go on; the exit code is then 1.
"""
import argparse
import collections
import glob
import gzip
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
DOWNLOAD = os.path.join(HERE, "download_gtdb.py")
RANKS = ("phylum", "class", "order", "family", "genus", "species")
COLUMNS = ("release", "version", "folder", "status", "genomes_in_taxonomy", "species", "genera", "families", "orders",
           "classes", "phyla", "bacteria_species", "archaea_species", "genomes_delivered", "genomes_missing",
           "simulation_species", "strains", "species_with_strains", "marker_files", "size_gb", "download_time", "log")


def say(message):
    print(f"[{time.strftime('%H:%M:%S')}] {message}", flush=True)


def clock(seconds):
    seconds = int(round(seconds))
    return f"{seconds // 3600}:{seconds // 60 % 60:02d}:{seconds % 60:02d}"


def folder_size(path):
    total = 0
    for root, _dirs, names in os.walk(path):
        for name in names:
            try:
                total += os.path.getsize(os.path.join(root, name))
            except OSError:
                pass
    return total


def taxonomy_counts(release_dir, number):
    """The clades of a release's taxonomy files: {"genomes": n, "species": n, ..., "bacteria_species": n,
    "archaea_species": n}; zeros without the files."""
    counts = collections.Counter()
    clades = collections.defaultdict(set)
    for mset in ("bac120", "ar53"):
        paths = glob.glob(os.path.join(release_dir, f"{mset}_taxonomy_r{number}.tsv*"))
        if not paths:
            continue
        opener = gzip.open if paths[0].endswith(".gz") else open
        with opener(paths[0], "rt") as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 2 or not f[1].startswith("d__"):
                    continue
                counts["genomes"] += 1
                names = f[1].split(";")
                for depth, rank in enumerate(RANKS, 1):
                    if len(names) > depth:
                        clades[rank].add(";".join(names[:depth + 1]))
                if len(names) == 7:
                    clades["species_" + names[0].removeprefix("d__").lower()].add(f[1])
    for rank in RANKS:
        counts[rank] = len(clades[rank])
    counts["bacteria_species"] = len(clades["species_bacteria"])
    counts["archaea_species"] = len(clades["species_archaea"])
    return counts


def summarize(folder, status, seconds, log):
    """The summary row of a release's inputs folder (download.json, the taxonomy, the genomes)."""
    row = {c: "" for c in COLUMNS}
    row.update(folder=folder, status=status, download_time=clock(seconds), log=log)
    try:
        with open(os.path.join(folder, "download.json")) as fh:
            state = json.load(fh)
    except (OSError, ValueError):
        state = {}
    release = state.get("release", {})
    number = release.get("number", "")
    row["release"], row["version"] = number, release.get("version", "")
    if number:
        counts = taxonomy_counts(os.path.join(folder, "release"), number)
        row.update(genomes_in_taxonomy=counts["genomes"], species=counts["species"], genera=counts["genus"],
                   families=counts["family"], orders=counts["order"], classes=counts["class"], phyla=counts["phylum"],
                   bacteria_species=counts["bacteria_species"], archaea_species=counts["archaea_species"])
        row["marker_files"] = sum(len(glob.glob(os.path.join(folder, "release", f"genomic_files_{kind}", "**", "*.fna*"),
                                                recursive=True)) for kind in ("reps", "all"))
    genomes = state.get("genomes", {})
    if genomes:
        row.update(genomes_delivered=genomes.get("delivered", ""), genomes_missing=genomes.get("missing", ""),
                   simulation_species=genomes.get("simulation_species", ""), strains=genomes.get("strains", ""),
                   species_with_strains=genomes.get("species_with_strains", ""))
    row["size_gb"] = f"{folder_size(folder) / 1e9:.1f}"
    return row


def write_summary(path, rows):
    with open(path + ".partial", "w", newline="\n") as fh:
        fh.write("\t".join(COLUMNS) + "\n")
        for row in rows:
            fh.write("\t".join(str(row[c]) for c in COLUMNS) + "\n")
    os.replace(path + ".partial", path)


def text_table(rows, columns):
    widths = [max(len(c), *(len(str(r[c])) for r in rows)) for c in columns]
    lines = ["  ".join(c.ljust(w) for c, w in zip(columns, widths)).rstrip()]
    lines += ["  ".join(str(r[c]).ljust(w) for c, w in zip(columns, widths)).rstrip() for r in rows]
    return "\n".join(lines)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-o", "--out", required=True, help="root folder: INPUTS/gtdb_r<release> per release")
    p.add_argument("--releases", required=True, help="comma-separated GTDB releases (207 or later; 214.1 for a point release)")
    p.add_argument("-t", "--threads", type=int, default=8, help="passed to download_gtdb.py (default 8)")
    args, passthrough = p.parse_known_args()
    releases = [r.strip().removeprefix("r") for r in args.releases.split(",") if r.strip()]
    if not releases:
        p.error("--releases: none given")
    os.makedirs(args.out, exist_ok=True)
    rows, failed = [], []
    for release in releases:
        folder = os.path.join(args.out, f"gtdb_r{release}")
        log = os.path.join(args.out, f"download_r{release}.log")
        command = [sys.executable, DOWNLOAD, "-o", folder, "--release", release, "-t", str(args.threads), *passthrough]
        say(f"GTDB r{release} into {folder} ({os.path.basename(log)})")
        began = time.time()
        with open(log, "a") as fh:
            fh.write(f"== {time.strftime('%Y-%m-%d %H:%M:%S')} {' '.join(command)}\n")
            fh.flush()
            result = subprocess.run(command, stdout=fh, stderr=subprocess.STDOUT)
        status = "ok" if result.returncode == 0 else f"failed ({result.returncode})"
        if result.returncode:
            failed.append(release)
            with open(log, errors="replace") as fh:
                tail = fh.read().splitlines()[-8:]
            say(f"    download_gtdb.py failed ({result.returncode}); the end of {log}:\n" + "\n".join("    " + t for t in tail))
        row = summarize(folder, status, time.time() - began, log)
        rows.append(row)
        say(f"    {status} in {row['download_time']}: {row['species']} species of {row['genomes_in_taxonomy']} genomes in the "
            f"taxonomy, {row['genomes_delivered'] or 0} genomes delivered ({row['genomes_missing'] or 0} missing), "
            f"{row['size_gb']} GB")
    summary = os.path.join(args.out, "download_summary.tsv")
    write_summary(summary, rows)
    print("\n" + text_table(rows, [c for c in COLUMNS if c not in ("folder", "log")]) + f"\n\nSummary: {summary}", flush=True)
    if failed:
        say(f"failed: r{', r'.join(failed)}")
        sys.exit(1)
    say(f"Inputs of {len(rows)} release{'s' if len(rows) > 1 else ''} in {args.out} (scripts/build_gtdb_releases.py --inputs {args.out})")


if __name__ == "__main__":
    main()
