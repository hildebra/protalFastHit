#!/usr/bin/env python3
"""Build and train protal databases, full and reduced, for several GTDB releases; no internet needed.

Usage: python3 scripts/build_gtdb_releases.py --inputs INPUTS --outdir OUT --protal build/protal \\
           --simulator build/simulate_metagenomes [-t 16] [--releases 220,226] [--variants full,n12]
           [--n-genes 12] [--scratch DIR] [--rerun] [build_gtdb_database.py options]

INPUTS is the folder of scripts/download_gtdb_releases.py (INPUTS/gtdb_r<release> per release; the
download phase, on a node with internet). For each release (--releases, default every release in
INPUTS) and each variant, in this order, the script runs scripts/build_gtdb_database.py into
OUT/r<release>_<variant>, with the options it does not know itself passed on (--samples,
--read-types, --features, ...):
  full   every marker gene, with --rank-genes: the genes ranked from its training database
         (OUT/r<release>_full/gene_ranking.tsv)
  n12    the --n-genes most distinctive genes (--genes-per-domain of them the best of each
         domain, so that archaea are covered), with --gene-ranking from the full variant's
         ranking when it is there, so that no ranking database is built
With --scratch, each database gets SCRATCH/r<release>_<variant>. A database whose
database.protal and model_logs/summary.txt are there is kept (--rerun builds it again; a
build_gtdb_database.py rerun itself skips what it completed). A failed database is reported and
the others go on; the exit code is then 1.

At the end OUT/build_summary.tsv (one line per database) and OUT/build_summary.txt are written
and printed: the release and protal version, the marker genes, the taxa the database holds
(species, genera, families, orders, classes, phyla; species per domain), the genomes simulated
from and the species held out, the size of database.protal, the time and peak memory of its
build and of the profiling of the training samples (the protal run of collect_training_data.py:
peak memory of the collector and protal together), the run's time, and for each read type's model
(pe, se, pb, ont) its F1, FP per sample, sensitivity and precision on the independent test set and
its F1 and FP per sample with species held out in training (model_logs/summary.txt).
"""
import argparse
import collections
import glob
import json
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
BUILD = os.path.join(HERE, "build_gtdb_database.py")
READ_TYPES = ("pe", "se", "pb", "ont")
RANKS = ("domain", "phylum", "class", "order", "family", "genus", "species")
COLUMNS = ["release", "variant", "status", "protal_version", "marker_genes", "genes", "species", "genera", "families",
           "orders", "classes", "phyla", "bacteria_species", "archaea_species", "genomes_simulated", "species_held_out",
           "database_gb", "build_time", "build_peak_gb", "profiling_peak_gb", "run_time"]
for _t in READ_TYPES:
    COLUMNS += [f"{_t}_test_F1", f"{_t}_test_FP_per_sample", f"{_t}_test_sensitivity", f"{_t}_test_precision",
                f"{_t}_heldout_F1", f"{_t}_heldout_FP_per_sample"]
COLUMNS += ["outdir", "log"]


def say(message):
    print(f"[{time.strftime('%H:%M:%S')}] {message}", flush=True)


def clock(seconds):
    seconds = int(round(seconds))
    return f"{seconds // 3600}:{seconds // 60 % 60:02d}:{seconds % 60:02d}"


def gigabytes(text):
    """'3.5 GB' or '512 MB' -> '3.5' (GB)."""
    m = re.match(r"([\d.]+) ([MG])B", text or "")
    if not m:
        return ""
    value = float(m.group(1)) / (1000 if m.group(2) == "M" else 1)
    return f"{value:.1f}" if value >= 0.1 else f"{value:.2f}"


def read_metadata(path):
    try:
        with open(path) as fh:
            return dict(line.rstrip("\n").split("\t", 1) for line in fh if "\t" in line)
    except OSError:
        return {}


def taxonomy_counts(path):
    """internal_taxonomy.dmp -> {rank: count} and the species per domain ("species_bacteria", ...)."""
    counts = collections.Counter()
    nodes = {}
    try:
        with open(path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 5 and f[0].isdigit():
                    nodes[int(f[0])] = (int(f[1]), f[3], f[4])
    except OSError:
        return counts
    for taxid, (parent, name, rank) in nodes.items():
        counts[rank] += 1
        if rank == "species":
            node, seen = taxid, set()
            while node in nodes and node not in seen:
                seen.add(node)
                up, up_name, up_rank = nodes[node]
                if up_rank == "domain":
                    counts["species_" + up_name.removeprefix("d__").lower()] += 1
                    break
                node = up
    return counts


def model_table(path):
    """model_logs/summary.txt -> {(read type, "test" | "heldout"): {column: value}}."""
    table = {}
    try:
        with open(path) as fh:
            lines = fh.read().splitlines()
    except OSError:
        return table
    header = None
    for line in lines:
        f = re.split(r"  +", line.rstrip())
        if f[:2] == ["read type", "evaluated on"]:
            header = f
        elif header and len(f) == len(header) and f[0] in READ_TYPES:
            which = "test" if f[1].startswith("independent") else "heldout"
            table[(f[0], which)] = dict(zip(header, f))
    return table


def run_numbers(log):
    """From a run's console log: the build's time and peak memory, the profiling's peak memory."""
    numbers = {"build_time": "", "build_peak_gb": "", "profiling_peak_gb": ""}
    try:
        with open(log, errors="replace") as fh:
            text = fh.read()
    except OSError:
        return numbers
    m = re.search(r"built protal_db(?: in the background)? in (\d+:\d\d:\d\d)(?:, peak memory ([\d.]+ [MG]B))?", text)
    if m:
        numbers["build_time"], numbers["build_peak_gb"] = m.group(1), gigabytes(m.group(2))
    m = re.search(r"collected in \d+:\d\d:\d\d, peak memory ([\d.]+ [MG]B)", text)
    if m:
        numbers["profiling_peak_gb"] = gigabytes(m.group(1))
    return numbers


def summarize(release, variant, outdir, status, seconds, log):
    row = {c: "" for c in COLUMNS}
    row.update(release=release, variant=variant, status=status, run_time=clock(seconds), outdir=outdir, log=log)
    db = os.path.join(outdir, "protal_db")
    metadata = read_metadata(os.path.join(db, "build_metadata.tsv"))
    row["protal_version"] = metadata.get("protal_version", "")
    genes = metadata.get("marker_genes", "")
    row["marker_genes"] = genes
    m = re.match(r"(\d+) of (\d+)", genes)
    row["genes"] = m.group(1) if m else ("all" if genes == "all" else "")
    row["genomes_simulated"] = metadata.get("genome_table", "")
    row["species_held_out"] = metadata.get("classifier_training_species_left_out", "")
    counts = taxonomy_counts(os.path.join(outdir, "internal_taxonomy.dmp"))
    if counts:
        row.update(species=counts["species"], genera=counts["genus"], families=counts["family"], orders=counts["order"],
                   classes=counts["class"], phyla=counts["phylum"], bacteria_species=counts["species_bacteria"],
                   archaea_species=counts["species_archaea"])
    if os.path.isfile(os.path.join(db, "database.protal")):
        row["database_gb"] = f"{os.path.getsize(os.path.join(db, 'database.protal')) / 1e9:.2f}"
    row.update(run_numbers(log))
    models = model_table(os.path.join(outdir, "model_logs", "summary.txt"))
    for t in READ_TYPES:
        test, held = models.get((t, "test"), {}), models.get((t, "heldout"), {})
        row[f"{t}_test_F1"] = test.get("F1", "")
        row[f"{t}_test_FP_per_sample"] = test.get("FP/sample", "")
        row[f"{t}_test_sensitivity"] = test.get("sensitivity", "")
        row[f"{t}_test_precision"] = test.get("precision", "")
        row[f"{t}_heldout_F1"] = held.get("F1", "")
        row[f"{t}_heldout_FP_per_sample"] = held.get("FP/sample", "")
    return row


def describe(row):
    """A readable block of a database's summary."""
    lines = [f"GTDB r{row['release']}, {row['variant']}: {row['status']}; protal {row['protal_version']}; {row['outdir']}",
             f"  marker genes: {row['marker_genes'] or '?'}",
             f"  taxa: {row['species']} species ({row['bacteria_species']} bacteria, {row['archaea_species']} archaea), "
             f"{row['genera']} genera, {row['families']} families, {row['orders']} orders, {row['classes']} classes, "
             f"{row['phyla']} phyla",
             f"  simulated from: {row['genomes_simulated'] or '?'}; species held out of the training database: "
             f"{row['species_held_out'] or '?'}",
             f"  database.protal: {row['database_gb'] or '?'} GB; its build {row['build_time'] or '?'}, peak memory "
             f"{row['build_peak_gb'] or '?'} GB; profiling the training samples: peak memory {row['profiling_peak_gb'] or '?'} GB; "
             f"the run {row['run_time']}"]
    for t in READ_TYPES:
        if row[f"{t}_test_F1"] or row[f"{t}_heldout_F1"]:
            lines.append(f"  {t}: independent test F1 {row[f'{t}_test_F1'] or '-'}, FP/sample {row[f'{t}_test_FP_per_sample'] or '-'}, "
                         f"sensitivity {row[f'{t}_test_sensitivity'] or '-'}, precision {row[f'{t}_test_precision'] or '-'}; "
                         f"species held out F1 {row[f'{t}_heldout_F1'] or '-'}, FP/sample {row[f'{t}_heldout_FP_per_sample'] or '-'}")
    return "\n".join(lines)


def write_summary(outdir, rows):
    tsv, txt = os.path.join(outdir, "build_summary.tsv"), os.path.join(outdir, "build_summary.txt")
    with open(tsv + ".partial", "w", newline="\n") as fh:
        fh.write("\t".join(COLUMNS) + "\n")
        for row in rows:
            fh.write("\t".join(str(row[c]) for c in COLUMNS) + "\n")
    os.replace(tsv + ".partial", tsv)
    text = "\n\n".join(describe(r) for r in rows) + "\n"
    with open(txt, "w", newline="\n") as fh:
        fh.write(text)
    return tsv, text


def releases_in(inputs):
    found = []
    for folder in sorted(glob.glob(os.path.join(inputs, "gtdb_r*"))):
        state = os.path.join(folder, "download.json")
        if os.path.isfile(state):
            found.append(os.path.basename(folder)[len("gtdb_r"):])
    return found


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--inputs", required=True, help="folder of download_gtdb_releases.py (INPUTS/gtdb_r<release>)")
    p.add_argument("--outdir", required=True, help="root folder: OUT/r<release>_<variant> per database")
    p.add_argument("--releases", help="comma-separated releases (default: every gtdb_r* in INPUTS with a download.json)")
    p.add_argument("--variants", default="full,n12", help="comma-separated: full, and n<N> for a reduced database "
                                                          "(default full,n12)")
    p.add_argument("--n-genes", type=int, default=12, help="the genes of the reduced variant (default 12)")
    p.add_argument("--genes-per-domain", type=int, help="of them, the best of each domain (default a third)")
    p.add_argument("--protal", default="protal", help="protal executable")
    p.add_argument("--simulator", default="simulate_metagenomes", help="simulate_metagenomes executable")
    p.add_argument("-t", "--threads", type=int, default=8)
    p.add_argument("--scratch", help="a fast local disk: SCRATCH/r<release>_<variant> for each database's samples")
    p.add_argument("--rerun", action="store_true", help="run build_gtdb_database.py again for databases already built")
    args, passthrough = p.parse_known_args()
    releases = [r.strip().removeprefix("r") for r in args.releases.split(",") if r.strip()] if args.releases \
        else releases_in(args.inputs)
    if not releases:
        p.error(f"no release: --releases, or gtdb_r* folders with a download.json in {args.inputs}")
    variants = [v.strip() for v in args.variants.split(",") if v.strip()]
    for v in variants:
        if v != "full" and not re.fullmatch(r"n\d+", v):
            p.error(f"--variants: {v} is neither full nor n<N>")
    variants.sort(key=lambda v: (v != "full", v))  # the full database first: its ranking serves the reduced one
    os.makedirs(args.outdir, exist_ok=True)
    rows, failed = [], []
    for release in releases:
        inputs = os.path.join(args.inputs, f"gtdb_r{release}")
        if not os.path.isfile(os.path.join(inputs, "download.json")):
            say(f"GTDB r{release}: {inputs} has no download.json (download_gtdb_releases.py first); skipped")
            failed.append(f"r{release}")
            continue
        for variant in variants:
            outdir = os.path.join(args.outdir, f"r{release}_{variant}")
            log = os.path.join(args.outdir, f"r{release}_{variant}.log")
            done = os.path.isfile(os.path.join(outdir, "protal_db", "database.protal")) and \
                os.path.isfile(os.path.join(outdir, "model_logs", "summary.txt"))
            command = [sys.executable, BUILD, "--inputs", inputs, "--outdir", outdir, "--protal", args.protal,
                       "--simulator", args.simulator, "-t", str(args.threads)]
            if args.scratch:
                command += ["--scratch", os.path.join(args.scratch, f"r{release}_{variant}")]
            if variant == "full":
                command += ["--rank-genes"]
            else:
                n = int(variant[1:]) if variant != f"n{args.n_genes}" else args.n_genes
                command += ["--n-genes", str(n)]
                if args.genes_per_domain is not None:
                    command += ["--genes-per-domain", str(args.genes_per_domain)]
                ranking = os.path.join(args.outdir, f"r{release}_full", "gene_ranking.tsv")
                if os.path.isfile(ranking):
                    command += ["--gene-ranking", ranking]
            command += passthrough
            began = time.time()
            if done and not args.rerun:
                say(f"GTDB r{release}, {variant}: {outdir} is built and trained; kept (--rerun builds it again)")
                status = "kept"
            else:
                say(f"GTDB r{release}, {variant}: building {outdir} ({os.path.basename(log)})")
                with open(log, "a") as fh:
                    fh.write(f"== {time.strftime('%Y-%m-%d %H:%M:%S')} {' '.join(command)}\n")
                    fh.flush()
                    result = subprocess.run(command, stdout=fh, stderr=subprocess.STDOUT)
                status = "ok" if result.returncode == 0 else f"failed ({result.returncode})"
                if result.returncode:
                    failed.append(f"r{release} {variant}")
                    with open(log, errors="replace") as fh:
                        tail = fh.read().splitlines()[-8:]
                    say(f"    build_gtdb_database.py failed ({result.returncode}); the end of {log}:\n" +
                        "\n".join("    " + t for t in tail))
            row = summarize(release, variant, outdir, status, time.time() - began, log)
            rows.append(row)
            say("    " + describe(row).replace("\n", "\n    "))
    tsv, text = write_summary(args.outdir, rows)
    print("\n" + text + f"\nSummary: {tsv}", flush=True)
    if failed:
        say(f"failed: {'; '.join(failed)}")
        sys.exit(1)


if __name__ == "__main__":
    main()
