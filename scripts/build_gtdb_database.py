#!/usr/bin/env python3
"""Build and train a complete protal database from an extracted GTDB release.

Usage: python3 scripts/build_gtdb_database.py --inputs INPUTS_DIR --outdir OUT_DIR
       python3 scripts/build_gtdb_database.py --gtdb GTDB_DIR --outdir OUT_DIR

--inputs is a folder of scripts/download_gtdb.py (run on a node with internet): the
release, the genomes to simulate from and the species pool, in one place that
serves every build of that release. It stands for --gtdb INPUTS/release
--extra-genomes INPUTS/genomes --simulate-species INPUTS/simulation_species.txt.
With --gtdb, the release must include taxonomy, marker gene FASTAs and extracted
whole genome FASTAs; a genome table can be supplied explicitly when the release
uses an unusual layout.

The presence model is trained on metagenomes simulated from these genomes, made
harder in two ways than simulating the database's own references:
- Strains. GTDB distributes whole genomes of the representatives only
  (genomic_files_reps), and with only those every simulated species is the
  database's own reference, closer to it than real strains are.
  scripts/download_gtdb.py picks other genomes of the species from GTDB's
  metadata and downloads them from NCBI, so that species are simulated from
  other strains, too. OUT_DIR/model_logs/genome_table.txt says how often.
- Species the database lacks. The samples are profiled against a training
  database (OUT_DIR/work/training_db; SCRATCH/training_db with --scratch) that
  leaves whole clades of every rank (--holdout-clades) and --holdout of the
  other species out
  (model_logs/heldout_species.txt): their reads land on relatives, as those of
  organisms GTDB lacks do in real samples, and the model learns to reject those
  relatives. The report gives false positive and false negative rates by rank.
  The finished database has all species and the model trained so. The training
  database costs a second index build and its disk space: it is built first,
  alone on the node, as the profiling waits for it; the finished database after
  it, at the idle scheduling class, on the cores the profiling leaves.

One model per read type (--read-types, default pe,se,pb,ont), trained in
parallel: paired-end reads, their first reads alone (single-end), and PacBio
and Nanopore reads of the same communities (all made by simulate_metagenomes). Besides the training
data, an independent test set of another design (--test-*: other depths,
community sizes, abundances and strain mixes) is profiled and scored by each
model: cross-validation on the training data cannot show what its design lacks.
Both collections simulate in the background once the training database is built
(they take minutes since simulate_metagenomes makes every read type itself, far
less than the profiling), at a lower priority than the profiling: as the
simulations go on, in protal runs of --profile-blocks GB of reads or more (at
most --profile-block-max), each design point's reads removed once profiled, so
that the samples need not all be on the disk at once (--profile-blocks 0: both
collections in one protal run once all are simulated). The largest design points
can be streamed into protal through named pipes instead of written
(--stream-above): by default only as many as the room on the samples' disk
requires, which the run estimates before its builds (the genome store, the
databases, the samples' SAMs and profiles, the reads) and says, beside the most it
took, at its end.
The design reaches the depths of real samples (2M and 10M read pairs, 1.5 and
6 Gb of long reads, a few samples each: DEPTH:SAMPLES), and each model gets a
knob curve over the sample's depth (--depth-knob-read-types): a deep sample
holds many more absent taxa with a few reads, and needs a higher threshold.

Four scenarios of real studies (--scenarios, scenarios.py: gut, soil, shallow
soil, 90% host reads) add samples of their own, each at a depth drawn around its
scenario's: their hold-in samples join the training data, their hold-out samples
the test set, and every model's report scores both. The models are gradient-boosted
trees on the default feature set, evaluated with rows, samples and species held
out (--features, --evaluation basic; the clades too with --evaluation full); --features auto lets each trainer choose its
set with species held out, and the run then says which set won and why.

OUT_DIR holds, each file once:
  protal_db/      the database: database.protal and build_metadata.tsv (what it was built from and with), and the
                  converter's gene2geneid.tsv and genome2tiid.tsv (gene ids and genomes by name)
  model_logs/     what tells whether the models are good, and what the run chose
  logs/           every step's log
  console.log     the console's lines
  work/           what a rerun reuses: the stage keys, the genome tables, the taxonomy, the training and test tables
                  (to retrain); without --scratch also the samples, the training database, the genome store and the
                  in-silico strains. Not needed to use the database
model_logs/ has summary.txt (TP, FP, TN, FN, sensitivity, specificity, precision and F1 of each model), each read type's
model (trained_model*.xml, as in the database) and training report (how it does on species and clades it was not trained
on, on the independent test set, false positive and false negative rates by rank, against the previous model and training
procedure), its numbers as JSON, the per-taxon predictions and calls, the threshold table, the parity check with
protal (parity.txt), the genome table summary, the species held out (heldout_species.txt, holdout.txt) and the clouds
that steered them; and what the model's conservation features rest on, on real genomes: gene_congeners.tsv (protal
--build: how each gene differs between congeners against within species), gene_incongruence.tsv (protal --build: every
near pair of gene copies across genera, and which copy is suspect, contamination or a transfer; the suspect ones go into
the database as suspect_copies.tsv and a run leaves their records out) and relatives_by_gene_conservation.txt
(trace_relatives.py: where the reads of the held-out species land, by the genes' factors). model_logs/error_reads/ tells
what each model's errors rest on in every sample of the training data and the test set (--error-reads, default all): per
read type, a table of every sample's error taxa and where their reads went, the non-hits among them (reads that seeded
on taxa but aligned nowhere, whose unmapped records protal writes for these samples; error_reads.py). With --share-logs
it also keeps the SAM records of those reads, the false positives' and the false negatives' in files of their own per
sample (a sample of each taxon's), and the run ends by packing OUT_DIR/<name>_share.tar.gz: console.log, logs/,
model_logs/ (without the models), the taxonomy and the training and test tables, and from the samples' disk the design
points' simulator logs and the in-silico strains' table, to copy off the cluster.

A reduced database holds a subset of the marker genes (--n-genes N: the N most
distinctive by prevalence x unique k-mer share, ranked by scripts/rank_genes.py
from a full build of the training database, or from --gene-ranking; --genes: the
genes named). The release is then converted whole into OUTDIR/work/converted and both
database folders are derived from it with the subset (their gene neighbours
counted over it); model_logs/gene_ranking.tsv and gene_subset.txt record the choice.
Each domain keeps its best genes in the subset (--genes-per-domain), so that
archaea are covered too. A run with every gene writes the ranking with
--rank-genes, for the reduced run to take (--gene-ranking); scripts/
build_gtdb_releases.py chains both for several releases.

The tools the run needs are checked before it starts, and protal and the simulator
must be of the source these scripts are at: of its version and (as their --version
says since 0.7.3) built from its commit, with nothing changed since in src/, lib/
or the build files (--no-binary-check runs them anyway). A run that stops (a failure,
SIGTERM, Ctrl-C) stops every command it started, and one that fails in the
background (the finished database's build) stops the run within seconds. A rerun
into the same OUTDIR resumes: the conversion and the two index builds are skipped
when their inputs are those of the run that completed them (OUTDIR/work/stages), and
the collector reuses the samples it simulated and profiled with the same database,
protal and design. An OUTDIR of the layout before 2026-10-08 (the logs, .stages and
the tables beside protal_db) is not resumed: its stages are done again.

Each stage writes to a log of its own in OUTDIR/logs. On the console, each line has the
time and how long the run has taken. A step says when it starts ("4/8 training
data (its log): ...") and, on an indented line, when it ends: how long it took, its
peak memory and a few numbers of what it made. With --progress-every, each stage
running also says every so many seconds how long it has run, the memory it takes
and the last line of its log.
"""
import argparse
import collections
import concurrent.futures
import contextlib
import csv
import fcntl
import glob
import gzip
import hashlib
import io
import json
import os
import random
import re
import shutil
import signal
import subprocess
import sys
import tarfile
import threading
import time

HERE = os.path.dirname(os.path.abspath(__file__))
CONVERTER = os.path.join(HERE, "mini_db", "gtdb_to_protal_db.py")
GENE_NEIGHBOURS = os.path.join(HERE, "mini_db", "gene_neighbours.py")
FOREIGN_RATES = os.path.join(HERE, "foreign_rates.py")
TRAINER = os.path.join(HERE, "machine_learning_cmdline.py")
COLLECTOR = os.path.join(HERE, "collect_training_data.py")
PARITY = os.path.join(HERE, "check_model_parity.py")
TRACE = os.path.join(HERE, "trace_relatives.py")
ERROR_READS = os.path.join(HERE, "error_reads.py")
ANCESTRY = os.path.join(HERE, "ancestry_sites.py")
RANKER = os.path.join(HERE, "rank_genes.py")
INSILICO = os.path.join(HERE, "insilico_strains.py")
# --insilico-ani when the conversion left no gene_positions.tsv (--no-gene-neighbours): no real strains to draw from.
INSILICO_FALLBACK_ANI = "97-99.5"
SOURCE = os.path.dirname(HERE)  # the checkout these scripts are part of
# OUTDIR's folders (the module's docstring): the database, the evaluation, every step's log, what a rerun reuses.
DATABASE, REPORTS, LOGS, WORK = "protal_db", "model_logs", "logs", "work"
# What protal and the simulator are built from (protal_commit.cmake marks a build of uncommitted changes to them).
BUILD_SOURCES = ("src", "lib", "CMakeLists.txt", "protal_config.h.in", "protal_commit.cmake")
ACCESSION = re.compile(r"(?:RS_|GB_)?(GC[AF]_\d{9}\.\d+)")
GENOME_CACHE = "genome_cache.tsv.gz"  # --genome-cache auto: in the --inputs folder (else OUTDIR/work), across builds
CONTIG_CACHE = "genome_contigs.tsv.gz"  # the reports' contig names (trace_relatives.genome_contigs), on the samples' disk
ANCESTRY_REFERENCE = "ancestry_reference.fna"  # the training database's reference.fna, a hard link kept for the report
PROTAL_CPU = "protal_cpu.tsv"  # collect_training_data.py's CPU table of its protal runs, in its -o folder
sys.path.insert(0, os.path.join(HERE, "mini_db"))
sys.path.insert(0, HERE)
import lineages  # noqa: E402
from gtdb_to_protal_db import (allele_genome, clear_build_outputs, full_reference_path, marker_files, normalize_accession,  # noqa: E402
                               read_gene_ids, read_gene_list, read_representatives,
                               REFERENCE_WRITTEN, remove_full_reference as remove_full_reference_files)
import rank_genes  # noqa: E402
import scenarios  # noqa: E402
from model_pmml import MODEL_FILES, write_placeholder  # noqa: E402
from model_features import DEFAULT_FEATURE_SET, FEATURE_SETS, feature_set_name  # noqa: E402
from collect_training_data import INSILICO_PREFIX, TABLES, clock, congener_spec, last_line, manifest_rows, read_clouds, simulation_bytes, simulation_state, units_of, parse_args as collector_args  # noqa: E402


def congener_text(spec):
    """A --congeners value (congener_spec) as the collector's option."""
    return spec[1] if spec[0] == "groups" else str(spec[1])


def stream_spec(text):
    """--stream-above: auto, or GB (a number, 0 or more)."""
    if text.strip().lower() == "auto":
        return "auto"
    try:
        value = float(text)
    except ValueError:
        raise argparse.ArgumentTypeError(f"auto or a number of GB, got {text!r}") from None
    if value < 0:
        raise argparse.ArgumentTypeError(f"auto or a number of GB, 0 or more, got {text!r}")
    return value

STARTED = time.time()
CONSOLE = {"file": None, "early": []}  # OUTDIR/console.log, and the lines said before it was opened
SAYING = threading.RLock()  # say(): one line at a time from the main thread and the reports' thread


def say(message):
    """Prints a message, its first line headed by the time and how long the run has taken (and adds it to
    console.log)."""
    line = f"[{time.strftime('%H:%M:%S')} +{clock(time.time() - STARTED)}] {message}"
    with SAYING:  # the reports' thread says its lines beside the main thread's (main(): reports beside the packaging)
        print(line, flush=True)
        console_log(line)


def console_log(text):
    """Adds the console's text to OUTDIR/console.log, or keeps it until the log is open."""
    if CONSOLE["file"] is None:
        CONSOLE["early"].append(text)
        return
    CONSOLE["file"].write(text + "\n")
    CONSOLE["file"].flush()


def open_console_log(path):
    """OUTDIR/console.log: what the run says on the console, after a line of when and how it was started (a rerun
    adds to it)."""
    fh = open(path, "a")
    fh.write(f"--- {time.strftime('%Y-%m-%d %H:%M:%S')}: {' '.join(sys.argv)}\n")
    CONSOLE["file"] = fh
    for line in CONSOLE["early"]:
        console_log(line)
    CONSOLE["early"].clear()


class Steps:
    """The steps of the run on the console: "3/8 what it does (its log)" when one starts, and indented lines
    of how it went (done())."""
    total, current, title = 0, 0, "starting"

    @classmethod
    def start(cls, text):
        cls.current += 1
        cls.title = text.split(" (")[0].split(":")[0]  # for CpuTimeline's rows
        say(f"{cls.current}/{cls.total} {text}")

    @staticmethod
    def done(text):
        say("    " + text)


def gigabytes(size):
    return f"{size / 1e9:.1f} GB" if size >= 1e9 else f"{size / 1e6:.0f} MB"


def group_memory(group):
    """The resident memory of the processes of a process group (a job and the commands it started), in
    bytes; None without /proc (not Linux)."""
    total, found = 0, False
    try:
        pids = [p for p in os.listdir("/proc") if p.isdigit()]
    except OSError:
        return None
    for pid in pids:
        try:
            with open(f"/proc/{pid}/stat") as fh:
                fields = fh.read().rsplit(")", 1)[1].split()  # the fields after the command's name
        except (OSError, IndexError):
            continue  # ended meanwhile
        if len(fields) > 21 and fields[2] == str(group):  # the process group, and the resident pages
            total += int(fields[21]) * os.sysconf("SC_PAGE_SIZE")
            found = True
    return total if found else None


def db_size(db):
    path = os.path.join(db, "database.protal")
    return f"; database.protal {gigabytes(os.path.getsize(path))}" if os.path.isfile(path) else ""


def tree_size(path):
    """The bytes of the files under a folder."""
    total = 0
    for root, _, files in os.walk(path):
        for name in files:
            try:
                total += os.lstat(os.path.join(root, name)).st_size
            except OSError:
                pass
    return total


class CpuTimeline(threading.Thread):
    """logs/cpu_timeline.tsv: every `every` seconds, how many cores the run kept busy since the last row (its own and its
    commands' CPU time: the CPU usage of the cgroup the script runs in, which a SLURM job's commands share with it; on a
    machine without one, the whole node's busy time from /proc/stat, said in the source column), the node's cores
    waiting for I/O (/proc/stat: the network file system's waits show there), the cores the run may use, the step on the
    console and the commands running. Written as it goes, not said on the console; what shows where the run leaves cores
    idle (docs/claude/2026-10-09-build-parallelism)."""

    def __init__(self, path, every=15.0):
        super().__init__(name="cpu timeline", daemon=True)
        self.path, self.every, self.stopping = path, every, threading.Event()
        self.usage_file, self.scale = self.cgroup_usage()

    @staticmethod
    def cgroup_usage():
        """(the file of this process's cgroup CPU usage, seconds per its unit), or (None, None) without one: cgroup v2's
        cpu.stat (usage_usec), else v1's cpuacct.usage (ns)."""
        try:
            with open("/proc/self/cgroup") as fh:
                lines = [line.rstrip("\n").split(":", 2) for line in fh]
        except OSError:
            return None, None
        for _, controllers, path in lines:
            if controllers == "":
                candidate = os.path.join("/sys/fs/cgroup", path.lstrip("/"), "cpu.stat")
                if os.path.isfile(candidate):
                    return candidate, 1e-6
            elif "cpuacct" in controllers.split(","):
                for mount in ("cpu,cpuacct", "cpuacct,cpu", "cpuacct"):
                    candidate = os.path.join("/sys/fs/cgroup", mount, path.lstrip("/"), "cpuacct.usage")
                    if os.path.isfile(candidate):
                        return candidate, 1e-9
        return None, None

    def cgroup_seconds(self):
        try:
            with open(self.usage_file) as fh:
                text = fh.read()
            if self.usage_file.endswith("cpu.stat"):
                value = next(int(line.split()[1]) for line in text.splitlines() if line.startswith("usage_usec"))
            else:
                value = int(text.split()[0])
            return value * self.scale
        except (OSError, ValueError, StopIteration):
            return None

    @staticmethod
    def node_seconds():
        """(busy, I/O wait) seconds of the whole node since boot, from /proc/stat's cpu line."""
        try:
            with open("/proc/stat") as fh:
                fields = [int(v) for v in fh.readline().split()[1:]]
            tick = os.sysconf("SC_CLK_TCK")
            busy = fields[0] + fields[1] + fields[2] + sum(fields[5:8])  # user, nice, system, irq, softirq, steal
            return busy / tick, fields[4] / tick
        except (OSError, ValueError, IndexError):
            return None, None

    def sample(self):
        own = self.cgroup_seconds() if self.usage_file else None
        busy, iowait = self.node_seconds()
        return time.time(), own, busy, iowait

    def run(self):
        try:
            allowed = len(os.sched_getaffinity(0))
        except (AttributeError, OSError):
            allowed = os.cpu_count() or 0
        before = self.sample()
        try:
            new = not os.path.isfile(self.path)
            with open(self.path, "a") as fh:
                if new:
                    fh.write("time\telapsed_s\tcores_busy\tsource\tiowait_cores\tallowed_cores\tstep\trunning\n")
                while True:
                    stopped = self.stopping.wait(self.every)
                    now = self.sample()
                    span = now[0] - before[0]
                    if span > 0:
                        cgroup = now[1] is not None and before[1] is not None
                        busy = (now[1] - before[1]) if cgroup else \
                            (now[2] - before[2]) if now[2] is not None and before[2] is not None else float("nan")
                        iowait = (now[3] - before[3]) if now[3] is not None and before[3] is not None else float("nan")
                        running = "; ".join(f"{job.label}{' (paused)' if job.paused else ''}" for job in list(Job.running))
                        fh.write(f"{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(now[0]))}\t{now[0] - STARTED:.0f}\t"
                                 f"{busy / span:.2f}\t{'cgroup' if cgroup else 'node'}\t{iowait / span:.2f}\t{allowed}\t"
                                 f"{Steps.current}/{Steps.total} {Steps.title}\t{running}\n")
                        fh.flush()
                    before = now
                    if stopped:
                        break
        except OSError:
            pass

    def stop(self):
        self.stopping.set()
        self.join(timeout=30)


class Scratch:
    """The space the run takes on --scratch: what its file system holds now less what it held at the start,
    and the most that reached. Looked at every few seconds while the script waits (one statvfs)."""

    def __init__(self, path):
        self.path, self.start, self.peak = path, shutil.disk_usage(path).used, 0

    def look(self):
        usage = shutil.disk_usage(self.path)
        self.peak = max(self.peak, usage.used - self.start)
        return usage

    def text(self):
        usage = self.look()
        return (f"{gigabytes(max(0, usage.used - self.start))} in use by the run (at most {gigabytes(self.peak)}), "
                f"{gigabytes(usage.free)} free")


# ---- room on the samples' disk ---------------------------------------------------------------------------------------
# What a run puts on the disk of the simulated samples (--scratch, or OUTDIR) besides the reads it writes, estimated
# before the builds and again once the training database is built, when --stream-above auto chooses from the free space
# what to stream into protal (docs/claude/2026-10-07-build-ordering).
STORE_BYTES_PER_BASE = 0.25  # the genome store: 2 bits a base, N runs and contig names (~50 GB at r226)
# The samples' SAMs, profiles and training dumps against their reads as the collector estimates them: r226 v15 kept
# ~62 GB at the end for 591 GB of reads, which the collector's constants put at 734 GB.
KEPT_SHARE = 0.1
DB_PER_REFERENCE = 1.5  # database.protal against its reference.fna: r226 v14's training database 21.8 GB from 15.4 GB
HOST_PER_GZ_BYTE = 3.4  # the host genome as plain sequence (scenarios.prepare_host) against its gzipped FASTA


def on_disk(path, root):
    """Whether `path` (or, before it is made, the nearest folder above it) is on the file system of `root`."""
    while not os.path.exists(path):
        parent = os.path.dirname(path)
        if parent == path:
            return False
        path = parent
    return os.stat(path).st_dev == os.stat(root).st_dev


def genome_bases(table):
    """The bases of a genome table's genomes (its fourth column, genome_length; each FASTA once)."""
    bases = {}
    with open(table) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 4 and fields[3].isdigit():
                bases[fields[2]] = int(fields[3])
    return sum(bases.values())


def planned_bases(table, share):
    """The bases of the genomes to simulate from before the in-silico strains are made (insilico_strains.py, in the
    background): the table's genomes, and `share` of the one-genome species' genomes again (a strain is about as long as
    its representative)."""
    bases, by_species = {}, collections.defaultdict(list)
    with open(table) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 4 and fields[3].isdigit():
                bases[fields[2]] = int(fields[3])
                by_species[fields[1].split(";")[-1]].append(int(fields[3]))
    return sum(bases.values()) + max(0.0, min(1.0, share)) * sum(g[0] for g in by_species.values() if len(g) == 1)


def simulation_sizes(collections_):
    """[(largest sample, reads)] in bytes of every simulation of the collections' commands, as the collector estimates
    them (collect_training_data.simulation_bytes): what --stream-above compares, and what each writes unless streamed."""
    sizes = []
    for _, command, _ in collections_:
        sizes += simulation_bytes(units_of(collector_args(command[2:]))[1]).values()
    return sizes


def stream_threshold(sizes, room):
    """The --stream-above, in bytes, that streams the fewest simulations while the reads of the others fit in `room` bytes
    at once: sizes [(largest sample, reads)], a simulation streamed when its largest sample is above the value. -> 0 when
    everything fits; else a value halfway between two simulations' largest samples (below every one: all streamed); None
    when nothing fits (room < 0)."""
    if room < 0:
        return None
    if sum(reads for _, reads in sizes) <= room:
        return 0
    values = sorted({largest for largest, _ in sizes}, reverse=True)
    for above, below in zip(values, values[1:]):  # those of `above` and larger streamed
        if sum(reads for largest, reads in sizes if largest <= below) <= room:
            return (above + below) / 2
    return values[-1] / 2


def disk_needs(samples_root, genome_store, sim_table, databases, hosts, sizes, keep_free):
    """What the run puts on the disk of `samples_root` besides the reads it writes, in bytes by part, those on another
    file system left out: the genome store's growth (STORE_BYTES_PER_BASE of sim_table's genomes, less what the store
    holds); each database still to be built (databases: [(name, folder)], DB_PER_REFERENCE of its reference.fna); the
    host genome of each collection that has not prepared it (hosts: [(FASTA, folder)]); the samples' SAMs and profiles
    (KEPT_SHARE of all the reads, sizes as simulation_sizes); and keep_free bytes. -> {part: bytes}."""
    needs = {}
    if genome_store and on_disk(genome_store, samples_root):
        held = tree_size(genome_store) if os.path.isdir(genome_store) else 0
        bases = genome_bases(sim_table) if isinstance(sim_table, str) else sim_table  # or the bases planned
        needs["genome store"] = max(0.0, STORE_BYTES_PER_BASE * bases - held)
    for name, folder in databases:
        reference = os.path.join(folder, "reference.fna")
        if on_disk(folder, samples_root) and os.path.isfile(reference):
            needs[name] = DB_PER_REFERENCE * os.path.getsize(reference)
    for fasta, folder in hosts:
        if fasta and os.path.isfile(fasta) and not os.path.isdir(folder) and on_disk(folder, samples_root):
            needs["host genome"] = needs.get("host genome", 0) + \
                os.path.getsize(fasta) * (HOST_PER_GZ_BYTE if fasta.endswith(".gz") else 1)
    needs["SAMs and profiles"] = KEPT_SHARE * sum(reads for _, reads in sizes)
    if keep_free > 0:
        needs["--keep-free"] = keep_free
    return needs


def needs_text(needs):
    return ", ".join(f"{part} ~{gigabytes(size)}" for part, size in needs.items() if size >= 1e6) or "nothing"


class Job:
    """A command run with its output to log, in a process group of its own: when the script stops, however
    it stops (a failure, an exception, SIGTERM, SIGINT, SIGHUP), the command is stopped with what it started
    (the collector's simulator and protal runs), so that a rerun does not race a build left running.
    While the script waits for one job, it looks at the others every few seconds: one that failed (the
    finished database's build in the background) stops the script then, not hours later. Every
    progress_every seconds (if not 0) it says how each job is doing."""
    running = []
    progress_every = 0
    scratch = None  # a Scratch with --scratch
    cpu_log = None  # logs/cpu_jobs.tsv: a row per command that ended (cpu_row), not said on the console

    def __init__(self, command, log, on_success=None, label=None, nice=0, idle=False, append=False):
        """nice: the command's niceness, more than the script's (the simulations, beside the profiling); idle: the idle
        scheduling class instead (lower_priority: the finished database's build, on the cores the others leave); append:
        its output added to the log (a step of several commands, one after the other), not replacing it."""
        os.makedirs(os.path.dirname(log), exist_ok=True)
        self.command, self.log, self.on_success, self.started = command, log, on_success, time.time()
        self.label = label or os.path.basename(log).removesuffix(".log")
        self.seconds = None  # set when it has ended
        self.peak = None  # the most memory it, or a command it ran, took, in bytes; set when it has ended
        self.cpu = None  # its CPU seconds (user and system), and those of the commands it ran and waited for; when ended
        self.nice, self.idle = nice, idle
        self.paused = False
        self.fh = open(log, "a" if append else "w")
        if append:
            self.fh.write(f"--- {' '.join(command)}\n")
            self.fh.flush()
        self.process = subprocess.Popen(command, stdout=self.fh, stderr=subprocess.STDOUT, start_new_session=True,
                                        preexec_fn=lower_priority(nice, idle))
        Job.running.append(self)

    def poll(self):
        """The command's exit code once it has ended, else None. It is reaped with wait4, which tells its peak
        memory, or that of the commands it ran if larger."""
        if self.process.returncode is None:
            try:
                pid, status, usage = os.wait4(self.process.pid, os.WNOHANG)
            except ChildProcessError:
                return self.process.poll()
            if pid:
                self.process.returncode = os.waitstatus_to_exitcode(status)
                self.peak = usage.ru_maxrss * (1 if sys.platform == "darwin" else 1024)
                self.cpu = usage.ru_utime + usage.ru_stime
        return self.process.returncode

    def ended(self, rc):
        self.fh.close()
        Job.running.remove(self)
        self.seconds = time.time() - self.started
        self.cpu_row(rc)
        if rc:
            stop(f"Command failed ({rc}); see {self.log}: {' '.join(self.command)}")
        if self.on_success:
            self.on_success()

    def finish(self):
        """Waits for the command; the job, once ended. Stops the script if it or another job failed."""
        checked = told = time.time()
        while self.seconds is None:
            rc = self.poll()
            if rc is not None:
                self.ended(rc)
                break
            time.sleep(0.2)
            now = time.time()
            if now - checked >= 1:  # a wait4 per job and one statvfs
                check_jobs()
                if Job.scratch:
                    Job.scratch.look()
                checked = now
            if 0 < Job.progress_every <= now - told:
                for job in Job.running:
                    say(job.status())
                if Job.scratch:
                    say("scratch " + Job.scratch.text())
                told = now
        return self

    def cpu_row(self, rc):
        """Adds the ended command's row to Job.cpu_log (with its header if new): when it started and how long it ran (s),
        its CPU seconds and the cores it kept busy on average (CPU / wall), its priority, peak memory, exit code, label and
        log. What tells a step that leaves cores idle from one that fills them (docs/claude/2026-10-09-build-parallelism)."""
        if not Job.cpu_log or self.cpu is None:
            return
        try:
            new = not os.path.isfile(Job.cpu_log)
            with open(Job.cpu_log, "a") as fh:
                if new:
                    fh.write("started_s\twall_s\tcpu_s\tcores_busy\tpriority\tpeak_gb\texit\tlabel\tlog\n")
                fh.write(f"{self.started - STARTED:.0f}\t{self.seconds:.1f}\t{self.cpu:.1f}\t"
                         f"{self.cpu / max(self.seconds, 1e-9):.2f}\t{'idle' if self.idle else f'nice {self.nice}'}\t"
                         f"{(self.peak or 0) / 1e9:.2f}\t{rc}\t{self.label}\t{os.path.basename(self.log)}\n")
        except OSError:
            pass

    def status(self):
        """How the command is doing: how long it has run, the memory it and the commands it started take, and
        the last line of its log."""
        memory = group_memory(self.process.pid)
        line = last_line(self.log, 200)
        return (f"{self.label}: {clock(time.time() - self.started)} so far" + (" (paused)" if self.paused else "") +
                (f", {gigabytes(memory)} in memory" if memory else "") +
                (f"; {os.path.basename(self.log)}: {line}" if line else ""))

    def took(self):
        """How long it took, and its peak memory."""
        return clock(self.seconds) + (f", peak memory {gigabytes(self.peak)}" if self.peak else "")

    def pause(self):
        """Stops the command and what it started (SIGSTOP to its group) until resume(): its memory stays, its cores and
        memory bandwidth go to what runs meanwhile. -> whether it was running."""
        if self.seconds is not None or self.paused or self.poll() is not None:
            return False
        try:
            os.killpg(self.process.pid, signal.SIGSTOP)
        except (ProcessLookupError, PermissionError):
            return False
        self.paused = True
        return True

    def resume(self):
        if self.paused:
            with contextlib.suppress(ProcessLookupError, PermissionError):
                os.killpg(self.process.pid, signal.SIGCONT)
            self.paused = False

    def kill(self):
        """Stops the command and every process of its group (a paused one continued, so that it takes the SIGTERM)."""
        for sig, wait in ((signal.SIGTERM, 60), (signal.SIGKILL, None)):
            try:
                os.killpg(self.process.pid, sig)
                if self.paused:
                    os.killpg(self.process.pid, signal.SIGCONT)
            except (ProcessLookupError, PermissionError):
                break
            try:
                self.process.wait(timeout=wait)
                break
            except subprocess.TimeoutExpired:
                continue
        self.fh.close()


def lower_priority(nice=0, idle=False):
    """A preexec_fn for a command run at a lower priority than the script: niceness raised by `nice`; idle, the idle
    scheduling class (SCHED_IDLE: its threads run only on cores no other thread wants, and give way at once), or the
    highest niceness where that class is missing. None for neither."""
    if not nice and not idle:
        return None

    def apply():
        if idle:
            try:
                os.sched_setscheduler(0, os.SCHED_IDLE, os.sched_param(0))
                return
            except (AttributeError, OSError):
                os.nice(19)
                return
        os.nice(nice)
    return apply


def check_jobs():
    """Ends the jobs whose command has ended: one that failed stops the script."""
    for job in list(Job.running):
        rc = job.poll()
        if rc is not None:
            job.ended(rc)


def run(command, log, on_success=None, label=None, append=False):
    return Job(command, log, on_success, label, append=append).finish()


def log_says(path, text):
    """Whether the log at `path` holds `text` (the line of a step that the run waits for)."""
    try:
        with open(path, errors="replace") as fh:
            return text in fh.read()
    except OSError:
        return False


def make_genome_table(gtdb, release, output, extra_dirs=(), species=None, threads=1):
    """The simulator's genome table (genome_rows, with the lengths: write_genome_table) at `output`. -> its genomes."""
    rows = genome_rows(gtdb, release, extra_dirs, species)
    write_genome_table(rows, output, threads)
    return len(rows)


GENOME_SUFFIXES = (".fna", ".fna.gz", ".fa", ".fa.gz", ".fasta", ".fasta.gz")  # by precedence for one accession


def genome_rows(gtdb, release, extra_dirs=(), species=None):
    """The rows [accession, lineage, FASTA path] of every genome FASTA of a GTDB species found in the release (and in
    extra_dirs), or only those of the species in the set `species`, by accession. Each folder is walked once (six
    recursive globs walked a network file system's 50,000 genomes six times); an accession with several files takes the
    first folder's, there the first suffix's of GENOME_SUFFIXES, then the first path in sorted order. Hidden files and
    folders are left out, as glob leaves them out."""
    taxonomy = {}
    for domain in ("bac120", "ar53"):
        for suffix in (".tsv", ".tsv.gz"):
            path = os.path.join(gtdb, f"{domain}_taxonomy_r{release}{suffix}")
            if os.path.isfile(path):
                opener = gzip.open if path.endswith(".gz") else open
                with opener(path, "rt") as fh:
                    for line in fh:
                        fields = line.rstrip("\n").split("\t")
                        if len(fields) >= 2:
                            match = ACCESSION.search(fields[0])
                            if match and (species is None or fields[1].split(";")[-1] in species):
                                taxonomy[match.group(1)] = fields[1]
                break
    genome_dirs = [os.path.join(gtdb, "genomic_files_all", f"gtdb_genomes_all_r{release}"),
                   os.path.join(gtdb, "genomic_files_reps", f"gtdb_genomes_reps_r{release}"), *extra_dirs]
    paths = {}  # accession -> (folder's rank, suffix's rank, path)
    for rank, root in enumerate(genome_dirs):
        if not os.path.isdir(root):
            continue
        found = {}
        for folder, subfolders, files in os.walk(root, followlinks=True):
            subfolders[:] = sorted(d for d in subfolders if not d.startswith("."))
            for name in files:
                suffix = next((i for i, s in enumerate(GENOME_SUFFIXES) if name.endswith(s)), None)
                m = ACCESSION.search(name) if suffix is not None and not name.startswith(".") else None
                if m and m.group(1) in taxonomy:
                    candidate = (suffix, os.path.join(folder, name))
                    if m.group(1) not in found or candidate < found[m.group(1)]:
                        found[m.group(1)] = candidate
        for accession, (suffix, path) in found.items():
            paths.setdefault(accession, (rank, suffix, path))
    if not paths:
        sys.exit("No extracted whole genome FASTAs found under genomic_files_all/gtdb_genomes_all_r" + release +
                 " or genomic_files_reps/gtdb_genomes_reps_r" + release + "; extract GTDB genome files or pass --genome-table")
    return [[a, taxonomy[a], os.path.abspath(paths[a][2])] for a in sorted(paths)]


class GenomeCache:
    """What a build reads every genome FASTA for, kept per input folder across builds (--genome-cache): its length
    (genome_length) and its contigs' names (the reports trace the reads' sources through them), by the FASTA's path, size
    and modification time. At r226 the pass over the ~50,000 genomes on the network file system took 8-22 min per build
    (docs/claude/2026-10-09-build-parallelism); a rebuild from the same --inputs reads none of them again. A gzipped
    TSV: path, size, mtime_ns, length, contig names (space-separated). Builds side by side merge what they add (save,
    under a lock)."""

    def __init__(self, path):
        self.path, self.facts, self.added = path, {}, {}
        if path and os.path.isfile(path):
            try:
                with self.lock(fcntl.LOCK_SH), gzip.open(path, "rt") as fh:
                    for line in fh:
                        f = line.rstrip("\n").split("\t")
                        if len(f) == 5 and f[3].isdigit():
                            self.facts[f[0]] = (f"{f[1]}:{f[2]}", int(f[3]), f[4].split(" ") if f[4] else [])
            except (OSError, EOFError, ValueError):  # a damaged cache is read again from the genomes
                self.facts = {}

    @contextlib.contextmanager
    def lock(self, how):
        with open(self.path + ".lock", "a") as fh:
            fcntl.flock(fh, how)
            try:
                yield
            finally:
                fcntl.flock(fh, fcntl.LOCK_UN)

    def get(self, path, stamp):
        """(length, contig names) of the FASTA at `path` if its stamp ("size:mtime_ns") is the one cached, else None."""
        fact = self.facts.get(path)
        return fact[1:] if fact and stamp and fact[0] == stamp else None

    def add(self, path, stamp, length, names):
        if stamp:
            self.facts[path] = self.added[path] = (stamp, length, names)

    def save(self):
        """Writes the cache with what this run added, merged with what the file holds now (another build may have added
        to it meanwhile)."""
        if not self.path or not self.added:
            return
        os.makedirs(os.path.dirname(os.path.abspath(self.path)), exist_ok=True)
        with self.lock(fcntl.LOCK_EX):
            merged = GenomeCache.__new__(GenomeCache)
            merged.path, merged.facts, merged.added = self.path, {}, {}
            if os.path.isfile(self.path):
                try:
                    with gzip.open(self.path, "rt") as fh:
                        for line in fh:
                            f = line.rstrip("\n").split("\t")
                            if len(f) == 5 and f[3].isdigit():
                                merged.facts[f[0]] = (f"{f[1]}:{f[2]}", int(f[3]), f[4].split(" ") if f[4] else [])
                except (OSError, EOFError, ValueError):
                    merged.facts = {}
            merged.facts.update(self.added)
            with gzip.open(self.path + ".partial", "wt", compresslevel=1) as fh:
                for path in sorted(merged.facts):
                    stamp, length, names = merged.facts[path]
                    size, mtime = stamp.split(":")
                    fh.write(f"{path}\t{size}\t{mtime}\t{length}\t{' '.join(names)}\n")
            os.replace(self.path + ".partial", self.path)
        self.added = {}


def file_stamps(paths, threads=1):
    """{path: "size:mtime_ns"} ("" for a file that cannot be read): the stamp genome_contigs' cache and GenomeCache key a
    FASTA by, taken on threads (one stat is a round trip to a network file system)."""
    def stamp(path):
        try:
            st = os.stat(path)
            return f"{st.st_size}:{st.st_mtime_ns}"
        except OSError:
            return ""
    paths = sorted(set(paths))
    with concurrent.futures.ThreadPoolExecutor(max(1, min(32, threads, len(paths) or 1))) as pool:
        return dict(zip(paths, pool.map(stamp, paths)))


def write_genome_table(rows, output, threads, cache=None, contigs=None):
    """Writes a genome table, rows [accession, taxonomy, FASTA path, ...], with each genome's length (genome_length) as
    the fourth column. A genome's length (and its contigs' names) comes from the GenomeCache `cache` when its FASTA is
    the one cached, else from one read of the FASTA (in processes: counting the letters holds the GIL), which also gives
    the names for the cache; a row that has its length (a given table's fourth column) keeps it. With `contigs` (the
    reports' contig cache, trace_relatives.genome_contigs), the names known are added to it, so that the reports read no
    genome again."""
    stamps = file_stamps((r[2] for r in rows), threads)
    facts = {}  # path -> (length, names or None)
    for r in rows:
        if r[2] in facts:
            continue
        cached = cache.get(r[2], stamps[r[2]]) if cache else None
        if cached:
            facts[r[2]] = cached
        elif len(r) >= 4 and r[3].isdigit():
            facts[r[2]] = (int(r[3]), None)
    todo = sorted({r[2] for r in rows if r[2] not in facts})
    if todo:
        with concurrent.futures.ProcessPoolExecutor(max(1, min(threads, len(todo)))) as pool:
            for path, fact in zip(todo, pool.map(genome_facts, todo, chunksize=8)):
                facts[path] = fact
                if cache:
                    cache.add(path, stamps[path], *fact)
    if cache:
        cache.save()
    with open(output + ".partial", "w") as fh:
        fh.writelines("\t".join(r[:3] + [str(facts[r[2]][0])]) + "\n" for r in rows)
    os.replace(output + ".partial", output)
    if contigs:
        import trace_relatives
        trace_relatives.add_to_contig_cache(contigs, [(path, stamps[path], names) for path, (_, names) in sorted(facts.items())
                                                      if names is not None and stamps[path]])


def table_rows(table):
    """The rows of a genome table that name a genome (its header and comments left out): [accession, lineage, FASTA
    path, ...]."""
    with open(table) as fh:
        return [line.rstrip("\n").split("\t") for line in fh if line.strip() and not line.startswith("#")]


def with_lengths(table, output, threads, cache=None, contigs=None):
    """The genome table to simulate from: `table` itself if every row has a genome length, else a copy with them
    (write_genome_table) at `output`. Rows that are no genome (a header, comments) are left out of the copy;
    the simulator reads a table of four columns without a header."""
    rows = table_rows(table)
    genomes = [r for r in rows if len(r) >= 3 and os.path.isfile(r[2])]
    if len(genomes) == len(rows) and all(len(r) >= 4 and r[3].isdigit() for r in rows):
        return table
    write_genome_table(genomes, output, threads, cache, contigs)
    return output


NON_LETTERS = bytes(b for b in range(256) if not (65 <= b <= 90 or 97 <= b <= 122))


def genome_length(path):
    """A genome's length as simulate_metagenomes counts it (read_genome_length): the letters of the lines that
    do not start with '>', of the gzip-decompressed file if its name ends in .gz."""
    return genome_facts(path)[0]


def genome_facts(path):
    """(length, contig names) of a genome FASTA, in one read: its length as genome_length counts it, and the first word
    of each header line (its contigs, which name the simulated reads; trace_relatives.contig_names)."""
    with open(path, "rb") as fh:
        data = fh.read()
    if path.endswith(".gz"):
        data = gzip.decompress(data)
    letters, names = len(data.translate(None, NON_LETTERS)), []

    def next_header(after):  # the start of the first header line at or after `after` (a line's start), -1 for none
        if data.startswith(b">", after):
            return after
        i = data.find(b"\n>", after)
        return i + 1 if i >= 0 else -1

    at = next_header(0)
    while at >= 0:
        end = data.find(b"\n", at)
        header = data[at:end] if end >= 0 else data[at:]
        letters -= len(header.translate(None, NON_LETTERS))
        words = header[1:].split(None, 1)
        if words:
            names.append(words[0].decode(errors="replace"))
        at = next_header(end + 1) if end >= 0 else -1
    return letters, names


def summarize_genome_table(path, reps, insilico_to_come=False):
    """What the simulations can draw: genomes and species by domain, and how often a simulated species is
    not the database's representative genome (the simulator picks a species, then one of its genomes): a real
    strain, or an in-silico strain (insilico_strains.py; INSILICO_PREFIX). The lines of genome_table.txt, the same
    in one line for the console, and a warning ("" for none): nearly every simulated species its reference, or most
    strains in-silico. insilico_to_come: the table before the in-silico strains (step 3), judged after them."""
    genomes, insilico, rep_count = collections.Counter(), collections.Counter(), collections.Counter()
    domain = {}
    with open(path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            lineage = next((f for f in fields if f.startswith("d__") and ";s__" in f), None)
            if lineage is None:
                continue
            species = lineage.split(";")[-1]
            genomes[species] += 1
            domain[species] = lineage.split(";")[0][3:]
            insilico[species] += fields[0].startswith(INSILICO_PREFIX)
            rep_count[species] += reps is not None and normalize_accession(fields[0]) in reps
    lines = [f"genome table {path}: {sum(genomes.values())} genomes of {len(genomes)} species"]
    for d in sorted(set(domain.values())):
        sp = [s for s in genomes if domain[s] == d]
        lines.append(f"  {d}: {len(sp)} species, {sum(genomes[s] for s in sp)} genomes")
    several = sum(1 for n in genomes.values() if n > 1)
    # A species' genomes other than its representative (without the representatives known: all but one), real or in silico.
    real = sum((n - insilico[s] - (rep_count[s] if reps is not None else 1)) / n for s, n in genomes.items())
    real = max(0.0, real) / max(1, len(genomes))
    made = sum(insilico[s] / n for s, n in genomes.items()) / max(1, len(genomes))
    other_strain = real + made
    split = f" (a real strain {100 * real:.1f}%, an in-silico strain {100 * made:.1f}%)" if made else ""
    lines.append(f"  {several} species have more than one genome; a simulated species is another genome than "
                 f"its representative {100 * other_strain:.1f}% of the time{split}")
    if reps is not None:
        lines.append(f"  {sum(rep_count.values())} of the genomes are species representatives (the database's references)")
    warning = ""
    if other_strain < 0.2 and not insilico_to_come:
        warning = ("WARNING: nearly all simulated species will be the database's own reference genome, closer to it "
                   "than real strains are. Add non-representative genomes (scripts/download_gtdb.py, then --inputs) to "
                   "train on real strain divergence.")
    elif made > real:
        warning = (f"WARNING: most simulated strains are in-silico ({100 * made:.1f}% of the simulated species against "
                   f"{100 * real:.1f}% real): the models learn strains mostly from in-silico ones, which they miss less "
                   "often than real ones (r226 v12: 2.5% against 6.1% with species held out; "
                   "docs/claude/2026-10-05-r226-v12-scenarios). Download more species with strains (download_gtdb.py "
                   "--species, then --inputs), or give fewer one-genome species an in-silico strain (--insilico-strains).")
    if warning:
        lines.append("  " + warning)
    elif insilico_to_come and other_strain < 0.2:
        lines.append("  the species with one genome get in-silico strains next (step 3)")
    domains = ", ".join(f"{d} {sum(1 for s in genomes if domain[s] == d)}" for d in sorted(set(domain.values())))
    brief = (f"{sum(genomes.values())} genomes of {len(genomes)} species ({domains}); a simulated species is another "
             f"genome than its representative {100 * other_strain:.1f}% of the time{split}")
    return lines, brief, warning


CLADE_RANKS = ("phylum", "class", "order", "family", "genus")


def parse_clades(text):
    """--holdout-clades "phylum:2,class:4,family:8" -> {"phylum": 2, ...}; "none" or "" -> {}."""
    clades = {}
    for part in (text or "").split(","):
        part = part.strip()
        if not part or part == "none":
            continue
        rank, _, count = part.partition(":")
        if rank not in CLADE_RANKS or not count.isdigit():
            sys.exit(f"--holdout-clades: expected RANK:COUNT with a rank of {', '.join(CLADE_RANKS)}, got {part!r}")
        clades[rank] = int(count)
    return clades


HOLDOUT_COMPLEX_DISTANCE = 0.01  # --holdout-complex-distance: a strain's distance on the marker genes
CLOUD_MAX_DISTANCE = 0.15  # protal's species_neighbours::kMaxDistance: farther congeners are not listed
CLOUD_BANDS = ((0.01, "within 0.01"), (0.02, "0.01-0.02"), (0.05, "0.02-0.05"), (CLOUD_MAX_DISTANCE, "0.05-0.15"))
CLOUD_BAND_NONE = "farther or none"


def cloud_band(distance):
    """The CLOUD_BANDS label of a distance to the nearest kept congener (None: none listed)."""
    if distance is None:
        return CLOUD_BAND_NONE
    return next((label for limit, label in CLOUD_BANDS if distance <= limit), CLOUD_BAND_NONE)


def species_complexes(clouds, distance):
    """{species: complex number} of the species within `distance` of a congener (read_clouds' distances, from either
    side): the connected components of those pairs, numbered from 1 in the order of their smallest species name.
    Species without such a congener are not listed."""
    parent = {}

    def find(s):
        while parent[s] != s:
            parent[s] = parent[parent[s]]
            s = parent[s]
        return s

    for species, congeners in clouds.items():
        for congener, d in congeners:
            if d <= distance:
                parent.setdefault(species, species)
                parent.setdefault(congener, congener)
                a, b = find(species), find(congener)
                if a != b:
                    parent[max(a, b)] = min(a, b)
    number = {root: i + 1 for i, root in enumerate(sorted({find(s) for s in parent}))}
    return {s: number[find(s)] for s in parent}


def nearest_kept(species, clouds, chosen):
    """The distance to the nearest congener of `species` not in `chosen` (kept in the training database), by the
    clouds; None without one listed (none within CLOUD_MAX_DISTANCE, or all held out)."""
    return min((d for c, d in clouds.get(species, ()) if c not in chosen), default=None)


def write_holdout(path, chosen, clouds=None, complex_distance=0):
    """heldout_species.txt: species, rank, clade; with clouds two more columns, the distance to the nearest kept
    congener ("-": none within CLOUD_MAX_DISTANCE) and the complex the species is held out with (c<number>; empty
    alone)."""
    complexes = species_complexes(clouds, complex_distance) if clouds and complex_distance > 0 else {}
    with open(path, "w") as fh:
        for s, (rank, clade) in sorted(chosen.items()):
            line = f"{s}\t{rank}\t{clade}"
            if clouds is not None:
                d = nearest_kept(s, clouds, chosen)
                line += f"\t{'-' if d is None else f'{d:.4f}'}\t{'c' + str(complexes[s]) if s in complexes else ''}"
            fh.write(line + "\n")


def read_holdout_details(path):
    """{species: (distance to the nearest kept congener or None, complex or "")} from the two columns write_holdout adds
    with clouds; {} for a file without them."""
    details = {}
    with open(path) as fh:
        for line in fh:
            fields = [f.strip() for f in line.rstrip("\n").split("\t")]
            if not fields[0] or fields[0].startswith("#") or len(fields) < 5:
                continue
            species = fields[0] if fields[0].startswith("s__") else "s__" + fields[0]
            details[species] = (None if fields[3] in ("", "-") else float(fields[3]), fields[4])
    return details


def species_neighbours_summary(log):
    """What protal --write_species_neighbours said (its "Species neighbours:" line, less the file it wrote)."""
    try:
        with open(log) as fh:
            lines = [line.strip() for line in fh if line.startswith("Species neighbours:")]
    except OSError:
        return "unknown (no log)"
    if not lines:
        return "unknown (no 'Species neighbours' line)"
    return lines[-1].removeprefix("Species neighbours:").strip().rsplit(": ", 1)[0]


def choose_holdout(genome_table, taxonomy, fraction, clades, max_share, seed, clouds=None, complex_distance=0):
    """The species a training database leaves out, as {species: (rank, clade)}; rank "species" for single
    species. First whole clades, clades[rank] of each rank from phylum down: drawn among those with at least
    two species the genome table can simulate (so that samples can have them), with at most max_share of the
    database's species, and in no clade drawn before. Then a random `fraction` of the species the genome
    table can simulate that no clade took, the same fraction in each domain.

    With `clouds` (read_clouds: each species' nearest congeners and their distances) and complex_distance > 0, species
    within that distance of a congener form complexes (species_complexes) that are held out or kept whole: the draw
    takes units, a complex (every species of it, also those the genome table cannot simulate) or a species alone, in
    random order until the fraction of the domain's species to simulate is reached, so that no held-out species leaves
    a near-identical twin in the training database (its reads would teach "absent" at the identity where a divergent
    strain teaches "present"). Without complexes the draw is the one of before 2026-10-08."""
    lineage_by_id, _ = lineages.from_taxonomy(taxonomy)
    db = {lin["species"]: lin for lin in lineage_by_id.values() if "species" in lin}
    complexes = species_complexes(clouds, complex_distance) if clouds and complex_distance > 0 else {}
    members_of = collections.defaultdict(set)
    for s, c in complexes.items():
        if s in db:
            members_of[c].add(s)
    pool_by_domain = collections.defaultdict(set)
    with open(genome_table) as fh:
        for line in fh:
            lineage = next((f for f in line.rstrip("\n").split("\t") if f.startswith("d__") and ";s__" in f), None)
            if lineage and lineage.split(";")[-1] in db:
                pool_by_domain[lineage.split(";")[0]].add(lineage.split(";")[-1])
    pool = set().union(*pool_by_domain.values()) if pool_by_domain else set()
    rng = random.Random(seed)
    chosen = {}
    for rank in CLADE_RANKS:
        if not clades.get(rank):
            continue
        members = collections.defaultdict(set)
        for species, lin in db.items():
            if rank in lin:
                members[lin[rank]].add(species)
        eligible = sorted(c for c, species in members.items() if len(species & pool) >= 2
                          and len(species) <= max_share * len(db) and not any(s in chosen for s in species))
        picked = rng.sample(eligible, min(clades[rank], len(eligible)))
        if len(picked) < clades[rank]:
            say(f"WARNING: only {len(eligible)} {rank} clades can be held out (two or more species to simulate, at most "
                f"{max_share:.0%} of the species); holding out {len(picked)} instead of {clades[rank]}")
        for clade in sorted(picked):
            for species in members[clade]:
                chosen[species] = (rank, clade)
    for domain in sorted(pool_by_domain):
        species = sorted(s for s in pool_by_domain[domain] if s not in chosen)
        if not complexes:
            for s in rng.sample(species, round(fraction * len(species))):
                chosen[s] = ("species", s)
            continue
        units, seen = [], set()
        for s in species:
            if s not in complexes:
                units.append([s])
            elif complexes[s] not in seen:
                seen.add(complexes[s])
                units.append(sorted(members_of[complexes[s]] - set(chosen)))
        rng.shuffle(units)
        target, taken = round(fraction * len(species)), 0
        for members in units:
            if taken >= target:
                break
            for s in members:
                chosen[s] = ("species", s)
            taken += sum(s in pool_by_domain[domain] for s in members)
    return chosen


def read_holdout(path):
    """heldout_species.txt: species, and optionally the rank it was held out at and the clade (a species
    alone: rank "species") -> {species: (rank, clade)}."""
    chosen = {}
    with open(path) as fh:
        for line in fh:
            fields = [f.strip() for f in line.rstrip("\n").split("\t")]
            if not fields[0] or fields[0].startswith("#"):
                continue
            species = fields[0] if fields[0].startswith("s__") else "s__" + fields[0]
            rank = fields[1] if len(fields) > 1 and fields[1] else "species"
            chosen[species] = (rank, fields[2] if len(fields) > 2 and fields[2] else species)
    return chosen


def pool_species(genome_table):
    """The species a genome table can simulate."""
    species = set()
    with open(genome_table) as fh:
        for line in fh:
            lineage = next((f for f in line.rstrip("\n").split("\t") if f.startswith("d__") and ";s__" in f), None)
            if lineage:
                species.add(lineage.split(";")[-1])
    return species


def describe_holdout(chosen, pool, details=None, complex_distance=0):
    """Lines saying what the training database leaves out; with details (read_holdout_details) also how near the
    species held out alone are to their nearest kept congener, and the complexes held out whole."""
    lines = []
    by_rank = collections.defaultdict(lambda: collections.defaultdict(set))
    for species, (rank, clade) in chosen.items():
        by_rank[rank][clade].add(species)
    for rank in (*CLADE_RANKS, "species"):
        if rank not in by_rank:
            continue
        species = set().union(*by_rank[rank].values())
        simulated = len(species & pool)
        if rank == "species":
            lines.append(f"  {len(species)} single species ({simulated} to simulate)")
        else:
            names = ", ".join(f"{c} ({len(s)} species, {len(s & pool)} to simulate)" for c, s in sorted(by_rank[rank].items()))
            lines.append(f"  {len(by_rank[rank])} {rank} clades, {len(species)} species ({simulated} to simulate): {names}")
    if details:
        alone = sorted(s for s, (rank, _) in chosen.items() if rank == "species")
        bands = collections.Counter(cloud_band(details[s][0]) if s in details else CLOUD_BAND_NONE for s in alone)
        lines.append("  nearest kept congener of the species held out alone (species_clouds.tsv): " +
                     ", ".join(f"{label} {bands[label]}" for label in (*(label for _, label in CLOUD_BANDS), CLOUD_BAND_NONE)))
        complexes = collections.defaultdict(set)
        for s in alone:
            if s in details and details[s][1]:
                complexes[details[s][1]].add(s)
        if complex_distance > 0:
            members = set().union(*complexes.values()) if complexes else set()
            lines.append(f"  species complexes (congeners within {complex_distance:g}) held out whole: {len(complexes)}, "
                         f"{len(members)} species ({len(members & pool)} to simulate); none torn")
    return lines


def holdout_brief(chosen):
    """What the training database leaves out, in a few words: the clades per rank and the species alone."""
    clades = collections.Counter(rank for rank, _ in set(chosen.values()) if rank != "species")
    in_clades = sum(1 for rank, _ in chosen.values() if rank != "species")
    alone = len(chosen) - in_clades
    parts = [f"{clades[r]} {r}" for r in CLADE_RANKS if clades.get(r)]
    text = f"{', '.join(parts)} clades ({in_clades} species)" if parts else ""
    return text + (" and " if text and alone else "") + (f"{alone} species alone" if alone else "")


def gene_conservation_summary(build_log):
    """What protal --build said of the genes' conservation factors (its "Gene conservation:" line, less the
    file it wrote), for build_metadata.tsv: the factors are in the database, so that the depth margin can be
    scaled by them (protal --gene_conservation db) without a rebuild."""
    try:
        with open(build_log) as fh:
            lines = [line.strip() for line in fh if line.startswith("Gene conservation:")]
    except OSError:
        return "unknown (no build log)"
    if not lines:
        return "none (this protal does not estimate them)"
    text = lines[-1].removeprefix("Gene conservation:").strip()
    return text.rsplit(": ", 1)[0] if text.endswith("gene_conservation.tsv") else text


# --maxnodes by default: the trainer's for boosting (machine_learning_cmdline.MODEL_DEFAULTS), and the forests' since the
# r226 v3 training; boosting's rounds (machine_learning_cmdline.GBM_ROUNDS).
BOOSTED_LEAVES, FOREST_LEAVES, TRAINER_ROUNDS = "63", "512,pb:128,ont:128", 250


def max_leaves(spec, read_type):
    """The leaves per tree of read_type's model from --maxnodes: a TYPE:N for it, else the bare N, else 256 (the
    trainer's default). ValueError when the list is not of N and TYPE:N."""
    default, given = 256, {}
    for item in (i.strip() for i in spec.split(",") if i.strip()):
        name, _, count = item.rpartition(":")
        if not count.isdigit() or (name and name not in TABLES):
            raise ValueError(f"--maxnodes: N or TYPE:N items ({', '.join(TABLES)}), got {item!r}")
        if name:
            given[name] = int(count)
        else:
            default = int(count)
    return given.get(read_type, default)


def depth_knob_types(args):
    """The read types whose models get knobs by sample depth (--depth-knob-read-types)."""
    return {t.strip() for t in args.depth_knob_read_types.split(",") if t.strip()}


def gene_congeners_summary(build_log):
    """What protal --build said of how the genes differ between congeners (its "Gene congeners:" line, less the
    file it wrote; gene_congeners.tsv beside the database and in model_logs/), for build_metadata.tsv."""
    try:
        with open(build_log) as fh:
            lines = [line.strip() for line in fh if line.startswith("Gene congeners:")]
    except OSError:
        return "unknown (no build log)"
    if not lines:
        return "none (this protal does not compare them)"
    text = lines[-1].removeprefix("Gene congeners:").strip()
    return text.rsplit(": ", 1)[0] if text.endswith("gene_congeners.tsv") else text


def suspect_copies_summary(build_log):
    """What protal --build said of the gene copies near-identical to another genus's (its "Suspect copies:" line, less
    the files it wrote; suspect_copies.tsv in the database, gene_incongruence.tsv beside it and in model_logs/), for
    build_metadata.tsv."""
    try:
        with open(build_log) as fh:
            lines = [line.strip() for line in fh if line.startswith("Suspect copies:")]
    except OSError:
        return "unknown (no build log)"
    if not lines:
        return "none (this protal does not look for them)"
    text = lines[-1].removeprefix("Suspect copies:").strip()
    return re.sub(r":? ?\S*(suspect_copies|gene_incongruence)\.tsv", "", text).strip()


def error_read_units(text, defs):
    """[(read type, scope)] of --error-reads: all (every read type; scope None: all its samples), none, or READ_TYPE,
    READ_TYPE:design (the design's samples) or READ_TYPE:SCENARIO, comma-separated; ValueError for an entry that is
    none of these, or of an unknown scenario."""
    out = []
    for entry in (e.strip() for e in (text or "").split(",")):
        if not entry or entry == "none":
            continue
        if entry == "all":
            out += [(kind, None) for kind in TABLES]
            continue
        kind, colon, scope = entry.partition(":")
        if kind not in TABLES or (colon and not scope):
            raise ValueError(f"{entry!r}: expected all, READ_TYPE, READ_TYPE:design or READ_TYPE:SCENARIO (read types "
                             f"{', '.join(TABLES)})")
        if scope and scope != "design" and scope not in defs:
            raise ValueError(f"{entry!r}: no scenario {scope!r} (scenarios: {', '.join(defs)})")
        out.append((kind, scope or None))
    return list(dict.fromkeys(out))


def reports(trace, units, prefixes, training, test, training_db, logs, steps, threads=1, contig_cache=None, sams=False,
            heldout=None, protal=None, taxonomy=None, work=None):
    """What the models' errors rest on, once the database is ready (neither feeds it), side by side, each on its share of
    the threads, each with its log in `steps` (OUTDIR/logs):
    - trace (heldout_species.txt, or None): model_logs/relatives_by_gene_conservation.txt (trace_relatives.py): where
      the paired-end reads of the species the training database lacks land, by the genes' conservation factors;
    - units (--error-reads): model_logs/error_reads/<read type>/ (error_reads.py, one run per read type, each on its share
      of error_reads.py's memory budget): a table of each model's false positives and false negatives in those samples
      and where their reads went, from the samples' SAMs, which protal wrote with an unmapped record for every read that
      seeded on taxa but aligned nowhere (collect_training_data.py --unmapped_reads); with sams (--share-logs) also the
      SAM records of those reads, the FP's and the FN's in files of their own;
    - then, with sams and the SAMs written: model_logs/ancestry_sites/<read type>.* (ancestry_reports): which side those
      reads take where the species differs from its congeners, from the training database's genes (protal --unpack_db
      into work, unless its reference.fna is still there) and its full reference's other genomes (if kept).
    The genomes' contig names, which both need for the paired-end samples, are read once before them into contig_cache.
    A failure is reported, and does not stop the build."""
    import error_reads as error_reads_script
    import trace_relatives as trace_script
    began = time.time()
    scopes = collections.defaultdict(list)
    for kind, scope in units:
        scopes[kind].append(scope)
    tasks, told = [], {}  # tasks: (what, command, log, error-read output folder or None)
    for kind, which in scopes.items():
        calls, out = prefixes[kind] + ".calls.tsv.gz", os.path.join(logs, "error_reads", kind)
        if not os.path.isfile(calls):
            told[kind] = f"{kind}: the trainer wrote no calls"
            continue
        command = [sys.executable, ERROR_READS, "--calls", calls, "--training", training, "--db", training_db,
                   "--samples", "all" if None in which else ",".join(which), "--read-type", kind, "--out", out,
                   "--sams", "FP,FN" if sams else "none"]
        command += ["--test", test] if test and os.path.isdir(test) else []
        command += ["--contig-cache", contig_cache] if contig_cache else []
        command += ["--heldout", heldout] if heldout and os.path.isfile(heldout) else []
        tasks.append((kind, command, os.path.join(steps, f"error_reads_{kind}.log"), out))
    relatives = os.path.join(logs, "relatives_by_gene_conservation")
    if trace:
        tasks.insert(0, ("trace", [sys.executable, TRACE, "--points", os.path.join(training, "points"), "--db", training_db,
                                   "--heldout", trace, "--out", relatives] +
                         (["--contig-cache", contig_cache] if contig_cache else []),
                         os.path.join(steps, "trace_relatives.log"), None))
    if not tasks:
        if told:
            say(f"    the reads of the models' errors: {'; '.join(told.values())}")
        return
    # Each report on all the threads: they differ in work (r226 v19, 16 threads each: pe 186 s, se 150 s, pb and ont
    # ~80 s), and a report that ends leaves its cores to the others, as a split by an estimate would not. Their memory
    # is split (error_reads.py keeps its workers within its share).
    share = max(1, threads)
    error_runs = sum(1 for task in tasks if task[3])
    budget = error_reads_script.memory_budget()
    memory = ["--memory", f"{budget / error_runs / 1e9:.3f}"] if error_runs and budget != float("inf") else []
    # The ancestry report's genes (training_db's reference.fna, packed by its build) unpacked meanwhile, if not kept.
    unpacking = start_unpack(training_db, protal, work, threads, steps) if sams and protal and work else None
    say(f"Reports of what the models' errors rest on ({', '.join(os.path.basename(t[2]) for t in tasks)}), side by side, "
        f"on the {share} thread{'s' if share > 1 else ''}")
    if contig_cache and (trace or any(kind in ("pe", "se") for kind in scopes)):
        fastas = []
        for collection in (training, test):
            for manifest in glob.glob(os.path.join(collection or "", "points", "*", "sim", "manifest.tsv")):
                fastas += [row.get("fasta_path") for row in manifest_rows(os.path.dirname(manifest))]
        try:
            contigs = trace_script.genome_contigs(fastas, threads, contig_cache)
            say(f"    the contig names of {len(contigs)} genomes in {clock(time.time() - began)} ({contig_cache})")
        except Exception as e:  # noqa: BLE001: each report reads what it misses itself
            say(f"    reading the genomes' contig names failed ({e}); each report reads them itself")
    running = []
    try:
        for what, command, log, out in tasks:
            command = command + ["--threads", str(share)] + (memory if out else [])
            fh = open(log, "w")
            fh.write(" ".join(command) + "\n")
            fh.flush()
            running.append((what, subprocess.Popen(command, stdout=fh, stderr=subprocess.STDOUT), log, out, fh))
        looked = time.time()
        while any(p.poll() is None for _, p, _, _, _ in running):
            time.sleep(0.5)
            if Job.scratch and time.time() - looked >= 5:
                Job.scratch.look()
                looked = time.time()
    finally:
        for _, p, _, _, fh in running:
            if p.poll() is None:
                p.kill()
                p.wait()
            fh.close()
    for what, p, log, out, _ in running:
        if not out:  # trace_relatives.py
            if p.returncode:
                say(f"    tracing the held-out species' reads failed ({p.returncode}; see {log}); the build went on")
                continue
            with open(relatives + ".txt") as fh:
                text = fh.read()
            first = text.splitlines()[0] if text.startswith("Not traced") else "model_logs/relatives_by_gene_conservation.txt"
            say(f"    the held-out species' reads by gene conservation: {first}")
        elif p.returncode:
            told[what] = f"{what}: failed ({p.returncode}; see {log})"
        else:
            told[what] = f"{what} {error_reads_summary(os.path.join(out, 'summary.tsv'))}"
    if told:
        say(f"    the reads of the models' errors (model_logs/error_reads, logs/error_reads_<read type>.log): "
            + "; ".join(told[kind] for kind in scopes if kind in told))
    say(f"    reported in {clock(time.time() - began)}")
    done = [what for what, p, _, out, _ in running if out and not p.returncode]
    if sams and done and protal and taxonomy and work:
        ancestry_reports(done, logs, training_db, protal, taxonomy, heldout, work, threads, steps, unpacking)
    elif unpacking:
        unpacking.kill()
        unpacking.wait()
        shutil.rmtree(os.path.join(work, "ancestry_files"), ignore_errors=True)


def kept_reference(training_db):
    """The training database's reference.fna where it is still there: the hard link the run keeps for the ancestry
    report (ANCESTRY_REFERENCE), or the file beside database.protal (a database built with --no_bundle); None if packed."""
    return next((p for p in (os.path.join(training_db, ANCESTRY_REFERENCE), os.path.join(training_db, "reference.fna"),
                             os.path.join(training_db, "reference.fna.zst")) if os.path.isfile(p)), None)


def start_unpack(training_db, protal, work, threads, steps):
    """protal --unpack_db of the training database into work/ancestry_files, for the ancestry report's genes, in the
    background (its Popen; logs/ancestry_unpack.log), or None when its reference.fna is kept (kept_reference)."""
    if kept_reference(training_db):
        return None
    unpacked = os.path.join(work, "ancestry_files")
    shutil.rmtree(unpacked, ignore_errors=True)
    fh = open(os.path.join(steps, "ancestry_unpack.log"), "w")
    command = [protal, "--unpack_db", "--db", os.path.join(training_db, "database.protal"), "--unpack_dir", unpacked,
               "-t", str(threads)]
    fh.write(" ".join(command) + "\n")
    fh.flush()
    process = subprocess.Popen(command, stdout=fh, stderr=subprocess.STDOUT)
    fh.close()  # the child holds its own descriptor
    return process


def ancestry_reports(kinds, logs, training_db, protal, taxonomy, heldout, work, threads, steps, unpacking=None):
    """ancestry_sites.py on every read type's error-read SAMs (model_logs/error_reads/<kind>) in one run, which reads the
    references once for all of them and analyses the read types side by side on `threads` (each read type's outputs as
    a run of its own would write them): model_logs/ancestry_sites/<kind>.{summary.txt,auc.tsv,taxa.tsv.gz,fragments.tsv.gz}
    and logs/ancestry_sites.log. The training database's genes come from its reference.fna if that is still there
    (kept_reference), else from protal --unpack_db into work/ancestry_files (removed afterwards; `unpacking`, started by
    start_unpack, or now); the species' other genomes from its full reference if it was kept (keep_full). A failure is
    reported, and does not stop the build."""
    began = time.time()
    reference = kept_reference(training_db)
    unpacked = None
    if not reference:
        unpacked = os.path.join(work, "ancestry_files")
        log = os.path.join(steps, "ancestry_unpack.log")
        if unpacking is None:
            unpacking = start_unpack(training_db, protal, work, threads, steps)
        code = unpacking.wait()
        if code or not os.path.isfile(os.path.join(unpacked, "reference.fna")):
            say(f"    the ancestry sites of the errors' reads: unpacking {os.path.basename(training_db)} failed ({code}; "
                f"see {log}); the build went on")
            shutil.rmtree(unpacked, ignore_errors=True)
            return
        reference = os.path.join(unpacked, "reference.fna")
    full = full_reference_path(training_db)
    out_dir = os.path.join(logs, "ancestry_sites")
    os.makedirs(out_dir, exist_ok=True)
    told = {}
    with_sams = []
    for kind in kinds:
        if glob.glob(os.path.join(logs, "error_reads", kind, "*", "*", "*.sam*")):
            with_sams.append(kind)
        else:
            told[kind] = f"{kind}: no SAMs"  # no errors, or none with a read
    if with_sams:
        command = [sys.executable, ANCESTRY, "--sams", *(os.path.join(logs, "error_reads", kind) for kind in with_sams),
                   "--out", *(os.path.join(out_dir, kind) for kind in with_sams), "--threads", str(max(1, threads)),
                   "--reference", reference, "--taxonomy", taxonomy]
        command += ["--heldout", heldout] if heldout and os.path.isfile(heldout) else []
        command += ["--full-reference", full] if full else []
        command += ["--allele-genome-share", f"{ALLELE_SHARE:g}"] if full and ALLELE_SHARE < 1 else []
        log = os.path.join(steps, "ancestry_sites.log")
        with open(log, "w") as fh:
            fh.write(" ".join(command) + "\n")
            fh.flush()
            code = subprocess.run(command, stdout=fh, stderr=subprocess.STDOUT).returncode
        for kind in with_sams:
            told[kind] = (f"{kind}: failed ({code}; see {log})" if code
                          else f"{kind}: {ancestry_summary(os.path.join(out_dir, kind + '.auc.tsv'))}")
    if unpacked:
        shutil.rmtree(unpacked, ignore_errors=True)
    told = [told[kind] for kind in kinds]
    say(f"    the ancestry sites of the errors' reads (model_logs/ancestry_sites, logs/ancestry_sites.log; "
        + (f"the species' alleles from {os.path.basename(full)}" if full else "the congener sites only, no full reference")
        + f"): {'; '.join(told)}; in {clock(time.time() - began)}")


def ancestry_summary(path):
    """ancestry_sites.py's auc.tsv in words: over the taxa with 10 or more sites, the AUC of plain identity and of the
    sites' signals for the false negatives' own reads against the false positives'."""
    if not os.path.isfile(path):
        return "no taxa to compare"
    rows = []
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["min_sites"] == "10" and row["identity_band"] == "all":
                rows.append(row)
    if not rows:
        return "fewer than 10 sites on every taxon, or taxa of one kind only"
    short = {"identity": "identity", "species_base_at_congener_sites": "species' base at the congener sites",
             "species_base_at_fixed_sites": "at the fixed sites", "fixed_site_identity": "fixed-site identity"}
    parts = [f"{short[r['signal']]} {float(r['auc']):.3f}" for r in rows if r["signal"] in short]
    return f"{rows[0]['taxa']} taxa ({rows[0]['fn']} FN) with 10 or more sites, AUC " + ", ".join(parts)


def error_reads_summary(path):
    """The samples, error taxa, fragments and size of error_reads.py's summary.tsv, in words."""
    if not os.path.isfile(path):
        return "no samples"
    c = collections.Counter()
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            c.update({"samples": 1, row["set"]: 1, "FP": int(row["FP"]), "FN": int(row["FN"]),
                      "unseen": int(row["unseen"]), "fragments": int(row["fragments"]),
                      "kept": int(row.get("sam_fragments") or 0), "bytes": int(row["sam_bytes"])})
    if not c["samples"]:
        return "no samples"
    return (f"{c['samples']} samples ({c['training']} training, {c['test']} test): {c['FP']} FP, {c['FN']} FN and "
            f"{c['unseen']} unseen taxa, {c['fragments']} fragments"
            + (f" ({c['kept']} in the SAMs, {gigabytes(c['bytes'])})" if c["kept"] else ""))


# A number of the tables with more than 9 significant digits (Python's repr of a double). 9 keep every float32 (what a
# forest compares, model_pmml.py) and are far beyond the features' precision, at three quarters of the text.
LONG_NUMBER = re.compile(rb"(?<![\w.+-])-?\d+\.\d{9,}(?:e[-+]?\d+)?(?![\w.])")


def shorten_table(job):
    """Copies a table with its numbers to 9 significant digits (LONG_NUMBER). -> its bytes before and after."""
    source, target = job
    def nine(m):
        return format(float(m.group()), ".9g").encode()
    with open(source, "rb") as fh, open(target, "wb") as out:
        for line in fh:
            out.write(LONG_NUMBER.sub(nine, line))
    return os.path.getsize(source), os.path.getsize(target)


GATHERED_MAX = 1 << 20  # a log gathered into the share archive longer than this keeps its first and last GATHERED_KEEP
GATHERED_KEEP = 1 << 18  # bytes


def natural_key(text):
    """Sorts names with numbers by their numbers (p2 before p10)."""
    return [int(x) if x.isdigit() else x for x in re.split(r"(\d+)", text)]


def concatenate_logs(paths, target, label):
    """Writes the files one after the other into target, each under a "==> label(path) <==" line; a file longer than
    GATHERED_MAX keeps its first and last GATHERED_KEEP bytes."""
    with open(target, "wb") as out:
        for path in paths:
            out.write(f"==> {label(path)} <==\n".encode())
            size = os.path.getsize(path)
            with open(path, "rb") as fh:
                if size <= GATHERED_MAX:
                    text = fh.read()
                else:
                    text = fh.read(GATHERED_KEEP) + f"\n[... {size - 2 * GATHERED_KEEP} bytes left out ...]\n".encode()
                    fh.seek(size - GATHERED_KEEP)
                    text += fh.read()
            out.write(text if text.endswith(b"\n") or not text else text + b"\n")


def share_archive(outdir, folders, threads=1, insilico_table=None):
    """--share-logs: OUTDIR/<name>_share.tar.gz, all under <name>/ where OUTDIR has them: console.log, logs/, model_logs/
    (the reports, predictions, calls and the error reads' tables and SAMs; not the models, which are in the database),
    protal_db/build_metadata.tsv, work/internal_taxonomy.dmp and each collection's tables (work/training/, work/test/;
    shortened by shorten_table). And from the samples' disk, which the end of a cluster job may clear: per collection
    (folders: [(name, its collect_training_data.py -o)]) its design points' simulator logs and parameters (simulate.log,
    design.log, stream*.log, run_params.tsv) as logs/simulations_<collection>.log, its protal runs' logs as
    logs/protal_runs_<collection>.log where the build did not gather them there (without --scratch), and
    insilico_strains.tsv (insilico_strains.py's per-strain table) as model_logs/insilico_strains.tsv. A failure is
    reported, and does not stop the build."""
    name = os.path.basename(os.path.normpath(outdir))
    target = os.path.join(outdir, name + "_share.tar.gz")
    began = time.time()
    shortened = os.path.join(outdir, WORK, "share_tables")
    try:
        shutil.rmtree(shortened, ignore_errors=True)
        jobs = []
        for which, folder in folders:
            for table in TABLES.values():
                if folder and os.path.isfile(os.path.join(folder, table)):
                    os.makedirs(os.path.join(shortened, which), exist_ok=True)
                    jobs.append((os.path.join(folder, table), os.path.join(shortened, which, table)))
        sizes = []
        if jobs:
            with concurrent.futures.ProcessPoolExecutor(max(1, min(threads, len(jobs)))) as pool:
                sizes = list(pool.map(shorten_table, jobs))
        gathered, simulator_logs = [], 0
        os.makedirs(shortened, exist_ok=True)
        for which, folder in folders:
            if not folder:
                continue
            points = os.path.join(folder, "points")
            logs = sorted(glob.glob(os.path.join(points, "*", "*.log")) +
                          glob.glob(os.path.join(points, "*", "*", "run_params.tsv")),
                          key=lambda p: natural_key(os.path.relpath(p, points)))
            if logs:
                path = os.path.join(shortened, f"simulations_{which}.log")
                concatenate_logs(logs, path, lambda p: os.path.relpath(p, points))
                gathered.append((path, f"{LOGS}/simulations_{which}.log"))
                simulator_logs += len(logs)
            cpu = os.path.join(folder, PROTAL_CPU)  # where the build did not copy it to logs/ (a stopped collection)
            if os.path.isfile(cpu) and not os.path.isfile(os.path.join(outdir, LOGS, f"protal_cpu_{which}.tsv")):
                gathered.append((cpu, f"{LOGS}/protal_cpu_{which}.tsv"))
            if not os.path.isfile(os.path.join(outdir, LOGS, f"protal_runs_{which}.log")):
                profile_all = os.path.join(folder, "profile_all")
                runs = sorted(glob.glob(os.path.join(profile_all, "**", "protal.log"), recursive=True),
                              key=lambda p: natural_key(os.path.relpath(p, profile_all)))
                if runs:
                    path = os.path.join(shortened, f"protal_runs_{which}.log")
                    concatenate_logs(runs, path, lambda p: os.path.relpath(os.path.dirname(p), profile_all))
                    gathered.append((path, f"{LOGS}/protal_runs_{which}.log"))
        if insilico_table and os.path.isfile(insilico_table):
            gathered.append((insilico_table, f"{REPORTS}/insilico_strains.tsv"))
        left_out = (".partial", ".xml", ".joblib")
        with tarfile.open(target + ".partial", "w:gz", compresslevel=6) as tar:
            for path in ("console.log", LOGS, REPORTS, os.path.join(DATABASE, "build_metadata.tsv"),
                         os.path.join(WORK, "internal_taxonomy.dmp")):
                if os.path.exists(os.path.join(outdir, path)):
                    tar.add(os.path.join(outdir, path), f"{name}/{path}",
                            filter=lambda t: None if t.name.endswith(left_out) else t)
            for _, table in jobs:
                tar.add(table, f"{name}/{WORK}/{os.path.relpath(table, shortened)}")
            for source, path in gathered:
                tar.add(source, f"{name}/{path}")
        os.replace(target + ".partial", target)
        before, after = sum(s[0] for s in sizes), sum(s[1] for s in sizes)
        say(f"Logs to share: {target}, {gigabytes(os.path.getsize(target))} in {clock(time.time() - began)} (the logs, "
            f"model_logs/ and {len(jobs)} tables, {gigabytes(before)} shortened to {gigabytes(after)}; from the samples' "
            f"disk {simulator_logs} simulator logs" + (" and the in-silico strains' table" if insilico_table and
                                                       os.path.isfile(insilico_table) else "") +
            f"); unpack with tar xzf {os.path.basename(target)}")
    except Exception as e:  # noqa: BLE001: the database is ready either way
        say(f"Packing the logs to share failed ({e}); the database is ready either way")
    finally:
        shutil.rmtree(shortened, ignore_errors=True)


def gene_neighbours_summary(build_log, what="Gene neighbours:"):
    """What protal --build said of the gene neighbours it packed (its "Gene neighbours:" line, or that of
    `what`: "Gene positions:"), for build_metadata.tsv."""
    try:
        with open(build_log) as fh:
            lines = [line.strip() for line in fh if line.startswith(what)]
    except OSError:
        return "unknown (no build log)"
    if not lines:
        return "none (this protal does not pack them)"
    return lines[-1].removeprefix(what).strip()


def genes_placed(log):
    """gene_neighbours.py's count of the genes it placed, from its log: "in N genomes (X exactly, Y by their
    k-mer trace)", or "" if the log has none."""
    try:
        with open(log) as fh:
            m = re.search(r"Genes placed: (\d+) exactly and (\d+) by their k-mer trace in the (\d+) genomes used", fh.read())
    except OSError:
        return ""
    return f" in {m.group(3)} genomes ({m.group(1)} exactly, {m.group(2)} by their k-mer trace)" if m else ""


def tool_versions(args):
    """protal's and the simulator's --version (their last line) and the scripts' commit ("(scripts changed since)" when
    the checkout's scripts differ from it), as build_metadata.tsv records them."""
    def output(command):
        try:
            return subprocess.run(command, capture_output=True, text=True, timeout=60).stdout.strip()
        except (OSError, subprocess.SubprocessError):
            return ""
    version = output([args.protal, "--version"]).splitlines()
    simulator_version = output([args.simulator, "--version"]).splitlines()
    commit = output(["git", "-C", HERE, "rev-parse", "HEAD"])
    if commit and output(["git", "-C", HERE, "status", "--porcelain", "--", "."]):
        commit += " (scripts changed since)"
    return {"protal_version": version[-1] if version else "unknown",
            "simulator_version": simulator_version[-1] if simulator_version else "unknown (older than 0.7.3)",
            "scripts_commit": commit or "unknown (not a git checkout)"}


def versions_at_end(started, ended):
    """build_metadata.tsv's versions: those the run started with (tool_versions at its start), each followed by what it
    was at the end where that differs, and the names of those that changed. A checkout pulled or rebuilt during a run of
    hours runs its later steps (the trainer, protal's runs) with the new scripts and binaries; read only at the end, as
    until 2026-10-06, the metadata named the new commit for a run that began with an older one (r226 v14, which trained
    500 rounds at 0.05 under 15006c2's commit; docs/claude/2026-10-06-r226-v14)."""
    versions, changed = {}, []
    for key, value in started.items():
        if ended.get(key, value) != value:
            versions[key] = f"{value}; at the end of the run: {ended[key]}"
            changed.append(key)
        else:
            versions[key] = value
    return versions, changed


def provenance(args, versions, release, genome_table, heldout, n_heldout, read_types, prefixes, genes="all",
               insilico="none"):
    """build_metadata.tsv: what the database was built from and with, so that two builds can be compared; versions:
    versions_at_end's."""
    genomes, species = 0, set()
    with open(genome_table) as fh:
        for line in fh:
            lineage = next((f for f in line.rstrip("\n").split("\t") if f.startswith("d__") and ";s__" in f), None)
            if lineage:
                genomes += 1
                species.add(lineage.split(";")[-1])
    clade_counts = collections.Counter(rank for rank, _ in set(read_holdout(heldout).values())) if n_heldout else {}
    rows = [("gtdb_release", f"r{release}"), ("built", time.strftime("%Y-%m-%d %H:%M:%S")),
            ("protal_version", versions["protal_version"]), ("protal_binary", args.protal),
            ("simulator_version", versions["simulator_version"]),
            ("scripts_commit", versions["scripts_commit"]), ("command", " ".join(sys.argv)),
            ("seed", args.seed), ("genome_table", f"{genomes} genomes of {len(species)} species"),
            ("insilico_strains", insilico),
            ("samples_streamed", getattr(args, "streaming", "")),
            ("marker_genes", genes),
            ("gene_conservation", gene_conservation_summary(os.path.join(args.outdir, LOGS, "index_and_package.log"))),
            ("gene_neighbours", gene_neighbours_summary(os.path.join(args.outdir, LOGS, "index_and_package.log"))),
            ("gene_positions", gene_neighbours_summary(os.path.join(args.outdir, LOGS, "index_and_package.log"),
                                                       "Gene positions:")),
            ("gene_congeners", gene_congeners_summary(os.path.join(args.outdir, LOGS, "index_and_package.log"))),
            ("suspect_copies", suspect_copies_summary(os.path.join(args.outdir, LOGS, "index_and_package.log"))),
            ("classifier_features", args.features), ("classifier_model", args.model),
            ("classifier_trees", args.ntree if args.model == "forest" else f"{args.rounds or TRAINER_ROUNDS} rounds"),
            ("classifier_max_leaves", ",".join(f"{t}:{max_leaves(args.maxnodes, t)}" for t in read_types)),
            ("classifier_evaluation", args.evaluation),
            ("classifier_previous_procedure",
             "compared" if args.previous_procedure and args.evaluation != "none" else "not compared"),
            ("classifier_training_species_left_out", n_heldout)]
    rows += [(f"classifier_training_{rank}_clades_left_out", clade_counts[rank]) for rank in CLADE_RANKS
             if clade_counts.get(rank)]
    rows += [("classifier_training_holdout_clouds", getattr(args, "holdout_clouds_note", "none"))]
    rows += [("classifier_training_samples", f"{args.samples} per design point"),
             ("classifier_training_design", f"read pairs {args.read_pairs}; read setups {args.read_setups}; "
                                            f"species per sample {args.species_per_sample}; strains "
                                            f"{args.strains_per_species or 'one'}; abundance {args.abundance or 'default'}; "
                                            f"congeners {congener_text(args.congeners)}"),
             ("classifier_read_types", ",".join(read_types)),
             ("classifier_depth_knobs", ",".join(t for t in read_types if t in depth_knob_types(args)) or "none"),
             ("classifier_call_mode", args.call_mode),
             ("classifier_scenarios", (", ".join(f"{name} {n_in} hold-in (training) and {n_out} hold-out (test) samples"
                                                 for name, (n_in, n_out) in args.scenario_samples_of.items())
                                       + f"; their rows weighted {args.scenario_weight:g}"
                                       + (f"; definitions {args.scenario_file}" if args.scenario_file else "")
                                       + (f"; host genome {args.host_genome}" if args.host_genome else "")
                                       + "".join(f"; {note}" for note in getattr(args, "scenario_notes", [])))
              if getattr(args, "scenario_samples_of", None) else
              "none" + "".join(f"; {note}" for note in getattr(args, "scenario_notes", [])))]
    for t in read_types:
        try:
            with open(prefixes[t] + ".metrics.json") as fh:
                metrics = json.load(fh)
        except (OSError, ValueError):
            continue
        if metrics.get("features_auto"):
            rows.append((f"model_{t}_features", f"{metrics['features_auto']['chosen']} (--features {args.features}: "
                                                f"{metrics['features_auto'].get('why', '')})"))
        scenario_rows = [m for m in metrics.get("scenarios", []) if m["scenario"] != "(design)" and
                         m["set"] in ("hold-out", "hold-in, species held out")]
        if scenario_rows:
            num = lambda v, unit="": "-" if v is None else f"{v:.4g}{unit}"
            rows.append((f"model_{t}_scenarios", "; ".join(
                f"{m['scenario']} {m['set']} F1 {num(m['F1'])}, FP rate {num(m['FP rate %'], '%')}, FN rate "
                f"{num(m['FN rate %'], '%')}" for m in scenario_rows)))
        species_cv = metrics.get("evaluation", {}).get("species", {})
        test = metrics.get("test", {}).get("this one", {})
        rows.append((f"model_{t}", f"species held out F1 {species_cv.get('F1')}, FP per sample "
                                   f"{species_cv.get('FP_per_sample')}; independent test F1 {test.get('F1')}, FP per "
                                   f"sample {test.get('FP_per_sample')}"))
        knobs = metrics.get("depth_knobs", {})
        curve = knobs.get("curve")
        if knobs.get("global_knob") is not None:  # one knob for every sample (a curve of one point)
            chosen = metrics.get("global_knob", {})
            rows.append((f"model_{t}_depth_knobs", f"none: {knobs.get('skipped')}; knob {knobs['global_knob']:g} for every "
                                                   f"sample, F1 {chosen.get('gain', 0):.4f} above 0.5 with species held out; "
                                                   f"independent test F1 at it {metrics.get('test_depth_knobs', {}).get('F1')}"))
        elif curve:
            rows.append((f"model_{t}_depth_knobs", "log10 fragments:knob " + ",".join(f"{x:.3f}:{k:g}" for x, k in curve) +
                         f"; independent test F1 at them {metrics.get('test_depth_knobs', {}).get('F1')}"))
        elif knobs.get("skipped"):
            rows.append((f"model_{t}_depth_knobs", "none: " + knobs["skipped"] + " (protal calls at --knob: no knob gained "
                                                   "0.002 with species held out)"))
        false_calls = metrics.get("false_calls")
        if false_calls:
            rows.append((f"model_{t}_false_calls", f"target {false_calls['fdr']} per sample; species held out F1 "
                                                   f"{false_calls['F1']}; independent test F1 at it "
                                                   f"{metrics.get('test_false_calls', {}).get('F1')}"))
    return rows


def model_scores(read_types, prefixes):
    """Each model's F1 at knob 0.5 with species held out in training and on the independent test set as protal calls
    (at the model's own knob or curve, if it has one), from its .metrics.json, in a few words."""
    scores, with_test = [], False
    for t in read_types:
        try:
            with open(prefixes[t] + ".metrics.json") as fh:
                metrics = json.load(fh)
        except (OSError, ValueError):
            scores.append(f"{t} -")
            continue
        cv = metrics.get("evaluation", {}).get("species", {}).get("F1")
        called = metrics.get("test_depth_knobs")
        test = (called or metrics.get("test", {}).get("this one", {})).get("F1")
        with_test |= test is not None
        scores.append(f"{t} " + ("-" if cv is None else f"{cv:.3f}") + ("" if test is None else f"/{test:.3f}") +
                      (f" ({called['knob']})" if called and called.get("knob") else ""))
    return f"F1 with species held out (knob 0.5){'/on the test set (as protal calls)' if with_test else ''}: {', '.join(scores)}"


def feature_choices(read_types, prefixes):
    """For each model whose trainer chose its feature set (--features auto), which set won and why (its F1 with
    species held out against the other sets', the rule, and the sets' F1 on the test sets), from its .metrics.json:
    [(read type, set, why)]."""
    out = []
    for t in read_types:
        try:
            with open(prefixes[t] + ".metrics.json") as fh:
                auto = json.load(fh).get("features_auto")
        except (OSError, ValueError):
            continue
        if auto:
            out.append((t, auto["chosen"], auto.get("why", "")))
    return out


def scenario_scores(read_types, prefixes, order=()):
    """The models' F1 per scenario on its hold-out samples (and its hold-in ones with species held out), from their
    .metrics.json (machine_learning_cmdline.py, section Scenarios), in a few words; the scenarios in `order` first, in
    that order."""
    found = collections.defaultdict(list)
    for t in read_types:
        try:
            with open(prefixes[t] + ".metrics.json") as fh:
                rows = json.load(fh).get("scenarios", [])
        except (OSError, ValueError):
            continue
        sets = {(r["scenario"], r["set"]): r for r in rows}
        for name in dict.fromkeys(r["scenario"] for r in rows if r["scenario"] != "(design)"):
            held, cv = sets.get((name, "hold-out")), sets.get((name, "hold-in, species held out"))
            f1 = lambda r: "-" if not r or r.get("F1") is None else f"{r['F1']:.3f}"
            found[name].append(f"{t} {f1(held)}/{f1(cv)}")
    names = [n for n in order if n in found] + [n for n in found if n not in order]
    return ("scenario F1 on the hold-out samples/hold-in with species held out: " +
            "; ".join(f"{name} {', '.join(found[name])}" for name in names)) if found else "no scenario scores"


def summary_lines(read_types, prefixes, db):
    """model_logs/summary.txt: for each read type's model, TP, FP, TN, FN and the rates, with species held
    out in training (cross-validation), on the independent test set and, with --scenarios, on each scenario's
    hold-out and hold-in samples, from its .metrics.json."""
    header = ("read type", "evaluated on", "knob", "taxa", "TP", "FP", "TN", "FN", "sensitivity", "specificity",
              "precision", "F1", "FP rate", "FN rate", "FP/sample")
    rows = []
    rate = lambda v: "-" if v is None else f"{v:.4f}"

    def row(t, label, knob, tp, fp, tn, fn, per_sample):
        return (t, label, knob, str(tp + fp + tn + fn), str(tp), str(fp), str(tn), str(fn),
                rate(tp / (tp + fn) if tp + fn else None), rate(tn / (tn + fp) if tn + fp else None),
                rate(tp / (tp + fp) if tp + fp else None), rate(2 * tp / (2 * tp + fp + fn) if tp + fp + fn else None),
                rate(fp / (fp + tn) if fp + tn else None), rate(fn / (tp + fn) if tp + fn else None),
                "-" if per_sample is None else f"{per_sample:.2f}")
    for t in read_types:
        try:
            with open(prefixes[t] + ".metrics.json") as fh:
                metrics = json.load(fh)
        except (OSError, ValueError):
            rows.append((t, "no metrics (training failed?)") + ("",) * (len(header) - 2))
            continue
        # How protal calls with the model: at its own knob or curve (the trainer's --depth-knobs), else at 0.5.
        called = metrics.get("test_depth_knobs") or {}
        label = called.get("knob") or ""
        own = "curve" if label.startswith("knob curve") else (label.removeprefix("knob ") or "0.5")
        species = metrics.get("evaluation", {}).get("species")
        if species:
            fp = species["FP"]
            rows.append(row(t, "species held out", "0.5", species["present"] - species["FN"], fp,
                            species["taxa"] - species["present"] - fp, species["FN"], species.get("FP_per_sample")))
        test = called if "TP" in called else metrics.get("test", {}).get("this one")
        if test:
            fp = test["FP"]
            rows.append(row(t, "independent test set", own if test is called else "0.5", test["present"] - test["FN"], fp,
                            test["taxa"] - test["present"] - fp, test["FN"], test.get("FP_per_sample")))
        # The scenarios: their hold-out samples, and their hold-in ones with species held out and in sample, as protal
        # calls (the trainer scores them at the model's knob).
        for m in metrics.get("scenarios", []):
            if m["scenario"] != "(design)" and m["set"] in ("hold-out", "hold-in, species held out", "hold-in, in sample"):
                rows.append(row(t, f"{m['scenario']}: {m['set']}", own, m["TP"], m["FP"], m["absent"] - m["FP"], m["FN"],
                                m["FP/sample"]))
    widths = [max(len(str(r[i])) for r in [header, *rows]) for i in range(len(header))]
    table = ["  ".join(str(v).ljust(w) for v, w in zip(r, widths)).rstrip() for r in [header, *rows]]
    choices = feature_choices(read_types, prefixes)
    if choices:  # --features auto: the set each trainer chose, and why
        table += ["", "Feature sets chosen (--features auto):"] + [f"  {t}: {name}: {why}" for t, name, why in choices]
    return [f"Presence models of {db}, taxa scored at the knob given (0.5, or the model's own, as protal calls by "
            "default): TP present and called, FP absent and called, TN "
            "absent and not called, FN present and not called; FP rate FP / (FP + TN), FN rate FN / (TP + FN). "
            "species held out: each taxon scored by models that did not see its species; independent test set: "
            "samples of another design, scored by the final model; <scenario>: hold-out: the scenario's samples of "
            "the test set (never trained on), hold-in: its training samples, scored with species held out, or in "
            "sample (by the final model, fitted on them). Details: trained_model*.report.txt", ""] + table


def composition_report(read_types, prefixes, training, test, heldout, logs, threads=1):
    """model_logs/composition_accuracy[_<read type>].tsv and composition_accuracy.txt (composition_accuracy.py): how well
    the samples' composition (the share of their reads the called species explain, the profile's unknown share "?", the
    average genome size, the species' genome sizes and depths) matches the simulator's truth, with the trained models'
    calls (PREFIX.calls.tsv.gz); composition_species[_<read type>].tsv.gz, every sample's species present or called with
    protal's depth and genome size against the truth (the profiles and manifests it rests on stay on the samples' disk).
    -> ({read type: one line of medians}, the summary's lines). A failure is reported, and does not stop the build."""
    import composition_accuracy
    briefs, text = {}, []
    for t in read_types:
        calls = prefixes[t] + ".calls.tsv.gz"
        if not os.path.isfile(calls):
            briefs[t] = f"no {os.path.basename(calls)}"
            continue
        suffix = "" if t == "pe" else "_" + t
        argv = ["--calls", calls, "--training", training, "--read-type", t, "--threads", str(threads),
                "--out", os.path.join(logs, f"composition_accuracy{suffix}.tsv"),
                "--species-out", os.path.join(logs, f"composition_species{suffix}.tsv.gz")]
        argv += ["--test", test] if test else []
        argv += ["--heldout", heldout] if heldout else []
        output = io.StringIO()
        try:
            with contextlib.redirect_stdout(output):
                lines = composition_accuracy.main(argv)
            briefs[t] = composition_accuracy.brief(lines)
        except BaseException as e:  # SystemExit too: the build goes on
            briefs[t] = f"failed: {e}"
        text += output.getvalue().splitlines() + [""]
    with open(os.path.join(logs, "composition_accuracy.txt"), "w") as fh:
        fh.write("\n".join(text) + "\n")
    return briefs, text


def remove_full_reference(folder):
    """Removes a folder's full reference (full_reference.fna.zst) once its build is done: the marker genes of
    every genome (86 GB raw at r226), which only that build reads. A rebuild converts the release again (the
    build packed reference.fna into database.protal). Says what it removed ("; ..."), or nothing."""
    path = full_reference_path(folder)
    if not path:
        return ""
    return f"; {os.path.basename(path)} removed ({gigabytes(remove_full_reference_files(folder))})"


def full_reference_fate(folder, keep, why="the ancestry report"):
    """remove_full_reference, or with keep the note that the full reference stays for `why`: the foreign scan
    (foreign_rates) reads its marker genes right after the build, the ancestry report (reports) the species' other
    genomes once the error reads are taken."""
    if not keep:
        return remove_full_reference(folder)
    path = full_reference_path(folder)
    return f"; {os.path.basename(path)} kept for {why} ({gigabytes(os.path.getsize(path))})" if path else ""


# The strain alleles' options of every protal --build of the run (--strain-alleles, --allele-genome-share,
# --index-alleles, off by default), set by main():
# both databases take their alleles from the same genomes, none the simulations draw strains from.
ALLELE_ARGS = []
# The share of the genomes the strain alleles come from (1: all), for the ancestry report's alleles (ancestry_reports).
ALLELE_SHARE = 1.0


def build_command(protal, db, threads, *extra):
    command = [protal, "--build", "--no_profile", "-t", str(threads), "--db", db,
               "--reference", os.path.join(db, "reference.fna"), *extra]
    if full_reference_path(db):
        command += ["--full_reference", full_reference_path(db), *ALLELE_ARGS]
    return command


def split_allele_genomes(table, reps, share, output):
    """The genome table to simulate from, without the genomes that give strain alleles: of each species' genomes other
    than its representative, those allele_genome() picks at `share` (protal --build --allele_genome_share takes its
    alleles from them), so that no simulated strain is its species' own allele. Without the representatives known
    (reps None) every genome is judged. Writes the kept rows to `output` in their order; -> (the table, rows removed).
    `table` itself, unchanged, if none is removed."""
    with open(table) as fh:
        lines = fh.readlines()
    kept = [line for line in lines if not is_allele_genome(line.rstrip("\n").split("\t"), reps, share)]
    if len(kept) == len(lines):
        return table, 0
    with open(output + ".partial", "w") as fh:
        fh.writelines(kept)
    os.replace(output + ".partial", output)
    return output, len(lines) - len(kept)


def is_allele_genome(fields, reps, share):
    """Whether a genome table's row (its fields) is a genome that gives strain alleles (split_allele_genomes)."""
    accession = normalize_accession(fields[0]) if len(fields) >= 3 and fields[0] else ""
    return bool(accession) and (reps is None or accession not in reps) and allele_genome(accession, share)


def stop_jobs():
    """Stops the jobs still running; one that already ended well still counts (its build is recorded)."""
    for job in list(Job.running):
        if job.poll() == 0:
            job.ended(0)
    for job in list(Job.running):
        job.kill()
    Job.running.clear()


def stop(message):
    stop_jobs()
    sys.exit(message)


def on_signal(signum, _frame):
    stop(f"Stopped by {signal.Signals(signum).name}, with the commands it ran")


# ---- resuming --------------------------------------------------------------------------------------------
# A rerun into the same OUTDIR (after a failure in training, say) skips the conversion and the two index
# builds (at GTDB scale an hour or two each) when their inputs are those of the run that completed them:
# OUTDIR/.stages/<stage>.json holds the key of those inputs. The collector resumes on its own.

class Stages:
    def __init__(self, folder):
        self.folder = folder
        os.makedirs(folder, exist_ok=True)

    def done(self, name, key):
        return self.stored(name) == json.loads(json.dumps(key))

    def stored(self, name):
        """The key a stage was marked with, None if it was not."""
        try:
            with open(os.path.join(self.folder, name + ".json")) as fh:
                return json.load(fh)
        except (OSError, ValueError):
            return None

    def mark(self, name, key):
        path = os.path.join(self.folder, name + ".json")
        with open(path + ".partial", "w") as fh:
            json.dump(key, fh, indent=1)
        os.replace(path + ".partial", path)

    def forget(self, name):
        if os.path.exists(os.path.join(self.folder, name + ".json")):
            os.remove(os.path.join(self.folder, name + ".json"))


def file_identity(path):
    """A file as part of a key: its real path, size and modification time."""
    real = os.path.realpath(shutil.which(path) or path)
    st = os.stat(real)
    return [real, st.st_size, st.st_mtime_ns]


def content_hash(path):
    with open(path, "rb") as fh:
        return hashlib.sha1(fh.read()).hexdigest()


def copy_file(source, target):
    """Copies a file given on the command line to where the run keeps it, unless it is that file already (e.g.
    --holdout-species OUTDIR/model_logs/heldout_species.txt of a rerun)."""
    if not (os.path.exists(target) and os.path.samefile(source, target)):
        shutil.copyfile(source, target)


def release_identity(gtdb, release):
    """The release's files the converter reads, as part of a key: the taxonomy and metadata files, and the
    marker gene folders (whose modification time changes with the files in them)."""
    items = []
    for folder in (gtdb, os.path.join(gtdb, "genomic_files_reps"), os.path.join(gtdb, "genomic_files_all"),
                   os.path.join(gtdb, "auxillary_files")):
        for name in sorted(os.listdir(folder)) if os.path.isdir(folder) else ():
            path = os.path.join(folder, name)
            if f"_r{release}" in name and (os.path.isfile(path) or "_marker_genes_" in name):
                st = os.stat(path)
                items.append([os.path.relpath(path, gtdb), st.st_size if os.path.isfile(path) else len(os.listdir(path)),
                              st.st_mtime_ns])
    return items


def source_version(source=SOURCE):
    """The protal version of the source these scripts are part of (CMakeLists.txt), or None outside a checkout."""
    try:
        with open(os.path.join(source, "CMakeLists.txt")) as fh:
            found = re.search(r"project\(protal VERSION ([0-9.]+)\)", fh.read())
    except OSError:
        return None
    return found.group(1) if found else None


def pbsim_model_found(name, models, pbsim):
    """Whether pbsim3's model `name` (a file, or NAME[.model] in `models` or pbsim3's data folder) is found, as
    collect_training_data.pbsim_model finds it."""
    if os.path.isfile(name):
        return True
    folders = [models] if models else []
    exe = shutil.which(pbsim)
    if exe:
        prefix = os.path.dirname(os.path.dirname(os.path.realpath(exe)))
        folders += sorted(glob.glob(os.path.join(prefix, "share", "pbsim*", "data"))) + \
            sorted(glob.glob(os.path.join(prefix, "share", "pbsim*"))) + [os.path.join(prefix, "data")]
    return any(os.path.isfile(os.path.join(folder, candidate)) for folder in folders for candidate in (name, name + ".model"))


def build_check(command, name, source=SOURCE):
    """Whether the binary `command` (protal or the simulator) was built from the source these scripts are at, from
    its --version (the version, and since 0.7.3 the commit it was built from): (problem, note), either of them None.
    A binary of another version, of a commit the checkout lacks, or of a commit whose source has changed since is a
    problem: with it the run would fail hours in (the r226 run of 2026-10-02 collected its training data with an
    older protal, whose dump lacked the trainer's features). One that cannot say is only noted."""
    rebuild = f"build it from {source} (docs/installation.md)"
    try:
        result = subprocess.run([command, "--version"], capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.SubprocessError) as e:
        return f"{name} ({command}) --version failed ({e}): {rebuild}", None
    lines = result.stdout.strip().splitlines()
    found = re.fullmatch(r"\S+ v([0-9.]+)(?: \(commit ([0-9a-f]{40})(, with uncommitted changes)?\))?", lines[-1]) \
        if lines else None
    if result.returncode or not found:
        said = (result.stdout + result.stderr).strip()[-200:]
        return f"{name} ({command}) --version gives no version ({result.returncode}: {said}), so it is older than " \
               f"these scripts: {rebuild}", None
    version, commit, modified = found.groups()
    wanted = source_version(source)
    if wanted and version != wanted:
        return f"{name} ({command}) is v{version}, these scripts v{wanted}: {rebuild}", None
    if not commit:
        return None, f"{name} ({command}) does not say which commit it was built from (not built in a git checkout), " \
                     "so it is not checked against the source"
    if modified:
        return None, f"{name} ({command}) was built with uncommitted changes to its source, so it is not checked " \
                     "against the source"

    def git(*arguments):
        return subprocess.run(["git", "-C", source, *arguments], capture_output=True, text=True)
    try:
        checkout = git("rev-parse", "--git-dir").returncode == 0
    except OSError:
        checkout = False
    if not checkout:
        return None, f"{source} is no git checkout (or there is no git), so {name} is not checked against the source"
    if git("cat-file", "-e", commit + "^{commit}").returncode:
        return f"{name} ({command}) was built from commit {commit[:10]}, which {source} does not have: {rebuild}", None
    changed = git("diff", "--quiet", commit, "--", *BUILD_SOURCES).returncode
    if changed == 1:
        commits = git("rev-list", "--count", f"{commit}..HEAD", "--", *BUILD_SOURCES).stdout.strip() or "?"
        since = [f"{commits} commit{'' if commits == '1' else 's'}"] if commits != "0" else []
        if git("diff", "--quiet", "HEAD", "--", *BUILD_SOURCES).returncode == 1:
            since.append("uncommitted changes")
        detail = " and ".join(since) if since else "the checkout is at another commit"
        return f"{name} ({command}) was built from commit {commit[:10]}, and src/, lib/ or the build files have " \
               f"changed since ({detail}): rebuild it (cmake --build <build dir> --target protal " \
               "simulate_metagenomes)", None
    if changed:
        return None, f"git cannot compare {source} with commit {commit[:10]}, so {name} is not checked against it"
    return None, None


def check_tools(args, read_types):
    """What the run needs later, checked before it starts: a missing tool, or protal or the simulator of an older
    source, would otherwise stop it after the conversion, the builds or the collection. The trainer's imports (a few
    seconds from an environment on a network file system) are checked beside the binaries."""
    problems, notes = [], []
    importing = subprocess.Popen([sys.executable, "-c", "import joblib, numpy, pandas, sklearn"], stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE, text=True)

    def executable(command, what, hint):
        if shutil.which(command) or (os.path.isfile(command) and os.access(command, os.X_OK)):
            return True
        problems.append(f"{what} ({command}) is not there: {hint}")
        return False

    def built_here(command, name):
        problem, note = build_check(command, name)
        if problem and args.no_binary_check:
            problem, note = None, problem + " (--no-binary-check: run anyway)"
        if problem:
            problems.append(problem + ", or pass --no-binary-check")
        if note:
            notes.append(note)

    if executable(args.protal, "protal", "build it (docs/installation.md), or pass --protal"):
        built_here(args.protal, "protal")
    if executable(args.simulator, "the simulator", "it is built with protal (target simulate_metagenomes), or pass "
                                                   "--simulator"):
        built_here(args.simulator, "the simulator")
    # The simulator makes the long reads itself; a qshmm setup's reads follow pbsim3's model file (pbsim3 installs it).
    for kind, setup in (("pb", args.pb_setup), ("ont", args.ont_setup)):
        if kind in read_types and setup.startswith("qshmm:") and len(setup.split(":")) > 1 and \
                not pbsim_model_found(setup.split(":")[1], args.pbsim_models, args.pbsim):
            problems.append(f"{kind} reads follow pbsim3's model {setup.split(':')[1]}, which is in neither "
                            f"--pbsim-models nor pbsim3's data folder (next to {args.pbsim}): install pbsim3 (e.g. "
                            "micromamba install -c conda-forge -c bioconda pbsim3), pass --pbsim-models, or leave "
                            f"{kind} out of --read-types")
    _, imports_said = importing.communicate()
    if importing.returncode:
        problems.append(f"{sys.executable} cannot import what the trainer needs ({(imports_said.strip().splitlines() or ['?'])[-1]}): "
                        "install scikit-learn, joblib, numpy and pandas (envs/protal-db-build.yaml), or run this "
                        "script with a Python that has them")
    if problems:
        sys.exit("Cannot start:\n  " + "\n  ".join(problems))
    for note in notes:
        say("Note: " + note)
    if not shutil.which("zstd"):
        say("Note: no zstd command (envs/protal-db-build.yaml has it): the conversion writes full_reference.fna "
            "uncompressed, 86 GB at r226")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--inputs", help="folder of scripts/download_gtdb.py: the release, the genomes to simulate from "
                                    "and the species pool (instead of --gtdb, --extra-genomes, --simulate-species)")
    p.add_argument("--gtdb", help="extracted GTDB release directory")
    p.add_argument("--outdir", required=True,
                   help="output root: the database in OUTDIR/protal_db, the evaluation in OUTDIR/model_logs, each step's "
                        "log in OUTDIR/logs and console.log, what a rerun reuses in OUTDIR/work")
    p.add_argument("--release", help="GTDB release number; detected from taxonomy filenames by default")
    p.add_argument("--genome-table", help="optional simulator table: accession, taxonomy, whole genome FASTA path "
                                          "(and genome length; without, the run writes a copy with the lengths, "
                                          "OUTDIR/work/genomes.tsv)")
    p.add_argument("--genome-cache", default="auto", metavar="FILE",
                   help="each genome FASTA's length and contig names, by its path, size and modification time, kept "
                        "across builds, so that a rebuild reads none of the genomes for them again (at r226 one pass over "
                        "~50,000 genomes on a network file system took 8-22 min). auto (the default): "
                        f"INPUTS/{GENOME_CACHE} with --inputs (where the folder is writable), else OUTDIR/work/{GENOME_CACHE}; "
                        "none: read every genome")
    p.add_argument("--extra-genomes", action="append", default=[],
                   help="folder of more whole genomes of GTDB species, found by the accession in their file names "
                        "(e.g. the NCBI genomes of download_gtdb.py); repeatable")
    p.add_argument("--holdout", type=float, default=0.3,
                   help="fraction of the species (of those the genome table can simulate and no clade of "
                        "--holdout-clades took) left out of a separate training database (default 0.3; 0.2 until "
                        "2026-10-03: more species missing from the database give the model more of the false "
                        "positives it must learn to reject, docs/claude/2026-10-03-false-positive-fixes): their "
                        "reads land on relatives, as those of species GTDB lacks do in real samples, and they are "
                        "the false positives the model must learn to reject. 0 with --holdout-clades none trains "
                        "on the database itself")
    p.add_argument("--holdout-clades", default="phylum:2,class:4,order:6,family:8,genus:12",
                   help="whole clades left out of the training database too, RANK:COUNT for ranks phylum, class, "
                        "order, family, genus (default phylum:2,class:4,order:6,family:8,genus:12; none for "
                        "species only): reads of organisms whose genus, family, order, class or phylum the "
                        "database lacks land on ever more distant relatives. The report gives false positive and "
                        "false negative rates by rank")
    p.add_argument("--holdout-max-share", type=float, default=0.02,
                   help="a held-out clade has at most this share of the database's species (default 0.02)")
    p.add_argument("--novel-clades-per-sample", type=int, default=1,
                   help="species of held-out clades in each sample, per rank held out (default 1; "
                        "collect_training_data.py --novel_clades)")
    p.add_argument("--holdout-species", help="file of the species to leave out (optionally with their rank and "
                                             "clade, as heldout_species.txt), instead of choosing them")
    p.add_argument("--holdout-complex-distance", type=float, default=HOLDOUT_COMPLEX_DISTANCE,
                   help="species within this distance of a congener on their references' marker genes (species_clouds.tsv, "
                        "protal --write_species_neighbours on the converted release) form a complex that the training "
                        f"database leaves out or keeps whole, never torn (default {HOLDOUT_COMPLEX_DISTANCE}: a strain's "
                        "distance; 0: species drawn one by one, as before 2026-10-08). Torn, a held-out twin's reads teach "
                        "'absent' at the identity where a divergent strain teaches 'present' (docs/claude/2026-10-08-r226-v17)")
    p.add_argument("--species-clouds", help="the species' nearest congeners to steer the hold-out by (species_neighbours.tsv's "
                                            "format, from an earlier build's database or protal --write_species_neighbours), "
                                            "instead of comparing the converted release's references here")
    p.add_argument("--one-build-at-a-time", action="store_true",
                   help="build the finished database after the model is trained, not while the training data are "
                        "collected (the default, which needs the memory of two builds, or of one build and the "
                        "collection's protal runs, at once)")
    p.add_argument("--training-db-level", type=int, default=3,
                   help="zstd level of the training database (default 3)")
    p.add_argument("--strain-alleles", type=int, default=4,
                   help="up to this many strain alleles per species and gene in both databases (protal --build --strain_alleles: "
                        "the other genomes' copies in the full reference as edits of the representative's, for the alignment "
                        "scores and the 'alleles' features; default 4, 0: none)")
    p.add_argument("--allele-genome-share", type=float, default=0.5,
                   help="the share of each species' genomes that give strain alleles, by a hash of the accession "
                        "(protal --allele_genome_share); the genome table loses those but the representatives, so that no "
                        "simulated strain is its species' own allele, and the ancestry report takes its alleles from the "
                        "same genomes (default 0.5)")
    p.add_argument("--index-alleles", type=float, default=1, metavar="DIVERGENCE",
                   help="the k-mers of the strain alleles at least this far from the representative's copy (edits per base) "
                        "go into both databases' indexes under the species as non-unique entries, so that a read of a deep "
                        "strain finds its species (protal --build --index_alleles; default 1: none, off since the r226 v22 "
                        "build, whose gain was too small for the entries; 0.01 indexes the alleles 1%% or more from the representative)")
    p.add_argument("--foreign-rates", action="store_true",
                   help="scan each database's full reference for the gene copies other species' reads reach "
                        "(scripts/foreign_rates.py) and store foreign_rates.tsv in it, for the 'foreign' features. Off by "
                        "default since the r226 v18 build: the scan tiles every species' marker genes alike, so the table "
                        "tells the models nothing of the simulation's species (the scan of the genomes at hand did, r226 "
                        "v17: docs/claude/2026-10-07-congener-gaps), but it added 0.001 of AUC where strains and novel "
                        "congeners overlap, for two scans of ~21 min, and the features are in no default set "
                        "(docs/claude/2026-10-08-r226-v18)")
    p.add_argument("--no-foreign-rates", action="store_true",
                   help="no scan (the default): the 'foreign' features are unknown (-1); overrides --foreign-rates")
    p.add_argument("--foreign-stride", type=int, default=250,
                   help="the scan's reads: 150 bases every this many bases of a gene copy (default 250: four reads of a "
                        "marker gene of 1 kb)")
    p.add_argument("--foreign-per-header", type=int, default=10,
                   help="the scan tiles at most this many copies of each species' gene, the full reference's first (default "
                        "10; 0: every genome's), so that a species with thousands of genomes gives no more reads than one "
                        "with ten")
    p.add_argument("--final-db-level", type=int, default=9,
                   help="zstd level of the finished database (default 9: at GTDB r226, protal's default 19 made the index "
                        "2.7%% smaller than level 3 for 21 more minutes of a 1:20 build; docs/claude/2026-10-02-r226-build-"
                        "evaluation)")
    p.add_argument("--n-genes", type=int, metavar="N",
                   help="build the database from the N most distinctive of the release's marker genes (a reduced "
                        "database: less memory, fewer hits per genome): the genes ranked by prevalence x unique "
                        "k-mer share (scripts/rank_genes.py) from a full build of the training database first "
                        "(model_logs/gene_ranking.tsv; --gene-ranking skips that build), the N best in "
                        "model_logs/gene_subset.txt; the training database and the finished database hold those genes "
                        "only, with their neighbours counted over them, and the models are trained on them")
    p.add_argument("--genes", metavar="LIST",
                   help="build the database from these marker genes instead: GTDB marker ids (PF00380.20, "
                        "TIGR00001; PF00380 matches any version) or protal gene ids (gene2geneid.tsv), "
                        "comma-separated or one per line in a file (first column, # comments)")
    p.add_argument("--gene-ranking", metavar="FILE",
                   help="with --n-genes: a ranking of this release's genes by scripts/rank_genes.py (from an earlier "
                        "full build, e.g. a run with --rank-genes), instead of building one")
    p.add_argument("--genes-per-domain", type=int, metavar="M",
                   help="with --n-genes: how many of the N genes are reserved for the best genes of each domain, "
                        "bacteria and archaea (default a third of N, at least 1): a gene of the archaeal marker set "
                        "alone scores low over all species, and without this a small subset would hold no gene "
                        "archaea have")
    p.add_argument("--rank-genes", action="store_true",
                   help="in a run with every gene: rank the genes from the training database once it is built "
                        "(unpacked on the scratch disk, scripts/rank_genes.py) into model_logs/gene_ranking.tsv, for a "
                        "reduced database of the same release (--n-genes N --gene-ranking OUTDIR/model_logs/gene_ranking.tsv, "
                        "which then builds no ranking database)")
    p.add_argument("--simulate-species",
                   help="file of the species to simulate from (e.g. simulation_species.txt of download_gtdb.py: "
                        "species with other strains, and some without): the simulator draws species uniformly, so "
                        "among all ~130,000 GTDB species the few with downloaded strains would hardly be drawn")
    p.add_argument("--protal", default="protal", help="protal executable")
    p.add_argument("--simulator", default="simulate_metagenomes", help="simulate_metagenomes executable")
    p.add_argument("--no-binary-check", action="store_true",
                   help="run protal and the simulator even when they were not built from the source these scripts are "
                        "at (another version, or a commit whose source has changed since): the check at the start "
                        "then only notes it")
    p.add_argument("-t", "--threads", type=int, default=8)
    p.add_argument("--samples", type=int, default=12,
                   help="samples per design point (default 12; on a GTDB-like world the model still improved "
                        "from 60 to 120 samples)")
    p.add_argument("--read-pairs", default="1000,2000,5000,20000,50000,100000,200000,500000,2000000:4,10000000:2,30000000:1",
                   help="read pairs per sample, one design point each, DEPTH:SAMPLES for other samples than --samples "
                        "(default 1000,2000,5000,20000,50000,100000,200000,500000,2000000:4,10000000:2,30000000:1: the "
                        "deepest point because a model with the sample's depth as a feature cannot extrapolate past the deepest "
                        "sample it saw, and real metagenomes reach it; without the shallowest, a model "
                        "missed 8%% of the present taxa of samples of 1000 read pairs; at GTDB r226 the absent taxa per "
                        "sample grew as depth^0.84 to 500,000, and real samples are often 5-50M; 2000, 50000 and "
                        "200000 since r226 v10, whose test set had 42%% of its paired-end false positives at 50,000 "
                        "read pairs and called 15%% of the absent taxa at 2000, depths between the points the model "
                        "was trained at, docs/claude/2026-10-04-r226-v10-evaluation)")
    p.add_argument("--read-setups", default="100:HS20:300:40,150:HSXt:350:50,250:MSv3:550:50",
                   help="LENGTH:INSTRUMENT:FRAGMENT_MEAN:FRAGMENT_SD, one design point each; INSTRUMENT: HS20, "
                        "HS25, HSXt (HiSeq X Ten), NovaSeq or MSv3 (MiSeq v3), simulate_metagenomes's models of them "
                        "(docs/databases.md#illumina-reads)")
    p.add_argument("--archaea", type=int, default=2)
    p.add_argument("--species-per-sample", default="20-200",
                   help="species per sample, drawn per sample (default 20-200: real gut samples hold 100-300 GTDB "
                        "species with a long tail of rare ones)")
    p.add_argument("--strains-per-species", default="0.5,0.2",
                   help="probabilities of a second, third, ... strain of a species in a sample (default 0.5,0.2, the "
                        "test set's; 0.3,0.1 before 2026-10-07): real samples often mix strains, which changes the "
                        "allele-frequency features. With fewer strains than the test set, the global knob was chosen on "
                        "samples that lost fewer strains to it: at r226 v15 pe's 0.80 lost the test design's samples "
                        "twice the present taxa per sample it lost the training design's, 90%% of them strains "
                        "(docs/claude/2026-10-07-r226-v15)")
    p.add_argument("--abundance", default="lognormal:1.3,2.0",
                   help="abundance model of the training samples: lognormal:SIGMA (continuous, no species below 1/1000 "
                        "of the median; until 2026-10-07 Poisson counts + 1, now poisson_lognormal:SIGMA), "
                        "powerlaw:ALPHA or negbin:R:P; lognormal:S1,S2,... gives a design point's samples the sigmas in "
                        "turn (default "
                        "lognormal:1.3,2.0: half the samples with the former sigma 1.3, half with the test set's 2.0, "
                        "so that the model's depth prior does not rest on one abundance distribution, "
                        "docs/claude/2026-10-03-false-positive-fixes)")
    p.add_argument("--read-types", default="pe,se,pb,ont",
                   help="the read types to train a model for (default pe,se,pb,ont): se from the paired-end "
                        "samples' first reads, pb and ont from long reads of the same communities (made by "
                        "simulate_metagenomes: HiFi reads by hifi_reads.py's model, Nanopore reads by pbsim3's qshmm "
                        "model). The models are trained in parallel; read types left out keep placeholders")
    p.add_argument("--long-read-bases", default="300000,1500000,6000000,30000000,150000000,1500000000:4,6000000000:2",
                   help="bases per long-read sample, one design point each, DEPTH:SAMPLES for other samples than "
                        "--long-read-samples (collect_training_data.py; real HiFi and Nanopore metagenomes are 5-30 Gb)")
    p.add_argument("--long-read-samples", type=int, default=36,
                   help="samples per long-read design point (default 36, at most the communities of the paired-end "
                        "points of its depth, 12 samples of 3 setups by default: at r226 the Nanopore model's learning "
                        "curve still fell at 126 samples)")
    p.add_argument("--pb-setup", default="hifi:15000:3000:3",
                   help="PacBio reads: hifi:LENGTH_MEAN:LENGTH_SD:Q_SD, HiFi reads by hifi_reads.py's model (made by "
                        "simulate_metagenomes), their quality by their length (Q50 at 5 kb to Q30 at 25 kb, Q20 at 50 kb) "
                        "and Q_SD around it (default), or a qshmm setup as --ont-setup's")
    p.add_argument("--ont-setup", default="qshmm:QSHMM-ONT-HQ:8000:6000:0.97:39/24/36",
                   help="Nanopore reads: qshmm:MODEL:LENGTH_MEAN:LENGTH_SD:ACCURACY_MEAN[:SUB/INS/DEL], pbsim3's "
                        "quality-score model (its MODEL.model file), as pbsim3 --strategy templ makes them, made by "
                        "simulate_metagenomes")
    p.add_argument("--pbsim", default="pbsim", help="pbsim3's executable, only to find its data folder of .model files "
                                                    "(pbsim3 itself is not run)")
    p.add_argument("--pbsim-models", help="folder of pbsim3's .model files (default: pbsim3's data folder)")
    p.add_argument("--insilico-strains", type=float, default=1.0, metavar="SHARE",
                   help="give this share of the species with one genome in the genome table an in-silico strain to be "
                        "simulated from, too (scripts/insilico_strains.py: a copy of the representative with "
                        "codon-aware substitutions, as far from it as the table's real strains are from theirs, each "
                        "marker gene by its conservation factor; default 1, 0: none). Without them one-genome species "
                        "were always simulated from the database's own reference, and a model given the GTDB cluster "
                        "sizes (+priors) learned to reject divergent reads on them (docs/claude/2026-10-04-r226-v10-"
                        "evaluation). They go to the samples' disk (SCRATCH/insilico_strains, OUTDIR/work/insilico_strains "
                        "without --scratch), the table simulated from to OUTDIR/work/genomes_simulated.tsv")
    p.add_argument("--insilico-ani", metavar="MIN-MAX",
                   help="draw the in-silico strains' genome ANI uniformly from MIN-MAX (e.g. 95-99) instead of the "
                        "real strains' marker divergence")
    p.add_argument("--insilico-multi", type=float, default=0.5, metavar="SHARE",
                   help="give this share of the species with real strains in the genome table an in-silico strain as "
                        "well (insilico_strains.py --multi-share; default 0.5, 0: none): a strain its species' alleles "
                        "do not reach, as a strain of a lineage GTDB has not sampled is. Trained on real strains alone, "
                        "half of whose relatives are the alleles, the models learn that a strain the alleles explain is "
                        "present and one they do not is not, and call the strains of lineages without an allele genome "
                        "at 0.3-0.6 against 0.98 (docs/claude/2026-10-09-ancestry-true-positive-test)")
    p.add_argument("--test-samples", type=int, default=4,
                   help="samples per design point of the independent test set (default 4; 0: none). The test set "
                        "has another design than the training data (--test-*), so that the report shows what "
                        "cross-validation on the training data cannot")
    p.add_argument("--test-read-pairs", default="500,2000,10000,50000,200000,1000000,5000000:2")
    p.add_argument("--test-species-per-sample", default="10-300")
    p.add_argument("--test-abundance", default="lognormal:2.0", help="(default lognormal:2.0: more uneven than training)")
    p.add_argument("--test-strains-per-species", default="0.5,0.2")
    p.add_argument("--test-long-read-bases", default="150000,1000000,5000000,25000000,250000000,3000000000:2")
    p.add_argument("--test-long-read-samples", type=int, default=8,
                   help="samples per long-read design point of the test set (default 8, at most the communities of "
                        "its paired-end points: with 4, the 22 samples of a long-read type could not tell its knob "
                        "curve from one knob)")
    p.add_argument("--scenarios", default=None,
                   help="scenarios to train on and score the models on, besides the design (default: all five; none "
                        "for none; scenarios.py): gut (150-1,000 "
                        "species, 5%% lacking from the database, Illumina PE 150 Q35 at 20M pairs, PacBio and Nanopore "
                        "at the same bases), moderate (1,000-5,000 species, 30%% lacking, Illumina at 10M pairs, Ultima "
                        "at 10M reads, PacBio and Nanopore at 3 Gb), soil (3,000-11,000 species, 60%% lacking, Ultima SE "
                        "300 Q25 at 20M reads, "
                        "Illumina, PacBio and Nanopore at the same bases), soil_shallow (soil communities at 5M "
                        "Illumina pairs, PacBio and Nanopore at 1.5 Gb), host (90%% human reads, 2-50 species of "
                        "power-law abundances, all four read types at 10M pairs or reads or 3 Gb), or those of "
                        "--scenario-file; comma-separated, NAME:N for N hold-in samples of one, all for every one. Each "
                        "scenario's hold-in samples (--scenario-samples) join the training data, as the design's do, "
                        "its hold-out samples (--scenario-test-samples, another seed) the test set, and every model's "
                        "report gives F1 and FP and FN rates on both (section Scenarios; model_logs/summary.txt). A "
                        "scenario larger than the genome table holds at its share of species the database lacks is "
                        "scaled down, and the run says so: soil's full size needs about 25,000 species to simulate from, "
                        "the download's default since 2026-10-05 (8,000 before). host needs the host genome (--host-genome, "
                        "or the download's); by default, without one, it is left out with a warning. Each sample's "
                        "depth is its scenario's times a factor from 1/8 to 2 (depth_range; 1/2 to 2 before "
                        "2026-10-07), its species count from the scenario's range, both log-uniform and stratified, "
                        "and its evenness one of four lognormal sigmas")
    p.add_argument("--scenario-file", help="JSON of scenarios by name, which add to or change the presets (scenarios.py)")
    p.add_argument("--scenario-samples", type=int, default=10,
                   help="hold-in samples per scenario, in the training data (default 10; 6 before 2026-10-07, 3 before "
                        "2026-10-06; 0: the scenarios are scored, not trained on). Sample-level features (the depth, "
                        "the sample's complexity) are learned from these few samples per scenario: with 6, the one at "
                        "the edge of its scenario's was scored like the design's samples (r226 v15)")
    p.add_argument("--scenario-test-samples", type=int, default=4,
                   help="hold-out samples per scenario, in the test set (default 4; 3 before 2026-10-07, 2 before "
                        "2026-10-06)")
    p.add_argument("--scenario-weight", type=float, default=0.25,
                   help="the sample weight of the scenarios' hold-in rows in the models, the design's 1 "
                        "(machine_learning_cmdline.py --scenario-weight; default 0.25: at r226 v12 their rows were 52-80%% "
                        "of the training rows, and at full weight they cost the design's test set 0.002-0.004 of F1, at "
                        "0.25 half to two thirds less with the soil scenarios' gain kept; "
                        "docs/claude/2026-10-05-r226-v12-scenarios)")
    p.add_argument("--host-genome",
                   help="FASTA (gzipped or not) of the host genome for scenarios with host reads (default: the human "
                        "genome download_gtdb.py fetched into --inputs)")
    p.add_argument("--error-reads", default="all",
                   help="the samples whose reads behind the models' errors are kept: protal writes them with an unmapped "
                        "record for every read that seeded on taxa but aligned nowhere (the non-hits; the profiles are "
                        "the same), and once the models are trained model_logs/error_reads/<read type>/ tells, in one "
                        "table of every sample's error taxa (taxa.tsv.gz), where the reads behind the model's false "
                        "positives and false negatives went, with their source genomes (error_reads.py). all (default: "
                        "every sample of the training data and the test set, every read type), none, or READ_TYPE, "
                        "READ_TYPE:design (the design's samples) or READ_TYPE:SCENARIO, comma-separated; read types and "
                        "scenarios not collected are left out. With --share-logs the reads' SAM records are kept too")
    p.add_argument("--share-logs", action="store_true",
                   help="keep the SAM records of the reads behind the models' errors (model_logs/error_reads/<read "
                        "type>/<set>/<design point>/: per sample <sample>.FP.sam.zst and <sample>.FN.sam.zst, at most 20 "
                        "fragments per taxon and reason, no qualities), measure which side those reads take where the "
                        "species differs from its congeners (model_logs/ancestry_sites/, ancestry_sites.py; the training "
                        "database's full reference is kept until then), and at the end pack OUTDIR/<OUTDIR's "
                        "name>_share.tar.gz: console.log, logs/, model_logs/ (without the models), the build's metadata, "
                        "the taxonomy and the training and test tables (their numbers to 9 significant digits), and from the "
                        "samples' disk the design points' simulator logs and parameters (logs/simulations_<collection>.log) "
                        "and the in-silico strains' table, to be copied off the cluster and read elsewhere")
    p.add_argument("--congeners", default="0.25:2-5", type=congener_spec,
                   help="relatives that share a sample, in the training data and the test set (collect_training_data.py "
                        "--congeners): SHARE:MIN-MAX, about SHARE of each sample's species in groups of MIN to MAX "
                        "species of one genus (default 0.25:2-5); N for N species of one genus per design point; 0 for "
                        "none. Real samples often hold congeners at very different abundances, uniform draws hardly "
                        "ever; and a model with the relatives features (--features normalized+adjacency+relatives) "
                        "trained on uniform draws learns that an abundant congener means a taxon is absent: on the "
                        "r226 tables such a model missed 69-73%% of the present species beside a congener 10 times as "
                        "abundant (docs/claude/2026-10-02-amplicon-denoising). On the benchmark world the groups left "
                        "the other models' test F1 within noise (docs/claude/2026-10-03-denoising-implementation)")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--model", choices=["gbm", "forest"], default="gbm",
                   help="the presence models (machine_learning_cmdline.py --model): gbm (default since 2026-10-06), "
                        "gradient-boosted trees, which beat the forest on the r226 v13 tables "
                        "(docs/claude/2026-10-06-r226-v13-soil); forest, a random forest, as before")
    p.add_argument("--ntree", type=int, default=64, help="a forest's trees (default 64)")
    p.add_argument("--rounds", type=int,
                   help=f"boosting's rounds (machine_learning_cmdline.py --rounds; default {TRAINER_ROUNDS})")
    p.add_argument("--maxnodes", default=None,
                   help=f"leaves per tree at most (machine_learning_cmdline.py --maxnodes), N or TYPE:N items, a bare N for "
                        f"the read types not named (default {BOOSTED_LEAVES} for boosting, {FOREST_LEAVES} for a forest: "
                        "at r226 512 leaves gave the short-read forests a lower log loss and fewer false positives than "
                        "256, 128 the long-read forests a lower log loss at the same F1)")
    p.add_argument("--no-gene-neighbours", action="store_true",
                   help="do not record which marker genes lie next to which in the genomes to simulate from "
                        "(gene_neighbours.py; protal then pairs no mates across neighbouring genes)")
    p.add_argument("--no-placeholder-models", action="store_true",
                   help="leave the se, pb and ont models out of the database (a run with such reads then stops "
                        "with an error) instead of placeholders that report no species until trained ones replace "
                        "them (placeholder_models.py)")
    p.add_argument("--evaluation", choices=["full", "basic", "none"], default="basic",
                   help="how much the trainer evaluates (machine_learning_cmdline.py --evaluation): basic (default since "
                        "2026-10-06), the models with rows, samples and species held out, which the summary and the "
                        "knob need; full also with clades held out and the studies (other feature sets, model size, "
                        "fewer samples), which multiply a boosted model's training several times")
    p.add_argument("--features", type=feature_set_name, default=DEFAULT_FEATURE_SET,
                   metavar="auto|" + "|".join(FEATURE_SETS[:2] + ("...",)),
                   help=f"the models' features (machine_learning_cmdline.py --features): {DEFAULT_FEATURE_SET} "
                        "(default since 2026-10-06: the set every model of the r226 v12 and v13 builds chose with auto, "
                        "and the reference's k-mer uniqueness); auto, each trainer chooses its set (below), which "
                        "doubles a boosted model's training with --evaluation basic; or feature groups joined by '+': "
                        "normalized, the normalised features; adjacency, the gene "
                        "neighbours'; distance, the four relative_* features that compare a taxon with its sample's "
                        "relatives by the distance of their references (they need --congeners groups; at r226 +0.004 "
                        "paired-end F1 at the knob curve, the other read types within noise, "
                        "docs/claude/2026-10-03-r226-v5-v6-training), or relatives, all the relatives features (better "
                        "ranking, no better at the knob curve there); depth, the sample's depth (with it the trainer "
                        "fits no knob curve and protal calls at --knob; on the r226 v5 tables false positives halved "
                        "at knob 0.5, docs/claude/2026-10-03-false-positive-anatomy); divergence, divergence by gene "
                        "conservation and codon position and lost mates; unfiltered, the reads before the MAPQ and "
                        "length filters and the failed candidates; ref, the reference's k-mer uniqueness in the "
                        "database; priors, what GTDB knows of the species "
                        "(docs/claude/2026-10-03-false-positive-fixes); all: every feature of the training dumps. "
                        "normalized+adjacency+distance is the set of protal 0.7.3's dumps, normalized+adjacency that "
                        "of older ones. auto: each model's set chosen by its trainer, the one of highest F1 with "
                        "species held out among those without the priors (auto+priors: with them), the default set "
                        "unless another is 0.002 better (machine_learning_cmdline.py choose_feature_set); the run says "
                        "which won and why, and each set's F1 on the test sets")
    p.add_argument("--call-mode", choices=["curve", "fdr"], default="curve",
                   help="how protal calls with the models by default: curve (default), the knob curve over the "
                        "sample's depth (--depth-knob-read-types); fdr, the highest-scoring taxa of each sample while "
                        "their expected share of false calls stays at the target the trainer chose "
                        "(machine_learning_cmdline.py --fdr-calls; protal --fdr), which follows the sample's number of "
                        "candidates at any depth (on the benchmark world 0.0005-0.004 F1 below the curve, "
                        "docs/claude/2026-10-03-denoising-implementation; at r226 0.001-0.007 below it, missing most "
                        "species of the shallowest samples, docs/claude/2026-10-03-r226-v5-v6-training). With fdr the "
                        "trainer reports both on the "
                        "test set")
    p.add_argument("--previous-procedure", action=argparse.BooleanOptionalAction, default=False,
                   help="the trainer also compares with its previous procedure (machine_learning_cmdline.py "
                        "--previous-procedure): for the first builds of a release; off by default, as it takes "
                        "most of the training time")
    p.add_argument("--depth-knob-read-types", default="pe,se,pb,ont",
                   help="read types (comma-separated) whose models also get a knob curve over the sample's depth "
                        "(machine_learning_cmdline.py --depth-knobs; default all four, '' for none): at GTDB r226 the best "
                        "threshold went from ~0.1 at 1,000 read pairs to ~0.9 at 500,000, and thresholds by depth raised "
                        "the test sets' F1 by 0.006-0.033 (docs/claude/2026-10-02-r226-build-evaluation). A model with "
                        "the sample's depth among its features (--features ... depth, the default) gets no curve: the "
                        "forest learns the depth itself and protal calls at --knob")
    p.add_argument("--progress-every", type=float, default=0,
                   help="seconds between status lines of the stages running, besides the lines of each step's start "
                        "and end: how long each has run, the memory it takes and the last line of its log (default "
                        "0: none)")
    p.add_argument("--scratch",
                   help="a fast local disk (a compute node's own) for the simulated samples: their reads, alignments, "
                        "profiles and the simulators' temporary files go to SCRATCH/training and SCRATCH/test, and "
                        "only the tables to OUTDIR/work; the training database is built there too (SCRATCH/training_db, "
                        "~22 GB at r226). A network file system (OUTDIR's, often) is slow at the many files the "
                        "simulators write and delete, and at writing a database; the converter spools the release's "
                        "marker genes there too. The run estimates before its builds what it needs there (the genome store "
                        "--genome-store, ~50 GB at r226; the training database, ~22 GB; the samples' SAMs and profiles; "
                        "their reads) and streams into protal what would not fit (--stream-above auto); its end says the "
                        "most it took against that estimate. A rerun reuses the samples and the store in the same SCRATCH")
    p.add_argument("--profile-blocks", type=float, default=20.0,
                   help="once the training database is built, profile the simulated samples as their simulations go "
                        "on, in protal runs of at least this many GB of reads (collect_training_data.py --follow; the "
                        "two collections' runs take turns), and remove each design point's reads once all its read types "
                        "are profiled (the SAMs, profiles and dumps are kept): the samples need not all be on the disk "
                        "at once. 0: profile both collections in one protal run once every sample is simulated, and "
                        "keep the reads (default 20)")
    p.add_argument("--keep-free", type=float, default=30.0,
                   help="with --profile-blocks: GB a simulation leaves free on the disk of the samples (--scratch or "
                        "OUTDIR/work), or it waits until profiled reads are removed (default 30)")
    p.add_argument("--stream-above", type=stream_spec, default="auto",
                   help="with --profile-blocks: a design point whose largest sample's reads would take more than this "
                        "many GB (compressed, estimated) is not written to the disk: protal reads its samples from named "
                        "pipes as simulate_metagenomes makes them, a sample at a time, in protal runs with the other "
                        "points streamed (collect_training_data.py --stream_above); smaller ones are written and profiled "
                        "in blocks. auto (the default since 2026-10-07; 2 before): chosen once the training database is "
                        "built, from the free space on the samples' disk less what the run needs there besides the reads "
                        "(the genome store, a database still to be built there, the samples' SAMs and profiles, "
                        "--keep-free): nothing streamed if every simulation's reads fit at once, else the largest "
                        "simulations, until the others' fit; the run stops before its builds if even streaming all "
                        "would not fit. A number: that many GB (0: none)")
    p.add_argument("--profile-block-max", type=float, default=200.0,
                   help="with --profile-blocks: the most GB of written reads one protal run takes (the files' size), at "
                        "least one design point's, so that a failed run costs at most that much profiling (default 200; "
                        "0: no limit). The streamed simulations whose communities are there share one run whatever their "
                        "size (until 2026-10-09 up to this cap): their reads are on no disk")
    p.add_argument("--profile-ahead", action=argparse.BooleanOptionalAction, default=True,
                   help="protal --profile_ahead in the collections' protal runs (the default; --no-profile-ahead: not): "
                        "each sample profiled while the next is aligned, on a quarter of the threads, so that a run's "
                        "profiling stage after its last alignment shrinks to what the worker had not started (at r226 v19 "
                        "the six runs' profiling stages took 15.5 min after their alignment; the profiles are the same)")
    p.add_argument("--genome-store", default="auto",
                   help="simulate_metagenomes's genome store for both collections (collect_training_data.py "
                        "--genome_store): each genome simulated is read and parsed from its FASTA once, written there "
                        "at 2 bits a base, and memory-mapped by every later sample, long-read round and simulation "
                        "instead of being inflated and parsed again (with the defaults at r226 ~1.5M genome reads of "
                        "~54k genomes, ~28 each); the same reads. auto (the default): SCRATCH/genome_store with "
                        "--scratch, OUTDIR/work/genome_store otherwise; or a folder; none: no store. It takes ~0.25 bytes a "
                        "base of the genomes simulated, ~50 GB at r226 besides the samples' space, and is kept for the "
                        "next build (a FASTA that changed is read again)")
    p.add_argument("--compressed-pipes", action="store_true",
                   help="with --stream-above: the streamed samples go through their named pipes compressed, as their "
                        "names say (collect_training_data.py --compressed_pipes); by default as plain FASTQ, which saves "
                        "compressing them only for protal to inflate them again")
    p.add_argument("--read-compression", choices=["zstd", "gzip"], default="zstd",
                   help="how the simulated samples' reads are written (collect_training_data.py --read_compression): "
                        "zstd (.fq.zst, the default: as small as BGZF or smaller, several times faster to write) or "
                        "gzip (.fq.gz); protal reads both")
    args = p.parse_args()
    args.no_foreign_rates = args.no_foreign_rates or not args.foreign_rates  # the scan only with --foreign-rates
    if args.strain_alleles < 0 or not 0 <= args.allele_genome_share <= 1 or args.index_alleles < 0:
        p.error("--strain-alleles is a count (0: none), --allele-genome-share a share from 0 to 1, --index-alleles a "
                "divergence of 0 or more (1 or more: none)")
    ALLELE_ARGS[:] = ["--strain_alleles", str(args.strain_alleles), "--allele_genome_share", f"{args.allele_genome_share:g}",
                      "--index_alleles", f"{args.index_alleles:g}"]
    global ALLELE_SHARE
    ALLELE_SHARE = args.allele_genome_share if args.strain_alleles > 0 else 1.0
    Job.progress_every = args.progress_every
    read_types = [t.strip() for t in args.read_types.split(",") if t.strip()]
    if not read_types or any(t not in TABLES for t in read_types):
        p.error(f"--read-types: a comma-separated list of {', '.join(TABLES)}, got {args.read_types!r}")
    if args.maxnodes is None:
        args.maxnodes = BOOSTED_LEAVES if args.model == "gbm" else FOREST_LEAVES
    try:
        max_leaves(args.maxnodes, "pe")
    except ValueError as e:
        p.error(str(e))
    if args.n_genes is not None and args.genes:
        p.error("give --n-genes or --genes, not both")
    if args.n_genes is not None and args.n_genes < 1:
        p.error("--n-genes: at least 1")
    if args.gene_ranking and args.n_genes is None:
        p.error("--gene-ranking goes with --n-genes")
    if args.gene_ranking and not os.path.isfile(args.gene_ranking):
        p.error(f"--gene-ranking: {args.gene_ranking} is not a file")
    if args.genes and not os.path.isfile(args.genes) and "," not in args.genes and not re.fullmatch(r"[\w.]+", args.genes):
        p.error(f"--genes: {args.genes} is neither a file nor a list of marker or gene ids")
    if args.genes_per_domain is not None and (args.n_genes is None or args.genes_per_domain < 0 or
                                              args.genes_per_domain * 2 > args.n_genes):
        p.error("--genes-per-domain goes with --n-genes, and twice it is at most N")
    if args.rank_genes and (args.n_genes is not None or args.genes):
        p.error("--rank-genes is for a run with every gene; a run with --n-genes ranks (or takes) the genes itself")
    if args.rank_genes and args.holdout <= 0 and args.holdout_clades.strip().lower() == "none" and not args.holdout_species:
        p.error("--rank-genes ranks from the training database: hold species out")
    # The scenarios (scenarios.py): their names, and what they need beyond the design, checked before anything runs.
    if args.scenario_samples < 0 or args.scenario_test_samples < 0 or args.scenario_weight < 0:
        p.error("--scenario-samples, --scenario-test-samples and --scenario-weight cannot be negative")
    scenarios_given = args.scenarios is not None  # by default all five, the host scenario only with a host genome
    if not scenarios_given:
        args.scenarios = ",".join(scenarios.PRESETS)
    try:
        scenario_defs = scenarios.definitions(args.scenario_file)
        # A scenario's hold-in samples (NAME:N in --scenarios, else --scenario-samples) and hold-out samples.
        chosen_scenarios = scenarios.selection(args.scenarios, args.scenario_samples, scenario_defs, keep_zero=True)
    except scenarios.ScenarioError as e:
        p.error(str(e))
    hold_in = {name: n for name, n in chosen_scenarios if n > 0}
    hold_out = {name: args.scenario_test_samples for name, _ in chosen_scenarios} if args.scenario_test_samples else {}
    selected = [name for name, _ in chosen_scenarios if name in hold_in or name in hold_out]
    args.scenario_samples_of = {name: (hold_in.get(name, 0), hold_out.get(name, 0)) for name in selected}  # provenance
    if selected and not any(r["type"] in read_types for name in selected for r in scenario_defs[name]["reads"]):
        p.error(f"--scenarios {args.scenarios}: none of their read types is among --read-types {args.read_types}")
    hosted = [name for name in selected if scenario_defs[name]["host_share"] > 0]
    with concurrent.futures.ThreadPoolExecutor(1) as beside:  # the versions read while the tools are checked
        versions = beside.submit(tool_versions, args)
        check_tools(args, read_types)
        args.versions_at_start = versions.result()  # build_metadata.tsv: what the run starts with (versions_at_end)
    state = {}
    if args.inputs:
        if args.gtdb:
            p.error("give --inputs or --gtdb, not both")
        state_file = os.path.join(args.inputs, "download.json")
        if not os.path.isfile(state_file):
            p.error(f"{args.inputs} has no download.json: is it a folder of download_gtdb.py?")
        with open(state_file) as fh:
            state = json.load(fh)
        args.gtdb = os.path.join(args.inputs, "release")
        args.release = args.release or state["release"]["number"]
        # The genomes and species to simulate shape the genome table this script makes, not a given one.
        if os.path.isdir(os.path.join(args.inputs, "genomes")) and not args.genome_table:
            args.extra_genomes.append(os.path.join(args.inputs, "genomes"))
        if not args.simulate_species and not args.genome_table and \
                os.path.isfile(os.path.join(args.inputs, "simulation_species.txt")):
            args.simulate_species = os.path.join(args.inputs, "simulation_species.txt")
    elif not args.gtdb:
        p.error("give --inputs (a folder of download_gtdb.py) or --gtdb")
    scenario_notes = []  # what the run changed of the scenarios asked for, for its console and build_metadata.tsv
    if hosted:
        if not args.host_genome and state.get("host", {}).get("path"):
            args.host_genome = os.path.join(args.inputs, state["host"]["path"])
        if args.host_genome and os.path.isfile(args.host_genome):
            args.host_genome = os.path.abspath(args.host_genome)
        elif scenarios_given:
            p.error(f"scenario {', '.join(hosted)}: its host reads need the host genome: --host-genome FASTA, or an "
                    "--inputs folder whose download fetched it (download_gtdb.py, --host_genome)")
        else:  # the default scenarios: the others run, and the run says why this one does not
            scenario_notes.append(f"{', '.join(hosted)} left out: no host genome (--host-genome, or rerun download_gtdb.py "
                                  "on its --inputs folder, which fetches only what is missing)")
            hold_in = {n: k for n, k in hold_in.items() if n not in hosted}
            hold_out = {n: k for n, k in hold_out.items() if n not in hosted}
            selected = [n for n in selected if n not in hosted]
            args.scenario_samples_of = {n: v for n, v in args.scenario_samples_of.items() if n not in hosted}
            args.host_genome = None
    try:
        args.error_units = [(kind, scope) for kind, scope in error_read_units(args.error_reads, scenario_defs)
                            if kind in read_types and (scope in (None, "design") or (
                                scope in selected and any(r["type"] == kind for r in scenario_defs[scope]["reads"])))]
    except ValueError as e:
        p.error(f"--error-reads: {e}")
    os.makedirs(args.outdir, exist_ok=True)
    open_console_log(os.path.join(args.outdir, "console.log"))
    # OUTDIR's folders (the module's docstring): the database, the evaluation (and what the run chose), every step's
    # log, and what a rerun reuses.
    db, logs, steps, work = (os.path.join(args.outdir, f) for f in (DATABASE, REPORTS, LOGS, WORK))
    for folder in (db, logs, steps, work):
        os.makedirs(folder, exist_ok=True)
    # The run's CPU use, in logs only (not on the console): each command's CPU time and the cores it kept busy
    # (cpu_jobs.tsv, Job.cpu_row), and how many cores the run kept busy over time (cpu_timeline.tsv, CpuTimeline); the
    # collections add their protal runs' (protal_cpu_<collection>.tsv). A rerun starts them anew.
    for name in ("cpu_jobs.tsv", "cpu_timeline.tsv"):
        with contextlib.suppress(OSError):
            os.remove(os.path.join(steps, name))
    Job.cpu_log = os.path.join(steps, "cpu_jobs.tsv")
    timeline = CpuTimeline(os.path.join(steps, "cpu_timeline.tsv"))
    timeline.start()

    def step_log(name):
        return os.path.join(steps, name)
    if args.release:
        release = str(args.release).removeprefix("r")
    else:
        found = set()
        for name in os.listdir(args.gtdb):
            m = re.match(r"(?:bac120|ar53)_taxonomy_r(\d+)", name)
            if m:
                found.add(m.group(1))
        if len(found) != 1:
            sys.exit(f"Cannot detect one GTDB release in {args.gtdb}; specify --release")
        release = next(iter(found))
    samples_root = work  # where the collections simulate and profile their samples
    if args.scratch:
        samples_root = os.path.abspath(args.scratch)
        os.makedirs(samples_root, exist_ok=True)
        Job.scratch = Scratch(samples_root)
    say(f"A protal database of GTDB r{release} in {args.outdir}, with {args.threads} threads; each step logs to a "
        f"file in {LOGS}/" + (f"; the simulated samples go to {samples_root} "
                            f"({gigabytes(shutil.disk_usage(samples_root).free)} free)" if args.scratch else "") +
        (f"; the stages running are reported every {clock(Job.progress_every)}" if Job.progress_every > 0 else ""))
    # The steps: genome table, release, marker genes (with --n-genes or --genes), databases, training data, test
    # set (if any), models, parity, packing.
    subset = args.n_genes is not None or bool(args.genes)
    has_test = args.test_samples > 0 or bool(hold_out)  # a test collection: the design's test set, the scenarios' hold-out
    clades = parse_clades(args.holdout_clades)
    steer_holdout = (args.holdout > 0 or bool(clades)) and not args.holdout_species and args.holdout_complex_distance > 0
    # The species clouds steer the hold-out and give the in-silico strains their congener sites.
    make_clouds = steer_holdout or (args.insilico_strains > 0 and not args.no_gene_neighbours)
    Steps.total = 8 + has_test + subset + (args.insilico_strains > 0) + (not args.no_foreign_rates) + make_clouds
    if args.genes:
        # The list is checked against the release's marker files before anything is converted (marker ids by
        # name, gene ids by their range; the ids themselves come from gene2geneid.tsv once it is there).
        markers = sorted({m for m, _ in marker_files(args.gtdb, "genomic_files_reps", "reps", release)})
        if markers:
            read_gene_list(args.genes, {m: i + 1 for i, m in enumerate(markers)})
    # A rerun skips what an earlier run into OUTDIR completed with the same inputs (see Stages). The gene neighbours have
    # a stage of their own (they depend on the genome table, the conversion does not), so that the release can be
    # converted before the genome table is made.
    stages = Stages(os.path.join(work, "stages"))
    convert_key = {"converter": content_hash(CONVERTER), "release": release, "gtdb": release_identity(args.gtdb, release),
                   "placeholders": not args.no_placeholder_models}
    final_base = {"convert": convert_key, "protal": file_identity(args.protal), "level": args.final_db_level, "alleles": ALLELE_ARGS}
    # --build packs the taxonomy into database.protal; the collector and the trainer read it (domains,
    # representative genomes). The training database has the same taxonomy.
    taxonomy = os.path.join(work, "internal_taxonomy.dmp")
    gene_table = os.path.join(work, "gene2geneid.tsv")  # the markers' gene ids, for --genes and the ranking
    neighbours_key = None  # set once the genome table is there

    def start_convert(into):
        """The converter into `into`, in the background (its Job)."""
        # With --scratch, the converter spools the marker genes on the node's disk (~5 GB compressed at r226).
        return Job([sys.executable, CONVERTER, "--gtdb", args.gtdb, "--outdir", into, "--release", release, "-t",
                    str(args.threads)] + (["--tmp", samples_root] if args.scratch else []),
                   step_log("convert.log"), label="converting the release")

    def find_neighbours(into):
        """gene_neighbours.py on the converted files in `into`, in the background (its Job). Which marker genes lie next to
        which in every genome to simulate from (their genes placed by their sequence or k-mer trace): gene_neighbours.tsv,
        the frequencies per clade, and gene_positions.tsv, which --build both pack; the training database's copy derives
        the frequencies anew from the positions, without the species it leaves out."""
        return Job([sys.executable, GENE_NEIGHBOURS, "--db", into, "--genome_table", genome_table, "-t", str(args.threads)],
                   step_log("gene_neighbours.log"), label="finding the genes' neighbours")

    def neighbours_done(job, into):
        """What the gene neighbours' job said; their stage marked for `into`."""
        stages.mark("gene_neighbours", {**neighbours_key, "into": os.path.basename(into)})
        return (f"the genes placed{genes_placed(step_log('gene_neighbours.log'))} and their neighbours counted in "
                f"{job.took()} (gene_neighbours.log)")

    def finish_convert(job, into):
        """The converted files of the release in `into` (the converter's job, started by start_convert), with the
        models of the read types but pe as placeholders (protal warns when it loads one), packed by --build like
        model_pe.xml, and the gene neighbours: gene_neighbours.py starts once the converter has written reference.fna and
        reference.map (REFERENCE_WRITTEN in convert.log), while it writes the full reference, which the neighbours do not
        read (1-1.5 min at r226). Says how long it took."""
        neighbours = None
        if not args.no_gene_neighbours:
            log = step_log("convert.log")
            while job.poll() is None and not log_says(log, REFERENCE_WRITTEN):
                check_jobs()
                time.sleep(0.5)
            if job.poll() in (None, 0):
                neighbours = find_neighbours(into)
        job.finish()
        shutil.copyfile(os.path.join(into, "internal_taxonomy.dmp"), taxonomy)
        shutil.copyfile(os.path.join(into, "gene2geneid.tsv"), gene_table)
        took = f"converted in {job.took()}"
        if neighbours:
            took += "; " + neighbours_done(neighbours.finish(), into)
        if not args.no_placeholder_models:
            for read_type in ("se", "pb", "ont"):
                write_placeholder(os.path.join(into, MODEL_FILES[read_type]), read_type)
        return took

    def convert(into):
        return finish_convert(start_convert(into), into)

    converted = None  # the folder with the converted files, once there
    # With a gene subset the whole release is converted here, and the database folders (the subset's genes,
    # with or without the species held out) are derived from it; without one, the finished database's folder
    # holds the conversion and its build consumes it (the folder serves the taxonomy alone otherwise).
    full = os.path.join(work, "converted")

    def convert_stage(folder):
        """The stage key of the whole release converted into `folder`: protal_db holds a gene subset's files after a
        reduced run, which a later run without the subset must not take for the whole release."""
        return {**convert_key, "into": os.path.basename(folder)}

    def converted_before(folder):
        """Whether `folder` holds the release converted by an earlier run (its stage marked, its files there)."""
        return stages.done("convert", convert_stage(folder)) and os.path.isfile(taxonomy) and \
            all(os.path.isfile(os.path.join(folder, f)) for f in ("reference.fna", "reference.map", "internal_taxonomy.dmp"))

    # The release is converted into protal_db in the background while the genome table is made (they need nothing of
    # each other; 4 and 8-22 min at r226), when it will be converted whatever the table holds: no gene subset, nothing
    # converted there before, and no finished database built from this release and protal (a conversion into protal_db
    # removes such a database's files; whether its gene neighbours still match is known only with the table).
    stored = stages.stored("protal_db")
    could_be_done = isinstance(stored, dict) and {k: v for k, v in stored.items() if k != "gene_neighbours"} == \
        json.loads(json.dumps(final_base)) and os.path.isfile(os.path.join(db, "database.protal"))
    early_convert = None
    if not subset and not could_be_done and not converted_before(db):
        Steps.start(f"converting GTDB r{release} (convert.log), in the background while the genome table is made")
        stages.forget("protal_db")
        early_convert = start_convert(db)

    # The genome table: the genomes that give strain alleles (split_allele_genomes) are no strains to simulate, and left
    # out first, so that only the genomes simulated are read for their lengths; the representatives stay. Their lengths
    # and contigs come from the genome cache (--genome-cache) where it has them, and their contigs go to the reports'
    # cache (genome_contigs.tsv.gz on the samples' disk), so that the reports read no genome again.
    genome_table = args.genome_table or os.path.join(work, "genomes.tsv")
    share = args.allele_genome_share if args.strain_alleles > 0 and args.allele_genome_share > 0 else 0
    if args.genome_cache == "none":
        cache = None
    elif args.genome_cache != "auto":
        cache = GenomeCache(os.path.abspath(args.genome_cache))
    elif args.inputs and os.access(args.inputs, os.W_OK):
        cache = GenomeCache(os.path.join(os.path.abspath(args.inputs), GENOME_CACHE))
    else:
        cache = GenomeCache(os.path.join(work, GENOME_CACHE))
    contig_cache = os.path.join(samples_root, CONTIG_CACHE)
    with concurrent.futures.ThreadPoolExecutor(1) as beside:  # the release's metadata read while the folders are walked
        reading_reps = beside.submit(read_representatives, args.gtdb, release)
        removed = 0
        if not args.genome_table:
            pool = None
            if args.simulate_species:
                with open(args.simulate_species) as fh:
                    pool = {n if n.startswith("s__") else "s__" + n
                            for n in (line.rstrip("\n").split("\t")[0].strip() for line in fh) if n and not n.startswith("#")}
            rows = genome_rows(args.gtdb, release, args.extra_genomes, pool)
            reps = reading_reps.result()
            kept = [r for r in rows if not (share and is_allele_genome(r, reps, share))]
            removed = len(rows) - len(kept)
            write_genome_table(kept, genome_table, args.threads, cache, contig_cache)
        elif args.extra_genomes or args.simulate_species:
            sys.exit("--extra-genomes and --simulate-species shape the genome table this script makes; "
                     "apply them to --genome-table instead")
        else:  # with the genomes' lengths, if it has none: the simulator would read every genome for them
            reps = reading_reps.result()
            rows = table_rows(args.genome_table)
            kept = [r for r in rows if not (share and is_allele_genome(r, reps, share))]
            removed = len(rows) - len(kept)
            if removed:
                genome_table = os.path.join(work, "genomes.tsv")
                write_genome_table([r for r in kept if len(r) >= 3 and os.path.isfile(r[2])], genome_table, args.threads,
                                   cache, contig_cache)
            else:
                genome_table = with_lengths(args.genome_table, os.path.join(work, "genomes.tsv"), args.threads, cache,
                                            contig_cache)
    allele_split = (f"; {removed} genomes left to the strain alleles, not simulated (--allele-genome-share "
                    f"{args.allele_genome_share:g})") if share else ""
    # The share of simulated species that are strains is judged after the in-silico strains (step 3), if any.
    summary, brief, warning = summarize_genome_table(genome_table, reps, insilico_to_come=args.insilico_strains > 0)
    if allele_split:
        summary.append(allele_split[2:])
    with open(os.path.join(logs, "genome_table.txt"), "w") as fh:
        fh.write("\n".join(summary) + "\n")
    Steps.start(f"genome table ({os.path.basename(genome_table)}, genome_table.txt): {brief}{allele_split}")
    if warning:
        say(warning)

    neighbours_key = None if args.no_gene_neighbours else \
        {"convert": convert_key, "script": content_hash(GENE_NEIGHBOURS), "genome_table": content_hash(genome_table)}
    final_key = {**final_base, "gene_neighbours": neighbours_key}
    final_done = stages.done("protal_db", final_key) and os.path.isfile(os.path.join(db, "database.protal"))

    def ensure_converted():
        """With a gene subset: the whole release in `full`, kept from an earlier run that stopped before the
        database folders were derived, or converted now."""
        nonlocal converted
        if converted:
            return
        if stages.done("convert", convert_stage(full)) and os.path.isfile(taxonomy) and os.path.isfile(gene_table) and \
                all(os.path.isfile(os.path.join(full, f)) for f in ("reference.fna", "reference.map", "internal_taxonomy.dmp")):
            clear_build_outputs(full)
            Steps.done(f"{os.path.basename(full)} holds the release converted by an earlier run; not converted again")
        else:
            Steps.done("the whole release: " + convert(full))
            stages.mark("convert", convert_stage(full))
        converted = full

    if subset:
        Steps.start(f"the release (convert.log): converted whole into {os.path.basename(full)}, which the database "
                    "folders of the gene subset are derived from")
        if not os.path.isfile(taxonomy) or not os.path.isfile(gene_table):
            ensure_converted()
        else:
            Steps.done("its taxonomy and gene ids are there from an earlier run; converted again only if a database "
                       "folder has to be derived")
    elif final_done:
        Steps.start(f"the release: {db} was built by an earlier run from the same release and protal; kept" +
                    remove_full_reference(db))  # left by an earlier version of this script
    elif early_convert is None and converted_before(db):
        # Converted by an earlier run that stopped before the build packed the files (not a gene subset's folder,
        # derived into it by a reduced run: that one marks the conversion of `full` only). Its gene neighbours are found
        # again if the genome table changed since.
        clear_build_outputs(db)
        converted = db
        Steps.start(f"the release: {db} holds the release converted by an earlier run; not converted again")
        if neighbours_key and not stages.done("gene_neighbours", {**neighbours_key, "into": os.path.basename(db)}):
            Steps.done("its gene neighbours again, for this genome table: " + neighbours_done(find_neighbours(db).finish(), db))
    else:
        if early_convert is None:
            Steps.start(f"converting GTDB r{release} (convert.log)")
            stages.forget("protal_db")
            early_convert = start_convert(db)
        Steps.done(f"r{release} " + finish_convert(early_convert, db))
        stages.mark("convert", convert_stage(db))
        converted = db
    if not subset and not os.path.isfile(taxonomy):
        Steps.done("its taxonomy: " + convert(full))
        converted = full

    # The species' nearest congeners by their references (species_clouds.tsv): protal --write_species_neighbours on the
    # converted release (what --build stores as species_neighbours.tsv, before any database exists), or --species-clouds
    # copied. The in-silico strains below take each marker gene's congener sites from the nearest congener, and the
    # hold-out keeps or leaves out whole the complexes they show.
    clouds_file = os.path.join(logs, "species_clouds.tsv")
    args.holdout_clouds_note = "none"  # build_metadata.tsv
    clouds = None
    if make_clouds:
        Steps.start("species clouds (species_clouds.tsv, species_clouds.log): every species' nearest congeners by their "
                    "references' marker genes" +
                    (f"; congeners within {args.holdout_complex_distance:g} form a complex the training database leaves "
                     "out or keeps whole" if steer_holdout else "") +
                    ("; the in-silico strains' congener sites" if args.insilico_strains > 0 else ""))
        clouds_key = {"convert": convert_key, "protal": final_key["protal"],
                      "given": content_hash(args.species_clouds) if args.species_clouds else None}
        if args.species_clouds:
            copy_file(args.species_clouds, clouds_file)
            stages.mark("species_clouds", clouds_key)
            Steps.done(f"copied from {args.species_clouds}")
        elif stages.done("species_clouds", clouds_key) and os.path.isfile(clouds_file):
            Steps.done("made by an earlier run from the same release; kept")
        else:
            stages.forget("species_clouds")
            if converted is None:
                if subset:
                    ensure_converted()
                else:
                    converted = full
                    Steps.done("the release again, for its species clouds (the finished database packed its references): " +
                               convert(converted))
            job = run([args.protal, "--write_species_neighbours", clouds_file, "--db", converted, "-t", str(args.threads)],
                      step_log("species_clouds.log"), label="comparing the species' references")
            stages.mark("species_clouds", clouds_key)
            Steps.done(f"made in {job.took()}: {species_neighbours_summary(step_log('species_clouds.log'))}")
        if steer_holdout:
            clouds = read_clouds(clouds_file, taxonomy)
            complexes = species_complexes(clouds, args.holdout_complex_distance)
            args.holdout_clouds_note = (f"{len(clouds)} species' congeners within {CLOUD_MAX_DISTANCE:g} "
                                        f"({'--species-clouds' if args.species_clouds else 'the converted release'}); "
                                        f"{len(set(complexes.values()))} complexes of {len(complexes)} species within "
                                        f"{args.holdout_complex_distance:g} of a congener, held out or kept whole")
            Steps.done(f"{len(clouds)} species compared; {len(set(complexes.values()))} complexes of {len(complexes)} "
                       f"species within {args.holdout_complex_distance:g} of a congener")

    # In-silico strains (insilico_strains.py): every species of the table with one genome gets a mutated copy of
    # its representative (codon-aware substitutions, the divergence of the table's real strains, the genes'
    # conservation factors, as many of its marker substitutions at the sites where its nearest congener differs, with
    # the congener's base, as the real strains have there: species_clouds.tsv), so that it is simulated from a strain
    # as often as a species with two genomes; and --insilico-multi of the species with real strains get one too, a
    # strain their alleles do not reach (half of a species' other genomes are its alleles, the other half its simulated
    # strains, so without these every simulated strain had a close allele). Without
    # them a model given GTDB's cluster sizes learns that a divergent read cloud on a one-genome species is a
    # relative the database lacks (docs/claude/2026-10-04-r226-v10-evaluation). The simulations draw from
    # genomes_simulated.tsv; the species to leave out and the gene neighbours come from the table itself. The strains'
    # FASTAs (one per species: 9,030 at r226 v17) go to the samples' disk, beside the genome store, not to OUTDIR.
    # Nothing but the simulations reads them (the hold-out and the scenarios take genomes.tsv, which has the same
    # species), so they are made in the background, at a lower priority, beside the hold-out, the training database's
    # files and build (2:24 at r226 v19); join_insilico() waits for them before the simulations start.
    sim_table, insilico_note = genome_table, "none (--insilico-strains 0)"
    insilico_job = None

    def insilico_summary(log):
        """The in-silico strains' lines, once made: what insilico_strains.py said, and the genome table simulated from."""
        nonlocal insilico_note
        insilico_note = last_line(log)
        Steps.done(insilico_note)
        summary, brief, warning = summarize_genome_table(sim_table, reps)
        with open(os.path.join(logs, "genome_table.txt"), "a") as fh:
            fh.write("with the in-silico strains (genomes_simulated.tsv, the genomes simulated from):\n" +
                     "\n".join(summary[1:]) + "\n" + insilico_note + "\n")
        Steps.done(f"simulated from: {brief}")
        if warning:
            say(warning)

    def join_insilico():
        """Waits for the in-silico strains made in the background, if they are, and says how it went."""
        nonlocal insilico_job
        if insilico_job is None:
            return
        job, insilico_job = insilico_job, None
        if job.seconds is None:
            Steps.done("waiting for the in-silico strains (insilico_strains.log), which the simulations draw from")
        job.finish()
        stages.mark("insilico", insilico_key)
        Steps.done(f"in-silico strains made in the background in {job.took()}" +
                   ("" if insilico_positions else f" (no gene_positions.tsv: ANI {args.insilico_ani or INSILICO_FALLBACK_ANI})"))
        insilico_summary(job.log)

    if args.insilico_strains > 0:
        sim_table = os.path.join(work, "genomes_simulated.tsv")
        strains = os.path.join(samples_root, "insilico_strains")
        log = step_log("insilico_strains.log")
        Steps.start("in-silico strains of the species with one genome" +
                    (f" and of {args.insilico_multi:g} of those with real strains" if args.insilico_multi > 0 else "") +
                    " (insilico_strains.log, genomes_simulated.tsv)")
        insilico_key = {"script": content_hash(INSILICO), "genome_table": content_hash(genome_table),
                        "convert": convert_key, "gene_neighbours": neighbours_key, "share": args.insilico_strains,
                        "multi": args.insilico_multi, "ani": args.insilico_ani, "seed": args.seed, "strains": strains,
                        "clouds": content_hash(clouds_file) if make_clouds and os.path.isfile(clouds_file) else None}
        # A new --scratch (the next job's node) has none of them: made again (the same strains, at the same seed).
        if stages.done("insilico", insilico_key) and os.path.isfile(sim_table) and \
                os.path.isfile(os.path.join(strains, "insilico_strains.tsv")):
            Steps.done("made by an earlier run from the same table; kept")
            insilico_summary(log)
        else:
            stages.forget("insilico")

            def positions_file():
                return next((p for p in (os.path.join(f, "gene_positions.tsv") for f in (converted, db, full) if f)
                             if os.path.isfile(p)), None)
            if positions_file() is None and not args.no_gene_neighbours:
                # The strains take each gene's divergence from the real strains' (gene_positions.tsv): without it they
                # would all get a uniform ANI. A finished database kept from an earlier run packed its positions, so
                # the release is converted again here (the training database's files are derived from it below).
                if subset:
                    ensure_converted()
                elif converted is None:
                    converted = full
                    Steps.done("the release again, for its gene positions (the finished database packed them): " +
                               convert(converted))
            positions = positions_file()
            command = [sys.executable, INSILICO, "--genome-table", genome_table, "--output", sim_table, "--out-dir",
                       strains, "--share", str(args.insilico_strains), "--multi-share", str(args.insilico_multi),
                       "--seed", str(args.seed), "-t", str(args.threads)]
            if positions:
                command += ["--positions", positions, "--taxonomy", taxonomy]
                if make_clouds and os.path.isfile(clouds_file):
                    command += ["--clouds", clouds_file]
            if args.insilico_ani or not positions:
                command += ["--ani", args.insilico_ani or INSILICO_FALLBACK_ANI]
            insilico_positions = positions
            insilico_job = Job(command, log, label="making in-silico strains", nice=10)
            Steps.done("made in the background, beside the hold-out and the training database; the simulations wait for "
                       "them")

    # The training database leaves some species out: the model then sees reads of species the database
    # lacks, which land on relatives, and reads of whole families, classes and phyla it lacks, which land on
    # distant ones. It is made from the converted files before --build packs them. heldout_species.txt: the
    # species, the rank they were held out at and the clade.
    heldout = os.path.join(logs, "heldout_species.txt")
    if args.holdout_species:
        copy_file(args.holdout_species, heldout)
    elif args.holdout > 0 or clades:
        chosen = choose_holdout(genome_table, taxonomy, args.holdout, clades, args.holdout_max_share, args.seed, clouds,
                                args.holdout_complex_distance)
        write_holdout(heldout, chosen, clouds, args.holdout_complex_distance)
    elif os.path.exists(heldout):
        os.remove(heldout)
    for note in scenario_notes:
        say(f"    WARNING: scenario {note}")
    if selected:
        # A scenario's share of species the database lacks comes from the species held out: enough of them, and of
        # the others, for its largest sample (scenarios.table_split), or the scenario is scaled down to what they hold
        # (scenarios.fit_species), as its collections do it.
        novel = read_holdout(heldout) if os.path.exists(heldout) else {}
        pool = pool_species(genome_table)  # the species of genomes_simulated.tsv, whose in-silico strains may still be made
        known, lacking = sum(s not in novel for s in pool), sum(s in novel for s in pool)
        notes = []
        for name in selected:
            try:
                d, scaled = scenarios.fit_species(scenario_defs[name], known, lacking)
            except scenarios.ScenarioError as e:
                sys.exit(f"scenario {name}: {e}")
            k, n = scenarios.table_split(known, lacking, d["novel_share"])
            notes.append(f"{name} {k + n} species ({n} the training database lacks) for samples of {d['species']}")
            if scaled:
                scenario_notes.append(f"{name} {scaled}")
                say(f"    WARNING: scenario {name} {scaled}")
        say(f"    scenarios: {', '.join(f'{name} {hold_in.get(name, 0)} hold-in and {hold_out.get(name, 0)} hold-out samples' for name in selected)}"
            f"; their genome tables: {', '.join(notes)}")
    args.scenario_notes = scenario_notes  # build_metadata.tsv
    # The marker genes of the databases: all of them, or a subset (--n-genes: the N most distinctive by
    # scripts/rank_genes.py, from a full build of the training database or --gene-ranking; --genes: the ones
    # named), listed in gene_subset.txt, which the folders are derived with.
    subset_file = os.path.join(logs, "gene_subset.txt")
    ranking_file = os.path.join(logs, "gene_ranking.tsv")
    genes_note = "all"  # build_metadata.tsv: which marker genes, chosen how
    if subset:
        gene_ids = read_gene_ids(gene_table)
        names = {gid: marker for marker, gid in gene_ids.items()}
        if args.genes:
            Steps.start("marker genes (gene_subset.txt): the ones listed (--genes)")
            chosen = read_gene_list(args.genes, gene_ids)
            if len(chosen) >= len(gene_ids):
                sys.exit(f"--genes names every one of the release's {len(gene_ids)} marker genes: no subset")
            with open(subset_file + ".partial", "w", newline="\n") as fh:
                fh.write(f"# {len(chosen)} of the {len(gene_ids)} marker genes of GTDB r{release} (--genes {args.genes})\n")
                for gid in chosen:
                    fh.write(f"# {names[gid]}\n{gid}\n")
            os.replace(subset_file + ".partial", subset_file)
            listed = ", ".join(names[g] for g in chosen)
            genes_note = f"{len(chosen)} of {len(gene_ids)} (--genes): {listed}"
            Steps.done(f"{len(chosen)} of the {len(gene_ids)} marker genes: {listed}")
        else:
            n = args.n_genes
            Steps.start(f"marker genes (gene_subset.txt): the {n} most distinctive by prevalence x unique k-mer share")
            if n >= len(gene_ids):
                sys.exit(f"--n-genes {n}: GTDB r{release} has {len(gene_ids)} marker genes; a subset has fewer")
            if args.gene_ranking:
                copy_file(args.gene_ranking, ranking_file)
                source = f"--gene-ranking {args.gene_ranking}"
                Steps.done(f"ranked by {args.gene_ranking} (copied to gene_ranking.tsv)")
            else:
                # A full build of the training database (every gene, the species held out left out; separate files,
                # as rank_genes.py reads them) ranks the genes: its unique k-mer table says how many of each gene's
                # k-mers name their species. The build costs what the training database's does; the folder goes
                # once ranked, the ranking stays for a rerun.
                ranking_key = {"convert": convert_key, "protal": final_key["protal"], "ranker": content_hash(RANKER),
                               "heldout": content_hash(heldout) if os.path.exists(heldout) else None}
                if stages.done("gene_ranking", ranking_key) and os.path.isfile(ranking_file):
                    Steps.done("ranked by an earlier run from a full build of the training database; kept (gene_ranking.tsv)")
                else:
                    stages.forget("gene_ranking")
                    ensure_converted()
                    ranking_db = os.path.join(samples_root, "ranking_db")
                    shutil.rmtree(ranking_db, ignore_errors=True)
                    command = [sys.executable, CONVERTER, "--from_db", converted, "--outdir", ranking_db, "-t", str(args.threads)]
                    if os.path.exists(heldout):
                        command += ["--exclude_species", heldout]
                    run(command, step_log("gene_ranking_files.log"), label="writing the ranking database's files")
                    job = run(build_command(args.protal, ranking_db, args.threads, "--compress_level", str(args.training_db_level),
                                            "--no_bundle"),
                              step_log("gene_ranking_build.log"), label="building the ranking database")
                    rows = rank_genes.rank(ranking_db, gene_table, taxonomy)
                    rank_genes.write_table(ranking_file, rows)
                    shutil.rmtree(ranking_db, ignore_errors=True)
                    stages.mark("gene_ranking", ranking_key)
                    Steps.done(f"ranked the {len(rows)} genes from a full build of the training database"
                               f"{' (every gene, ' + str(len(read_holdout(heldout))) + ' species left out)' if os.path.exists(heldout) else ''}"
                               f", built in {job.took()} (gene_ranking_build.log): gene_ranking.tsv")
                source = "a full build of the training database"
            rows = rank_genes.read_ranking(ranking_file)
            unknown = [str(r["gene_id"]) for r in rows if r["gene_id"] not in names]
            if unknown:
                sys.exit(f"{ranking_file} ranks genes this release lacks ({', '.join(unknown[:5])}): a ranking of another release?")
            if n > len(rows):
                sys.exit(f"--n-genes {n}: {ranking_file} ranks {len(rows)} genes")
            # Each domain keeps its best genes (a gene of the archaeal marker set alone scores low over all species).
            chosen = rank_genes.select(rows, n, args.genes_per_domain)
            rank_genes.write_subset(subset_file, chosen, len(rows), source)
            listed = ", ".join(r["marker"] for r in chosen)
            covered = rank_genes.coverage(chosen)
            genes_note = (f"{n} of {len(gene_ids)}, the most distinctive by prevalence x unique k-mer share "
                          f"(scripts/rank_genes.py, {args.genes_per_domain or rank_genes.per_domain_default(n)} per domain, "
                          f"from {source}): {listed}" +
                          (f"; in half the species or more of {', '.join(f'{d} {c}' for d, c in covered.items())}" if covered else ""))
            Steps.done(f"the {n} best of {len(rows)}: {rank_genes.describe(chosen)}")
        final_key["genes"] = content_hash(subset_file)
        final_done = stages.done("protal_db", final_key) and os.path.isfile(os.path.join(db, "database.protal"))
        if final_done:
            Steps.done(f"{os.path.basename(db)} was built by an earlier run from these genes; kept" + remove_full_reference(db))
        else:
            ensure_converted()
            stages.forget("protal_db")
            job = run([sys.executable, CONVERTER, "--from_db", converted, "--genes", subset_file, "--outdir", db,
                       "-t", str(args.threads)], step_log("protal_db_files.log"),
                      label=f"deriving {os.path.basename(db)}'s files")
            Steps.done(f"{os.path.basename(db)}'s files derived for these genes in {job.took()} (protal_db_files.log)")

    training_db, training_done = db, False
    n_heldout, files_took = 0, ""
    keep_full = False  # the training database's full reference kept for the ancestry report (reports)
    if os.path.exists(heldout):
        chosen = read_holdout(heldout)
        n_heldout = len(chosen)
        with open(os.path.join(logs, "holdout.txt"), "w") as fh:
            fh.write("\n".join([f"training database: {n_heldout} species left out ({heldout})"] +
                               describe_holdout(chosen, pool_species(genome_table), read_holdout_details(heldout),
                                                args.holdout_complex_distance if clouds else 0)) + "\n")
        # Read only by the collections and the parity check: on --scratch, its build writes and they load it from
        # the node's disk (writing database.protal to a network file system was 6 of the 15 min of an r226 build).
        training_db = os.path.join(samples_root, "training_db")
        # Its full reference stays for the foreign scan right after its build (foreign_rates) and, with --share-logs, for
        # the ancestry report (ancestry_sites.py reads the species' other genomes from it, after the error reads).
        keep_full = args.share_logs or not args.no_foreign_rates
        Steps.start(f"training database ({os.path.basename(training_db)}, training_db_index.log): {n_heldout} species left "
                    f"out, {holdout_brief(chosen)} (model_logs/holdout.txt)")
        training_key = {"convert": convert_key, "gene_neighbours": neighbours_key, "heldout": content_hash(heldout),
                        "protal": final_key["protal"], "level": args.training_db_level, "alleles": ALLELE_ARGS}
        if subset:
            training_key["genes"] = final_key["genes"]
        # The key is kept in the training database's folder too: on a --scratch that several OUTDIRs share, another
        # build may have replaced that folder (other species held out) since this OUTDIR's stage was marked.
        training_stamp = Stages(training_db)
        training_done = stages.done("training_db", training_key) and training_stamp.done("built_for", training_key) and \
            os.path.isfile(os.path.join(training_db, "database.protal"))
        if training_done:
            Steps.done(f"{training_db} was built by an earlier run with the same species left out"
                       f"{' and the same genes' if subset else ''}; kept" + full_reference_fate(training_db, keep_full))
        else:
            stages.forget("training_db")
            training_stamp.forget("built_for")
            if subset:
                ensure_converted()
            elif converted is None:  # the finished database's build consumed them
                converted = full
                Steps.done("the release again (the finished database's build took its files): " + convert(converted))
            job = run([sys.executable, CONVERTER, "--from_db", converted, "--exclude_species", heldout, "--outdir",
                       training_db, "-t", str(args.threads)] + (["--genes", subset_file] if subset else []),
                      step_log("training_db.log"),
                      label="leaving the species out")
            files_took = f"; its files written in {job.took()} (training_db.log)"
            if args.share_logs:
                # The ancestry report (--share-logs) reads the training database's reference.fna, which its build packs
                # and removes: a hard link keeps it, so that the report needs no protal --unpack_db (1:59 at r226 v19).
                kept = os.path.join(training_db, ANCESTRY_REFERENCE)
                with contextlib.suppress(OSError):
                    os.remove(kept)
                with contextlib.suppress(OSError):
                    os.link(os.path.join(training_db, "reference.fna"), kept)
    if converted and converted != db:
        join_insilico()  # they may read the gene positions of the release converted again, which goes now
        shutil.rmtree(converted, ignore_errors=True)

    # The finished database is needed only for --add_model at the end: it is built in the background from the
    # start, while the training database is built and the training data collected, unless one build at a time.
    final_log = step_log("index_and_package.log")
    final_build = None
    final_level = ("--compress_level", str(args.final_db_level))

    def built(name, job, extra=""):
        """The line of a build that ended."""
        return (f"built {os.path.basename(name)}{' in the background' if job is final_build else ''} in "
                f"{job.took()}{db_size(name)}{extra}")

    def built_final():
        stages.mark("protal_db", final_key)

    def built_final_in_background():  # says so when the script notices, whatever step it is at
        built_final()
        Steps.done(built(db, final_build, full_reference_fate(db, scan, "the foreign scan")))

    # Training data of every read type (pe, se from its first reads, pb and ont from long reads of the same
    # communities), then an independent test set of another design, both profiled against the training database.
    def collect_command(out, samples, read_pairs, species, abundance, strains, long_bases, long_samples, seed,
                        scenario_samples=None):
        """The collector's command; scenario_samples: {scenario: samples} collected beside the design."""
        command = [sys.executable, COLLECTOR, "--db", training_db, "--genome_table", sim_table, "-o", out,
                   "--protal", args.protal, "--simulator", args.simulator, "--samples", str(samples),
                   "--long_read_samples", str(long_samples),
                   "--read_pairs", read_pairs, "--read_setups", args.read_setups, "--archaea", str(args.archaea),
                   "--species_per_sample", species, "--seed", str(seed), "-t", str(args.threads),
                   "--taxonomy", taxonomy, "--congeners", congener_text(args.congeners), "--read_types", ",".join(read_types),
                   "--long_read_bases", long_bases, "--pb_setup", args.pb_setup, "--ont_setup", args.ont_setup,
                   "--pbsim", args.pbsim, "--read_compression", args.read_compression]
        command += ["--abundance", abundance] if abundance else []
        command += ["--strains_per_species", strains] if strains else []
        command += ["--pbsim_models", args.pbsim_models] if args.pbsim_models else []
        if n_heldout:
            command += ["--novel_species", heldout, "--novel_clades", str(args.novel_clades_per_sample)]
            if os.path.isfile(clouds_file):  # meta_novel_distance: how far an absent taxon's held-out congener is
                command += ["--species_clouds", clouds_file]
        if scenario_samples:
            command += ["--scenarios", ",".join(f"{name}:{n}" for name, n in scenario_samples.items())]
            command += ["--scenario_file", os.path.abspath(args.scenario_file)] if args.scenario_file else []
            command += ["--host_genome", args.host_genome] if args.host_genome else []
        if args.error_units:
            command += ["--unmapped_reads", ",".join(kind + (f":{scope}" if scope else "") for kind, scope in args.error_units)]
        if genome_store:  # both collections share it
            command += ["--genome_store", genome_store]
        if args.compressed_pipes:
            command += ["--compressed_pipes"]
        if args.profile_ahead:
            command += ["--profile_ahead"]
        return command

    genome_store = None  # --genome-store
    if args.genome_store == "auto":
        genome_store = os.path.join(samples_root, "genome_store")
    elif args.genome_store != "none":
        genome_store = os.path.abspath(args.genome_store)
    training = os.path.join(samples_root, "training")
    test = os.path.join(samples_root, "test")
    collections_ = [("training data", collect_command(training, args.samples, args.read_pairs, args.species_per_sample,
                                                      args.abundance, args.strains_per_species, args.long_read_bases,
                                                      args.long_read_samples, args.seed, hold_in),
                     step_log("training_data.log"))]
    if has_test:  # the scenarios' hold-out samples are of another seed, as the test set's design is
        collections_.append(("independent test set" if args.test_samples > 0 else "scenarios' hold-out samples",
                             collect_command(test, args.test_samples, args.test_read_pairs, args.test_species_per_sample,
                                             args.test_abundance, args.test_strains_per_species,
                                             args.test_long_read_bases, args.test_long_read_samples, args.seed + 1000,
                                             hold_out),
                             step_log("test_data.log")))

    if training_db == db:
        Steps.start(f"database ({os.path.basename(db)}, index_and_package.log): no species left out" +
                    ("; built by an earlier run, kept" if final_done else ""))
    # What the run needs on the samples' disk besides the reads (disk_needs), and what the reads take, said before the
    # builds; with --stream-above auto the run stops here if even streaming every simulation would not fit. What to
    # stream is chosen once the training database is built, from the space free then.
    sizes = simulation_sizes(collections_)
    reads = sum(r for _, r in sizes)
    hosts = [(o.host_genome, os.path.join(o.out, "host")) for o in (collector_args(c[2:]) for _, c, _ in collections_)]
    keep_free = args.keep_free * 1e9 if args.profile_blocks > 0 else 0
    to_build = ([] if final_done else [("database" if training_db == db else "finished database", db)]) + \
        ([("training database", training_db)] if training_db != db and not training_done else [])
    # The in-silico strains may still be made: their genomes counted at their representatives' lengths.
    store_bases = planned_bases(genome_table, args.insilico_strains) if insilico_job is not None else sim_table
    needs = disk_needs(samples_root, genome_store, store_bases, to_build, hosts, sizes, keep_free)
    free = shutil.disk_usage(samples_root).free
    say(f"    room on {samples_root}: {gigabytes(free)} free; besides the reads the run needs ~"
        f"{gigabytes(sum(needs.values()))} there ({needs_text(needs)}), and the reads of its {len(sizes)} simulations take "
        f"~{gigabytes(reads)} (the collector's estimate, ~1.25 times what r226 v15 wrote)" +
        ("; what to stream into protal is chosen once the training database is built (--stream-above auto)"
         if args.profile_blocks > 0 and args.stream_above == "auto" else ""))
    if args.profile_blocks > 0 and args.stream_above == "auto" and free < sum(needs.values()):
        stop(f"Not enough room on {samples_root}: {gigabytes(free)} free, and the run needs ~{gigabytes(sum(needs.values()))} "
             f"there even with every simulation streamed into protal ({needs_text(needs)}): give it a larger disk "
             "(--scratch), fewer samples or scenarios or a lower --keep-free, or --stream-above GB to run anyway")
    if args.profile_blocks <= 0 and free < sum(needs.values()) + reads:
        say(f"    WARNING: with --profile-blocks 0 every read is on the disk at once: ~{gigabytes(sum(needs.values()) + reads)} "
            f"needed, {gigabytes(free)} free")

    # The training database gates the profiling, which is the build's longest stage now that the simulations take minutes
    # (docs/claude/2026-10-07-build-ordering): it is built alone, with every core and no simulation writing to its disk.
    # The finished database is needed only for --add_model at the end: built after it, at the idle scheduling class,
    # on the cores the profiling leaves (and paused while the models are trained).
    # The scan of a database's full reference for the gene copies other species' reads reach (scripts/foreign_rates.py):
    # every species' marker genes tiled alike (a few genomes each), aligned against the database, the table stored in it
    # (--add_tables). The training database right after its build, from its full reference, which lacks the held-out
    # species (so the table knows nothing of them); the finished database after its own build, from its full reference
    # (every species), before its models go in. A database's full reference goes once its scan is done (the training
    # database's with --share-logs after the ancestry report). The finished database's table waits in work/foreign_rates/
    # for a rerun, not beside database.protal, which holds it. Only with --foreign-rates (r226 v18: 0.001 of AUC for ~21 min).
    scan = not args.no_foreign_rates

    def foreign_rates(against):
        table = os.path.join(against if against != db else os.path.join(work, "foreign_rates"), "foreign_rates.tsv")
        stage = "foreign_rates" if against == db else "foreign_rates_training"
        key = {"database": final_key if against == db else training_key, "stride": args.foreign_stride,
               "per_header": args.foreign_per_header, "script": content_hash(FOREIGN_RATES)}
        log = step_log(stage + ".log")
        Steps.start(f"the gene copies' foreign reads of {os.path.basename(against)} ({os.path.basename(log)}): reads every "
                    f"{args.foreign_stride} bases of {'every' if args.foreign_per_header == 0 else 'at most ' + str(args.foreign_per_header)} "
                    f"cop{'y' if args.foreign_per_header == 1 else 'ies'} of each species' gene in its full reference, aligned "
                    f"against it")
        if stages.done(stage, key) and os.path.isfile(table):
            Steps.done(f"{table} was made by an earlier run for the same database; kept")
            return table
        full = full_reference_path(against)
        if not full:
            Steps.done(f"left out: {os.path.basename(against)} has no full reference (removed by an earlier run?); the "
                       "'foreign' features are unknown (-1)")
            return None
        if not os.path.isfile(taxonomy):
            Steps.done(f"left out: no taxonomy at {taxonomy}; the 'foreign' features are unknown (-1)")
            return None
        stages.forget(stage)
        command = [sys.executable, FOREIGN_RATES, "--db", against, "--full-reference", full, "--taxonomy", taxonomy,
                   "--out", table, "--stride", str(args.foreign_stride), "--per-header", str(args.foreign_per_header),
                   "-t", str(args.threads), "--protal", args.protal, "--simulate", args.simulator,
                   "--workdir", os.path.join(samples_root, "foreign_rates_scan")]
        os.makedirs(os.path.dirname(table), exist_ok=True)
        job = run(command, log, label="scanning the full reference")
        said = last_line(log)
        run([args.protal, "--add_tables", table, "--db", against, "-t", str(args.threads)], log,
            lambda: stages.mark(stage, key), f"storing the table in {os.path.basename(against)}", append=True)
        Steps.done(f"{said}; in {job.took()}; stored in {os.path.basename(against)}")
        return table

    def scanned(folder, keep=False):
        """The scan of a database just built, then its full reference removed (unless keep); the line to add to the
        build's."""
        if not scan:
            return remove_full_reference(folder)
        foreign_rates(folder)
        return "" if keep else remove_full_reference(folder)

    if final_done:
        pass
    elif training_db == db:
        job = run(build_command(args.protal, db, args.threads, *final_level), final_log, built_final,
                  f"building {os.path.basename(db)}")
        Steps.done(built(db, job, full_reference_fate(db, scan, "the foreign scan")))
        if scan:
            removed = scanned(db)
            if removed:
                Steps.done(f"{os.path.basename(db)}'s full reference, read by the scan{removed}")
    if training_db != db and not training_done:
        # Read only for the training samples and the parity check: zstd level 3 packs it in a fraction of the
        # time of level 19 (which half of a build spent on), and loads as fast.
        job = run(build_command(args.protal, training_db, args.threads, "--compress_level", str(args.training_db_level)),
                  step_log("training_db_index.log"),
                  lambda: (stages.mark("training_db", training_key), training_stamp.mark("built_for", training_key)),
                  f"building {os.path.basename(training_db)}")
        Steps.done(built(training_db, job, full_reference_fate(training_db, keep_full, "the foreign scan" if scan else
                                                                "the ancestry report") + files_took))
    if scan and training_db != db:
        removed = scanned(training_db, keep=args.share_logs)
        if removed:
            Steps.done(f"{os.path.basename(training_db)}'s full reference, read by the scan{removed}")
    if not final_done and training_db != db:
        if args.one_build_at_a_time:
            Steps.done(f"{os.path.basename(db)} is built after the models (--one-build-at-a-time)")
        else:
            final_build = Job(build_command(args.protal, db, args.threads, *final_level), final_log,
                              built_final_in_background, label=f"building {os.path.basename(db)} in the background",
                              idle=True)
            Steps.done(f"building {os.path.basename(db)} meanwhile, in the background at the idle scheduling class: on "
                       "the cores the profiling leaves (index_and_package.log)")

    # What to stream into protal (--stream-above): with auto, the fewest simulations whose streaming leaves the others'
    # reads room on the disk at once, besides what disk_needs counts now that the training database is built. The
    # simulations draw from the in-silico strains: they are made by now (or waited for).
    join_insilico()
    free = shutil.disk_usage(samples_root).free
    needs = disk_needs(samples_root, genome_store, sim_table,
                       [("finished database", db)] if not final_done and training_db != db else [], hosts, sizes,
                       keep_free)
    ancestry_reference = os.path.join(training_db, ANCESTRY_REFERENCE)  # kept for the ancestry report until the end
    if os.path.isfile(ancestry_reference) and on_disk(training_db, samples_root):
        needs["the ancestry report's reference.fna"] = os.path.getsize(ancestry_reference)
    room = free - sum(needs.values())
    auto = args.stream_above == "auto"
    if args.profile_blocks <= 0:
        stream_above = 0.0
    elif auto:
        threshold = stream_threshold(sizes, room)
        if threshold is None:
            stop(f"Not enough room on {samples_root}: {gigabytes(free)} free, and the run needs ~"
                 f"{gigabytes(sum(needs.values()))} there even with every simulation streamed into protal "
                 f"({needs_text(needs)}): give it a larger disk (--scratch), fewer samples or scenarios or a lower "
                 "--keep-free, or --stream-above GB to run anyway")
        stream_above = threshold / 1e9
    else:
        stream_above = args.stream_above
    streamed = [r for largest, r in sizes if stream_above > 0 and largest > stream_above * 1e9]
    written = reads - sum(streamed)
    chosen = (f"{len(streamed)} of the {len(sizes)} simulations streamed into protal (a sample above "
              f"{stream_above:.3g} GB; ~{gigabytes(sum(streamed))} of reads), ~{gigabytes(written)} written"
              if streamed else f"none of the {len(sizes)} simulations streamed: their reads written (~{gigabytes(written)})")
    how = "--profile-blocks 0" if args.profile_blocks <= 0 else ("--stream-above auto" if auto else
                                                                 f"--stream-above {stream_above:g}")
    say(f"    {how}: {chosen}; {gigabytes(free)} free on {samples_root}, ~{gigabytes(sum(needs.values()))} of it for the "
        f"run besides the reads ({needs_text(needs)})")
    if not auto and args.profile_blocks > 0 and written > room:
        say(f"    WARNING: the reads written (~{gigabytes(written)}) do not fit at once in the room left "
            f"(~{gigabytes(max(0, room))}): the simulations will wait for reads to be profiled and removed")
    args.streaming = f"{how}: {chosen}; {gigabytes(free)} free on the samples' disk once the training database was built"
    disk_estimate = sum(needs.values()) + written  # at most, all at once; said again beside the peak at the end

    # The simulations: both collections from here on, in the background at a lower priority than the profiling
    # (collect_training_data.py --simulate_only), profiled as they go (follow(): the reads profiled removed, so that the
    # simulations keep --keep-free GB free and wait for that while the disk is full), or once all are simulated
    # (--profile-blocks 0, collect()).
    simulations = {}
    for what, command, log in collections_:
        options = ["--min_free", f"{args.keep_free:g}"] if args.profile_blocks > 0 and args.keep_free > 0 else []
        if stream_above > 0:  # the follower streams them (the same value for both)
            options += ["--stream_above", repr(stream_above)]
        simulations[what] = Job(command + ["--simulate_only"] + options, log.removesuffix(".log") + "_simulation.log",
                                lambda what=what: Steps.done(f"simulated the {what} in the background in "
                                                             f"{simulations[what].took()}"),
                                label=f"simulating the {what}", nice=10)
    Steps.done(f"simulating the {' and the '.join(simulations)} in the background "
               f"({', '.join(os.path.basename(j.log) for j in simulations.values())})")
    if args.rank_genes:
        # The genes ranked from the training database (every gene, the species held out left out), unpacked on
        # the scratch disk for rank_genes.py: a reduced database of this release takes the ranking
        # (--gene-ranking) instead of building one.
        if training_db == db:
            sys.exit("--rank-genes ranks from the training database: hold species out")
        ranking_key = {"training": training_key, "ranker": content_hash(RANKER)}
        if stages.done("gene_ranking", ranking_key) and os.path.isfile(ranking_file):
            Steps.done("the genes ranked by an earlier run from the training database; kept (gene_ranking.tsv)")
        else:
            stages.forget("gene_ranking")
            unpacked = os.path.join(samples_root, "ranking_files")
            shutil.rmtree(unpacked, ignore_errors=True)
            job = run([args.protal, "--unpack_db", "--db", os.path.join(training_db, "database.protal"), "--unpack_dir",
                       unpacked, "-t", str(args.threads)], step_log("gene_ranking.log"),
                      label="unpacking the training database to rank its genes")
            # gene_congeners.tsv is a report of the build beside database.protal, not a member of it.
            congeners = os.path.join(training_db, "gene_congeners.tsv")
            if os.path.isfile(congeners):
                shutil.copy(congeners, unpacked)
            rows = rank_genes.rank(unpacked, gene_table, taxonomy)
            rank_genes.write_table(ranking_file, rows)
            shutil.rmtree(unpacked, ignore_errors=True)
            stages.mark("gene_ranking", ranking_key)
            twelve = rank_genes.select(rows, min(12, len(rows)))
            Steps.done(f"ranked the {len(rows)} genes of {os.path.basename(training_db)} in {job.took()} (--rank-genes, "
                       f"gene_ranking.log): gene_ranking.tsv; the {len(twelve)} a reduced database would take: "
                       + rank_genes.describe(twelve))

    def announce(what, command, log):
        """The step of a collection: what it makes. -> the collector's options."""
        opts = collector_args(command[2:])
        samples, in_scenarios = collections.Counter(), collections.Counter()
        for unit in units_of(opts)[1]:
            samples[unit["type"]] += unit["samples"]
            in_scenarios[unit["type"]] += unit["samples"] if unit.get("scenario") else 0
        Steps.start(f"{what} ({os.path.basename(log)}): {', '.join(f'{n} {t}' for t, n in samples.items())} samples" +
                    (f" (of them in the scenarios {opts.scenarios}: "
                     f"{', '.join(f'{n} {t}' for t, n in in_scenarios.items() if n)})" if +in_scenarios else ""))
        return opts

    def report(opts, log, job):
        """What a collection's tables hold (present and absent taxa per read type) and how much space its samples
        take. With --scratch, its tables are copied to OUTDIR/work/<collection>."""
        type_of = {name: t for t, name in TABLES.items()}
        counts = []
        with open(log, errors="replace") as fh:
            for line in fh:
                m = re.match(r"\d+ taxa in (\S+): (\d+) present, (\d+) absent", line)
                if m:
                    counts.append(f"{type_of.get(os.path.basename(m.group(1)), m.group(1))} {m.group(2)}/{m.group(3)}")
        Steps.done(f"collected in {job.took()}; taxa present/absent: {', '.join(counts) or 'none'}; "
                   f"{gigabytes(tree_size(opts.out))} in {opts.out}" + (f"; scratch: {Job.scratch.text()}" if Job.scratch else ""))
        if args.scratch:
            keep = os.path.join(work, os.path.basename(opts.out))
            os.makedirs(keep, exist_ok=True)
            for table in TABLES.values():
                if os.path.isfile(os.path.join(opts.out, table)):
                    shutil.copy(os.path.join(opts.out, table), keep)
            # Each protal run's log (its stage timers: how long the alignment and the profiling stage took), which a
            # node's disk cleared after the job would lose, one after the other in logs/protal_runs_<collection>.log.
            def run_name(path):
                name = os.path.relpath(os.path.dirname(path), os.path.join(opts.out, "profile_all"))
                return "all" if name == "." else name.replace(os.sep, "_")
            runs = sorted(glob.glob(os.path.join(opts.out, "profile_all", "**", "protal.log"), recursive=True),
                          key=lambda p: [int(x) if x.isdigit() else x for x in re.split(r"(\d+)", run_name(p))])
            if runs:
                with open(step_log(f"protal_runs_{os.path.basename(opts.out)}.log"), "w") as out:
                    for path in runs:
                        out.write(f"==> {run_name(path)} <==\n")
                        with open(path, errors="replace") as fh:
                            shutil.copyfileobj(fh, out)
        # The collector's CPU table of its protal runs (each run, and protal's own stages and samples: misc/cpu.tsv) as
        # logs/protal_cpu_<collection>.tsv, beside the run's other CPU logs.
        cpu = os.path.join(opts.out, PROTAL_CPU)
        if os.path.isfile(cpu):
            shutil.copyfile(cpu, step_log(f"protal_cpu_{os.path.basename(opts.out)}.tsv"))

    def collect(what, command, log):
        """Runs the collector once the collection's simulations are done (in the background), so that it profiles
        them (both collections in one protal run)."""
        opts = announce(what, command, log)
        if simulations[what].seconds is None:
            Steps.done(f"waiting for its simulations in the background ({os.path.basename(simulations[what].log)})")
            simulations[what].finish()  # its line comes from its on_success
        if what == collections_[0][0] and len(collections_) > 1:
            # The other collection's samples go into this collection's protal run, which loads the database once:
            # the other collector maps them (--prepare_profiling), and finds them profiled when it runs.
            other, other_command, other_log = collections_[1]
            if simulations[other].seconds is None:
                Steps.done(f"waiting for the {other}'s simulations too, to profile both in one protal run")
                simulations[other].finish()
            run(other_command + ["--prepare_profiling"], other_log.removesuffix(".log") + "_map.log",
                label=f"mapping the {other}'s samples")
            other_map = os.path.join(collector_args(other_command[2:]).out, "profile_all", "samples.map")
            if os.path.isfile(other_map):
                command = command + ["--also_profile", other_map]
                Steps.done(f"profiling the {other}'s samples in the same protal run")
        report(opts, log, run(command, log, label=f"collecting {what}"))

    def follow_all():
        """--profile-blocks: a collector per collection profiles its samples as they are simulated
        (collect_training_data.py --follow), in protal runs that take turns (a lock on the scratch disk), removing the
        reads profiled, then writes the collection's tables."""
        lock = os.path.join(samples_root, "protal.lock")
        followers = {}
        for what, command, log in collections_:
            out = collector_args(command[2:]).out
            # Its simulations must have said that they run and have prepared what they write before simulating (or
            # have ended): a follower that finds neither simulates what is left itself, which it must not do beside
            # them.
            while simulations[what].seconds is None and not (simulation_state(out) or {}).get("prepared"):
                check_jobs()
                time.sleep(0.5)
            stream = ["--stream_above", repr(stream_above)] if stream_above > 0 else []
            followers[what] = Job(command + ["--follow", "--profile_block", f"{args.profile_blocks:g}", "--protal_lock", lock,
                                             "--profile_block_max", f"{args.profile_block_max:g}"] + stream,
                                  log, label=f"profiling the {what}")
        for what, command, log in collections_:
            opts = announce(what, command, log)
            Steps.done(f"profiled as its simulations go on, in protal runs of {args.profile_blocks:g} GB of reads or "
                       "more" + (f", at most {args.profile_block_max:g}" if args.profile_block_max > 0 else "") +
                       " (taking turns with the other collection's); each design point's reads removed once profiled")
            report(opts, log, followers[what].finish())

    if args.profile_blocks > 0:
        follow_all()
    else:
        for what, command, log in collections_:
            collect(what, command, log)

    # One model per read type, trained in parallel, its outputs in model_logs/ (the model as it goes into the database,
    # the report, the numbers, the predictions and calls).
    prefixes = {t: os.path.join(logs, "trained_model" + ("" if t == "pe" else "_" + t)) for t in read_types}
    trainer_threads = max(1, args.threads // len(read_types))
    models = f"the {', '.join(read_types)} model{'s' if len(read_types) > 1 else ''}"
    # Boosting's OpenMP threads wait on each other at every step: nothing else runs beside the trainers. The finished
    # database's build, if still running, is paused (its memory kept) and goes on after them.
    check_jobs()
    paused = final_build is not None and final_build.pause()
    Steps.start(f"training {models} (classifier_training*.log){' in parallel' if len(read_types) > 1 else ''}, "
                f"{trainer_threads} thread{'s' if trainer_threads > 1 else ''} each")
    if paused:
        Steps.done(f"{os.path.basename(db)}'s build in the background paused meanwhile")
    trainers = {}
    for t in read_types:
        command = [sys.executable, TRAINER, "--truth-file", os.path.join(training, TABLES[t]),
                   "--output-prefix", prefixes[t], "--features", args.features, "--model", args.model,
                   "--ntree", str(args.ntree), "--maxnodes", str(max_leaves(args.maxnodes, t)), "--seed", str(args.seed),
                   "--threads", str(trainer_threads), "--scenario-weight", str(args.scenario_weight),
                   "--taxonomy", taxonomy, "--evaluation", args.evaluation]
        if args.rounds:
            command += ["--rounds", str(args.rounds)]
        if args.previous_procedure:
            command += ["--previous-procedure"]
        if t in depth_knob_types(args):
            command += ["--depth-knobs"]
        if args.call_mode == "fdr":
            command += ["--fdr-calls"]
        if has_test and os.path.isfile(os.path.join(test, TABLES[t])):
            command += ["--test-file", os.path.join(test, TABLES[t])]
        trainers[t] = Job(command, step_log("classifier_training" + ("" if t == "pe" else "_" + t) + ".log"),
                          label=f"training the {t} model")
    seconds = max(job.finish().seconds for job in trainers.values())
    if paused:
        final_build.resume()
    each = f" ({', '.join(f'{t} {clock(job.seconds)}' for t, job in trainers.items())})" if len(trainers) > 1 else ""
    Steps.done(f"trained in {clock(seconds)}{each}; {model_scores(read_types, prefixes)}")
    for t, name, why in feature_choices(read_types, prefixes):
        Steps.done(f"{t} model's features: {name}: {why}")
    if selected:
        Steps.done(scenario_scores(read_types, prefixes, selected))
    # The reports of what the models' errors rest on need only the trainers' calls, the samples' SAMs and the training
    # database (neither feeds the database): they start now, in a thread beside the parity check, the composition, the
    # wait for the finished database's build and --add_model (r226 v18 waited 21 min for that build before them).
    trace = heldout if "pe" in read_types and training_db != db and os.path.isfile(heldout) else None

    def report_all():
        try:
            reports(trace, args.error_units, prefixes, training, test if has_test else None, training_db, logs, steps,
                    args.threads, os.path.join(samples_root, CONTIG_CACHE), args.share_logs,
                    heldout if training_db != db else None, args.protal, taxonomy, samples_root)
        except Exception as e:  # noqa: BLE001: the database is ready either way
            say(f"    the reports failed ({e}); the database is ready either way")
    reporting = threading.Thread(target=report_all, name="reports", daemon=True)
    reporting.start()
    # protal must score as the trainer does, and compute the features as it did during collection: the read types'
    # checks side by side (protal --profile_only: no index loaded), their logs then joined in parity.log.
    Steps.start("checking that protal scores the models as the trainer does (parity.log)")
    began = time.time()
    parity_log = step_log("parity.log")
    checks = [Job([sys.executable, PARITY, "--db", training_db, "--model", prefixes[t] + ".xml", "--training", training,
                   "--read_type", t, "--protal", args.protal, "-t", str(max(1, args.threads // len(read_types)))],
                  step_log(f"parity_{t}.log"), label=f"checking the {t} model's parity with protal") for t in read_types]
    for job in checks:
        job.finish()
    concatenate_logs([job.log for job in checks], parity_log, lambda p: os.path.basename(p).removesuffix(".log"))
    for job in checks:
        os.remove(job.log)
    Steps.done(f"{', '.join(read_types)}: the same probabilities and features, checked in {clock(time.time() - began)}")
    # model_logs/parity.txt: each read type's result (check_model_parity.py's parity.txt), which names its model.
    results = []
    for t in read_types:
        path = os.path.join(training, "parity" if t == "pe" else "parity_" + t, "parity.txt")
        if os.path.isfile(path):
            with open(path) as fh:
                results.append(fh.read().rstrip("\n"))
    with open(os.path.join(logs, "parity.txt"), "w") as fh:
        fh.write("\n".join(results) + "\n")
    # The samples' composition against the simulator's truth, with the models' calls: how far the explained share, the
    # unknown share and the genome sizes can be trusted (composition_accuracy.py).
    Steps.start("checking the samples' composition against the truth (model_logs/composition_accuracy*.tsv)")
    began = time.time()
    composition_briefs, composition_text = composition_report(read_types, prefixes, training, test if has_test else None,
                                                              heldout if training_db != db else None, logs, args.threads)
    for t, line in composition_briefs.items():
        Steps.done(f"{t}: {line}")
    Steps.done(f"checked in {clock(time.time() - began)}; details: model_logs/composition_accuracy.txt")
    # The models go into the database, while the reports of what their errors rest on run (report_all, reports()).
    Steps.start(f"adding {models} to {os.path.basename(db)} (final_package.log)")
    if final_build is not None:
        if final_build.seconds is None:
            Steps.done(f"waiting for {os.path.basename(db)}'s build in the background (index_and_package.log)")
        final_build.finish()  # its line comes from built_final_in_background
    elif training_db != db and not final_done:
        Steps.done(f"building {os.path.basename(db)} first (index_and_package.log)")
        job = run(build_command(args.protal, db, args.threads, *final_level), final_log, built_final,
                  f"building {os.path.basename(db)}")
        Steps.done(built(db, job, full_reference_fate(db, scan, "the foreign scan")))
    if scan and training_db != db:
        # The finished database's own scan (every species, from its full reference), before the models go in: --add_model
        # then replaces them in place at the file's end.
        removed = scanned(db)
        if removed:
            Steps.done(f"{os.path.basename(db)}'s full reference, read by the scan{removed}")
    # The trained models replace the shipped one and the placeholders in database.protal; --add_model checks
    # each and replaces them in place at the end of the file (seconds; without placeholders, a read type's
    # model is a new member and the ~20 GB file is rewritten once instead).
    began = time.time()
    run([args.protal, "--add_model", ",".join(prefixes[t] + ".xml" for t in read_types), "--read_type", ",".join(read_types),
         "--db", db, "-t", str(args.threads)], step_log("final_package.log"), label=f"adding {models}")
    Steps.done(f"added in {clock(time.time() - began)}{db_size(db)}")
    versions, changed = versions_at_end(args.versions_at_start, tool_versions(args))
    with open(os.path.join(db, "build_metadata.tsv"), "w") as fh:
        fh.write("".join(f"{k}\t{v}\n" for k, v in provenance(args, versions, release, genome_table, heldout, n_heldout,
                                                                  read_types, prefixes, genes_note, insilico_note)))
    if changed:
        say(f"Warning: {', '.join(k.replace('_', ' ') for k in changed)} changed during the run, so its later steps ran "
            "other scripts or binaries than its first (build_metadata.tsv has both): rebuild the database if that matters")
    # The build's reports beside database.protal (not members of it) go to model_logs/: protal_db/ is the database.
    for name in ("gene_congeners.tsv", "gene_incongruence.tsv"):
        if os.path.isfile(os.path.join(db, name)):
            os.replace(os.path.join(db, name), os.path.join(logs, name))
    summary = summary_lines(read_types, prefixes, db) + ["", *composition_text]
    with open(os.path.join(logs, "summary.txt"), "w") as fh:
        fh.write("\n".join(summary) + "\n")
    with SAYING:  # in one piece, whatever the reports' thread says meanwhile
        print("\n" + "\n".join(summary) + "\n", flush=True)
        console_log("\n" + "\n".join(summary) + "\n")
    say(f"Ready protal database: {db}{db_size(db)}" +
        (f" (marker genes: {genes_note.split(',')[0].split(' (')[0]}, gene_subset.txt)" if subset else "") +
        f"; model evaluation: {logs} (start with trained_model.report.txt, and trained_model_<read type>.report.txt)")
    reporting.join()  # started after the training (report_all)
    with contextlib.suppress(OSError):
        os.remove(os.path.join(training_db, ANCESTRY_REFERENCE))  # read by the ancestry report; 15 GB at r226
    if keep_full:
        if removed := remove_full_reference(training_db):
            say(f"    the training database's full reference, kept for the ancestry report{removed}")
    if Job.scratch:
        Job.scratch.look()
        parts = [("the genome store", genome_store), ("the training database", training_db if training_db != db else None),
                 ("the training samples", training), ("the test samples", test if has_test else None)]
        held = ", ".join(f"{what} {gigabytes(tree_size(path))}" for what, path in parts
                         if path and os.path.isdir(path) and on_disk(path, samples_root))
        say(f"The run took at most {gigabytes(Job.scratch.peak)} on {samples_root} (estimated at most "
            f"~{gigabytes(disk_estimate)} with what it streamed, every read written at once); there now: {held} "
            "(the samples left for a rerun)")
    if genome_store and os.path.isdir(genome_store):
        say(f"The genome store {genome_store} holds {gigabytes(tree_size(genome_store))}, kept for the next build "
            "(remove it to free the space; --genome-store none builds without one)")
    timeline.stop()  # its last row, before the archive packs logs/
    if args.share_logs:
        share_archive(args.outdir, [("training", training), ("test", test if has_test else None)], args.threads,
                      os.path.join(samples_root, "insilico_strains", "insilico_strains.tsv"))
    files = {name: sum(len(names) for _, _, names in os.walk(os.path.join(args.outdir, name)))
             for name in (DATABASE, REPORTS, LOGS, WORK)}
    say(f"In {args.outdir}: {DATABASE}/ the database ({files[DATABASE]} files), {REPORTS}/ the evaluation "
        f"({files[REPORTS]}), {LOGS}/ and console.log every step's log ({files[LOGS]}), {WORK}/ what a rerun reuses "
        f"({files[WORK]}; not needed to use the database)")


if __name__ == "__main__":
    for sig in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP):
        signal.signal(sig, on_signal)
    try:
        main()
    finally:  # an uncaught exception, too
        stop_jobs()
