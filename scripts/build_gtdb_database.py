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
  other strains, too. OUT_DIR/genome_table.txt says how often.
- Species the database lacks. The samples are profiled against a training
  database (OUT_DIR/training_db; SCRATCH/training_db with --scratch) that
  leaves whole clades of every rank (--holdout-clades) and --holdout of the
  other species out
  (OUT_DIR/heldout_species.txt): their reads land on relatives, as those of
  organisms GTDB lacks do in real samples, and the model learns to reject those
  relatives. The report gives false positive and false negative rates by rank.
  The finished database has all species and the model trained so. The training
  database costs a second index build (while the first one runs) and its disk
  space.

One model per read type (--read-types, default pe,se,pb,ont), trained in
parallel: paired-end reads, their first reads alone (single-end), and PacBio
and Nanopore reads of the same communities (all made by simulate_metagenomes). Besides the training
data, an independent test set of another design (--test-*: other depths,
community sizes, abundances and strain mixes) is profiled and scored by each
model: cross-validation on the training data cannot show what its design lacks.
The simulations need no database: both collections simulate in the background,
at a lower priority than the builds, from the moment the species to leave out
are chosen, and profile their samples once the training database is built: as
the simulations go on, in protal runs of --profile-blocks GB of reads, each
design point's reads removed once profiled, so that the samples need not all
be on the disk at once (--profile-blocks 0: both collections in one protal run
once all are simulated).
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

OUT_DIR/model_logs/ collects what tells whether the models are good: summary.txt
(TP, FP, TN, FN, sensitivity, specificity, precision and F1 of each model), each read
type's training report (how it does on species and clades it was not trained on,
on the independent test set, false positive and false negative rates by rank,
against the previous model and training procedure), its numbers as JSON, the
per-taxon predictions, the threshold table, the parity check with protal, the
genome table summary and build_metadata.tsv (what the database was built from
and with); and what the model's conservation features rest on, on real genomes:
gene_congeners.tsv (protal --build: how each gene differs between congeners
against within species), gene_incongruence.tsv (protal --build: every near pair of
gene copies across genera, and which copy is suspect, contamination or a transfer; the
suspect ones go into the database as suspect_copies.tsv and a run leaves their records
out) and relatives_by_gene_conservation.txt (trace_relatives.py: where the reads of the
held-out species land, by the genes' factors). model_logs/error_reads/ keeps the reads behind each model's errors in
every sample of the training data and the test set (--error-reads, default all): per sample, the SAM records of the
reads on its false positives and of its false negatives' reads wherever they went, the non-hits among them (reads that
seeded on taxa but aligned nowhere, whose unmapped records protal writes for these samples), each with its source
genome, and a table of the error taxa (error_reads.py).

A reduced database holds a subset of the marker genes (--n-genes N: the N most
distinctive by prevalence x unique k-mer share, ranked by scripts/rank_genes.py
from a full build of the training database, or from --gene-ranking; --genes: the
genes named). The release is then converted whole into OUTDIR/.converted and both
database folders are derived from it with the subset (their gene neighbours
counted over it); OUTDIR/gene_ranking.tsv and gene_subset.txt record the choice.
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
when their inputs are those of the run that completed them (OUTDIR/.stages), and
the collector reuses the samples it simulated and profiled with the same database,
protal and design.

Each stage writes to a log of its own in OUTDIR. On the console, each line has the
time and how long the run has taken. A step says when it starts ("4/8 training
data (its log): ...") and, on an indented line, when it ends: how long it took, its
peak memory and a few numbers of what it made. With --progress-every, each stage
running also says every so many seconds how long it has run, the memory it takes
and the last line of its log.
"""
import argparse
import collections
import concurrent.futures
import csv
import glob
import gzip
import hashlib
import json
import os
import random
import re
import shutil
import signal
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
CONVERTER = os.path.join(HERE, "mini_db", "gtdb_to_protal_db.py")
GENE_NEIGHBOURS = os.path.join(HERE, "mini_db", "gene_neighbours.py")
TRAINER = os.path.join(HERE, "machine_learning_cmdline.py")
COLLECTOR = os.path.join(HERE, "collect_training_data.py")
PARITY = os.path.join(HERE, "check_model_parity.py")
TRACE = os.path.join(HERE, "trace_relatives.py")
ERROR_READS = os.path.join(HERE, "error_reads.py")
RANKER = os.path.join(HERE, "rank_genes.py")
INSILICO = os.path.join(HERE, "insilico_strains.py")
# --insilico-ani when the conversion left no gene_positions.tsv (--no-gene-neighbours): no real strains to draw from.
INSILICO_FALLBACK_ANI = "97-99.5"
SOURCE = os.path.dirname(HERE)  # the checkout these scripts are part of
# What protal and the simulator are built from (protal_commit.cmake marks a build of uncommitted changes to them).
BUILD_SOURCES = ("src", "lib", "CMakeLists.txt", "protal_config.h.in", "protal_commit.cmake")
ACCESSION = re.compile(r"(?:RS_|GB_)?(GC[AF]_\d{9}\.\d+)")
sys.path.insert(0, os.path.join(HERE, "mini_db"))
sys.path.insert(0, HERE)
import lineages  # noqa: E402
from gtdb_to_protal_db import (clear_build_outputs, full_reference_path, marker_files, normalize_accession,  # noqa: E402
                               read_gene_ids, read_gene_list, read_representatives,
                               remove_full_reference as remove_full_reference_files)
import rank_genes  # noqa: E402
import scenarios  # noqa: E402
from model_pmml import MODEL_FILES, write_placeholder  # noqa: E402
from model_features import DEFAULT_FEATURE_SET, FEATURE_SETS, feature_set_name  # noqa: E402
from collect_training_data import INSILICO_PREFIX, TABLES, clock, congener_spec, last_line, simulation_state, units_of, parse_args as collector_args  # noqa: E402


def congener_text(spec):
    """A --congeners value (congener_spec) as the collector's option."""
    return spec[1] if spec[0] == "groups" else str(spec[1])

STARTED = time.time()


def say(message):
    """Prints a message, its first line headed by the time and how long the run has taken."""
    print(f"[{time.strftime('%H:%M:%S')} +{clock(time.time() - STARTED)}] {message}", flush=True)


class Steps:
    """The steps of the run on the console: "3/8 what it does (its log)" when one starts, and indented lines
    of how it went (done())."""
    total, current = 0, 0

    @classmethod
    def start(cls, text):
        cls.current += 1
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

    def __init__(self, command, log, on_success=None, label=None, nice=0):
        """nice: the command's niceness, more than the script's (the simulations, beside the builds)."""
        os.makedirs(os.path.dirname(log), exist_ok=True)
        self.command, self.log, self.on_success, self.started = command, log, on_success, time.time()
        self.label = label or os.path.basename(log).removesuffix(".log")
        self.seconds = None  # set when it has ended
        self.peak = None  # the most memory it, or a command it ran, took, in bytes; set when it has ended
        self.fh = open(log, "w")
        self.process = subprocess.Popen(command, stdout=self.fh, stderr=subprocess.STDOUT, start_new_session=True,
                                        preexec_fn=(lambda: os.nice(nice)) if nice else None)
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
        return self.process.returncode

    def ended(self, rc):
        self.fh.close()
        Job.running.remove(self)
        self.seconds = time.time() - self.started
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
            if now - checked >= 5:
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

    def status(self):
        """How the command is doing: how long it has run, the memory it and the commands it started take, and
        the last line of its log."""
        memory = group_memory(self.process.pid)
        line = last_line(self.log, 200)
        return (f"{self.label}: {clock(time.time() - self.started)} so far" +
                (f", {gigabytes(memory)} in memory" if memory else "") +
                (f"; {os.path.basename(self.log)}: {line}" if line else ""))

    def took(self):
        """How long it took, and its peak memory."""
        return clock(self.seconds) + (f", peak memory {gigabytes(self.peak)}" if self.peak else "")

    def kill(self):
        """Stops the command and every process of its group."""
        for sig, wait in ((signal.SIGTERM, 60), (signal.SIGKILL, None)):
            try:
                os.killpg(self.process.pid, sig)
            except (ProcessLookupError, PermissionError):
                break
            try:
                self.process.wait(timeout=wait)
                break
            except subprocess.TimeoutExpired:
                continue
        self.fh.close()


def check_jobs():
    """Ends the jobs whose command has ended: one that failed stops the script."""
    for job in list(Job.running):
        rc = job.poll()
        if rc is not None:
            job.ended(rc)


def run(command, log, on_success=None, label=None):
    return Job(command, log, on_success, label).finish()


def make_genome_table(gtdb, release, output, extra_dirs=(), species=None, threads=1):
    """The simulator's genome table: every genome FASTA of a GTDB species found in the release (and in
    extra_dirs), or only those of the species in the set `species`. A fourth column holds each genome's length
    (genome_length): without it simulate_metagenomes reads every genome of the table for its length, once per
    run, that is once per design point (~20 ms per genome; 20,000 genomes at r226). Lengths of an earlier table
    at `output` are kept for the FASTAs not changed since it was written."""
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
    paths = {}
    for root in genome_dirs:
        if not os.path.isdir(root):
            continue
        for pattern in ("**/*.fna", "**/*.fna.gz", "**/*.fa", "**/*.fa.gz", "**/*.fasta", "**/*.fasta.gz"):
            for path in glob.iglob(os.path.join(root, pattern), recursive=True):
                m = ACCESSION.search(os.path.basename(path))
                if m and m.group(1) in taxonomy:
                    paths.setdefault(m.group(1), path)
    if not paths:
        sys.exit("No extracted whole genome FASTAs found under genomic_files_all/gtdb_genomes_all_r" + release +
                 " or genomic_files_reps/gtdb_genomes_reps_r" + release + "; extract GTDB genome files or pass --genome-table")
    rows = [[a, taxonomy[a], os.path.abspath(paths[a])] for a in sorted(paths)]
    write_genome_table(rows, output, threads)
    return len(paths)


def write_genome_table(rows, output, threads):
    """Writes a genome table, rows [accession, taxonomy, FASTA path, ...], with each genome's length
    (genome_length) as the fourth column; lengths of an earlier table at `output` are kept for the FASTAs not
    changed since it was written."""
    known = {}  # path -> length, from an earlier table, for FASTAs older than it
    if os.path.isfile(output):
        written = os.path.getmtime(output)
        with open(output) as fh:
            for fields in (line.rstrip("\n").split("\t") for line in fh):
                if len(fields) >= 4 and fields[3].isdigit() and os.path.isfile(fields[2]) and \
                        os.path.getmtime(fields[2]) < written:
                    known[fields[2]] = int(fields[3])
    todo = sorted({r[2] for r in rows if r[2] not in known})
    if todo:  # in processes: counting the letters holds the GIL
        with concurrent.futures.ProcessPoolExecutor(max(1, min(threads, len(todo)))) as pool:
            known.update(zip(todo, pool.map(genome_length, todo, chunksize=8)))
    with open(output + ".partial", "w") as fh:
        fh.writelines("\t".join(r[:3] + [str(known[r[2]])]) + "\n" for r in rows)
    os.replace(output + ".partial", output)


def with_lengths(table, output, threads):
    """The genome table to simulate from: `table` itself if every row has a genome length, else a copy with them
    (write_genome_table) at `output`. Rows that are no genome (a header, comments) are left out of the copy;
    the simulator reads a table of four columns without a header."""
    with open(table) as fh:
        rows = [line.rstrip("\n").split("\t") for line in fh if line.strip() and not line.startswith("#")]
    genomes = [r for r in rows if len(r) >= 3 and os.path.isfile(r[2])]
    if len(genomes) == len(rows) and all(len(r) >= 4 and r[3].isdigit() for r in rows):
        return table
    write_genome_table(genomes, output, threads)
    return output


NON_LETTERS = bytes(b for b in range(256) if not (65 <= b <= 90 or 97 <= b <= 122))


def genome_length(path):
    """A genome's length as simulate_metagenomes counts it (read_genome_length): the letters of the lines that
    do not start with '>', of the gzip-decompressed file if its name ends in .gz."""
    with open(path, "rb") as fh:
        data = fh.read()
    if path.endswith(".gz"):
        data = gzip.decompress(data)
    letters = len(data.translate(None, NON_LETTERS))
    return letters - sum(len(h.translate(None, NON_LETTERS)) for h in re.findall(rb"^>[^\n]*", data, re.M))


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


def choose_holdout(genome_table, taxonomy, fraction, clades, max_share, seed):
    """The species a training database leaves out, as {species: (rank, clade)}; rank "species" for single
    species. First whole clades, clades[rank] of each rank from phylum down: drawn among those with at least
    two species the genome table can simulate (so that samples can have them), with at most max_share of the
    database's species, and in no clade drawn before. Then a random `fraction` of the species the genome
    table can simulate that no clade took, the same fraction in each domain."""
    lineage_by_id, _ = lineages.from_taxonomy(taxonomy)
    db = {lin["species"]: lin for lin in lineage_by_id.values() if "species" in lin}
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
        for s in rng.sample(species, round(fraction * len(species))):
            chosen[s] = ("species", s)
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


def describe_holdout(chosen, pool):
    """Lines saying what the training database leaves out."""
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


def trace_relatives(training, training_db, heldout, logs, outdir, threads=1, contig_cache=None):
    """model_logs/relatives_by_gene_conservation.txt (trace_relatives.py): where the paired-end reads of the species the
    training database lacks land, by the genes' conservation factors, on real genomes. A failure is reported, and does
    not stop the build: the models do not depend on it."""
    log = os.path.join(outdir, "trace_relatives.log")
    out = os.path.join(logs, "relatives_by_gene_conservation")
    began = time.time()
    with open(log, "w") as fh:
        rc = subprocess.run([sys.executable, TRACE, "--points", os.path.join(training, "points"), "--db", training_db,
                             "--heldout", heldout, "--out", out, "--threads", str(threads)] +
                            (["--contig-cache", contig_cache] if contig_cache else []),
                            stdout=fh, stderr=subprocess.STDOUT).returncode
    if rc:
        say(f"    tracing the held-out species' reads failed ({rc}; see {log}); the build goes on")
        return
    with open(out + ".txt") as fh:
        text = fh.read()
    first = text.splitlines()[0] if text.startswith("Not traced") else "model_logs/relatives_by_gene_conservation.txt"
    say(f"    the held-out species' reads by gene conservation, in {clock(time.time() - began)}: {first}")


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


def error_reads(units, prefixes, training, test, training_db, logs, outdir, threads=1, contig_cache=None):
    """model_logs/error_reads/<read type>/ (error_reads.py): the SAM records of the reads behind each model's false
    positives and false negatives in the samples of --error-reads, which protal wrote with an unmapped record for every
    read that seeded on taxa but aligned nowhere (collect_training_data.py --unmapped_reads), and a table of their error
    taxa. A failure is reported, and does not stop the build: the models do not depend on it."""
    log = os.path.join(outdir, "error_reads.log")
    began, told = time.time(), []
    scopes = collections.defaultdict(list)
    for kind, scope in units:
        scopes[kind].append(scope)
    with open(log, "w") as fh:
        for kind, which in scopes.items():
            calls, out = prefixes[kind] + ".calls.tsv.gz", os.path.join(logs, "error_reads", kind)
            if not os.path.isfile(calls):
                told.append(f"{kind}: the trainer wrote no calls")
                continue
            command = [sys.executable, ERROR_READS, "--calls", calls, "--training", training, "--db", training_db,
                       "--samples", "all" if None in which else ",".join(which), "--read-type", kind, "--out", out,
                       "--threads", str(threads)]
            command += ["--test", test] if test and os.path.isdir(test) else []
            command += ["--contig-cache", contig_cache] if contig_cache else []
            fh.write(" ".join(command) + "\n")
            fh.flush()
            rc = subprocess.run(command, stdout=fh, stderr=subprocess.STDOUT).returncode
            if rc:
                say(f"    taking the reads of the {kind} model's errors failed ({rc}; see {log}); the build goes on")
                continue
            told.append(f"{kind} {error_reads_summary(os.path.join(out, 'summary.tsv'))}")
    if told:
        say(f"    the reads of the models' errors (model_logs/error_reads, in {clock(time.time() - began)}): "
            + "; ".join(told))


def error_reads_summary(path):
    """The samples, error taxa, fragments and size of error_reads.py's summary.tsv, in words."""
    if not os.path.isfile(path):
        return "no samples"
    c = collections.Counter()
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            c.update({"samples": 1, row["set"]: 1, "FP": int(row["FP"]), "FN": int(row["FN"]),
                      "unseen": int(row["unseen"]), "fragments": int(row["fragments"]), "bytes": int(row["sam_bytes"])})
    if not c["samples"]:
        return "no samples"
    return (f"{c['samples']} samples ({c['training']} training, {c['test']} test): {c['FP']} FP, {c['FN']} FN and "
            f"{c['unseen']} unseen taxa, {c['fragments']} fragments, {gigabytes(c['bytes'])}")


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
            ("marker_genes", genes),
            ("gene_conservation", gene_conservation_summary(os.path.join(args.outdir, "index_and_package.log"))),
            ("gene_neighbours", gene_neighbours_summary(os.path.join(args.outdir, "index_and_package.log"))),
            ("gene_positions", gene_neighbours_summary(os.path.join(args.outdir, "index_and_package.log"),
                                                       "Gene positions:")),
            ("gene_congeners", gene_congeners_summary(os.path.join(args.outdir, "index_and_package.log"))),
            ("suspect_copies", suspect_copies_summary(os.path.join(args.outdir, "index_and_package.log"))),
            ("classifier_features", args.features), ("classifier_model", args.model),
            ("classifier_trees", args.ntree if args.model == "forest" else f"{args.rounds or TRAINER_ROUNDS} rounds"),
            ("classifier_max_leaves", ",".join(f"{t}:{max_leaves(args.maxnodes, t)}" for t in read_types)),
            ("classifier_evaluation", args.evaluation),
            ("classifier_previous_procedure",
             "compared" if args.previous_procedure and args.evaluation != "none" else "not compared"),
            ("classifier_training_species_left_out", n_heldout)]
    rows += [(f"classifier_training_{rank}_clades_left_out", clade_counts[rank]) for rank in CLADE_RANKS
             if clade_counts.get(rank)]
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


def remove_full_reference(folder):
    """Removes a folder's full reference (full_reference.fna.zst) once its build is done: the marker genes of
    every genome (86 GB raw at r226), which only that build reads. A rebuild converts the release again (the
    build packed reference.fna into database.protal). Says what it removed ("; ..."), or nothing."""
    path = full_reference_path(folder)
    if not path:
        return ""
    return f"; {os.path.basename(path)} removed ({gigabytes(remove_full_reference_files(folder))})"


def build_command(protal, db, threads, *extra):
    command = [protal, "--build", "--no_profile", "-t", str(threads), "--db", db,
               "--reference", os.path.join(db, "reference.fna"), *extra]
    if full_reference_path(db):
        command += ["--full_reference", full_reference_path(db)]
    return command


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
        try:
            with open(os.path.join(self.folder, name + ".json")) as fh:
                return json.load(fh) == json.loads(json.dumps(key))
        except (OSError, ValueError):
            return False

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
    source, would otherwise stop it after the conversion, the builds or the collection."""
    problems, notes = [], []

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
    imports = subprocess.run([sys.executable, "-c", "import joblib, numpy, pandas, sklearn"], capture_output=True, text=True)
    if imports.returncode:
        problems.append(f"{sys.executable} cannot import what the trainer needs ({imports.stderr.strip().splitlines()[-1]}): "
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
    p.add_argument("--outdir", required=True, help="output root; database is written to OUTDIR/protal_db")
    p.add_argument("--release", help="GTDB release number; detected from taxonomy filenames by default")
    p.add_argument("--genome-table", help="optional simulator table: accession, taxonomy, whole genome FASTA path "
                                          "(and genome length; without, the run writes a copy with the lengths, "
                                          "OUTDIR/genomes.tsv)")
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
    p.add_argument("--one-build-at-a-time", action="store_true",
                   help="build the finished database after the model is trained, not while the training data are "
                        "collected (the default, which needs the memory of two builds, or of one build and the "
                        "collection's protal runs, at once)")
    p.add_argument("--training-db-level", type=int, default=3,
                   help="zstd level of the training database (default 3)")
    p.add_argument("--final-db-level", type=int, default=9,
                   help="zstd level of the finished database (default 9: at GTDB r226, protal's default 19 made the index "
                        "2.7%% smaller than level 3 for 21 more minutes of a 1:20 build; docs/claude/2026-10-02-r226-build-"
                        "evaluation)")
    p.add_argument("--n-genes", type=int, metavar="N",
                   help="build the database from the N most distinctive of the release's marker genes (a reduced "
                        "database: less memory, fewer hits per genome): the genes ranked by prevalence x unique "
                        "k-mer share (scripts/rank_genes.py) from a full build of the training database first "
                        "(OUTDIR/gene_ranking.tsv; --gene-ranking skips that build), the N best in "
                        "OUTDIR/gene_subset.txt; the training database and the finished database hold those genes "
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
                        "(unpacked on the scratch disk, scripts/rank_genes.py) into OUTDIR/gene_ranking.tsv, for a "
                        "reduced database of the same release (--n-genes N --gene-ranking OUTDIR/gene_ranking.tsv, "
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
    p.add_argument("--strains-per-species", default="0.3,0.1",
                   help="probabilities of a second, third, ... strain of a species in a sample (default 0.3,0.1): "
                        "real samples often mix strains, which changes the allele-frequency features")
    p.add_argument("--abundance", default="lognormal:1.3,2.0",
                   help="abundance model of the training samples: lognormal:SIGMA, powerlaw:ALPHA or negbin:R:P; "
                        "lognormal:S1,S2,... gives a design point's samples the sigmas in turn (default "
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
                        "evaluation). They go to OUTDIR/insilico_strains, the table simulated from to "
                        "OUTDIR/genomes_simulated.tsv")
    p.add_argument("--insilico-ani", metavar="MIN-MAX",
                   help="draw the in-silico strains' genome ANI uniformly from MIN-MAX (e.g. 95-99) instead of the "
                        "real strains' marker divergence")
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
                   help="scenarios to train on and score the models on, besides the design (default: all four; none "
                        "for none; scenarios.py): gut (~400 "
                        "species, 5%% lacking from the database, Illumina PE 150 Q35 at 20M pairs, PacBio and Nanopore "
                        "at the same bases), soil (~10,000 species, 60%% lacking, Ultima SE 300 Q25 at 20M reads, "
                        "Illumina, PacBio and Nanopore at the same bases), soil_shallow (the soil communities at 5M "
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
                        "depth is its scenario's times a factor from 1/2 to 2 (depth_spread)")
    p.add_argument("--scenario-file", help="JSON of scenarios by name, which add to or change the presets (scenarios.py)")
    p.add_argument("--scenario-samples", type=int, default=6,
                   help="hold-in samples per scenario, in the training data (default 6; 3 before 2026-10-06; 0: the "
                        "scenarios are scored, not trained on)")
    p.add_argument("--scenario-test-samples", type=int, default=3,
                   help="hold-out samples per scenario, in the test set (default 3; 2 before 2026-10-06)")
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
                        "the same), and once the models are trained model_logs/error_reads/ keeps, per sample, the SAM "
                        "records of the reads behind the model's false positives and false negatives, with their source "
                        "genomes (error_reads.py). all (default: every sample of the training data and the test set, "
                        "every read type), none, or READ_TYPE, READ_TYPE:design (the design's samples) or "
                        "READ_TYPE:SCENARIO, comma-separated; read types and scenarios not collected are left out")
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
                        "only the tables to OUTDIR; the training database is built there too (SCRATCH/training_db, "
                        "~22 GB at r226). A network file system (OUTDIR's, often) is slow at the many files the "
                        "simulators write and delete, and at writing a database; the converter spools the release's "
                        "marker genes there too. With the defaults the r226 run took up to 120 GB there, without the training database "
                        "(give it 175 GB, docs/databases.md); a rerun reuses the samples in the same SCRATCH")
    p.add_argument("--profile-blocks", type=float, default=20.0,
                   help="once the training database is built, profile the simulated samples as their simulations go "
                        "on, in protal runs of at least this many GB of reads (collect_training_data.py --follow; the "
                        "two collections' runs take turns), and remove each design point's reads once all its read types "
                        "are profiled (the SAMs, profiles and dumps are kept): the samples need not all be on the disk "
                        "at once. 0: profile both collections in one protal run once every sample is simulated, and "
                        "keep the reads (default 20)")
    p.add_argument("--keep-free", type=float, default=30.0,
                   help="with --profile-blocks: GB a simulation leaves free on the disk of the samples (--scratch or "
                        "OUTDIR), or it waits until profiled reads are removed (default 30)")
    p.add_argument("--stream-above", type=float, default=2.0,
                   help="with --profile-blocks: a design point whose largest sample's reads would take more than this "
                        "many GB (compressed, estimated) is not written to the disk: protal reads its samples from named "
                        "pipes as simulate_metagenomes makes them, in a protal run of its own, a sample at a time "
                        "(collect_training_data.py --stream_above); smaller ones are written and profiled in blocks "
                        "(default 2; 0: none)")
    p.add_argument("--read-compression", choices=["zstd", "gzip"], default="zstd",
                   help="how the simulated samples' reads are written (collect_training_data.py --read_compression): "
                        "zstd (.fq.zst, the default: as small as BGZF or smaller, several times faster to write) or "
                        "gzip (.fq.gz); protal reads both")
    args = p.parse_args()
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
    scenarios_given = args.scenarios is not None  # by default all four, the host scenario only with a host genome
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
    check_tools(args, read_types)
    args.versions_at_start = tool_versions(args)  # build_metadata.tsv: what the run starts with (versions_at_end)
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
    db = os.path.join(args.outdir, "protal_db")
    os.makedirs(db, exist_ok=True)
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
    samples_root = args.outdir  # where the collections simulate and profile their samples
    if args.scratch:
        samples_root = os.path.abspath(args.scratch)
        os.makedirs(samples_root, exist_ok=True)
        Job.scratch = Scratch(samples_root)
    say(f"A protal database of GTDB r{release} in {args.outdir}, with {args.threads} threads; each step logs to a "
        "file there" + (f"; the simulated samples go to {samples_root} "
                            f"({gigabytes(shutil.disk_usage(samples_root).free)} free)" if args.scratch else "") +
        (f"; the stages running are reported every {clock(Job.progress_every)}" if Job.progress_every > 0 else ""))
    # The steps: genome table, release, marker genes (with --n-genes or --genes), databases, training data, test
    # set (if any), models, parity, packing.
    subset = args.n_genes is not None or bool(args.genes)
    has_test = args.test_samples > 0 or bool(hold_out)  # a test collection: the design's test set, the scenarios' hold-out
    Steps.total = 7 + has_test + subset + (args.insilico_strains > 0)
    if args.genes:
        # The list is checked against the release's marker files before anything is converted (marker ids by
        # name, gene ids by their range; the ids themselves come from gene2geneid.tsv once it is there).
        markers = sorted({m for m, _ in marker_files(args.gtdb, "genomic_files_reps", "reps", release)})
        if markers:
            read_gene_list(args.genes, {m: i + 1 for i, m in enumerate(markers)})
    genome_table = args.genome_table or os.path.join(args.outdir, "genomes.tsv")
    if not args.genome_table:
        pool = None
        if args.simulate_species:
            with open(args.simulate_species) as fh:
                pool = {n if n.startswith("s__") else "s__" + n
                        for n in (line.rstrip("\n").split("\t")[0].strip() for line in fh) if n and not n.startswith("#")}
        make_genome_table(args.gtdb, release, genome_table, args.extra_genomes, pool, args.threads)
    elif args.extra_genomes or args.simulate_species:
        sys.exit("--extra-genomes and --simulate-species shape the genome table this script makes; "
                 "apply them to --genome-table instead")
    else:  # with the genomes' lengths, if it has none: the simulator would read every genome for them
        genome_table = with_lengths(args.genome_table, os.path.join(args.outdir, "genomes.tsv"), args.threads)
    logs = os.path.join(args.outdir, "model_logs")
    os.makedirs(logs, exist_ok=True)
    reps = read_representatives(args.gtdb, release)
    # The share of simulated species that are strains is judged after the in-silico strains (step 3), if any.
    summary, brief, warning = summarize_genome_table(genome_table, reps, insilico_to_come=args.insilico_strains > 0)
    with open(os.path.join(args.outdir, "genome_table.txt"), "w") as fh:
        fh.write("\n".join(summary) + "\n")
    Steps.start(f"genome table ({os.path.basename(genome_table)}, genome_table.txt): {brief}")
    if warning:
        say(warning)

    # A rerun skips what an earlier run into OUTDIR completed with the same inputs (see Stages).
    stages = Stages(os.path.join(args.outdir, ".stages"))
    convert_key = {"converter": content_hash(CONVERTER), "release": release, "gtdb": release_identity(args.gtdb, release),
                   "placeholders": not args.no_placeholder_models,
                   "gene_neighbours": None if args.no_gene_neighbours else [content_hash(GENE_NEIGHBOURS), content_hash(genome_table)]}
    final_key = {"convert": convert_key, "protal": file_identity(args.protal), "level": args.final_db_level}
    final_done = stages.done("protal_db", final_key) and os.path.isfile(os.path.join(db, "database.protal"))
    # --build packs the taxonomy into database.protal; the collector and the trainer read it (domains,
    # representative genomes). The training database has the same taxonomy.
    taxonomy = os.path.join(args.outdir, "internal_taxonomy.dmp")
    gene_table = os.path.join(args.outdir, "gene2geneid.tsv")  # the markers' gene ids, for --genes and the ranking

    def convert(into):
        """The converted files of the release in `into`, with the models of the read types but pe as
        placeholders (protal warns when it loads one), packed by --build like model_pe.xml. Says how long it
        took."""
        # With --scratch, the converter spools the marker genes on the node's disk (~5 GB compressed at r226).
        job = run([sys.executable, CONVERTER, "--gtdb", args.gtdb, "--outdir", into, "--release", release, "-t",
                   str(args.threads)] + (["--tmp", samples_root] if args.scratch else []),
                  os.path.join(args.outdir, "convert.log"), label="converting the release")
        shutil.copyfile(os.path.join(into, "internal_taxonomy.dmp"), taxonomy)
        shutil.copyfile(os.path.join(into, "gene2geneid.tsv"), gene_table)
        took = f"converted in {job.took()}"
        if not args.no_gene_neighbours:
            # Which marker genes lie next to which in every genome to simulate from (their genes placed by
            # their sequence or k-mer trace): gene_neighbours.tsv, the frequencies per clade, and
            # gene_positions.tsv, which --build both pack; the training database's copy derives the frequencies
            # anew from the positions, without the species it leaves out.
            log = os.path.join(args.outdir, "gene_neighbours.log")
            job = run([sys.executable, GENE_NEIGHBOURS, "--db", into, "--genome_table", genome_table, "-t", str(args.threads)],
                      log, label="finding the genes' neighbours")
            took += f"; the genes placed{genes_placed(log)} and their neighbours counted in {job.took()} (gene_neighbours.log)"
        if not args.no_placeholder_models:
            for read_type in ("se", "pb", "ont"):
                write_placeholder(os.path.join(into, MODEL_FILES[read_type]), read_type)
        return took

    converted = None  # the folder with the converted files, once there
    # With a gene subset the whole release is converted here, and the database folders (the subset's genes,
    # with or without the species held out) are derived from it; without one, the finished database's folder
    # holds the conversion and its build consumes it (the folder serves the taxonomy alone otherwise).
    full = os.path.join(args.outdir, ".converted")

    def ensure_converted():
        """With a gene subset: the whole release in `full`, kept from an earlier run that stopped before the
        database folders were derived, or converted now."""
        nonlocal converted
        if converted:
            return
        if stages.done("convert", convert_key) and os.path.isfile(taxonomy) and os.path.isfile(gene_table) and \
                all(os.path.isfile(os.path.join(full, f)) for f in ("reference.fna", "reference.map", "internal_taxonomy.dmp")):
            clear_build_outputs(full)
            Steps.done(f"{os.path.basename(full)} holds the release converted by an earlier run; not converted again")
        else:
            Steps.done("the whole release: " + convert(full))
            stages.mark("convert", convert_key)
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
    elif stages.done("convert", convert_key) and os.path.isfile(taxonomy) and \
            all(os.path.isfile(os.path.join(db, f)) for f in ("reference.fna", "reference.map", "internal_taxonomy.dmp")):
        # Converted by an earlier run that stopped before the build packed the files.
        clear_build_outputs(db)
        converted = db
        Steps.start(f"the release: {db} holds the release converted by an earlier run; not converted again")
    else:
        Steps.start(f"converting GTDB r{release} (convert.log)")
        stages.forget("protal_db")
        Steps.done(convert(db))
        stages.mark("convert", convert_key)
        converted = db
    if not subset and not os.path.isfile(taxonomy):
        Steps.done("its taxonomy: " + convert(full))
        converted = full

    # In-silico strains (insilico_strains.py): every species of the table with one genome gets a mutated copy of
    # its representative (codon-aware substitutions, the divergence of the table's real strains, the genes'
    # conservation factors), so that it is simulated from a strain as often as a species with two genomes. Without
    # them a model given GTDB's cluster sizes learns that a divergent read cloud on a one-genome species is a
    # relative the database lacks (docs/claude/2026-10-04-r226-v10-evaluation). The simulations draw from
    # genomes_simulated.tsv; the species to leave out and the gene neighbours come from the table itself.
    sim_table, insilico_note = genome_table, "none (--insilico-strains 0)"
    if args.insilico_strains > 0:
        sim_table = os.path.join(args.outdir, "genomes_simulated.tsv")
        log = os.path.join(args.outdir, "insilico_strains.log")
        Steps.start("in-silico strains of the species with one genome (insilico_strains.log, genomes_simulated.tsv)")
        insilico_key = {"script": content_hash(INSILICO), "genome_table": content_hash(genome_table),
                        "convert": convert_key, "share": args.insilico_strains, "ani": args.insilico_ani,
                        "seed": args.seed}
        if stages.done("insilico", insilico_key) and os.path.isfile(sim_table):
            Steps.done("made by an earlier run from the same table; kept")
        else:
            stages.forget("insilico")

            def positions_file():
                return next((p for p in (os.path.join(f, "gene_positions.tsv") for f in (converted, db, full) if f)
                             if os.path.isfile(p)), None)
            if positions_file() is None and subset and not args.no_gene_neighbours:
                ensure_converted()
            positions = positions_file()
            command = [sys.executable, INSILICO, "--genome-table", genome_table, "--output", sim_table, "--out-dir",
                       os.path.join(args.outdir, "insilico_strains"), "--share", str(args.insilico_strains),
                       "--seed", str(args.seed), "-t", str(args.threads)]
            if positions:
                command += ["--positions", positions, "--taxonomy", taxonomy]
            if args.insilico_ani or not positions:
                command += ["--ani", args.insilico_ani or INSILICO_FALLBACK_ANI]
            job = run(command, log, label="making in-silico strains")
            stages.mark("insilico", insilico_key)
            Steps.done(f"made in {job.took()}" + ("" if positions else
                                                    f" (no gene_positions.tsv: ANI {args.insilico_ani or INSILICO_FALLBACK_ANI})"))
        insilico_note = last_line(log)
        Steps.done(insilico_note)
        summary, brief, warning = summarize_genome_table(sim_table, reps)
        with open(os.path.join(args.outdir, "genome_table.txt"), "a") as fh:
            fh.write("with the in-silico strains (genomes_simulated.tsv, the genomes simulated from):\n" +
                     "\n".join(summary[1:]) + "\n" + insilico_note + "\n")
        Steps.done(f"simulated from: {brief}")
        if warning:
            say(warning)

    # The training database leaves some species out: the model then sees reads of species the database
    # lacks, which land on relatives, and reads of whole families, classes and phyla it lacks, which land on
    # distant ones. It is made from the converted files before --build packs them. heldout_species.txt: the
    # species, the rank they were held out at and the clade.
    heldout = os.path.join(args.outdir, "heldout_species.txt")
    clades = parse_clades(args.holdout_clades)
    if args.holdout_species:
        shutil.copyfile(args.holdout_species, heldout)
    elif args.holdout > 0 or clades:
        chosen = choose_holdout(genome_table, taxonomy, args.holdout, clades, args.holdout_max_share, args.seed)
        with open(heldout, "w") as fh:
            fh.write("".join(f"{s}\t{rank}\t{clade}\n" for s, (rank, clade) in sorted(chosen.items())))
    elif os.path.exists(heldout):
        os.remove(heldout)
    for note in scenario_notes:
        say(f"    WARNING: scenario {note}")
    if selected:
        # A scenario's share of species the database lacks comes from the species held out: enough of them, and of
        # the others, for its largest sample (scenarios.table_split), or the scenario is scaled down to what they hold
        # (scenarios.fit_species), as its collections do it.
        novel = read_holdout(heldout) if os.path.exists(heldout) else {}
        pool = pool_species(sim_table)
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
    subset_file = os.path.join(args.outdir, "gene_subset.txt")
    ranking_file = os.path.join(args.outdir, "gene_ranking.tsv")
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
                shutil.copyfile(args.gene_ranking, ranking_file)
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
                    run(command, os.path.join(args.outdir, "gene_ranking_files.log"), label="writing the ranking database's files")
                    job = run(build_command(args.protal, ranking_db, args.threads, "--compress_level", str(args.training_db_level),
                                            "--no_bundle"),
                              os.path.join(args.outdir, "gene_ranking_build.log"), label="building the ranking database")
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
        shutil.copy(subset_file, logs)
        if os.path.isfile(ranking_file):
            shutil.copy(ranking_file, logs)
        final_key["genes"] = content_hash(subset_file)
        final_done = stages.done("protal_db", final_key) and os.path.isfile(os.path.join(db, "database.protal"))
        if final_done:
            Steps.done(f"{os.path.basename(db)} was built by an earlier run from these genes; kept" + remove_full_reference(db))
        else:
            ensure_converted()
            stages.forget("protal_db")
            job = run([sys.executable, CONVERTER, "--from_db", converted, "--genes", subset_file, "--outdir", db,
                       "-t", str(args.threads)], os.path.join(args.outdir, "protal_db_files.log"),
                      label=f"deriving {os.path.basename(db)}'s files")
            Steps.done(f"{os.path.basename(db)}'s files derived for these genes in {job.took()} (protal_db_files.log)")

    training_db, training_done = db, False
    n_heldout, files_took = 0, ""
    if os.path.exists(heldout):
        chosen = read_holdout(heldout)
        n_heldout = len(chosen)
        with open(os.path.join(logs, "holdout.txt"), "w") as fh:
            fh.write("\n".join([f"training database: {n_heldout} species left out ({heldout})"] +
                               describe_holdout(chosen, pool_species(genome_table))) + "\n")
        # Read only by the collections and the parity check: on --scratch, its build writes and they load it from
        # the node's disk (writing database.protal to a network file system was 6 of the 15 min of an r226 build).
        training_db = os.path.join(samples_root, "training_db")
        Steps.start(f"training database ({os.path.basename(training_db)}, training_db_index.log): {n_heldout} species left "
                    f"out, {holdout_brief(chosen)} (model_logs/holdout.txt)")
        training_key = {"convert": convert_key, "heldout": content_hash(heldout), "protal": final_key["protal"],
                        "level": args.training_db_level}
        if subset:
            training_key["genes"] = final_key["genes"]
        training_done = stages.done("training_db", training_key) and \
            os.path.isfile(os.path.join(training_db, "database.protal"))
        if training_done:
            Steps.done(f"{training_db} was built by an earlier run with the same species left out"
                       f"{' and the same genes' if subset else ''}; kept" + remove_full_reference(training_db))
        else:
            stages.forget("training_db")
            if subset:
                ensure_converted()
            elif converted is None:  # the finished database's build consumed them
                converted = full
                Steps.done("the release again (the finished database's build took its files): " + convert(converted))
            job = run([sys.executable, CONVERTER, "--from_db", converted, "--exclude_species", heldout, "--outdir",
                       training_db, "-t", str(args.threads)] + (["--genes", subset_file] if subset else []),
                      os.path.join(args.outdir, "training_db.log"),
                      label="leaving the species out")
            files_took = f"; its files written in {job.took()} (training_db.log)"
    if converted and converted != db:
        shutil.rmtree(converted, ignore_errors=True)

    # The finished database is needed only for --add_model at the end: it is built in the background from the
    # start, while the training database is built and the training data collected, unless one build at a time.
    final_log = os.path.join(args.outdir, "index_and_package.log")
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
        Steps.done(built(db, final_build, remove_full_reference(db)))

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
        if scenario_samples:
            command += ["--scenarios", ",".join(f"{name}:{n}" for name, n in scenario_samples.items())]
            command += ["--scenario_file", os.path.abspath(args.scenario_file)] if args.scenario_file else []
            command += ["--host_genome", args.host_genome] if args.host_genome else []
        if args.error_units:
            command += ["--unmapped_reads", ",".join(kind + (f":{scope}" if scope else "") for kind, scope in args.error_units)]
        return command

    training = os.path.join(samples_root, "training")
    test = os.path.join(samples_root, "test")
    collections_ = [("training data", collect_command(training, args.samples, args.read_pairs, args.species_per_sample,
                                                      args.abundance, args.strains_per_species, args.long_read_bases,
                                                      args.long_read_samples, args.seed, hold_in),
                     os.path.join(args.outdir, "training_data.log"))]
    if has_test:  # the scenarios' hold-out samples are of another seed, as the test set's design is
        collections_.append(("independent test set" if args.test_samples > 0 else "scenarios' hold-out samples",
                             collect_command(test, args.test_samples, args.test_read_pairs, args.test_species_per_sample,
                                             args.test_abundance, args.test_strains_per_species,
                                             args.test_long_read_bases, args.test_long_read_samples, args.seed + 1000,
                                             hold_out),
                             os.path.join(args.outdir, "test_data.log")))

    if training_db == db:
        Steps.start(f"database ({os.path.basename(db)}, index_and_package.log): no species left out" +
                    ("; built by an earlier run, kept" if final_done else ""))
    # The simulations need neither database, only the genome table and the species left out: both collections
    # simulate from here on, in the background and at a lower priority than the builds (collect_training_data.py
    # --simulate_only), and profile what they simulated once the training database is there (collect(), or with
    # --profile-blocks follow(): as they simulate, the reads profiled removed, so the simulations keep --keep-free GB
    # free and wait for that while the disk is full).
    simulations = {}
    for what, command, log in collections_:
        keep_free = ["--min_free", f"{args.keep_free:g}"] if args.profile_blocks > 0 and args.keep_free > 0 else []
        if args.profile_blocks > 0 and args.stream_above > 0:  # the follower streams them (the same value for both)
            keep_free += ["--stream_above", f"{args.stream_above:g}"]
        simulations[what] = Job(command + ["--simulate_only"] + keep_free, log.removesuffix(".log") + "_simulation.log",
                                lambda what=what: Steps.done(f"simulated the {what} in the background in "
                                                             f"{simulations[what].took()}"),
                                label=f"simulating the {what}", nice=10)
    Steps.done(f"simulating the {' and the '.join(simulations)} meanwhile, in the background "
               f"({', '.join(os.path.basename(j.log) for j in simulations.values())})")
    if final_done:
        pass
    elif training_db != db and not args.one_build_at_a_time:
        final_build = Job(build_command(args.protal, db, args.threads, *final_level), final_log, built_final_in_background,
                          label=f"building {os.path.basename(db)} in the background")
        Steps.done(f"building {os.path.basename(db)} meanwhile, in the background (index_and_package.log)")
    elif training_db == db:
        job = run(build_command(args.protal, db, args.threads, *final_level), final_log, built_final,
                  f"building {os.path.basename(db)}")
        Steps.done(built(db, job, remove_full_reference(db)))
    else:
        Steps.done(f"{os.path.basename(db)} is built after the models (--one-build-at-a-time)")
    if training_db != db and not training_done:
        # Read only for the training samples and the parity check: zstd level 3 packs it in a fraction of the
        # time of level 19 (which half of a build spent on), and loads as fast.
        job = run(build_command(args.protal, training_db, args.threads, "--compress_level", str(args.training_db_level)),
                  os.path.join(args.outdir, "training_db_index.log"), lambda: stages.mark("training_db", training_key),
                  f"building {os.path.basename(training_db)}")
        Steps.done(built(training_db, job, remove_full_reference(training_db) + files_took))
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
                       unpacked, "-t", str(args.threads)], os.path.join(args.outdir, "gene_ranking.log"),
                      label="unpacking the training database to rank its genes")
            # gene_congeners.tsv is a report of the build beside database.protal, not a member of it.
            congeners = os.path.join(training_db, "gene_congeners.tsv")
            if os.path.isfile(congeners):
                shutil.copy(congeners, unpacked)
            rows = rank_genes.rank(unpacked, gene_table, taxonomy)
            rank_genes.write_table(ranking_file, rows)
            shutil.copy(ranking_file, logs)
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
        take. With --scratch, its tables are copied to OUTDIR."""
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
            keep = os.path.join(args.outdir, os.path.basename(opts.out))
            os.makedirs(keep, exist_ok=True)
            for table in TABLES.values():
                if os.path.isfile(os.path.join(opts.out, table)):
                    shutil.copy(os.path.join(opts.out, table), keep)

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
            stream = ["--stream_above", f"{args.stream_above:g}"] if args.stream_above > 0 else []
            followers[what] = Job(command + ["--follow", "--profile_block", f"{args.profile_blocks:g}", "--protal_lock", lock]
                                  + stream,
                                  log, label=f"profiling the {what}")
        for what, command, log in collections_:
            opts = announce(what, command, log)
            Steps.done(f"profiled as its simulations go on, in protal runs of {args.profile_blocks:g} GB of reads or "
                       "more (taking turns with the other collection's); each design point's reads removed once profiled")
            report(opts, log, followers[what].finish())

    if args.profile_blocks > 0:
        follow_all()
    else:
        for what, command, log in collections_:
            collect(what, command, log)

    # One model per read type, trained in parallel.
    prefixes = {t: os.path.join(args.outdir, "trained_model" + ("" if t == "pe" else "_" + t)) for t in read_types}
    trainer_threads = max(1, args.threads // len(read_types))
    models = f"the {', '.join(read_types)} model{'s' if len(read_types) > 1 else ''}"
    Steps.start(f"training {models} (classifier_training*.log){' in parallel' if len(read_types) > 1 else ''}, "
                f"{trainer_threads} thread{'s' if trainer_threads > 1 else ''} each")
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
        trainers[t] = Job(command, os.path.join(args.outdir, "classifier_training" + ("" if t == "pe" else "_" + t) + ".log"),
                          label=f"training the {t} model")
    seconds = max(job.finish().seconds for job in trainers.values())
    each = f" ({', '.join(f'{t} {clock(job.seconds)}' for t, job in trainers.items())})" if len(trainers) > 1 else ""
    Steps.done(f"trained in {clock(seconds)}{each}; {model_scores(read_types, prefixes)}")
    for t, name, why in feature_choices(read_types, prefixes):
        Steps.done(f"{t} model's features: {name}: {why}")
    if selected:
        Steps.done(scenario_scores(read_types, prefixes, selected))
    # protal must score as the trainer does, and compute the features as it did during collection.
    Steps.start("checking that protal scores the models as the trainer does (parity*.log)")
    began = time.time()
    for t in read_types:
        run([sys.executable, PARITY, "--db", training_db, "--model", prefixes[t] + ".xml", "--training", training,
             "--read_type", t, "--protal", args.protal, "-t", str(args.threads)],
            os.path.join(args.outdir, "parity" + ("" if t == "pe" else "_" + t) + ".log"),
            label=f"checking the {t} model's parity with protal")
    Steps.done(f"{', '.join(read_types)}: the same probabilities and features, checked in {clock(time.time() - began)}")
    for t in read_types:
        prefix = prefixes[t]
        parity = os.path.join(training, "parity" if t == "pe" else "parity_" + t, "parity.txt")
        for name in (prefix + ".report.txt", prefix + ".metrics.json", prefix + ".thresholds.tsv", prefix + ".varimp.tsv",
                     prefix + ".predictions.tsv.gz", prefix + ".test_predictions.tsv.gz",
                     prefix + ".scenario_predictions.tsv.gz", prefix + ".calls.tsv.gz"):
            if os.path.isfile(name):
                shutil.copy(name, logs)
        if os.path.isfile(parity):
            shutil.copy(parity, os.path.join(logs, "parity.txt" if t == "pe" else f"parity_{t}.txt"))
    for name in [os.path.join(args.outdir, n) for n in ("training_data_simulation.log", "training_data.log",
                                                         "test_data_simulation.log", "test_data.log", "genome_table.txt")] + [heldout]:
        if os.path.isfile(name):
            shutil.copy(name, logs)
    # The genomes' contig names, read once for both (beside the samples).
    contig_cache = os.path.join(samples_root, "genome_contigs.tsv.gz")
    if "pe" in read_types and training_db != db and os.path.isfile(heldout):
        trace_relatives(training, training_db, heldout, logs, args.outdir, args.threads, contig_cache)
    if args.error_units:
        error_reads(args.error_units, prefixes, training, test if has_test else None, training_db, logs, args.outdir,
                    args.threads, contig_cache)
    Steps.start(f"adding {models} to {os.path.basename(db)} (final_package.log)")
    if final_build is not None:
        if final_build.seconds is None:
            Steps.done(f"waiting for {os.path.basename(db)}'s build in the background (index_and_package.log)")
        final_build.finish()  # its line comes from built_final_in_background
    elif training_db != db and not final_done:
        Steps.done(f"building {os.path.basename(db)} first (index_and_package.log)")
        job = run(build_command(args.protal, db, args.threads, *final_level), final_log, built_final,
                  f"building {os.path.basename(db)}")
        Steps.done(built(db, job, remove_full_reference(db)))
    # The trained models replace the shipped one and the placeholders in database.protal; --add_model checks
    # each and replaces them in place at the end of the file (seconds; without placeholders, a read type's
    # model is a new member and the ~20 GB file is rewritten once instead).
    began = time.time()
    run([args.protal, "--add_model", ",".join(prefixes[t] + ".xml" for t in read_types), "--read_type", ",".join(read_types),
         "--db", db, "-t", str(args.threads)], os.path.join(args.outdir, "final_package.log"), label=f"adding {models}")
    Steps.done(f"added in {clock(time.time() - began)}{db_size(db)}")
    versions, changed = versions_at_end(args.versions_at_start, tool_versions(args))
    with open(os.path.join(db, "build_metadata.tsv"), "w") as fh:
        fh.write("".join(f"{k}\t{v}\n" for k, v in provenance(args, versions, release, genome_table, heldout, n_heldout,
                                                                  read_types, prefixes, genes_note, insilico_note)))
    if changed:
        say(f"Warning: {', '.join(k.replace('_', ' ') for k in changed)} changed during the run, so its later steps ran "
            "other scripts or binaries than its first (build_metadata.tsv has both): rebuild the database if that matters")
    shutil.copy(os.path.join(db, "build_metadata.tsv"), logs)
    if os.path.isfile(os.path.join(db, "gene_congeners.tsv")):
        shutil.copy(os.path.join(db, "gene_congeners.tsv"), logs)
    if os.path.isfile(os.path.join(db, "gene_incongruence.tsv")):
        shutil.copy(os.path.join(db, "gene_incongruence.tsv"), logs)
    summary = summary_lines(read_types, prefixes, db)
    with open(os.path.join(logs, "summary.txt"), "w") as fh:
        fh.write("\n".join(summary) + "\n")
    print("\n" + "\n".join(summary) + "\n", flush=True)
    if Job.scratch:
        Job.scratch.look()
        say(f"The run took at most {gigabytes(Job.scratch.peak)} on {samples_root}; the simulated samples there "
            f"({gigabytes(tree_size(training) + tree_size(test))}) are left for a rerun")
    say(f"Ready protal database: {db}{db_size(db)}" +
        (f" (marker genes: {genes_note.split(',')[0].split(' (')[0]}, gene_subset.txt)" if subset else "") +
        f"; model evaluation: {logs} (start with trained_model.report.txt, and trained_model_<read type>.report.txt)")


if __name__ == "__main__":
    for sig in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP):
        signal.signal(sig, on_signal)
    try:
        main()
    finally:  # an uncaught exception, too
        stop_jobs()
