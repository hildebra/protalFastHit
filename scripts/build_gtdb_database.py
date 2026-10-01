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
  database (OUT_DIR/training_db) that leaves whole clades of every rank
  (--holdout-clades) and --holdout of the other species out
  (OUT_DIR/heldout_species.txt): their reads land on relatives, as those of
  organisms GTDB lacks do in real samples, and the model learns to reject those
  relatives. The report gives false positive and false negative rates by rank.
  The finished database has all species and the model trained so. The training
  database costs a second index build (while the first one runs) and its disk
  space.

One model per read type (--read-types, default pe,se,pb,ont), trained in
parallel: paired-end reads (ART), their first reads alone (single-end), and
PacBio and Nanopore reads of the same communities (pbsim3). Besides the training
data, an independent test set of another design (--test-*: other depths,
community sizes, abundances and strain mixes) is profiled and scored by each
model: cross-validation on the training data cannot show what its design lacks.

OUT_DIR/model_logs/ collects what tells whether the models are good: summary.txt
(TP, FP, TN, FN, sensitivity, specificity, precision and F1 of each model), each read
type's training report (how it does on species and clades it was not trained on,
on the independent test set, false positive and false negative rates by rank,
against the previous model and training procedure), its numbers as JSON, the
per-taxon predictions, the threshold table, the parity check with protal, the
genome table summary and build_metadata.tsv (what the database was built from
and with).

The tools the run needs are checked before it starts. A run that stops (a failure,
SIGTERM, Ctrl-C) stops every command it started, and one that fails in the
background (the finished database's build) stops the run within seconds. A rerun
into the same OUTDIR resumes: the conversion and the two index builds are skipped
when their inputs are those of the run that completed them (OUTDIR/.stages), and
the collector reuses the samples it simulated and profiled with the same database,
protal and design.

Each stage writes to a log of its own in OUTDIR. On the console, each line has the
time and how long the run has taken; a stage says when it starts and when it ends
(how long it took, its peak memory), and every --progress-every seconds each stage
running says how long it has run, the memory it takes and the last line of its log.
"""
import argparse
import collections
import glob
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
TRAINER = os.path.join(HERE, "random_forest_cmdline.py")
COLLECTOR = os.path.join(HERE, "collect_training_data.py")
PARITY = os.path.join(HERE, "check_model_parity.py")
# Read types whose models get knobs by sample depth (random_forest_cmdline.py --depth-knobs): long reads, where they
# raised the test F1 (docs/claude/2026-10-01-f1-opportunities); short reads gained nothing.
DEPTH_KNOB_READ_TYPES = ("pb", "ont")
ACCESSION = re.compile(r"(?:RS_|GB_)?(GC[AF]_\d{9}\.\d+)")
sys.path.insert(0, os.path.join(HERE, "mini_db"))
sys.path.insert(0, HERE)
import lineages  # noqa: E402
from gtdb_to_protal_db import (clear_build_outputs, full_reference_path, normalize_accession,  # noqa: E402
                               read_representatives, remove_full_reference as remove_full_reference_files)
from model_pmml import MODEL_FILES, write_placeholder  # noqa: E402
from collect_training_data import TABLES, clock, last_line, units_of, parse_args as collector_args  # noqa: E402

STARTED = time.time()


def say(message):
    """Prints a message, its first line headed by the time and how long the run has taken."""
    print(f"[{time.strftime('%H:%M:%S')} +{clock(time.time() - STARTED)}] {message}", flush=True)


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
        return (f"{self.path}: {gigabytes(max(0, usage.used - self.start))} in use by the run (at most "
                f"{gigabytes(self.peak)} so far), {gigabytes(usage.free)} free")


class Job:
    """A command run with its output to log, in a process group of its own: when the script stops, however
    it stops (a failure, an exception, SIGTERM, SIGINT, SIGHUP), the command is stopped with what it started
    (the collector's simulator, ART and protal runs), so that a rerun does not race a build left running.
    While the script waits for one job, it looks at the others every few seconds: one that failed (the
    finished database's build in the background) stops the script then, not hours later. Every
    progress_every seconds it says how each job is doing."""
    running = []
    progress_every = 600
    scratch = None  # a Scratch with --scratch

    def __init__(self, command, log, on_success=None, label=None):
        os.makedirs(os.path.dirname(log), exist_ok=True)
        self.command, self.log, self.on_success, self.started = command, log, on_success, time.time()
        self.label = label or os.path.basename(log).removesuffix(".log")
        self.seconds = None  # set when it has ended
        self.peak = None  # the most memory it, or a command it ran, took, in bytes; set when it has ended
        self.fh = open(log, "w")
        self.process = subprocess.Popen(command, stdout=self.fh, stderr=subprocess.STDOUT, start_new_session=True)
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


def make_genome_table(gtdb, release, output, extra_dirs=(), species=None):
    """The simulator's genome table: every genome FASTA of a GTDB species found in the release (and in
    extra_dirs), or only those of the species in the set `species`."""
    taxonomy = {}
    for domain in ("bac120", "ar53"):
        for suffix in (".tsv", ".tsv.gz"):
            path = os.path.join(gtdb, f"{domain}_taxonomy_r{release}{suffix}")
            if os.path.isfile(path):
                import gzip
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
    with open(output, "w") as fh:
        for accession in sorted(paths):
            fh.write(f"{accession}\t{taxonomy[accession]}\t{os.path.abspath(paths[accession])}\n")
    return len(paths)


def summarize_genome_table(path, reps):
    """What the simulations can draw: genomes and species by domain, and how often a simulated species is
    not the database's representative genome (the simulator picks a species, then one of its genomes)."""
    genomes = collections.Counter()
    domain = {}
    rep_genomes = 0
    with open(path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            lineage = next((f for f in fields if f.startswith("d__") and ";s__" in f), None)
            if lineage is None:
                continue
            species = lineage.split(";")[-1]
            genomes[species] += 1
            domain[species] = lineage.split(";")[0][3:]
            rep_genomes += reps is not None and normalize_accession(fields[0]) in reps
    lines = [f"genome table {path}: {sum(genomes.values())} genomes of {len(genomes)} species"]
    for d in sorted(set(domain.values())):
        sp = [s for s in genomes if domain[s] == d]
        lines.append(f"  {d}: {len(sp)} species, {sum(genomes[s] for s in sp)} genomes")
    several = sum(1 for n in genomes.values() if n > 1)
    other_strain = sum((n - 1) / n for n in genomes.values()) / max(1, len(genomes))
    lines.append(f"  {several} species have more than one genome; a simulated species is another genome than "
                 f"its representative {100 * other_strain:.1f}% of the time")
    if reps is not None:
        lines.append(f"  {rep_genomes} of the genomes are species representatives (the database's references)")
    if other_strain < 0.2:
        lines.append("  WARNING: nearly all simulated species will be the database's own reference genome, closer to it "
                     "than real strains are. Add non-representative genomes (scripts/download_gtdb.py, then "
                     "--inputs) to train on real strain divergence.")
    return lines


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


def gene_neighbours_summary(build_log):
    """What protal --build said of the gene neighbours it packed (its "Gene neighbours:" line), for
    build_metadata.tsv."""
    try:
        with open(build_log) as fh:
            lines = [line.strip() for line in fh if line.startswith("Gene neighbours:")]
    except OSError:
        return "unknown (no build log)"
    if not lines:
        return "none (this protal does not pack them)"
    return lines[-1].removeprefix("Gene neighbours:").strip()


def provenance(args, release, genome_table, heldout, n_heldout, read_types, prefixes):
    """build_metadata.tsv: what the database was built from and with, so that two builds can be compared."""
    def output(command):
        try:
            return subprocess.run(command, capture_output=True, text=True, timeout=60).stdout.strip()
        except (OSError, subprocess.SubprocessError):
            return ""
    version = output([args.protal, "--version"]).splitlines()
    commit = output(["git", "-C", HERE, "rev-parse", "HEAD"])
    if commit and output(["git", "-C", HERE, "status", "--porcelain", "--", "."]):
        commit += " (scripts changed since)"
    genomes, species = 0, set()
    with open(genome_table) as fh:
        for line in fh:
            lineage = next((f for f in line.rstrip("\n").split("\t") if f.startswith("d__") and ";s__" in f), None)
            if lineage:
                genomes += 1
                species.add(lineage.split(";")[-1])
    clade_counts = collections.Counter(rank for rank, _ in set(read_holdout(heldout).values())) if n_heldout else {}
    rows = [("gtdb_release", f"r{release}"), ("built", time.strftime("%Y-%m-%d %H:%M:%S")),
            ("protal_version", version[-1] if version else "unknown"), ("protal_binary", args.protal),
            ("scripts_commit", commit or "unknown (not a git checkout)"), ("command", " ".join(sys.argv)),
            ("seed", args.seed), ("genome_table", f"{genomes} genomes of {len(species)} species"),
            ("gene_conservation", gene_conservation_summary(os.path.join(args.outdir, "index_and_package.log"))),
            ("gene_neighbours", gene_neighbours_summary(os.path.join(args.outdir, "index_and_package.log"))),
            ("classifier_features", "normalized"), ("classifier_trees", args.ntree),
            ("classifier_max_leaves", args.maxnodes), ("classifier_training_species_left_out", n_heldout)]
    rows += [(f"classifier_training_{rank}_clades_left_out", clade_counts[rank]) for rank in CLADE_RANKS
             if clade_counts.get(rank)]
    rows += [("classifier_training_samples", f"{args.samples} per design point"),
             ("classifier_training_design", f"read pairs {args.read_pairs}; read setups {args.read_setups}; "
                                            f"species per sample {args.species_per_sample}; strains "
                                            f"{args.strains_per_species or 'one'}; abundance {args.abundance or 'default'}"),
             ("classifier_read_types", ",".join(read_types)),
             ("classifier_depth_knobs", ",".join(t for t in read_types if t in DEPTH_KNOB_READ_TYPES) or "none")]
    for t in read_types:
        try:
            with open(prefixes[t] + ".metrics.json") as fh:
                metrics = json.load(fh)
        except (OSError, ValueError):
            continue
        species_cv = metrics.get("evaluation", {}).get("species", {})
        test = metrics.get("test", {}).get("this one", {})
        rows.append((f"model_{t}", f"species held out F1 {species_cv.get('F1')}, FP per sample "
                                   f"{species_cv.get('FP_per_sample')}; independent test F1 {test.get('F1')}, FP per "
                                   f"sample {test.get('FP_per_sample')}"))
        knobs = metrics.get("depth_knobs", {}).get("knobs")
        if knobs:
            rows.append((f"model_{t}_depth_knobs", ",".join(f"{b}:{k:g}" for b, k in sorted(knobs.items())) +
                         f"; independent test F1 at them {metrics.get('test_depth_knobs', {}).get('F1')}"))
    return rows


def summary_lines(read_types, prefixes, db):
    """model_logs/summary.txt: for each read type's model, TP, FP, TN, FN and the rates, with species held
    out in training (cross-validation) and on the independent test set, from its .metrics.json."""
    header = ("read type", "evaluated on", "taxa", "TP", "FP", "TN", "FN", "sensitivity", "specificity",
              "precision", "F1", "FP/sample")
    rows = []
    for t in read_types:
        try:
            with open(prefixes[t] + ".metrics.json") as fh:
                metrics = json.load(fh)
        except (OSError, ValueError):
            rows.append((t, "no metrics (training failed?)") + ("",) * (len(header) - 2))
            continue
        for label, m in (("species held out", metrics.get("evaluation", {}).get("species")),
                         ("independent test set", metrics.get("test", {}).get("this one"))):
            if not m:
                continue
            tp, fn, fp = m["present"] - m["FN"], m["FN"], m["FP"]
            tn = m["taxa"] - m["present"] - fp
            rate = lambda v: "-" if v is None else f"{v:.4f}"
            rows.append((t, label, str(m["taxa"]), str(tp), str(fp), str(tn), str(fn), rate(m.get("sensitivity")),
                         rate(tn / (tn + fp) if tn + fp else None), rate(m.get("precision")), rate(m.get("F1")),
                         "-" if m.get("FP_per_sample") is None else f"{m['FP_per_sample']:.2f}"))
    widths = [max(len(str(r[i])) for r in [header, *rows]) for i in range(len(header))]
    table = ["  ".join(str(v).ljust(w) for v, w in zip(r, widths)).rstrip() for r in [header, *rows]]
    return [f"Presence models of {db}, taxa scored at knob 0.5: TP present and called, FP absent and called, TN "
            "absent and not called, FN present and not called. species held out: each taxon scored by forests "
            "that did not see its species; independent test set: samples of another design, scored by the "
            "final model. Details: trained_model*.report.txt", ""] + table


def remove_full_reference(folder):
    """Removes a folder's full reference (full_reference.fna.zst) once its build is done: the marker genes of
    every genome (86 GB raw at r226), which only that build reads. A rebuild converts the release again (the
    build packed reference.fna into database.protal)."""
    path = full_reference_path(folder)
    if path:
        say(f"Removed {path} ({gigabytes(remove_full_reference_files(folder))}): only the build reads it")


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
    for folder in (gtdb, os.path.join(gtdb, "genomic_files_reps"), os.path.join(gtdb, "genomic_files_all")):
        for name in sorted(os.listdir(folder)) if os.path.isdir(folder) else ():
            path = os.path.join(folder, name)
            if f"_r{release}" in name and (os.path.isfile(path) or "_marker_genes_" in name):
                st = os.stat(path)
                items.append([os.path.relpath(path, gtdb), st.st_size if os.path.isfile(path) else len(os.listdir(path)),
                              st.st_mtime_ns])
    return items


def check_tools(args, read_types):
    """What the run needs later, checked before it starts: a missing tool would otherwise stop it after the
    conversion, the builds or the collection."""
    problems = []

    def executable(command, what, hint):
        if shutil.which(command) or (os.path.isfile(command) and os.access(command, os.X_OK)):
            return True
        problems.append(f"{what} ({command}) is not there: {hint}")
        return False

    if executable(args.protal, "protal", "build it (docs/installation.md), or pass --protal"):
        version = subprocess.run([args.protal, "--version"], capture_output=True, text=True)
        if version.returncode:
            problems.append(f"{args.protal} --version failed ({version.returncode}): {version.stdout[-300:]}{version.stderr[-300:]}")
    executable(args.simulator, "the simulator", "it is built with protal (target simulate_metagenomes), or pass --simulator")
    executable("art_illumina", "ART", "the simulator simulates the Illumina reads with it: install ART (conda: art, "
                                      "envs/protal-db-build.yaml)")
    if any(t in ("pb", "ont") for t in read_types):
        executable(args.pbsim, "pbsim3", "pb and ont reads are simulated with it: install it (e.g. micromamba install "
                                         "-c conda-forge -c bioconda pbsim3), pass --pbsim, or leave pb and ont out of "
                                         "--read-types")
    imports = subprocess.run([sys.executable, "-c", "import joblib, numpy, pandas, sklearn"], capture_output=True, text=True)
    if imports.returncode:
        problems.append(f"{sys.executable} cannot import what the trainer needs ({imports.stderr.strip().splitlines()[-1]}): "
                        "install scikit-learn, joblib, numpy and pandas (envs/protal-db-build.yaml), or run this "
                        "script with a Python that has them")
    if problems:
        sys.exit("Cannot start:\n  " + "\n  ".join(problems))
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
    p.add_argument("--genome-table", help="optional simulator table: accession, taxonomy, whole genome FASTA path")
    p.add_argument("--extra-genomes", action="append", default=[],
                   help="folder of more whole genomes of GTDB species, found by the accession in their file names "
                        "(e.g. the NCBI genomes of download_gtdb.py); repeatable")
    p.add_argument("--holdout", type=float, default=0.2,
                   help="fraction of the species (of those the genome table can simulate and no clade of "
                        "--holdout-clades took) left out of a separate training database (default 0.2): their "
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
                   help="zstd level of the training database (default 3; the finished database uses protal's default)")
    p.add_argument("--simulate-species",
                   help="file of the species to simulate from (e.g. simulation_species.txt of download_gtdb.py: "
                        "species with other strains, and some without): the simulator draws species uniformly, so "
                        "among all ~130,000 GTDB species the few with downloaded strains would hardly be drawn")
    p.add_argument("--protal", default="protal", help="protal executable")
    p.add_argument("--simulator", default="simulate_metagenomes", help="simulate_metagenomes executable")
    p.add_argument("-t", "--threads", type=int, default=8)
    p.add_argument("--samples", type=int, default=12,
                   help="samples per design point (default 12; on a GTDB-like world the model still improved "
                        "from 60 to 120 samples)")
    p.add_argument("--read-pairs", default="1000,5000,20000,100000,500000",
                   help="read pairs per sample, one design point each (default 1000,5000,20000,100000,500000: "
                        "without the shallowest, a model missed 8%% of the present taxa of samples of 1000 read "
                        "pairs)")
    p.add_argument("--read-setups", default="100:HS20:300:40,150:HSXt:350:50,250:MSv3:550:50",
                   help="LENGTH:ART_PROFILE:FRAGMENT_MEAN:FRAGMENT_SD, one design point each (HSXt: HiSeq X, the "
                        "closest of ART's profiles to NovaSeq; file=R1.txt+R2.txt: profiles art_profiler_illumina "
                        "made from real reads)")
    p.add_argument("--archaea", type=int, default=2)
    p.add_argument("--species-per-sample", default="20-200",
                   help="species per sample, drawn per sample (default 20-200: real gut samples hold 100-300 GTDB "
                        "species with a long tail of rare ones)")
    p.add_argument("--strains-per-species", default="0.3,0.1",
                   help="probabilities of a second, third, ... strain of a species in a sample (default 0.3,0.1): "
                        "real samples often mix strains, which changes the allele-frequency features")
    p.add_argument("--abundance", default="",
                   help="abundance model of the training samples: lognormal:SIGMA, powerlaw:ALPHA or negbin:R:P "
                        "(default: the simulator's, Poisson-lognormal with sigma 1.3)")
    p.add_argument("--read-types", default="pe,se,pb,ont",
                   help="the read types to train a model for (default pe,se,pb,ont): se from the paired-end "
                        "samples' first reads, pb and ont from long reads of the same communities (pbsim3). The "
                        "models are trained in parallel; read types left out keep placeholders")
    p.add_argument("--long-read-bases", default="300000,1500000,6000000,30000000,150000000",
                   help="bases per long-read sample, one design point each (collect_training_data.py)")
    p.add_argument("--pb-setup", default="errhmm:ERRHMM-SEQUEL:15000:3000:0.999",
                   help="pbsim3 METHOD:MODEL:LENGTH_MEAN:LENGTH_SD:ACCURACY_MEAN of PacBio reads")
    p.add_argument("--ont-setup", default="qshmm:QSHMM-ONT-HQ:8000:6000:0.97:39/24/36",
                   help="pbsim3 METHOD:MODEL:LENGTH_MEAN:LENGTH_SD:ACCURACY_MEAN of Nanopore reads")
    p.add_argument("--pbsim", default="pbsim", help="pbsim3 executable, for pb and ont")
    p.add_argument("--pbsim-models", help="folder of pbsim3's .model files (default: found next to the executable)")
    p.add_argument("--test-samples", type=int, default=4,
                   help="samples per design point of the independent test set (default 4; 0: none). The test set "
                        "has another design than the training data (--test-*), so that the report shows what "
                        "cross-validation on the training data cannot")
    p.add_argument("--test-read-pairs", default="500,2000,10000,50000,200000,1000000")
    p.add_argument("--test-species-per-sample", default="10-300")
    p.add_argument("--test-abundance", default="lognormal:2.0", help="(default lognormal:2.0: more uneven than training)")
    p.add_argument("--test-strains-per-species", default="0.5,0.2")
    p.add_argument("--test-long-read-bases", default="150000,1000000,5000000,25000000,250000000")
    p.add_argument("--congeners", type=int, default=0,
                   help="species of one genus in every sample of a design point (collect_training_data.py "
                        "--congeners): relatives that share a sample, as in real samples")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--ntree", type=int, default=64)
    p.add_argument("--maxnodes", type=int, default=128)
    p.add_argument("--no-gene-neighbours", action="store_true",
                   help="do not record which marker genes lie next to which in the genomes to simulate from "
                        "(gene_neighbours.py; protal then pairs no mates across neighbouring genes)")
    p.add_argument("--no-placeholder-models", action="store_true",
                   help="leave the se, pb and ont models out of the database (a run with such reads then stops "
                        "with an error) instead of placeholders that report no species until trained ones replace "
                        "them (placeholder_models.py)")
    p.add_argument("--evaluation", choices=["full", "basic", "none"], default="full",
                   help="how much the trainer evaluates (random_forest_cmdline.py --evaluation)")
    p.add_argument("--progress-every", type=float, default=600,
                   help="seconds between the status lines of the stages running: how long each has run, the memory "
                        "it takes and the last line of its log (default 600; 0: none)")
    p.add_argument("--scratch",
                   help="a fast local disk (a compute node's own) for the simulated samples: their reads, alignments, "
                        "profiles and the simulators' temporary files go to SCRATCH/training and SCRATCH/test, and "
                        "only the tables to OUTDIR. A network file system (OUTDIR's, often) is slow at the many "
                        "files the simulators write and delete. With the defaults the run takes about 60 GB there at "
                        "its peak (docs/building-a-database.md); a rerun reuses the samples in the same SCRATCH")
    args = p.parse_args()
    Job.progress_every = args.progress_every
    read_types = [t.strip() for t in args.read_types.split(",") if t.strip()]
    if not read_types or any(t not in TABLES for t in read_types):
        p.error(f"--read-types: a comma-separated list of {', '.join(TABLES)}, got {args.read_types!r}")
    check_tools(args, read_types)
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
    say(f"A protal database of GTDB r{release} in {args.outdir}, with {args.threads} threads" +
        (f"; the stages running are reported every {clock(Job.progress_every)}" if Job.progress_every > 0 else ""))
    samples_root = args.outdir  # where the collections simulate and profile their samples
    if args.scratch:
        samples_root = os.path.abspath(args.scratch)
        os.makedirs(samples_root, exist_ok=True)
        Job.scratch = Scratch(samples_root)
        say(f"The simulated samples go to {samples_root}: {gigabytes(shutil.disk_usage(samples_root).free)} free")
    genome_table = args.genome_table or os.path.join(args.outdir, "genomes.tsv")
    if not args.genome_table:
        pool = None
        if args.simulate_species:
            with open(args.simulate_species) as fh:
                pool = {n if n.startswith("s__") else "s__" + n
                        for n in (line.rstrip("\n").split("\t")[0].strip() for line in fh) if n and not n.startswith("#")}
        count = make_genome_table(args.gtdb, release, genome_table, args.extra_genomes, pool)
        say(f"Found {count} genome FASTAs for training")
    elif args.extra_genomes or args.simulate_species:
        sys.exit("--extra-genomes and --simulate-species shape the genome table this script makes; "
                 "apply them to --genome-table instead")
    logs = os.path.join(args.outdir, "model_logs")
    os.makedirs(logs, exist_ok=True)
    summary = summarize_genome_table(genome_table, read_representatives(args.gtdb, release))
    with open(os.path.join(args.outdir, "genome_table.txt"), "w") as fh:
        fh.write("\n".join(summary) + "\n")
    say("\n".join(summary))

    # A rerun skips what an earlier run into OUTDIR completed with the same inputs (see Stages).
    stages = Stages(os.path.join(args.outdir, ".stages"))
    convert_key = {"converter": content_hash(CONVERTER), "release": release, "gtdb": release_identity(args.gtdb, release),
                   "placeholders": not args.no_placeholder_models,
                   "gene_neighbours": None if args.no_gene_neighbours else [content_hash(GENE_NEIGHBOURS), content_hash(genome_table)]}
    final_key = {"convert": convert_key, "protal": file_identity(args.protal)}
    final_done = stages.done("protal_db", final_key) and os.path.isfile(os.path.join(db, "database.protal"))
    # --build packs the taxonomy into database.protal; the collector and the trainer read it (domains,
    # representative genomes). The training database has the same taxonomy.
    taxonomy = os.path.join(args.outdir, "internal_taxonomy.dmp")

    def convert(into):
        """The converted files of the release in `into`, with the models of the read types but pe as
        placeholders (protal warns when it loads one), packed by --build like model_pe.xml."""
        log = os.path.join(args.outdir, "convert.log")
        say(f"Converting GTDB r{release} into {into} ({log})")
        job = run([sys.executable, CONVERTER, "--gtdb", args.gtdb, "--outdir", into, "--release", release, "-t",
                   str(args.threads)], log, label="converting the release")
        say(f"Converted GTDB r{release} in {job.took()}")
        shutil.copyfile(os.path.join(into, "internal_taxonomy.dmp"), taxonomy)
        if not args.no_gene_neighbours:
            # Which marker genes lie next to which, from the representatives among the genomes to simulate from
            # (gene_neighbours.tsv, which --build packs; a copy without some species keeps it, per clade).
            log = os.path.join(args.outdir, "gene_neighbours.log")
            job = run([sys.executable, GENE_NEIGHBOURS, "--db", into, "--genome_table", genome_table, "-t", str(args.threads)],
                      log, label="finding the genes' neighbours")
            say(f"Found the genes' neighbours in the genomes in {job.took()} ({log})")
        if not args.no_placeholder_models:
            for read_type in ("se", "pb", "ont"):
                write_placeholder(os.path.join(into, MODEL_FILES[read_type]), read_type)

    converted = None  # the folder with the converted files, once there
    if final_done:
        say(f"{db} was built by an earlier run from the same release and protal; kept")
        remove_full_reference(db)  # left by an earlier version of this script
    elif stages.done("convert", convert_key) and os.path.isfile(taxonomy) and \
            all(os.path.isfile(os.path.join(db, f)) for f in ("reference.fna", "reference.map", "internal_taxonomy.dmp")):
        # Converted by an earlier run that stopped before the build packed the files.
        clear_build_outputs(db)
        converted = db
        say(f"{db} holds the release converted by an earlier run; not converted again")
    else:
        stages.forget("protal_db")
        convert(db)
        stages.mark("convert", convert_key)
        converted = db
    if not os.path.isfile(taxonomy):
        convert(os.path.join(args.outdir, ".converted"))
        converted = os.path.join(args.outdir, ".converted")

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
    training_db, training_done = db, False
    n_heldout, holdout_lines = 0, []
    if os.path.exists(heldout):
        chosen = read_holdout(heldout)
        n_heldout = len(chosen)
        holdout_lines = [f"training database: {n_heldout} species left out ({heldout})"] + \
            describe_holdout(chosen, pool_species(genome_table))
        with open(os.path.join(logs, "holdout.txt"), "w") as fh:
            fh.write("\n".join(holdout_lines) + "\n")
        say("\n".join(holdout_lines))
        training_db = os.path.join(args.outdir, "training_db")
        training_key = {"convert": convert_key, "heldout": content_hash(heldout), "protal": final_key["protal"],
                        "level": args.training_db_level}
        training_done = stages.done("training_db", training_key) and \
            os.path.isfile(os.path.join(training_db, "database.protal"))
        if training_done:
            say(f"{training_db} was built by an earlier run with the same species left out; kept")
            remove_full_reference(training_db)
        else:
            stages.forget("training_db")
            if converted is None:  # the finished database's build consumed them
                converted = os.path.join(args.outdir, ".converted")
                convert(converted)
            job = run([sys.executable, CONVERTER, "--from_db", converted, "--exclude_species", heldout, "--outdir",
                       training_db, "-t", str(args.threads)], os.path.join(args.outdir, "training_db.log"),
                      label="leaving the species out")
            say(f"Wrote the files of {training_db}, without the species left out, in {job.took()}")
    if converted and converted != db:
        shutil.rmtree(converted, ignore_errors=True)

    # The finished database is needed only for --add_model at the end: it is built in the background from the
    # start, while the training database is built and the training data collected, unless one build at a time.
    final_log = os.path.join(args.outdir, "index_and_package.log")
    final_build = None

    def built_final():
        stages.mark("protal_db", final_key)
        remove_full_reference(db)

    def built_training():
        stages.mark("training_db", training_key)
        remove_full_reference(training_db)

    if final_done:
        pass
    elif training_db != db and not args.one_build_at_a_time:
        final_build = Job(build_command(args.protal, db, args.threads), final_log, built_final,
                          label=f"building {os.path.basename(db)} in the background")
        say(f"Building {db} in the background ({final_log})")
    elif training_db == db:
        say(f"Building {db} ({final_log})")
        job = run(build_command(args.protal, db, args.threads), final_log, built_final, f"building {os.path.basename(db)}")
        say(f"Built {db} in {job.took()}{db_size(db)}")
    if training_db != db and not training_done:
        # Read only for the training samples and the parity check: zstd level 3 packs it in a fraction of the
        # time of level 19 (which half of a build spent on), and loads as fast.
        log = os.path.join(args.outdir, "training_db_index.log")
        say(f"Building {training_db} ({log})")
        job = run(build_command(args.protal, training_db, args.threads, "--compress_level", str(args.training_db_level)),
                  log, built_training, f"building {os.path.basename(training_db)}")
        say(f"Built {training_db} in {job.took()}{db_size(training_db)}")

    # Training data of every read type (pe, se from its first reads, pb and ont from long reads of the same
    # communities), then an independent test set of another design, both against the training database.
    def collect_command(out, samples, read_pairs, species, abundance, strains, long_bases, seed):
        command = [sys.executable, COLLECTOR, "--db", training_db, "--genome_table", genome_table, "-o", out,
                   "--protal", args.protal, "--simulator", args.simulator, "--samples", str(samples),
                   "--read_pairs", read_pairs, "--read_setups", args.read_setups, "--archaea", str(args.archaea),
                   "--species_per_sample", species, "--seed", str(seed), "-t", str(args.threads),
                   "--taxonomy", taxonomy, "--congeners", str(args.congeners), "--read_types", ",".join(read_types),
                   "--long_read_bases", long_bases, "--pb_setup", args.pb_setup, "--ont_setup", args.ont_setup,
                   "--pbsim", args.pbsim]
        command += ["--abundance", abundance] if abundance else []
        command += ["--strains_per_species", strains] if strains else []
        command += ["--pbsim_models", args.pbsim_models] if args.pbsim_models else []
        if n_heldout:
            command += ["--novel_species", heldout, "--novel_clades", str(args.novel_clades_per_sample)]
        return command

    def collect(what, command, log):
        """Runs the collector, saying before what it makes and after what its tables hold and how much space its
        samples take. With --scratch, its tables are copied to OUTDIR."""
        opts = collector_args(command[2:])
        units = collections.Counter(unit["type"] for unit in units_of(opts)[1])
        say(f"Collecting {what}: {', '.join(f'{n * opts.samples} {t}' for t, n in units.items())} samples, "
            f"{opts.samples} per design point ({log})")
        job = run(command, log, label=f"collecting {what}")
        with open(log, errors="replace") as fh:
            tables = [line.rstrip() for line in fh if re.match(r"\d+ taxa in ", line)]
        say(f"Collected {what} in {job.took()}; {gigabytes(tree_size(opts.out))} in {opts.out}" +
            (f"; scratch {Job.scratch.text()}" if Job.scratch else "") + "".join(f"\n  {line}" for line in tables))
        if args.scratch:
            keep = os.path.join(args.outdir, os.path.basename(opts.out))
            os.makedirs(keep, exist_ok=True)
            for table in TABLES.values():
                if os.path.isfile(os.path.join(opts.out, table)):
                    shutil.copy(os.path.join(opts.out, table), keep)

    training = os.path.join(samples_root, "training")
    collect(f"the training data ({', '.join(read_types)})",
            collect_command(training, args.samples, args.read_pairs, args.species_per_sample, args.abundance,
                            args.strains_per_species, args.long_read_bases, args.seed),
            os.path.join(args.outdir, "training_data.log"))
    test = os.path.join(samples_root, "test")
    if args.test_samples > 0:
        collect("the independent test set",
                collect_command(test, args.test_samples, args.test_read_pairs, args.test_species_per_sample,
                                args.test_abundance, args.test_strains_per_species, args.test_long_read_bases,
                                args.seed + 1000), os.path.join(args.outdir, "test_data.log"))

    # One model per read type, trained in parallel.
    prefixes = {t: os.path.join(args.outdir, "trained_model" + ("" if t == "pe" else "_" + t)) for t in read_types}
    trainer_threads = max(1, args.threads // len(read_types))
    models = f"the {', '.join(read_types)} model{'s' if len(read_types) > 1 else ''}"
    say(f"Training {models}{' in parallel' if len(read_types) > 1 else ''}, {trainer_threads} "
        f"thread{'s' if trainer_threads > 1 else ''} each ({os.path.join(args.outdir, 'classifier_training*.log')})")
    trainers = {}
    for t in read_types:
        command = [sys.executable, TRAINER, "--truth-file", os.path.join(training, TABLES[t]),
                   "--output-prefix", prefixes[t], "--features", "normalized", "--ntree", str(args.ntree),
                   "--maxnodes", str(args.maxnodes), "--seed", str(args.seed), "--threads", str(trainer_threads),
                   "--taxonomy", taxonomy, "--evaluation", args.evaluation]
        if t in DEPTH_KNOB_READ_TYPES:
            command += ["--depth-knobs"]
        if args.test_samples > 0 and os.path.isfile(os.path.join(test, TABLES[t])):
            command += ["--test-file", os.path.join(test, TABLES[t])]
        trainers[t] = Job(command, os.path.join(args.outdir, "classifier_training" + ("" if t == "pe" else "_" + t) + ".log"),
                          label=f"training the {t} model")
    seconds = max(job.finish().seconds for job in trainers.values())
    if len(trainers) > 1:
        say(f"Trained {models} in parallel in {clock(seconds)}" +
            "".join(f"\n  {t}: {job.took()}" for t, job in trainers.items()))
    else:
        say(f"Trained {models} in {trainers[read_types[0]].took()}")
    # protal must score as the trainer does, and compute the features as it did during collection.
    for t in read_types:
        job = run([sys.executable, PARITY, "--db", training_db, "--model", prefixes[t] + ".xml", "--training", training,
                   "--read_type", t, "--protal", args.protal, "-t", str(args.threads)],
                  os.path.join(args.outdir, "parity" + ("" if t == "pe" else "_" + t) + ".log"),
                  label=f"checking the {t} model's parity with protal")
        say(f"protal scores the {t} model as the trainer does ({job.took()})")
    for t in read_types:
        prefix = prefixes[t]
        parity = os.path.join(training, "parity" if t == "pe" else "parity_" + t, "parity.txt")
        for name in (prefix + ".report.txt", prefix + ".metrics.json", prefix + ".thresholds.tsv", prefix + ".varimp.tsv",
                     prefix + ".predictions.tsv.gz", prefix + ".test_predictions.tsv.gz"):
            if os.path.isfile(name):
                shutil.copy(name, logs)
        if os.path.isfile(parity):
            shutil.copy(parity, os.path.join(logs, "parity.txt" if t == "pe" else f"parity_{t}.txt"))
    for name in (os.path.join(args.outdir, "training_data.log"), os.path.join(args.outdir, "test_data.log"),
                 os.path.join(args.outdir, "genome_table.txt"), heldout):
        if os.path.isfile(name):
            shutil.copy(name, logs)
    if final_build is not None:
        if final_build.seconds is None:
            say(f"Waiting for the build of {db} in the background ({final_log})")
        say(f"Built {db} in the background in {final_build.finish().took()}{db_size(db)}")
    elif training_db != db and not final_done:
        say(f"Building {db} ({final_log})")
        job = run(build_command(args.protal, db, args.threads), final_log, built_final, f"building {os.path.basename(db)}")
        say(f"Built {db} in {job.took()}{db_size(db)}")
    # The trained models replace the shipped one and the placeholders in database.protal; --add_model checks
    # each and copies the other parts as they are.
    say(f"Adding {models} to {db}")
    for t in read_types:
        run([args.protal, "--add_model", prefixes[t] + ".xml", "--read_type", t, "--db", db, "-t", str(args.threads)],
            os.path.join(args.outdir, "final_package" + ("" if t == "pe" else "_" + t) + ".log"),
            label=f"adding the {t} model")
    with open(os.path.join(db, "build_metadata.tsv"), "w") as fh:
        fh.write("".join(f"{k}\t{v}\n" for k, v in provenance(args, release, genome_table, heldout, n_heldout,
                                                                  read_types, prefixes)))
    shutil.copy(os.path.join(db, "build_metadata.tsv"), logs)
    summary = summary_lines(read_types, prefixes, db)
    with open(os.path.join(logs, "summary.txt"), "w") as fh:
        fh.write("\n".join(summary) + "\n")
    say("\n".join(summary))
    if Job.scratch:
        Job.scratch.look()
        say(f"The run took at most {gigabytes(Job.scratch.peak)} on {samples_root}; the simulated samples there "
            f"({gigabytes(tree_size(training) + tree_size(test))}) are left for a rerun")
    say(f"Ready protal database: {db}{db_size(db)}\nModel evaluation: {logs} (start with trained_model.report.txt, and "
        "trained_model_<read type>.report.txt)")


if __name__ == "__main__":
    for sig in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP):
        signal.signal(sig, on_signal)
    try:
        main()
    finally:  # an uncaught exception, too
        stop_jobs()
