#!/usr/bin/env python3
"""test_protal_e2e.py - end-to-end checks of protal and simulate_metagenomes on the mini database.

Simulates reads from the database's reference genes (sequencing errors only, fixed seeds) and runs the real
binaries, checking what the unit tests cannot: exit codes, output files, SAM records, strain MSAs, reruns, that
failures are reported and stay with their sample, a profile's truth counts and abundances, that reads of nothing in
the database give an empty profile, and a gradient-boosted model scored end to end. examples/mini_db/run.sh checks
accuracy on reads of whole genomes.

  just e2e                                  # builds the mini DB and the binaries first
  PROTAL_TEST_DB=data/mini_db/protal_db PROTAL=build/protal SIMULATE=build/simulate_metagenomes \
      PROTAL_TESTS_REQUIRED=1 python3 -m unittest -v tests/e2e/test_protal_e2e.py

The tests are written for the mini database of scripts/mini_db/build_mini_db.sh: they name its three species
(Mockella alpha and beta, two congeners, and Fakibacter gamma), their taxids and genes. They need Linux, the zstd CLI
(or Python 3.14), and for a test or two numpy.

PROTAL_TEST_DB         protal database (only read): the single file database.protal (or its folder), or separate raw
                       or zstd-compressed files (index.prx.zst, reference.fna.zst). Tests that read the database's
                       files get them unpacked (protal --unpack_db) into a temporary folder. Without it the tests that
                       need it are skipped; VersionTest, SimulatorTest, QcmsaContractTest, BuildIndexTest and
                       GeneNeighboursTest need none.
PROTAL                 protal binary (default: build/protal)
SIMULATE               simulate_metagenomes binary (default: build/simulate_metagenomes)
PROTAL_TESTS_REQUIRED  1: a missing prerequisite (the database, a binary, the zstd CLI, numpy) fails
                       the tests that need it instead of skipping them, as CI wants
PROTAL_TEST_KEEP       set: keep the temporary folders

What it costs: about 1.5 minutes on 4 cores, 3 GB of memory and of /tmp at a time. Whatever the reference, every
protal run that loads an index, and every build, handles its fixed-size key map (~3 GB raw): a second or more each.
So one paired-end run aligns the module's samples (baseline()), and tests that change only the profiling profile its
SAMs again (--profile_only, or the SAMs copied to where a run looks for them), which loads no index; the read types
share one map run (read_type_run()); builds compress fast (FAST_BUILD), and four builds remain.
"""

import csv
import filecmp
import functools
import glob
import gzip
import hashlib
import io
import math
import os
import random
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import threading
import time
import unittest

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
DB = os.environ.get("PROTAL_TEST_DB", "")
DB = os.path.abspath(DB) if DB else ""  # protal runs in temporary folders
PROTAL = os.path.abspath(os.environ.get("PROTAL", os.path.join(ROOT, "build", "protal")))
SIMULATE = os.path.abspath(os.environ.get("SIMULATE", os.path.join(ROOT, "build", "simulate_metagenomes")))
QCMSA = os.path.join(ROOT, "scripts", "qcmsa.py")
REQUIRED = os.environ.get("PROTAL_TESTS_REQUIRED", "") not in ("", "0")
KEEP = bool(os.environ.get("PROTAL_TEST_KEEP"))
# The test builds' compression: fast. Whatever the reference, a build compresses the index's fixed-size key map (~3 GB
# raw), which takes seconds at --build's default level 19 on 2 threads and a fraction of that at level 1 on 4.
FAST_BUILD = ("--compress_level", "1", "-t", "4")
READS = None  # directory with the simulated reads, set up once per module (require_database)
FILES = DB     # the database's separate files: DB, or DB unpacked if it is a single file
UNPACKED = None


class MissingPrerequisite(RuntimeError):
    """A prerequisite of a test is missing, and PROTAL_TESTS_REQUIRED=1 asks for every test to run."""


def unavailable(reason):
    """The exception for a missing prerequisite: a skip, or an error with PROTAL_TESTS_REQUIRED=1."""
    if REQUIRED:
        return MissingPrerequisite(f"{reason} (PROTAL_TESTS_REQUIRED=1: missing prerequisites fail)")
    return unittest.SkipTest(reason)


def require_zstd(purpose):
    if not shutil.which("zstd"):
        raise unavailable(f"the zstd CLI is needed {purpose}")


def require_binary(path, name, variable):
    if not os.access(path, os.X_OK):
        raise unavailable(f"{name} not found at {path} (set {variable})")


def scratch(prefix):
    """A temporary folder, removed when the module's tests are done (kept with PROTAL_TEST_KEEP)."""
    path = tempfile.mkdtemp(prefix=prefix)
    if not KEEP:
        unittest.addModuleCleanup(shutil.rmtree, path, ignore_errors=True)
    return path


def single_file(db):
    """The single-file database db is or holds (protal's precedence: separate files first), else None."""
    if os.path.isfile(db):
        return db
    bundle = os.path.join(db, "database.protal")
    has_index = any(os.path.exists(os.path.join(db, f)) for f in ("index.prx", "index.prx.zst"))
    return bundle if os.path.exists(bundle) and not has_index else None


def db_file(name):
    """Path of a database file as protal picks it: <name> if present, else <name>.zst."""
    raw = os.path.join(FILES, name)
    return raw if os.path.exists(raw) or not os.path.exists(raw + ".zst") else raw + ".zst"


def db_content(name):
    """Content of a database file (decompressed if the database holds <name>.zst)."""
    path = db_file(name)
    if path.endswith(".zst"):
        return subprocess.run(["zstd", "-dc", path], check=True, stdout=subprocess.PIPE).stdout
    with open(path, "rb") as fh:
        return fh.read()


def symlink_db(folder, skip=(), extra=None):
    """A database in `folder` of symlinks to the test database's files, without those in `skip` (database.protal
    always: separate files go first anyway) and with `extra` ({name: path}) linked in."""
    os.makedirs(folder, exist_ok=True)
    for f in glob.glob(os.path.join(FILES, "*")):
        name = os.path.basename(f)
        if name not in skip and name != "database.protal" and name not in (extra or {}):
            os.symlink(f, os.path.join(folder, name))
    for name, path in (extra or {}).items():
        os.symlink(path, os.path.join(folder, name))
    return folder


_database = None  # the outcome of set_up_database: None before it ran, True, or the exception it raised


def require_database():
    """The database's files and the simulated reads, set up once for the module (set_up_database); every class that
    needs them gets the same skip or error."""
    global _database
    if _database is None:
        try:
            set_up_database()
            _database = True
        except Exception as e:  # unittest.SkipTest too
            _database = e
    if _database is not True:
        raise type(_database)(*_database.args)


def set_up_database():
    global READS, FILES, UNPACKED
    if not DB or not os.path.exists(DB):
        raise unavailable("set PROTAL_TEST_DB to a protal database (just mini-db builds data/mini_db/protal_db)")
    require_binary(PROTAL, "protal", "PROTAL")
    if not shutil.which("zstd") and sys.version_info < (3, 14):
        raise unavailable("the zstd CLI (or Python 3.14) is needed to read protal's default .sam.zst files")
    bundle = single_file(DB)
    if bundle:
        UNPACKED = scratch("protal_e2e_db_")
        rc, log = run(UNPACKED, "--unpack_db", "--db", bundle, "--unpack_dir", UNPACKED, "-t", "4")
        if rc != 0:
            raise RuntimeError("protal --unpack_db failed:\n" + log[-3000:])
        FILES = UNPACKED
    index = db_file("index.prx")
    if not os.path.isfile(index) or os.path.getsize(index) == 0:
        raise unavailable(f"{DB} holds no index: set PROTAL_TEST_DB to a protal database")
    if db_file("reference.fna").endswith(".zst"):
        require_zstd("to read reference.fna.zst")
    READS = scratch("protal_e2e_reads_")
    simulate_reads("sa", pairs_per_gene=12, seed=1)
    simulate_reads("sb", pairs_per_gene=12, seed=2)
    simulate_reads("sr", pairs_per_gene=12, seed=3, random_r2=True)
    # Two species at 12 pairs per gene and one error-free pair of Mockella alpha, too few to call it: the sample in
    # which the knob decides a call. (Three pairs were, until the k-mers of 15-base cores with one value in the index
    # counted as unique, 2026-10-06: in the mini database they are 85% of the unique k-mers, and the 0.6.0a forest it
    # ships scores three perfect pairs 0.55, one pair 0.34.)
    alpha = species()["s__Mockella alpha"]
    simulate_reads("thin", pairs_per_gene=12, seed=6, taxa=set(species().values()) - {alpha},
                   extra_pairs=few_pairs(alpha, n=1))
    simulate_foreign_reads("foreign", seed=7)


def revcomp(seq):
    return seq[::-1].translate(str.maketrans("ACGTN", "TGCAN"))


def reference_genes():
    """[(name, sequence)] of the database's reference, names <taxid>_<gene id>."""
    return list(_reference_genes(FILES))


@functools.lru_cache(maxsize=None)
def _reference_genes(files):
    genes, name = [], None
    for line in db_content("reference.fna").decode().splitlines():
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
        elif line:
            genes.append((name, line.upper()))
    return tuple(genes)


def taxonomy():
    """internal_taxonomy.dmp as {taxid: row} (row: the file's columns by name)."""
    with open(db_file("internal_taxonomy.dmp")) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = [dict(zip(header, line.rstrip("\n").split("\t"))) for line in fh if line.strip()]
    return {row["id"]: row for row in rows}


def species():
    """{name: taxid} of the species that have reference genes (taxids as text)."""
    with_genes = {name.split("_")[0] for name, _ in reference_genes()}
    return {row["name"]: taxid for taxid, row in taxonomy().items() if row["rank"] == "species" and taxid in with_genes}


def few_pairs(taxid, n=3):
    """n error-free read pairs of a species, one from each of its first n genes of 300 bp or more."""
    genes = [seq for name, seq in reference_genes() if name.split("_")[0] == taxid and len(seq) >= 300][:n]
    return [(gene[:100], revcomp(gene[200:300])) for gene in genes]


def simulate_reads(prefix, pairs_per_gene, seed, random_r2=False, pairs_per_kb=None, taxa=None, extra_pairs=()):
    """Write <prefix>_R1.fq/_R2.fq in READS: 220-320 bp fragments, 100 bp reads, both orientations,
    0.5% substitutions and a few '#' (Q2) bases. With random_r2 the second mate cannot align.
    With pairs_per_kb, genes get pairs in proportion to their length (even depth), shorter genes
    included, instead of pairs_per_gene each. taxa (taxids as text) limits the genes to those species';
    extra_pairs ([(read1, read2)]) are appended as they are (quality I). Returns the number of pairs."""
    rng = random.Random(seed)

    def mutate(seq):
        bases, qual = list(seq), ["I"] * len(seq)
        for i, b in enumerate(bases):
            if rng.random() < 0.005:
                bases[i] = rng.choice([c for c in "ACGT" if c != b])
            if rng.random() < 0.01:
                qual[i] = "#"
        return "".join(bases), "".join(qual)

    n = 0
    with open(os.path.join(READS, f"{prefix}_R1.fq"), "w") as r1, open(os.path.join(READS, f"{prefix}_R2.fq"), "w") as r2:
        for name, gene in reference_genes():
            if taxa is not None and name.split("_")[0] not in taxa:
                continue
            if pairs_per_kb is None:
                if len(gene) < 320:
                    continue
                pairs, shortest, longest = pairs_per_gene, 220, 320
            else:
                pairs = round(pairs_per_kb * len(gene) / 1000)
                shortest, longest = min(220, len(gene)), min(320, len(gene))
            for _ in range(pairs):
                flen = rng.randint(shortest, longest)
                start = rng.randint(0, len(gene) - flen)
                frag = gene[start:start + flen]
                if rng.random() < 0.5:
                    frag = revcomp(frag)
                s1, q1 = mutate(frag[:100])
                s2, q2 = mutate(revcomp(frag[-100:]))
                if random_r2:
                    s2 = "".join(rng.choice("ACGT") for _ in range(100))
                n += 1
                r1.write(f"@{prefix}.{n}/1\n{s1}\n+\n{q1}\n")
                r2.write(f"@{prefix}.{n}/2\n{s2}\n+\n{q2}\n")
        for s1, s2 in extra_pairs:
            n += 1
            r1.write(f"@{prefix}.{n}/1\n{s1}\n+\n{'I' * len(s1)}\n")
            r2.write(f"@{prefix}.{n}/2\n{s2}\n+\n{'I' * len(s2)}\n")
    return n


def simulate_foreign_reads(prefix, seed, pairs=150):
    """Reads of nothing in the database: `pairs` pairs of random sequence, and as many from the reference genes with a
    quarter of their bases substituted (a genome far from every one in the database: a few seeds, no alignment)."""
    rng = random.Random(seed)
    genes = [seq for _, seq in reference_genes() if len(seq) >= 320]
    with open(os.path.join(READS, f"{prefix}_R1.fq"), "w") as r1, open(os.path.join(READS, f"{prefix}_R2.fq"), "w") as r2:
        for i in range(1, 2 * pairs + 1):
            if i <= pairs:
                frag = "".join(rng.choice("ACGT") for _ in range(300))
            else:
                gene = rng.choice(genes)
                start = rng.randint(0, len(gene) - 300)
                frag = "".join(rng.choice([c for c in "ACGT" if c != b]) if rng.random() < 0.25 else b
                               for b in gene[start:start + 300])
            r1.write(f"@{prefix}.{i}/1\n{frag[:100]}\n+\n{'I' * 100}\n")
            r2.write(f"@{prefix}.{i}/2\n{revcomp(frag[-100:])}\n+\n{'I' * 100}\n")


def head_reads(src, dst, pairs):
    """The first `pairs` records of a FASTQ file, written to dst."""
    with open(src) as fh:
        lines = fh.readlines()[:4 * pairs]
    with open(dst, "w") as fh:
        fh.writelines(lines)
    return dst


def run(cwd, *args, binary=None, timeout=900):
    """Run protal (or `binary`) in cwd; returns (exit code, combined output)."""
    proc = subprocess.run(["timeout", str(timeout), binary or PROTAL, *args], cwd=cwd,
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, errors="replace")
    return proc.returncode, proc.stdout


def python(*args):
    """Run a Python script (sys.executable); its output, or RuntimeError with it if it fails."""
    proc = subprocess.run([sys.executable, *args], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if proc.returncode != 0:
        raise RuntimeError(os.path.basename(args[0]) + " failed:\n" + proc.stdout[-3000:])
    return proc.stdout


def reads(*prefixes):
    """-1/-2/--prefix arguments for the given simulated samples."""
    return ["-1", ",".join(os.path.join(READS, f"{p}_R1.fq") for p in prefixes),
            "-2", ",".join(os.path.join(READS, f"{p}_R2.fq") for p in prefixes),
            "--prefix", ",".join(prefixes)]


def single_reads(*prefixes):
    """-1/--prefix arguments: the first mates of the given simulated samples, as single-end reads."""
    return ["-1", ",".join(os.path.join(READS, f"{p}_R1.fq") for p in prefixes), "--prefix", ",".join(prefixes)]


def profile_only(cwd, out, sams, *extra, prefixes=None, db=None, threads=2, strains=False):
    """protal --profile_only of the SAM files `sams` into `out` (--no_strains unless `strains`)."""
    args = ["--db", db or DB, "--profile_only", ",".join(sams), "-o", out, "-t", str(threads), "--no_qcmsa"]
    if prefixes:
        args += ["--prefix", ",".join(prefixes)]
    if not strains:
        args.append("--no_strains")
    return run(cwd, *args, *extra)


def sam_path(path):
    """The SAM protal wrote as `path` (a name ending in .sam), in whichever format: path, or path with
    .zst (the default for names protal picks) or .gz appended."""
    return next((p for p in (path, path + ".zst", path + ".gz") if os.path.exists(p)), path)


def find_sams(pattern):
    """glob for SAMs in any format: pattern (ending in .sam), and with .zst or .gz appended."""
    return sorted(glob.glob(pattern) + glob.glob(pattern + ".zst") + glob.glob(pattern + ".gz"))


def sam_text(path):
    """The text of a SAM (plain, .gz or .zst; see sam_path)."""
    path = sam_path(path)
    if path.endswith(".gz"):
        with gzip.open(path, "rt") as fh:
            return fh.read()
    if path.endswith(".zst"):
        try:
            from compression import zstd  # Python 3.14+
            with zstd.open(path, "rt") as fh:
                return fh.read()
        except ImportError:
            return subprocess.run(["zstd", "-dcq", path], check=True, stdout=subprocess.PIPE, text=True).stdout
    with open(path) as fh:
        return fh.read()


def open_sam(path):
    return io.StringIO(sam_text(path))


def sam_records(path):
    return [line.split("\t") for line in sam_text(path).splitlines() if not line.startswith("@")]


def read_text(path):
    with open(path) as fh:
        return fh.read()


def read_table(path):
    """A tab-separated file with a header line, as (header, rows)."""
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        return header, [line.rstrip("\n").split("\t") for line in fh]


def read_dicts(path):
    """A tab-separated file with a header line, as a list of dicts."""
    header, rows = read_table(path)
    return [dict(zip(header, row)) for row in rows]


def profile_species(text):
    """The species names a .profile reports (rep genome, lineage, abundance per line)."""
    return sorted(line.split("\t")[1].split(";")[-1] for line in text.splitlines() if line.strip())


def digest(path):
    sha = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 24), b""):
            sha.update(block)
    return sha.hexdigest()


def write_tiny_db(db, genes, taxonomy_text):
    """A database folder of genes ({(taxid, gene id): sequence}) for --build: reference.fna, reference.map and
    internal_taxonomy.dmp."""
    os.makedirs(db, exist_ok=True)
    with open(os.path.join(db, "reference.fna"), "w") as fna, open(os.path.join(db, "reference.map"), "w") as mp:
        offset = 0
        for (taxid, gene), seq in genes.items():
            header = f">{taxid}_{gene}\n"
            fna.write(header + seq + "\n")
            mp.write(f"{taxid}\t{gene}\t{offset + len(header)}\t{offset + len(header) + len(seq)}\n")
            offset += len(header) + len(seq) + 1
    with open(os.path.join(db, "internal_taxonomy.dmp"), "w") as fh:
        fh.write(taxonomy_text)
    return db


def remove_indexes(db):
    """Remove a test build's index (~3 GB raw whatever the reference: a fixed-size key map)."""
    for index in ("index.prx", "index.prx.zst"):
        if os.path.exists(os.path.join(db, index)):
            os.remove(os.path.join(db, index))


class Baseline:
    """The module's paired-end run: protal with its defaults (-1/-2/-o, strain MSAs; --no_qcmsa, and a truth file per
    sample) over sa and sb (12 pairs on every gene of 320 bp or more), sr (sa's read1 mates with random read2 mates),
    thin (Mockella beta and Fakibacter gamma, and one pair of Mockella alpha) and foreign (reads of nothing in the
    database). Tests that change only the profiling profile its SAMs again."""

    SAMPLES = ("sa", "sb", "sr", "thin", "foreign")

    def __init__(self):
        self.work = scratch("protal_e2e_baseline_")
        self.out = os.path.join(self.work, "out")
        taxa = species()
        genus = next(t for t, row in taxonomy().items() if row["rank"] == "genus")
        truth = {"sa": list(taxa.values()), "sb": list(taxa.values()), "sr": list(taxa.values()),
                 "thin": list(taxa.values()), "foreign": [genus]}  # foreign: a genus, in the taxonomy but no species
        files = []
        for sample in self.SAMPLES:
            files.append(os.path.join(self.work, f"{sample}.truth.tsv"))
            with open(files[-1], "w") as fh:
                fh.write("".join(f"{t}\n" for t in truth[sample]))
        self.rc, self.log = run(self.work, "--db", DB, *reads(*self.SAMPLES), "-o", "out", "-t", "4", "--no_qcmsa",
                                "--profile_truth", ",".join(files))
        if self.rc != 0:
            raise RuntimeError(f"the baseline run exited {self.rc}:\n" + self.log[-3000:])

    def path(self, *parts):
        return os.path.join(self.out, *parts)

    def sam(self, sample):
        return sam_path(self.path(f"{sample}.sam"))

    def text(self, name):
        return read_text(self.path(name))

    def profiles(self, sample):
        """{file name: text} of a sample's profile files (.profile, .profile.log, .profile.gene.log, .profile.genes.log)."""
        return {name: self.text(name) for name in (f"{sample}.profile", f"{sample}.profile.log",
                                                    f"{sample}.profile.gene.log", f"{sample}.profile.genes.log")}

    def place_sams(self, folder, samples, names=None):
        """Copy the samples' SAMs into folder (as names[i].sam.zst), where a run with -o folder (or a map's SAM folder)
        takes them instead of aligning the reads."""
        os.makedirs(folder, exist_ok=True)
        for sample, name in zip(samples, names or samples):
            shutil.copy(self.sam(sample), os.path.join(folder, f"{name}.sam.zst"))
        return folder


_baseline = None


def baseline():
    """The module's Baseline, run once on first use."""
    global _baseline
    if _baseline is None:
        require_database()
        try:
            _baseline = Baseline()
        except Exception as e:
            _baseline = e
    if isinstance(_baseline, Exception):
        raise RuntimeError(f"no baseline run: {_baseline}")
    return _baseline


class WorkDir(unittest.TestCase):
    """A test class with its own scratch directory, removed after the class (unless PROTAL_TEST_KEEP is set), also when
    its setUpClass fails."""

    @classmethod
    def setUpClass(cls):
        cls.work = tempfile.mkdtemp(prefix=f"protal_e2e_{cls.__name__}_")
        if not KEEP:
            cls.addClassCleanup(shutil.rmtree, cls.work, ignore_errors=True)

    def path(self, *parts):
        return os.path.join(self.work, *parts)


class ProtalTest(WorkDir):
    """A test class that runs protal without the test database."""

    @classmethod
    def setUpClass(cls):
        require_binary(PROTAL, "protal", "PROTAL")
        super().setUpClass()


class DbTest(WorkDir):
    """A test class that needs the test database and the simulated reads."""

    @classmethod
    def setUpClass(cls):
        require_database()
        super().setUpClass()


class CompleteRunTest(DbTest):
    """The baseline run over five samples (two normal ones, sr whose read2 mates cannot align, thin, foreign): its
    files and records."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.base = baseline()

    def test_exit_code(self):
        self.assertEqual(self.base.rc, 0, self.base.log[-3000:])

    def test_outputs(self):
        base = self.base
        for sample in base.SAMPLES:
            self.assertTrue(os.path.isfile(base.path(f"{sample}.sam.zst")), "SAMs protal names are zstd-compressed")
            self.assertTrue(os.path.isfile(base.path(f"{sample}.profile")), sample)
            self.assertTrue(os.path.isfile(base.path("misc", f"{sample}_runtime.tsv")), sample)
        self.assertTrue(sam_records(base.sam("sa")))
        self.assertEqual(glob.glob(base.path("**", "*.partial"), recursive=True), [], "no .partial SAM left")
        # strains/ and misc/ are made in -1/-2/-o mode too (they once were not). Every species is reported by at least
        # two samples and gets an MSA, with its partition file.
        msas = sorted(os.path.basename(f) for f in glob.glob(base.path("strains", "*.raw.msa.fna")))
        self.assertEqual(msas, sorted(name.replace(" ", "_") + ".raw.msa.fna" for name in species()))
        for msa in msas:
            self.assertTrue(os.path.isfile(base.path("strains", msa.replace(".raw.msa.fna", ".raw.partition.txt"))), msa)

    def test_read_names_and_pairs(self):
        records = sam_records(self.base.sam("sa"))
        self.assertTrue(records)
        bad = [r[0] for r in records if not (r[0].startswith("sa.") and r[0][3:].isdigit())]
        self.assertEqual(bad[:5], [], "QNAME is the read id without its /1 /2 suffix, nothing more")
        self.assertTrue(any(int(r[1]) & 0x2 for r in records), "proper pairs are flagged")

    def test_read1_only_pairs_are_written(self):
        records = [r for r in sam_records(self.base.sam("sr")) if not int(r[1]) & 0x100]
        read1_only = [r for r in records if int(r[1]) & 0x40 and int(r[1]) & 0x8]
        with open(os.path.join(READS, "sr_R1.fq")) as fh:
            total = sum(1 for _ in fh) // 4
        self.assertGreater(len(read1_only), total // 2, "primary read1 records with mate unmapped (0x8)")

    def test_reference_calls_are_not_blanked(self):
        bases = ns = 0
        for msa in glob.glob(self.base.path("strains", "*.raw.msa.fna")):
            with open(msa) as fh:
                for line in fh:
                    if not line.startswith(">"):
                        seq = line.strip()
                        bases += len(seq)
                        ns += seq.count("N")
        self.assertGreater(bases, 0)
        self.assertLess(ns / bases, 0.005, f"N fraction {100 * ns / bases:.3f}% with error-only reads")

        stats = glob.glob(self.base.path("strains", "*.snp_stats.tsv"))
        self.assertTrue(stats)
        retained = 0
        for path in stats:
            retained += sum(int(row["refs_retained"]) for row in read_dicts(path))
        self.assertGreater(retained, 0, "reference calls retained at variant positions")

    def test_partitions_are_one_based_and_cover_the_msa(self):
        parts = glob.glob(self.base.path("strains", "*.raw.partition.txt"))
        self.assertEqual(len(parts), len(species()), "one partition file per species' MSA")
        for part in parts:
            with open(part) as fh:
                ranges = [tuple(int(x) for x in line.split("=")[1].split("-")) for line in fh if line.strip()]
            with open(part.replace(".raw.partition.txt", ".raw.msa.fna")) as fh:
                length = len([line for line in fh if not line.startswith(">")][0].strip())
            self.assertTrue(ranges, part)
            self.assertEqual(ranges[0][0], 1, part)
            self.assertEqual(ranges[-1][1], length, part)
            for (_, end), (start, _) in zip(ranges, ranges[1:]):
                self.assertEqual(start, end + 1, part)


class AccuracyTest(DbTest):
    """Is the baseline's profile right? Its truth counts, abundances against the pairs simulated per species, and reads of
    nothing in the database. examples/mini_db/run.sh is the deeper check: reads of genomes, strains other than the
    references, abundances by cell."""

    # The largest difference allowed between a species' abundance and the one its simulated pairs give.
    ABUNDANCE_TOLERANCE = 0.01

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.base = baseline()

    def test_truth_counts(self):
        # sa, sb and sr hold the three species, thin two of them and one pair of the third (too few to call it),
        # foreign none: its truth names a genus, a taxon the database has no genes of.
        expected = {"sa": "TP 3, FP 0, FN 0 (and 0", "sb": "TP 3, FP 0, FN 0 (and 0", "sr": "TP 3, FP 0, FN 0 (and 0",
                    "thin": "TP 2, FP 0, FN 1 (and 0", "foreign": "TP 0, FP 0, FN 0 (and 1"}
        for sample, counts in expected.items():
            self.assertIn(f"Sample {sample}: {counts} true species not in the database)", self.base.log)

    def test_reads_of_nothing_in_the_database_give_an_empty_profile(self):
        self.assertEqual(self.base.rc, 0)
        self.assertEqual(self.base.text("foreign.profile"), "")
        self.assertIn("No taxon passes the model in sample foreign", self.base.log)
        _, rows = read_table(self.base.path("foreign.profile.log"))
        self.assertEqual([row for row in rows if row[0] != "0"], [], "no taxon called")
        self.assertEqual([r for r in sam_records(self.base.sam("foreign")) if not int(r[1]) & 0x4 and int(r[4]) >= 4], [],
                         "no read aligns with a MAPQ the profiler takes")

    def expected_abundances(self):
        """{species name: relative abundance} that sa's and sb's pairs give: 12 pairs on each gene of 320 bp or more,
        each pair 200 bases of the gene (the mates never overlap), so a gene's depth is 2400 / its length. A species'
        depth is the median over the genes with reads (all hit genes and enough depth: BlendedDepth's weight is 1);
        abundances are the depths' shares."""
        lengths = {}
        for name, seq in reference_genes():
            if len(seq) >= 320:
                lengths.setdefault(name.split("_")[0], []).append(len(seq))
        depths = {}
        for taxid, values in lengths.items():
            per_gene = sorted(2400 / length for length in values)
            mid = len(per_gene) // 2
            depths[taxid] = per_gene[mid] if len(per_gene) % 2 else (per_gene[mid - 1] + per_gene[mid]) / 2
        total = sum(depths.values())
        return {name: depths[taxid] / total for name, taxid in species().items()}

    def test_profile_format_and_abundances(self):
        tax = taxonomy()
        expected = self.expected_abundances()
        for sample in ("sa", "sb"):
            rows = [line.split("\t") for line in self.base.text(f"{sample}.profile").splitlines()]
            self.assertEqual(len(rows), 3, sample)
            found = {}
            for rep_genome, lineage, abundance in rows:  # three fields per line
                name = lineage.split(";")[-1]
                self.assertTrue(lineage.startswith("d__"), lineage)
                self.assertEqual(rep_genome, tax[species()[name]]["rep_genome"], lineage)
                found[name] = float(abundance)
            self.assertAlmostEqual(sum(found.values()), 1, delta=1e-5)
            for name, share in expected.items():
                self.assertLess(abs(found[name] - share), self.ABUNDANCE_TOLERANCE,
                                f"{sample}, {name}: abundance {found[name]:.4f}, its pairs give {share:.4f}")


class OutputFilesTest(DbTest):
    """One sample from a map, with a sample ID and a truth file in the map, and --taxon_statistics: what the profile's
    companion files hold. The SAM is the baseline's (a map without a SAM column looks for <OUTPUT_DIR>/<PREFIX>.sam.zst):
    no read is aligned again."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        base = baseline()
        truth = os.path.join(cls.work, "truth.tsv")
        taxa = species()
        genus = next(t for t, row in taxonomy().items() if row["rank"] == "genus")
        with open(truth, "w") as fh:
            fh.write("".join(f"{t}\n" for t in taxa.values()) + f"{genus}\n")  # a genus: in the taxonomy, not in the database
        cls.sample_map = os.path.join(cls.work, "samples.map")
        with open(cls.sample_map, "w") as fh:
            fh.write(f"#OUTPUT_DIR\t{os.path.join(cls.work, 'out')}\n#INPUT_DIR\t{READS}\n")
            fh.write("#SAMPLEID\tPREFIX\tFIRST\tSECOND\tPROFILE_TRUTH\n")
            fh.write(f"sample_a\tpa\tsa_R1.fq\tsa_R2.fq\t{truth}\n")
        base.place_sams(os.path.join(cls.work, "out"), ["sa"], ["pa"])
        cls.rc, cls.log = run(cls.work, "--db", DB, "--map", cls.sample_map, "-t", "2", "--no_qcmsa", "--taxon_statistics")

    def test_exit_code(self):
        self.assertEqual(self.rc, 0, self.log[-3000:])
        self.assertIn("All alignments are present", self.log)

    def test_the_sample_id_names_the_sample(self):
        header, rows = read_table(self.path("out", "pa.profile.gene.log"))
        self.assertEqual(header[0], "Sample")
        self.assertTrue(rows)
        self.assertEqual({row[0] for row in rows}, {"sample_a"})
        self.assertGreater(max(int(row[header.index("CoverageSum")]) for row in rows), 0)

    def test_truth_counts_name_the_sample(self):
        self.assertIn("Sample sample_a: TP 3, FP 0, FN 0 (and 1 true species not in the database)", self.log)
        self.assertNotIn("Truth: 0", self.log)
        self.assertEqual(read_text(self.path("out", "pa.profile")), baseline().text("sa.profile"))

    def test_statistics_for_a_single_sample(self):
        # The per-taxon files are written only with --taxon_statistics (since 0.7.6), then with one sample, too.
        self.assertTrue(os.path.isdir(baseline().path("misc")))
        self.assertEqual(glob.glob(baseline().path("misc", "*.statistics.tsv")), [], "not written by default")
        stats = glob.glob(self.path("out", "misc", "*.statistics.tsv"))
        self.assertEqual(len(stats), len(species()), "one per taxon with reads")
        _, rows = read_table(stats[0])
        self.assertEqual([row[0] for row in rows], ["sample_a"])

    def test_profile_log_columns(self):
        header, rows = read_table(self.path("out", "pa.profile.log"))
        self.assertEqual(header[:2], ["Predicted", "Probability"])
        self.assertTrue(rows)
        for row in rows:
            self.assertEqual(len(row), len(header))
            self.assertGreaterEqual(float(row[header.index("VCov")]), 0)
            self.assertNotIn("VCOV: -1", row[header.index("Summary")])
            if row[0] == "0":
                self.assertEqual(float(row[header.index("Abundance")]), 0)

    def test_genes_log_abundances(self):
        header, rows = read_table(self.path("out", "pa.profile.genes.log"))
        self.assertTrue(rows)
        called = 0
        for row in rows:
            value = float(row[header.index("TaxAbundance")])
            self.assertTrue(0 <= value <= 1, row)
            if row[0] == "0":
                self.assertEqual(value, 0, "a rejected taxon has no abundance")
            else:
                called += 1
        self.assertGreater(called, 0)
        self.assertTrue(any(int(row[header.index("UniqueMerReads")]) > 0 for row in rows))
        self.assertGreater(max(int(row[header.index("TotalReads")]) for row in rows), 10)
        self.assertTrue(all(float(row[header.index("MAPQ")]) <= 255 for row in rows), "MAPQ is a mean, not a sum")


class MateAssignmentTest(DbTest):
    """Fragments whose best alignment is mate 2's alone, and fragments over two genes, reach the SAM (one run, two
    samples). The run has no -o: its outputs go to the folder it runs in (such runs once aborted on making the folder
    "")."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        rng = random.Random(21)
        mate2 = []
        for _, gene in reference_genes():
            if len(gene) >= 320:
                start = rng.randint(0, len(gene) - 100)
                mate2.append(("".join(rng.choice("ACGT") for _ in range(100)), revcomp(gene[start:start + 100])))
        genes = {}
        for name, seq in reference_genes():
            taxid, gene = (int(x) for x in name.split("_")[:2])
            genes[(taxid, gene)] = seq
        junction, cls.expected = [], []
        for (taxid, gene), seq in sorted(genes.items()):
            nxt = genes.get((taxid, gene + 1))
            if nxt is not None and len(seq) >= 100 and len(nxt) >= 100:
                junction.append((seq[-100:], revcomp(nxt[:100])))
                cls.expected.append((f"{taxid}_{gene}", f"{taxid}_{gene + 1}"))
        cls.pairs = {"mate2": mate2, "junction": junction}
        for name, pairs in cls.pairs.items():
            with open(os.path.join(cls.work, f"{name}_R1.fq"), "w") as r1, open(os.path.join(cls.work, f"{name}_R2.fq"), "w") as r2:
                for i, (s1, s2) in enumerate(pairs, 1):
                    r1.write(f"@{name}.{i}/1\n{s1}\n+\n{'I' * len(s1)}\n")
                    r2.write(f"@{name}.{i}/2\n{s2}\n+\n{'I' * len(s2)}\n")
        os.makedirs(os.path.join(cls.work, "out"))
        cls.rc, cls.log = run(os.path.join(cls.work, "out"), "--db", DB,
                              "-1", ",".join(os.path.join(cls.work, f"{n}_R1.fq") for n in cls.pairs),
                              "-2", ",".join(os.path.join(cls.work, f"{n}_R2.fq") for n in cls.pairs),
                              "--prefix", ",".join(cls.pairs), "-t", "2", "--no_qcmsa", "--no_strains")

    def test_outputs_in_the_current_folder(self):
        self.assertEqual(self.rc, 0, self.log[-3000:])
        for name in self.pairs:
            self.assertTrue(find_sams(self.path("out", f"{name}.sam")), name)
            self.assertTrue(os.path.isfile(self.path("out", f"{name}.profile")), name)
        self.assertTrue(os.path.isdir(self.path("out", "misc")))

    def primary(self, name):
        self.assertEqual(self.rc, 0, self.log[-3000:])
        return [r for r in sam_records(find_sams(self.path("out", f"{name}.sam"))[0]) if not int(r[1]) & 0x100]

    def test_pairs_where_only_mate2_aligns(self):
        read2_only = [r for r in self.primary("mate2") if int(r[1]) & 0x80 and int(r[1]) & 0x8]
        pairs = len(self.pairs["mate2"])
        self.assertGreater(len(read2_only), 0.9 * pairs, f"{len(read2_only)} of {pairs} written")

    def test_pairs_over_two_genes(self):
        by_read = {}
        for r in self.primary("junction"):
            by_read.setdefault(r[0], []).append(r)
        good = 0
        for i, (gene_a, gene_b) in enumerate(self.expected, 1):
            recs = by_read.get(f"junction.{i}", [])
            mates = {("1" if int(r[1]) & 0x40 else "2"): r for r in recs}
            if (len(recs) == 2 and mates.get("1", [None, None, None])[2] == gene_a and
                    mates.get("2", [None, None, None])[2] == gene_b and all(int(r[4]) >= 4 for r in recs)):
                good += 1
        pairs = len(self.pairs["junction"])
        self.assertGreater(good, 0.9 * pairs, f"{good} of {pairs} fragments written with both mates")


class MsaSampleSelectionTest(DbTest):
    """A species' MSA takes only the samples in which the model accepts the species: thin's one pair of Mockella alpha
    does not call it, so thin is left out of its MSA (with --msa_min_hcov 0 coverage would not keep it out) and is in the
    MSAs of the two species it reports."""

    def test_rejected_samples_are_left_out(self):
        samples = ["sa", "sb", "thin"]
        baseline().place_sams(self.path("out"), samples)
        rc, log = run(self.work, "--db", DB, *reads(*samples), "-o", "out", "-t", "2", "--no_qcmsa", "--msa_min_hcov", "0")
        self.assertEqual(rc, 0, log[-3000:])
        self.assertIn("All alignments are present", log)
        self.assertNotIn("Mockella alpha", read_text(self.path("out", "thin.profile")), "one pair does not call the species")
        for name in species():
            with open(self.path("out", "strains", name.replace(" ", "_") + ".raw.msa.fna")) as fh:
                names = [line[1:].strip() for line in fh if line.startswith(">")]
            self.assertIn("sa", names, name)
            self.assertIn("sb", names, name)
            if name == "s__Mockella alpha":
                self.assertNotIn("thin", names)
            else:
                self.assertIn("thin", names, name)


class StrainEdgeCaseTest(DbTest):
    def test_species_without_msa_genes(self):
        # No position reaches --msa_min_depth, so no species has MSA columns (this used to segfault).
        # Two samples: MSAs are built only across samples.
        baseline().place_sams(self.path("out"), ["sa", "sb"])
        rc, log = run(self.work, "--db", DB, *reads("sa", "sb"), "-o", "out", "-t", "2", "--no_qcmsa",
                      "--msa_min_depth", "100000", "--msa_min_hcov", "0")
        self.assertEqual(rc, 0, log[-3000:])
        self.assertIn("has enough coverage for an MSA", log)
        self.assertEqual(glob.glob(self.path("out", "strains", "*.raw.msa.fna")), [])


class MSAKnobTest(DbTest):
    """A species' MSA holds the samples whose profile reports it; --msa_knob sets another threshold. Each test profiles
    the baseline's SAMs of sa and sb into an output folder of its own."""

    def strains_run(self, out, *extra):
        baseline().place_sams(self.path(out), ["sa", "sb"])
        rc, log = run(self.work, "--db", DB, *reads("sa", "sb"), "-o", out, "-t", "2", "--no_qcmsa", "--msa_min_hcov", "0",
                      *extra)
        self.assertEqual(rc, 0, log[-3000:])
        return log

    def calls(self, out, sample):
        """{species: (reported, probability)} of a sample's profile, species spelled as in species.tsv."""
        rows = read_dicts(self.path(out, f"{sample}.profile.log"))
        return {r["Name"].replace(" ", "_"): (r["Predicted"] == "1", float(r["Probability"])) for r in rows}

    def species_list(self, out):
        return {r["species"]: int(r["samples"]) for r in read_dicts(self.path(out, "strains", "species.tsv"))}

    def test_msa_samples_mirror_the_profiles(self):
        self.strains_run("out")
        calls = {s: self.calls("out", s) for s in ("sa", "sb")}
        listed = self.species_list("out")
        both = [sp for sp in calls["sa"] if calls["sa"][sp][0] and calls["sb"].get(sp, (False, 0))[0]]
        self.assertTrue(both, "species reported by both samples")
        for species_name in set(calls["sa"]) | set(calls["sb"]):
            reported = sum(calls[s].get(species_name, (False, 0))[0] for s in calls)
            if reported >= 2:
                self.assertEqual(listed.get(species_name), reported, species_name)
            else:
                self.assertNotIn(species_name, listed, "an MSA needs 2 samples that report the species")

    def test_msa_knob_admits_unreported_species(self):
        # No species scores --knob 1, yet --msa_knob 0 builds the MSAs of all species with reads in both
        # samples. At 2.4x on every gene, the species' reads are strong evidence: those the profiles
        # leave out are listed in unreported_species.tsv.
        log = self.strains_run("out_knob", "--knob", "1", "--msa_knob", "0")
        calls = {s: self.calls("out_knob", s) for s in ("sa", "sb")}
        listed = self.species_list("out_knob")
        self.assertTrue(set(calls["sa"]) & set(calls["sb"]))
        for species_name in set(calls["sa"]) & set(calls["sb"]):
            self.assertEqual(listed.get(species_name), 2, species_name)
        rows = read_dicts(self.path("out_knob", "misc", "unreported_species.tsv"))
        self.assertTrue(rows, log[-3000:])
        for r in rows:
            reported, probability = calls[r["sample"]][r["species"].replace(" ", "_")]
            self.assertFalse(reported, r)
            self.assertLess(probability, 1)
            self.assertEqual(r["passes_msa_knob"], "yes")
        self.assertIn("are not reported, although their own reads are strong evidence", log)


class CombiningRunsTest(DbTest):
    """Strain MSAs over the SAMs of several runs (--profile_only, which loads no index): a pattern that protal expands,
    the arguments an unquoted one becomes, and SAMs of the same name, named by their folders, all give the MSAs of the
    same SAMs listed with --prefix; without -o the strain outputs go to the folder protal runs in. The baseline's sa and
    sb are study1's samples, its thin is study2's sa."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        base = baseline()
        cls.sams = {}
        for study, sample, source in (("study1", "sa", "sa"), ("study1", "sb", "sb"), ("study2", "sa", "thin")):
            folder = os.path.join(cls.work, study, "alignments")
            os.makedirs(folder, exist_ok=True)
            cls.sams[f"{study}_{sample}"] = shutil.copy(base.sam(source), os.path.join(folder, f"{sample}.sam.zst"))
        cls.rc, cls.log = profile_only(cls.work, os.path.join(cls.work, "listed"), list(cls.sams.values()),
                                       prefixes=list(cls.sams), strains=True)

    def strain_files(self, folder):
        """{file name: text} of the raw MSAs, partitions and meta tables in folder/strains."""
        names = [os.path.basename(f) for pattern in ("*.raw.msa.fna", "*.raw.partition.txt", "*.meta.tsv")
                 for f in glob.glob(self.path(folder, "strains", pattern))]
        return {name: read_text(self.path(folder, "strains", name)) for name in names}

    def test_the_listed_samples(self):
        self.assertEqual(self.rc, 0, self.log[-3000:])
        files = self.strain_files("listed")
        self.assertTrue(any(name.endswith(".raw.msa.fna") for name in files))
        for name, text in files.items():
            if name.endswith(".raw.msa.fna"):
                rows = [line[1:].strip() for line in text.splitlines() if line.startswith(">")]
                self.assertIn("study1_sa", rows, name)
                self.assertIn("study1_sb", rows, name)

    def test_a_pattern_without_outdir(self):
        # Two SAMs named sa: every sample is named by its study folder. Without -o the profiles go next to the SAMs and
        # the strain outputs into the folder protal runs in.
        os.makedirs(self.path("here"))
        rc, log = run(self.path("here"), "--db", DB, "--profile_only", os.path.join(self.work, "study*", "alignments", "*.sam.zst"),
                      "-t", "2", "--no_qcmsa")
        self.assertEqual(rc, 0, log[-3000:])
        self.assertIn("matches 3 SAM file(s)", log)
        self.assertIn("named by the folders in which their paths differ: study1_sa", log)
        self.assertEqual(self.strain_files("here"), self.strain_files("listed"))
        self.assertTrue(os.path.isfile(self.path("study2", "alignments", "sa.profile")))

    def test_an_unquoted_pattern(self):
        # The shell expands it: the option takes the first SAM, the others are arguments (once ignored without a word).
        rc, log = run(self.work, "--db", DB, "--profile_only", *self.sams.values(), "-o", "unquoted", "-t", "2", "--no_qcmsa")
        self.assertEqual(rc, 0, log[-3000:])
        self.assertEqual(self.strain_files("unquoted"), self.strain_files("listed"))

    def test_a_folder_per_sample(self):
        for name in ("sa", "sb"):
            os.makedirs(self.path("per_sample", name))
            shutil.copy(self.sams[f"study1_{name}"], self.path("per_sample", name, "aln.sam.zst"))
        rc, log = run(self.work, "--db", DB, "--profile_only", self.path("per_sample", "*", "aln.sam.zst"), "-o", "per_sample_out",
                      "-t", "2", "--no_qcmsa")
        self.assertEqual(rc, 0, log[-3000:])
        msas = glob.glob(self.path("per_sample_out", "strains", "*.raw.msa.fna"))
        self.assertTrue(msas, log[-3000:])
        for msa in msas:
            rows = [line[1:].strip() for line in read_text(msa).splitlines() if line.startswith(">")]
            self.assertEqual(rows[1:], ["sa", "sb"], msa)

    def test_spilled_evidence(self):
        # --strain_spill: the samples' strain evidence goes to files, read back per species: the same strain outputs, and
        # the files are gone at the end. A folder protal cannot write to stops it before it starts.
        rc, log = run(self.work, "--db", DB, "--profile_only", ",".join(self.sams.values()), "--prefix", ",".join(self.sams),
                      "-o", "spilled", "-t", "2", "--no_qcmsa", "--strain_spill", self.path("spill"))
        self.assertEqual(rc, 0, log[-3000:])
        self.assertEqual(self.strain_files("spilled"), self.strain_files("listed"))
        self.assertEqual(os.listdir(self.path("spill")), [])
        with open(self.path("a_file"), "w") as fh:
            fh.write("not a folder\n")
        rc, log = run(self.work, "--db", DB, "--profile_only", ",".join(self.sams.values()), "--prefix", ",".join(self.sams),
                      "-o", "spill_refused", "-t", "2", "--no_qcmsa", "--strain_spill", self.path("a_file"))
        self.assertNotIn(rc, (0, 1), log[-3000:])
        self.assertIn("is not a folder protal can write to", log)

    def test_earlier_results_stop_protal_unless_forced(self):
        # Profiling into a folder with an earlier run's results stops before it starts; --force writes them again,
        # replacing the profiling rows of misc/<sample>_runtime.tsv (a rerun once added a second set).
        os.makedirs(self.path("again"))
        args = ["--db", DB, "--profile_only", ",".join(self.sams.values()), "--prefix", ",".join(self.sams),
                "-o", "again", "-t", "2", "--no_qcmsa"]
        rc, log = run(self.work, *args)
        self.assertEqual(rc, 0, log[-3000:])
        _, runtime = read_table(self.path("again", "misc", "study1_sa_runtime.tsv"))
        self.assertTrue(runtime)
        rc, log = run(self.work, *args)
        self.assertNotIn(rc, (0, 1), log[-3000:])
        self.assertIn("--profile_only would overwrite the results of an earlier run", log)
        rc, log = run(self.work, *args, "--force")
        self.assertEqual(rc, 0, log[-3000:])
        _, rerun = read_table(self.path("again", "misc", "study1_sa_runtime.tsv"))
        self.assertEqual([row[0] for row in rerun], [row[0] for row in runtime])

    def test_patterns_without_sams_stop_protal(self):
        rc, log = run(self.work, "--db", DB, "--profile_only", "nothing/*.sam.zst,study1/alignments,study1/alignments/*.err",
                      "-o", "nothing", "-t", "2", "--no_qcmsa")
        self.assertNotIn(rc, (0, 1), log[-3000:])
        self.assertIn("nothing/*.sam.zst matches no SAM file", log)
        self.assertIn("study1/alignments is a folder", log)
        self.assertIn("study1/alignments/*.err matches no SAM file", log)
        self.assertFalse(os.path.exists(self.path("nothing", "strains", "species.tsv")))

    def test_a_stray_argument_stops_an_alignment_run(self):
        # -1 *_R1.fq unquoted would align the first file only.
        rc, log = run(self.work, "--db", DB, "-1", os.path.join(READS, "sa_R1.fq"), os.path.join(READS, "sb_R1.fq"),
                      "--prefix", "sa", "-o", "stray")
        self.assertEqual(rc, 2, log[-3000:])
        self.assertIn("unexpected argument(s): " + os.path.join(READS, "sb_R1.fq"), log)


class LowCoverageAbundanceTest(DbTest):
    """The depth estimate, and so relative abundances, stays proportional at low coverage."""

    @staticmethod
    def taxon_depths(genes_log):
        return {row["TaxID"]: float(row["TaxVCOV"]) for row in read_dicts(genes_log)}

    def test_subsampled_depths_scale_with_the_read_count(self):
        # About 3x on every gene; the subsamples have about 0.12x and 0.36x.
        simulate_reads("even", 0, seed=4, pairs_per_kb=15)
        rng = random.Random(5)
        with open(os.path.join(READS, "even_R1.fq")) as f1, open(os.path.join(READS, "even_R2.fq")) as f2:
            r1, r2 = f1.readlines(), f2.readlines()
        fractions = {"low": 0.04, "mid": 0.12}
        for name, fraction in fractions.items():
            keep = [i for i in range(len(r1) // 4) if rng.random() < fraction]
            for lines, mate in ((r1, 1), (r2, 2)):
                with open(self.path(f"{name}_R{mate}.fq"), "w") as fh:
                    for i in keep:
                        fh.writelines(lines[4 * i:4 * i + 4])
            fractions[name] = len(keep) / (len(r1) // 4)
        rc, log = run(self.work, "--db", DB,
                      "-1", ",".join([os.path.join(READS, "even_R1.fq")] + [self.path(f"{n}_R1.fq") for n in fractions]),
                      "-2", ",".join([os.path.join(READS, "even_R2.fq")] + [self.path(f"{n}_R2.fq") for n in fractions]),
                      "--prefix", "full," + ",".join(fractions), "-o", "out", "-t", "3", "--no_strains")
        self.assertEqual(rc, 0, log[-3000:])
        full = self.taxon_depths(self.path("out", "full.profile.genes.log"))
        self.assertEqual(sorted(full), sorted(species().values()))
        for name, fraction in fractions.items():
            depths = self.taxon_depths(self.path("out", f"{name}.profile.genes.log"))
            self.assertEqual(sorted(depths), sorted(full), f"{name}: every species has reads")
            for taxid, depth in depths.items():
                ratio = depth / (fraction * full[taxid])
                self.assertLess(abs(ratio - 1), 0.3, f"{name} ({fraction:.3f} of the reads), taxon {taxid}: "
                                                     f"depth {depth:.4f} vs {fraction * full[taxid]:.4f} expected")


class GeneConservationTest(DbTest):
    """--build estimates how fast each gene diverges within species from the other genomes' copies in
    --full_reference (gene_conservation.tsv) and stores it in the database; queries take it for the model's
    conservation features, scale the depth identity margin by it with --gene_conservation db, and keep the same margin
    on every gene by default."""

    def build(self, name, rates):
        """A database of 4 species of one genus with 12 genes of 600 bp each, each species 5% times rates[i] from
        their ancestor at gene i; in full_reference, 2 other genomes per species whose gene i differs from the
        representative's at 2% times rates[i]. (Without other copies there are no factors: BuildIndexTest.)"""
        rng = random.Random(21)

        def mutate(seq, rate):
            return "".join(rng.choice([c for c in "ACGT" if c != b]) if rng.random() < rate else b for b in seq)

        ancestor = {gene: "".join(rng.choice("ACGT") for _ in range(600)) for gene in range(1, len(rates) + 1)}
        genes = {(taxid, gene): mutate(ancestor[gene], 0.05 * rates[gene - 1])
                 for taxid in range(1, 5) for gene in range(1, len(rates) + 1)}
        db = write_tiny_db(self.path(name), genes,
                           "id\tparent_id\texternal_id\tname\trank\tlevel\trep_genome\n5\t5\t0\troot\tno rank\t0\t\n" +
                           "6\t5\t0\tg__Genus\tgenus\t6\t\n" +
                           "".join(f"{t}\t6\t0\ts__Genus species{t}\tspecies\t7\tGCF_{t}\n" for t in range(1, 5)))
        full = os.path.join(db, "full_reference.fna")
        with open(full, "w") as fh:
            for (taxid, gene), seq in genes.items():
                fh.write(f">{taxid}_{gene}\n{seq}\n")
                for _ in range(2):
                    fh.write(f">{taxid}_{gene}\n{mutate(seq, 0.02 * rates[gene - 1])}\n")
        rc, log = run(self.work, "--build", "--no_bundle", "--no_profile", *FAST_BUILD, "--db", db,
                      "--reference", os.path.join(db, "reference.fna"), "--full_reference", full)
        self.assertEqual(rc, 0, log[-3000:])
        remove_indexes(db)
        return db, log

    def test_build_estimates_the_genes_factors(self):
        rates = [0.3] * 6 + [1.7] * 6
        db, log = self.build("db", rates)
        self.assertIn("Gene conservation: factors", log)
        with open(os.path.join(db, "gene_conservation.tsv")) as fh:
            header = fh.readline()
            factors = {int(f[0]): float(f[1]) for f in (line.split("\t") for line in fh)}
        self.assertEqual(header, "geneid\tfactor\tspecies\n")
        self.assertEqual(sorted(factors), list(range(1, 13)))
        slow = [factors[g] for g in range(1, 7)]
        fast = [factors[g] for g in range(7, 13)]
        self.assertLess(max(slow), 1, factors)
        self.assertGreater(min(fast), 1, factors)
        # Between the 4 congeners, too, the slow genes differ least (gene_congeners.tsv, a report beside the database).
        self.assertRegex(log, r"Gene congeners: 6 pairs of species of 1 genera \(4 species\); the genes' divergence between "
                              r"congeners correlates 0\.[5-9]\d with their factors \(Spearman, 12 genes\)")
        rows = read_dicts(os.path.join(db, "gene_congeners.tsv"))
        between = {int(r["geneid"]): float(r["between_factor"]) for r in rows}
        self.assertLess(max(between[g] for g in range(1, 7)), min(between[g] for g in range(7, 13)), between)
        self.assertEqual({r["species"] for r in rows}, {"4"})

    def test_queries_scale_the_margin_only_when_asked(self):
        if not os.path.exists(db_file("gene_conservation.tsv")):
            raise unavailable("the test database has no gene_conservation.tsv (built by an earlier protal)")
        base = baseline()
        ones = self.path("ones.tsv")
        with open(db_file("gene_conservation.tsv")) as src, open(ones, "w") as dst:
            dst.write(src.readline())
            dst.writelines(line.split("\t")[0] + "\t1\n" for line in src)
        same = "for the conservation features; the depth identity margin is the same on every gene"
        self.assertIn(same, base.log)  # the default
        profiles = {}
        for name, extra, expected in (("scaled", ["--gene_conservation", "db"], "they scale the depth identity margin per gene"),
                                      ("none", ["--gene_conservation", "none"], same),
                                      ("ones", ["--gene_conservation", ones], "Gene conservation: factors 1-1 for")):
            rc, log = profile_only(self.work, self.path(name), [base.sam("sa")], *extra, prefixes=["sa"])
            self.assertEqual(rc, 0, log[-3000:])
            self.assertIn(expected, log)
            profiles[name] = read_text(self.path(name, "sa.profile"))
        # The reads are the reference genes' with 0.5% errors: far above any gene's threshold.
        self.assertEqual(profiles["scaled"], profiles["none"])
        self.assertEqual(profiles["ones"], profiles["none"])
        self.assertEqual(base.text("sa.profile"), profiles["none"])

        rc, log = profile_only(self.work, self.path("missing"), [base.sam("sa")], "--gene_conservation", "no_such.tsv",
                               prefixes=["sa"])
        self.assertNotEqual(rc, 0, log[-3000:])
        self.assertIn("--gene_conservation does not exist: no_such.tsv", log)


class GeneNeighboursTest(WorkDir):
    """A database of a synthetic release whose markers lie in operon-like clusters (simulate_gtdb_release.py
    --operons), with the gene neighbours of its genomes (gene_neighbours.py: the per-clade frequencies and where
    each gene lies in each genome), which --build checks and packs; read pairs drawn from the genomes span
    neighbouring genes. protal pairs mates across them (a proper pair on two references) and gives the adjacency
    features, from database.protal alone; --no_gene_neighbours does neither. (Needs protal, not the test database.)"""

    @classmethod
    def setUpClass(cls):
        require_binary(PROTAL, "protal", "PROTAL")
        super().setUpClass()
        scripts = os.path.join(ROOT, "scripts", "mini_db")
        cls.gtdb, cls.db = os.path.join(cls.work, "gtdb"), os.path.join(cls.work, "db")
        python(os.path.join(scripts, "simulate_gtdb_release.py"), "--outdir", cls.gtdb, "--operons", "--genome_length",
               "200000", "--genomes_per_species", "1", "--contigs", "1")
        python(os.path.join(scripts, "gtdb_to_protal_db.py"), "--gtdb", cls.gtdb, "--outdir", cls.db)
        genomes = os.path.join(cls.gtdb, "simulation", "genomes.tsv")
        python(os.path.join(scripts, "gene_neighbours.py"), "--db", cls.db, "--genome_table", genomes)
        cls.table = read_text(os.path.join(cls.db, "gene_neighbours.tsv"))
        cls.positions = read_text(os.path.join(cls.db, "gene_positions.tsv"))
        cls.inputs = os.path.join(cls.work, "inputs")  # the converted files, for builds of their own
        shutil.copytree(cls.db, cls.inputs)
        rc, cls.build_log = run(cls.work, "--build", "--no_profile", *FAST_BUILD, "--db", cls.db,
                                "--reference", os.path.join(cls.db, "reference.fna"))
        if rc != 0:
            raise RuntimeError("protal --build failed:\n" + cls.build_log[-3000:])
        with open(genomes) as fh:
            accessions = [row["accession"] for row in csv.DictReader(fh, delimiter="\t")]
        community = os.path.join(cls.work, "community.tsv")
        with open(community, "w") as fh:
            fh.write("accession\trelative_abundance\n")
            fh.writelines(f"{a}\t1\n" for a in accessions)
        python(os.path.join(scripts, "simulate_reads.py"), "--genomes", genomes, "--community", community,
               "--out_prefix", os.path.join(cls.work, "s"), "--pairs", "8000", "--seed", "5")
        with open(os.path.join(cls.inputs, "internal_taxonomy.dmp")) as fh:  # --build packed the database's copy
            next(fh)
            cls.species = [r[0] for r in (line.split("\t") for line in fh) if r[4] == "species"]

    def profile(self, name, *extra, db=None):
        truth = self.path("truth.tsv")
        with open(truth, "w") as fh:
            fh.write("\n".join(self.species) + "\n")
        rc, log = run(self.work, "--db", db or self.db, "-1", self.path("s_R1.fq"), "-2", self.path("s_R2.fq"), "--prefix", "s",
                      "-o", name, "-t", "2", "--no_strains", "--profile_truth", truth, *extra)
        self.assertEqual(rc, 0, log[-3000:])
        records = sam_records(find_sams(self.path(name, "s.sam"))[0])
        return log, records, read_dicts(self.path(name, "s.profile.truth_annotated"))

    @staticmethod
    def across(records):
        """The primary records of proper pairs whose mates are on two genes."""
        return [r for r in records if int(r[1]) & 0x2 and not int(r[1]) & 0x100 and r[6] not in ("=", "*")]

    def test_the_build_checks_and_packs_the_neighbours(self):
        self.assertIn("Gene neighbours: ", self.build_log)
        self.assertIn("stored in the database", self.build_log)
        # Each species' representative is one sequence: read as circular.
        n = len(self.species)
        self.assertIn(f"in {n} genomes of {n} species, {n} read as circular (gene_positions.tsv), stored in the database, "
                      "not read by queries", self.build_log)
        unpacked = self.path("unpacked")
        rc, log = run(self.work, "--unpack_db", "--db", self.db, "--unpack_dir", unpacked, "-t", "2")
        self.assertEqual(rc, 0, log[-3000:])
        self.assertEqual(read_text(os.path.join(unpacked, "gene_neighbours.tsv")), self.table)
        self.assertEqual(read_text(os.path.join(unpacked, "gene_positions.tsv")), self.positions)

    def test_bad_tables_stop_the_build(self):
        line = self.positions.splitlines()[-1].split("\t")
        for name, table, row, problem in (
                ("bad_positions", "gene_positions.tsv", "\t".join(line[:5] + ["999"] + line[6:]),
                 f"species {line[1]} has no gene 999 in the database"),
                ("bad_neighbours", "gene_neighbours.tsv",
                 self.table.splitlines()[-1].split("\t", 1)[0] + "\t999\t3\t0\t0\t1\t1\t0\t0\t0", "gene 999 is not in the database")):
            with self.subTest(table):
                db = self.path(name)
                shutil.copytree(self.inputs, db)
                with open(os.path.join(db, table), "a") as fh:
                    fh.write(row + "\n")
                rc, log = run(self.work, "--build", "--no_profile", "--no_bundle", "-t", "2", "--db", db,
                              "--reference", os.path.join(db, "reference.fna"))
                self.assertEqual(rc, 8, log[-3000:])
                self.assertIn(problem, log)
                remove_indexes(db)

    def test_mates_pair_across_neighbouring_genes(self):
        log, records, rows = self.profile("with")
        self.assertRegex(log, r"Gene neighbours: \d+ rules of \d+ clades")
        paired = int(re.search(r"Gene neighbours: (\d+) fragments paired across two neighbouring genes", log).group(1))
        self.assertGreater(paired, 50)
        across = self.across(records)
        self.assertGreater(len(across), paired)  # both mates of nearly every such fragment (a few fail the identity filter)
        self.assertLessEqual(len(across), 2 * paired)
        # Each such pair's genes face each other in its species' clade.
        adjacent = set()
        for line in self.table.splitlines():
            f = line.split("\t")
            if f[0].isdigit() and f[3] != "0":
                adjacent.add((f[1], f[3]))
        for r in across:
            taxid, gene = r[2].split("_")
            other_taxid, other = r[6].split("_")
            self.assertEqual(taxid, other_taxid)
            self.assertIn((gene, other), adjacent)
            self.assertEqual(r[8], "0")  # TLEN on two references
        present = [row for row in rows if row["truth"] == "1"]
        self.assertTrue(present)
        self.assertTrue(all(float(row["adjacent_expected_share"]) > 0.9 for row in present), present)

        # The database file alone, in a folder of its own: the same, to the last digit (the profile does not depend on
        # the threads or the order of the reads).
        alone = self.path("bundle_only")
        os.makedirs(alone)
        shutil.copy(os.path.join(self.db, "database.protal"), alone)
        log_alone, records_alone, rows_alone = self.profile("alone", db=alone)
        self.assertEqual(os.listdir(alone), ["database.protal"])  # nothing unpacked or written next to it
        self.assertIn(f"Gene neighbours: {paired} fragments paired across two neighbouring genes", log_alone)
        self.assertEqual(rows_alone, rows)
        self.assertEqual(sorted(records_alone), sorted(records))

        log, records, rows = self.profile("without", "--no_gene_neighbours")
        self.assertIn("Gene neighbours: not used (--no_gene_neighbours)", log)
        self.assertEqual(self.across(records), [])
        self.assertTrue(all(float(row["adjacent_expected_share"]) == 0 for row in rows))


class ModelContractTest(DbTest):
    """The training dump (the baseline's truth annotation) holds the features the model is scored with and every group
    the trainer can choose, and --no_strains changes no profile."""

    def test_truth_annotation_has_the_model_features(self):
        rows = read_dicts(baseline().path("sa.profile.truth_annotated"))
        self.assertTrue(rows)
        header = list(rows[0])
        with open(db_file("model_pe.xml")) as fh:
            fields = set(re.findall(r'<DataField name="([^"]+)"', fh.read())) - {"truth"}
        self.assertEqual(sorted(fields - set(header)), [], "every model input is in the training dump")
        sys.path.insert(0, os.path.join(ROOT, "scripts"))
        import model_features
        self.assertEqual([c for c in model_features.feature_set_columns(model_features.DEFAULT_FEATURE_SET)
                          if c not in header], [], "the trainer's default features are dumped")
        for group, columns in model_features.FEATURE_GROUPS.items():
            self.assertEqual([c for c in columns if c not in header], [], f"the {group} features are dumped")
        self.assertEqual(model_features.feature_columns(header, model_features.DEFAULT_FEATURE_SET),
                         model_features.feature_set_columns(model_features.DEFAULT_FEATURE_SET))
        conservation = os.path.exists(db_file("gene_conservation.tsv"))
        for row in rows:
            # The reference genes' reads with 0.5% substitutions and 1% Q2 bases: about what their qualities explain.
            self.assertLess(abs(float(row["excess_median"])), 0.02, row)
            self.assertLessEqual(float(row["excess_high_share"]), 0.5, row)
            share = float(row["conserved_hit_share"])
            self.assertTrue(0 <= share <= 1 if conservation else share == 0.5, row)
            for prefix, counts in (("RAF", "AF"), ("RA", "A")):
                total = sum(float(row[f"{counts}{i}"]) for i in range(5))
                for i in range(5):
                    expected = float(row[f"{counts}{i}"]) / total if total else 0
                    self.assertAlmostEqual(float(row[f"{prefix}{i}"]), expected, places=9, msg=f"{prefix}{i}")
            self.assertEqual(row["truth"], "1")
        # The sample's complexity: the same in every row of the sample, whose taxa (each with fragments) are the rows;
        # the reads have 0.5% substitutions.
        from model_features import COMPLEXITY_FEATURES
        for name in COMPLEXITY_FEATURES:
            self.assertEqual(len({row[name] for row in rows}), 1, name)
        self.assertAlmostEqual(float(rows[0]["sample_log_taxa"]), math.log10(len(rows)), places=12)
        self.assertTrue(0 <= float(rows[0]["sample_low_identity"]) <= 1, rows[0])
        self.assertTrue(0.98 < float(rows[0]["sample_identity"]) <= 1, rows[0])

    def test_no_strains_changes_no_profile(self):
        # The baseline built strain MSAs over its samples; the same SAMs profiled with --no_strains.
        base = baseline()
        base.place_sams(self.path("out"), ["sa", "sb"])
        rc, log = run(self.work, "--db", DB, *reads("sa", "sb"), "-o", "out", "-t", "2", "--no_qcmsa", "--no_strains")
        self.assertEqual(rc, 0, log[-3000:])
        self.assertIn("All alignments are present", log)
        self.assertFalse(os.path.exists(self.path("out", "strains")))
        for sample in ("sa", "sb"):
            for name, text in base.profiles(sample).items():
                self.assertEqual(read_text(self.path("out", name)), text, name)


def model_with_header(model_path, path, extensions):
    """The database's model with `extensions` (XML text) added to its header, written to path."""
    with open(model_path) as fh:
        xml = fh.read()
    xml, n = re.subn(r"(<Header\b[^>]*[^/]>)", lambda m: m.group(1) + "\n  " + extensions, xml, count=1)
    if n != 1:
        raise AssertionError(f"{model_path} has no header")
    with open(path, "w") as fh:
        fh.write(xml)
    return path


class KnobSample(DbTest):
    """Tests of how a sample's calls are made, on the baseline's sample thin: Mockella beta and Fakibacter gamma at 12
    pairs per gene, which the model reports, and one pair of Mockella alpha, which it rejects at the default knob
    0.5 but reports at knob 0. So the profile tells which knob was applied."""

    def profile(self, name, *extra):
        rc, log = profile_only(self.work, self.path(name), [baseline().sam("thin")], *extra, prefixes=["thin"])
        self.assertEqual(rc, 0, log[-3000:])
        return log, read_text(self.path(name, "thin.profile"))

    def fails(self, name, model_extensions, code, problem=None, *extra):
        """protal stops with `code` on the model with these header extensions, before profiling, saying `problem`
        (None: the message is a unit test's)."""
        model = model_with_header(db_file("model_pe.xml"), self.path(name + ".xml"), model_extensions)
        rc, log = profile_only(self.work, self.path(name), [baseline().sam("thin")], "--model", model, *extra,
                               prefixes=["thin"])
        self.assertEqual(rc, code, log[-3000:])
        self.assertFalse(os.path.exists(self.path(name, "thin.profile")), "no profile may be written")
        if problem is not None:
            self.assertIn(problem, log)

    def assert_knob_zero(self, profile):
        self.assertEqual(profile_species(profile), sorted(species()), "knob 0 reports Mockella alpha's one pair")

    def assert_default(self, profile):
        self.assertEqual(profile, baseline().text("thin.profile"))
        self.assertEqual(profile_species(profile), sorted(set(species()) - {"s__Mockella alpha"}))


class DepthKnobsTest(KnobSample):
    """A model with knobs by sample depth (machine_learning_cmdline.py --depth-knobs, in its header): a sample's taxa are
    reported at the knob of its depth bin unless --knob is given."""

    def model(self, name, value, extension="protal_depth_knobs", more=""):
        return model_with_header(db_file("model_pe.xml"), self.path(name),
                                 f'<Extension name="{extension}" value="{value}"/>' + more)

    def test_the_samples_depth_knob_unless_knob_is_given(self):
        every_bin = ",".join(f"{b}:0" for b in range(2, 7))
        model = self.model("knobs.xml", every_bin)
        log, by_depth = self.profile("by_depth", "--model", model)
        self.assertIn("knobs by sample depth (bin b: 10^b to 10^(b+1) fragments, 2 also fewer, 6 also more): "
                      "2: 0, 3: 0, 4: 0, 5: 0, 6: 0; other depths --knob 0.5", log)
        self.assertRegex(log, r"Sample thin: \d+ fragments, knob 0 \(the model's for depth bin [2-6]\)")
        self.assert_knob_zero(by_depth)
        _, at_zero = self.profile("at_zero", "--knob", "0")
        self.assertEqual(by_depth, at_zero)

        log, given = self.profile("given", "--model", model, "--knob", "0.5")
        self.assertIn("; not used, --knob is given", log)
        # No knob chosen for the sample ("Sample thin: N fragments, knob ..."); its alignment counts start with
        # "Sample thin: " too.
        self.assertNotRegex(log, r"Sample thin: \d+ fragments, ")
        self.assert_default(given)

    def test_a_knob_curve(self):
        # The trainer's knob curve (since 0.7.3): read at the sample's depth, linear between its points and the ends'
        # beyond them, so that every depth has a knob of the model's.
        curve = "protal_depth_knob_curve"
        model = self.model("curve.xml", "0.000:0,12.000:0", curve)
        log, curved = self.profile("curved", "--model", model)
        self.assertIn("knobs by sample depth (log10 of the sample's fragments: knob; linear between, the ends' beyond): "
                      "0: 0, 12: 0", log)
        self.assertRegex(log, r"Sample thin: \d+ fragments, knob 0 \(the model's for that depth\)")
        self.assert_knob_zero(curved)
        log, given = self.profile("curve_given", "--model", model, "--knob", "0.5")
        self.assertIn("; not used, --knob is given", log)
        self.assert_default(given)

    def test_malformed_depth_knobs(self):
        # The malformed values' messages are tested by ModelFeatures.MalformedDepthKnob*; here protal stops on them.
        for name, value, extension, more, problem in (
                ("bad", "2:0.3,9:0.4", "protal_depth_knobs", "", None),
                ("curve_bad", "3:0.2,2:0.4", "protal_depth_knob_curve", "", None),
                ("both", "2:0.2", "protal_depth_knob_curve", '<Extension name="protal_depth_knobs" value="2:0.3"/>',
                 "it has depth knobs twice")):
            with self.subTest(name):
                self.fails(name, f'<Extension name="{extension}" value="{value}"/>' + more, 2, problem)


class FalseCallsTest(KnobSample):
    """A model with calibrated calls (machine_learning_cmdline.py --fdr-calls, in its header): with --fdr F a sample
    reports its highest-scoring taxa while their expected share of false calls stays at F; without --fdr (or with
    --fdr 0 or --knob) the calls are not used and the knob curve or --knob applies."""

    @staticmethod
    def calls(curve="0:0,1:1", prior="0.5", fdr="0.05"):
        return (f'<Extension name="protal_calibration" value="{curve}"/><Extension name="protal_prior" value="{prior}"/>'
                f'<Extension name="protal_fdr" value="{fdr}"/>')

    def test_calls_only_with_fdr(self):
        model = model_with_header(db_file("model_pe.xml"), self.path("calls.xml"), self.calls(fdr="0.000001"))
        # Without --fdr the model's calibrated calls are not used, whatever its target: the profile is the default's.
        log, none = self.profile("none", "--model", model)
        self.assertIn("calls at an expected share of false calls of 1e-06 (calibrated, 2 points; training prior 0.5); "
                      "not used (--fdr F would use them)", log)
        self.assertNotIn("expected share of false calls of at most", log)
        self.assert_default(none)
        log, strict = self.profile("strict", "--model", model, "--fdr", "0.000001")
        self.assertIn("; at --fdr 1e-06, the depth knobs are not used", log)
        self.assertRegex(log, r"Sample thin: \d+ fragments, 0 taxa at an expected share of false calls of at most 1e-06")
        self.assertEqual(strict.strip(), "", "a target no taxon meets calls none")
        log, generous = self.profile("generous", "--model", model, "--fdr", "0.9")
        self.assertIn("; at --fdr 0.9", log)
        self.assertRegex(log, r"Sample thin: \d+ fragments, 3 taxa at an expected share of false calls of at most 0.9")
        self.assert_knob_zero(generous)
        for name, extra, said in (("off", ["--fdr", "0"], "; not used (--fdr F would use them)"),
                                  ("knob", ["--knob", "0.5"], "; not used, --knob is given")):
            log, profile = self.profile(name, "--model", model, *extra)
            self.assertIn(said, log)
            self.assert_default(profile)

    def test_fdr_needs_calibrated_calls_and_no_knob(self):
        thin = baseline().sam("thin")
        for name, extra, code, problem in (("uncalibrated", ["--fdr", "0.1"], 2, "--fdr needs a model with a calibration"),
                                           ("both", ["--fdr", "0.1", "--knob", "0.5"], None, "--fdr and --knob exclude each other"),
                                           ("range", ["--fdr", "1"], None, "--fdr must be at least 0 and below 1")):
            with self.subTest(name):
                rc, log = profile_only(self.work, self.path(name), [thin], *extra, prefixes=["thin"])
                if code is None:
                    self.assertNotEqual(rc, 0, log[-3000:])
                else:
                    self.assertEqual(rc, code, log[-3000:])
                self.assertIn(problem, log)

    def test_malformed_calibration(self):
        for name, extensions, problem in (("decreasing", self.calls(curve="0:0.5,1:0.4"), "its calibration is malformed"),
                                          ("partial", '<Extension name="protal_fdr" value="0.1"/>', "need all of protal_calibration"),
                                          ("prior", self.calls(prior="1"), "its prior is malformed")):
            with self.subTest(name):
                self.fails(name, extensions, 2, problem)


def gbm_node(value=0.0, feature=0, threshold=0.0, left=0, right=0, leaf=False, missing_left=True):
    """A node of a tree as HistGradientBoostingClassifier stores it (its predictor's nodes): x <= threshold goes left."""
    return {"value": value, "is_leaf": leaf, "feature_idx": feature, "num_threshold": threshold, "left": left,
            "right": right, "missing_go_to_left": missing_left}


class GradientBoostedModelTest(DbTest):
    """The default model type since 1750475: a gradient-boosted model as scripts/model_pmml.py's write_boosted exports it
    (a modelChain of regression trees into a logit RegressionModel), loaded by cPMML from a database
    (tests/e2e/data/model_gbm_small.xml as its model_pe.xml) and scored on the baseline's SAMs of sa and thin. The
    trees are set by hand, their leaves sums of powers of two: a taxon's probability is the logistic of the baseline
    plus the leaves its features reach, computed here from the features the dump (the truth annotation) gives. Every
    leaf is reached by one of the six taxa (the thresholds lie between the mini database's values), and the profiles
    report exactly the taxa at the knob (0.5) or above: some, not all."""

    MODEL = os.path.join(DATA, "model_gbm_small.xml")
    FEATURES = ["fragments", "hit_gene_fraction", "identity", "sample_log_fragments"]
    BASELINE = -0.5
    # Per round the nodes, as sklearn's predictors hold them; the first's leaves without the baseline.
    TREES = [
        # Few fragments (thin's Mockella alpha: 1), or reads on fewer of the genes (Mockella beta: 108 of 115 genes,
        # the others 109 of 116 and 110 of 117): absent.
        [gbm_node(feature=0, threshold=10.0, left=1, right=2), gbm_node(-2.0, leaf=True),
         gbm_node(feature=1, threshold=0.9394, left=3, right=4, missing_left=False), gbm_node(-0.5, leaf=True),
         gbm_node(1.5, leaf=True)],
        # Reads with 0.5% substitutions against reads without errors (thin's one pair of Mockella alpha).
        [gbm_node(feature=2, threshold=0.999, left=1, right=2, missing_left=False), gbm_node(0.5, leaf=True),
         gbm_node(-1.0, leaf=True)],
        # The sample's depth: thin (10^3.42 fragments) against sa (10^3.59).
        [gbm_node(feature=3, threshold=3.5, left=1, right=2), gbm_node(0.25, leaf=True), gbm_node(-0.25, leaf=True)],
        # A round whose best split gained nothing: one leaf.
        [gbm_node(0.125, leaf=True)],
    ]
    ANNOTATION = "tests/e2e/test_protal_e2e.py GradientBoostedModelTest: trees set by hand, not trained"

    @classmethod
    def leaf(cls, tree, row):
        """The index of the leaf of `tree` that a dump row's features reach."""
        i = 0
        while not tree[i]["is_leaf"]:
            value = float(row[cls.FEATURES[tree[i]["feature_idx"]]])
            i = tree[i]["left"] if value <= tree[i]["num_threshold"] else tree[i]["right"]
        return i

    @classmethod
    def probability(cls, row):
        raw = cls.BASELINE + sum(tree[cls.leaf(tree, row)]["value"] for tree in cls.TREES)
        return 1 / (1 + math.exp(-raw))

    @classmethod
    def exported(cls, path):
        """Write the trees with model_pmml.write_boosted, from a stand-in of a fitted HistGradientBoostingClassifier."""
        sys.path.insert(0, os.path.join(ROOT, "scripts"))
        import model_pmml

        class Predictor:
            def __init__(self, nodes):
                self.nodes = nodes

        class Model:
            classes_ = [0, 1]
            n_features_in_ = len(cls.FEATURES)
            n_trees_per_iteration_ = 1
            is_categorical_ = None
            _baseline_prediction = [[cls.BASELINE]]
            _predictors = [[Predictor(tree)] for tree in cls.TREES]

        model_pmml.write_boosted(Model(), cls.FEATURES, path, [cls.ANNOTATION])
        return path

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        base = baseline()
        cls.db = symlink_db(os.path.join(cls.work, "db"), extra={"model_pe.xml": cls.MODEL})
        truths = []
        for sample in ("sa", "thin"):
            truths.append(os.path.join(cls.work, f"{sample}.truth.tsv"))
            with open(truths[-1], "w") as fh:
                fh.write("".join(f"{t}\n" for t in species().values()))
        cls.rc, cls.log = profile_only(cls.work, os.path.join(cls.work, "out"), [base.sam("sa"), base.sam("thin")],
                                       "--profile_truth", ",".join(truths), prefixes=["sa", "thin"], db=cls.db)

    def rows(self, sample):
        return read_dicts(self.path("out", f"{sample}.profile.truth_annotated"))

    def test_the_fixture_is_what_the_exporter_writes(self):
        try:
            import numpy  # noqa: F401  (model_pmml needs it)
        except ImportError:
            raise unavailable("numpy is needed to run scripts/model_pmml.py")
        self.assertEqual(read_text(self.exported(self.path("exported.xml"))), read_text(self.MODEL),
                         "regenerate the fixture: GradientBoostedModelTest.exported(path)")

    def test_protal_scores_the_hand_computed_probabilities(self):
        self.assertEqual(self.rc, 0, self.log[-3000:])
        self.assertIn("Model of paired-end reads: " + os.path.join(self.db, "model_pe.xml"), self.log)
        called = rejected = 0
        reached = set()
        for sample in ("sa", "thin"):
            rows = self.rows(sample)
            self.assertEqual(len(rows), len(species()), sample)
            for row in rows:
                expected = self.probability(row)
                reached |= {(t, self.leaf(tree, row)) for t, tree in enumerate(self.TREES)}
                self.assertAlmostEqual(float(row["probability"]), expected, delta=1e-12, msg=f"{sample}, {row['taxon_name']}")
                self.assertEqual(row["prediction"], "1" if expected >= 0.5 else "0", f"{sample}, {row['taxon_name']}")
                called += expected >= 0.5
                rejected += expected < 0.5
            reported = {row["taxon_name"] for row in rows if row["prediction"] == "1"}
            self.assertEqual(profile_species(read_text(self.path("out", f"{sample}.profile"))), sorted(reported), sample)
        leaves = {(t, i) for t, tree in enumerate(self.TREES) for i, node in enumerate(tree) if node["is_leaf"]}
        self.assertEqual(sorted(leaves - reached), [], "leaves no taxon reaches: set the thresholds between the taxa's values")
        self.assertGreater(called, 0)
        self.assertGreater(rejected, 0)

    def test_the_trainers_scorer_agrees(self):
        # scripts/model_pmml.py's PmmlBoosted scores a model as cPMML does (the trainer's knob and calibration use it).
        try:
            import numpy as np
        except ImportError:
            raise unavailable("numpy is needed to run scripts/model_pmml.py")
        sys.path.insert(0, os.path.join(ROOT, "scripts"))
        import model_pmml
        model = model_pmml.load_model(self.MODEL)
        self.assertEqual(model.features, self.FEATURES)
        rows = self.rows("sa") + self.rows("thin")
        X = np.array([[float(row[f]) for f in self.FEATURES] for row in rows])
        np.testing.assert_allclose(model.predict(X), [float(row["probability"]) for row in rows], rtol=0, atol=1e-12)


class BuildIndexTest(ProtalTest):
    """--build of a few genes: it checks every k-mer against the full reference, also those whose core occurs once in the
    index; ambiguous bases (N, IUPAC codes) do not go into the index: k-mers whose window holds one are left out of the
    counting, the placing and the uniqueness check (AmbiguousKmers.* test the windows); a full reference with no other
    copies of the genes gives no gene conservation factors. And the passes that count and place the k-mers in -t
    threads (by key range, batches applied in reference order) give the index and unique_kmers.tsv of the one-thread
    passes (--serial_index_passes): the genes' k-mers in 1 KB batches, those of taxon 2's genes 3-6 in other batches
    than taxon 1's copies. Both builds compress the index alike (FAST_BUILD), so the same index gives the same
    index.prx.zst."""

    TAXONOMY = ("id\tparent_id\texternal_id\tname\trank\tlevel\trep_genome\n"
                "3\t3\t0\troot\tno rank\t0\t\n"
                "1\t3\t0\ts__Alpha one\tspecies\t7\tGCF_1\n"
                "2\t3\t0\ts__Beta two\tspecies\t7\tGCF_2\n")
    CODES = {3: "N", 4: "Y", 5: "R", 6: "k"}  # gene id: the ambiguous base in taxon 1's copy

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        rng = random.Random(8)
        genes = {(1, 1): None, (1, 2): None, (2, 1): None}
        for key in genes:
            genes[key] = "".join(rng.choice("ACGT") for _ in range(900))
        # Genes 3-6: taxon 2 has a clean copy, taxon 1 the same with one ambiguous base in the middle.
        for gene, code in cls.CODES.items():
            seq = "".join(rng.choice("ACGT") for _ in range(3000))
            genes[(2, gene)] = seq
            genes[(1, gene)] = seq[:1500] + code + seq[1501:]
        # Another genome of taxon 2 carries taxon 1's gene 2 unchanged.
        full = os.path.join(cls.work, "full_reference.fna")
        cls.logs, cls.outputs = {}, {}
        # --no_bundle: unique_kmers.tsv stays a file of its own.
        for name, extra in (("db", ["--index_batch_kb", "1"]), ("serial", ["--serial_index_passes"])):
            db = write_tiny_db(os.path.join(cls.work, name), genes, cls.TAXONOMY)
            with open(os.path.join(db, "reference.fna")) as src, open(full, "w") as dst:
                dst.write(src.read() + ">2_7\n" + genes[(1, 2)] + "\n")
            rc, cls.logs[name] = run(cls.work, "--build", "--no_bundle", "--no_profile", *FAST_BUILD, *extra, "--db", db,
                                     "--reference", os.path.join(db, "reference.fna"), "--full_reference", full)
            if rc != 0:
                raise RuntimeError("protal --build failed:\n" + cls.logs[name][-3000:])
            cls.outputs[name] = {f: digest(os.path.join(db, f)) for f in ("index.prx.zst", "unique_kmers.tsv")}
        cls.db = os.path.join(cls.work, "db")
        cls.uniques = {}  # (taxid, gene): (short unique, long unique, all k-mers)
        with open(os.path.join(cls.db, "unique_kmers.tsv")) as fh:
            for line in fh:
                f = line.split("\t")
                cls.uniques[(int(f[0]), int(f[1]))] = (int(f[2]), int(f[4]), int(f[8]))

    def test_the_index_is_the_same_in_any_threads_and_batches(self):
        self.assertEqual(self.outputs["db"], self.outputs["serial"])

    def test_a_gene_another_taxon_carries_is_not_unique(self):
        self.assertGreater(self.uniques[(1, 1)][0], 0, self.uniques)
        self.assertGreater(self.uniques[(2, 1)][0], 0, self.uniques)
        self.assertGreater(self.uniques[(1, 2)][2], 0, self.uniques)
        self.assertEqual(self.uniques[(1, 2)][:2], (0, 0), "no k-mer of gene 1_2 is unique to taxon 1")

    def test_windows_with_an_ambiguous_base_are_not_indexed(self):
        for gene, code in self.CODES.items():
            with self.subTest(code):
                # The 31 windows over the base are gone, so a fifth of them are fewer syncmers of the gene.
                lost = self.uniques[(2, gene)][2] - self.uniques[(1, gene)][2]
                self.assertGreater(lost, 0)
                self.assertLessEqual(lost, 31)
                # Taxon 1's copy shares every k-mer it has with taxon 2's; taxon 2's clean copy has exactly the k-mers
                # over the base to itself.
                self.assertEqual(self.uniques[(1, gene)][:2], (0, 0))
                self.assertEqual(self.uniques[(2, gene)][0], lost)

    def test_no_other_copies_give_no_conservation_factors(self):
        # The full reference holds the representatives' genes themselves and a copy of a gene the database lacks.
        self.assertIn("Gene conservation: no factors", self.logs["db"])
        self.assertFalse(os.path.exists(os.path.join(self.db, "gene_conservation.tsv")))


class QcmsaContractTest(WorkDir):
    """qcmsa counts only the samples that are in the MSA."""

    META_HEADER = ("sample\tgene_id\tvertical_coverage\tcounts_vcov1\tcounts_vcov2\tmulti_allelic\tfiltered\t"
                   "multi_rate_vcov1\tfiltered_rate_vcov1\tmulti_rate_vcov2\tfiltered_rate_vcov2\tmedian_vcov\t"
                   "hcov\tgene_length\tmean_vcov_nonzero\tmedian_vcov_nonzero\n")

    def qcmsa(self, *args):
        return run(self.work, QCMSA, *args, binary=sys.executable)

    def test_samples_missing_from_the_msa_do_not_count(self):
        with open(self.path("x.raw.msa.fna"), "w") as fh:
            fh.write(">x_reference\nACGTACGTAC\n>s1\nACGTACGTAC\n>s2\nACGTACGTAC\n")
        with open(self.path("x.raw.partition.txt"), "w") as fh:
            fh.write("DNA, gene1 = 1-10\n")
        with open(self.path("x.meta.tsv"), "w") as fh:
            fh.write(self.META_HEADER)
            for sample in ("s1", "s2", "s3", "s4"):  # s3 and s4 are not in the MSA
                fh.write(f"{sample}\t1\t5\t10\t10\t0\t0\t0\t0\t0\t0\t5\t1\t10\t5\t5\n")
        args = [self.path("x.raw.msa.fna"), self.path("x.raw.partition.txt"), self.path("x.meta.tsv")]

        rc, log = self.qcmsa(*args, "--prefix", self.path("two"), "--gene-min-samples", "2")
        self.assertEqual(rc, 0, log)
        self.assertIn("Loaded meta: 2 samples", log)
        self.assertFalse(os.path.exists(self.path("two.msa.fna")), "2 samples are not more than 2")

        rc, log = self.qcmsa(*args, "--prefix", self.path("one"), "--gene-min-samples", "1")
        self.assertEqual(rc, 0, log)
        self.assertTrue(os.path.exists(self.path("one.msa.fna")), log)

    def write_species(self, name, rows):
        """An MSA of one 8-column gene (rows: [(name, sequence)]) with meta rows for its samples."""
        with open(self.path(name + ".raw.msa.fna"), "w") as fh:
            fh.writelines(f">{n}\n{s}\n" for n, s in rows)
        with open(self.path(name + ".raw.partition.txt"), "w") as fh:
            fh.write("DNA, gene1 = 1-8\n")
        with open(self.path(name + ".meta.tsv"), "w") as fh:
            fh.write(self.META_HEADER)
            for sample, _ in rows[1:]:
                fh.write(f"{sample}\t1\t5\t8\t8\t0\t0\t0\t0\t0\t0\t5\t1\t8\t5\t5\n")
        return [self.path(name + ".raw.msa.fna"), self.path(name + ".raw.partition.txt"), self.path(name + ".meta.tsv")]

    def read_msa(self, path):
        with open(path) as fh:
            lines = fh.read().split()
        return dict(zip((n[1:] for n in lines[0::2]), lines[1::2]))

    def test_site_cleanup(self):
        # Columns: 1 constant; 2 only s1 differs (s1's own mutation); 3 only the reference differs;
        # 4 A and R (an ambiguity, so constant); 5 an insertion column with only N, '-' and IUPAC
        # codes; 6 s1 and s2 differ; 7 and 8 constant.
        args = self.write_species("y", [("y_reference", "AAGA-AAA"), ("s1", "ACAA-CAA"), ("s2", "AAARNCAA"),
                                        ("s3", "AAAA-AAA"), ("s4", "AAAAYAAA")])
        rc, log = self.qcmsa(*args, "--prefix", self.path("default"))
        self.assertEqual(rc, 0, log)
        msa = self.read_msa(self.path("default.msa.fna"))
        self.assertEqual(msa["s1"], "ACAACAA", "every column but the one without a base; singletons kept")
        self.assertEqual(msa["y_reference"], "AAGAAAA")

        rc, log = self.qcmsa(*args, "--prefix", self.path("parsimony"), "--min-parsimony-samples", "2")
        self.assertEqual(rc, 0, log)
        # Column 2 (one sample differs) goes; column 3 stays: the reference row is not a sample.
        self.assertEqual(self.read_msa(self.path("parsimony.msa.fna"))["s1"], "AAACAA")

        rc, log = self.qcmsa(*args, "--prefix", self.path("variable"), "--discard-constant")
        self.assertEqual(rc, 0, log)
        self.assertEqual(self.read_msa(self.path("variable.msa.fna"))["s1"], "CAC", "A and R count as constant")

    def test_duplicate_names_stop_qcmsa(self):
        args = self.write_species("z", [("z_reference", "AAAAAAAA"), ("s1", "ACAAACAA"), ("s1", "AAAAAAAA")])
        rc, log = self.qcmsa(*args, "--prefix", self.path("dup"))
        self.assertNotEqual(rc, 0, log)
        self.assertIn("names 1 sequence(s) more than once (s1)", log)

    def test_coverage_gate_reads_the_msa(self):
        # The meta says every cell is fully covered; s3's row writes 2 of the gene's 8 positions, below
        # --gene-min-hcov 0.3, so its cell is gap-filled.
        args = self.write_species("c", [("c_reference", "AAAAAAAA"), ("s1", "ACAAAAAA"), ("s2", "AAAAAAAA"),
                                        ("s3", "AC------")])
        rc, log = self.qcmsa(*args, "--prefix", self.path("c"))
        self.assertEqual(rc, 0, log)
        self.assertEqual(self.read_msa(self.path("c.msa.fna"))["s3"], "--------", log)
        rc, log = self.qcmsa(*args, "--prefix", self.path("c2"), "--gene-min-hcov", "0.2")
        self.assertEqual(rc, 0, log)
        self.assertEqual(self.read_msa(self.path("c2.msa.fna"))["s3"], "AC------", log)

    def test_multi_allelic_filter_judges_rates(self):
        # 8 samples, 4 genes of 1000 positions with >= 2 reads. s1-s6 hold no IUPAC code; s7, deep,
        # holds 1 per gene (0.1%, noise); s8, a mixture, 10 per gene (1%). By counts of multi-allelic
        # genes s7 and s8 are alike; by their rates, only s8 is above the 0.2% floor.
        samples = [f"s{i}" for i in range(1, 9)]
        multi = {"s7": 1, "s8": 10}
        with open(self.path("m.raw.msa.fna"), "w") as fh:
            fh.write(">m_reference\n" + "A" * 32 + "\n")
            for i, s in enumerate(samples):
                fh.write(f">{s}\n" + ("A" * 7 + "ACGT"[i % 4]) * 4 + "\n")
        with open(self.path("m.raw.partition.txt"), "w") as fh:
            fh.writelines(f"DNA, gene{g} = {8 * g - 7}-{8 * g}\n" for g in range(1, 5))
        with open(self.path("m.meta.tsv"), "w") as fh:
            fh.write(self.META_HEADER)
            for s in samples:
                m = multi.get(s, 0)
                for g in range(1, 5):
                    fh.write(f"{s}\t{g}\t20\t1000\t1000\t{m}\t0\t{m / 1000}\t0\t{m / 1000}\t0\t20\t1\t1000\t20\t20\n")
        args = [self.path("m.raw.msa.fna"), self.path("m.raw.partition.txt"), self.path("m.meta.tsv")]
        rc, log = self.qcmsa(*args, "--prefix", self.path("m"))
        self.assertEqual(rc, 0, log)
        msa = self.read_msa(self.path("m.msa.fna"))
        self.assertNotIn("s8", msa, log)
        self.assertEqual(sorted(msa), sorted(["m_reference"] + samples[:7]), log)
        self.assertEqual(len(msa["s7"]), 32, "no gene removed or masked")
        self.assertIn("sample_filtered\ts8\t4\tmulti-allelic rate 0.0100 > 0.0020", read_text(self.path("m.qcmsa_summary.tsv")))

        # With a floor of 0.05%, the fence (0.0625%) decides, and s7 goes too.
        rc, log = self.qcmsa(*args, "--prefix", self.path("floor"), "--mrate2-min-rate", "0.0005")
        self.assertEqual(rc, 0, log)
        self.assertEqual(sorted(self.read_msa(self.path("floor.msa.fna"))), sorted(["m_reference"] + samples[:6]))

    def test_no_msa_leaves_no_stale_output(self):
        args = self.write_species("w", [("w_reference", "AAAAAAAA"), ("s1", "ACAAACAA"), ("s2", "AAAAAAAA")])
        rc, log = self.qcmsa(*args, "--prefix", self.path("w"))
        self.assertEqual(rc, 0, log)
        self.assertTrue(os.path.exists(self.path("w.msa.fna")), log)
        rc, log = self.qcmsa(*args, "--prefix", self.path("w"), "--gene-min-samples", "5")
        self.assertEqual(rc, 0, log)
        self.assertFalse(os.path.exists(self.path("w.msa.fna")), "an earlier run's MSA is removed")
        self.assertFalse(os.path.exists(self.path("w.partition.txt")))
        self.assertIn("status\tno_msa\t\tevery gene was filtered", read_text(self.path("w.qcmsa_summary.tsv")))


class MapUtilsTest(DbTest):
    """protal_map_utils resolves relative map paths as protal does."""

    def test_relative_paths_resolve_like_protal(self):
        os.makedirs(self.path("maps"))
        os.makedirs(self.path("reads"))
        maps = {}
        for sample in ("sa", "sb"):
            for mate in (1, 2):  # a few reads: protal must find them, not profile them
                head_reads(os.path.join(READS, f"{sample}_R{mate}.fq"), self.path("reads", f"{sample}_R{mate}.fq"), 100)
            maps[sample] = self.path("maps", f"{sample}.map")
            with open(maps[sample], "w") as fh:
                fh.write("#OUTPUT_DIR\tout\n#SAM_OUTPUT_DIR\taln\n#SAMPLEID\tPREFIX\tFIRST\tSECOND\tSAM\n")
                fh.write(f"{sample}\t{sample}\treads/{sample}_R1.fq\treads/{sample}_R2.fq\t{sample}.sam\n")
        sample_map = maps["sa"]
        tool = os.path.join(ROOT, "scripts", "protal_map_utils")

        # Read paths are relative to the directory protal runs in, not to the map's directory.
        rc, log = run(self.work, "validate", "--map", sample_map, binary=tool)
        self.assertEqual(rc, 0, log)

        rc, log = run(self.work, "--db", DB, "--map", sample_map, "-t", "2", "--no_qcmsa", "--no_profile")
        self.assertEqual(rc, 0, log[-3000:])
        self.assertTrue(sam_records(self.path("out", "aln", "sa.sam")),
                        "protal found the reads and wrote OUTPUT_DIR/SAM_OUTPUT_DIR")

        rc, merged = run(self.work, "merge", "--map", maps["sa"], maps["sb"], binary=tool)
        self.assertEqual(rc, 0, merged)
        lines = merged.splitlines()
        variables = dict(line.split("\t", 1) for line in lines if line.startswith("#") and not line.startswith("#SAMPLEID"))
        header = next(line for line in lines if line.startswith("#SAMPLEID")).split("\t")
        rows = [dict(zip(header, line.split("\t"))) for line in lines if line and not line.startswith("#")]
        self.assertEqual(len(rows), 2)
        for row in rows:
            first = os.path.join(variables.get("#INPUT_DIR", self.work), row["FIRST"])
            self.assertEqual(os.path.realpath(first), os.path.realpath(self.path("reads", row["#SAMPLEID"] + "_R1.fq")))
        # The SAM sa's run wrote stays where it is, by its absolute path, so that a run of the merged map profiles it;
        # sb has none yet, and keeps the name its map gives it, in the merged map's SAM folder.
        sams = {row["#SAMPLEID"]: row["SAM"] for row in rows}
        self.assertEqual(os.path.realpath(sams["sa"]), os.path.realpath(self.path("out", "aln", "sa.sam")))
        self.assertEqual(sams["sb"], "sb.sam")

        def merged_sams(*options):
            rc, merged = run(self.work, "merge", "--map", maps["sa"], maps["sb"], *options, binary=tool)
            self.assertEqual(rc, 0, merged)
            return {line.split("\t")[0]: line.split("\t")[header.index("SAM")] for line in merged.splitlines()
                    if line and not line.startswith("#")}

        # With --use-sampleid merge names the new SAMs itself; the ending chooses protal's output format. --new-sams
        # gives every sample a new one, as merge did up to 0.7.8.
        for option, ending in ((None, ".sam.zst"), ("--zstd", ".sam.zst"), ("--gzip", ".sam.gz"), ("--nogzip", ".sam")):
            formats = [option] if option else []
            self.assertEqual(merged_sams("--use-sampleid", "--new-sams", *formats), {"sa": "sa" + ending, "sb": "sb" + ending}, option)
            self.assertEqual(merged_sams("--use-sampleid", *formats)["sb"], "sb" + ending, option)
        self.assertEqual(merged_sams("--new-sams"), {"sa": "sa.sam", "sb": "sb.sam"})
        rc, log = run(self.work, "merge", "--map", maps["sa"], maps["sb"], "--use-sampleid", "--nogzip", "--zstd", binary=tool)
        self.assertNotEqual(rc, 0, "one format only")

        # Once sb has its SAM too, a run of the merged map aligns nothing: it profiles both runs' SAMs where they are.
        with open(self.path("out", "aln", "sb.sam"), "w") as fh:
            fh.write(sam_text(baseline().sam("sb")))
        with open(self.path("merged.map"), "w") as fh:
            fh.write(run(self.work, "merge", "--map", maps["sa"], maps["sb"], "--out", "merged", binary=tool)[1])
        rc, log = run(self.work, "validate", "--map", self.path("merged.map"), binary=tool)
        self.assertEqual(rc, 0, log)
        rc, log = run(self.work, "--db", DB, "--map", self.path("merged.map"), "-t", "2", "--no_qcmsa", "--no_profile")
        self.assertEqual(rc, 0, log[-3000:])
        self.assertIn("All alignments are present", log)
        self.assertFalse(os.path.exists(self.path("merged", "alignments", "sa.sam")))

    def test_generate_names_samples_without_mate_numbers(self):
        # The sample ID of a pair is its files' name without the mate number, an R before it and Illumina's chunk
        # number (x_R1.fq once gave x_R, Illumina's x_S1_L001_R1_001 the whole name); .fq.zst files are found too.
        tool = os.path.join(ROOT, "scripts", "protal_map_utils")
        expected = {"a_R1.fq": "a", "b_1.fastq.gz": "b", "c.R1.fq.gz": "c", "d_S1_L001_R1_001.fastq.gz": "d_S1_L001",
                    "e_R1_trimmed.fq": "e_trimmed", "f_1.fq.zst": "f"}
        os.makedirs(self.path("names"))
        for first in expected:
            second = first.replace("R1", "R2") if "R1" in first else first.replace("_1.", "_2.")
            for name in (first, second):
                open(self.path("names", name), "w").close()
        rc, out = run(self.work, "generate", "--input", self.path("names"), binary=tool)
        self.assertEqual(rc, 0, out)
        rows = [line.split("\t") for line in out.splitlines() if line and not line.startswith("#")]
        self.assertEqual(sorted(row[0] for row in rows), sorted(expected.values()), out)

    def test_merging_single_end_maps(self):
        # '-' as SECOND (single-end reads) stays '-', and a map without a SECOND column merges too.
        tool = os.path.join(ROOT, "scripts", "protal_map_utils")
        for name, header, row in (("a", "#SAMPLEID\tFIRST\tSECOND\tPREFIX", "x\treads/x.fq\t-\tx"),
                                  ("b", "#SAMPLEID\tFIRST\tSECOND\tPREFIX", "y\treads/y.fq\t-\ty"),
                                  ("c", "#SAMPLEID\tFIRST\tPREFIX", "x\treads/x.fq\tx"),
                                  ("d", "#SAMPLEID\tFIRST\tPREFIX", "y\treads/y.fq\ty")):
            with open(self.path(f"{name}.map"), "w") as fh:
                fh.write(f"#OUTPUT_DIR\tout_{name}\n{header}\n{row}\n")
        rc, merged = run(self.work, "merge", "--map", self.path("a.map"), self.path("b.map"), "--out", "m", binary=tool)
        self.assertEqual(rc, 0, merged)
        rows = [line.split("\t") for line in merged.splitlines() if line and not line.startswith("#")]
        self.assertEqual([row[2] for row in rows], ["-", "-"], merged)
        rc, merged = run(self.work, "merge", "--map", self.path("c.map"), self.path("d.map"), "--out", "m", binary=tool)
        self.assertEqual(rc, 0, merged)

    def test_sample_ids_that_cannot_name_files(self):
        # #SAMPLEID names files (misc/<sample>_runtime.tsv) and MSA rows: every unsafe one is reported as the map is read.
        sam = sam_path(baseline().path("sa.sam"))
        ids = ["a/b", "c:d", "-e", "..", "f\x01g", "x" * 201, "fine"]
        with open(self.path("bad.map"), "w", encoding="utf-8") as fh:
            fh.write(f"#OUTPUT_DIR\tbad\n#SAMPLEID\tFIRST\tSECOND\tSAM\tPREFIX\n")
            for i, sample_id in enumerate(ids, 1):
                fh.write(f"{sample_id}\tgone_1.fq\tgone_2.fq\t{sam}\tp{i}\n")
        rc, log = run(self.work, "--db", DB, "--map", self.path("bad.map"), "-t", "2", "--no_qcmsa")
        self.assertEqual(rc, 9, log[-3000:])
        for line, reason in ((3, "contains '/'"), (4, "contains ':'"), (5, "starts with '-'"), (6, "names a folder"),
                             (7, "control character"), (8, "201 bytes long")):
            self.assertIn(f"Line {line}: sample ID", log)
            self.assertIn(reason, log)
        self.assertNotIn("Line 9:", log)


class RerunTest(DbTest):
    """Existing SAM files are reused (left as they are), and --no_profile is honoured when all of them exist."""

    def test_reruns(self):
        # The files' times are set an hour back: a file written again gets a time of now.
        baseline().place_sams(self.path("out"), ["sa", "sb"])
        args = ["--db", DB, *reads("sa", "sb"), "-o", "out", "-t", "4", "--no_qcmsa"]
        sam = self.path("out", "sa.sam.zst")
        with open(sam, "rb") as fh:
            before = fh.read()
        past = int(time.time()) - 3600
        os.utime(sam, (past, past))

        rc, log = run(self.work, *args)
        self.assertEqual(rc, 0, log[-3000:])
        self.assertIn("All alignments are present", log)
        self.assertEqual(os.path.getmtime(sam), past, "SAM untouched by the rerun")
        with open(sam, "rb") as fh:
            self.assertEqual(fh.read(), before)
        self.assertEqual(read_text(self.path("out", "sa.profile")), baseline().text("sa.profile"))

        os.utime(self.path("out", "sa.profile"), (past, past))
        rc, log = run(self.work, *args, "--no_profile")
        self.assertEqual(rc, 0, log[-3000:])
        self.assertEqual(os.path.getmtime(self.path("out", "sa.profile")), past, "--no_profile leaves the profiles alone")


class FailureTest(DbTest):
    """Failures that do not stop the run are summarised and make protal exit 1."""

    def test_mismatched_read_counts(self):
        with open(os.path.join(READS, "sa_R1.fq")) as fh:
            r1 = fh.readlines()[:400]   # 100 records
        with open(os.path.join(READS, "sa_R2.fq")) as fh:
            r2 = fh.readlines()[:396]   # 99 records
        with open(self.path("bad_R1.fq"), "w") as fh:
            fh.writelines(r1)
        with open(self.path("bad_R2.fq"), "w") as fh:
            fh.writelines(r2)
        baseline().place_sams(self.path("out_bad"), ["sa"], ["good"])
        rc, log = run(self.work, "--db", DB, "-1", f"{READS}/sa_R1.fq,{self.path('bad_R1.fq')}",
                      "-2", f"{READS}/sa_R2.fq,{self.path('bad_R2.fq')}", "--prefix", "good,bad",
                      "-o", "out_bad", "-t", "2", "--no_qcmsa")
        self.assertEqual(rc, 1, log[-3000:])
        self.assertRegex(log, r"protal finished with \d+ error")
        self.assertEqual(read_text(self.path("out_bad", "good.profile")), baseline().text("sa.profile"),
                         "the other sample is still profiled")
        self.assertFalse(find_sams(self.path("out_bad", "bad.sam")))

    def test_missing_qcmsa(self):
        # Two samples, so that there are strain MSAs for qcmsa to filter.
        baseline().place_sams(self.path("out_noqc"), ["sa", "sb"])
        rc, log = run(self.work, "--db", DB, *reads("sa", "sb"), "-o", "out_noqc", "-t", "2",
                      "--qcmsa_script", self.path("no", "such", "qcmsa"))
        self.assertEqual(rc, 1, log[-3000:])
        self.assertIn("qcmsa not found", log)


class BadInputTest(DbTest):
    """Read files protal cannot use stop it before aligning (exit 30) or fail their sample (exit 1) with the reason,
    instead of giving an empty profile and exit 0, and the other samples of the run are profiled; pipes and
    compressed files are read as plain files are. One run of the readable variants, one of the broken ones."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        r1, r2 = os.path.join(READS, "sa_R1.fq"), os.path.join(READS, "sa_R2.fq")
        with open(r1, "rb") as fh:
            data1 = fh.read()
        with open(r2, "rb") as fh:
            data2 = fh.read()
        # Readable: sa's reads compressed with zstd, two frames one after the other (cat a.zst b.zst), as files (zst)
        # and from pipes (zpipe); from pipes, gzipped and plain (pipe, as from process substitution, -1 <(zcat
        # a.fq.gz)); with mates of other names (names); and empty files (empty).
        cls.samples = {}
        cls.feeds = []
        cls.zstd = shutil.which("zstd") is not None
        if cls.zstd:
            frames = []
            for data in (data1, data2):
                half = data.index(b"\n@", len(data) // 2) + 1
                frames.append(b"".join(subprocess.run(["zstd", "-q", "-c", "-3"], input=part, check=True,
                                                      stdout=subprocess.PIPE).stdout for part in (data[:half], data[half:])))
            cls.samples["zst"] = [cls.write("z_R1.fq.zst", frames[0]), cls.write("z_R2.fq.zst", frames[1])]
            cls.samples["zpipe"] = [cls.fifo("zpipe_R1.fq.zst", frames[0]), cls.fifo("zpipe_R2.fq.zst", frames[1])]
        cls.samples["pipe"] = [cls.fifo("pipe_R1.fq.gz", gzip.compress(data1)), cls.fifo("pipe_R2.fq", data2)]
        lines = data2.decode().splitlines(keepends=True)
        cls.renamed = len(lines) // 4
        for i in range(0, len(lines), 4):
            lines[i] = lines[i].replace("@sa.", "@other.", 1)
        cls.samples["names"] = [r1, cls.write("renamed_R2.fq", "".join(lines).encode())]
        cls.samples["empty"] = [cls.write("empty_R1.fq", b""), cls.write("empty_R2.fq", b"")]
        writers = [threading.Thread(target=cls.feed, args=feed, daemon=True) for feed in cls.feeds]
        for writer in writers:
            writer.start()
        cls.rc, cls.log = cls.sample_run("out", cls.samples)
        for writer in writers:
            writer.join(timeout=10)

        # Broken, beside a good sample (sa's reads): bzip2-compressed reads, which protal does not read (they go through
        # a pipe: <(bzcat ...)); corrupt zstd; a damaged gzip member.
        half = len(data1) // 2
        cls.broken = {
            "bz2": [cls.write(f"b_R{m}.fq.bz2", b"BZh91AY&SY" + random.Random(m).randbytes(4000)) for m in (1, 2)],
            "zstd": [cls.write(f"c_R{m}.fq.zst", b"\x28\xb5\x2f\xfd" + random.Random(m + 2).randbytes(4000)) for m in (1, 2)],
            "gzip": [cls.write("damaged_R1.fq.gz", gzip.compress(data1[:half], mtime=0) + b"\x1fX" +
                               gzip.compress(data1[half:], mtime=0)[2:]), r2],
            "good": [r1, r2]}
        cls.broken_rc, cls.broken_log = cls.sample_run("out_broken", cls.broken)

    @classmethod
    def write(cls, name, data):
        path = os.path.join(cls.work, name)
        with open(path, "wb") as fh:
            fh.write(data)
        return path

    @classmethod
    def fifo(cls, name, data):
        path = os.path.join(cls.work, name)
        os.mkfifo(path)
        cls.feeds.append((data, path))
        return path

    @staticmethod
    def feed(data, fifo):
        with open(fifo, "wb") as fout:
            fout.write(data)

    @classmethod
    def sample_run(cls, out, samples):
        return run(cls.work, "--db", DB, "-1", ",".join(f for f, _ in samples.values()),
                   "-2", ",".join(s for _, s in samples.values()), "--prefix", ",".join(samples), "-o", out, "-t", "2",
                   "--no_qcmsa", "--no_strains", timeout=300)

    def test_compressed_reads_and_pipes_profile_as_the_plain_files(self):
        self.assertEqual(self.rc, 0, self.log[-3000:])
        for sample in ("zst", "zpipe", "pipe"):
            if sample in self.samples:
                with self.subTest(sample):
                    self.assertEqual(read_text(self.path("out", f"{sample}.profile")), baseline().text("sa.profile"))
        if not self.zstd:
            raise unavailable("the zstd command is needed to write zstd-compressed reads")

    def test_mates_of_other_names_are_warned_about(self):
        self.assertRegex(self.log, rf"Warning: sample names: the mates of {self.renamed} read pair\(s\) have different names "
                                   r"\(e\.g\. sa\.(\d+)/1, other\.\1/2\)")

    def test_empty_read_files_are_warned_about(self):
        self.assertIn("Warning: sample empty has no reads", self.log)

    def test_unusable_reads_fail_their_sample(self):
        # Which sample failed, and that it alone did (the reasons' wording: SeqReader's unit tests).
        self.assertEqual(self.broken_rc, 1, self.broken_log[-3000:])
        self.assertIn("Reading the FASTQ files of sample bz2 failed", self.broken_log)
        self.assertIn("The FASTQ files of sample zstd are truncated or corrupt", self.broken_log)
        self.assertIn("The FASTQ files of sample gzip are truncated or corrupt", self.broken_log)
        for sample in ("bz2", "zstd", "gzip"):
            self.assertFalse(glob.glob(self.path("out_broken", f"{sample}.sam*")), f"{sample}: no SAM may be written")
        self.assertEqual(read_text(self.path("out_broken", "good.profile")), baseline().text("sa.profile"),
                         "the good sample is profiled as on its own")

    def test_an_unreadable_read_file_stops_protal(self):
        r1, r2 = os.path.join(READS, "sa_R1.fq"), os.path.join(READS, "sa_R2.fq")
        rc, log = run(self.work, "--db", DB, "-1", self.work, "-2", r2, "--prefix", "s", "-o", "out_dir")
        self.assertEqual(rc, 30, log[-3000:])
        self.assertIn(f"-1 file cannot be read: {self.work}", log)
        if os.geteuid() != 0:  # root reads any file
            locked = self.path("locked_R2.fq")
            shutil.copy(r2, locked)
            os.chmod(locked, 0)
            rc, log = run(self.work, "--db", DB, "-1", r1, "-2", locked, "--prefix", "s", "-o", "out_locked")
            self.assertEqual(rc, 30, log[-3000:])
            self.assertIn(f"-2 file cannot be read: {locked}", log)
        for out in ("out_dir", "out_locked"):
            self.assertFalse(glob.glob(self.path(out, "*.sam*")), "no read may be aligned")


def is_seekable(path):
    """True if a zstd file ends with a seek table (zstd seekable format, as protal writes)."""
    with open(path, "rb") as fh:
        fh.seek(-4, os.SEEK_END)
        return fh.read(4) == b"\xb1\xea\x92\x8f"


class CompressedDatabaseTest(DbTest):
    """Raw, seekable (--compress_db --no_bundle, loaded in parallel), single-frame (zstd CLI) and
    single-file (--compress_db: database.protal) copies of the database give identical results,
    with one thread or several. The raw index (~3 GB, a fixed-size key map) lives only while setUpClass makes the
    copies and the raw run; its digest stands for it after."""

    KINDS = ("raw", "seekable", "single", "bundle")
    FILES = ("index.prx", "reference.fna", "reference.map", "internal_taxonomy.dmp", "unique_kmers.tsv", "model_pe.xml")
    BIG = ("index.prx", "reference.fna")

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        require_zstd("to make raw and compressed copies of the database")
        cls.dbs = {kind: os.path.join(cls.work, f"{kind}_db") for kind in cls.KINDS}
        for d in cls.dbs.values():
            os.mkdir(d)
        for f in glob.glob(os.path.join(FILES, "*")):
            name = os.path.basename(f)
            if not name.startswith(cls.BIG) and name != "database.protal":
                for d in cls.dbs.values():
                    os.symlink(f, os.path.join(d, name))
        # A raw copy: --decompress_db on symlinks to the database's files (a compressed index is in
        # protal's column format, which zstd -d does not turn back into index.prx).
        for name in cls.BIG:
            path = db_file(name)
            os.symlink(path, os.path.join(cls.dbs["raw"], os.path.basename(path)))
        cls.decompress_rc, cls.decompress_log = run(cls.work, "--decompress_db", "--db", cls.dbs["raw"], "-t", "4")
        if cls.decompress_rc != 0:
            raise RuntimeError("protal --decompress_db failed:\n" + cls.decompress_log[-3000:])
        raw_index = os.path.join(cls.dbs["raw"], "index.prx")
        cls.raw_digest = digest(raw_index)
        cls.raw_files = {name: (os.path.isfile(os.path.join(cls.dbs["raw"], name)),
                                os.path.lexists(os.path.join(cls.dbs["raw"], name + ".zst"))) for name in cls.BIG}
        for name in cls.BIG:
            raw = os.path.join(cls.dbs["raw"], name)
            # One frame, no seek table (-f: raw may be a symlink, which the zstd CLI skips otherwise).
            subprocess.run(["zstd", "-q", "-f", "-1", "-T0", raw, "-o", os.path.join(cls.dbs["single"], name + ".zst")],
                           check=True)
            os.symlink(raw, os.path.join(cls.dbs["seekable"], name))
            os.symlink(raw, os.path.join(cls.dbs["bundle"], name))
        # --compress_db --no_bundle turns the symlinked raw files into .zst files (removing the
        # links): the index in the column format, the reference seekable. Without --no_bundle, it
        # packs all database files into database.protal.
        small = ("--compress_level", "3", "--compress_frame_mb", "1")
        cls.compress_rc, cls.compress_log = run(cls.work, "--compress_db", "--no_bundle", "--db", cls.dbs["seekable"],
                                                "-t", "4", *small)
        cls.bundle_rc, cls.bundle_log = run(cls.work, "--compress_db", "--db", cls.dbs["bundle"], "-t", "4", *small)
        cls.bundle = os.path.join(cls.dbs["bundle"], "database.protal")
        cls.raw_kept = {name: os.path.exists(os.path.join(cls.dbs["raw"], name)) for name in cls.BIG}
        cls.expected = cls.result(cls.dbs["raw"], "out_raw_1", 1)
        os.remove(raw_index)

    @classmethod
    def result(cls, db, out, threads, *extra):
        """Sorted SAM records, profile and log of sample sa on db."""
        rc, log = run(cls.work, "--db", db, *reads("sa"), "-o", out, "-t", str(threads), "--no_qcmsa", *extra)
        if rc != 0:
            raise AssertionError(f"protal exited {rc} on {db}:\n" + log[-3000:])
        with open_sam(find_sams(os.path.join(cls.work, out, "sa.sam"))[0]) as sam:
            records = sorted(line for line in sam if not line.startswith("@"))
        return records, read_text(os.path.join(cls.work, out, "sa.profile")), log

    def test_decompress_db(self):
        for name in self.BIG:
            self.assertEqual(self.raw_files[name], (True, False), f"raw {name}, no {name}.zst")
        with open(os.path.join(self.dbs["raw"], "reference.map"), "rb") as fh:
            ends = [int(line.split()[3]) for line in fh]
        self.assertGreaterEqual(os.path.getsize(os.path.join(self.dbs["raw"], "reference.fna")), max(ends))

    def test_compress_db(self):
        self.assertEqual(self.compress_rc, 0, self.compress_log[-3000:])
        for name in self.BIG:
            self.assertFalse(os.path.lexists(os.path.join(self.dbs["seekable"], name)), f"{name} replaced")
            self.assertTrue(is_seekable(os.path.join(self.dbs["seekable"], name + ".zst")), f"{name}.zst is seekable")
            self.assertFalse(is_seekable(os.path.join(self.dbs["single"], name + ".zst")))
            self.assertTrue(self.raw_kept[name], "the raw files stay")
        self.assertFalse(os.path.exists(os.path.join(self.dbs["seekable"], "database.protal")))
        head = subprocess.run(["zstd", "-dc", os.path.join(self.dbs["seekable"], "index.prx.zst")],
                              stdout=subprocess.PIPE).stdout[:8]
        self.assertEqual(head, b"PRXSPLT1", "the index is in the column format")
        rc, log = run(self.work, "--compress_db", "--no_bundle", "--db", self.dbs["seekable"], "-t", "2")
        self.assertEqual(rc, 0, log[-3000:])
        self.assertIn("already in the column format; kept", log)
        self.assertIn("already seekable; kept", log)

    def test_single_file_db(self):
        """--compress_db packs the files into database.protal and removes them."""
        self.assertEqual(self.bundle_rc, 0, self.bundle_log[-3000:])
        self.assertTrue(is_seekable(self.bundle))
        for name in self.FILES:
            for variant in (name, name + ".zst"):
                self.assertFalse(os.path.lexists(os.path.join(self.dbs["bundle"], variant)), f"{variant} is packed")
        head = subprocess.run(["zstd", "-dc", self.bundle], stdout=subprocess.PIPE).stdout[:8]
        self.assertEqual(head, b"PROTALDB", "the file starts with its directory")
        # A single file with a current gene table (--compress_db packs one in) is left as it is.
        with open(self.bundle, "rb") as fh:
            before = fh.read()
        rc, log = run(self.work, "--compress_db", "--db", self.dbs["bundle"], "-t", "2")
        self.assertEqual(rc, 0, log[-3000:])
        self.assertIn("is a single-file database with a current gene_table.bin; kept", log)
        with open(self.bundle, "rb") as fh:
            self.assertEqual(fh.read(), before)
        rc, log = run(self.work, "--unpack_db", "--db", self.dbs["seekable"])
        self.assertEqual(rc, 30, log[-3000:])
        self.assertIn("--unpack_db needs a single-file database", log)

    def test_round_trip_is_byte_identical(self):
        """--decompress_db of the column-format index, and of the single file, gives exactly the raw files."""
        for kind in ("seekable", "bundle"):
            with self.subTest(kind):
                db = self.path(f"round_trip_{kind}_db")
                os.mkdir(db)
                for f in glob.glob(os.path.join(self.dbs[kind], "*")):
                    os.symlink(os.path.realpath(f), os.path.join(db, os.path.basename(f)))
                rc, log = run(self.work, "--decompress_db", "--db", db, "-t", "4")
                self.assertEqual(rc, 0, log[-3000:])
                self.assertFalse(os.path.lexists(os.path.join(db, "database.protal")))
                self.assertEqual(digest(os.path.join(db, "index.prx")), self.raw_digest, f"{kind}: index.prx")
                for name in self.FILES[1:]:
                    self.assertTrue(filecmp.cmp(os.path.join(db, name), os.path.join(self.dbs["raw"], name), shallow=False),
                                    f"{kind}: {name}")
                shutil.rmtree(db)  # a raw index.prx takes ~3 GB

    def test_unpack_db(self):
        """--unpack_db writes the files (index.prx.zst, reference.fna raw) and keeps database.protal."""
        out = self.path("unpacked")
        rc, log = run(self.work, "--unpack_db", "--db", self.bundle, "--unpack_dir", out, "-t", "2")
        self.assertEqual(rc, 0, log[-3000:])
        self.assertTrue(os.path.exists(self.bundle))
        for name in self.FILES:
            path = os.path.join(out, name + ".zst" if name == "index.prx" else name)
            self.assertTrue(os.path.isfile(path), path)
            if name != "index.prx":
                self.assertTrue(filecmp.cmp(path, os.path.join(self.dbs["raw"], name), shallow=False), name)
        self.assertTrue(is_seekable(os.path.join(out, "index.prx.zst")))

    def test_identical_results(self):
        # The raw index with one thread, against the baseline (the test database, 4 threads) and the other kinds (in
        # frames of 1 MB, loaded on several threads). The single file named directly (the baseline names its folder).
        base = baseline()
        with open_sam(base.sam("sa")) as sam:
            self.assertEqual(self.expected[0], sorted(line for line in sam if not line.startswith("@")))
        self.assertEqual(self.expected[1], base.text("sa.profile"))
        for kind, db, threads in (("seekable", self.dbs["seekable"], 4), ("single", self.dbs["single"], 2),
                                  ("bundle", self.bundle, 4)):
            with self.subTest(kind):
                sam, profile, log = self.result(db, f"out_{kind}_{threads}", threads)
                where = f"index.prx in {self.bundle}" if kind == "bundle" else os.path.join(self.dbs[kind], "index.prx")
                self.assertIn("Load index " + where, log)
                self.assertEqual(sam, self.expected[0], f"SAM differs: {kind} database, {threads} threads")
                self.assertEqual(profile, self.expected[1], f"profile differs: {kind} database, {threads} threads")

    def test_corrupt_seekable_index(self):
        bad_db = self.path("bad_seekable_db")
        os.mkdir(bad_db)
        for f in glob.glob(os.path.join(self.dbs["seekable"], "*")):
            if os.path.basename(f) != "index.prx.zst":
                os.symlink(os.path.realpath(f), os.path.join(bad_db, os.path.basename(f)))
        with open(os.path.join(self.dbs["seekable"], "index.prx.zst"), "rb") as fh:
            data = bytearray(fh.read())
        data[len(data) // 3] ^= 0x5A  # inside some frame; the seek table at the end is intact
        with open(os.path.join(bad_db, "index.prx.zst"), "wb") as fh:
            fh.write(data)
        rc, log = run(self.work, "--db", bad_db, *reads("sa"), "-o", "out_bad", "-t", "4", "--no_qcmsa")
        self.assertEqual(rc, 8, log[-3000:])
        self.assertRegex(log, r"Invalid index .*(frame|chunk) \d+ of \d+.*truncated or corrupt")

    def test_corrupt_single_file(self):
        with open(self.bundle, "rb") as fh:
            data = bytearray(fh.read())
        cut = self.path("cut.protal")
        with open(cut, "wb") as fh:
            fh.write(data[:-20])
        rc, log = run(self.work, "--db", cut, *reads("sa"), "-o", "out_cut", "-t", "4", "--no_qcmsa")
        self.assertEqual(rc, 30, log[-3000:])
        self.assertIn("seek table at its end is missing", log)
        self.assertNotIn("does not exist", log)
        data[len(data) // 3] ^= 0x5A  # inside the index's frames
        bad = self.path("bad.protal")
        with open(bad, "wb") as fh:
            fh.write(data)
        rc, log = run(self.work, "--db", bad, *reads("sa"), "-o", "out_bad_file", "-t", "4", "--no_qcmsa")
        self.assertEqual(rc, 8, log[-3000:])
        self.assertRegex(log, r"Invalid index index.prx in .*bad.protal: (frame|chunk) \d+ of \d+")

    def test_no_preload_needs_a_raw_reference(self):
        rc, log = run(self.work, "--db", self.dbs["seekable"], *reads("sa"), "-o", "out_lazy", "-t", "1", "--no_qcmsa",
                      "--preload_genomes_off")
        self.assertNotEqual(rc, 0)
        self.assertIn("--preload_genomes_off needs an uncompressed reference", log)

    def test_no_preload_on_a_single_file_says_how_to_unpack(self):
        """The message gives the commands that unpack the database and rerun protal; both work as given."""
        db = self.path("lazy_db")
        os.mkdir(db)
        os.symlink(os.path.realpath(self.bundle), os.path.join(db, "database.protal"))
        args = ["--db", db, *reads("sa"), "-o", "out_lazy_file", "-t", "2", "--no_qcmsa", "--preload_genomes_off"]
        rc, log = run(self.work, *args)
        self.assertEqual(rc, 30, log[-3000:])
        lines = log.splitlines()
        start = next(i for i, line in enumerate(lines) if "--preload_genomes_off reads genes one by one" in line)
        self.assertIn("is a single-file database", lines[start])
        unpack, rerun = lines[start + 1].strip(), lines[start + 3].strip()
        self.assertEqual(shlex.split(unpack), [PROTAL, "--unpack_db", "--db", os.path.join(db, "database.protal"), "-t", "2"])
        self.assertEqual(shlex.split(rerun), [PROTAL, *args])
        for command in (unpack, rerun):
            proc = subprocess.run(shlex.split(command), cwd=self.work, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                  text=True, errors="replace")
            self.assertEqual(proc.returncode, 0, proc.stdout[-3000:])
        self.assertIn(f"separate files in {db}", proc.stdout)
        self.assertNotIn("Preload genomes took", proc.stdout)
        with open_sam(find_sams(self.path("out_lazy_file", "sa.sam"))[0]) as sam:
            self.assertEqual(sorted(line for line in sam if not line.startswith("@")), self.expected[0])
        self.assertEqual(read_text(self.path("out_lazy_file", "sa.profile")), self.expected[1])


class ReadTypeModelTest(DbTest):
    """A database holds one presence model per read type (--read_type); --add_model stores one. Each test works on a
    copy of its own of the test database's single file (about 1 MB). The model added is the small gradient-boosted one
    of GradientBoostedModelTest: --add_model compresses what it stores at zstd level 19, which takes seconds for the
    mini database's 10 MB random forest."""

    MODEL = os.path.join(DATA, "model_gbm_small.xml")

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        # The test database's single file, or its files packed.
        cls.source = single_file(DB)
        if not cls.source:
            packed = symlink_db(os.path.join(cls.work, "source"))
            rc, log = run(cls.work, "--compress_db", "--db", packed, "--compress_level", "1", "-t", "4")
            if rc != 0:
                raise RuntimeError("protal --compress_db failed:\n" + log[-3000:])
            cls.source = os.path.join(packed, "database.protal")
        # sa's profile with the small model, of the paired-end reads and of the single-end ones (read_type_run's sa).
        cls.expected = {}
        for name, sam, model in (("pe", baseline().sam("sa"), "--model"), ("se", read_type_run().sam("sa"), "--model_se")):
            out = os.path.join(cls.work, f"expected_{name}")
            rc, log = profile_only(cls.work, out, [sam], model, cls.MODEL, prefixes=["sa"])
            if rc != 0:
                raise RuntimeError("protal --profile_only failed:\n" + log[-3000:])
            cls.expected[name] = read_text(os.path.join(out, "sa.profile"))

    def copy(self, name):
        """A copy of the test database's single file: (its folder, the file)."""
        os.mkdir(self.path(name))
        shutil.copy(self.source, self.path(name, "database.protal"))
        return self.path(name), self.path(name, "database.protal")

    def profile_only(self, db, out, *extra, sam=None):
        return profile_only(self.work, self.path(out), [sam or baseline().sam("sa")], *extra, prefixes=["sa"], db=db)

    def test_add_model_for_a_read_type(self):
        db, bundle = self.copy("db_add")
        rc, log = self.profile_only(db, "out_pe")
        self.assertEqual(rc, 0, log[-3000:])
        self.assertIn("Model of paired-end reads: model_pe.xml in " + bundle, log)
        self.assertEqual(read_text(self.path("out_pe", "sa.profile")), baseline().text("sa.profile"))
        rc, log = self.profile_only(db, "out_pb", "--read_type", "pb")
        self.assertEqual(rc, 30, log[-3000:])
        self.assertIn("no model for --read_type pb (PacBio reads): model_PB.xml in", log)
        self.assertIn("--add_model MODEL.xml --read_type pb --db " + bundle, log)

        # A model stored as the single-end one: paired-end reads profiled as single-end ones with it, as with the model
        # given as a file to paired-end reads.
        rc, log = run(self.work, "--add_model", self.MODEL, "--read_type", "se", "--db", db, "-t", "2")
        self.assertEqual(rc, 0, log[-3000:])
        self.assertIn(f"Stored {self.MODEL} as model_se.xml in {bundle}", log)
        self.assertRegex(log, r"Models for read types: pe, se\b")
        rc, log = self.profile_only(db, "out_se", "--read_type", "se")
        self.assertEqual(rc, 0, log[-3000:])
        self.assertIn("Model of single-end reads: model_se.xml in " + bundle, log)
        self.assertIn("holds paired-end reads; profiled as single-end reads (se, --read_type or READ_TYPE)", log)
        self.assertEqual(read_text(self.path("out_se", "sa.profile")), self.expected["pe"])
        rc, log = run(self.work, "--unpack_db", "--db", bundle, "--unpack_dir", self.path("unpacked"))
        self.assertEqual(rc, 0, log[-3000:])
        self.assertTrue(filecmp.cmp(self.path("unpacked", "model_se.xml"), self.MODEL, shallow=False))
        # Single-end reads take it from the single file, as from a file given.
        rc, log = self.profile_only(bundle, "out_single_end", sam=read_type_run().sam("sa"))
        self.assertEqual(rc, 0, log[-3000:])
        self.assertIn("Model of single-end reads: model_se.xml in " + bundle, log)
        self.assertEqual(read_text(self.path("out_single_end", "sa.profile")), self.expected["se"])

    def test_placeholder_model(self):
        # scripts/placeholder_models.py fills a read type's slot until a trained model replaces it: it reports
        # no species (--knob 0: every taxon with reads), and protal warns whenever it loads it.
        db, bundle = self.copy("db_ont")
        python(os.path.join(ROOT, "scripts", "placeholder_models.py"), "-o", self.path("placeholders"), "--read_types", "ont")
        placeholder = self.path("placeholders", "model_ONT.xml")
        warning = "is a placeholder, not a trained model"
        rc, log = run(self.work, "--add_model", placeholder, "--read_type", "ont", "--db", db, "-t", "2")
        self.assertEqual(rc, 0, log[-3000:])
        self.assertIn(warning, log)
        self.assertRegex(log, r"Models for read types: [^\n]*\bont\b")
        rc, log = self.profile_only(db, "out_ont", "--read_type", "ont")
        self.assertEqual(rc, 0, log[-3000:])
        self.assertIn("Model of ONT reads: model_ONT.xml in " + bundle, log)
        self.assertIn(warning, log)
        self.assertEqual(read_text(self.path("out_ont", "sa.profile")), "", "a placeholder reports no species")
        rc, log = self.profile_only(db, "out_ont_all", "--read_type", "ont", "--knob", "0")
        self.assertEqual(rc, 0, log[-3000:])
        self.assertEqual(profile_species(read_text(self.path("out_ont_all", "sa.profile"))), sorted(species()))
        rc, log = self.profile_only(db, "out_pe_again")
        self.assertEqual(rc, 0, log[-3000:])
        self.assertNotIn("placeholder", log, "the paired-end model is not one")

    def test_several_models_at_once(self):
        # --add_model A,B --read_type se,pb: every model checked first, the database rewritten once.
        db, bundle = self.copy("db_several")
        model = self.MODEL
        rc, log = run(self.work, "--add_model", f"{model},{model}", "--read_type", "se,pb", "--db", db, "-t", "2")
        self.assertEqual(rc, 0, log[-3000:])
        self.assertIn(f"Stored {model} as model_se.xml in {bundle}", log)
        self.assertIn(f"Stored {model} as model_PB.xml in {bundle}", log)
        self.assertEqual(log.count("Models for read types:"), 1, "one rewrite")
        self.assertRegex(log, r"Models for read types: [^\n]*\bse\b[^\n]*\bpb\b")
        rc, log = run(self.work, "--unpack_db", "--db", bundle, "--unpack_dir", self.path("unpacked_two"))
        self.assertEqual(rc, 0, log[-3000:])
        for member in ("model_se.xml", "model_PB.xml"):
            self.assertTrue(filecmp.cmp(self.path("unpacked_two", member), model, shallow=False), member)
        # Their read types one each, and every model usable, else nothing is written.
        with open(bundle, "rb") as fh:
            before = fh.read()
        bad = self.path("bad_second.xml")
        with open(bad, "w") as fh:
            fh.write("<PMML>\n")
        for models, types, code, problem in ((f"{model},{model}", "se", 30, "--read_type must give as many"),
                                             (f"{model},{model}", "se,se", 30, "--read_type names se twice"),
                                             (model, "se,pb", 30, "--read_type gives 2 read types for --add_model's 1 model"),
                                             (f"{model},{bad}", "se,ont", 2, "Cannot load the model " + bad),
                                             (bad, "ont", 2, "Cannot load the model")):
            with self.subTest(types=types, models=models):
                rc, log = run(self.work, "--add_model", models, "--read_type", types, "--db", db)
                self.assertEqual(rc, code, log[-3000:])
                self.assertIn(problem, log)
                with open(bundle, "rb") as fh:
                    self.assertEqual(fh.read(), before, "the database is unchanged")

    def test_models_are_replaced_in_place(self):
        # The models are the last members of database.protal (--build puts them there): --add_model replaces one there
        # in place, the members before it not rewritten. The same model again, at the level the database was packed at
        # (--build's default), gives the same bytes.
        db, bundle = self.copy("db_in_place")
        rc, log = run(self.work, "--add_model", self.MODEL, "--read_type", "pe", "--db", db, "-t", "2")
        self.assertEqual(rc, 0, log[-3000:])
        self.assertRegex(log, r"Replace the last \d+ of \d+ members of \S+ in place")
        self.assertFalse(os.path.exists(bundle + ".journal"))
        with open(bundle, "rb") as fh:
            replaced = fh.read()
        rc, log = run(self.work, "--add_model", self.MODEL, "--read_type", "pe", "--db", db, "-t", "2")
        self.assertEqual(rc, 0, log[-3000:])
        self.assertRegex(log, r"Replace the last \d+ of \d+ members of \S+ in place")
        with open(bundle, "rb") as fh:
            self.assertEqual(fh.read(), replaced)
        rc, log = self.profile_only(db, "out_pe_in_place")
        self.assertEqual(rc, 0, log[-3000:])
        self.assertIn("Model of paired-end reads: model_pe.xml in " + bundle, log)
        self.assertEqual(read_text(self.path("out_pe_in_place", "sa.profile")), self.expected["pe"])

    def test_read_type_checks(self):
        rc, log = run(self.work, "--db", DB, *reads("sa"), "-o", "out_x", "--read_type", "nanopore")
        self.assertEqual(rc, 30, log[-3000:])
        self.assertIn("is 'nanopore' (--read_type or READ_TYPE): give one of pe, se, pb, ont", log)
        rc, log = run(self.work, "--db", DB, *reads("sa"), "-o", "out_y", "--read_type", "se")
        self.assertEqual(rc, 30, log[-3000:])
        self.assertIn("has single-end reads (se), which come in one file, but it has a second read file", log)

    def test_older_database_with_model_xml(self):
        """model.xml of a database from before read types serves paired-end reads."""
        db = symlink_db(self.path("old_db"), skip=("model_pe.xml", "model_se.xml", "model_PB.xml", "model_ONT.xml"),
                        extra={"model.xml": db_file("model_pe.xml")})
        rc, log = profile_only(self.work, self.path("out_old"), [baseline().sam("sa")], prefixes=["sa"], db=db)
        self.assertEqual(rc, 0, log[-3000:])
        self.assertIn("Model of paired-end reads: " + os.path.join(db, "model.xml"), log)
        self.assertEqual(read_text(self.path("out_old", "sa.profile")), baseline().text("sa.profile"))


class FailFastTest(DbTest):
    """Problems with the database or the inputs stop protal before any read is aligned: the exit code, that no SAM was
    written, and the message where no unit test checks it."""

    def db_copy(self, name, replace=None, drop=()):
        """A database of symlinks to DB, with the files in `replace` ({name: bytes}) written instead
        and the files in `drop` left out."""
        replace = replace or {}
        db = symlink_db(self.path(name), skip=set(replace) | set(drop))
        for base, content in replace.items():
            with open(os.path.join(db, base), "wb") as fh:
                fh.write(content)
        return db

    def query(self, db, out, *extra, samples=("sa",)):
        rc, log = run(self.work, "--db", db, *reads(*samples), "-o", out, "-t", "1", "--no_qcmsa", *extra)
        self.assertFalse(glob.glob(self.path(out, "*.sam*")), "no read may be aligned")
        return rc, log

    def test_reference_changed_since_build(self):
        db = self.db_copy("db_ref", {"reference.fna": db_content("reference.fna") + b">9_1\nACGT\n"})
        rc, log = self.query(db, "out_ref")
        if "no reference fingerprint" in log:
            raise unavailable("the test database's index predates the reference fingerprint")
        self.assertEqual(rc, 8, log[-3000:])
        self.assertIn("built against a different reference", log)

    def test_malformed_map(self):
        db = self.db_copy("db_map", {"reference.map": db_content("reference.map") + b"1\t999\t5\n"})
        rc, log = self.query(db, "out_map")
        self.assertEqual(rc, 8, log[-3000:])
        self.assertIn("Invalid reference map", log)

    def test_missing_db(self):
        """A --db path with nothing at it is reported as missing, relative to where protal runs."""
        for args in (reads("sa") + ["-o", "out_nodb", "--no_qcmsa"], ["--unpack_db"], ["--compress_db"]):
            with self.subTest(args[0]):
                rc, log = run(self.work, "--db", "no/such_db", *args)
                self.assertEqual(rc, 30, log[-3000:])
                self.assertIn("--db no/such_db does not exist (relative to the working directory "
                              f"{os.path.realpath(self.work)})", log)
                self.assertNotIn("holds separate files", log)
                self.assertNotIn("Sequence file does not exist", log)

    def test_malformed_gene_conservation(self):
        # Read by every query: the factors give the model's conservation features.
        db = self.db_copy("db_conservation", {"gene_conservation.tsv": b"geneid\tfactor\tspecies\n1\tfast\t3\n"})
        rc, log = self.query(db, "out_conservation")
        self.assertEqual(rc, 8, log[-3000:])
        self.assertIn("Invalid gene conservation factors", log)

    def test_missing_model_and_unique_kmers(self):
        db = self.db_copy("db_files", drop=("model_pe.xml", "unique_kmers.tsv"))
        rc, log = self.query(db, "out_files")
        self.assertEqual(rc, 30, log[-3000:])
        self.assertIn("The database has no model for --read_type pe", log)
        self.assertIn("Unique k-mer file does not exist", log)

    def test_unusable_models(self):
        model = db_content("model_pe.xml").decode()
        # Not a model; an input protal does not compute; a model predicting other labels than TRUE/FALSE.
        unknown = model.replace("<MiningSchema>", '<MiningSchema>\n<MiningField name="moon_phase"/>', 1)
        unknown = re.sub(r"(<DataDictionary[^>]*>)", r'\1\n<DataField name="moon_phase" optype="continuous" dataType="double"/>',
                         unknown, count=1)
        labels = model.replace('value="TRUE"', 'value="present"').replace('score="TRUE"', 'score="present"')
        # The last two problems' messages: ModelFeatures.TheModelMustFitProtal.
        for name, text, message in (("db_model", "<PMML>\n", "Cannot load the model"), ("db_unknown", unknown, "moon_phase"),
                                    ("db_labels", labels, None)):
            with self.subTest(name):
                db = self.db_copy(name, {"model_pe.xml": text.encode()})
                rc, log = self.query(db, "out_" + name)
                self.assertEqual(rc, 2, log[-3000:])
                if message:
                    self.assertIn(message, log)

    def test_a_read_type_without_a_model_stops_before_aligning(self):
        # The test database has a model for paired-end reads only.
        rc, log = run(self.work, "--db", DB, *single_reads("sa"), "-o", "out_se", "-t", "1", "--no_qcmsa")
        self.assertEqual(rc, 30, log[-3000:])
        self.assertIn("no model for --read_type se (single-end reads)", log)
        self.assertFalse(glob.glob(self.path("out_se", "*.sam*")), "no read may be aligned")

    def test_knob_is_a_probability(self):
        rc, log = self.query(DB, "out_knob", "--knob", "1.5")
        self.assertNotEqual(rc, 0, log[-3000:])
        self.assertIn("--knob must be between 0 and 1", log)

    def test_one_truth_file_per_sample(self):
        truth = self.path("truth.tsv")
        with open(truth, "w") as fh:
            fh.write("s__Mockella alpha\n")
        rc, log = self.query(DB, "out_truth", "--profile_truth", truth, samples=("sa", "sb"))
        self.assertEqual(rc, 30, log[-3000:])
        self.assertIn("must name one file per sample", log)

    def test_map_row_without_a_profile_cell(self):
        sample_map = self.path("samples.map")
        with open(sample_map, "w") as fh:
            fh.write(f"#OUTPUT_DIR\t{self.path('out_map_rows')}\n#SAMPLEID\tPREFIX\tFIRST\tSECOND\tPROFILE\n")
            fh.write(f"sa\tsa\t{READS}/sa_R1.fq\t{READS}/sa_R2.fq\tsa.profile\n")
            fh.write(f"sb\tsb\t{READS}/sb_R1.fq\t{READS}/sb_R2.fq\n")
        rc, log = run(self.work, "--db", DB, "--map", sample_map, "-t", "1", "--no_qcmsa")
        self.assertEqual(rc, 9, log[-3000:])  # the message: SampleMap.RejectsRowsWithMissingOrEmptyCells
        self.assertFalse(glob.glob(self.path("out_map_rows", "**", "*.sam*"), recursive=True))

    def test_benchmark_needs_reads_named_by_gene(self):
        rc, log = self.query(DB, "out_bench", "--benchmark_alignment")
        self.assertEqual(rc, 2, log[-3000:])
        self.assertIn("--benchmark_alignment needs reads named <taxid>_<gene id>", log)

    def test_build_rejects_a_reference_the_map_does_not_describe(self):
        db = self.path("build_db")
        os.mkdir(db)
        for f in ("reference.fna", "reference.map", "internal_taxonomy.dmp"):
            with open(os.path.join(db, f), "wb") as fh:
                fh.write(db_content(f))
        bad = self.path("bad_reference.fna")
        with open(bad, "wb") as fh:
            fh.write(db_content("reference.fna") + b">9_1\nACGTACGT\n>unnamed\nACGTACGT\n")
        rc, log = run(self.work, "--build", "--no_profile", "-t", "1", "--db", db,
                      "--reference", bad, "--full_reference", bad)
        self.assertEqual(rc, 8, log[-3000:])
        self.assertIn("--reference record >9_1: not in", log)
        self.assertIn("--reference record >unnamed: header is not <taxid>_<gene id>", log)
        self.assertRegex(log, r"2 of \d+ --reference records do not match")
        for written in ("index.prx", "index.prx.zst", "database.protal"):
            self.assertFalse(os.path.exists(os.path.join(db, written)), written)


class SamInputTest(DbTest):
    """--profile_only reads SAM files as other tools may leave them (the baseline's SAM of sa, edited)."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        base = baseline()
        lines = sam_text(base.sam("sa")).splitlines()
        cls.header = [line for line in lines if line.startswith("@")]
        cls.records = [line for line in lines if not line.startswith("@")]
        cls.profile_text = base.text("sa.profile")
        if not cls.records or not cls.profile_text.strip():
            raise RuntimeError("the baseline's sa has no records or reports no taxon")

    def write_sam(self, name, lines, final_newline=True):
        sam = self.path(f"{name}.sam")
        with open(sam, "w", newline="") as fh:
            fh.write("\n".join(lines) + ("\n" if final_newline else ""))
        return sam

    def profile_only(self, *sams):
        return run(self.work, "--db", DB, "--profile_only", ",".join(sams), "-o", self.path("out_" + os.path.basename(sams[0])),
                   "-t", "1", "--no_qcmsa")

    def test_edited_sam_profiles_like_the_original(self):
        half = len(self.records) // 2
        unmapped = "x.1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII"
        foreign = "\t".join(["y.1", "0", "chr1"] + self.records[0].split("\t")[3:])
        edited = (self.header + self.records[:half] + ["", unmapped, foreign] +
                  [r + "\r" for r in self.records[half:half + 10]] + self.records[half + 10:])
        sam = self.write_sam("edited", edited, final_newline=False)
        rc, log = self.profile_only(sam)
        self.assertEqual(rc, 0, log[-3000:])
        # The added unmapped record, and those protal wrote itself for reads that seeded on taxa but aligned nowhere.
        # (SamReader.SkipsHeadersBlankLinesAndUnusableRecords counts the other kinds of skipped records.)
        own_unmapped = sum(1 for r in self.records if int(r.split("\t")[1]) & 0x4)
        self.assertIn(f"skipped {own_unmapped + 1} record(s): unmapped", log)
        self.assertEqual(read_text(self.path("out_edited.sam", "edited.profile")), self.profile_text)

    def test_header_only_sam_gets_an_empty_profile(self):
        sam = self.write_sam("header_only", self.header)
        # A sample without records that do not fit gets no misc/<sample>.err, and an earlier run's is removed.
        os.makedirs(self.path("out_header_only.sam", "misc"))
        with open(self.path("out_header_only.sam", "misc", "header_only.err"), "w") as fh:
            fh.write("an earlier run's\n")
        rc, log = self.profile_only(sam)
        self.assertEqual(rc, 0, log[-3000:])
        self.assertEqual(read_text(self.path("out_header_only.sam", "header_only.profile")), "")
        self.assertFalse(os.path.exists(self.path("out_header_only.sam", "misc", "header_only.err")))

    def test_records_that_do_not_fit_go_to_misc(self):
        # A record whose bases differ from its gene under an M (a SAM aligned against another database) is left out and
        # listed in misc/<sample>.err of the run, not beside the SAM, which may be in another run's or a read-only folder.
        i = next(k for k, r in enumerate(self.records)
                 if not int(r.split("\t")[1]) & 0x4 and re.match(r"\d+M", r.split("\t")[5]) and r.split("\t")[9][0] in "ACGT")
        fields = self.records[i].split("\t")
        fields[9] = ("C" if fields[9][0] != "C" else "G") + fields[9][1:]
        sam = self.write_sam("misfit", self.header + self.records[:i] + ["\t".join(fields)] + self.records[i + 1:])
        rc, log = self.profile_only(sam)
        self.assertEqual(rc, 0, log[-3000:])
        err = self.path("out_misfit.sam", "misc", "misfit.err")
        self.assertIn("reads have an alignment that does not fit the database", log)
        self.assertIn("listed in " + err, log)
        self.assertIn(fields[0], read_text(err))
        self.assertFalse(os.path.exists(sam + ".err"))

    def test_unreadable_sam_fails_only_its_sample(self):
        good = self.write_sam("good", self.header + self.records)
        broken = self.write_sam("broken", self.header + self.records[:50] + ["sa.9\t0\t1_1"])
        rc, log = self.profile_only(good, broken)
        self.assertEqual(rc, 1, log[-3000:])
        self.assertRegex(log, r"Cannot read the SAM file of sample \S+ \(.*broken\.sam\)")
        self.assertEqual(read_text(self.path("out_good.sam", "good.profile")), self.profile_text)

    def test_truncated_gzip_sam_fails_its_sample(self):
        sam = self.write_sam("cut", self.header + self.records)
        with open(sam, "rb") as fh, gzip.open(sam + ".gz", "wb") as gz:
            gz.write(fh.read())
        os.remove(sam)
        size = os.path.getsize(sam + ".gz")
        with open(sam + ".gz", "r+b") as fh:
            fh.truncate(size // 2)
        rc, log = self.profile_only(sam + ".gz")
        self.assertEqual(rc, 1, log[-3000:])
        self.assertRegex(log, r"Cannot read the SAM file of sample cut \(.*cut\.sam\.gz\): the file is truncated or corrupt")

    def test_sam_of_another_database_fails_its_sample(self):
        sq = next(i for i, line in enumerate(self.header) if line.startswith("@SQ"))
        name, length = re.match(r"@SQ\tSN:(\S+)\tLN:(\d+)", self.header[sq]).groups()
        header = list(self.header)
        header[sq] = f"@SQ\tSN:{name}\tLN:{int(length) + 7}"
        sam = self.write_sam("other_db", header + self.records)
        rc, log = self.profile_only(sam)
        self.assertEqual(rc, 1, log[-3000:])  # the message: FromSam.RejectsASamAlignedAgainstAnotherDatabase
        self.assertRegex(log, r"Cannot read the SAM file of sample \S+ \(.*other_db\.sam\)")

    def test_samples_cannot_share_an_output_file(self):
        sam = self.write_sam("twice", self.header + self.records)
        rc, log = self.profile_only(sam, sam)
        self.assertNotEqual(rc, 0, log[-3000:])
        self.assertIn("samples 1 and 2 would both use the SAM file", log)
        self.assertFalse(glob.glob(self.path("out_twice.sam", "*.profile")))


class CompressedSamOutputTest(DbTest):
    """A map's SAM names choose the format: .sam, .sam.gz (BGZF, compressed while aligning) or
    .sam.zst (seekable zstd). Headers list the genes that records name; --full_sam_header all. (A part of sa's reads.)"""

    FORMATS = {"plain": "sam", "gz": "sam.gz", "zst": "sam.zst"}

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        for mate in (1, 2):
            head_reads(os.path.join(READS, f"sa_R{mate}.fq"), os.path.join(cls.work, f"part_R{mate}.fq"), 1500)
        cls.sample_map = os.path.join(cls.work, "samples.map")
        with open(cls.sample_map, "w") as fh:
            fh.write(f"#OUTPUT_DIR\t{os.path.join(cls.work, 'out')}\n#INPUT_DIR\t{cls.work}\n")
            fh.write("#SAMPLEID\tPREFIX\tFIRST\tSECOND\tSAM\n")
            # Names of their own: a rerun takes an existing x.sam for x.sam.gz or x.sam.zst.
            for name, ext in cls.FORMATS.items():
                fh.write(f"s_{name}\ts_{name}\tpart_R1.fq\tpart_R2.fq\ts_{name}.{ext}\n")
        cls.rc, cls.log = run(cls.work, "--db", DB, "--map", cls.sample_map, "-t", "1", "--no_qcmsa", "--no_strains")

    def sam(self, name):
        return self.path("out", "alignments", f"s_{name}.{self.FORMATS[name]}")

    def text(self, name):
        return sam_text(self.sam(name))

    def test_exit_code(self):
        self.assertEqual(self.rc, 0, self.log[-3000:])
        self.assertNotIn("pigz", self.log)

    def test_every_format_holds_the_same_sam(self):
        plain = self.text("plain")
        self.assertTrue(any(not line.startswith("@") for line in plain.splitlines()))
        for name in ("gz", "zst"):
            self.assertEqual(self.text(name), plain, name)  # one thread: the same records in the same order
        self.assertEqual(glob.glob(self.path("out", "**", "*partial*"), recursive=True), [], "no temporary files left")

    def test_gzip_is_bgzf_with_its_end_of_file_block(self):
        with open(self.sam("gz"), "rb") as fh:
            data = fh.read()
        self.assertEqual(data[:4], b"\x1f\x8b\x08\x04")
        self.assertEqual(data[12:16], b"BC\x02\x00")
        self.assertEqual(data[-28:], bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000"))

    def test_the_header_lists_the_genes_the_records_name(self):
        lines = self.text("plain").splitlines()
        listed = [line.split("\t")[1][3:] for line in lines if line.startswith("@SQ")]
        named = set()
        for line in lines:
            if not line.startswith("@"):
                fields = line.split("\t")
                if fields[2] == "*":  # an unmapped record (a read that seeded on taxa but aligned nowhere): no gene
                    continue
                named.add(fields[2])
                if fields[6] not in ("=", "*"):
                    named.add(fields[6])
        self.assertEqual(sorted(listed), sorted(named))
        self.assertEqual(len(listed), len(set(listed)))
        self.assertTrue(set(listed) <= {name for name, _ in reference_genes()})

    def test_profiles_do_not_depend_on_the_format(self):
        profiles = {}
        for name in self.FORMATS:
            profiles[name] = read_text(glob.glob(self.path("out", "**", f"s_{name}.profile"), recursive=True)[0])
        self.assertTrue(profiles["plain"].strip())
        self.assertEqual(profiles["gz"], profiles["plain"])
        self.assertEqual(profiles["zst"], profiles["plain"])

    def test_a_rerun_leaves_the_compressed_sams_alone(self):
        before = {}
        for name in self.FORMATS:
            with open(self.sam(name), "rb") as fh:
                before[name] = (fh.read(), os.path.getmtime(self.sam(name)))
        rc, log = run(self.work, "--db", DB, "--map", self.sample_map, "-t", "1", "--no_qcmsa", "--no_strains")
        self.assertEqual(rc, 0, log[-3000:])
        self.assertIn("All alignments are present", log)
        for name in self.FORMATS:
            with open(self.sam(name), "rb") as fh:
                self.assertEqual((fh.read(), os.path.getmtime(self.sam(name))), before[name], name)

    def test_sam_format_names_the_sams_protal_picks(self):
        # The baseline's -1/-2 run shows the default, .sam.zst; --sam_format picks another. The plain one with
        # --full_sam_header, which lists every gene.
        reads_part = ["-1", self.path("part_R1.fq"), "-2", self.path("part_R2.fq"), "--prefix", "sa"]
        for fmt, name, extra in (("gz", "sa.sam.gz", []), ("sam", "sa.sam", ["--full_sam_header"])):
            rc, log = run(self.work, "--db", DB, *reads_part, "-o", f"out_{fmt}", "-t", "1", "--no_qcmsa", "--no_profile",
                          "--sam_format", fmt, *extra)
            self.assertEqual(rc, 0, log[-3000:])
            self.assertEqual(os.listdir(self.path(f"out_{fmt}")).count(name), 1, os.listdir(self.path(f"out_{fmt}")))
        with open(self.path("out_sam", "sa.sam")) as fh:
            full = [line for line in fh if line.startswith("@SQ")]
        self.assertEqual(len(full), len(reference_genes()))
        listed = [line for line in self.text("plain").splitlines() if line.startswith("@SQ")]
        self.assertTrue(set(line + "\n" for line in listed) <= set(full))
        rc, log = run(self.work, "--db", DB, *reads_part, "-o", "out_bad", "--sam_format", "bam")
        self.assertEqual(rc, 2, log[-3000:])
        self.assertIn("--sam_format must be zst, gz or sam", log)


class QcmsaTest(DbTest):
    """The post-filter runs on the strain MSAs (of the baseline's SAMs), and --qcmsa_args reaches it intact."""

    def test_filtered_msa(self):
        # With its defaults, qcmsa filters an MSA of three samples (a gene needs two).
        samples = ("sa", "sb", "sr")
        baseline().place_sams(self.path("out"), samples)
        rc, log = run(self.work, "--db", DB, *reads(*samples), "-o", "out", "-t", "4", "--qcmsa_script", QCMSA)
        self.assertEqual(rc, 0, log[-3000:])
        filtered = [f for f in glob.glob(self.path("out", "strains", "*.msa.fna")) if not f.endswith(".raw.msa.fna")]
        self.assertTrue(filtered, "qcmsa wrote filtered MSAs")
        with open(self.path("out", "strains", "species.tsv")) as fh:
            listed = [line.rstrip("\n").split("\t") for line in fh][1:]
        self.assertEqual(sorted(row[4] for row in listed if row[4] != "-"), sorted(os.path.basename(f) for f in filtered))

        # A rerun in which qcmsa keeps nothing reports it and leaves no filtered MSA of the first run behind.
        rc, log = run(self.work, "--db", DB, *reads(*samples), "-o", "out", "-t", "4",
                      "--qcmsa_script", QCMSA, "--qcmsa_args", "--gene-min-samples 100")
        self.assertEqual(rc, 0, log[-3000:])
        self.assertIn("qcmsa kept no gene or sample", log)
        self.assertFalse([f for f in filtered if os.path.exists(f)], "stale filtered MSAs")
        self.assertEqual([f for f in glob.glob(self.path("out", "strains", "*.msa.fna")) if not f.endswith(".raw.msa.fna")], [])
        with open(self.path("out", "strains", "species.tsv")) as fh:
            self.assertTrue(all(line.rstrip("\n").split("\t")[4] == "-" for line in list(fh)[1:]))

    def test_sample_ids_must_be_unique(self):
        sample_map = self.path("dup.map")
        with open(sample_map, "w") as fh:
            fh.write(f"#OUTPUT_DIR\t{self.path('out_dup')}\n#SAMPLEID\tPREFIX\tFIRST\tSECOND\n")
            fh.write(f"same\tsa\t{READS}/sa_R1.fq\t{READS}/sa_R2.fq\n")
            fh.write(f"same\tsb\t{READS}/sb_R1.fq\t{READS}/sb_R2.fq\n")
        rc, log = run(self.work, "--db", DB, "--map", sample_map, "-t", "1", "--no_qcmsa")
        self.assertNotEqual(rc, 0, log[-3000:])
        self.assertIn("share the sample ID 'same'", log)
        self.assertFalse(glob.glob(self.path("out_dup", "**", "*.sam*"), recursive=True))


def simulate_long_reads(path, reads, seed, genes_per_read=(3, 8), long_read=0, error=0.001, indels=0.6,
                        read_name="m64001_000000/{}/ccs", quality="I"):
    """Write PacBio-like reads to path: reference genes (either strand) between random stretches of
    0.5-3 kb, with 0.1% errors (substitutions and 1 bp indels, indels 60% of them), 3-8 genes per
    read; half of the reads start inside a gene and half end inside one (30-70% of it on the read).
    With long_read, one more read of that length with a gene every 5 kb. error, indels, read_name and
    quality make ONT-like reads. Returns per read its genes as (name, start, end) on the read."""
    rng = random.Random(seed)
    genes = [(name, seq) for name, seq in reference_genes() if len(seq) >= 300]
    substitution, insertion = error * (1 - indels), error * (1 - indels / 2)

    def mutate(seq):
        out = []
        for b in seq:
            r = rng.random()
            if r < substitution:
                out.append(rng.choice([c for c in "ACGT" if c != b]))
            elif r < insertion:
                out.append(b + rng.choice("ACGT"))
            elif r >= error:
                out.append(b)
        return "".join(out)

    def spacer(n):
        return "".join(rng.choice("ACGT") for _ in range(n))

    def read_of(pieces):
        seq, placed = "", []
        for kind, name, piece in pieces:
            piece = mutate(piece)
            if kind == "gene":
                placed.append((name, len(seq), len(seq) + len(piece)))
            seq += piece
        return seq, placed

    truth = []
    with open(path, "w") as fh:
        layouts = []

        def oriented_gene():
            name, seq = rng.choice(genes)
            return name, seq if rng.random() < 0.5 else revcomp(seq)

        for _ in range(reads):
            pieces = [("spacer", None, spacer(rng.randint(500, 3000)))]
            for _ in range(rng.randint(*genes_per_read)):
                pieces.append(("gene", *oriented_gene()))
                pieces.append(("spacer", None, spacer(rng.randint(500, 3000))))
            if rng.random() < 0.5:  # the read starts inside a gene
                name, seq = oriented_gene()
                pieces[0] = ("gene", name, seq[-int(len(seq) * rng.uniform(0.3, 0.7)):])
            if rng.random() < 0.5:  # ... and ends inside one
                name, seq = oriented_gene()
                pieces[-1] = ("gene", name, seq[:int(len(seq) * rng.uniform(0.3, 0.7))])
            layouts.append(pieces)
        if long_read:
            pieces, length = [], 0
            while length < long_read:
                name, seq = rng.choice(genes)
                gap = max(0, 5000 - len(seq))
                pieces += [("spacer", None, spacer(gap)), ("gene", name, seq if rng.random() < 0.5 else revcomp(seq))]
                length += gap + len(seq)
            layouts.append(pieces)
        for i, pieces in enumerate(layouts, 1):
            seq, placed = read_of(pieces)
            truth.append(placed)
            fh.write(f"@{read_name.format(i)}\n{seq}\n+\n{quality * len(seq)}\n")
    return truth


def representative_records(records):
    """Per read name, its records that are not secondary (the primary and the supplementary ones)."""
    by_read = {}
    for r in records:
        if not int(r[1]) & 0x100:
            by_read.setdefault(r[0], []).append(r)
    return by_read


def cigar_ops(cigar):
    return [(int(n), op) for n, op in re.findall(r"(\d+)([MIDNSHPX=])", cigar)]


def fastq_records(path):
    """[(name, sequence)] of a FASTQ file."""
    with open(path) as fh:
        lines = fh.read().splitlines()
    return [(lines[i][1:], lines[i + 1]) for i in range(0, len(lines), 4)]


def find_one(folder, name):
    """The one file `name` below folder."""
    matches = glob.glob(os.path.join(glob.escape(folder), "**", name), recursive=True)
    if len(matches) != 1:
        raise AssertionError(f"expected one {name} in {folder}: {matches}")
    return matches[0]


_read_type_db = None


def read_type_db():
    """A database of symlinks to the test database's files in which its paired-end model also serves single-end,
    PacBio and ONT reads (model_se.xml, model_PB.xml, model_ONT.xml): the test database has none, so the profiles of
    these read types test the plumbing, not a model."""
    global _read_type_db
    if _read_type_db is None:
        require_database()
        model = db_file("model_pe.xml")
        _read_type_db = symlink_db(os.path.join(scratch("protal_e2e_read_type_db_"), "db"),
                                   extra={"model_se.xml": model, "model_PB.xml": model, "model_ONT.xml": model})
    return _read_type_db


class ReadTypeRun:
    """One map run (strain MSAs on, --no_qcmsa) of samples of every read type on read_type_db(): pe (sa's pairs), sa and
    sb (their first mates, single-end; sa's read type found from its reads, sb's given), la and lb (PacBio-like reads of
    several genes each, la with a read of 150 kb), oa and ob (ONT-like: 2% errors, most of them 1 bp indels, Q17; oa
    with a read of 150 kb)."""

    def __init__(self):
        self.db = read_type_db()
        self.work = scratch("protal_e2e_read_types_")
        self.out = os.path.join(self.work, "out")
        self.long = {p: os.path.join(self.work, f"{p}.fq") for p in ("la", "lb", "oa", "ob")}
        ont = dict(error=0.02, indels=0.7, read_name="ont_read_{}", quality="2")
        self.truth = {"la": simulate_long_reads(self.long["la"], 40, seed=11, long_read=150000),
                      "lb": simulate_long_reads(self.long["lb"], 40, seed=12),
                      "oa": simulate_long_reads(self.long["oa"], 40, seed=21, long_read=150000, **ont),
                      "ob": simulate_long_reads(self.long["ob"], 40, seed=22, **ont)}
        rows = [("pe", os.path.join(READS, "sa_R1.fq"), os.path.join(READS, "sa_R2.fq"), "-"),
                ("sa", os.path.join(READS, "sa_R1.fq"), "-", "-"), ("sb", os.path.join(READS, "sb_R1.fq"), "-", "se"),
                ("la", self.long["la"], "-", "PB"), ("lb", self.long["lb"], "-", "pb"),
                ("oa", self.long["oa"], "-", "ont"), ("ob", self.long["ob"], "-", "ONT")]
        sample_map = os.path.join(self.work, "samples.map")
        with open(sample_map, "w") as fh:
            fh.write(f"#OUTPUT_DIR\t{self.out}\n#SAM_OUTPUT_DIR\t.\n#SAMPLEID\tPREFIX\tFIRST\tSECOND\tREAD_TYPE\n")
            fh.writelines(f"{name}\t{name}\t{first}\t{second}\t{read_type}\n" for name, first, second, read_type in rows)
        self.rc, self.log = run(self.work, "--db", self.db, "--map", sample_map, "-t", "4", "--no_qcmsa")
        if self.rc != 0:
            raise RuntimeError(f"the read types' run exited {self.rc}:\n" + self.log[-3000:])

    def sam(self, sample):
        return find_one(self.out, f"{sample}.sam.zst")

    def text(self, name):
        return read_text(find_one(self.out, name))

    def msa_rows(self, species_name):
        """The row names of a species' strain MSA."""
        with open(find_one(self.out, species_name.replace(" ", "_") + ".raw.msa.fna")) as fh:
            return [line[1:].strip() for line in fh if line.startswith(">")]


_read_type_run = None


def read_type_run():
    """The module's ReadTypeRun, run once on first use."""
    global _read_type_run
    if _read_type_run is None:
        require_database()
        try:
            _read_type_run = ReadTypeRun()
        except Exception as e:
            _read_type_run = e
    if isinstance(_read_type_run, Exception):
        raise RuntimeError(f"no read types' run: {_read_type_run}")
    return _read_type_run


class ReadTypeTest(DbTest):
    """A test class on the read types' run (read_type_run())."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.runs = read_type_run()

    def species_called(self, sample):
        return profile_species(self.runs.text(f"{sample}.profile"))

    def species_on_reads(self, sample):
        """The species of the genes simulate_long_reads placed on a long-read sample's reads."""
        names = {taxid: name for name, taxid in species().items()}
        return sorted({names[gene.split("_")[0]] for placed in self.runs.truth[sample] for gene, _, _ in placed})


class SingleEndTest(ReadTypeTest):
    """Single-end reads: the first mates of the simulated samples sa and sb."""

    def reads_by_name(self, prefix):
        return {name.rsplit("/", 1)[0]: seq for name, seq in fastq_records(os.path.join(READS, f"{prefix}_R1.fq"))}

    def test_exit_code_and_model(self):
        self.assertIn("Align the single-end reads of sample sa", self.runs.log)
        self.assertIn("Align the single-end reads of sample sb", self.runs.log)
        self.assertIn("Model of single-end reads: " + os.path.join(self.runs.db, "model_se.xml"), self.runs.log)

    def test_records_are_unpaired_reads(self):
        all_records = sam_records(self.runs.sam("sa"))
        # Unmapped records (flag 4, no gene, no sequence) stand for reads that seeded on taxa but aligned nowhere;
        # they carry the taxa as a ZF tag and nothing else.
        unmapped = [r for r in all_records if int(r[1]) & 0x4]
        for r in unmapped:
            self.assertEqual((r[2], r[3], r[5], r[9], r[10]), ("*", "0", "*", "*", "*"), r)
            self.assertTrue(any(f.startswith("ZF:Z:") and f[5:] for f in r[11:]), f"an unmapped record names its failed candidates: {r}")
        records = [r for r in all_records if not int(r[1]) & 0x4]
        self.assertTrue(records)
        self.assertEqual([r[1] for r in records if int(r[1]) & 0xCD], [], "no pair, mate or unmapped flags")
        self.assertEqual(len({r[0] for r in records}), len(records), "one record per read (-m 1)")
        reads_of = self.reads_by_name("sa")
        for r in records:
            read = reads_of[r[0]]  # QNAME is the read id without its /1
            self.assertEqual(r[9], revcomp(read) if int(r[1]) & 0x10 else read, "SEQ in reference orientation")
            self.assertEqual((r[6], r[7], r[8]), ("*", "0", "0"))
        self.assertTrue(any(int(r[1]) & 0x10 for r in records) and any(not int(r[1]) & 0x10 for r in records))
        self.assertTrue(any(int(r[4]) >= 4 for r in records), "reads with a MAPQ the profiler takes")

    def test_profiles_and_strains(self):
        for sample in ("sa", "sb"):
            self.assertEqual(self.species_called(sample), sorted(species()), f"{sample}: its three species")
            self.assertTrue(os.path.exists(find_one(self.runs.out, f"{sample}_runtime.tsv")))
        for name in species():  # species in both samples get MSAs
            rows = self.runs.msa_rows(name)
            self.assertIn("sa", rows, name)
            self.assertIn("sb", rows, name)


class PacBioTest(ReadTypeTest):
    """PacBio-like long reads of several genes each (la, lb)."""

    def read_seqs(self, prefix):
        records = fastq_records(self.runs.long[prefix])
        return [seq for _, seq in records], [name for name, _ in records]

    def test_exit_code_and_model(self):
        log = self.runs.log
        self.assertIn("Align the PacBio reads of sample la", log)
        self.assertIn("Model of PacBio reads: " + os.path.join(self.runs.db, "model_PB.xml"), log)
        self.assertIn("1 read(s) longer than 65000 bp were seeded in chunks", log)
        self.assertRegex(log, r"\d+ gene hits fit several taxa \(MAPQ < 4\); \d+ gene hits took their read's consensus taxon "
                              r"\(another best hit or a higher MAPQ\), and \d+ had no hit of it or a clearly better one of "
                              r"another taxon \(written with MAPQ 0\)")
        self.assertIn("@CO\tprotal read type: pb\n", sam_text(self.runs.sam("la")))

    def test_records_hold_their_aligned_bases(self):
        seqs, names = self.read_seqs("la")
        read_of = dict(zip(names, seqs))
        records = [r for r in sam_records(self.runs.sam("la")) if not int(r[1]) & 0x4]  # not the unmapped ones (ZF)
        self.assertTrue(records)
        for r in records:
            flag, ops, read = int(r[1]), cigar_ops(r[5]), read_of[r[0]]
            clip_left = ops[0][0] if ops[0][1] == "H" else 0
            clip_right = ops[-1][0] if ops[-1][1] == "H" else 0
            self.assertEqual(clip_left + len(r[9]) + clip_right, len(read), r[:6])
            oriented = revcomp(read) if flag & 0x10 else read
            self.assertEqual(r[9], oriented[clip_left:clip_left + len(r[9])], r[:6])
            self.assertEqual((r[6], r[7], r[8]), ("*", "0", "0"))
            self.assertFalse(flag & 0xCD, "no pair, mate or unmapped flags")
        for name, reps in representative_records(records).items():
            primaries = [r for r in reps if not int(r[1]) & 0x800]
            self.assertEqual(len(primaries), 1, f"one primary record per read: {name}")

    def test_every_gene_is_found_once(self):
        for prefix in ("la", "lb"):
            _, names = self.read_seqs(prefix)
            by_read = representative_records(sam_records(self.runs.sam(prefix)))
            found = missed = twice = 0
            for name, placed in zip(names, self.runs.truth[prefix]):
                reps = by_read.get(name, [])
                for gene, start, end in placed:
                    hits = [r for r in reps if r[2] == gene and sum(n for n, op in cigar_ops(r[5]) if op in "MX=D") >= 0.95 * (end - start)]
                    found += len(hits) >= 1
                    missed += not hits
                    twice += len(hits) > sum(1 for g, _, _ in placed if g == gene)
            self.assertEqual(twice, 0, f"{prefix}: no gene counted twice")
            self.assertGreaterEqual(found / (found + missed), 0.95, f"{prefix}: {found} genes found, {missed} missed")

    def test_genes_across_chunk_boundaries(self):
        _, names = self.read_seqs("la")
        long_name, placed = names[-1], self.runs.truth["la"][-1]
        reps = representative_records(sam_records(self.runs.sam("la"))).get(long_name, [])
        found = [g for g, _, _ in placed if any(r[2] == g for r in reps)]
        self.assertGreaterEqual(len(found), 0.95 * len(placed))
        self.assertLessEqual(len(reps), len(placed), "each gene of the 150 kb read at most once")

    def test_profiles(self):
        for sample in ("la", "lb"):
            self.assertEqual(self.species_called(sample), self.species_on_reads(sample), sample)
        self.assertTrue(os.path.exists(find_one(self.runs.out, "la_runtime.tsv")))

    def test_a_rerun_takes_the_read_type_of_the_sam_it_reuses(self):
        os.makedirs(self.path("rerun"))
        shutil.copy(self.runs.sam("la"), self.path("rerun", "la.sam.zst"))
        rc, log = run(self.work, "--db", self.runs.db, "-1", self.runs.long["la"], "--prefix", "la", "-o", "rerun", "-t", "2",
                      "--no_qcmsa")
        self.assertEqual(rc, 0, log[-3000:])
        self.assertIn("Note: Sample la: PacBio reads", log)  # their names say so, as the SAM's header
        self.assertIn("Model of PacBio reads: " + os.path.join(self.runs.db, "model_PB.xml"), log)
        self.assertEqual(read_text(self.path("rerun", "la.profile")), self.runs.text("la.profile"))
        # A read type given wins, as with --profile_only.
        model = db_file("model_pe.xml")
        rc, log = run(self.work, "--db", self.runs.db, "-1", self.runs.long["la"], "--prefix", "la", "-o", "rerun", "-t", "2",
                      "--no_qcmsa", "--read_type", "ont", "--model_ont", model)
        self.assertEqual(rc, 0, log[-3000:])
        self.assertIn("holds PacBio reads; profiled as ONT reads (ont, --read_type or READ_TYPE; --force aligns them again)", log)
        self.assertIn("Model of ONT reads: " + model, log)

    def test_long_reads_given_as_short_ones_stop(self):
        rc, log = run(self.work, "--db", self.runs.db, "-1", self.runs.long["lb"], "--read_type", "se", "-o", "out_short",
                      "-t", "1", "--no_qcmsa")
        self.assertEqual(rc, 30, log[-3000:])
        self.assertIn("too long for short reads: give --read_type pb", log)


class OntTest(ReadTypeTest):
    """ONT-like long reads (oa, ob)."""

    def test_exit_code_and_model(self):
        log = self.runs.log
        self.assertIn("Align the ONT reads of sample oa (-a 0.85)", log)
        self.assertIn("Model of ONT reads: " + os.path.join(self.runs.db, "model_ONT.xml"), log)
        self.assertIn("@CO\tprotal read type: ont\n", sam_text(self.runs.sam("oa")))
        # The options summary shows the values the ONT reads get.
        self.assertIn("max score ani:       0.900000 (ONT reads: 0.850000)\n", log)
        self.assertIn("snp min af:          0.150000 (ONT reads: 0.200000)\n", log)
        self.assertIn("x-drop:              1000 (short reads; long reads: none)\n", log)

    def test_long_reads_are_aligned_without_x_drop(self):
        # --x_drop is for short reads: over the gene-long windows of long reads it lost alignments.
        rc, log = run(self.work, "--db", self.runs.db, "-1", self.runs.long["oa"], "--prefix", "oa", "--read_type", "ont",
                      "-o", "out_xdrop", "-t", "4", "--no_qcmsa", "--no_profile", "--x_drop", "50")
        self.assertEqual(rc, 0, log[-3000:])
        self.assertEqual(sorted(sam_records(self.path("out_xdrop", "oa.sam"))), sorted(sam_records(self.runs.sam("oa"))))

    def test_every_gene_is_found_once(self):
        for prefix in ("oa", "ob"):
            by_read = representative_records(sam_records(self.runs.sam(prefix)))
            found = missed = twice = 0
            for (name, _), placed in zip(fastq_records(self.runs.long[prefix]), self.runs.truth[prefix]):
                reps = by_read.get(name, [])
                for gene, start, end in placed:
                    hits = [r for r in reps if r[2] == gene and sum(n for n, op in cigar_ops(r[5]) if op in "MX=D") >= 0.9 * (end - start)]
                    found += len(hits) >= 1
                    missed += not hits
                    twice += len(hits) > sum(1 for g, _, _ in placed if g == gene)
            self.assertEqual(twice, 0, f"{prefix}: no gene counted twice")
            self.assertGreaterEqual(found / (found + missed), 0.9, f"{prefix}: {found} genes found, {missed} missed")

    def test_profiles(self):
        # The stand-in model (the paired-end one) is not fit for 2% errors: which species it reports is not the point
        # here, but that every species on the reads reaches it and none other is reported.
        for sample in ("oa", "ob"):
            on_reads = self.species_on_reads(sample)
            scored = sorted(row["Name"] for row in read_dicts(find_one(self.runs.out, f"{sample}.profile.log")))
            self.assertEqual(scored, on_reads, sample)
            called = self.species_called(sample)
            self.assertTrue(called, sample)
            self.assertLessEqual(set(called), set(on_reads), sample)


class ReadTypesTest(ReadTypeTest):
    """What every read type does alike: one map of all of them, --profile_only taking the model of the SAM's reads, FASTA
    reads, and the read type found from the reads."""

    TYPES = {"sa": ("se", "single-end reads"), "la": ("pb", "PacBio reads"), "oa": ("ont", "ONT reads")}

    def test_a_map_mixes_every_read_type(self):
        log = self.runs.log
        self.assertIn("1 paired-end, 2 single-end, 2 PacBio, 2 ONT sample(s)", log)
        self.assertIn("Model of paired-end reads: ", log)
        self.assertIn("Align the paired-end reads of sample pe (-a 0.9)", log)
        self.assertIn("Align the ONT reads of sample ob (-a 0.85)", log)
        self.assertNotIn("minimum allele frequencies for", log)
        paired = sam_records(self.runs.sam("pe"))
        self.assertTrue(paired and all(int(r[1]) & 0x1 for r in paired))
        for sample, (read_type, _) in self.TYPES.items():
            records = sam_records(self.runs.sam(sample))
            self.assertTrue(records and not any(int(r[1]) & 0x1 for r in records), sample)
            if read_type != "se":
                self.assertIn(f"@CO\tprotal read type: {read_type}\n", sam_text(self.runs.sam(sample)))
        self.assertEqual(self.species_called("pe"), sorted(species()))

    def test_profile_only_takes_the_model_of_the_sams_reads(self):
        # The test database has no model for these read types: the SAM's records (unpaired) and its header (@CO) ask
        # for one.
        sams = {sample: self.runs.sam(sample) for sample in self.TYPES}
        for sample, (read_type, name) in self.TYPES.items():
            with self.subTest(read_type):
                rc, log = profile_only(self.work, self.path(f"missing_{read_type}"), [sams[sample]], prefixes=[sample])
                self.assertEqual(rc, 30, log[-3000:])
                self.assertIn(f"no model for --read_type {read_type} ({name})", log)
                self.assertIn(f"--model_{read_type}", log)
        model = db_file("model_pe.xml")
        rc, log = profile_only(self.work, self.path("po"), list(sams.values()), "--model_se", model, "--model_pb", model,
                               "--model_ont", model, prefixes=list(sams))
        self.assertEqual(rc, 0, log[-3000:])
        for sample, (_, name) in self.TYPES.items():
            self.assertIn(f"Model of {name}: {model}", log)
            self.assertEqual(read_text(self.path("po", f"{sample}.profile")), self.runs.text(f"{sample}.profile"), sample)

    def test_fasta_reads(self):
        # FASTA reads get a quality: Q30 (?) short and PacBio reads, Q18 (3) ONT reads; -a given replaces the read
        # types' identities.
        sample_map = self.path("fasta.map")
        expected = {"safa": ("se", os.path.join(READS, "sa_R1.fq"), "?"), "lbfa": ("pb", self.runs.long["lb"], "?"),
                    "obfa": ("ont", self.runs.long["ob"], "3")}
        with open(sample_map, "w") as fh:
            fh.write(f"#OUTPUT_DIR\t{self.path('out_fasta')}\n#SAM_OUTPUT_DIR\t.\n#SAMPLEID\tPREFIX\tFIRST\tSECOND\tREAD_TYPE\n")
            for sample, (read_type, fastq, _) in expected.items():
                with open(self.path(f"{sample}.fa"), "w") as fa:
                    fa.writelines(f">{name.rsplit('/', 1)[0] if read_type == 'se' else name}\n{seq}\n"
                                  for name, seq in fastq_records(fastq))
                fh.write(f"{sample}\t{sample}\t{self.path(sample + '.fa')}\t-\t{read_type}\n")
        rc, log = run(self.work, "--db", self.runs.db, "--map", sample_map, "-a", "0.8", "-t", "4", "--no_qcmsa", "--no_strains")
        self.assertEqual(rc, 0, log[-3000:])
        self.assertIn("Align the ONT reads of sample obfa (-a 0.8)", log)
        for sample, (read_type, _, quality) in expected.items():
            with self.subTest(read_type):
                records = sam_records(find_one(self.path("out_fasta"), f"{sample}.sam.zst"))
                self.assertTrue(records)
                self.assertEqual({q for r in records for q in r[10]}, {quality})
                _, rows = read_table(find_one(self.path("out_fasta"), f"{sample}.profile.log"))
                self.assertTrue(rows, "the profiler takes the records")

    def test_read_types_by_the_reads_and_prefixes_by_the_files(self):
        # Without --read_type: PacBio names (movie/ZMW) make PacBio reads; other names by their quality, Q17 here (ONT
        # reads); short ones single-end reads. Without --prefix the read files name the samples.
        ont = self.path("ont_like.fq")
        simulate_long_reads(ont, 10, seed=31, error=0.02, indels=0.7, read_name="read_{}", quality="2")
        short = head_reads(os.path.join(READS, "sa_R1.fq"), self.path("sa_R1.fq"), 200)
        rc, log = run(self.work, "--db", self.runs.db, "-1", ",".join([self.runs.long["lb"], ont, short]), "-o", "out_auto",
                      "-t", "2", "--no_qcmsa", "--no_profile")
        self.assertEqual(rc, 0, log[-3000:])
        self.assertIn("Note: Sample lb: PacBio reads (the first ", log)
        self.assertIn("named as PacBio names reads", log)
        self.assertIn("@CO\tprotal read type: pb\n", sam_text(self.path("out_auto", "lb.sam")))
        self.assertIn("Note: Sample ont_like: ONT reads (the first 10 reads", log)
        self.assertIn("median read quality Q17.0 (PacBio from Q25, ONT below)", log)
        self.assertIn("@CO\tprotal read type: ont\n", sam_text(self.path("out_auto", "ont_like.sam")))
        self.assertIn("Align the single-end reads of sample sa_R1", log)
        self.assertTrue(sam_records(self.path("out_auto", "sa_R1.sam")))

    def test_a_long_read_after_the_first_100_fails_its_sample(self):
        # The check before aligning sees the first 100 reads; the reader stops at any later long read.
        head = head_reads(os.path.join(READS, "sa_R1.fq"), self.path("head.fq"), 150)
        late = self.path("late_long.fq")
        with open(late, "w") as fh:
            fh.write(read_text(head) + "@long\n" + "ACGT" * 500 + "\n+\n" + "I" * 2000 + "\n" + read_text(head))
        rc, log = run(self.work, "--db", self.runs.db, "-1", late, "--prefix", "late", "-o", "out_late", "-t", "2", "--no_qcmsa")
        self.assertEqual(rc, 1, log[-3000:])
        self.assertIn("read long has 2000 bp, too long for short reads; give --read_type pb or ont", log)
        self.assertFalse(glob.glob(self.path("out_late", "*.sam*")))


class PhasingTest(DbTest):
    """A long-read sample of two strains of a species gets a strain MSA row per strain (Haplotypes.h). The species'
    genes laid out as a genome (300-800 bp apart, in the database's order); strain 2 has 1.5% of the bases of every
    gene substituted. Sample pa holds strain 1, pm strains 1 and 2 at 70:30: PacBio-like reads of 6-10 kb at 0.1%
    errors, about 25x. The MSA has pa, pm_hap1 (strain 1) and pm_hap2 (strain 2); with --no_phasing (the same SAMs) pa
    and pm."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.db = read_type_db()
        rng = random.Random(23)
        genes = [(name, seq) for name, seq in reference_genes() if name.startswith("1_")]
        cls.strains = [dict(genes), {}]
        for name, seq in genes:
            cls.strains[1][name] = "".join(rng.choice([c for c in "ACGT" if c != b]) if rng.random() < 0.015 else b
                                           for b in seq)
        spacers = [("".join(rng.choice("ACGT") for _ in range(rng.randint(300, 800)))) for _ in genes]
        genomes = []
        for strain in cls.strains:
            genome = ""
            for (name, _), spacer in zip(genes, spacers):
                genome += spacer + strain[name]
            genomes.append(genome)

        def sample(path, shares, seed):
            r = random.Random(seed)
            length = len(genomes[0])
            reads = int(25 * length / 8000)
            with open(path, "w") as fh:
                for i in range(reads):
                    genome = genomes[0] if r.random() < shares[0] else genomes[1]
                    size = r.randint(6000, 10000)
                    start = r.randint(0, max(0, length - size))
                    seq = genome[start:start + size]
                    seq = "".join(r.choice([c for c in "ACGT" if c != b]) if r.random() < 0.001 else b for b in seq)
                    if r.random() < 0.5:
                        seq = revcomp(seq)
                    fh.write(f"@m64001_000000/{i + 1}/ccs\n{seq}\n+\n{'I' * len(seq)}\n")

        cls.reads = {"pa": os.path.join(cls.work, "pa.fq"), "pm": os.path.join(cls.work, "pm.fq")}
        sample(cls.reads["pa"], (1.0, 0.0), 1)
        sample(cls.reads["pm"], (0.7, 0.3), 2)
        common = ["--db", cls.db, "-1", ",".join(cls.reads.values()), "--prefix", "pa,pm", "--read_type", "pb", "-t", "4",
                  "--no_qcmsa", "--msa_knob", "0"]
        cls.rc, cls.log = run(cls.work, *common, "-o", "out")
        os.makedirs(os.path.join(cls.work, "off"))
        for prefix in ("pa", "pm"):
            if os.path.exists(os.path.join(cls.work, "out", f"{prefix}.sam.zst")):
                shutil.copy(os.path.join(cls.work, "out", f"{prefix}.sam.zst"), os.path.join(cls.work, "off"))
        cls.rc_off, cls.log_off = run(cls.work, *common, "-o", "off", "--no_phasing")
        # The same SAMs once more, each sample's strain evidence (its long reads' phase records too) spilled to a file
        # and read back per species (--strain_spill).
        os.makedirs(os.path.join(cls.work, "spilled"))
        for prefix in ("pa", "pm"):
            if os.path.exists(os.path.join(cls.work, "out", f"{prefix}.sam.zst")):
                shutil.copy(os.path.join(cls.work, "out", f"{prefix}.sam.zst"), os.path.join(cls.work, "spilled"))
        cls.rc_spill, cls.log_spill = run(cls.work, *common, "-o", "spilled", "--strain_spill", os.path.join(cls.work, "spill"))

    def test_spilled_evidence_gives_the_same_rows(self):
        self.assertEqual(self.rc_spill, 0, self.log_spill[-3000:])
        names = sorted(os.path.basename(p) for pattern in ("*.raw.msa.fna", "*.raw.partition.txt", "*.meta.tsv", "*.haplotypes.tsv")
                       for p in glob.glob(self.path("out", "strains", pattern)))
        self.assertTrue(any(name.endswith(".haplotypes.tsv") for name in names), names)
        for name in names:
            self.assertEqual(read_text(self.path("spilled", "strains", name)), read_text(self.path("out", "strains", name)), name)
        self.assertEqual(os.listdir(self.path("spill")), [], "the spill files are removed")

    def msa(self, out):
        paths = []
        for path in glob.glob(self.path(out, "strains", "*.raw.msa.fna")):
            if ">pa\n" in read_text(path):
                paths.append(path)
        self.assertEqual(len(paths), 1, paths)
        rows, name = {}, None
        for line in read_text(paths[0]).splitlines():
            if line.startswith(">"):
                name = line[1:]
                rows[name] = ""
            else:
                rows[name] += line
        return paths[0], rows

    def test_each_strain_gets_its_row(self):
        self.assertEqual(self.rc, 0, self.log[-3000:])
        self.assertRegex(self.log, r"pm: \d+ long reads at \d+ blocks of multi-allelic sites .*; 2 strain rows \(shares 0\.[67]")
        path, rows = self.msa("out")
        self.assertEqual(sorted(n for n in rows if not n.endswith("_reference")), ["pa", "pm_hap1", "pm_hap2"])
        reference = next(seq for n, seq in rows.items() if n.endswith("_reference"))
        # The MSA's columns of reference bases in gene order (the partitions), against the strains' genes.
        partitions = [line for line in read_text(path.replace(".raw.msa.fna", ".raw.partition.txt")).splitlines() if line]
        strain_of = {"pm_hap1": 0, "pm_hap2": 1}
        matched = {r: [0, 0] for r in strain_of}  # bases of the row's own strain, of the other one, at differing sites
        for part in partitions:
            gene, span = re.match(r"DNA, gene(\d+) = (\d+)-(\d+)", part).group(1), re.match(r".* = (\d+)-(\d+)", part).groups()
            columns = [i for i in range(int(span[0]) - 1, int(span[1])) if reference[i] != "-"]
            name = f"1_{gene}"
            if name not in self.strains[0]:
                continue
            for k, i in enumerate(columns):
                a, b = self.strains[0][name][k], self.strains[1][name][k]
                if a == b:
                    continue
                for row, s in strain_of.items():
                    c = rows[row][i]
                    own, other = (a, b) if s == 0 else (b, a)
                    matched[row][0] += c == own
                    matched[row][1] += c == other
        for row, (own, other) in matched.items():
            self.assertGreater(own, 0, row)
            self.assertLessEqual(other, 0.02 * own, (row, own, other))
        lines = read_dicts(path.replace(".raw.msa.fna", ".haplotypes.tsv"))
        self.assertTrue(lines)
        self.assertTrue(all(line["sample"] == "pm" for line in lines))
        self.assertTrue(any(line["phased"] == "yes" for line in lines))
        meta = read_table(path.replace(".raw.msa.fna", ".meta.tsv"))[1]
        self.assertEqual({row[0] for row in meta}, {"pa", "pm_hap1", "pm_hap2"})

    def test_one_row_per_sample_without_phasing(self):
        self.assertEqual(self.rc_off, 0, self.log_off[-3000:])
        self.assertIn("All alignments are present", self.log_off)
        path, rows = self.msa("off")
        self.assertEqual(sorted(n for n in rows if not n.endswith("_reference")), ["pa", "pm"])
        self.assertGreater(sum(c in "RYSWKMBDHV" for c in rows["pm"]), 0)  # the mixture's IUPAC codes
        self.assertFalse(os.path.exists(path.replace(".raw.msa.fna", ".haplotypes.tsv")))


class VersionTest(unittest.TestCase):
    """--version of protal and the simulator: the version of CMakeLists.txt, and the commit they were built from when
    built in a git checkout (build_gtdb_database.py checks the binaries it runs by it). No database needed."""

    def test_both_say_the_version_and_commit(self):
        version = re.search(r"project\(protal VERSION ([0-9.]+)\)", read_text(os.path.join(ROOT, "CMakeLists.txt"))).group(1)
        for binary, name, variable in ((PROTAL, "protal", "PROTAL"), (SIMULATE, "simulate_metagenomes", "SIMULATE")):
            with self.subTest(name):
                require_binary(binary, name, variable)
                said = subprocess.run([binary, "--version"], capture_output=True, text=True, check=True).stdout
                self.assertRegex(said, rf"^{name} v{re.escape(version)}"
                                       r"( \(commit [0-9a-f]{40}(, with uncommitted changes)?\))?\n$")


class SimulatorTest(WorkDir):
    """simulate_metagenomes on genomes of its own (no database needed)."""

    @classmethod
    def setUpClass(cls):
        require_binary(SIMULATE, "simulate_metagenomes", "SIMULATE")
        super().setUpClass()

    def test_too_few_read_pairs_fails_fast(self):
        with open(self.path("genomes.tsv"), "w") as table:
            for sp in range(1, 4):
                for strain in "ab":
                    fasta = self.path(f"g{sp}{strain}.fa")
                    with open(fasta, "w") as fh:
                        fh.write(">c1\n" + "".join("ACGT"[(i * 7 + sp) % 4] for i in range(3000)) + "\n")
                    table.write(f"g{sp}{strain}\td__B;p__P;c__C;o__O;f__F;g__G;s__G sp{sp}\t{fasta}\n")
        # Used to spin forever when there are fewer read pairs than species.
        rc, log = run(self.work, "--genome_table", "genomes.tsv", "--test", "--seed", "1",
                      "--total_read_pairs", "2", "--species_per_sample", "3", "--output_dir", "sim",
                      binary=SIMULATE, timeout=60)
        self.assertEqual(rc, 1, log)
        self.assertIn("increase --total_read_pairs", log)

    def test_taxon_quota_counts(self):
        # --taxon d__A:2 asks for two species of d__A per sample, the rest drawn from all species. It used
        # to fill every sample with d__A species (the quota's copy was never counted down).
        with open(self.path("quota_genomes.tsv"), "w") as table:
            for domain, count in (("A", 6), ("B", 30)):
                for sp in range(count):
                    fasta = self.path(f"{domain}{sp}.fa")
                    with open(fasta, "w") as fh:
                        fh.write(">c1\n" + "".join("ACGT"[(i * 7 + sp) % 4] for i in range(3000)) + "\n")
                    table.write(f"{domain}{sp}\td__{domain};p__P{domain};c__C{domain};o__O{domain};f__F{domain};"
                                f"g__G{domain};s__G{domain} sp{sp}\t{fasta}\n")
        rc, log = run(self.work, "--genome_table", "quota_genomes.tsv", "--test", "--seed", "1", "--samples", "8",
                      "--total_read_pairs", "1000", "--species_per_sample", "6", "--taxon", "d__A:2",
                      "--output_dir", "sim_quota", binary=SIMULATE, timeout=60)
        self.assertEqual(rc, 0, log)
        domains = {}
        for row in read_dicts(self.path("sim_quota", "manifest.tsv")):
            domains.setdefault(row["sample"], []).append(row["taxonomy"].split(";")[0])
        self.assertEqual(len(domains), 8)
        for sample, found in domains.items():
            self.assertEqual(len(found), 6, sample)
            self.assertGreaterEqual(found.count("d__A"), 2, sample)
        # The four species drawn at random are d__A 4 times in 34, so about 2.5 d__A per sample, not 6.
        self.assertLess(sum(f.count("d__A") for f in domains.values()) / len(domains), 3.5)

    def test_species_counts_in_turn(self):
        # --species_per_sample 2,5,3: the samples' species counts in turn, cyclically (collect_training_data.py spreads
        # a scenario's richness over its samples itself); a count of 0 is refused.
        with open(self.path("count_genomes.tsv"), "w") as table:
            for sp in range(8):
                fasta = self.path(f"n{sp}.fa")
                with open(fasta, "w") as fh:
                    fh.write(">c1\n" + "".join("ACGT"[(i * 7 + sp) % 4] for i in range(3000)) + "\n")
                table.write(f"n{sp}\td__B;p__P;c__C;o__O;f__F;g__G;s__G n{sp}\t{fasta}\n")
        rc, log = run(self.work, "--genome_table", "count_genomes.tsv", "--test", "--seed", "1", "--samples", "4",
                      "--total_read_pairs", "1000", "--species_per_sample", "2,5,3", "--output_dir", "sim_counts",
                      binary=SIMULATE, timeout=60)
        self.assertEqual(rc, 0, log)
        species = {}
        for row in read_dicts(self.path("sim_counts", "manifest.tsv")):
            species.setdefault(row["sample"], set()).add(row["taxonomy"])
        self.assertEqual([len(species[s]) for s in sorted(species, key=lambda s: int(s.rsplit("_", 1)[1]))], [2, 5, 3, 2])
        rc, log = run(self.work, "--genome_table", "count_genomes.tsv", "--test", "--seed", "1", "--samples", "2",
                      "--total_read_pairs", "1000", "--species_per_sample", "2,0", "--output_dir", "sim_zero",
                      binary=SIMULATE, timeout=60)
        self.assertNotEqual(rc, 0, log)
        self.assertIn("--species_per_sample takes", log)

    def test_samples_on_threads_are_the_same(self):
        # -t: samples written side by side (their designs and the reads' seeds drawn first, in order), each sample's
        # reads made in process (IlluminaSimulator) and BGZF-compressed in the genomes' order: the same files byte for
        # byte on 1 and 3 threads, and no temporary files left.
        with open(self.path("art_genomes.tsv"), "w") as table:
            for sp in range(6):
                fasta = self.path(f"t{sp}.fa.gz")
                rng = random.Random(sp)
                with gzip.open(fasta, "wt") as fh:
                    fh.write(">c1\n" + "".join(rng.choice("ACGT") for _ in range(6000)) + "\n")
                table.write(f"t{sp}\td__B;p__P;c__C;o__O;f__F;g__G;s__G sp{sp}\t{fasta}\n")
        common = ["--genome_table", "art_genomes.tsv", "--seed", "4", "--samples", "4", "--total_read_pairs", "600",
                  "--species_per_sample", "3", "--read_length", "100", "--sequencer", "HS20", "--fragment_mean", "300",
                  "--fragment_stdev", "30"]
        outputs = {}
        for threads in ("1", "3"):
            out = f"sim_t{threads}"
            rc, log = run(self.work, *common, "-t", threads, "--output_dir", out, binary=SIMULATE, timeout=300)
            self.assertEqual(rc, 0, log[-3000:])
            files = sorted(os.listdir(self.path(out, "reads")))
            self.assertEqual(files, [f"sample_{i}_R{r}.fq.gz" for i in range(1, 5) for r in (1, 2)])
            self.assertEqual([f for f in os.listdir(self.path(out)) if f.endswith("_tmp")], [], "temporary files left")
            outputs[threads] = {}
            for f in files:
                with open(self.path(out, "reads", f), "rb") as fh:
                    outputs[threads][f] = fh.read()
            outputs[threads]["manifest"] = re.sub(re.escape(out), "OUT", read_text(self.path(out, "manifest.tsv"))).encode()
        self.assertEqual(outputs["1"], outputs["3"])
        with gzip.open(self.path("sim_t1", "reads", "sample_1_R1.fq.gz"), "rt") as fh:
            self.assertEqual(sum(1 for _ in fh) // 4, 600)
        # --reads_compression zstd: .fq.zst files of the same reads (the training data collector's default), named so in
        # the protal map; the same bytes on any number of threads.
        require_zstd("to read zstd-compressed reads")
        zstd = {}
        for threads in ("1", "3"):
            out = f"sim_zstd_t{threads}"
            rc, log = run(self.work, *common, "-t", threads, "--output_dir", out, "--reads_compression", "zstd",
                          "--protal_metafile", self.path(out, "p"), binary=SIMULATE, timeout=300)
            self.assertEqual(rc, 0, log[-3000:])
            files = sorted(os.listdir(self.path(out, "reads")))
            self.assertEqual(files, [f"sample_{i}_R{r}.fq.zst" for i in range(1, 5) for r in (1, 2)])
            zstd[threads] = {}
            for f in files:
                with open(self.path(out, "reads", f), "rb") as fh:
                    zstd[threads][f] = fh.read()
            self.assertIn("sample_1\tsample_1_R1.fq.zst\tsample_1_R2.fq.zst\tsample_1.sam.zst\tsample_1\t",
                          read_text(self.path(out, "protal.meta")))
        self.assertEqual(zstd["1"], zstd["3"])
        for name, data in zstd["1"].items():
            self.assertEqual(data[:4], b"\x28\xb5\x2f\xfd", name)
            plain = subprocess.run(["zstd", "-dcq"], input=data, check=True, stdout=subprocess.PIPE).stdout
            self.assertEqual(plain, gzip.decompress(outputs["1"][name.replace(".zst", ".gz")]), name)


if __name__ == "__main__":
    sys.exit(unittest.main())
