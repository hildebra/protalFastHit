#!/usr/bin/env python3
"""test_protal_e2e.py - end-to-end checks of protal and simulate_metagenomes on a small database.

Simulates paired reads from the database's reference genes (sequencing errors only, fixed seed)
and runs the real binaries, checking what the unit tests cannot: exit codes, output files, SAM
records, strain MSAs, reruns, and that failures are reported.

  PROTAL_TEST_DB=data/mini_db/protal_db python3 -m unittest -v tests/e2e/test_protal_e2e.py
  just e2e                                  # builds the mini DB first

PROTAL_TEST_DB  protal database, e.g. from scripts/mini_db/build_mini_db.sh (required; only read):
                the single file database.protal (or its folder), or separate raw or zstd-compressed
                files (index.prx.zst, reference.fna.zst). Tests that read the database's files
                get them unpacked (protal --unpack_db) into a temporary folder. The zstd CLI is
                needed for a compressed database.
PROTAL          protal binary (default: build/protal)
SIMULATE        simulate_metagenomes binary (default: build/simulate_metagenomes; optional)
"""

import csv
import filecmp
import glob
import gzip
import os
import random
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
import unittest

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
DB = os.environ.get("PROTAL_TEST_DB", "")
DB = os.path.abspath(DB) if DB else ""  # protal runs in temporary folders
PROTAL = os.path.abspath(os.environ.get("PROTAL", os.path.join(ROOT, "build", "protal")))
SIMULATE = os.path.abspath(os.environ.get("SIMULATE", os.path.join(ROOT, "build", "simulate_metagenomes")))
QCMSA = os.path.join(ROOT, "scripts", "qcmsa.py")
READS = None  # directory with the simulated reads, set up once per module
FILES = DB     # the database's separate files: DB, or DB unpacked if it is a single file
UNPACKED = None


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


def setUpModule():
    global READS, FILES, UNPACKED
    if not DB or not os.path.exists(DB):
        raise unittest.SkipTest("set PROTAL_TEST_DB to a protal database (just mini-db builds data/mini_db/protal_db)")
    if not os.access(PROTAL, os.X_OK):
        raise unittest.SkipTest(f"protal binary not found at {PROTAL} (set PROTAL)")
    bundle = single_file(DB)
    if bundle:
        UNPACKED = tempfile.mkdtemp(prefix="protal_e2e_db_")
        rc, log = run(UNPACKED, "--unpack_db", "--db", bundle, "--unpack_dir", UNPACKED, "-t", "4")
        if rc != 0:
            raise RuntimeError("protal --unpack_db failed:\n" + log[-3000:])
        FILES = UNPACKED
    index = db_file("index.prx")
    if not os.path.isfile(index) or os.path.getsize(index) == 0:
        raise unittest.SkipTest("set PROTAL_TEST_DB to a protal database (just mini-db builds data/mini_db/protal_db)")
    if db_file("reference.fna").endswith(".zst") and not shutil.which("zstd"):
        raise unittest.SkipTest("the zstd CLI is needed to read reference.fna.zst")
    READS = tempfile.mkdtemp(prefix="protal_e2e_reads_")
    simulate_reads("sa", pairs_per_gene=12, seed=1)
    simulate_reads("sb", pairs_per_gene=12, seed=2)
    simulate_reads("sr", pairs_per_gene=12, seed=3, random_r2=True)


def tearDownModule():
    for d in (READS, UNPACKED):
        if d:
            shutil.rmtree(d, ignore_errors=True)


def revcomp(seq):
    return seq[::-1].translate(str.maketrans("ACGTN", "TGCAN"))


def reference_genes():
    genes, name = [], None
    path = db_file("reference.fna")
    if path.endswith(".zst"):
        lines = subprocess.run(["zstd", "-dc", path], check=True, stdout=subprocess.PIPE, text=True).stdout.splitlines()
    else:
        with open(path) as fh:
            lines = fh.read().splitlines()
    for line in lines:
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
        elif line:
            genes.append((name, line.upper()))
    return genes


def simulate_reads(prefix, pairs_per_gene, seed, random_r2=False, pairs_per_kb=None):
    """Write <prefix>_R1.fq/_R2.fq in READS: 220-320 bp fragments, 100 bp reads, both orientations,
    0.5% substitutions and a few '#' (Q2) bases. With random_r2 the second mate cannot align.
    With pairs_per_kb, genes get pairs in proportion to their length (even depth), shorter genes
    included, instead of pairs_per_gene each."""
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
        for _, gene in reference_genes():
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
    return n


def run(cwd, *args, binary=None, timeout=900):
    """Run protal (or `binary`) in cwd; returns (exit code, combined output)."""
    proc = subprocess.run(["timeout", str(timeout), binary or PROTAL, *args], cwd=cwd,
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, errors="replace")
    return proc.returncode, proc.stdout


def reads(*prefixes):
    """-1/-2/--prefix arguments for the given simulated samples."""
    return ["-1", ",".join(os.path.join(READS, f"{p}_R1.fq") for p in prefixes),
            "-2", ",".join(os.path.join(READS, f"{p}_R2.fq") for p in prefixes),
            "--prefix", ",".join(prefixes)]


def single_reads(*prefixes):
    """-1/--prefix arguments: the first mates of the given simulated samples, as single-end reads."""
    return ["-1", ",".join(os.path.join(READS, f"{p}_R1.fq") for p in prefixes), "--prefix", ",".join(prefixes)]


def sam_records(path):
    with open(path) as fh:
        return [line.rstrip("\n").split("\t") for line in fh if not line.startswith("@")]


class WorkDir(unittest.TestCase):
    """A test class with its own scratch directory."""

    @classmethod
    def setUpClass(cls):
        cls.work = tempfile.mkdtemp(prefix=f"protal_e2e_{cls.__name__}_")

    @classmethod
    def tearDownClass(cls):
        if not os.environ.get("PROTAL_TEST_KEEP"):
            shutil.rmtree(cls.work, ignore_errors=True)

    def path(self, *parts):
        return os.path.join(self.work, *parts)


class CompleteRunTest(WorkDir):
    """One run over two normal samples and one whose read2 mates cannot align."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.rc, cls.log = run(cls.work, "--db", DB, *reads("sa", "sb", "sr"), "-o", "out", "-t", "4", "--no_qcmsa")

    def test_exit_code(self):
        self.assertEqual(self.rc, 0, self.log[-3000:])

    def test_outputs(self):
        self.assertTrue(glob.glob(self.path("out", "sa*.sam")))
        self.assertTrue(glob.glob(self.path("out", "sa*.profile")))
        self.assertEqual(glob.glob(self.path("out", "**", "*.partial"), recursive=True), [], "no .partial SAM left")
        # strains/ and misc/ were never created in -1/-2/-o mode
        self.assertTrue(os.path.isdir(self.path("out", "strains")))
        self.assertTrue(os.path.isdir(self.path("out", "misc")))
        self.assertTrue(glob.glob(self.path("out", "strains", "*.raw.msa.fna")))

    def test_read_names_and_pairs(self):
        records = sam_records(glob.glob(self.path("out", "sa*.sam"))[0])
        self.assertTrue(records)
        bad = [r[0] for r in records if not (r[0].startswith("sa.") and r[0][3:].isdigit())]
        self.assertEqual(bad[:5], [], "QNAME is the read id without its /1 /2 suffix, nothing more")
        self.assertTrue(any(int(r[1]) & 0x2 for r in records), "proper pairs are flagged")

    def test_read1_only_pairs_are_written(self):
        records = [r for r in sam_records(glob.glob(self.path("out", "sr*.sam"))[0]) if not int(r[1]) & 0x100]
        read1_only = [r for r in records if int(r[1]) & 0x40 and int(r[1]) & 0x8]
        with open(os.path.join(READS, "sr_R1.fq")) as fh:
            total = sum(1 for _ in fh) // 4
        self.assertGreater(len(read1_only), total // 2, "primary read1 records with mate unmapped (0x8)")

    def test_reference_calls_are_not_blanked(self):
        bases = ns = 0
        for msa in glob.glob(self.path("out", "strains", "*.raw.msa.fna")):
            with open(msa) as fh:
                for line in fh:
                    if not line.startswith(">"):
                        seq = line.strip()
                        bases += len(seq)
                        ns += seq.count("N")
        self.assertGreater(bases, 0)
        self.assertLess(ns / bases, 0.005, f"N fraction {100 * ns / bases:.3f}% with error-only reads")

        stats = glob.glob(self.path("out", "strains", "*.snp_stats.tsv"))[0]
        with open(stats) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            column = header.index("refs_retained")
            retained = sum(int(line.split("\t")[column]) for line in fh)
        self.assertGreater(retained, 0, "reference calls retained at variant positions")


    def test_partitions_are_one_based_and_cover_the_msa(self):
        for part in glob.glob(self.path("out", "strains", "*.raw.partition.txt")):
            with open(part) as fh:
                ranges = [tuple(int(x) for x in line.split("=")[1].split("-")) for line in fh if line.strip()]
            with open(part.replace(".raw.partition.txt", ".raw.msa.fna")) as fh:
                length = len([line for line in fh if not line.startswith(">")][0].strip())
            self.assertEqual(ranges[0][0], 1, part)
            self.assertEqual(ranges[-1][1], length, part)
            for (_, end), (start, _) in zip(ranges, ranges[1:]):
                self.assertEqual(start, end + 1, part)


def read_table(path):
    """A tab-separated file with a header line, as (header, rows)."""
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        return header, [line.rstrip("\n").split("\t") for line in fh]


class OutputFilesTest(WorkDir):
    """One sample from a map, with a truth file: what the profile's companion files hold."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        truth = os.path.join(cls.work, "truth.tsv")
        with open(truth, "w") as fh:
            fh.write("1\n2\n3\n14\n")  # 14 is a genus: in the taxonomy, not in the database
        sample_map = os.path.join(cls.work, "samples.map")
        with open(sample_map, "w") as fh:
            fh.write(f"#OUTPUT_DIR\t{os.path.join(cls.work, 'out')}\n#INPUT_DIR\t{READS}\n")
            fh.write("#SAMPLEID\tPREFIX\tFIRST\tSECOND\tPROFILE_TRUTH\n")
            fh.write(f"sample_a\tpa\tsa_R1.fq\tsa_R2.fq\t{truth}\n")
        cls.rc, cls.log = run(cls.work, "--db", DB, "--map", sample_map, "-t", "2", "--no_qcmsa")

    def test_exit_code(self):
        self.assertEqual(self.rc, 0, self.log[-3000:])

    def test_the_sample_id_names_the_sample(self):
        header, rows = read_table(self.path("out", "pa.profile.gene.log"))
        self.assertEqual(header[0], "Sample")
        self.assertTrue(rows)
        self.assertEqual({row[0] for row in rows}, {"sample_a"})
        self.assertGreater(max(int(row[header.index("CoverageSum")]) for row in rows), 0)

    def test_truth_counts_name_the_sample(self):
        self.assertRegex(self.log, r"Sample sample_a: TP \d+, FP \d+, FN \d+ \(and 1 true species not in the database\)")
        self.assertNotIn("Truth: 0", self.log)

    def test_statistics_for_a_single_sample(self):
        stats = glob.glob(self.path("out", "misc", "*.statistics.tsv"))
        self.assertTrue(stats, "written with one sample, too")
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
        called = abundance = 0
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


class MateAssignmentTest(WorkDir):
    """Fragments whose best alignment is mate 2's alone, and fragments over two genes, reach the SAM."""

    def align(self, name, pairs):
        with open(self.path(f"{name}_R1.fq"), "w") as r1, open(self.path(f"{name}_R2.fq"), "w") as r2:
            for i, (s1, s2) in enumerate(pairs, 1):
                r1.write(f"@{name}.{i}/1\n{s1}\n+\n{'I' * len(s1)}\n")
                r2.write(f"@{name}.{i}/2\n{s2}\n+\n{'I' * len(s2)}\n")
        rc, log = run(self.work, "--db", DB, "-1", self.path(f"{name}_R1.fq"), "-2", self.path(f"{name}_R2.fq"),
                      "--prefix", name, "-o", "out", "-t", "2", "--no_qcmsa", "--no_strains")
        self.assertEqual(rc, 0, log[-3000:])
        records = sam_records(glob.glob(self.path("out", f"{name}*.sam"))[0])
        return [r for r in records if not int(r[1]) & 0x100]

    def test_pairs_where_only_mate2_aligns(self):
        rng = random.Random(21)
        pairs = []
        for _, gene in reference_genes():
            if len(gene) >= 320:
                start = rng.randint(0, len(gene) - 100)
                pairs.append(("".join(rng.choice("ACGT") for _ in range(100)), revcomp(gene[start:start + 100])))
        records = self.align("mate2", pairs)
        read2_only = [r for r in records if int(r[1]) & 0x80 and int(r[1]) & 0x8]
        self.assertGreater(len(read2_only), 0.9 * len(pairs), f"{len(read2_only)} of {len(pairs)} written")

    def test_pairs_over_two_genes(self):
        genes = {}
        for name, seq in reference_genes():
            taxid, gene = (int(x) for x in name.split("_")[:2])
            genes[(taxid, gene)] = seq
        pairs, expected = [], []
        for (taxid, gene), seq in sorted(genes.items()):
            nxt = genes.get((taxid, gene + 1))
            if nxt is not None and len(seq) >= 100 and len(nxt) >= 100:
                pairs.append((seq[-100:], revcomp(nxt[:100])))
                expected.append((f"{taxid}_{gene}", f"{taxid}_{gene + 1}"))
        records = self.align("junction", pairs)
        by_read = {}
        for r in records:
            by_read.setdefault(r[0], []).append(r)
        good = 0
        for i, (gene_a, gene_b) in enumerate(expected, 1):
            recs = by_read.get(f"junction.{i}", [])
            mates = {("1" if int(r[1]) & 0x40 else "2"): r for r in recs}
            if (len(recs) == 2 and mates.get("1", [None, None, None])[2] == gene_a and
                    mates.get("2", [None, None, None])[2] == gene_b and all(int(r[4]) >= 4 for r in recs)):
                good += 1
        self.assertGreater(good, 0.9 * len(pairs), f"{good} of {len(pairs)} fragments written with both mates")


class MsaSampleSelectionTest(WorkDir):
    """A species' MSA takes only the samples in which the model accepts the species."""

    def test_rejected_samples_are_left_out(self):
        # A sample with three read pairs of Mockella alpha, too few to call it.
        taxid = None
        with open(db_file("internal_taxonomy.dmp")) as fh:
            for line in fh:
                f = line.split("\t")
                if f[3] == "s__Mockella alpha":
                    taxid = f[0]
        genes = [seq for name, seq in reference_genes() if name.split("_")[0] == taxid and len(seq) >= 300][:3]
        with open(self.path("few_R1.fq"), "w") as r1, open(self.path("few_R2.fq"), "w") as r2:
            for i, gene in enumerate(genes, 1):
                r1.write(f"@few.{i}/1\n{gene[:100]}\n+\n{'I' * 100}\n")
                r2.write(f"@few.{i}/2\n{revcomp(gene[200:300])}\n+\n{'I' * 100}\n")
        first = ",".join([os.path.join(READS, "sa_R1.fq"), os.path.join(READS, "sb_R1.fq"), self.path("few_R1.fq")])
        second = ",".join([os.path.join(READS, "sa_R2.fq"), os.path.join(READS, "sb_R2.fq"), self.path("few_R2.fq")])
        rc, log = run(self.work, "--db", DB, "-1", first, "-2", second, "--prefix", "sa,sb,few", "-o", "out",
                      "-t", "2", "--no_qcmsa", "--msa_min_hcov", "0")
        self.assertEqual(rc, 0, log[-3000:])
        with open(self.path("out", "few.profile")) as fh:
            self.assertNotIn("Mockella alpha", fh.read(), "the three pairs do not call the species")
        with open(self.path("out", "strains", "s__Mockella_alpha.raw.msa.fna")) as fh:
            names = [line[1:].strip() for line in fh if line.startswith(">")]
        self.assertIn("sa", names)
        self.assertIn("sb", names)
        self.assertNotIn("few", names)


class StrainEdgeCaseTest(WorkDir):
    def test_species_without_msa_genes(self):
        # No position reaches --msa_min_depth, so no species has MSA columns (this used to segfault).
        # Two samples: MSAs are built only across samples.
        rc, log = run(self.work, "--db", DB, *reads("sa", "sb"), "-o", "out", "-t", "2", "--no_qcmsa",
                      "--msa_min_depth", "100000", "--msa_min_hcov", "0")
        self.assertEqual(rc, 0, log[-3000:])
        self.assertIn("has enough coverage for an MSA", log)
        self.assertEqual(glob.glob(self.path("out", "strains", "*.raw.msa.fna")), [])


class MSAKnobTest(WorkDir):
    """A species' MSA holds the samples whose profile reports it; --msa_knob sets another threshold."""

    def calls(self, sample):
        """{species: (reported, probability)} of a sample's profile, species spelled as in species.tsv."""
        with open(self.path("out", f"{sample}.profile.log")) as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        return {r["Name"].replace(" ", "_"): (r["Predicted"] == "1", float(r["Probability"])) for r in rows}

    def species_list(self):
        with open(self.path("out", "strains", "species.tsv")) as fh:
            return {r["species"]: int(r["samples"]) for r in csv.DictReader(fh, delimiter="\t")}

    def test_msa_samples_mirror_the_profiles(self):
        rc, log = run(self.work, "--db", DB, *reads("sa", "sb"), "-o", "out", "-t", "2", "--no_qcmsa",
                      "--msa_min_hcov", "0")
        self.assertEqual(rc, 0, log[-3000:])
        calls = {s: self.calls(s) for s in ("sa", "sb")}
        listed = self.species_list()
        for species in set(calls["sa"]) | set(calls["sb"]):
            reported = sum(calls[s].get(species, (False, 0))[0] for s in calls)
            if reported >= 2:
                self.assertEqual(listed.get(species), reported, species)
            else:
                self.assertNotIn(species, listed, "an MSA needs 2 samples that report the species")

    def test_msa_knob_admits_unreported_species(self):
        # No species scores --knob 1, yet --msa_knob 0 builds the MSAs of all species with reads in both
        # samples. At 2.4x on every gene, the species' reads are strong evidence: those the profiles
        # leave out are listed in unreported_species.tsv.
        rc, log = run(self.work, "--db", DB, *reads("sa", "sb"), "-o", "out", "-t", "2", "--no_qcmsa",
                      "--msa_min_hcov", "0", "--knob", "1", "--msa_knob", "0")
        self.assertEqual(rc, 0, log[-3000:])
        calls = {s: self.calls(s) for s in ("sa", "sb")}
        listed = self.species_list()
        for species in set(calls["sa"]) & set(calls["sb"]):
            self.assertEqual(listed.get(species), 2, species)
        with open(self.path("out", "misc", "unreported_species.tsv")) as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        self.assertTrue(rows, log[-3000:])
        for r in rows:
            reported, probability = calls[r["sample"]][r["species"].replace(" ", "_")]
            self.assertFalse(reported, r)
            self.assertLess(probability, 1)
            self.assertEqual(r["passes_msa_knob"], "yes")
        self.assertIn("are not reported, although their own reads are strong evidence", log)


class LowCoverageAbundanceTest(WorkDir):
    """The depth estimate, and so relative abundances, stays proportional at low coverage."""

    @staticmethod
    def taxon_depths(genes_log):
        depths = {}
        with open(genes_log) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            taxid, vcov = header.index("TaxID"), header.index("TaxVCOV")
            for line in fh:
                fields = line.rstrip("\n").split("\t")
                depths[fields[taxid]] = float(fields[vcov])
        return depths

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
        self.assertTrue(full)
        for name, fraction in fractions.items():
            depths = self.taxon_depths(self.path("out", f"{name}.profile.genes.log"))
            for taxid, depth in depths.items():
                ratio = depth / (fraction * full[taxid])
                self.assertLess(abs(ratio - 1), 0.3, f"{name} ({fraction:.3f} of the reads), taxon {taxid}: "
                                                     f"depth {depth:.4f} vs {fraction * full[taxid]:.4f} expected")


class ModelContractTest(WorkDir):
    """The training dump holds the features the model is scored with, and --no_strains changes no profile."""

    def test_truth_annotation_has_the_model_features(self):
        truth = self.path("truth.tsv")
        with open(truth, "w") as fh:
            fh.write("1\n2\n3\n")
        rc, log = run(self.work, "--db", DB, *reads("sa"), "-o", "out", "-t", "2", "--no_strains",
                      "--profile_truth", truth)
        self.assertEqual(rc, 0, log[-3000:])
        with open(self.path("out", "sa.profile.truth_annotated")) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            rows = [dict(zip(header, line.rstrip("\n").split("\t"))) for line in fh]
        with open(db_file("model_pe.xml")) as fh:
            fields = set(re.findall(r'<DataField name="([^"]+)"', fh.read())) - {"truth"}
        self.assertEqual(sorted(fields - set(header)), [], "every model input is in the training dump")
        self.assertTrue(rows)
        for row in rows:
            for prefix, counts in (("RAF", "AF"), ("RA", "A")):
                total = sum(float(row[f"{counts}{i}"]) for i in range(5))
                for i in range(5):
                    expected = float(row[f"{counts}{i}"]) / total if total else 0
                    self.assertAlmostEqual(float(row[f"{prefix}{i}"]), expected, places=9, msg=f"{prefix}{i}")
            self.assertEqual(row["truth"], "1")

    def test_no_strains_changes_no_profile(self):
        outputs = {}
        for name, extra in (("strains", []), ("no_strains", ["--no_strains"])):
            rc, log = run(self.work, "--db", DB, *reads("sa", "sb"), "-o", name, "-t", "2", "--no_qcmsa", *extra)
            self.assertEqual(rc, 0, log[-3000:])
            outputs[name] = {}
            for f in glob.glob(self.path(name, "*.profile*")):
                with open(f) as fh:
                    outputs[name][os.path.basename(f)] = fh.read()
        self.assertTrue(outputs["strains"])
        self.assertEqual(outputs["strains"], outputs["no_strains"])


class BuildUniquenessTest(WorkDir):
    """--build checks every k-mer against the full reference, also those whose core occurs once in the index."""

    def test_a_gene_another_taxon_carries_is_not_unique(self):
        rng = random.Random(8)
        genes = {(1, 1): None, (1, 2): None, (2, 1): None}
        for key in genes:
            genes[key] = "".join(rng.choice("ACGT") for _ in range(900))
        db = self.path("db")
        os.mkdir(db)
        with open(os.path.join(db, "reference.fna"), "w") as fna, open(os.path.join(db, "reference.map"), "w") as mp:
            offset = 0
            for (taxid, gene), seq in genes.items():
                header = f">{taxid}_{gene}\n"
                fna.write(header + seq + "\n")
                mp.write(f"{taxid}\t{gene}\t{offset + len(header)}\t{offset + len(header) + len(seq)}\n")
                offset += len(header) + len(seq) + 1
        with open(os.path.join(db, "internal_taxonomy.dmp"), "w") as fh:
            fh.write("id\tparent_id\texternal_id\tname\trank\tlevel\trep_genome\n"
                     "3\t3\t0\troot\tno rank\t0\t\n"
                     "1\t3\t0\ts__Alpha one\tspecies\t7\tGCF_1\n"
                     "2\t3\t0\ts__Beta two\tspecies\t7\tGCF_2\n")
        # Another genome of taxon 2 carries taxon 1's gene 2 unchanged.
        full = self.path("full_reference.fna")
        with open(os.path.join(db, "reference.fna")) as src, open(full, "w") as dst:
            dst.write(src.read() + ">2_7\n" + genes[(1, 2)] + "\n")

        # --no_bundle: unique_kmers.tsv stays a file of its own.
        rc, log = run(self.work, "--build", "--no_bundle", "--no_profile", "-t", "1", "--db", db,
                      "--reference", os.path.join(db, "reference.fna"), "--full_reference", full)
        self.assertEqual(rc, 0, log[-3000:])
        uniques = {}
        with open(os.path.join(db, "unique_kmers.tsv")) as fh:
            for line in fh:
                f = line.split("\t")
                uniques[(int(f[0]), int(f[1]))] = (int(f[2]), int(f[4]), int(f[8]))
        self.assertGreater(uniques[(1, 1)][0], 0, uniques)
        self.assertGreater(uniques[(2, 1)][0], 0, uniques)
        self.assertGreater(uniques[(1, 2)][2], 0, uniques)
        self.assertEqual(uniques[(1, 2)][:2], (0, 0), "no k-mer of gene 1_2 is unique to taxon 1")
        for index in ("index.prx", "index.prx.zst"):  # 3 GB raw; --build compresses by default
            if os.path.exists(os.path.join(db, index)):
                os.remove(os.path.join(db, index))


class QcmsaContractTest(WorkDir):
    """qcmsa counts only the samples that are in the MSA."""

    META_HEADER = ("sample\tgene_id\tvertical_coverage\tcounts_vcov1\tcounts_vcov2\tmulti_allelic\tfiltered\t"
                   "multi_rate_vcov1\tfiltered_rate_vcov1\tmulti_rate_vcov2\tfiltered_rate_vcov2\tmedian_vcov\t"
                   "hcov\tgene_length\tmean_vcov_nonzero\tmedian_vcov_nonzero\n")

    def test_samples_missing_from_the_msa_do_not_count(self):
        with open(self.path("x.raw.msa.fna"), "w") as fh:
            fh.write(">x_reference\nACGTACGTAC\n>s1\nACGTACGTAC\n>s2\nACGTACGTAC\n")
        with open(self.path("x.raw.partition.txt"), "w") as fh:
            fh.write("DNA, gene1 = 1-10\n")
        with open(self.path("x.meta.tsv"), "w") as fh:
            fh.write(self.META_HEADER)
            for sample in ("s1", "s2", "s3", "s4"):  # s3 and s4 are not in the MSA
                fh.write(f"{sample}\t1\t5\t10\t10\t0\t0\t0\t0\t0\t0\t5\t1\t10\t5\t5\n")
        args = [QCMSA, self.path("x.raw.msa.fna"), self.path("x.raw.partition.txt"), self.path("x.meta.tsv")]

        rc, log = run(self.work, *args, "--prefix", self.path("two"), "--gene-min-samples", "2", binary="python3")
        self.assertEqual(rc, 0, log)
        self.assertIn("Loaded meta: 2 samples", log)
        self.assertFalse(os.path.exists(self.path("two.msa.fna")), "2 samples are not more than 2")

        rc, log = run(self.work, *args, "--prefix", self.path("one"), "--gene-min-samples", "1", binary="python3")
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
        return [QCMSA, self.path(name + ".raw.msa.fna"), self.path(name + ".raw.partition.txt"),
                self.path(name + ".meta.tsv")]

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
        rc, log = run(self.work, *args, "--prefix", self.path("default"), binary="python3")
        self.assertEqual(rc, 0, log)
        msa = self.read_msa(self.path("default.msa.fna"))
        self.assertEqual(msa["s1"], "ACAACAA", "every column but the one without a base; singletons kept")
        self.assertEqual(msa["y_reference"], "AAGAAAA")

        rc, log = run(self.work, *args, "--prefix", self.path("parsimony"), "--min-parsimony-samples", "2", binary="python3")
        self.assertEqual(rc, 0, log)
        # Column 2 (one sample differs) goes; column 3 stays: the reference row is not a sample.
        self.assertEqual(self.read_msa(self.path("parsimony.msa.fna"))["s1"], "AAACAA")

        rc, log = run(self.work, *args, "--prefix", self.path("variable"), "--discard-constant", binary="python3")
        self.assertEqual(rc, 0, log)
        self.assertEqual(self.read_msa(self.path("variable.msa.fna"))["s1"], "CAC", "A and R count as constant")

    def test_duplicate_names_stop_qcmsa(self):
        args = self.write_species("z", [("z_reference", "AAAAAAAA"), ("s1", "ACAAACAA"), ("s1", "AAAAAAAA")])
        rc, log = run(self.work, *args, "--prefix", self.path("dup"), binary="python3")
        self.assertNotEqual(rc, 0, log)
        self.assertIn("names 1 sequence(s) more than once (s1)", log)

    def test_coverage_gate_reads_the_msa(self):
        # The meta says every cell is fully covered; s3's row writes 2 of the gene's 8 positions, below
        # --gene-min-hcov 0.3, so its cell is gap-filled.
        args = self.write_species("c", [("c_reference", "AAAAAAAA"), ("s1", "ACAAAAAA"), ("s2", "AAAAAAAA"),
                                        ("s3", "AC------")])
        rc, log = run(self.work, *args, "--prefix", self.path("c"), binary="python3")
        self.assertEqual(rc, 0, log)
        self.assertEqual(self.read_msa(self.path("c.msa.fna"))["s3"], "--------", log)
        rc, log = run(self.work, *args, "--prefix", self.path("c2"), "--gene-min-hcov", "0.2", binary="python3")
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
        args = [QCMSA, self.path("m.raw.msa.fna"), self.path("m.raw.partition.txt"), self.path("m.meta.tsv")]
        rc, log = run(self.work, *args, "--prefix", self.path("m"), binary="python3")
        self.assertEqual(rc, 0, log)
        msa = self.read_msa(self.path("m.msa.fna"))
        self.assertNotIn("s8", msa, log)
        self.assertEqual(sorted(msa), sorted(["m_reference"] + samples[:7]), log)
        self.assertEqual(len(msa["s7"]), 32, "no gene removed or masked")
        with open(self.path("m.qcmsa_summary.tsv")) as fh:
            summary = fh.read()
        self.assertIn("sample_filtered\ts8\t4\tmulti-allelic rate 0.0100 > 0.0020", summary)

        # With a floor of 0.05%, the fence (0.0625%) decides, and s7 goes too.
        rc, log = run(self.work, *args, "--prefix", self.path("floor"), "--mrate2-min-rate", "0.0005", binary="python3")
        self.assertEqual(rc, 0, log)
        self.assertEqual(sorted(self.read_msa(self.path("floor.msa.fna"))), sorted(["m_reference"] + samples[:6]))

    def test_no_msa_leaves_no_stale_output(self):
        args = self.write_species("w", [("w_reference", "AAAAAAAA"), ("s1", "ACAAACAA"), ("s2", "AAAAAAAA")])
        rc, log = run(self.work, *args, "--prefix", self.path("w"), binary="python3")
        self.assertEqual(rc, 0, log)
        self.assertTrue(os.path.exists(self.path("w.msa.fna")), log)
        rc, log = run(self.work, *args, "--prefix", self.path("w"), "--gene-min-samples", "5", binary="python3")
        self.assertEqual(rc, 0, log)
        self.assertFalse(os.path.exists(self.path("w.msa.fna")), "an earlier run's MSA is removed")
        self.assertFalse(os.path.exists(self.path("w.partition.txt")))
        with open(self.path("w.qcmsa_summary.tsv")) as fh:
            self.assertIn("status\tno_msa\t\tevery gene was filtered", fh.read())


class MapUtilsTest(WorkDir):
    """protal_map_utils resolves relative map paths as protal does."""

    def test_relative_paths_resolve_like_protal(self):
        os.makedirs(self.path("maps"))
        os.makedirs(self.path("reads"))
        maps = {}
        for sample in ("sa", "sb"):
            for mate in (1, 2):
                shutil.copy(os.path.join(READS, f"{sample}_R{mate}.fq"), self.path("reads", f"{sample}_R{mate}.fq"))
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
        self.assertTrue(os.path.isfile(self.path("out", "aln", "sa.sam")), "protal found the reads and wrote OUTPUT_DIR/SAM_OUTPUT_DIR")

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


class LauncherTest(WorkDir):
    """The protal launcher runs the binaries installed next to it before any on $PATH."""

    def stub(self, directory, name, text):
        os.makedirs(directory, exist_ok=True)
        path = os.path.join(directory, name)
        with open(path, "w") as fh:
            fh.write(f'#!/bin/sh\necho {text} "$@"\n')
        os.chmod(path, 0o755)
        return path

    def test_own_install_wins_over_path(self):
        install, other = self.path("install"), self.path("other")
        os.makedirs(install)
        launcher = os.path.join(install, "protal")
        shutil.copy(os.path.join(ROOT, "protal_launcher"), launcher)
        os.chmod(launcher, 0o755)
        self.stub(install, "protal_baseline", "own-baseline")
        self.stub(other, "protal_avx2", "other-avx2")
        self.stub(other, "protal_baseline", "other-baseline")
        env = dict(os.environ, PATH=other + os.pathsep + os.environ["PATH"])
        env.pop("PROTAL_NO_AVX2", None)

        out = subprocess.run([launcher, "--x", "a b"], env=env, stdout=subprocess.PIPE, text=True)
        self.assertEqual(out.stdout.strip(), "own-baseline --x a b")

        # Without a binary of its own, it falls back to $PATH.
        os.remove(os.path.join(install, "protal_baseline"))
        env["PROTAL_NO_AVX2"] = "1"
        out = subprocess.run([launcher, "--x"], env=env, stdout=subprocess.PIPE, text=True)
        self.assertEqual(out.stdout.strip(), "other-baseline --x")


class RerunTest(WorkDir):
    """Existing SAM files are reused, and --no_profile is honoured when all of them exist."""

    def test_reruns(self):
        args = ["--db", DB, *reads("sa", "sb"), "-o", "out", "-t", "4", "--no_qcmsa"]
        rc, log = run(self.work, *args)
        self.assertEqual(rc, 0, log[-3000:])
        sam = glob.glob(self.path("out", "sa*.sam"))[0]
        sam_mtime = os.path.getmtime(sam)
        time.sleep(1)

        rc, log = run(self.work, *args)
        self.assertEqual(rc, 0)
        self.assertIn("All alignments are present", log)
        self.assertEqual(os.path.getmtime(sam), sam_mtime, "SAM untouched by the rerun")

        profile = glob.glob(self.path("out", "sa*.profile"))[0]
        profile_mtime = os.path.getmtime(profile)
        time.sleep(1)
        rc, log = run(self.work, *args, "--no_profile")
        self.assertEqual(rc, 0)
        self.assertEqual(os.path.getmtime(profile), profile_mtime, "--no_profile leaves the profiles alone")


class FailureTest(WorkDir):
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
        rc, log = run(self.work, "--db", DB, "-1", f"{READS}/sa_R1.fq,{self.path('bad_R1.fq')}",
                      "-2", f"{READS}/sa_R2.fq,{self.path('bad_R2.fq')}", "--prefix", "good,bad",
                      "-o", "out_bad", "-t", "2", "--no_qcmsa")
        self.assertEqual(rc, 1, log[-3000:])
        self.assertRegex(log, r"protal finished with \d+ error")
        self.assertTrue(glob.glob(self.path("out_bad", "good*.profile")), "the other sample is still profiled")

    def test_missing_qcmsa(self):
        # Two samples, so that there are strain MSAs for qcmsa to filter.
        rc, log = run(self.work, "--db", DB, *reads("sa", "sb"), "-o", "out_noqc", "-t", "2",
                      "--qcmsa_script", self.path("no", "such", "qcmsa"))
        self.assertEqual(rc, 1)
        self.assertIn("qcmsa not found", log)

    def test_truncated_index(self):
        bad_db = self.path("bad_db")
        os.mkdir(bad_db)
        index = db_file("index.prx")
        for f in glob.glob(os.path.join(FILES, "*")):
            if not os.path.basename(f).startswith(("index.prx", "database.protal")):
                os.symlink(f, os.path.join(bad_db, os.path.basename(f)))
        with open(index, "rb") as src, open(os.path.join(bad_db, os.path.basename(index)), "wb") as dst:
            dst.write(src.read(min(1 << 20, os.path.getsize(index) // 2)))
        rc, log = run(self.work, "--db", bad_db, *reads("sa"), "-o", "out_idx", "-t", "1", "--no_qcmsa")
        self.assertEqual(rc, 8)
        self.assertRegex(log, r"Invalid index .*truncated or corrupt")


def is_seekable(path):
    """True if a zstd file ends with a seek table (zstd seekable format, as protal writes)."""
    with open(path, "rb") as fh:
        fh.seek(-4, os.SEEK_END)
        return fh.read(4) == b"\xb1\xea\x92\x8f"


class CompressedDatabaseTest(WorkDir):
    """Raw, seekable (--compress_db --no_bundle, loaded in parallel), single-frame (zstd CLI) and
    single-file (--compress_db: database.protal) copies of the database give identical results,
    with one thread or several."""

    KINDS = ("raw", "seekable", "single", "bundle")
    FILES = ("index.prx", "reference.fna", "reference.map", "internal_taxonomy.dmp", "unique_kmers.tsv", "model_pe.xml")

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        if not shutil.which("zstd"):
            raise unittest.SkipTest("the zstd CLI is needed to make raw and compressed copies of the database")
        cls.dbs = {kind: os.path.join(cls.work, f"{kind}_db") for kind in cls.KINDS}
        for d in cls.dbs.values():
            os.mkdir(d)
        big = ("index.prx", "reference.fna")
        for f in glob.glob(os.path.join(FILES, "*")):
            name = os.path.basename(f)
            if not name.startswith(big) and name != "database.protal":
                for d in cls.dbs.values():
                    os.symlink(f, os.path.join(d, name))
        # A raw copy: --decompress_db on symlinks to the database's files (a compressed index is in
        # protal's column format, which zstd -d does not turn back into index.prx).
        for name in big:
            path = db_file(name)
            os.symlink(path, os.path.join(cls.dbs["raw"], os.path.basename(path)))
        cls.decompress_rc, cls.decompress_log = run(cls.work, "--decompress_db", "--db", cls.dbs["raw"], "-t", "4")
        for name in big:
            raw = os.path.join(cls.dbs["raw"], name)
            # -f: raw may be a symlink, which the zstd CLI skips otherwise.
            subprocess.run(["zstd", "-q", "-f", "-3", "--long=27", raw, "-o", os.path.join(cls.dbs["single"], name + ".zst")],
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

    def result(self, db, out, threads, *extra):
        """Sorted SAM records and profile of sample sa on db."""
        rc, log = run(self.work, "--db", db, *reads("sa"), "-o", out, "-t", str(threads), "--no_qcmsa", *extra)
        self.assertEqual(rc, 0, log[-3000:])
        with open(glob.glob(self.path(out, "sa*.sam"))[0]) as sam, open(glob.glob(self.path(out, "sa*.profile"))[0]) as prof:
            return sorted(line for line in sam if not line.startswith("@")), prof.read(), log

    def test_decompress_db(self):
        self.assertEqual(self.decompress_rc, 0, self.decompress_log[-3000:])
        for name in ("index.prx", "reference.fna"):
            self.assertTrue(os.path.isfile(os.path.join(self.dbs["raw"], name)), f"raw {name}")
            self.assertFalse(os.path.lexists(os.path.join(self.dbs["raw"], name + ".zst")))
        with open(os.path.join(self.dbs["raw"], "reference.map"), "rb") as fh:
            ends = [int(line.split()[3]) for line in fh]
        self.assertGreaterEqual(os.path.getsize(os.path.join(self.dbs["raw"], "reference.fna")), max(ends))

    def test_compress_db(self):
        self.assertEqual(self.compress_rc, 0, self.compress_log[-3000:])
        for name in ("index.prx", "reference.fna"):
            self.assertFalse(os.path.lexists(os.path.join(self.dbs["seekable"], name)), f"{name} replaced")
            self.assertTrue(is_seekable(os.path.join(self.dbs["seekable"], name + ".zst")), f"{name}.zst is seekable")
            self.assertFalse(is_seekable(os.path.join(self.dbs["single"], name + ".zst")))
            self.assertTrue(os.path.exists(os.path.join(self.dbs["raw"], name)), "the raw files stay")
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
        rc, log = run(self.work, "--compress_db", "--db", self.dbs["bundle"], "-t", "2")
        self.assertEqual(rc, 0, log[-3000:])
        self.assertIn("already a single-file database; kept", log)
        rc, log = run(self.work, "--unpack_db", "--db", self.dbs["seekable"])
        self.assertEqual(rc, 30, log[-3000:])
        self.assertIn("--unpack_db needs a single-file database", log)

    def test_round_trip_is_byte_identical(self):
        """--decompress_db of the column-format index, and of the single file, gives exactly the raw files."""
        for kind in ("seekable", "bundle"):
            db = self.path(f"round_trip_{kind}_db")
            os.mkdir(db)
            for f in glob.glob(os.path.join(self.dbs[kind], "*")):
                os.symlink(os.path.realpath(f), os.path.join(db, os.path.basename(f)))
            rc, log = run(self.work, "--decompress_db", "--db", db, "-t", "4")
            self.assertEqual(rc, 0, log[-3000:])
            self.assertFalse(os.path.lexists(os.path.join(db, "database.protal")))
            for name in self.FILES:
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
        expected = self.result(self.dbs["raw"], "out_raw_1", 1)
        for kind in self.KINDS:
            for threads in (1, 4):
                if (kind, threads) == ("raw", 1):
                    continue
                sam, profile, log = self.result(self.dbs[kind], f"out_{kind}_{threads}", threads)
                where = f"index.prx in {self.bundle}" if kind == "bundle" else os.path.join(self.dbs[kind], "index.prx")
                self.assertIn("Load index " + where, log)
                self.assertEqual(sam, expected[0], f"SAM differs: {kind} database, {threads} threads")
                self.assertEqual(profile, expected[1], f"profile differs: {kind} database, {threads} threads")
        # The single file named directly.
        sam, profile, log = self.result(self.bundle, "out_bundle_file", 4)
        self.assertEqual((sam, profile), expected[:2])

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
        expected = self.result(self.dbs["raw"], "out_lazy_expected", 2)
        with open(glob.glob(self.path("out_lazy_file", "sa*.sam"))[0]) as sam:
            self.assertEqual(sorted(line for line in sam if not line.startswith("@")), expected[0])
        with open(glob.glob(self.path("out_lazy_file", "sa*.profile"))[0]) as prof:
            self.assertEqual(prof.read(), expected[1])

class ReadTypeModelTest(WorkDir):
    """A database holds one presence model per read type (--read_type); --add_model stores one."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        # A single-file database of our own, which --add_model may rewrite.
        cls.db = os.path.join(cls.work, "db")
        os.mkdir(cls.db)
        for f in glob.glob(os.path.join(FILES, "*")):
            if os.path.basename(f) != "database.protal":
                os.symlink(f, os.path.join(cls.db, os.path.basename(f)))
        cls.pack_rc, cls.pack_log = run(cls.work, "--compress_db", "--db", cls.db, "-t", "4", "--compress_level", "3")
        cls.bundle = os.path.join(cls.db, "database.protal")
        cls.pe_rc, cls.pe_log = run(cls.work, "--db", cls.db, *reads("sa"), "-o", "out_pe", "-t", "2", "--no_qcmsa")

    def profile_only(self, out, *extra):
        sam = glob.glob(self.path("out_pe", "sa*.sam"))[0]
        return run(self.work, "--db", self.db, "--profile_only", sam, "--prefix", "sa", "-o", out, "-t", "2",
                   "--no_qcmsa", *extra)

    def test_add_model_for_a_read_type(self):
        self.assertEqual(self.pack_rc, 0, self.pack_log[-3000:])
        self.assertEqual(self.pe_rc, 0, self.pe_log[-3000:])
        self.assertIn("Model of paired-end reads: model_pe.xml in " + self.bundle, self.pe_log)
        rc, log = self.profile_only("out_pb", "--read_type", "pb")
        self.assertEqual(rc, 30, log[-3000:])
        self.assertIn("no model for --read_type pb (PacBio reads): model_PB.xml in", log)
        self.assertIn("--add_model MODEL.xml --read_type pb --db " + self.bundle, log)

        # The paired-end model stored as the single-end one: the same profile.
        rc, log = run(self.work, "--add_model", db_file("model_pe.xml"), "--read_type", "se", "--db", self.db, "-t", "2")
        self.assertEqual(rc, 0, log[-3000:])
        self.assertIn("Stored", log)
        self.assertRegex(log, r"Models for read types: pe, se\b")
        rc, log = self.profile_only("out_se", "--read_type", "se")
        self.assertEqual(rc, 0, log[-3000:])
        self.assertIn("Model of single-end reads: model_se.xml in " + self.bundle, log)
        self.assertIn("holds paired-end reads; profiled as single-end reads (se, --read_type or READ_TYPE)", log)
        with open(glob.glob(self.path("out_pe", "sa*.profile"))[0]) as a, open(glob.glob(self.path("out_se", "sa*.profile"))[0]) as b:
            self.assertEqual(a.read(), b.read())
        rc, log = run(self.work, "--unpack_db", "--db", self.bundle, "--unpack_dir", self.path("unpacked"))
        self.assertEqual(rc, 0, log[-3000:])
        self.assertTrue(filecmp.cmp(self.path("unpacked", "model_se.xml"), db_file("model_pe.xml"), shallow=False))

    def test_placeholder_model(self):
        # scripts/placeholder_models.py fills a read type's slot until a trained model replaces it: it reports
        # no species (--knob 0: every taxon with reads), and protal warns whenever it loads it.
        subprocess.run([sys.executable, os.path.join(ROOT, "scripts", "placeholder_models.py"), "-o",
                        self.path("placeholders"), "--read_types", "ont"], check=True, capture_output=True)
        placeholder = self.path("placeholders", "model_ONT.xml")
        warning = "is a placeholder, not a trained model"
        rc, log = run(self.work, "--add_model", placeholder, "--read_type", "ont", "--db", self.db, "-t", "2")
        self.assertEqual(rc, 0, log[-3000:])
        self.assertIn(warning, log)
        self.assertRegex(log, r"Models for read types: [^\n]*\bont\b")
        rc, log = self.profile_only("out_ont", "--read_type", "ont")
        self.assertEqual(rc, 0, log[-3000:])
        self.assertIn("Model of ONT reads: model_ONT.xml in " + self.bundle, log)
        self.assertIn(warning, log)
        with open(glob.glob(self.path("out_ont", "sa*.profile"))[0]) as fh:
            self.assertNotIn("s__", fh.read(), "a placeholder reports no species")
        rc, log = self.profile_only("out_ont_all", "--read_type", "ont", "--knob", "0")
        self.assertEqual(rc, 0, log[-3000:])
        with open(glob.glob(self.path("out_ont_all", "sa*.profile"))[0]) as all_taxa, \
                open(glob.glob(self.path("out_pe", "sa*.profile"))[0]) as pe:
            reported = {line.split("\t")[1] for line in all_taxa if "s__" in line}
            self.assertTrue({line.split("\t")[1] for line in pe if "s__" in line} <= reported)
        rc, log = self.profile_only("out_pe_again")
        self.assertNotIn("placeholder", log, "the paired-end model is not one")

    def test_an_unusable_model_is_not_added(self):
        with open(self.bundle, "rb") as fh:
            before = fh.read()
        bad = self.path("bad.xml")
        with open(bad, "w") as fh:
            fh.write("<PMML>\n")
        rc, log = run(self.work, "--add_model", bad, "--read_type", "ont", "--db", self.db)
        self.assertEqual(rc, 2, log[-3000:])
        self.assertIn("Cannot load the model", log)
        with open(self.bundle, "rb") as fh:
            self.assertEqual(fh.read(), before, "the database is unchanged")

    def test_read_type_checks(self):
        rc, log = run(self.work, "--db", self.db, *reads("sa"), "-o", "out_x", "--read_type", "nanopore")
        self.assertEqual(rc, 30, log[-3000:])
        self.assertIn("is 'nanopore' (--read_type or READ_TYPE): give one of pe, se, pb, ont", log)
        rc, log = run(self.work, "--db", self.db, *reads("sa"), "-o", "out_y", "--read_type", "se")
        self.assertEqual(rc, 30, log[-3000:])
        self.assertIn("has single-end reads (se), which come in one file, but it has a second read file", log)

    def test_older_database_with_model_xml(self):
        """model.xml of a database from before read types serves paired-end reads."""
        db = self.path("old_db")
        os.mkdir(db)
        for f in glob.glob(os.path.join(FILES, "*")):
            name = os.path.basename(f)
            if name.startswith("model_") or name == "database.protal":
                continue
            os.symlink(f, os.path.join(db, name))
        os.symlink(db_file("model_pe.xml"), os.path.join(db, "model.xml"))
        rc, log = run(self.work, "--db", db, *reads("sa"), "-o", "out_old", "-t", "2", "--no_qcmsa")
        self.assertEqual(rc, 0, log[-3000:])
        self.assertIn("Model of paired-end reads: " + os.path.join(db, "model.xml"), log)


class FailFastTest(WorkDir):
    """Problems with the database or the inputs stop protal before any read is aligned."""

    def db_copy(self, name, replace=None, drop=()):
        """A database of symlinks to DB, with the files in `replace` ({name: bytes}) written instead
        and the files in `drop` left out."""
        replace = replace or {}
        db = self.path(name)
        os.mkdir(db)
        for f in glob.glob(os.path.join(FILES, "*")):
            base = os.path.basename(f)
            if base not in replace and base not in drop and base != "database.protal":
                os.symlink(f, os.path.join(db, base))
        for base, content in replace.items():
            with open(os.path.join(db, base), "wb") as fh:
                fh.write(content)
        return db

    @staticmethod
    def db_file(name):
        """Content of a database file (decompressed if the database holds <name>.zst)."""
        path = db_file(name)
        if path.endswith(".zst"):
            return subprocess.run(["zstd", "-dc", path], check=True, stdout=subprocess.PIPE).stdout
        with open(path, "rb") as fh:
            return fh.read()

    def query(self, db, out, *extra, samples=("sa",)):
        rc, log = run(self.work, "--db", db, *reads(*samples), "-o", out, "-t", "1", "--no_qcmsa", *extra)
        self.assertFalse(glob.glob(self.path(out, "*.sam*")), "no read may be aligned")
        return rc, log

    def test_reference_changed_since_build(self):
        db = self.db_copy("db_ref", {"reference.fna": self.db_file("reference.fna") + b">9_1\nACGT\n"})
        rc, log = self.query(db, "out_ref")
        if "no reference fingerprint" in log:
            self.skipTest("the test database's index predates the reference fingerprint")
        self.assertEqual(rc, 8, log[-3000:])
        self.assertIn("index.prx was built against a different reference", log)

    def test_malformed_map(self):
        db = self.db_copy("db_map", {"reference.map": self.db_file("reference.map") + b"1\t999\t5\n"})
        rc, log = self.query(db, "out_map")
        self.assertEqual(rc, 8, log[-3000:])
        self.assertRegex(log, r"Invalid reference map .*expected 4 tab-separated columns, found 3")

    def test_missing_db(self):
        """A --db path with nothing at it is reported as missing, relative to where protal runs."""
        for args in (reads("sa") + ["-o", "out_nodb", "--no_qcmsa"], ["--unpack_db"], ["--compress_db"]):
            rc, log = run(self.work, "--db", "no/such_db", *args)
            self.assertEqual(rc, 30, log[-3000:])
            self.assertIn("--db no/such_db does not exist (relative to the working directory "
                          f"{os.path.realpath(self.work)})", log)
            self.assertNotIn("holds separate files", log)
            self.assertNotIn("Sequence file does not exist", log)

    def test_missing_model_and_unique_kmers(self):
        db = self.db_copy("db_files", drop=("model_pe.xml", "unique_kmers.tsv"))
        rc, log = self.query(db, "out_files")
        self.assertEqual(rc, 30, log[-3000:])
        self.assertIn("The database has no model for --read_type pe", log)
        self.assertIn("Unique k-mer file does not exist", log)

    def test_corrupt_model(self):
        db = self.db_copy("db_model", {"model_pe.xml": b"<PMML>\n"})
        rc, log = self.query(db, "out_model")
        self.assertEqual(rc, 2, log[-3000:])
        self.assertIn("Cannot load the model", log)

    def test_model_protal_cannot_feed(self):
        model = self.db_file("model_pe.xml").decode()
        # An input protal does not compute, and a model predicting other labels than TRUE/FALSE.
        unknown = model.replace("<MiningSchema>", '<MiningSchema>\n<MiningField name="moon_phase"/>', 1)
        unknown = re.sub(r"(<DataDictionary[^>]*>)", r'\1\n<DataField name="moon_phase" optype="continuous" dataType="double"/>',
                         unknown, count=1)
        labels = model.replace('value="TRUE"', 'value="present"').replace('score="TRUE"', 'score="present"')
        for name, text, message in (("db_unknown", unknown, "input(s) protal does not compute: moon_phase"),
                                    ("db_labels", labels, "has no value TRUE")):
            db = self.db_copy(name, {"model_pe.xml": text.encode()})
            rc, log = self.query(db, "out_" + name)
            self.assertEqual(rc, 2, log[-3000:])
            self.assertIn("Cannot use the model", log)
            self.assertIn(message, log)

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
        self.assertIn("must name one file per sample: 1 given for 2 samples", log)

    def test_map_row_without_a_profile_cell(self):
        sample_map = self.path("samples.map")
        with open(sample_map, "w") as fh:
            fh.write(f"#OUTPUT_DIR\t{self.path('out_map_rows')}\n#SAMPLEID\tPREFIX\tFIRST\tSECOND\tPROFILE\n")
            fh.write(f"sa\tsa\t{READS}/sa_R1.fq\t{READS}/sa_R2.fq\tsa.profile\n")
            fh.write(f"sb\tsb\t{READS}/sb_R1.fq\t{READS}/sb_R2.fq\n")
        rc, log = run(self.work, "--db", DB, "--map", sample_map, "-t", "1", "--no_qcmsa")
        self.assertEqual(rc, 9, log[-3000:])
        self.assertIn("Line 4: no value in column 5 (PROFILE)", log)
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
                fh.write(self.db_file(f))
        bad = self.path("bad_reference.fna")
        with open(bad, "wb") as fh:
            fh.write(self.db_file("reference.fna") + b">9_1\nACGTACGT\n>unnamed\nACGTACGT\n")
        rc, log = run(self.work, "--build", "--no_profile", "-t", "1", "--db", db,
                      "--reference", bad, "--full_reference", bad)
        self.assertEqual(rc, 8, log[-3000:])
        self.assertIn("--reference record >9_1: not in", log)
        self.assertIn("--reference record >unnamed: header is not <taxid>_<gene id>", log)
        self.assertRegex(log, r"2 of \d+ --reference records do not match")
        for written in ("index.prx", "index.prx.zst", "database.protal"):
            self.assertFalse(os.path.exists(os.path.join(db, written)), written)


class SamInputTest(WorkDir):
    """--profile_only reads SAM files as other tools may leave them."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        rc, log = run(cls.work, "--db", DB, *reads("sa"), "-o", "out", "-t", "2", "--no_qcmsa")
        assert rc == 0, log[-3000:]
        cls.sam = glob.glob(os.path.join(cls.work, "out", "sa*.sam"))[0]
        with open(cls.sam) as fh:
            lines = fh.read().splitlines()
        cls.header = [line for line in lines if line.startswith("@")]
        cls.records = [line for line in lines if not line.startswith("@")]
        with open(os.path.join(cls.work, "out", "sa.profile")) as fh:
            cls.profile_text = fh.read()
        assert cls.profile_text.strip(), "the reference profile lists taxa"

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
        self.assertIn("skipped 1 record(s): unmapped", log)
        self.assertIn("skipped 1 record(s): reference is not a protal gene", log)
        with open(self.path("out_edited.sam", "edited.profile")) as fh:
            self.assertEqual(fh.read(), self.profile_text)

    def test_header_only_sam_gets_an_empty_profile(self):
        sam = self.write_sam("header_only", self.header)
        rc, log = self.profile_only(sam)
        self.assertEqual(rc, 0, log[-3000:])
        self.assertIn("contains no usable alignments", log)
        self.assertTrue(os.path.isfile(self.path("out_header_only.sam", "header_only.profile")))

    def test_unreadable_sam_fails_only_its_sample(self):
        good = self.write_sam("good", self.header + self.records)
        broken = self.write_sam("broken", self.header + self.records[:50] + ["sa.9\t0\t1_1"])
        rc, log = self.profile_only(good, broken)
        self.assertEqual(rc, 1, log[-3000:])
        self.assertRegex(log, r"Cannot read the SAM file of sample \S+ \(.*broken\.sam\): line \d+: expected at least 11")
        with open(self.path("out_good.sam", "good.profile")) as fh:
            self.assertEqual(fh.read(), self.profile_text)

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
        self.assertEqual(rc, 1, log[-3000:])
        self.assertIn(f"gene {name} is {int(length) + 7} bp in the SAM header (@SQ) but {length} bp in the database", log)

    def test_samples_cannot_share_an_output_file(self):
        sam = self.write_sam("twice", self.header + self.records)
        rc, log = self.profile_only(sam, sam)
        self.assertNotEqual(rc, 0, log[-3000:])
        self.assertIn("samples 1 and 2 would both use the SAM file", log)
        self.assertFalse(glob.glob(self.path("out_twice.sam", "*.profile")))


class QcmsaTest(WorkDir):
    """The post-filter runs, and --qcmsa_args reaches it intact."""

    def test_no_filtered_msa_is_reported(self):
        rc, log = run(self.work, "--db", DB, *reads("sa", "sb"), "-o", "out_none", "-t", "4",
                      "--qcmsa_script", QCMSA, "--qcmsa_args", "--gene-min-samples 100")
        self.assertEqual(rc, 0, log[-3000:])
        self.assertIn("qcmsa kept no gene or sample", log)
        filtered = [f for f in glob.glob(self.path("out_none", "strains", "*.msa.fna")) if not f.endswith(".raw.msa.fna")]
        self.assertEqual(filtered, [])

    def test_filtered_msa(self):
        # With its defaults, qcmsa filters an MSA of three samples (a gene needs two).
        rc, log = run(self.work, "--db", DB, *reads("sa", "sb", "sr"), "-o", "out", "-t", "4", "--qcmsa_script", QCMSA)
        self.assertEqual(rc, 0, log[-3000:])
        filtered = [f for f in glob.glob(self.path("out", "strains", "*.msa.fna")) if not f.endswith(".raw.msa.fna")]
        self.assertTrue(filtered, "qcmsa wrote filtered MSAs")
        with open(self.path("out", "strains", "species.tsv")) as fh:
            listed = [line.rstrip("\n").split("\t") for line in fh][1:]
        self.assertEqual(sorted(row[4] for row in listed if row[4] != "-"), sorted(os.path.basename(f) for f in filtered))

        # A rerun in which qcmsa keeps nothing leaves no filtered MSA of the first run behind.
        rc, log = run(self.work, "--db", DB, *reads("sa", "sb", "sr"), "-o", "out", "-t", "4",
                      "--qcmsa_script", QCMSA, "--qcmsa_args", "--gene-min-samples 100")
        self.assertEqual(rc, 0, log[-3000:])
        self.assertFalse([f for f in filtered if os.path.exists(f)], "stale filtered MSAs")
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


class SingleEndTest(WorkDir):
    """Single-end reads (the first mates of the simulated samples), profiled with the single-end model.
    The test database has none, so a copy of it gets its paired-end model as model_se.xml: the
    profiles then test the plumbing, not the model."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.db = os.path.join(cls.work, "se_db")
        os.mkdir(cls.db)
        for f in glob.glob(os.path.join(FILES, "*")):
            if os.path.basename(f) != "database.protal":
                os.symlink(f, os.path.join(cls.db, os.path.basename(f)))
        cls.model = db_file("model_pe.xml")
        os.symlink(cls.model, os.path.join(cls.db, "model_se.xml"))
        cls.rc, cls.log = run(cls.work, "--db", cls.db, *single_reads("sa", "sb"), "-o", "out", "-t", "4", "--no_qcmsa")

    def reads_by_name(self, prefix):
        with open(os.path.join(READS, f"{prefix}_R1.fq")) as fh:
            lines = fh.read().splitlines()
        return {lines[i][1:].rsplit("/", 1)[0]: lines[i + 1] for i in range(0, len(lines), 4)}

    def test_exit_code_and_model(self):
        self.assertEqual(self.rc, 0, self.log[-3000:])
        self.assertIn("Align the single-end reads of sample sa", self.log)
        self.assertIn("Model of single-end reads: " + os.path.join(self.db, "model_se.xml"), self.log)
        self.assertNotIn("Model of paired-end reads", self.log)

    def test_records_are_unpaired_reads(self):
        records = sam_records(self.path("out", "sa.sam"))
        self.assertTrue(records)
        self.assertEqual([r[1] for r in records if int(r[1]) & 0xCD], [], "no pair, mate or unmapped flags")
        self.assertEqual(len({r[0] for r in records}), len(records), "one record per read (-m 1)")
        reads = self.reads_by_name("sa")
        for r in records:
            read = reads[r[0]]  # QNAME is the read id without its /1
            self.assertEqual(r[9], revcomp(read) if int(r[1]) & 0x10 else read, "SEQ in reference orientation")
            self.assertEqual((r[6], r[7], r[8]), ("*", "0", "0"))
        self.assertTrue(any(int(r[1]) & 0x10 for r in records) and any(not int(r[1]) & 0x10 for r in records))
        self.assertTrue(any(int(r[4]) >= 4 for r in records), "reads with a MAPQ the profiler takes")

    def test_profiles_and_strains(self):
        header, rows = read_table(self.path("out", "sa.profile.log"))
        self.assertTrue(rows)
        self.assertTrue(any(row[0] == "1" for row in rows), "a species passes the model")
        self.assertTrue(glob.glob(self.path("out", "strains", "*.raw.msa.fna")), "species in both samples get MSAs")
        self.assertTrue(os.path.exists(self.path("out", "misc", "sa_runtime.tsv")))

    def test_profile_only_takes_the_model_of_the_sams_reads(self):
        sam = self.path("out", "sa.sam")
        # The test database has no model_se.xml: the SAM's unpaired records ask for one.
        rc, log = run(self.work, "--db", DB, "--profile_only", sam, "-o", self.path("po_missing"), "-t", "2", "--no_qcmsa")
        self.assertEqual(rc, 30, log[-3000:])
        self.assertIn("no model for --read_type se (single-end reads)", log)
        rc, log = run(self.work, "--db", DB, "--profile_only", sam, "--model_se", self.model, "-o", self.path("po"),
                      "-t", "2", "--no_qcmsa")
        self.assertEqual(rc, 0, log[-3000:])
        self.assertIn("Model of single-end reads: " + self.model, log)
        with open(self.path("po", "sa.profile")) as again, open(self.path("out", "sa.profile")) as first:
            self.assertEqual(again.read(), first.read())

    def test_a_missing_single_end_model_stops_before_aligning(self):
        rc, log = run(self.work, "--db", DB, *single_reads("sa"), "-o", "out_nomodel", "-t", "1", "--no_qcmsa")
        self.assertEqual(rc, 30, log[-3000:])
        self.assertIn("no model for --read_type se (single-end reads)", log)
        self.assertIn("--model_se", log)
        self.assertFalse(glob.glob(self.path("out_nomodel", "*.sam*")))

    def test_a_map_mixes_paired_and_single_end_samples(self):
        sample_map = self.path("mixed.map")
        with open(sample_map, "w") as fh:
            fh.write(f"#OUTPUT_DIR\t{self.path('out_mixed')}\n#INPUT_DIR\t{READS}\n#SAMPLEID\tPREFIX\tFIRST\tSECOND\n")
            fh.write("pe\tpe\tsa_R1.fq\tsa_R2.fq\nse\tse\tsb_R1.fq\t-\n")
        rc, log = run(self.work, "--db", self.db, "--map", sample_map, "-t", "2", "--no_qcmsa")
        self.assertEqual(rc, 0, log[-3000:])
        self.assertIn("Model of paired-end reads: ", log)
        self.assertIn("Model of single-end reads: ", log)
        self.assertIn("1 paired-end, 1 single-end, 0 PacBio, 0 ONT sample(s)", log)
        paired = sam_records(self.path("out_mixed", "pe.sam"))
        single = sam_records(self.path("out_mixed", "se.sam"))
        self.assertTrue(paired and all(int(r[1]) & 0x1 for r in paired))
        self.assertTrue(single and not any(int(r[1]) & 0x1 for r in single))
        self.assertTrue(os.path.exists(self.path("out_mixed", "se.profile")))

    def test_fasta_reads_get_q30(self):
        fasta = self.path("sa.fa")
        with open(fasta, "w") as fh:
            for name, seq in self.reads_by_name("sa").items():
                fh.write(f">{name}\n{seq}\n")
        rc, log = run(self.work, "--db", self.db, "-1", fasta, "--prefix", "safa", "-o", "out_fasta", "-t", "2", "--no_qcmsa")
        self.assertEqual(rc, 0, log[-3000:])
        records = sam_records(self.path("out_fasta", "safa.sam"))
        self.assertTrue(records)
        self.assertEqual({q for r in records for q in r[10]}, {"?"}, "Q30 for every base")
        _, rows = read_table(self.path("out_fasta", "safa.profile.log"))
        self.assertTrue(rows, "the profiler takes the records")

    def test_the_prefix_comes_from_the_read_file(self):
        rc, log = run(self.work, "--db", self.db, "-1", os.path.join(READS, "sa_R1.fq"), "-o", "out_prefix", "-t", "1",
                      "--no_profile")
        self.assertEqual(rc, 0, log[-3000:])
        self.assertTrue(os.path.exists(self.path("out_prefix", "sa_R1.sam")))

    def test_a_single_file_database_holds_model_se(self):
        if not os.path.exists(os.path.join(FILES, "index.prx.zst")):
            self.skipTest("packing a raw index into database.protal takes long")
        db = self.path("se_bundle")
        os.mkdir(db)
        for f in glob.glob(os.path.join(self.db, "*")):
            os.symlink(os.path.realpath(f), os.path.join(db, os.path.basename(f)))
        rc, log = run(self.work, "--compress_db", "--db", db, "-t", "2", "--compress_level", "3")
        self.assertEqual(rc, 0, log[-3000:])
        self.assertFalse(os.path.lexists(os.path.join(db, "model_se.xml")), "model_se.xml is packed")
        bundle = os.path.join(db, "database.protal")
        rc, log = run(self.work, "--db", bundle, *single_reads("sa"), "-o", "out_bundle", "-t", "2", "--no_qcmsa")
        self.assertEqual(rc, 0, log[-3000:])
        self.assertIn(f"Model of single-end reads: model_se.xml in {bundle}", log)
        with open(self.path("out_bundle", "sa.profile")) as packed, open(self.path("out", "sa.profile")) as files:
            self.assertEqual(packed.read(), files.read())


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


class PacBioTest(WorkDir):
    """PacBio-like long reads of several genes each, profiled with the PacBio model. The test database
    has none, so a copy of it gets its paired-end model as model_PB.xml: the profiles test the
    plumbing, not the model."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.db = os.path.join(cls.work, "pacbio_db")
        os.mkdir(cls.db)
        for f in glob.glob(os.path.join(FILES, "*")):
            if os.path.basename(f) != "database.protal":
                os.symlink(f, os.path.join(cls.db, os.path.basename(f)))
        cls.model = db_file("model_pe.xml")
        os.symlink(cls.model, os.path.join(cls.db, "model_PB.xml"))
        cls.reads = {p: os.path.join(cls.work, f"{p}.fq") for p in ("la", "lb")}
        cls.truth = {"la": simulate_long_reads(cls.reads["la"], 40, seed=11, long_read=150000),
                     "lb": simulate_long_reads(cls.reads["lb"], 40, seed=12)}
        cls.rc, cls.log = run(cls.work, "--db", cls.db, "-1", ",".join(cls.reads.values()), "--prefix", "la,lb",
                              "--read_type", "pb", "-o", "out", "-t", "4", "--no_qcmsa")

    def read_seqs(self, prefix):
        with open(self.reads[prefix]) as fh:
            lines = fh.read().splitlines()
        return [lines[i + 1] for i in range(0, len(lines), 4)], [lines[i][1:] for i in range(0, len(lines), 4)]

    def test_exit_code_and_model(self):
        self.assertEqual(self.rc, 0, self.log[-3000:])
        self.assertIn("Align the PacBio reads of sample la", self.log)
        self.assertIn("Model of PacBio reads: " + os.path.join(self.db, "model_PB.xml"), self.log)
        self.assertIn("1 read(s) longer than 65000 bp were seeded in chunks", self.log)
        self.assertRegex(self.log, r"\d+ of \d+ gene hits that fit several taxa \(MAPQ < 4\) were settled by their read's other genes")
        with open(self.path("out", "la.sam")) as fh:
            self.assertIn("@CO\tprotal read type: pb\n", fh.read())

    def test_records_hold_their_aligned_bases(self):
        seqs, names = self.read_seqs("la")
        read_of = dict(zip(names, seqs))
        records = sam_records(self.path("out", "la.sam"))
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
            by_read = representative_records(sam_records(self.path("out", f"{prefix}.sam")))
            found = missed = twice = 0
            for name, placed in zip(names, self.truth[prefix]):
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
        long_name, placed = names[-1], self.truth["la"][-1]
        reps = representative_records(sam_records(self.path("out", "la.sam"))).get(long_name, [])
        found = [g for g, _, _ in placed if any(r[2] == g for r in reps)]
        self.assertGreaterEqual(len(found), 0.95 * len(placed))
        self.assertLessEqual(len(reps), len(placed), "each gene of the 150 kb read at most once")

    def test_profiles(self):
        header, rows = read_table(self.path("out", "la.profile.log"))
        self.assertTrue(rows)
        self.assertTrue(os.path.exists(self.path("out", "misc", "la_runtime.tsv")))

    def test_profile_only_takes_the_pacbio_model(self):
        sam = self.path("out", "la.sam")
        rc, log = run(self.work, "--db", DB, "--profile_only", sam, "-o", self.path("po_missing"), "-t", "2", "--no_qcmsa")
        self.assertEqual(rc, 30, log[-3000:])
        self.assertIn("no model for --read_type pb (PacBio reads)", log)
        self.assertIn("--model_pb", log)
        rc, log = run(self.work, "--db", DB, "--profile_only", sam, "--model_pb", self.model, "-o", self.path("po"),
                      "-t", "2", "--no_qcmsa")
        self.assertEqual(rc, 0, log[-3000:])
        with open(self.path("po", "la.profile")) as again, open(self.path("out", "la.profile")) as first:
            self.assertEqual(again.read(), first.read())

    def test_long_reads_given_as_short_ones_stop(self):
        rc, log = run(self.work, "--db", self.db, "-1", self.reads["lb"], "-o", "out_short", "-t", "1", "--no_qcmsa")
        self.assertEqual(rc, 30, log[-3000:])
        self.assertIn("too long for short reads: give --read_type pb", log)

    def test_fasta_reads_get_q30(self):
        seqs, names = self.read_seqs("lb")
        fasta = self.path("lb.fa")
        with open(fasta, "w") as fh:
            for name, seq in zip(names, seqs):
                fh.write(f">{name}\n{seq}\n")
        rc, log = run(self.work, "--db", self.db, "-1", fasta, "--read_type", "pb", "--prefix", "lbfa", "-o", "out_fasta",
                      "-t", "2", "--no_qcmsa")
        self.assertEqual(rc, 0, log[-3000:])
        records = sam_records(self.path("out_fasta", "lbfa.sam"))
        self.assertTrue(records)
        self.assertEqual({q for r in records for q in r[10]}, {"?"}, "Q30 for every base")
        _, rows = read_table(self.path("out_fasta", "lbfa.profile.log"))
        self.assertTrue(rows, "the profiler takes the records")

    def test_a_map_names_the_read_type(self):
        sample_map = self.path("typed.map")
        with open(sample_map, "w") as fh:
            fh.write(f"#OUTPUT_DIR\t{self.path('out_map')}\n#SAMPLEID\tPREFIX\tFIRST\tSECOND\tREAD_TYPE\n")
            fh.write(f"pe\tpe\t{READS}/sa_R1.fq\t{READS}/sa_R2.fq\t-\n")
            fh.write(f"lb\tlb\t{self.reads['lb']}\t-\tPB\n")
        rc, log = run(self.work, "--db", self.db, "--map", sample_map, "-t", "2", "--no_qcmsa")
        self.assertEqual(rc, 0, log[-3000:])
        self.assertIn("1 paired-end, 0 single-end, 1 PacBio, 0 ONT sample(s)", log)
        self.assertTrue(all(int(r[1]) & 0x1 for r in sam_records(self.path("out_map", "pe.sam"))))
        with open(self.path("out_map", "lb.sam")) as fh:
            self.assertIn("@CO\tprotal read type: pb\n", fh.read())


class OntTest(WorkDir):
    """ONT-like long reads (2% errors, most of them 1 bp indels, Q17) profiled with the ONT model, as
    in PacBioTest a stand-in: the test database's paired-end model as model_ONT.xml."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.db = os.path.join(cls.work, "ont_db")
        os.mkdir(cls.db)
        for f in glob.glob(os.path.join(FILES, "*")):
            if os.path.basename(f) != "database.protal":
                os.symlink(f, os.path.join(cls.db, os.path.basename(f)))
        cls.model = db_file("model_pe.xml")
        os.symlink(cls.model, os.path.join(cls.db, "model_ONT.xml"))
        cls.reads = {p: os.path.join(cls.work, f"{p}.fq") for p in ("oa", "ob")}
        ont = dict(error=0.02, indels=0.7, read_name="ont_read_{}", quality="2")
        cls.truth = {"oa": simulate_long_reads(cls.reads["oa"], 40, seed=21, long_read=150000, **ont),
                     "ob": simulate_long_reads(cls.reads["ob"], 40, seed=22, **ont)}
        cls.rc, cls.log = run(cls.work, "--db", cls.db, "-1", ",".join(cls.reads.values()), "--prefix", "oa,ob",
                              "--read_type", "ont", "-o", "out", "-t", "4", "--no_qcmsa")

    def read_names(self, prefix):
        with open(self.reads[prefix]) as fh:
            return [line[1:] for i, line in enumerate(fh.read().splitlines()) if i % 4 == 0]

    def test_exit_code_and_model(self):
        self.assertEqual(self.rc, 0, self.log[-3000:])
        self.assertIn("Align the ONT reads of sample oa (-a 0.85)", self.log)
        self.assertIn("Model of ONT reads: " + os.path.join(self.db, "model_ONT.xml"), self.log)
        self.assertIn("1 read(s) longer than 65000 bp were seeded in chunks", self.log)
        with open(self.path("out", "oa.sam")) as fh:
            self.assertIn("@CO\tprotal read type: ont\n", fh.read())

    def test_every_gene_is_found_once(self):
        for prefix in ("oa", "ob"):
            by_read = representative_records(sam_records(self.path("out", f"{prefix}.sam")))
            found = missed = twice = 0
            for name, placed in zip(self.read_names(prefix), self.truth[prefix]):
                reps = by_read.get(name, [])
                for gene, start, end in placed:
                    hits = [r for r in reps if r[2] == gene and sum(n for n, op in cigar_ops(r[5]) if op in "MX=D") >= 0.9 * (end - start)]
                    found += len(hits) >= 1
                    missed += not hits
                    twice += len(hits) > sum(1 for g, _, _ in placed if g == gene)
            self.assertEqual(twice, 0, f"{prefix}: no gene counted twice")
            self.assertGreaterEqual(found / (found + missed), 0.9, f"{prefix}: {found} genes found, {missed} missed")

    def test_profile_only_takes_the_ont_model(self):
        sam = self.path("out", "oa.sam")
        rc, log = run(self.work, "--db", DB, "--profile_only", sam, "-o", self.path("po_missing"), "-t", "2", "--no_qcmsa")
        self.assertEqual(rc, 30, log[-3000:])
        self.assertIn("no model for --read_type ont (ONT reads)", log)
        self.assertIn("--model_ont", log)
        rc, log = run(self.work, "--db", DB, "--profile_only", sam, "--model_ont", self.model, "-o", self.path("po"),
                      "-t", "2", "--no_qcmsa")
        self.assertEqual(rc, 0, log[-3000:])
        with open(self.path("po", "oa.profile")) as again, open(self.path("out", "oa.profile")) as first:
            self.assertEqual(again.read(), first.read())

    def test_fasta_reads_get_q18_and_a_given_identity(self):
        fasta = self.path("ob.fa")
        with open(self.reads["ob"]) as fq, open(fasta, "w") as fh:
            lines = fq.read().splitlines()
            for i in range(0, len(lines), 4):
                fh.write(f">{lines[i][1:]}\n{lines[i + 1]}\n")
        rc, log = run(self.work, "--db", self.db, "-1", fasta, "--read_type", "ont", "-a", "0.8", "--prefix", "obfa",
                      "-o", "out_fasta", "-t", "2", "--no_qcmsa")
        self.assertEqual(rc, 0, log[-3000:])
        self.assertIn("Align the ONT reads of sample obfa (-a 0.8)", log)
        records = sam_records(self.path("out_fasta", "obfa.sam"))
        self.assertTrue(records)
        self.assertEqual({q for r in records for q in r[10]}, {"3"}, "Q18 for every base")

    def test_a_map_mixes_paired_end_and_ont_samples(self):
        sample_map = self.path("typed.map")
        with open(sample_map, "w") as fh:
            fh.write(f"#OUTPUT_DIR\t{self.path('out_map')}\n#SAMPLEID\tPREFIX\tFIRST\tSECOND\tREAD_TYPE\n")
            fh.write(f"pe\tpe\t{READS}/sa_R1.fq\t{READS}/sa_R2.fq\t-\n")
            fh.write(f"ob\tob\t{self.reads['ob']}\t-\tont\n")
        rc, log = run(self.work, "--db", self.db, "--map", sample_map, "-t", "2", "--no_qcmsa")
        self.assertEqual(rc, 0, log[-3000:])
        self.assertIn("1 paired-end, 0 single-end, 0 PacBio, 1 ONT sample(s)", log)
        self.assertIn("Align the paired-end reads of sample pe (-a 0.9)", log)
        self.assertIn("Align the ONT reads of sample ob (-a 0.85)", log)
        self.assertNotIn("minimum allele frequencies for", log)
        with open(self.path("out_map", "ob.sam")) as fh:
            self.assertIn("@CO\tprotal read type: ont\n", fh.read())


class SimulatorTest(WorkDir):
    def test_too_few_read_pairs_fails_fast(self):
        if not os.access(SIMULATE, os.X_OK):
            self.skipTest(f"simulate_metagenomes not found at {SIMULATE}")
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
        if not os.access(SIMULATE, os.X_OK):
            self.skipTest(f"simulate_metagenomes not found at {SIMULATE}")
        with open(self.path("genomes.tsv"), "w") as table:
            for domain, count in (("A", 6), ("B", 30)):
                for sp in range(count):
                    fasta = self.path(f"{domain}{sp}.fa")
                    with open(fasta, "w") as fh:
                        fh.write(">c1\n" + "".join("ACGT"[(i * 7 + sp) % 4] for i in range(3000)) + "\n")
                    table.write(f"{domain}{sp}\td__{domain};p__P{domain};c__C{domain};o__O{domain};f__F{domain};"
                                f"g__G{domain};s__G{domain} sp{sp}\t{fasta}\n")
        rc, log = run(self.work, "--genome_table", "genomes.tsv", "--test", "--seed", "1", "--samples", "8",
                      "--total_read_pairs", "1000", "--species_per_sample", "6", "--taxon", "d__A:2",
                      "--output_dir", "sim", binary=SIMULATE, timeout=60)
        self.assertEqual(rc, 0, log)
        domains = {}
        with open(self.path("sim", "manifest.tsv")) as fh:
            header = next(fh).rstrip("\n").split("\t")
            for line in fh:
                row = dict(zip(header, line.rstrip("\n").split("\t")))
                domains.setdefault(row["sample"], []).append(row["taxonomy"].split(";")[0])
        self.assertEqual(len(domains), 8)
        for sample, found in domains.items():
            self.assertEqual(len(found), 6, sample)
            self.assertGreaterEqual(found.count("d__A"), 2, sample)
        # The four species drawn at random are d__A 4 times in 34, so about 2.5 d__A per sample, not 6.
        self.assertLess(sum(f.count("d__A") for f in domains.values()) / len(domains), 3.5)


if __name__ == "__main__":
    sys.exit(unittest.main())
