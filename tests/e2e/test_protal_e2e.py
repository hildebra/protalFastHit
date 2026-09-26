#!/usr/bin/env python3
"""test_protal_e2e.py - end-to-end checks of protal and simulate_metagenomes on a small database.

Simulates paired reads from the database's reference genes (sequencing errors only, fixed seed)
and runs the real binaries, checking what the unit tests cannot: exit codes, output files, SAM
records, strain MSAs, reruns, and that failures are reported.

  PROTAL_TEST_DB=data/mini_db/protal_db python3 -m unittest -v tests/e2e/test_protal_e2e.py
  just e2e                                  # builds the mini DB first

PROTAL_TEST_DB  protal database, e.g. from scripts/mini_db/build_mini_db.sh (required; only read).
                Raw or zstd-compressed (index.prx.zst, reference.fna.zst); the zstd CLI is
                needed for a compressed one.
PROTAL          protal binary (default: build/protal)
SIMULATE        simulate_metagenomes binary (default: build/simulate_metagenomes; optional)
"""

import glob
import os
import random
import re
import shutil
import subprocess
import sys
import tempfile
import time
import unittest

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
DB = os.environ.get("PROTAL_TEST_DB", "")
PROTAL = os.path.abspath(os.environ.get("PROTAL", os.path.join(ROOT, "build", "protal")))
SIMULATE = os.path.abspath(os.environ.get("SIMULATE", os.path.join(ROOT, "build", "simulate_metagenomes")))
QCMSA = os.path.join(ROOT, "scripts", "qcmsa.py")
READS = None  # directory with the simulated reads, set up once per module


def db_file(name):
    """Path of a database file as protal picks it: <name> if present, else <name>.zst."""
    raw = os.path.join(DB, name)
    return raw if os.path.exists(raw) or not os.path.exists(raw + ".zst") else raw + ".zst"


def setUpModule():
    global READS
    index = db_file("index.prx")
    if not DB or not os.path.isfile(index) or os.path.getsize(index) == 0:
        raise unittest.SkipTest("set PROTAL_TEST_DB to a protal database (just mini-db builds data/mini_db/protal_db)")
    if db_file("reference.fna").endswith(".zst") and not shutil.which("zstd"):
        raise unittest.SkipTest("the zstd CLI is needed to read reference.fna.zst")
    if not os.access(PROTAL, os.X_OK):
        raise unittest.SkipTest(f"protal binary not found at {PROTAL} (set PROTAL)")
    READS = tempfile.mkdtemp(prefix="protal_e2e_reads_")
    simulate_reads("sa", pairs_per_gene=12, seed=1)
    simulate_reads("sb", pairs_per_gene=12, seed=2)
    simulate_reads("sr", pairs_per_gene=12, seed=3, random_r2=True)


def tearDownModule():
    if READS:
        shutil.rmtree(READS, ignore_errors=True)


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


class StrainEdgeCaseTest(WorkDir):
    def test_species_without_msa_genes(self):
        # No gene reaches --snp_min_cov, so no species has MSA columns (this used to segfault).
        # Two samples: MSAs are built only across samples.
        rc, log = run(self.work, "--db", DB, *reads("sa", "sb"), "-o", "out", "-t", "2", "--no_qcmsa",
                      "--snp_min_cov", "100000", "--msa_min_hcov", "0")
        self.assertEqual(rc, 0, log[-3000:])
        self.assertIn("has enough coverage for an MSA", log)
        self.assertEqual(glob.glob(self.path("out", "strains", "*.raw.msa.fna")), [])


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
        with open(os.path.join(DB, "model.xml")) as fh:
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
        for f in glob.glob(os.path.join(DB, "*")):
            if not os.path.basename(f).startswith("index.prx"):
                os.symlink(f, os.path.join(bad_db, os.path.basename(f)))
        with open(index, "rb") as src, open(os.path.join(bad_db, os.path.basename(index)), "wb") as dst:
            dst.write(src.read(min(1 << 20, os.path.getsize(index) // 2)))
        rc, log = run(self.work, "--db", bad_db, *reads("sa"), "-o", "out_idx", "-t", "1", "--no_qcmsa")
        self.assertEqual(rc, 8)
        self.assertRegex(log, r"Invalid index .*truncated or corrupt")


class CompressedDatabaseTest(WorkDir):
    """A raw and a zstd-compressed copy of the database give byte-identical results."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        if not shutil.which("zstd"):
            raise unittest.SkipTest("the zstd CLI is needed to make a raw and a compressed copy of the database")
        cls.dbs = {"raw": os.path.join(cls.work, "raw_db"), "zst": os.path.join(cls.work, "zst_db")}
        for d in cls.dbs.values():
            os.mkdir(d)
        big = ("index.prx", "reference.fna")
        for f in glob.glob(os.path.join(DB, "*")):
            name = os.path.basename(f)
            if not name.startswith(big):
                for d in cls.dbs.values():
                    os.symlink(f, os.path.join(d, name))
        for name in big:
            path = db_file(name)
            if path.endswith(".zst"):
                os.symlink(path, os.path.join(cls.dbs["zst"], name + ".zst"))
                subprocess.run(["zstd", "-q", "-d", path, "-o", os.path.join(cls.dbs["raw"], name)], check=True)
            else:
                os.symlink(path, os.path.join(cls.dbs["raw"], name))
                subprocess.run(["zstd", "-q", "-3", "--long=27", path, "-o", os.path.join(cls.dbs["zst"], name + ".zst")],
                               check=True)

    def test_identical_results(self):
        outputs = {}
        for kind, db in self.dbs.items():
            rc, log = run(self.work, "--db", db, *reads("sa"), "-o", f"out_{kind}", "-t", "1", "--no_qcmsa")
            self.assertEqual(rc, 0, log[-3000:])
            self.assertIn("index.prx.zst" if kind == "zst" else "index.prx", log)
            with open(glob.glob(self.path(f"out_{kind}", "sa*.sam"))[0], "rb") as sam, \
                 open(glob.glob(self.path(f"out_{kind}", "sa*.profile"))[0], "rb") as profile:
                outputs[kind] = (sam.read(), profile.read())
        self.assertEqual(outputs["raw"][0], outputs["zst"][0], "SAM differs between raw and compressed database")
        self.assertEqual(outputs["raw"][1], outputs["zst"][1], "profile differs between raw and compressed database")

    def test_no_preload_needs_a_raw_reference(self):
        rc, log = run(self.work, "--db", self.dbs["zst"], *reads("sa"), "-o", "out_lazy", "-t", "1", "--no_qcmsa",
                      "--preload_genomes_off")
        self.assertNotEqual(rc, 0)
        self.assertIn("--preload_genomes_off needs an uncompressed reference", log)


class FailFastTest(WorkDir):
    """Problems with the database or the inputs stop protal before any read is aligned."""

    def db_copy(self, name, replace=None, drop=()):
        """A database of symlinks to DB, with the files in `replace` ({name: bytes}) written instead
        and the files in `drop` left out."""
        replace = replace or {}
        db = self.path(name)
        os.mkdir(db)
        for f in glob.glob(os.path.join(DB, "*")):
            base = os.path.basename(f)
            if base not in replace and base not in drop:
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

    def test_missing_model_and_unique_kmers(self):
        db = self.db_copy("db_files", drop=("model.xml", "unique_kmers.tsv"))
        rc, log = self.query(db, "out_files")
        self.assertEqual(rc, 30, log[-3000:])
        self.assertIn("Model file does not exist", log)
        self.assertIn("Unique k-mer file does not exist", log)

    def test_corrupt_model(self):
        db = self.db_copy("db_model", {"model.xml": b"<PMML>\n"})
        rc, log = self.query(db, "out_model")
        self.assertEqual(rc, 2, log[-3000:])
        self.assertIn("Cannot load the model", log)

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
        self.assertFalse(os.path.exists(os.path.join(db, "index.prx")))
        self.assertFalse(os.path.exists(os.path.join(db, "index.prx.zst")))


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
        with open(self.path("edited.profile")) as fh:
            self.assertEqual(fh.read(), self.profile_text)

    def test_header_only_sam_gets_an_empty_profile(self):
        sam = self.write_sam("header_only", self.header)
        rc, log = self.profile_only(sam)
        self.assertEqual(rc, 0, log[-3000:])
        self.assertIn("contains no usable alignments", log)
        self.assertTrue(os.path.isfile(self.path("header_only.profile")))

    def test_unreadable_sam_fails_only_its_sample(self):
        good = self.write_sam("good", self.header + self.records)
        broken = self.write_sam("broken", self.header + self.records[:50] + ["sa.9\t0\t1_1"])
        rc, log = self.profile_only(good, broken)
        self.assertEqual(rc, 1, log[-3000:])
        self.assertRegex(log, r"Cannot read the SAM file of sample \S+ \(.*broken\.sam\): line \d+: expected at least 11")
        with open(self.path("good.profile")) as fh:
            self.assertEqual(fh.read(), self.profile_text)


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
        # qcmsa keeps a gene only if MORE than --gene-min-samples samples pass (default 3), which
        # three samples never do; relax it to exercise the output path.
        rc, log = run(self.work, "--db", DB, *reads("sa", "sb", "sr"), "-o", "out", "-t", "4",
                      "--qcmsa_script", QCMSA, "--qcmsa_args", "--gene-min-samples 1")
        self.assertEqual(rc, 0, log[-3000:])
        filtered = [f for f in glob.glob(self.path("out", "strains", "*.msa.fna")) if not f.endswith(".raw.msa.fna")]
        self.assertTrue(filtered, "qcmsa wrote filtered MSAs")


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


if __name__ == "__main__":
    sys.exit(unittest.main())
