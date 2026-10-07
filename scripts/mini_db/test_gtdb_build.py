#!/usr/bin/env python3
"""test_gtdb_build.py - checks for build_gtdb_database.py, build_gtdb_releases.py, rank_genes.py and the reduced
(gene subset) database, short of a build: the build script's options and functions, its check of the binaries' commit
(git, cmake), build_gtdb_releases.py with a stand-in build script, and the converter's gene subsets; with $PROTAL (a
protal binary) protal --build of a gene subset and rank_genes.py on a full build. A missing prerequisite skips its
test, or fails it with PROTAL_TESTS_REQUIRED=1 (scripts/prerequisites.py). The build end to end: test_gtdb_pipeline.py.

  PROTAL=build/protal python3 -m unittest scripts/mini_db/test_gtdb_build.py
"""

import collections
import contextlib
import csv
import gzip
import io
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
import unittest
import unittest.mock

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mini_db_fixtures import (CONVERT, HERE, RANK_GENES, check_reference_map, full_reference, map_rows,  # noqa: E402
                              operon_release, read_table, run, same_files, shared_release)
import prerequisites  # noqa: E402

import build_gtdb_database as build  # noqa: E402
import build_gtdb_releases as releases  # noqa: E402

class BuildOptionsTest(unittest.TestCase):
    def test_leaves_per_read_type(self):
        # build_gtdb_database.py --maxnodes: N for the read types not named, TYPE:N for one; 256 without either.
        self.assertEqual([build.max_leaves("512,pb:128,ont:128", t) for t in ("pe", "se", "pb", "ont")],
                         [512, 512, 128, 128])
        self.assertEqual(build.max_leaves("64", "ont"), 64)
        self.assertEqual(build.max_leaves("se:32", "se"), 32)
        self.assertEqual(build.max_leaves("se:32", "pe"), 256)
        for bad in ("512,xx:4", "pb:", "many"):
            with self.assertRaises(ValueError):
                build.max_leaves(bad, "pe")

    def test_versions_are_those_the_run_started_with(self):
        # build_metadata.tsv records the versions read at the run's start, with what they were at its end if a pull or
        # a rebuild changed them meanwhile (r226 v14 recorded the end's commit alone for a run that began with another).
        started = {"protal_version": "protal v0.7.8 (commit aaa)", "scripts_commit": "aaa"}
        self.assertEqual(build.versions_at_end(started, dict(started)), (started, []))
        versions, changed = build.versions_at_end(started, {"protal_version": "protal v0.7.8 (commit aaa)",
                                                            "scripts_commit": "bbb"})
        self.assertEqual(changed, ["scripts_commit"])
        self.assertEqual(versions, {"protal_version": "protal v0.7.8 (commit aaa)",
                                    "scripts_commit": "aaa; at the end of the run: bbb"})

    def test_stream_above_from_the_room(self):
        # --stream-above auto: the fewest simulations streamed, the largest samples first, for the reads of the others to
        # fit in the room at once; between two simulations' largest samples, so that equal ones go together.
        sizes = [(10, 30), (8, 16), (5, 10), (5, 5), (1, 2)]  # (largest sample, all reads), 63 together
        self.assertEqual(build.stream_threshold(sizes, 63), 0)  # all fit: none streamed
        self.assertEqual(build.stream_threshold(sizes, 40), 9)  # the one of 10 streamed: 33 left
        self.assertEqual(build.stream_threshold(sizes, 17), 6.5)  # those of 10 and 8: 17 left
        self.assertEqual(build.stream_threshold(sizes, 16), 3)  # both of 5 too: 2 left
        self.assertEqual(build.stream_threshold(sizes, 1), 0.5)  # every one
        self.assertIsNone(build.stream_threshold(sizes, -1))  # not even then
        self.assertEqual(build.stream_threshold([], 0), 0)
        self.assertEqual((build.stream_spec("auto"), build.stream_spec("AUTO"), build.stream_spec("2.5")),
                         ("auto", "auto", 2.5))
        for bad in ("-1", "lots"):
            with self.assertRaises(build.argparse.ArgumentTypeError):
                build.stream_spec(bad)

    def test_room_the_run_needs(self):
        # What the run puts on the samples' disk besides the reads: the genome store's growth (0.25 bytes a base of the
        # genomes simulated, each FASTA once, less what the store holds), a database to be built (1.5 times its
        # reference.fna; one without its files left out), the host genome of a collection that has not prepared it
        # (3.4 times a gzipped FASTA), the SAMs and profiles (a tenth of the reads) and --keep-free.
        with tempfile.TemporaryDirectory() as root:
            store = os.path.join(root, "genome_store")
            os.makedirs(store)
            with open(os.path.join(store, "a.g2b"), "wb") as fh:
                fh.write(b"x" * 100)
            table = os.path.join(root, "genomes_simulated.tsv")
            with open(table, "w") as fh:
                fh.write("a\td__B;s__x\t/g/a.fna\t1000\nb\td__B;s__y\t/g/b.fna\t3000\nb2\td__B;s__y\t/g/b.fna\t3000\n")
            training_db = os.path.join(root, "training_db")
            os.makedirs(training_db)
            with open(os.path.join(training_db, "reference.fna"), "wb") as fh:
                fh.write(b"A" * 200)
            host = os.path.join(root, "host.fna.gz")
            with open(host, "wb") as fh:
                fh.write(b"x" * 50)
            os.makedirs(os.path.join(root, "training", "host"))  # prepared by an earlier run
            needs = build.disk_needs(root, store, table, [("training database", training_db),
                                                          ("finished database", os.path.join(root, "nothing_yet"))],
                                     [(host, os.path.join(root, "training", "host")),
                                      (host, os.path.join(root, "test", "host"))], [(5, 1000), (2, 600)], 30)
            expected = {"genome store": 900, "training database": 300, "host genome": 170, "SAMs and profiles": 160,
                        "--keep-free": 30}
            self.assertEqual(set(needs), set(expected))
            for part, size in expected.items():
                self.assertAlmostEqual(needs[part], size, msg=part)
            self.assertTrue(build.on_disk(os.path.join(root, "not", "made"), root))
            self.assertNotIn("--keep-free", build.disk_needs(root, None, table, [], [], [], 0))

    @unittest.skipUnless(hasattr(os, "SCHED_IDLE") and os.path.isdir("/proc"), "Linux: SCHED_IDLE and /proc")
    def test_a_job_at_idle_priority_paused(self):
        # The finished database's build: a Job at the idle scheduling class, paused (SIGSTOP to its group) while the
        # models are trained and continued after them; one killed while paused still ends.
        with tempfile.TemporaryDirectory() as root:
            log = os.path.join(root, "job.log")
            job = build.Job([sys.executable, "-c", "import os, time; print(os.sched_getscheduler(0) == os.SCHED_IDLE, "
                                                   "flush=True); time.sleep(60)"], log, idle=True)

            def state():
                with open(f"/proc/{job.process.pid}/stat") as fh:
                    return fh.read().rsplit(")", 1)[1].split()[0]

            def said():
                with open(log) as fh:
                    return fh.read().strip()
            try:
                deadline = time.time() + 30
                while not said() and time.time() < deadline:
                    time.sleep(0.05)
                self.assertEqual(said(), "True")
                self.assertTrue(job.pause())
                deadline = time.time() + 10
                while state() != "T" and time.time() < deadline:
                    time.sleep(0.05)
                self.assertEqual(state(), "T")
                self.assertIn("(paused)", job.status())
                job.resume()
                deadline = time.time() + 10
                while state() == "T" and time.time() < deadline:
                    time.sleep(0.05)
                self.assertNotEqual(state(), "T")
                self.assertTrue(job.pause())
                began = time.time()
                job.kill()
                self.assertLess(time.time() - began, 30)  # SIGCONT after the SIGTERM: not the 60 s SIGKILL wait
                self.assertIsNotNone(job.process.poll())
            finally:
                if job.process.poll() is None:
                    job.process.kill()
                    job.process.wait()
                if job in build.Job.running:
                    build.Job.running.remove(job)

    def test_error_reads(self):
        # --error-reads: all (every read type's samples), none, READ_TYPE, READ_TYPE:design or READ_TYPE:SCENARIO.
        import scenarios
        defs = scenarios.definitions()
        self.assertEqual(build.error_read_units("all", defs), [("pe", None), ("se", None), ("pb", None), ("ont", None)])
        self.assertEqual(build.error_read_units("pe:soil,pb:design,se,se", defs), [("pe", "soil"), ("pb", "design"),
                                                                                  ("se", None)])
        self.assertEqual(build.error_read_units("none", defs), [])
        for bad in ("pe:mars", "soil", "pe:"):
            with self.assertRaises(ValueError):
                build.error_read_units(bad, defs)


class BuildFunctionsTest(unittest.TestCase):
    """build_gtdb_database.py's parts on their own: the genome table's summary and lengths, and what the build log
    says of the genes' conservation."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def test_genome_table_summary(self):
        # build_gtdb_database.py's genome table summary: how often a simulated species is a real strain or an in-silico
        # one, and its warnings (judged after the in-silico strains, when they come).
        path = os.path.join(self.tmp.name, "summary_genomes.tsv")

        def table(species):  # {name: (real strains, in-silico strains)}, each species with its representative
            reps = set()
            with open(path, "w") as fh:
                for i, (name, (real, made)) in enumerate(species.items()):
                    lineage = f"d__Bacteria;p__P;c__C;o__O;f__F;g__G;s__{name}"
                    fh.write(f"GCF_{i:09d}.1\t{lineage}\t/x.fna\t100\n")
                    reps.add(f"GCF_{i:09d}.1")
                    for j in range(real):
                        fh.write(f"GCA_{i:06d}{j:03d}.1\t{lineage}\t/x.fna\t100\n")
                    for j in range(made):  # insilico_strains.strain_name: the accession with '_' for '.'
                        fh.write(f"insilico_GCF_{i:09d}_{j + 1}\t{lineage}\t/x.fna\t100\n")
            return reps
        # Representatives only: the old warning, unless in-silico strains come next.
        reps = table({f"A{i}": (0, 0) for i in range(10)})
        lines, brief, warning = build.summarize_genome_table(path, reps)
        self.assertIn("nearly all simulated species will be the database's own reference", warning)
        lines, brief, warning = build.summarize_genome_table(path, reps, insilico_to_come=True)
        self.assertEqual(warning, "")
        self.assertIn("  the species with one genome get in-silico strains next (step 3)", lines)
        # Two species of two real strains, eight one-genome species with an in-silico strain: 13.3% real, 40% in silico.
        reps = table({**{f"R{i}": (2, 0) for i in range(2)}, **{f"M{i}": (0, 1) for i in range(8)}})
        lines, brief, warning = build.summarize_genome_table(path, reps)
        self.assertIn("another genome than its representative 53.3% of the time (a real strain 13.3%, an in-silico strain "
                      "40.0%)", brief)
        self.assertIn("WARNING: most simulated strains are in-silico (40.0% of the simulated species against 13.3% real)",
                      warning)
        # Mostly real strains: no warning.
        reps = table({**{f"R{i}": (2, 0) for i in range(8)}, **{f"M{i}": (0, 1) for i in range(2)}})
        self.assertEqual(build.summarize_genome_table(path, reps)[2], "")

    def test_gene_conservation_in_the_build_metadata(self):
        # build_metadata.tsv records what protal --build said of the genes' conservation factors.
        log = os.path.join(self.tmp.name, "index_and_package.log")
        with open(log, "w") as fh:
            fh.write("Uniqueness check took 1s\nGene conservation: factors 0.26-3.3 for 168 genes, from 120 species "
                     "(65559 copies compared): /data/db/gene_conservation.tsv\nGene conservation took 1s\n")
        self.assertEqual(build.gene_conservation_summary(log),
                         "factors 0.26-3.3 for 168 genes, from 120 species (65559 copies compared)")
        with open(log, "w") as fh:
            fh.write("Gene conservation: no factors (0 species with other genomes' copies of their genes in x.fna, 0 of "
                     "them with enough genes that differ from the representative's): every gene keeps the whole margin\n")
        self.assertTrue(build.gene_conservation_summary(log).startswith("no factors (0 species"))
        with open(log, "w") as fh:
            fh.write("Run build took 5s\n")
        self.assertEqual(build.gene_conservation_summary(log), "none (this protal does not estimate them)")
        self.assertEqual(build.gene_conservation_summary(log + ".missing"), "unknown (no build log)")

    def test_genome_table_lengths(self):
        # build_gtdb_database.py gives the simulator each genome's length (it would read every genome for it at
        # each design point): letters outside header lines, gzipped or not; a given table without lengths gets a
        # copy with them, one with them is taken as it is.
        root = os.path.join(self.tmp.name, "lengths")
        os.makedirs(root)
        plain, zipped = os.path.join(root, "a.fna"), os.path.join(root, "b.fna.gz")
        text = ">a one 123 ACGT\r\nACGTNNacgt\r\n\r\nRYK-*\n>second\nAC GT\n"  # 10 + 3 + 4 letters
        with open(plain, "w", newline="") as fh:
            fh.write(text)
        with gzip.open(zipped, "wt", newline="") as fh:
            fh.write(text * 2)
        self.assertEqual(build.genome_length(plain), 17)
        self.assertEqual(build.genome_length(zipped), 34)
        given = os.path.join(root, "given.tsv")
        with open(given, "w") as fh:
            fh.write(f"name\ttaxonomy\tfasta_path\nGA\td__B;s__A\t{plain}\nGB\td__B;s__B\t{zipped}\n")
        copy = build.with_lengths(given, os.path.join(root, "genomes.tsv"), 2)
        with open(copy) as fh:
            self.assertEqual(fh.read(), f"GA\td__B;s__A\t{plain}\t17\nGB\td__B;s__B\t{zipped}\t34\n")
        self.assertEqual(build.with_lengths(copy, os.path.join(root, "other.tsv"), 2), copy)


class CladeHoldoutTest(unittest.TestCase):
    """The species build_gtdb_database.py leaves out of the training database (whole clades, then single species),
    and the lineages the training scripts read, on the default synthetic release (shared_release)."""

    def test_lineages_and_clade_holdout(self):
        # The training scripts read lineages from internal_taxonomy.dmp; build_gtdb_database.py holds out whole
        # clades and then single species.
        import collect_training_data as collect
        import lineages
        gtdb, db = shared_release()
        taxonomy = os.path.join(db, "internal_taxonomy.dmp")
        by_id, ids = lineages.from_taxonomy(taxonomy)
        with open(os.path.join(db, "genome2tiid.tsv")) as fh:
            for accession, taxid, _rep, lineage in (line.rstrip("\n").split("\t") for line in fh):
                self.assertEqual(by_id[taxid], lineages.from_string(lineage), accession)
        table = os.path.join(gtdb, "simulation", "genomes.tsv")
        pool = build.pool_species(table)
        species = {lin["species"]: lin for lin in by_id.values() if "species" in lin}
        families = {}
        for name, lin in species.items():
            families.setdefault(lin["family"], set()).add(name)
        eligible = {f for f, members in families.items() if len(members & pool) >= 2}
        self.assertTrue(eligible, "the synthetic release needs a family with two species")
        chosen = build.choose_holdout(table, taxonomy, 0.5, {"family": 1}, 1.0, 3)
        self.assertEqual(chosen, build.choose_holdout(table, taxonomy, 0.5, {"family": 1}, 1.0, 3))
        held_family = {clade for rank, clade in chosen.values() if rank == "family"}
        self.assertEqual(len(held_family), 1)
        family = held_family.pop()
        self.assertIn(family, eligible)
        self.assertEqual({s for s, (rank, _) in chosen.items() if rank == "family"}, families[family])
        # Then that share of the species left, of each domain, rounded: of the one left (Fakibacter gamma), none at
        # 0.5, itself at 1.
        self.assertEqual(chosen, {"s__Mockella alpha": ("family", "f__Simulaceae"),
                                  "s__Mockella beta": ("family", "f__Simulaceae")})
        self.assertEqual(build.choose_holdout(table, taxonomy, 1.0, {"family": 1}, 1.0, 3),
                         {**chosen, "s__Fakibacter gamma": ("species", "s__Fakibacter gamma")})
        # a clade holding more than the share, or inside one taken before, is not drawn
        self.assertEqual(build.choose_holdout(table, taxonomy, 0, {"family": 1}, 0.0, 3), {})
        phylum_first = build.choose_holdout(table, taxonomy, 0, {"phylum": 1, "family": 5}, 1.0, 3)
        phylum = next(c for r, c in phylum_first.values() if r == "phylum")
        for s, (rank, clade) in phylum_first.items():
            if rank == "family":
                self.assertNotEqual(species[s]["phylum"], phylum)
        # each rank alone takes a whole clade of that rank; the synthetic release has one clade with two
        # species, from phylum to genus, so with every rank asked for, the first (phylum) takes it
        for rank in build.CLADE_RANKS:
            alone = build.choose_holdout(table, taxonomy, 0, {rank: 1}, 1.0, 3)
            taken = {c for r, c in alone.values()}
            self.assertEqual({r for r, _ in alone.values()}, {rank})
            self.assertEqual(len(taken), 1)
            clade = taken.pop()
            self.assertEqual(set(alone), {s for s, lin in species.items() if lin.get(rank) == clade})
            self.assertEqual(collect.novel_clades(alone, table, 1), {rank: [clade]})
        every = build.choose_holdout(table, taxonomy, 0, build.parse_clades("phylum:1,class:1,order:1,family:1,genus:2"), 1.0, 3)
        self.assertEqual({r for r, _ in every.values()}, {"phylum"})
        # heldout_species.txt round trip, and the collector's reading of it
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "heldout.txt")
            with open(path, "w") as fh:
                fh.write("".join(f"{s}\t{r}\t{c}\n" for s, (r, c) in sorted(chosen.items())) + "s__Alone one\n")
            self.assertEqual(build.read_holdout(path), {**chosen, "s__Alone one": ("species", "s__Alone one")})
            self.assertEqual(collect.read_novel(path), build.read_holdout(path))
        self.assertEqual(build.parse_clades("phylum:2,class:4"), {"phylum": 2, "class": 4})
        self.assertEqual(build.parse_clades("none"), {})


class BinaryCheckTest(unittest.TestCase):
    """The commit a build records for --version (protal_commit.cmake), and build_gtdb_database.py's check at its start
    that protal and the simulator were built from the source its scripts are at, on a throwaway checkout."""

    def setUp(self):
        if not prerequisites.on_path("git"):
            prerequisites.missing("no git")
        self.tmp = tempfile.TemporaryDirectory()
        self.source = os.path.join(self.tmp.name, "checkout")
        os.makedirs(os.path.join(self.source, "src"))
        self.write("CMakeLists.txt", "project(protal VERSION 0.7.3)\n")
        self.write("src/a.h", "int a;\n")
        self.git("init", "-q")
        self.commit("one")
        self.built = self.git("rev-parse", "HEAD")

    def tearDown(self):
        self.tmp.cleanup()

    def write(self, path, text, mode="w"):
        with open(os.path.join(self.source, path), mode) as fh:
            fh.write(text)

    def git(self, *arguments):
        return subprocess.run(["git", "-C", self.source, "-c", "user.name=test", "-c", "user.email=test@example.org",
                               *arguments], check=True, capture_output=True, text=True).stdout.strip()

    def commit(self, message):
        self.git("add", "-A")
        self.git("commit", "-q", "-m", message)

    def check(self, said, code=0):
        # A stand-in binary that says `said` to --version.
        binary = os.path.join(self.tmp.name, "protal")
        with open(binary, "w") as fh:
            fh.write(f"#!/bin/sh\necho '{said}'\nexit {code}\n")
        os.chmod(binary, 0o755)
        return build.build_check(binary, "protal", self.source)

    def test_the_commit_in_the_build(self):
        if not prerequisites.on_path("cmake"):
            prerequisites.missing("no cmake")
        header = os.path.join(self.tmp.name, "protal_commit.h")

        def written():
            subprocess.run(["cmake", f"-DSOURCE_DIR={self.source}", f"-DOUTPUT={header}", "-P",
                            os.path.join(HERE, "..", "..", "protal_commit.cmake")], check=True, capture_output=True)
            with open(header) as fh:
                return fh.read()

        def expected(commit, modified):
            return f'#pragma once\n#define PROTAL_GIT_COMMIT "{commit}"\n#define PROTAL_GIT_MODIFIED {modified}\n'

        self.assertEqual(written(), expected(self.built, 0))
        self.write("README.md", "docs\n")  # not what the binaries are built from
        self.assertEqual(written(), expected(self.built, 0))
        self.write("src/b.h", "int b;\n")  # a new file in src/ is
        self.assertEqual(written(), expected(self.built, 1))
        shutil.rmtree(os.path.join(self.source, ".git"))
        self.assertEqual(written(), expected("", 0))

    def test_binaries_of_another_source_stop_the_run(self):
        self.assertEqual(self.check(f"protal v0.7.3 (commit {self.built})"), (None, None))
        # Commits that leave src/, lib/ and the build files alone keep the binary current.
        self.write("README.md", "docs\n")
        self.commit("docs")
        self.assertEqual(self.check(f"protal v0.7.3 (commit {self.built})"), (None, None))
        problem, _ = self.check("protal v0.7.1")
        self.assertIn("is v0.7.1, these scripts v0.7.3", problem)
        problem, _ = self.check("unknown option(s): --version", 1)  # simulate_metagenomes before 0.7.3
        self.assertIn("--version gives no version (1: unknown option", problem)
        problem, _ = self.check(f"protal v0.7.3 (commit {'0' * 40})")
        self.assertIn(f"was built from commit 0000000000, which {self.source} does not have", problem)
        self.write("src/a.h", "int a2;\n")
        problem, _ = self.check(f"protal v0.7.3 (commit {self.built})")
        self.assertIn(f"was built from commit {self.built[:10]}, and src/, lib/ or the build files have changed since "
                      "(uncommitted changes): rebuild it", problem)
        self.commit("two")
        problem, _ = self.check(f"protal v0.7.3 (commit {self.built})")
        self.assertIn("have changed since (1 commit): rebuild it", problem)
        self.assertIsNone(self.check(f"protal v0.7.3 (commit {self.git('rev-parse', 'HEAD')})")[0])

    def test_binaries_that_cannot_say_are_noted(self):
        problem, note = self.check("protal v0.7.3")
        self.assertIsNone(problem)
        self.assertIn("does not say which commit it was built from", note)
        problem, note = self.check(f"protal v0.7.3 (commit {self.built}, with uncommitted changes)")
        self.assertIsNone(problem)
        self.assertIn("was built with uncommitted changes to its source", note)
        shutil.rmtree(os.path.join(self.source, ".git"))
        problem, note = self.check(f"protal v0.7.3 (commit {self.built})")
        self.assertIsNone(problem)
        self.assertIn("is no git checkout", note)


class GeneSubsetTest(unittest.TestCase):
    """A database of a subset of the marker genes, on the operon release (operon_release). The converter's --genes
    (with --from_db, and with --gtdb the same files) keeps the genes named, by GTDB marker id or protal gene id, with
    their ids in reference.fna, reference.map, the full reference and gene2geneid.tsv, and counts the gene neighbours
    anew over them: a gene whose neighbour in the genome is left out faces the nearest gene kept, as a read of the
    reduced database meets it. With $PROTAL: protal --build of such a folder lists only its genes in unique_kmers.tsv
    and gene_conservation.tsv; rank_genes.py ranks a full build's genes and writes a gene list that the converter and
    --build_gene_subset take; and --build_gene_subset refuses a folder whose neighbours were counted over every
    gene."""

    MAX_GAP = 3000

    @classmethod
    def setUpClass(cls):
        cls.gtdb, cls.db, _ = operon_release()
        cls.tmp = tempfile.TemporaryDirectory()
        from gtdb_to_protal_db import read_gene_ids
        cls.gene_ids = read_gene_ids(os.path.join(cls.db, "gene2geneid.tsv"))
        cls.markers = {gid: marker for marker, gid in cls.gene_ids.items()}
        # The subset: every other gene along the longest contig of one representative (the genes in between are
        # left out, so their neighbours in the subset's tables are the next genes kept): the first two by their
        # marker id, the third by the id without its version, the rest by protal gene id.
        with open(os.path.join(cls.db, "internal_taxonomy.dmp")) as fh:
            next(fh)
            cls.nodes = {int(r[0]): r for r in (l.rstrip("\n").split("\t") for l in fh)}
        cls.taxid, cls.rep = next((t, r[6]) for t, r in cls.nodes.items() if r[3] == "s__Mockella s0")
        by_contig = collections.defaultdict(list)
        for r in read_table(os.path.join(cls.db, "gene_positions.tsv")):
            if r["accession"] == cls.rep:
                by_contig[r["contig"]].append((int(r["start"]), int(r["end"]), int(r["gene"])))
        along = sorted(max(by_contig.values(), key=len))
        # Three genes in a row within the neighbour gap: the first and the third are kept, the middle one left
        # out, so that in the subset the first faces the third; four more genes from elsewhere on the contig.
        first = next(i for i in range(len(along) - 2)
                     if along[i + 1][0] - along[i][1] <= cls.MAX_GAP and along[i + 2][0] - along[i + 1][1] <= cls.MAX_GAP
                     and along[i + 2][0] - along[i][1] <= cls.MAX_GAP)
        others = [k for k in range(len(along)) if k not in (first, first + 1, first + 2)]
        indices = [first, first + 2] + others[::max(1, len(others) // 4)][:4]
        cls.subset = sorted(along[k][2] for k in indices)
        tokens = [cls.markers[g] for g in cls.subset[:2]] + [cls.markers[cls.subset[2]].split(".")[0]] + \
            [str(g) for g in cls.subset[3:]]
        cls.spec = ",".join(tokens)
        cls.copy = os.path.join(cls.tmp.name, "subset")
        cls.derived = subprocess.run([sys.executable, CONVERT, "--from_db", cls.db, "--genes", cls.spec, "--outdir", cls.copy],
                                     check=True, capture_output=True, text=True).stderr

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    @staticmethod
    def records(data):
        """[(gene id, header, sequence)] of a FASTA's bytes, one sequence line per record."""
        lines = data.split(b"\n")
        return [(int(lines[i][1:].split(b"_")[1]), lines[i], lines[i + 1]) for i in range(0, len(lines) - 1, 2)]

    def test_the_copy_holds_the_genes_named_with_their_ids(self):
        self.assertEqual(len(self.subset), 6)
        rows = check_reference_map(self, self.copy)
        self.assertEqual(sorted({r[1] for r in rows}), self.subset)
        self.assertEqual([r[:2] for r in rows], [r[:2] for r in map_rows(self.db) if r[1] in self.subset])
        with open(os.path.join(self.copy, "reference.fna"), "rb") as fh:
            data = fh.read()
        with open(os.path.join(self.db, "reference.fna"), "rb") as fh:
            self.assertEqual(self.records(data), [r for r in self.records(fh.read()) if r[0] in self.subset])
        self.assertEqual(self.records(full_reference(self.copy)),
                         [r for r in self.records(full_reference(self.db)) if r[0] in self.subset])
        with open(os.path.join(self.copy, "gene2geneid.tsv")) as fh:
            kept = {m: int(g) for m, g in (line.rstrip("\n").split("\t") for line in fh)}
        self.assertEqual(sorted(kept.values()), self.subset)
        self.assertEqual(kept, {m: g for m, g in self.gene_ids.items() if g in self.subset})
        same_files(self, self.db, self.copy, "internal_taxonomy.dmp", "species_priors.tsv", "genome2tiid.tsv")
        # The genomes with a gene kept: the archaea have none of the bacterial genes chosen here.
        self.assertRegex(self.derived, r"gene_neighbours\.tsv derived from the 1[2-6] genomes of [6-8] species kept: \d+ lines "
                                       r"\(over the 6 genes kept\)")
        self.assertRegex(self.derived, r"without 0 species, 6 genes kept: \d+ representative sequences kept")

    def test_a_direct_conversion_gives_the_same_files(self):
        direct = os.path.join(self.tmp.name, "direct")
        run(CONVERT, "--gtdb", self.gtdb, "--outdir", direct, "--genes", self.spec)
        same_files(self, direct, self.copy, "reference.fna", "reference.map", "gene2geneid.tsv", "internal_taxonomy.dmp",
                   "species_priors.tsv")
        self.assertEqual(full_reference(direct), full_reference(self.copy))
        # A file lists the genes too (its first column), with or without the species left out.
        listed = os.path.join(self.tmp.name, "genes.txt")
        with open(listed, "w") as fh:
            fh.write("# the subset\n" + "".join(f"{t}\tcomment\n" for t in self.spec.split(",")))
        excluded = os.path.join(self.tmp.name, "excluded.txt")
        with open(excluded, "w") as fh:
            fh.write("s__Otherella s1\n")
        both = os.path.join(self.tmp.name, "both")
        run(CONVERT, "--from_db", self.db, "--genes", listed, "--exclude_species", excluded, "--outdir", both)
        rows = map_rows(both)
        self.assertEqual(sorted({r[1] for r in rows}), self.subset)
        self.assertEqual(len({r[0] for r in rows}), len({r[0] for r in map_rows(self.copy)}) - 1)
        same_files(self, both, self.copy, "gene2geneid.tsv")
        for bad in ("PF99999.1", "0", str(max(self.gene_ids.values()) + 1), "#"):
            with self.assertRaises(subprocess.CalledProcessError, msg=bad):
                run(CONVERT, "--from_db", self.db, "--genes", bad, "--outdir", os.path.join(self.tmp.name, "bad"))

    def test_neighbours_are_counted_over_the_genes_kept(self):
        subset = set(self.subset)
        kept = read_table(os.path.join(self.copy, "gene_neighbours.tsv"))
        self.assertTrue(kept)
        self.assertTrue(all(int(r["gene"]) in subset and int(r["partner"]) in subset | {0} for r in kept))
        positions = read_table(os.path.join(self.copy, "gene_positions.tsv"))
        self.assertEqual({int(r["gene"]) for r in positions}, subset)
        self.assertEqual({r["accession"] for r in positions},
                         {r["accession"] for r in read_table(os.path.join(self.db, "gene_positions.tsv")) if int(r["gene"]) in subset})
        # What the representative's genome gives among the genes kept (neighbour_ends on its placements) is in
        # the copy's table, under the species or one of its clades; a gene whose neighbour was left out faces the
        # next gene kept, farther away, where the full table had the gene left out.
        import gene_neighbours as gn
        ancestors, node = set(), self.taxid
        while node not in ancestors:
            ancestors.add(node)
            node = int(self.nodes[node][1])
        species, settings = gn.read_positions(os.path.join(self.db, "gene_positions.tsv"), genes=subset)
        max_gap = int(settings.get("max_gap", self.MAX_GAP))
        expected = gn.neighbour_ends(species[self.taxid][self.rep], max_gap)
        lines = {(int(r["clade"]), int(r["gene"]), int(r["end"]), int(r["partner"]), int(r["partner_end"])) for r in kept}
        for gene, end, partner, partner_end, gap in expected:
            self.assertTrue(any((clade, gene, end, partner, partner_end) in lines for clade in ancestors),
                            f"gene {gene} end {end} facing {partner} (end {partner_end}, {gap} bases)")
        whole, _ = gn.read_positions(os.path.join(self.db, "gene_positions.tsv"))
        before = {(g, e): (p, gap) for g, e, p, _, gap in gn.neighbour_ends(whole[self.taxid][self.rep], max_gap)}
        bridged = [(g, e, p, gap) for g, e, p, _, gap in expected if p and before[(g, e)][0] not in subset]
        self.assertTrue(bridged, "no gene end faces a gene kept where a gene left out was")
        self.assertTrue(all(gap > before[(g, e)][1] for g, e, p, gap in bridged))

    def test_protal_builds_the_subset_and_ranks_a_full_build(self):
        protal = os.environ.get("PROTAL", "")
        if not prerequisites.executable(protal):
            prerequisites.missing("$PROTAL names no protal binary")
        from gtdb_to_protal_db import full_reference_path

        def build_db(folder, *extra):  # zstd level 1 without long-distance matching: the database's content is the same
            return subprocess.run([protal, "--build", "--no_profile", "-t", "2", "--no_bundle", "--compress_level", "1",
                                   "--compress_window_log", "0", "--db", folder, "--reference",
                                   os.path.join(folder, "reference.fna"), "--full_reference", full_reference_path(folder),
                                   *extra], capture_output=True, text=True)

        def genes_of(folder, name):
            with open(os.path.join(folder, name)) as fh:
                return sorted({int(line.split("\t")[1]) for line in fh if line[0].isdigit() and line.split("\t")[1].isdigit()})

        subset = os.path.join(self.tmp.name, "subset_built")
        shutil.copytree(self.copy, subset)
        result = build_db(subset)
        self.assertEqual(result.returncode, 0, result.stdout[-2000:] + result.stderr[-2000:])
        self.assertEqual(genes_of(subset, "unique_kmers.tsv"), self.subset)
        if os.path.isfile(os.path.join(subset, "gene_conservation.tsv")):  # too few genes per species for factors here
            self.assertTrue(set(genes_of(subset, "gene_conservation.tsv")) <= set(self.subset))
        else:
            self.assertIn("Gene conservation: no factors", result.stdout)
        self.assertRegex(result.stdout, r"Gene neighbours: \d+ rules of \d+ clades from 1[2-6] genomes")
        # Every gene, built and ranked: a table of all genes with the domains' columns, the 3 chosen as a gene
        # list, one of them reserved for archaea (a gene of their marker set alone is rare over all species).
        every = os.path.join(self.tmp.name, "every_gene")
        run(CONVERT, "--from_db", self.db, "--outdir", every)
        same_files(self, every, self.db, "reference.map")
        result = build_db(every)
        self.assertEqual(result.returncode, 0, result.stdout[-2000:] + result.stderr[-2000:])
        table, listed = os.path.join(self.tmp.name, "ranking.tsv"), os.path.join(self.tmp.name, "best.txt")
        ranking = subprocess.run([sys.executable, RANK_GENES, "--db", every, "-o", table, "--top", "3", "--subset", listed],
                                 check=True, capture_output=True, text=True)
        self.assertRegex(ranking.stderr, r"best\.txt: 3 genes: \S+, \S+, \S+; scores [0-9.]+ down to [0-9.]+; in half the "
                                         r"species or more of: bacteria \d, archaea \d")
        rows = read_table(table)
        self.assertEqual(len(rows), len({r[1] for r in map_rows(every)}))  # the reference's genes
        self.assertEqual([int(r["rank"]) for r in rows], list(range(1, len(rows) + 1)))
        scores = [float(r["score"]) for r in rows]
        self.assertEqual(scores, sorted(scores, reverse=True))
        self.assertTrue(all(0 <= float(r[c]) <= 1 for r in rows
                            for c in ("score", "prevalence", "unique_share", "bacteria_prevalence", "archaea_prevalence")))
        self.assertTrue(all(self.gene_ids[r["marker"]] == int(r["gene_id"]) for r in rows))
        self.assertTrue(all(int(r["species"]) <= 8 and float(r["mean_length"]) > 0 for r in rows))
        self.assertGreater(scores[0], 0)
        archaeal = {int(r["gene_id"]) for r in rows if float(r["archaea_prevalence"]) >= 0.5}
        bacterial = {int(r["gene_id"]) for r in rows if float(r["bacteria_prevalence"]) >= 0.5}
        self.assertTrue(archaeal - bacterial, "genes of the archaeal marker set alone")
        self.assertTrue(all(float(r["score"]) < 0.5 for r in rows if int(r["gene_id"]) in archaeal - bacterial))
        with open(listed) as fh:
            lines = [line.rstrip("\n") for line in fh]
        best = [int(line) for line in lines if not line.startswith("#")]
        self.assertEqual(len(best), 3)
        self.assertTrue(set(best) & archaeal, "a gene archaea have")
        self.assertTrue(set(best) & bacterial, "a gene bacteria have")
        self.assertEqual(best, sorted(best, key=lambda g: next(int(r["rank"]) for r in rows if int(r["gene_id"]) == g)))
        self.assertEqual(sum(line.startswith("# rank ") for line in lines), 3)
        self.assertTrue(any("chosen for archaea" in line for line in lines))
        # Without the domains' share, the 3 best by the overall score are bacterial genes only.
        plain = os.path.join(self.tmp.name, "plain.txt")
        subprocess.run([sys.executable, RANK_GENES, "--db", every, "--top", "3", "--per-domain", "0", "--subset", plain],
                       check=True, capture_output=True, text=True)
        with open(plain) as fh:
            self.assertEqual([int(line) for line in fh if not line.startswith("#")], [int(r["gene_id"]) for r in rows[:3]])
        self.assertFalse({int(r["gene_id"]) for r in rows[:3]} & (archaeal - bacterial))
        from_list = os.path.join(self.tmp.name, "from_list")
        run(CONVERT, "--from_db", self.db, "--genes", listed, "--outdir", from_list)
        self.assertEqual(sorted({r[1] for r in map_rows(from_list)}), sorted(best))
        # protal --build_gene_subset: refused with a neighbours table counted over every gene; without one, the
        # database gets the listed genes' k-mers, rows and factors only, the other genes stay in reference.fna.
        with_neighbours = os.path.join(self.tmp.name, "with_neighbours")
        run(CONVERT, "--from_db", self.db, "--outdir", with_neighbours)
        result = build_db(with_neighbours, "--build_gene_subset", listed)
        self.assertEqual(result.returncode, 8)
        self.assertIn("Cannot build with --build_gene_subset: the folder has gene_neighbours.tsv", result.stderr)
        for name in ("gene_neighbours.tsv", "gene_positions.tsv"):
            os.remove(os.path.join(with_neighbours, name))
        result = build_db(with_neighbours, "--build_gene_subset", listed)
        self.assertEqual(result.returncode, 0, result.stdout[-2000:] + result.stderr[-2000:])
        self.assertEqual(genes_of(with_neighbours, "unique_kmers.tsv"), sorted(best))
        if os.path.isfile(os.path.join(with_neighbours, "gene_conservation.tsv")):
            self.assertTrue(set(genes_of(with_neighbours, "gene_conservation.tsv")) <= set(best))
        self.assertEqual(map_rows(with_neighbours), map_rows(self.db))  # every gene is still in the reference


# A stand-in for build_gtdb_database.py: each call is a line of JSON (its arguments) in $FAKE_BUILD_LOG; it fails for an
# OUTDIR named in $FAKE_BUILD_FAIL; else it writes what build_gtdb_releases.py reads of a build (the database, its
# metadata, the taxonomy, the models' summary, the ranking with --rank-genes) and the console lines of a build's times.
FAKE_BUILD = r'''#!/usr/bin/env python3
import json, os, sys
args = sys.argv[1:]
get = lambda k: args[args.index(k) + 1] if k in args else None
out = get("--outdir")
with open(os.environ["FAKE_BUILD_LOG"], "a") as fh:
    fh.write(json.dumps(args) + "\n")
if os.path.basename(out) in os.environ.get("FAKE_BUILD_FAIL", "").split(","):
    sys.exit("the stand-in fails as asked")
db = os.path.join(out, "protal_db")
os.makedirs(db, exist_ok=True)
os.makedirs(os.path.join(out, "model_logs"), exist_ok=True)
with open(os.path.join(db, "database.protal"), "wb") as fh:
    fh.write(b"\0" * 12000000)
n = get("--n-genes")
genes = "all" if n is None else f"{n} of 120, the most distinctive (1 per domain, from --gene-ranking): A, B, C"
with open(os.path.join(db, "build_metadata.tsv"), "w") as fh:
    fh.write(f"protal_version\tprotal v0.7.8 (commit 1234)\nmarker_genes\t{genes}\ngenome_table\t30 genomes of 12 species\n"
             "classifier_training_species_left_out\t4\n")
nodes = [(1, 1, "root", "no rank"), (2, 1, "d__Bacteria", "domain"), (3, 1, "d__Archaea", "domain"),
         (4, 2, "p__A", "phylum"), (5, 3, "p__B", "phylum"), (6, 4, "c__A", "class"), (7, 5, "c__B", "class"),
         (8, 6, "o__A", "order"), (9, 7, "o__B", "order"), (10, 8, "f__A", "family"), (11, 9, "f__B", "family"),
         (12, 10, "g__A", "genus"), (13, 11, "g__B", "genus"), (14, 12, "s__A a", "species"),
         (15, 12, "s__A b", "species"), (16, 13, "s__B a", "species")]
with open(os.path.join(out, "internal_taxonomy.dmp"), "w") as fh:
    fh.write("taxid\tparent\tlevel\tname\trank\n" + "".join(f"{t}\t{p}\t-\t{name}\t{rank}\n" for t, p, name, rank in nodes))
with open(os.path.join(out, "model_logs", "summary.txt"), "w") as fh:
    fh.write("read type  evaluated on           knob  taxa  TP  FP  TN  FN  sensitivity  specificity  precision  F1      FP rate  FN rate  FP/sample\n"
             "pe         species held out       0.5   100   9   1   88  2   0.8182       0.9888       0.9000     0.8571  0.0112   0.1818   0.50\n"
             "pe         independent test set   0.5   50    5   0   44  1   0.8333       1.0000       1.0000     0.9091  0.0000   0.1667   0.00\n")
if "--rank-genes" in args:
    with open(os.path.join(out, "gene_ranking.tsv"), "w") as fh:
        fh.write("rank\tgene_id\tmarker\tscore\n1\t7\tA\t0.9\n")
print("[00:00:01]     built protal_db in the background in 0:01:02, peak memory 1.5 GB; full_reference.fna removed")
print("[00:00:02]     collected in 0:00:30, peak memory 512 MB; taxa present/absent: pe 9/91")
'''


class BuildReleasesTest(unittest.TestCase):
    """build_gtdb_releases.py with a stand-in build script: the full database of each release first (with --rank-genes),
    then the reduced ones from its ranking, each in OUT/r<release>_<variant> with its options and those it does not
    know passed on; the summary of every database; a database built before is kept; a failure is reported, the
    others built."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.inputs = os.path.join(self.tmp.name, "inputs")
        os.makedirs(os.path.join(self.inputs, "gtdb_r226"))
        with open(os.path.join(self.inputs, "gtdb_r226", "download.json"), "w") as fh:
            json.dump({"release": {"number": "226"}}, fh)
        os.makedirs(os.path.join(self.inputs, "gtdb_r220"))  # no download.json: not a release
        self.out = os.path.join(self.tmp.name, "dbs")
        self.calls = os.path.join(self.tmp.name, "calls.jsonl")
        fake = os.path.join(self.tmp.name, "build_gtdb_database.py")
        with open(fake, "w") as fh:
            fh.write(FAKE_BUILD)
        self.addCleanup(setattr, releases, "BUILD", releases.BUILD)
        releases.BUILD = fake

    def releases(self, *extra, fail=""):
        """build_gtdb_releases.main() with the stand-in -> (exit code, console)."""
        argv = ["build_gtdb_releases.py", "--inputs", self.inputs, "--outdir", self.out, "--variants", "n3,full",
                "--n-genes", "3", "--protal", "/p/protal", "--simulator", "/p/sim", "-t", "2",
                "--scratch", os.path.join(self.tmp.name, "scratch"), "--samples", "2", *extra]
        console = io.StringIO()
        code = 0
        environment = {"FAKE_BUILD_LOG": self.calls, "FAKE_BUILD_FAIL": fail}
        saved = {k: os.environ.get(k) for k in environment}
        os.environ.update(environment)
        try:
            with contextlib.redirect_stdout(console), unittest.mock.patch.object(sys, "argv", argv):
                releases.main()
        except SystemExit as stop:
            code = stop.code
        finally:
            for k, v in saved.items():
                if v is None:
                    os.environ.pop(k, None)
                else:
                    os.environ[k] = v
        return code, console.getvalue()

    def builds(self):
        """The stand-in's calls: [{option: value}], in order."""
        if not os.path.exists(self.calls):
            return []
        with open(self.calls) as fh:
            calls = [json.loads(line) for line in fh]
        return [{a: (c[i + 1] if i + 1 < len(c) and not c[i + 1].startswith("--") else True)
                 for i, a in enumerate(c) if a.startswith("-")} for c in calls]

    def summary(self):
        with open(os.path.join(self.out, "build_summary.tsv")) as fh:
            return list(csv.DictReader(fh, delimiter="\t"))

    def test_full_first_then_the_reduced_from_its_ranking(self):
        code, console = self.releases()
        self.assertIn(code, (0, None), console)
        full, reduced = self.builds()
        self.assertEqual((full["--outdir"], reduced["--outdir"]),
                         (os.path.join(self.out, "r226_full"), os.path.join(self.out, "r226_n3")))
        self.assertEqual(full["--inputs"], os.path.join(self.inputs, "gtdb_r226"))
        self.assertTrue(full.get("--rank-genes"))
        self.assertNotIn("--n-genes", full)
        self.assertEqual((reduced["--n-genes"], reduced["--gene-ranking"]),
                         ("3", os.path.join(self.out, "r226_full", "gene_ranking.tsv")))
        self.assertNotIn("--rank-genes", reduced)
        for call, variant in ((full, "full"), (reduced, "n3")):
            self.assertEqual(call["--scratch"], os.path.join(self.tmp.name, "scratch", f"r226_{variant}"))
            self.assertEqual((call["--samples"], call["--protal"], call["-t"]), ("2", "/p/protal", "2"))  # passed on
        self.assertTrue(os.path.isfile(os.path.join(self.out, "r226_full.log")))
        rows = self.summary()
        self.assertEqual([(r["release"], r["variant"], r["status"]) for r in rows], [("226", "full", "ok"), ("226", "n3", "ok")])
        self.assertEqual([r["genes"] for r in rows], ["all", "3"])
        for r in rows:
            self.assertEqual((r["protal_version"], r["species"], r["genera"], r["phyla"], r["bacteria_species"],
                              r["archaea_species"]), ("protal v0.7.8 (commit 1234)", "3", "2", "2", "2", "1"))
            self.assertEqual((r["genomes_simulated"], r["species_held_out"], r["database_gb"]),
                             ("30 genomes of 12 species", "4", "0.01"))
            self.assertEqual((r["build_time"], r["build_peak_gb"], r["profiling_peak_gb"]), ("0:01:02", "1.5", "0.5"))
            self.assertEqual((r["pe_test_F1"], r["pe_test_FP_per_sample"], r["pe_heldout_F1"], r["se_test_F1"]),
                             ("0.9091", "0.00", "0.8571", ""))
        with open(os.path.join(self.out, "build_summary.txt")) as fh:
            text = fh.read()
        self.assertRegex(text, r"GTDB r226, full: ok; (protal )?protal v0\.7\.8 \(commit 1234\)")
        self.assertIn("  pe: independent test F1 0.9091, FP/sample 0.00", text)
        self.assertIn("Summary: " + os.path.join(self.out, "build_summary.tsv"), console)
        # A rerun keeps both; --rerun builds them again.
        code, console = self.releases()
        self.assertIn(code, (0, None), console)
        self.assertEqual(len(self.builds()), 2)
        self.assertEqual([r["status"] for r in self.summary()], ["kept", "kept"])
        self.releases("--rerun")
        self.assertEqual(len(self.builds()), 4)

    def test_a_failed_database_is_reported_and_the_others_built(self):
        code, console = self.releases(fail="r226_full")
        self.assertEqual(code, 1)
        self.assertEqual(len(self.builds()), 2)
        self.assertNotIn("--gene-ranking", self.builds()[1])  # no ranking from the failed full database
        self.assertEqual([r["status"] for r in self.summary()], ["failed (1)", "ok"])
        self.assertIn("failed: r226 full", console)
        # A release without download.json is skipped, and the run fails.
        code, console = self.releases("--releases", "220,226")
        self.assertEqual(code, 1)
        self.assertIn("failed: r220", console)


if __name__ == "__main__":
    unittest.main()
