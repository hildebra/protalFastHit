#!/usr/bin/env python3
"""test_gtdb_pipeline.py - build_gtdb_database.py end to end (GtdbBuildTest), on a synthetic GTDB-like release served
by stand-ins of GTDB's mirror and NCBI: the protal databases built, the models trained and packed, a rerun, a reduced
database, another seed, the samples profiled as they are simulated, scenarios, failures and SIGTERM.

Needs $PROTAL, $SIMULATE (simulate_metagenomes) and a Python with scikit-learn, joblib, numpy and pandas
($PROTAL_TRAIN_PYTHON, default this one); without them skipped, or failed with PROTAL_TESTS_REQUIRED=1
(scripts/prerequisites.py). A few minutes. The build script's parts on their own: test_gtdb_build.py.

  PROTAL=build/protal SIMULATE=build/simulate_metagenomes PROTAL_TRAIN_PYTHON=python3 \
      python3 -m unittest scripts/mini_db/test_gtdb_pipeline.py
"""

import collections
import csv
import glob
import gzip
import json
import os
import random
import re
import shutil
import signal
import subprocess
import sys
import tempfile
import time
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mini_db_fixtures import BUILD, DOWNLOAD, LINEAGES, SIMULATE, FakeNcbi, run  # noqa: E402
import prerequisites  # noqa: E402

import build_gtdb_database as build  # noqa: E402
import build_gtdb_releases as releases  # noqa: E402


class GtdbBuildTest(unittest.TestCase):
    """build_gtdb_database.py end to end, on a synthetic GTDB-like release of 60 species downloaded from stand-ins of
    GTDB's mirror and NCBI. Each pipeline run serves every check of what it does:
      test_a  the database built with every gene, its pe and se models trained, the genes ranked from the training
              database (--rank-genes; full_build, which the other tests share); a rerun builds nothing; a reduced
              database of the 3 best genes ranked from a full build of the training database, the same ranking as
              --rank-genes's, both databases built again from the samples already simulated, its models with the
              relatives features and calls at a target share of false calls
      test_b  a reduced database ranked by a given table (--gene-ranking), whose finished database's build fails in the
              background: the run stops at once
      test_c  a copy of the full build with another seed: the finished database kept, the training database built
              again, the samples simulated again; stopped by SIGTERM, the run leaves no command running
      test_f  scenarios with host reads, the feature sets chosen by the trainers, the reads behind the models' errors
      test_g  the full build's samples profiled as they are simulated: the same tables
    Needs $PROTAL, $SIMULATE and a Python with scikit-learn ($PROTAL_TRAIN_PYTHON, default this one)."""

    full = None  # full_build()'s run, once made

    @classmethod
    def setUpClass(cls):
        cls.python = os.environ.get("PROTAL_TRAIN_PYTHON", sys.executable)
        missing = []
        if not prerequisites.executable(os.environ.get("PROTAL", "")):
            missing.append("$PROTAL (a protal binary)")
        if not prerequisites.executable(os.environ.get("SIMULATE", "")):
            missing.append("$SIMULATE (simulate_metagenomes)")
        try:
            trainer_ok = subprocess.run([cls.python, "-c", "import joblib, numpy, pandas, sklearn"],
                                        capture_output=True).returncode == 0
        except OSError:
            trainer_ok = False
        if not trainer_ok:
            missing.append(f"scikit-learn, joblib, numpy and pandas in {cls.python} ($PROTAL_TRAIN_PYTHON)")
        if missing:
            prerequisites.missing("needs " + ", ".join(missing))
        cls.tmp = tempfile.TemporaryDirectory()
        cls.full = None
        work = cls.tmp.name
        lineages = os.path.join(work, "lineages.txt")
        with open(lineages, "w") as fh:
            subprocess.run([sys.executable, LINEAGES, "--species", "60", "--archaea", "0.1", "--seed", "1"], stdout=fh,
                           check=True)
        gtdb = os.path.join(work, "gtdb")
        run(SIMULATE, "--outdir", gtdb, "--lineages", lineages, "--genomes_per_species", "3", "--genome_length", "40000",
            "--strain_divergence", "0.002-0.015", "--species_divergence", "0.015-0.04", "--seed", "3")
        ncbi = FakeNcbi(gtdb, os.path.join(work, "served"))
        cls.inputs = os.path.join(work, "inputs")
        try:
            subprocess.run([sys.executable, DOWNLOAD, "-o", cls.inputs, *ncbi.download_options, "--species", "15",
                            "--per_species", "2", "--rep_only_species", "10", "--batch", "8", "-t", "2", "--host_genome",
                            "none"], env=ncbi.env, check=True, capture_output=True)
        finally:
            ncbi.stop()

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    @classmethod
    def build(cls, out, *extra, protal=None, wait=True, scenarios=False, error_reads="none", foreign=True):
        """A run of build_gtdb_database.py into tmp/out; its console is also kept as tmp/out.log. The default scenarios
        (scenarios.py) have the depths of real studies, so the runs but test_f's (which defines small ones) have none;
        only test_f keeps the reads of the models' errors; zstd level 1 for both databases, whose content is the same at
        any level. The runs scan the full references for the foreign rates (--foreign-rates, off by default), so that
        their tables compare; test_f has the default, no scan."""
        command = [cls.python, BUILD, "--inputs", cls.inputs, "--outdir", os.path.join(cls.tmp.name, out),
                   "--protal", protal or os.environ["PROTAL"], "--simulator", os.environ["SIMULATE"], "-t", "2",
                   "--samples", "2", "--read-pairs", "1000,4000", "--read-setups", "100:HS20:300:40",
                   "--species-per-sample", "6-8", "--archaea", "1", "--holdout-max-share", "0.2",
                   "--holdout-clades", "family:1,genus:1", "--read-types", "pe,se", "--test-samples", "1",
                   "--test-read-pairs", "2000", "--ntree", "16", "--rounds", "40", "--evaluation", "basic",
                   "--progress-every", "5", "--training-db-level", "1", "--final-db-level", "1",
                   "--error-reads", error_reads, *([] if scenarios else ["--scenarios", "none"]),
                   *(["--foreign-rates"] if foreign else []), *extra]
        if not wait:
            return subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=3000)
        with open(os.path.join(cls.tmp.name, out + ".log"), "a") as fh:
            fh.write(result.stdout)
        return result

    @classmethod
    def full_build(cls):
        """The build of every gene into out (its samples on the scratch disk tmp/scratch, both collections profiled in one
        protal run, the genes ranked: --rank-genes), made once for the class -> its run. Its tables are copied to
        full_tables (test_g compares its own with them), and out and scratch as they are to full_copy (test_c runs
        another seed there): test_a builds out again."""
        if cls.full is None:
            scratch = os.path.join(cls.tmp.name, "scratch")
            cls.full = cls.build("out", "--scratch", scratch, "--profile-blocks", "0", "--rank-genes")
            if cls.full.returncode == 0:
                for collection in ("training", "test"):
                    os.makedirs(os.path.join(cls.tmp.name, "full_tables", collection))
                    for table in ("training_data.tsv", "training_data_se.tsv"):
                        shutil.copy(os.path.join(cls.tmp.name, "out", "work", collection, table),
                                    os.path.join(cls.tmp.name, "full_tables", collection))
                shutil.copytree(os.path.join(cls.tmp.name, "out"), os.path.join(cls.tmp.name, "full_copy", "out"),
                                symlinks=True)
                shutil.copytree(scratch, os.path.join(cls.tmp.name, "full_copy", "scratch"), symlinks=True)
        return cls.full

    def text(self, *path):
        with open(os.path.join(self.tmp.name, *path)) as fh:
            return fh.read()

    def metadata(self, out):
        with open(os.path.join(self.tmp.name, out, "protal_db", "build_metadata.tsv")) as fh:
            return dict(line.rstrip("\n").split("\t", 1) for line in fh)

    def samples(self, out, collection, table):
        """{scenario ("" for the design): sample names} of a collection's table (OUTDIR/work/<collection>), and its rows."""
        with open(os.path.join(self.tmp.name, out, "work", collection, table)) as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        samples = collections.defaultdict(set)
        for row in rows:
            samples[row.get("meta_scenario", "")].add(row["meta_sample"])
        return samples, rows

    @staticmethod
    def gene_list(path):
        """The gene ids of a gene_subset.txt, in its order."""
        with open(path) as fh:
            return [int(line) for line in fh if line.strip() and not line.startswith("#")]

    def summary_row(self, out, variant):
        """build_gtdb_releases.py's summary of a run into out (BuildReleasesTest runs the script itself)."""
        return releases.summarize("226", variant, os.path.join(self.tmp.name, out), "ok", 60,
                                  os.path.join(self.tmp.name, out + ".log"))

    def test_a_build_rerun_and_reduced_database(self):
        first = self.full_build()
        self.assertEqual(first.returncode, 0, first.stdout[-3000:])
        out = os.path.join(self.tmp.name, "out")
        scratch = ("--scratch", os.path.join(self.tmp.name, "scratch"), "--profile-blocks", "0")
        for path in ("protal_db/database.protal", "model_logs/summary.txt",
                     "work/stages/convert.json", "work/stages/protal_db.json", "work/stages/training_db.json"):
            self.assertTrue(os.path.isfile(os.path.join(out, path)), path)
        # OUTDIR holds the database, the evaluation, the logs and what a rerun reuses, each file once: protal_db/ the
        # database and what names its genes and genomes (the build's reports beside it moved to model_logs/, the foreign
        # scan's table to work/), the models' files in model_logs/ only, the in-silico strains on the scratch disk.
        self.assertEqual(sorted(os.listdir(out)), ["console.log", "logs", "model_logs", "protal_db", "work"])
        self.assertEqual(sorted(os.listdir(os.path.join(out, "protal_db"))),
                         ["build_metadata.tsv", "database.protal", "gene2geneid.tsv", "genome2tiid.tsv"])
        for name in ("trained_model.xml", "trained_model_se.report.txt", "trained_model.calls.tsv.gz", "parity.txt",
                     "genome_table.txt", "heldout_species.txt", "gene_incongruence.tsv"):
            self.assertTrue(os.path.isfile(os.path.join(out, "model_logs", name)), name)
        self.assertEqual(glob.glob(os.path.join(out, "model_logs", "*.log")) +
                         glob.glob(os.path.join(out, "model_logs", "build_metadata.tsv")), [])
        self.assertIn("(se reads)", self.text("out", "model_logs", "parity.txt"))  # every read type's check in one file
        self.assertTrue(os.path.isfile(os.path.join(self.tmp.name, "scratch", "insilico_strains", "insilico_strains.tsv")))
        self.assertFalse(os.path.exists(os.path.join(out, "work", "insilico_strains")))
        self.assertRegex(first.stdout, r"In \S+: protal_db/ the database \(4 files\), model_logs/ the evaluation \(\d+\)")
        # The training database, read only by the collections and the parity check, is built on the scratch disk.
        self.assertTrue(os.path.isfile(os.path.join(self.tmp.name, "scratch", "training_db", "database.protal")))
        self.assertFalse(os.path.exists(os.path.join(out, "work", "training_db")))
        # full_reference.fna, which only the builds read, is gone once they are done.
        for path in ("out/protal_db/full_reference.fna", "scratch/training_db/full_reference.fna"):
            for name in (path, path + ".zst"):
                self.assertFalse(os.path.exists(os.path.join(self.tmp.name, name)), name)
        # The genes' conservation factors are in the database, and their summary in build_metadata.tsv.
        metadata = self.metadata("out")
        self.assertRegex(metadata["gene_conservation"], r"^factors [0-9.]+-[0-9.]+ for \d+ genes, from \d+ species")
        self.assertEqual(metadata["classifier_previous_procedure"], "not compared")  # without --previous-procedure
        # The versions the run started with; nothing changed them during it.
        self.assertRegex(metadata["protal_version"], r"^protal v[0-9.]+")
        self.assertNotIn("at the end of the run", metadata["protal_version"] + metadata["scripts_commit"])
        self.assertNotIn("changed during the run", first.stdout)
        # The models train on the default feature set, evaluated with samples, species and clades held out (the build's
        # defaults: --features DEFAULT_FEATURE_SET, --evaluation basic); test_f has each trainer choose (--features auto);
        # the reduced database below trains the relatives features and the calls at a target share of false calls.
        import model_features
        self.assertEqual((metadata["classifier_features"], metadata["classifier_evaluation"]),
                         (model_features.DEFAULT_FEATURE_SET, "basic"))
        self.assertNotIn("model_pe_features", metadata)
        self.assertNotIn("Feature sets chosen", self.text("out", "model_logs", "summary.txt"))
        # The samples' composition against the simulator's truth, with the models' calls (composition_accuracy.py): a
        # table per read type of the training and test samples; the species' genome sizes are the simulated genomes'
        # (the release's metadata gives their lengths), the summary goes into summary.txt and the console.
        for t, name in (("pe", "composition_accuracy.tsv"), ("se", "composition_accuracy_se.tsv")):
            with open(os.path.join(out, "model_logs", name)) as fh:
                rows = list(csv.DictReader(fh, delimiter="\t"))
            self.assertEqual({row["set"] for row in rows}, {"training", "test"}, t)
            self.assertEqual({row["read_type"] for row in rows}, {t})
            sizes = [float(row["size_ratio"]) for row in rows if row["size_ratio"] != "NA"]
            self.assertTrue(sizes, t)
            self.assertLess(abs(sorted(sizes)[len(sizes) // 2] - 1), 0.05, t)
            self.assertTrue(all(row["unknown"] != "NA" and 0 <= float(row["unknown"]) <= 1 for row in rows), t)
        self.assertIn("Composition of the pe samples against their truth", self.text("out", "model_logs", "summary.txt"))
        self.assertRegex(first.stdout, r"\d+/\d+ checking the samples' composition against the truth")
        self.assertRegex(first.stdout, r"pe: median errors over \d+ test samples without a host: explained share [-+]")
        self.assertIn(f"--features {model_features.DEFAULT_FEATURE_SET} --model gbm",
                      self.text("out", "logs", "classifier_training_se.log"))
        self.assertEqual(metadata["classifier_scenarios"], "none")
        self.assertIn("gene copies", metadata["suspect_copies"])  # the build looked for suspect copies
        # The gene copies' gaps to their congeners' copies (--build), in the databases the training samples were profiled
        # with. The scan of each database's full reference for the gene copies other species' reads reach (--foreign-rates,
        # which build() passes; test_f has the default, no scan): the training database's right after its build, from its
        # full reference, which lacks the held-out species; the finished database's before its models went in; every copy
        # listed; the full references gone after.
        self.assertRegex(self.text("out", "logs", "index_and_package.log"), r"Congener gaps: \d+ gene copies of \d+ species")
        # The strain alleles (--strain-alleles 4, --allele-genome-share 0.5): both builds took them from the genomes of the
        # allele share, which the genome table lost (none of them simulated).
        for log in ("index_and_package.log", "training_db_index.log"):
            self.assertRegex(self.text("out", "logs", log), r"Strain alleles: (\d+ alleles \(\d+ edits\) of \d+ gene copies of \d+ "
                                                            r"species, up to 4 each|none) \(of \d+ full-reference copies")
            self.assertIn("outside --allele_genome_share 0.5", self.text("out", "logs", log))
        self.assertRegex(first.stdout, r"genome table \(genomes\.tsv, genome_table\.txt\): .*; \d+ genomes left to the strain "
                                       r"alleles, not simulated \(--allele-genome-share 0\.5\)")
        alleles = re.search(r"Strain alleles: (\d+) alleles", self.text("out", "logs", "training_db_index.log"))
        for stage, folder in (("foreign_rates_training", os.path.join(self.tmp.name, "scratch", "training_db")),
                              ("foreign_rates", os.path.join(out, "work", "foreign_rates"))):
            self.assertRegex(self.text("out", "logs", stage + ".log"),
                             r"Tiles: \d+ reads of 150 bases every 250 bases of \d+ of \d+ records \(\d+ headers, at most 10 "
                             r"records each\) of .*full_reference\.fna")
            self.assertRegex(self.text("out", "logs", stage + ".log"),
                             r"Foreign rates: \d+ reads with a record, \d+ counted \(MAPQ >= 4\), 0 of an unknown taxon, from "
                             r"\d+ gene copies of \d+ species' genes \(at most 10 each\); (\d+) gene copies listed, \d+ reached")
            self.assertTrue(os.path.isfile(os.path.join(out, "work", "stages", stage + ".json")), stage)
            self.assertTrue(os.path.isfile(os.path.join(folder, "foreign_rates.tsv")), stage)
            # --add_tables in the same log, after the scan's lines
            self.assertRegex(self.text("out", "logs", stage + ".log"), r"\n--- \S+ --add_tables (.|\n)*gene copies of")
        training_table = self.text("scratch", "training_db", "foreign_rates.tsv").splitlines()
        final_table = self.text("out", "work", "foreign_rates", "foreign_rates.tsv").splitlines()
        self.assertLess(len(training_table), len(final_table))  # the held-out species' copies only in the finished one
        self.assertTrue(all(":" in line for line in final_table[2:]))
        self.assertIn("the foreign scan", self.text("out", "console.log"))
        # The hold-out steered by the species' clouds (species_clouds.tsv, protal --write_species_neighbours on the
        # converted release): heldout_species.txt with the distance to the nearest kept congener and the complex, and
        # holdout.txt with the bands, all in model_logs; the clouds also in the collector's command (meta_novel_distance).
        self.assertRegex(self.text("out", "logs", "species_clouds.log"), r"Species neighbours: \d+ pairs of congeners compared")
        self.assertEqual(self.text("out", "model_logs", "species_clouds.tsv").splitlines()[0], "taxid\tneighbours")
        held = [line.split("\t") for line in self.text("out", "model_logs", "heldout_species.txt").splitlines()]
        self.assertTrue(held and all(len(f) == 5 for f in held), held)
        self.assertTrue(all(f[3] == "-" or float(f[3]) <= 0.15 for f in held), held)
        self.assertIn("nearest kept congener of the species held out alone", self.text("out", "model_logs", "holdout.txt"))
        self.assertIn("species complexes (congeners within 0.01) held out whole", self.text("out", "model_logs", "holdout.txt"))
        self.assertRegex(metadata["classifier_training_holdout_clouds"], r"^\d+ species' congeners within 0\.15 \(the converted "
                                                                         r"release\); \d+ complexes of \d+ species within 0\.01")
        self.assertRegex(self.text("out", "logs", "training_data.log"), r"species clouds: \d+ species' nearest congeners")
        self.assertIn("; congeners 0.25:2-5", metadata["classifier_training_design"])
        commands = []
        for path in glob.glob(os.path.join(self.tmp.name, "scratch", "**", "run_params.tsv"), recursive=True):
            with open(path) as fh:
                commands.append(fh.read())
        self.assertTrue(commands, "the simulators' run_params.tsv")  # beside the samples, in --scratch
        self.assertTrue(all("--congener_groups 0.25:2-5" in c for c in commands))
        # The species with one genome (10 of the download are representatives only) got in-silico strains, and the
        # collections simulated from them, too.
        self.assertRegex(metadata["insilico_strains"], r"^\d+ in-silico strains of the \d+ species with one genome")
        self.assertTrue(any(line.startswith("insilico_")
                            for line in self.text("out", "work", "genomes_simulated.tsv").splitlines()))
        self.assertFalse(any(line.startswith("insilico_") for line in self.text("out", "work", "genomes.tsv").splitlines()))
        self.assertIn("with the in-silico strains", self.text("out", "model_logs", "genome_table.txt"))
        self.assertTrue(all("genomes_simulated.tsv" in c for c in commands))
        self.assertEqual(metadata["classifier_call_mode"], "curve")
        self.assertNotIn("model_pe_false_calls", metadata)
        self.assertNotIn("--fdr-calls", self.text("out", "logs", "classifier_training.log"))
        # What the conservation features rest on, on the release's genomes: how the genes differ between congeners
        # (protal --build), and where the held-out species' reads land (trace_relatives.py).
        self.assertRegex(metadata["gene_congeners"], r"^\d+ pairs of species of \d+ genera")
        # The gene neighbours of every genome to simulate from: the frequencies and the positions in the database,
        # the training database's frequencies derived anew without the species it leaves out.
        self.assertRegex(metadata["gene_neighbours"], r"^\d+ rules of \d+ clades from \d+ genomes")
        self.assertRegex(metadata["gene_positions"], r"^\d+ genes \(\d+ placed by their k-mer trace\) in \d+ genomes of "
                                                     r"\d+ species, \d+ read as circular")
        self.assertTrue(os.path.isfile(os.path.join(out, "logs", "gene_neighbours.log")))
        self.assertRegex(self.text("out", "logs", "training_db.log"), r"gene_neighbours\.tsv derived from the \d+ genomes of "
                                                             r"\d+ species kept")
        for name in ("gene_congeners.tsv", "relatives_by_gene_conservation.txt"):
            self.assertTrue(os.path.isfile(os.path.join(out, "model_logs", name)), name)
        # The held-out species' reads were traced to genes (through their genomes' contigs; at r226 v10 none was).
        traced = self.text("out", "model_logs", "relatives_by_gene_conservation.txt")
        self.assertNotIn("Not traced", traced)
        self.assertRegex(traced, r"over [1-9]\d* genes")
        # The training data: 2 samples at each of the 2 depths, pe and se, each table with present and absent taxa.
        for table in ("training_data.tsv", "training_data_se.tsv"):
            samples, rows = self.samples("out", "training", table)
            self.assertEqual(len(samples[""]), 4, table)
            self.assertEqual({r["truth"] for r in rows}, {"0", "1"}, table)
            # The per-copy tables were loaded for the profiling: the features are known (not -1) where reads landed, the
            # foreign rates among them (the scan's table in the training database); meta_novel_distance is a distance
            # or empty.
            for feature in ("gap_informative_share", "untried_candidate_rate", "foreign_scanned_share"):
                self.assertTrue(any(float(r[feature]) >= 0 for r in rows), f"{table}: {feature}")
            # The strain alleles' features: known wherever the training database has the table.
            self.assertEqual(any(float(r["allele_copy_share"]) >= 0 for r in rows), alleles is not None, table)
            self.assertTrue(all(r["meta_novel_distance"] == "" or 0 <= float(r["meta_novel_distance"]) <= 0.15 for r in rows), table)
        # A stage running for a while (5 s here, --progress-every) says how it is doing.
        self.assertRegex(first.stdout, r": \d+:\d\d:\d\d so far")
        # The models go into the database in one rewrite, each with its knob curve over depth (the trainer's
        # --depth-knobs for every read type), gradient-boosted trees of up to 63 leaves (--model, --maxnodes); the
        # converter logs its steps' times.
        self.assertEqual(glob.glob(os.path.join(out, "logs", "final_package_*.log")), [])
        self.assertEqual(self.text("out", "logs", "final_package.log").count("Models for read types:"), 1)
        self.assertEqual(metadata["classifier_depth_knobs"], "pe,se")
        self.assertEqual((metadata["classifier_model"], metadata["classifier_trees"]), ("gbm", "40 rounds"))
        self.assertEqual(metadata["classifier_max_leaves"], "pe:63,se:63")
        self.assertIn("--model gbm --ntree 16 --maxnodes 63", self.text("out", "logs", "classifier_training_se.log"))
        self.assertIn("gradient-boosted trees: 40 rounds", self.text("out", "logs", "classifier_training_se.log"))
        self.assertRegex(self.text("out", "logs", "convert.log"), r"spooled the representatives' marker genes \(\d+ species\): [\d.]+ s")
        self.assertRegex(self.text("out", "logs", "convert.log"), r"joined them into full_reference\.fna(\.zst)?: [\d.]+ s")
        # The training database is built alone; then the finished database in the background at the idle scheduling
        # class, and both collections' simulations, what to stream chosen from the room on the samples' disk (here all
        # of it is room: nothing streamed). The models go into the database before the reports.
        self.assertLess(first.stdout.index("built training_db in"),
                        first.stdout.index("building protal_db meanwhile, in the background at the idle scheduling class"))
        self.assertLess(first.stdout.index("built training_db in"),
                        first.stdout.index("simulating the training data and the independent test set in the background"))
        self.assertRegex(first.stdout, r"room on \S+scratch: [\d.]+ [MG]B free; besides the reads the run needs")
        self.assertRegex(first.stdout, r"--profile-blocks 0: none of the \d+ simulations streamed")
        self.assertLess(first.stdout.index("Ready protal database"), first.stdout.index("Reports of what the models' errors"))
        self.assertTrue(self.text("out", "logs", "protal_runs_training.log").startswith("==> all <==\n"),
                        "the protal run's log copied off the scratch disk")
        simulation = self.text("out", "logs", "training_data_simulation.log")
        self.assertRegex(simulation, r"2 of 2 paired-end design points, \d+:\d\d:\d\d in all")
        self.assertIn("profiling left to a run without --simulate_only", simulation)
        self.assertNotIn("protal profiled", simulation)
        collection = self.text("out", "logs", "training_data.log")
        self.assertNotRegex(collection, "simulating")
        # The test set's samples (1 pe, 1 se) are profiled in the training data's protal run: the database loads once.
        self.assertIn("and 2 samples of another collection in one protal run", collection)
        self.assertRegex(collection, r"protal profiled 10 samples in \d+:\d\d:\d\d")
        self.assertNotIn("protal profiled", self.text("out", "logs", "test_data.log"))
        self.assertRegex(self.text("out", "logs", "test_data.log"), r"\d+ taxa in \S+training_data\.tsv: \d+ present")
        # The genome table has each genome's length, as the simulator counts it when the table has none.
        with open(os.path.join(out, "work", "genomes.tsv")) as fh:
            table = [line.rstrip("\n").split("\t") for line in fh]
        self.assertTrue(table and all(len(f) == 4 and f[3].isdigit() for f in table))
        three = os.path.join(self.tmp.name, "three_columns.tsv")
        with open(three, "w") as fh:
            fh.write("".join("\t".join(f[:3]) + "\n" for f in table))
        subprocess.run([os.environ["SIMULATE"], "--genome_table", three, "-o", os.path.join(self.tmp.name, "lengths"),
                        "-n", "3", "--total_read_pairs", "100", "--species_per_sample", "6", "--seed", "4", "--test"],
                       check=True, capture_output=True)
        with open(os.path.join(self.tmp.name, "lengths", "manifest.tsv")) as fh:
            counted = {row["genome"]: row["genome_length"] for row in csv.DictReader(fh, delimiter="\t")}
        self.assertTrue(counted)
        self.assertEqual(counted, {f[0]: f[3] for f in table if f[0] in counted})
        self.assertTrue(os.path.isfile(os.path.join(out, "work", "training", "training_data_se.tsv")))
        self.assertFalse(os.path.exists(os.path.join(out, "work", "training", "points")))
        self.assertTrue(os.path.isdir(os.path.join(self.tmp.name, "scratch", "training", "points")))
        self.assertRegex(first.stdout, r"The run took at most [\d.]+ [MG]B on \S+scratch")
        # --rank-genes: the genes of the training database ranked (unpacked on the scratch disk, then removed).
        with open(os.path.join(out, "model_logs", "gene_ranking.tsv")) as fh:
            ranking_text = fh.read()
        ranking = [line.split("\t") for line in ranking_text.splitlines()]
        self.assertEqual(ranking[0][:4], ["rank", "gene_id", "marker", "score"])
        self.assertEqual(len(ranking) - 1, len(build.read_gene_ids(os.path.join(out, "work", "gene2geneid.tsv"))))
        self.assertFalse(os.path.exists(os.path.join(self.tmp.name, "scratch", "ranking_files")))
        # build_gtdb_releases.py's summary of the build: the release's species (not only those simulated from), the
        # numbers it reads from the console, the models' scores.
        full = self.summary_row("out", "full")
        self.assertEqual((full["genes"], full["species"]), ("all", 60))
        self.assertRegex(full["protal_version"], r"^protal v[0-9.]+")
        self.assertTrue(full["genera"] <= full["species"] and full["phyla"] >= 1, full)
        self.assertEqual(full["bacteria_species"] + full["archaea_species"], 60)
        self.assertGreater(full["archaea_species"], 0)
        self.assertRegex(full["genomes_simulated"], r"^\d+ genomes of \d+ species$")
        self.assertRegex(str(full["species_held_out"]), r"^\d+$")
        self.assertRegex(full["database_gb"], r"^\d+\.\d\d$")
        self.assertRegex(full["build_time"], r"^\d+:\d\d:\d\d$")
        self.assertRegex(full["build_peak_gb"], r"^\d+\.\d+$")
        self.assertRegex(full["profiling_peak_gb"], r"^\d+\.\d+$")
        for t in ("pe", "se"):
            self.assertRegex(full[f"{t}_test_F1"], r"^[0-9.]+$|^-$")
            self.assertRegex(full[f"{t}_heldout_F1"], r"^[0-9.]+$")
            self.assertRegex(full[f"{t}_heldout_FP_per_sample"], r"^[0-9.]+$")
        self.assertEqual((full["pb_test_F1"], full["ont_test_F1"]), ("", ""))
        described = releases.describe(full)
        self.assertRegex(described, r"taxa: 60 species \(\d+ bacteria, \d+ archaea\), \d+ genera, \d+ families, \d+ orders, "
                                    r"\d+ classes, \d+ phyla")
        self.assertRegex(described, r"\n  pe: independent test F1 [0-9.-]+, FP/sample [0-9.-]+, sensitivity [0-9.-]+, "
                                    r"precision [0-9.-]+; species held out F1 [0-9.]+, FP/sample [0-9.]+")

        # A rerun (here without the models' evaluation) converts and builds nothing: the stages, the training database,
        # the in-silico strains and the ranking stay as they were (the finished database gets its models again), and the
        # collector reuses its samples and dumps.
        training_db = os.path.join(self.tmp.name, "scratch", "training_db", "database.protal")
        kept = [os.path.join(out, "work", "stages", f"{stage}.json") for stage in ("convert", "protal_db", "training_db")] + \
            [training_db, os.path.join(out, "work", "genomes_simulated.tsv"),
             os.path.join(out, "model_logs", "gene_ranking.tsv"),
             os.path.join(self.tmp.name, "scratch", "insilico_strains", "insilico_strains.tsv")]
        before = {path: os.stat(path).st_mtime_ns for path in kept}
        again = self.build("out", *scratch, "--rank-genes", "--evaluation", "none")
        self.assertEqual(again.returncode, 0, again.stdout[-3000:])
        self.assertEqual({path: os.stat(path).st_mtime_ns for path in kept}, before)
        self.assertNotRegex(again.stdout, r"built (protal|training)_db")
        self.assertNotRegex(self.text("out", "logs", "training_data_simulation.log"), "simulating")
        self.assertNotRegex(self.text("out", "logs", "training_data.log"), "simulating|profiling")

        # A reduced database: the 3 best genes ranked from a full build of the training database (--n-genes 3 without
        # --gene-ranking; test_b gives one), one of them a gene archaea have. The ranking is --rank-genes's above:
        # the release ranked by build_gtdb_releases.py's two runs gets the same genes as by one. Both database folders
        # are derived from the release converted anew, both databases built again from the samples already simulated,
        # the whole conversion and the ranking database gone at the end. Its models with the relatives features and calls
        # at a target share of false calls: trained, checked for parity with protal, and in the database.
        built = [os.path.join(out, "work", "stages", "protal_db.json"), training_db]  # marked once built, and the build
        reduced = self.build("out", *scratch, "--n-genes", "3", "--features", "normalized+adjacency+relatives",
                             "--call-mode", "fdr")
        self.assertEqual(reduced.returncode, 0, reduced.stdout[-3000:])
        self.assertTrue(all(os.stat(path).st_mtime_ns != before[path] for path in built), "both databases built again")
        self.assertNotRegex(self.text("out", "logs", "training_data_simulation.log"), "simulating")
        self.assertFalse(os.path.exists(os.path.join(out, "work", "converted")))
        self.assertFalse(os.path.exists(os.path.join(self.tmp.name, "scratch", "ranking_db")))
        for path in ("logs/gene_ranking_build.log", "logs/gene_ranking_files.log", "work/stages/gene_ranking.json",
                     "model_logs/gene_ranking.tsv", "model_logs/gene_subset.txt"):
            self.assertTrue(os.path.isfile(os.path.join(out, path)), path)
        self.assertNotIn("genes kept", self.text("out", "logs", "gene_ranking_files.log"))  # every gene, the species left out
        # The same ranking as --rank-genes's, the congeners' columns too (gene_congeners.tsv, beside the training
        # database's database.protal, is copied to the unpacked folder that --rank-genes ranks).
        self.assertEqual(self.text("out", "model_logs", "gene_ranking.tsv"), ranking_text)
        self.assertRegex(ranking[1][ranking[0].index("between_factor")], r"^[0-9.]+$")  # the best gene's, a bacterial one
        subset_text = self.text("out", "model_logs", "gene_subset.txt")
        self.assertIn("ranked from a full build of the training database", subset_text.splitlines()[0])
        chosen = self.gene_list(os.path.join(out, "model_logs", "gene_subset.txt"))
        self.assertEqual(len(chosen), 3)
        self.assertIn(int(ranking[1][1]), chosen)  # the best overall, and the best of each domain
        header = ranking[0]
        archaeal = {int(r[1]) for r in ranking[1:] if float(r[header.index("archaea_prevalence")]) >= 0.5}
        self.assertTrue(set(chosen) & archaeal, "no gene archaea have among the three chosen")
        self.assertIn("chosen for archaea", subset_text)
        metadata = self.metadata("out")
        self.assertRegex(metadata["marker_genes"], r"^3 of \d+, the most distinctive by prevalence x unique k-mer share "
                                                   r"\(scripts/rank_genes\.py, 1 per domain, from a full build of the training "
                                                   r"database\): \S+, \S+, \S+; in half the species or more of bacteria \d, archaea \d$")
        self.assertRegex(self.text("out", "logs", "training_db.log"), r"derived from the \d+ genomes of \d+ species kept: "
                                                             r"\d+ lines \(over the 3 genes kept\)")
        self.assertIn(", 3 genes kept: ", self.text("out", "logs", "protal_db_files.log"))
        self.assertEqual(metadata["classifier_features"], "normalized+adjacency+relatives")
        self.assertEqual(metadata["classifier_call_mode"], "fdr")
        for t in ("pe", "se"):
            self.assertRegex(metadata[f"model_{t}_false_calls"], r"^target [0-9.]+ per sample; species held out F1 [0-9.]+")
            log = self.text("out", "logs", "classifier_training" + ("" if t == "pe" else "_se") + ".log")
            self.assertIn("--fdr-calls", log)
            self.assertIn("em_own_share", log)  # among the model's features
        row = self.summary_row("out", "n3")
        self.assertEqual(row["genes"], "3")
        self.assertIn("1 per domain", row["marker_genes"])
        self.assertLess(float(row["database_gb"]), float(full["database_gb"]) + 0.01)

        # The training database's folder on the scratch disk says what it was built from (built_for.json, as
        # work/stages/training_db.json; the rerun above kept the database by it). Another OUTDIR's build on the same scratch
        # disk (other species held out) leaves its own there: this OUTDIR's next run then builds the training database
        # again instead of taking that one. Stopped once it has decided.
        stamp = os.path.join(self.tmp.name, "scratch", "training_db", "built_for.json")
        with open(stamp) as fh:
            self.assertEqual(json.load(fh), json.loads(self.text("out", "work", "stages", "training_db.json")))
        with open(stamp, "w") as fh:
            json.dump({"heldout": "another OUTDIR's species"}, fh)
        console = os.path.join(out, "console.log")
        seen = os.path.getsize(console)
        rerun = self.build("out", *scratch, "--n-genes", "3", "--features", "normalized+adjacency+relatives",
                           "--call-mode", "fdr", wait=False)
        said, deadline = "", time.time() + 600
        while rerun.poll() is None and os.path.exists(stamp) and time.time() < deadline and \
                "with the same species left out" not in said:
            time.sleep(0.2)
            with open(console) as fh:
                fh.seek(seen)
                said = fh.read()
        if rerun.poll() is None:
            rerun.send_signal(signal.SIGTERM)
        rerun.communicate(timeout=120)
        self.assertFalse(os.path.exists(stamp), "the training database of another build was kept")
        self.assertNotIn("with the same species left out", said)

    def test_b_a_failed_background_build_stops_the_run(self):
        # A reduced database ranked by a given table (--n-genes 3 --gene-ranking, as build_gtdb_releases.py builds its
        # reduced variants), here of four genes without the domains' columns: the table copied, the subset its 3 best,
        # no ranking database built. The finished database's build, started in the background once the training database
        # is built, then fails at once: the run stops within seconds, before any model is trained (at this scale the
        # collection may have ended by then; at GTDB scale it takes hours).
        given = os.path.join(self.tmp.name, "given_ranking.tsv")
        with open(given, "w") as fh:
            fh.write("rank\tgene_id\tmarker\tscore\n" +
                     "".join(f"{rank}\t{gene}\tgiven{gene}\t{1 - rank / 10:.1f}\n" for rank, gene in enumerate((9, 4, 6, 2), 1)))
        failing = os.path.join(self.tmp.name, "failing_protal")
        with open(failing, "w") as fh:
            fh.write(f'#!/bin/sh\ncase "$*" in *--build*protal_db*) exit 3;; esac\nexec {os.environ["PROTAL"]} "$@"\n')
        os.chmod(failing, 0o755)
        started = time.time()
        result = self.build("out_fail", "--insilico-strains", "0", "--no-gene-neighbours", "--n-genes", "3",
                            "--gene-ranking", given, protal=failing)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Command failed (3)", result.stdout)
        self.assertIn("index_and_package.log", result.stdout)
        out = os.path.join(self.tmp.name, "out_fail")
        self.assertNotRegex(result.stdout, r"\d+/\d+ training the ")
        self.assertEqual(glob.glob(os.path.join(out, "model_logs", "trained_model*.xml")), [])
        self.assertLess(time.time() - started, 600)
        self.assertEqual(self.text("out_fail", "model_logs", "gene_ranking.tsv"), self.text("given_ranking.tsv"))
        self.assertEqual(self.gene_list(os.path.join(out, "model_logs", "gene_subset.txt")), [9, 4, 6])
        self.assertIn(f"ranked from --gene-ranking {given}", self.text("out_fail", "model_logs", "gene_subset.txt"))
        for name in ("logs/gene_ranking_build.log", "logs/gene_ranking_files.log", "work/stages/gene_ranking.json",
                     "work/ranking_db"):
            self.assertFalse(os.path.exists(os.path.join(out, name)), name)
        # A run without the gene subset into the same folder: protal_db holds the subset's files, derived from the release
        # converted whole into work/converted, which the stopped run marked. They must not pass for the whole release (the
        # finished database would have the subset's genes): the release is converted again. Stopped once it says.
        console = os.path.join(out, "console.log")
        seen = os.path.getsize(console)
        full = self.build("out_fail", "--insilico-strains", "0", "--no-gene-neighbours", wait=False)
        said, deadline = "", time.time() + 600
        while full.poll() is None and time.time() < deadline and \
                not any(s in said for s in ("converting GTDB", "holds the release converted by an earlier run")):
            time.sleep(0.5)
            with open(console) as fh:
                fh.seek(seen)
                said = fh.read()
        if full.poll() is None:
            full.send_signal(signal.SIGTERM)
        full.communicate(timeout=120)
        self.assertIn("converting GTDB", said)
        self.assertNotIn("holds the release converted by an earlier run", said)

    def test_c_another_seed_then_sigterm(self):
        # The full build's copy (full_copy) run with another seed: other species held out, so the finished database is
        # kept and the training database built again (its stage marked anew), from the release converted anew, and
        # every design point simulated again rather than mixed into the tables. Stopped by SIGTERM in its collection,
        # the run leaves no command running.
        first = self.full_build()
        self.assertEqual(first.returncode, 0, first.stdout[-3000:])
        out = os.path.join(self.tmp.name, "full_copy", "out")
        stages = os.path.join(out, "work", "stages")
        final = [os.path.join(out, "protal_db", "database.protal"), os.path.join(stages, "protal_db.json")]
        before = {path: os.stat(path).st_mtime_ns for path in final + [os.path.join(stages, "training_db.json")]}
        with open(os.path.join(out, "model_logs", "heldout_species.txt")) as fh:
            held_out = fh.read()
        run_ = self.build(os.path.join("full_copy", "out"), "--scratch", os.path.join(self.tmp.name, "full_copy", "scratch"),
                          "--profile-blocks", "0", "--rank-genes", "--seed", "2", wait=False)
        simulation = os.path.join(out, "logs", "training_data_simulation.log")

        def resimulated():
            try:
                with open(simulation) as fh:
                    return "was simulated from other inputs (or by an older collector): simulating it again" in fh.read()
            except OSError:
                return False

        def training_db_built():
            path = os.path.join(stages, "training_db.json")
            return os.path.exists(path) and os.stat(path).st_mtime_ns != before[path]
        deadline = time.time() + 900
        while not (resimulated() and training_db_built()) and run_.poll() is None and time.time() < deadline:
            time.sleep(0.5)
        self.assertIsNone(run_.poll(), "the run ended before its collection")
        self.assertTrue(resimulated(), "the design points were not simulated again")
        self.assertTrue(training_db_built(), "the training database was not built again")
        self.assertEqual({path: os.stat(path).st_mtime_ns for path in final}, {path: before[path] for path in final})
        with open(os.path.join(out, "model_logs", "heldout_species.txt")) as fh:
            self.assertNotEqual(fh.read(), held_out)
        time.sleep(1)  # the simulations go on in the background, the run ranks the genes or waits for them
        run_.send_signal(signal.SIGTERM)
        output, _ = run_.communicate(timeout=120)
        self.assertNotEqual(run_.returncode, 0)
        self.assertIn("Stopped by SIGTERM", output)
        self.assertNotRegex(output, r"built protal_db")
        # The in-silico strains of the new seed take each gene's divergence from the real strains' (gene_positions.tsv),
        # which the finished database kept packed: the release is converted again for them, not left to a uniform ANI.
        self.assertIn("the release again, for its gene positions", output)
        self.assertNotIn("no gene_positions.tsv", output)
        left = []
        for pid in os.listdir("/proc"):
            try:
                with open(f"/proc/{pid}/cmdline", "rb") as fh:
                    if os.path.join(self.tmp.name, "full_copy").encode() in fh.read():
                        left.append(pid)
            except OSError:
                pass
        self.assertEqual(left, [], "commands of the stopped run are still running")

    def scenario_inputs(self):
        """test_f's small scenarios (definitions of the presets' names) and host genome. -> their paths."""
        work = os.path.join(self.tmp.name, "scenario_inputs")
        definitions, host = os.path.join(work, "scenarios.json"), os.path.join(work, "host.fna.gz")
        if os.path.isfile(definitions):
            return definitions, host
        os.makedirs(work, exist_ok=True)
        rng = random.Random(9)
        with gzip.open(host, "wt") as fh:
            for c in (1, 2):
                seq = "".join(rng.choice("ACGT") for _ in range(150000))
                fh.write(f">chr{c}\n" + "\n".join(seq[i:i + 80] for i in range(0, len(seq), 80)) + "\n")
        illumina = {"type": "pe", "length": 100, "profile": "HS20", "fragment_mean": 300, "fragment_sd": 40, "quality": 30}
        ultima = {"type": "se", "setup": "ultima:300:40:25:2"}
        small = {"abundance": "lognormal:1.5", "strains": "0.3", "congeners": "0", "host_share": 0}
        with open(definitions, "w") as fh:
            json.dump({"gut": {**small, "species": "5-6", "novel_share": 0.3,
                               "reads": [{**illumina, "depth": 3000}, {**ultima, "depth": 1000}]},
                       "moderate": {**small, "species": "10-15", "novel_share": 0.3,
                                    "reads": [{**illumina, "depth": 1500}]},
                       "soil": {**small, "species": "30-40", "novel_share": 0.6, "reads": [{**illumina, "depth": 2000}]},
                       "soil_shallow": {**small, "species": "5", "novel_share": 0.3, "reads": [{**illumina, "depth": 1000}]},
                       # host at one depth, so that its host reads can be counted; the others drawn around theirs
                       "host": {**small, "species": "2-3", "novel_share": 0.3, "abundance": "powerlaw:1.0", "strains": "",
                                "host_share": 0.9, "depth_spread": 1,
                                "reads": [{**illumina, "depth": 20000}, {**ultima, "depth": 5000}]}}, fh)
        return definitions, host

    def test_f_scenarios(self):
        # The default scenarios (gut, moderate, soil, soil_shallow, host; here made small by a --scenario-file of their
        # names),
        # one with 90% host reads: their hold-in samples in the training data, their hold-out samples in the test set,
        # each read type's report and the summary scoring both; soil, larger than the genome table holds at 60% held
        # out, scaled down; the feature sets chosen by the trainers (--features auto, not the default) and why; the
        # reads behind the models' errors (--error-reads all, the build's default), their SAMs and the archive to share
        # (--share-logs). Both collections profiled in one protal run, their reads kept (--profile-blocks 0). No scan of
        # the full references for the foreign rates (the default; test_a has --foreign-rates): no table, the features
        # unknown, and --features auto still tries the candidate set with them.
        import compressed
        definitions, host = self.scenario_inputs()
        scratch = os.path.join(self.tmp.name, "scenario_scratch")
        result = self.build("scenarios", "--scenario-file", definitions, "--scenario-samples", "2",
                            "--scenario-test-samples", "1", "--host-genome", host, "--scratch", scratch,
                            "--profile-blocks", "0", "--features", "auto", "--share-logs",
                            scenarios=True, error_reads="all", foreign=False)
        self.assertEqual(result.returncode, 0, result.stdout[-3000:])
        for stage in ("foreign_rates_training", "foreign_rates"):
            self.assertFalse(os.path.exists(os.path.join(self.tmp.name, "scenarios", "logs", stage + ".log")), stage)
            self.assertFalse(os.path.exists(os.path.join(self.tmp.name, "scenarios", "work", "stages", stage + ".json")),
                             stage)
        self.assertFalse(os.path.exists(os.path.join(scratch, "training_db", "foreign_rates.tsv")))
        _, rows = self.samples("scenarios", "training", "training_data.tsv")
        self.assertTrue(all(float(r["foreign_scanned_share"]) == -1 for r in rows))
        simulation = self.text("scenarios", "logs", "training_data_simulation.log")
        self.assertRegex(simulation, r"scenario gut \(2 samples\): \d+ species \(\d+ the database lacks, [\d.]+%; \d+ it "
                                     r"has\) for samples of 5-6; Illumina reads at Q30 \(HS20, --mean_quality\)")
        self.assertIn("sc_host_pe_p20000 simulated (2 samples)", simulation)
        # A host sample: the community's 2,000 read pairs and the host's 18,000; its Ultima reads, 90% of them host's.
        points = os.path.join(scratch, "training", "points")
        host_reads = glob.glob(os.path.join(points, "sc_host_pe_p20000", "sim", "reads", "*_R1.fq.zst"))
        self.assertEqual(len(host_reads), 2)  # zstd by default (--read-compression)
        # Every read file the build simulated is zstd (the build's default --read-compression: the design's paired-end
        # points by simulate_metagenomes --reads_compression zstd, the scenarios' with their host reads in the same run,
        # Ultima and long reads by simulate_metagenomes --long_samples), and protal profiled them all (the tables below).
        for collection in ("training", "test"):
            simulated = glob.glob(os.path.join(scratch, collection, "points", "*", "sim", "reads", "*.fq*"))
            self.assertTrue(any("/rl100_p" in p.replace(os.sep, "/") for p in simulated), collection)  # the design's
            for path in simulated:
                self.assertTrue(path.endswith(".fq.zst"), path)
                with open(path, "rb") as fh:
                    self.assertEqual(fh.read(4), compressed.ZSTD_MAGIC, path)
        for reads in host_reads:
            names = compressed.read_text(reads).splitlines()[0::4]
            self.assertEqual(sum(n.startswith("@h_") for n in names), 18000)  # exactly the host's
            self.assertEqual(len(names), 20000)  # and exactly the community's 2,000
        lines = compressed.read_text(glob.glob(os.path.join(points, "sc_host_se_ultima_r5000", "sim", "reads",
                                                            "*.fq.zst"))[0]).splitlines()
        lengths = [len(r) for r in lines[1::4]]
        self.assertAlmostEqual(sum(lengths) / len(lengths), 300, delta=15)
        quality = [ord(c) - 33 for q in lines[3::4] for c in q]
        self.assertAlmostEqual(sum(quality) / len(quality), 25, delta=1.5)
        # Its tables: the scenarios' rows beside the design's, in the training data (which the models are fitted on)
        # and in the test set; the training data's design samples as without scenarios (4 pe, 4 se).
        for collection, n in (("training", 2), ("test", 1)):
            for table, names in (("training_data.tsv", ("gut", "moderate", "soil", "soil_shallow", "host")),
                                 ("training_data_se.tsv", ("gut", "host"))):
                samples, _ = self.samples("scenarios", collection, table)
                self.assertEqual({k: len(v) for k, v in samples.items() if k}, {name: n for name in names})
                self.assertEqual(len(samples[""]), 4 if collection == "training" else 1)
        # Each model's report scores the scenarios, hold-in and hold-out, and every candidate feature set on every test
        # set; and so does the summary.
        for t, names in (("", ("gut", "moderate", "soil", "soil_shallow", "host")), ("_se", ("gut", "host"))):
            report = self.text("scenarios", "model_logs", f"trained_model{t}.report.txt")
            self.assertIn("## Scenarios: hold-in and hold-out samples", report)
            self.assertIn("## Feature set chosen (--features auto, species held out)", report)
            self.assertRegex(report, r"scenario host: hold-out F1 [0-9.-]+ \(FP rate [0-9.-]+%, FN rate [0-9.-]+%, "
                                     r"1 samples\); hold-in with species held out F1")
            with open(os.path.join(self.tmp.name, "scenarios", "model_logs", f"trained_model{t}.metrics.json")) as fh:
                metrics = json.load(fh)
            sets = {(r["scenario"], r["set"]) for r in metrics["scenarios"]}
            for name in names:
                self.assertTrue({(name, "hold-out"), (name, "hold-in, in sample"),
                                 (name, "hold-in, species held out")} <= sets, sets)
            auto = metrics["features_auto"]
            self.assertIn(auto["chosen"], [r["features"] for r in auto["candidates"]])
            self.assertTrue(any("foreign" in r["features"].split("+") for r in auto["candidates"]))
            self.assertEqual(set(auto["held_out"]), {"test set"} | {f"{name} hold-out" for name in names})
            self.assertTrue(all(f"F1 {name} hold-out" in r for r in auto["candidates"] for name in names))
        summary = self.text("scenarios", "model_logs", "summary.txt")
        for label, count in (("gut: hold-out", 2), ("gut: hold-in, species held out", 2), ("host: hold-in, in sample", 2),
                             ("soil: hold-out", 1), ("soil_shallow: hold-out", 1)):
            self.assertEqual(summary.count(label), count, label)
        self.assertIn("FP rate  FN rate", summary)
        self.assertRegex(summary, r"Feature sets chosen \(--features auto\):\n  pe: normalized\S*: F1 ")
        metadata = self.metadata("scenarios")
        self.assertTrue(metadata["classifier_scenarios"].startswith(
            "gut 2 hold-in (training) and 1 hold-out (test) samples, moderate 2 hold-in (training) and 1 hold-out"))
        self.assertIn("; soil scaled from 30-40 to ", metadata["classifier_scenarios"])
        self.assertRegex(metadata["model_pe_features"], r"^normalized\S* \(--features auto: F1 ")
        self.assertRegex(metadata["model_se_scenarios"], r"gut hold-out F1 ([0-9.]+|-), FP rate ([0-9.]+%|-)")
        # The reads behind each model's errors: every sample profiled with an unmapped record for each read that seeded
        # but aligned nowhere (no counts in the SAM header), and per sample of the training data and the test set, of
        # each read type, the records of the errors' reads, each with its source.
        unmapped = 0
        for collection in ("training", "test"):
            for sam in glob.glob(os.path.join(scratch, collection, "points", "*", "protal*", "alignments", "*.sam*")):
                if not sam.endswith((".sam", ".sam.gz", ".sam.zst")):
                    continue
                lines = compressed.read_text(sam).splitlines()
                unmapped += sum(not line.startswith("@") and line.split("\t")[1] == "4" for line in lines)
                self.assertFalse(any(line.startswith("@CO\tprotal failed candidates") for line in lines), sam)
        self.assertGreater(unmapped, 0)
        logs = os.path.join(self.tmp.name, "scenarios", "model_logs")
        for t in ("pe", "se"):
            # The errors are those of the model's calls (the trainer's calls.tsv.gz, copied to model_logs): every sample
            # of both tables.
            expected, samples_of_calls = collections.Counter(), set()
            with gzip.open(os.path.join(logs, f"trained_model{'' if t == 'pe' else '_' + t}.calls.tsv.gz"), "rt") as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    samples_of_calls.add((row["set"], row["meta_sample"]))
                    if row["call"] != row["truth"]:
                        expected[(row["set"], row["meta_sample"], "FP" if row["call"] == "1" else "FN")] += 1
            errors = os.path.join(logs, "error_reads", t)
            with open(os.path.join(errors, "summary.tsv")) as fh:
                samples = list(csv.DictReader(fh, delimiter="\t"))
            self.assertEqual({(r["set"], r["sample"]) for r in samples}, samples_of_calls, t)
            self.assertEqual({r["scenario"] for r in samples},
                             {"design", "gut", "host"} | ({"moderate", "soil", "soil_shallow"} if t == "pe" else set()), t)
            # Every sample's error taxa in one table (no file per sample), with its design point and scenario.
            with gzip.open(os.path.join(errors, "taxa.tsv.gz"), "rt") as fh:
                taxa = list(csv.DictReader(fh, delimiter="\t"))
            self.assertEqual(glob.glob(os.path.join(errors, "*", "*", "*.taxa.tsv")), [])
            point_of = {(r["set"], r["sample"]): (r["point"], r["scenario"]) for r in samples}
            self.assertTrue(all(point_of[(x["set"], x["sample"])] == (x["point"], x["scenario"]) for x in taxa), t)
            for r in samples:
                for kind in ("FP", "FN"):
                    self.assertEqual(int(r[kind]), expected[(r["set"], r["sample"], kind)], (t, r["sample"], kind))
                # The records of the FP's and the FN's reads in files of their own (no unseen ones), QUAL left out,
                # each taken for a reason of its file's kind.
                written, kinds = 0, set()
                for path in filter(None, r["sams"].split(",")):
                    kind = path.rsplit(".", 3)[1]
                    self.assertIn(kind, ("FP", "FN"), path)
                    kinds.add(kind)
                    lines = compressed.read_text(os.path.join(errors, path)).splitlines()
                    records = [line.split("\t") for line in lines if not line.startswith("@")]
                    written += len(records)
                    self.assertTrue(all(f[10] == "*" for f in records), path)
                    tags = [{t[:2]: t[5:] for t in f[11:] if t[:2] in ("xg", "xs", "xe")} for f in records]
                    self.assertTrue(all(len(t) == 3 for t in tags), path)
                    if kind == "FP":
                        self.assertTrue(all("FP:" in t["xe"] for t in tags), path)
                    self.assertTrue(any(line.startswith("@CO\terror_reads.py: ") for line in lines))
                self.assertEqual(written, int(r["sam_records"]), (t, r["sample"]))
                self.assertLessEqual(int(r["sam_fragments"]), int(r["fragments"]))
                if int(r["FP"]):
                    self.assertIn("FP", kinds, (t, r["sample"]))  # an FP taxon has a row: reads on it
                # Every read's source genome known (a paired-end read is named after its contig, a drawn read after
                # its genome's place), but the host's paired-end reads.
                if r["scenario"] != "host" or t != "pe":
                    self.assertEqual(r["unknown_source"], "0", (t, r["sample"]))
                errors_of = collections.Counter(x["error"] for x in taxa
                                                if (x["set"], x["sample"]) == (r["set"], r["sample"]))
                self.assertEqual((errors_of["FP"], errors_of["FN"], errors_of["unseen"]),
                                 (int(r["FP"]), int(r["FN"]), int(r["unseen"])))
            self.assertGreater(sum(int(r["records"]) for r in samples), 0, t)
            if any(int(r["FP"]) + int(r["FN"]) for r in samples):
                self.assertGreater(sum(int(r["sam_records"]) for r in samples), 0, t)
            # Which side the errors' reads take where the species differs from its congeners (ancestry_sites.py, after
            # the error reads, from the training database's genes and its full reference, kept until then): a summary,
            # the AUC table and the counts per taxon and record.
            ancestry = os.path.join(logs, "ancestry_sites", t)
            if not any(int(r["sam_records"]) for r in samples):
                continue
            summary = self.text("scenarios", "model_logs", "ancestry_sites", f"{t}.summary.txt")
            self.assertRegex(summary, r"^\d+ SAMs of \d+ samples, \d+ counted records", t)
            self.assertIn("alleles: ", summary)  # the full reference was there
            self.assertIn("Per taxon, FN own (1) against FP genus (0)", summary)
            with open(ancestry + ".auc.tsv") as fh:
                self.assertEqual(fh.readline().rstrip("\n").split("\t"), ["min_sites", "identity_band", "signal", "taxa", "fn", "auc"])
            for suffix in (".taxa.tsv.gz", ".fragments.tsv.gz"):
                with gzip.open(ancestry + suffix, "rt") as fh:
                    self.assertTrue(fh.readline().startswith("sample\t"), suffix)
            self.assertTrue(os.path.isfile(os.path.join(self.tmp.name, "scenarios", "logs", f"ancestry_sites_{t}.log")))
        console = self.text("scenarios", "console.log")
        self.assertIn("the ancestry sites of the errors' reads (model_logs/ancestry_sites", console)
        self.assertIn("full_reference.fna", console)  # kept for the report, then removed
        self.assertFalse(glob.glob(os.path.join(scratch, "training_db", "full_reference.fna*")))
        self.assertFalse(os.path.exists(os.path.join(scratch, "ancestry_files")))
        # The archive to share (--share-logs), laid out as OUTDIR: the console's lines and the other logs, model_logs/
        # with the error reads' SAMs but not the models, the build's metadata, the taxonomy, and the tables, their numbers
        # to 9 significant digits (the same rows and values to 1e-8).
        import math
        import tarfile
        with tarfile.open(os.path.join(self.tmp.name, "scenarios", "scenarios_share.tar.gz")) as tar:
            members = set(tar.getnames())
            for name in ("console.log", "logs/training_data.log", "model_logs/summary.txt",
                         "model_logs/error_reads/pe/summary.tsv", "model_logs/error_reads/pe/taxa.tsv.gz",
                         "model_logs/ancestry_sites/pe.summary.txt", "model_logs/ancestry_sites/pe.auc.tsv",
                         "protal_db/build_metadata.tsv", "work/internal_taxonomy.dmp", "work/training/training_data.tsv",
                         "work/training/training_data_se.tsv", "work/test/training_data.tsv"):
                self.assertIn("scenarios/" + name, members)
            self.assertTrue(any(m.endswith(".FP.sam.zst") or m.endswith(".FN.sam.zst") for m in members))
            self.assertFalse(any(m.endswith((".partial", ".xml", ".joblib")) or "share_tables" in m for m in members))
            self.assertFalse(any(m.startswith("scenarios/protal_db/database") for m in members))
            self.assertIn("Ready protal database", tar.extractfile("scenarios/console.log").read().decode())
            shared = tar.extractfile("scenarios/work/training/training_data.tsv").read().decode().splitlines()
        with open(os.path.join(scratch, "training", "training_data.tsv")) as fh:
            original = fh.read().splitlines()
        self.assertEqual(len(shared), len(original))
        for a, b in zip(original, shared):
            for x, y in zip(a.split("\t"), b.split("\t"), strict=True):
                if x != y:
                    self.assertTrue(math.isclose(float(x), float(y), rel_tol=1e-8), (x, y))
                    self.assertLessEqual(len(y), len(x))

        # Without the host genome, the default scenarios run without host and say so; asked for, host stops the build
        # before anything runs.
        default = self.build("scenarios_default_nohost", "--scenario-file", definitions, scenarios=True, wait=False)
        seen = ""
        for line in default.stdout:
            seen += line
            if "left out: no host genome" in line:
                break
        default.send_signal(signal.SIGTERM)
        default.communicate(timeout=120)
        self.assertIn("WARNING: scenario host left out: no host genome (--host-genome, or rerun download_gtdb.py", seen)
        stopped = self.build("scenarios_nohost", "--scenarios", "host", "--scenario-file", definitions, scenarios=True)
        self.assertNotEqual(stopped.returncode, 0)
        self.assertIn("its host reads need the host genome", stopped.stdout)

    def test_g_profiled_as_simulated(self):
        # The full build's samples profiled as they are simulated (--profile-blocks, here a few kB: a protal run as soon
        # as anything is simulated), the two collections' runs taking turns: the same tables as full_build()'s one
        # protal run; every read removed once profiled (the SAMs, profiles and dumps kept). Its models are not evaluated:
        # what the run profiled is in the tables.
        first = self.full_build()
        self.assertEqual(first.returncode, 0, first.stdout[-3000:])
        scratch = os.path.join(self.tmp.name, "follow_scratch")
        result = self.build("follow", "--scratch", scratch, "--profile-blocks", "0.00001", "--keep-free", "1",
                            "--evaluation", "none")
        self.assertEqual(result.returncode, 0, result.stdout[-3000:])
        self.assertIn("profiled as its simulations go on, in protal runs of 1e-05 GB of reads or more", result.stdout)
        for collection in ("training", "test"):
            for table in ("training_data.tsv", "training_data_se.tsv"):
                self.assertEqual(self.text("follow", "work", collection, table), self.text("full_tables", collection, table),
                                 f"{collection}/{table}")
            points = os.path.join(scratch, collection, "points")
            self.assertEqual(glob.glob(os.path.join(points, "*", "sim", "reads", "*.fq.*")), [])
            self.assertTrue(glob.glob(os.path.join(points, "*", "sim", "reads_removed.txt")))
            self.assertTrue(glob.glob(os.path.join(points, "*", "protal*", "alignments", "*.sam*")))
        log = self.text("follow", "logs", "training_data.log")
        self.assertRegex(log, r"protal run 1: \d+ design points, [\d.]+ GB of reads; the simulations (go on|have ended)")
        self.assertRegex(log, r"every design point is profiled, in \d+ protal runs?")
        self.assertIn("keeping 1 GB free on", self.text("follow", "logs", "training_data_simulation.log"))

    def test_h_streamed(self):
        # Every design point streamed (--stream-above, here a few bytes): protal reads each point's samples from named
        # pipes as simulate_metagenomes makes them (pe, then se from a run that writes only read 1s), as plain FASTQ, a
        # protal run per point, no read on the disk; the genomes from a genome store on the scratch disk
        # (--genome-store auto): the same tables as full_build()'s one protal run.
        first = self.full_build()
        self.assertEqual(first.returncode, 0, first.stdout[-3000:])
        scratch = os.path.join(self.tmp.name, "stream_scratch")
        result = self.build("stream", "--scratch", scratch, "--profile-blocks", "20", "--stream-above", "0.000000001",
                            "--evaluation", "none", "--genome-store", "auto")
        self.assertEqual(result.returncode, 0, result.stdout[-3000:])
        self.assertTrue(glob.glob(os.path.join(scratch, "genome_store", "*.g2b")), "the genomes in the store")
        for collection in ("training", "test"):
            for table in ("training_data.tsv", "training_data_se.tsv"):
                self.assertEqual(self.text("stream", "work", collection, table), self.text("full_tables", collection, table),
                                 f"{collection}/{table}")
            points = os.path.join(scratch, collection, "points")
            self.assertEqual(glob.glob(os.path.join(points, "*", "sim", "reads", "*")), [], "no reads, no pipes left")
            self.assertEqual(glob.glob(os.path.join(points, "*", "stream_*", "")), [], "the streams' folders removed")
            self.assertTrue(glob.glob(os.path.join(points, "*", "stream_pe.log")), "the simulators' logs kept")
            self.assertTrue(glob.glob(os.path.join(points, "*", "streamed.json")))
            self.assertTrue(glob.glob(os.path.join(points, "*", "protal*", "alignments", "*.sam*")))
        log = self.text("stream", "logs", "training_data.log")
        self.assertRegex(log, r"\d+ simulations streamed into protal \(a sample above 1e-09 GB\)")
        # The streamed points whose communities are there share a protal run (--profile-block-max, 200 GB): here both.
        self.assertRegex(log, r"protal run 1: 2 simulations streamed \(rl100_p1000, rl100_p4000; 4 design points, \d+ "
                              r"samples, read from named pipes")
        self.assertRegex(log, r"every design point is profiled, in 1 protal run\b")
        self.assertIn("--stream-above 1e-09: 3 of the 3 simulations streamed into protal", result.stdout)


if __name__ == "__main__":
    unittest.main()
