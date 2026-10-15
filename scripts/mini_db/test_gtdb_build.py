#!/usr/bin/env python3
"""test_gtdb_build.py - checks for build_gtdb_database.py, build_gtdb_releases.py, rank_genes.py and the reduced
(gene subset) database, short of a build: the build script's options and functions, its check of the binaries' commit
(git, cmake), build_gtdb_releases.py with a stand-in build script, and the converter's gene subsets; with $PROTAL (a
protal binary) protal --build of a gene subset and rank_genes.py on a full build. A missing prerequisite skips its
test, or fails it with PROTAL_TESTS_REQUIRED=1 (scripts/prerequisites.py). The build end to end: test_gtdb_pipeline.py.

  PROTAL=build/protal python3 -m unittest scripts/mini_db/test_gtdb_build.py
"""

import argparse
import collections
import contextlib
import csv
import gzip
import io
import json
import os
import re
import shutil
import subprocess
import sys
import tarfile
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
import model_pmml  # noqa: E402

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

    def test_kinds_of_models(self):
        # --model: the kinds to train, the first the database's default; --maxnodes for it, the others their kind's
        # leaves. The finished database's placeholders: the models at other gates, then every model of the other kinds
        # (model_pe_forest.xml, model_pe_c1_forest.xml); refreshing them removes an earlier build's of other gates and
        # kinds and puts a build_metadata.tsv where there is none (packed last, replaced in place at the end).
        self.assertEqual(build.model_kinds_spec("gbm,forest"), ["gbm", "forest"])
        self.assertEqual(build.model_kinds_spec("forest,"), ["forest"])
        for bad in ("", "gbm,gbm", "rf", "gbm,rf"):
            with self.assertRaises(argparse.ArgumentTypeError):
                build.model_kinds_spec(bad)
        args = argparse.Namespace(model="gbm", maxnodes="31")
        self.assertEqual((build.kind_leaves(args, "gbm"), build.kind_leaves(args, "forest")), ("31", build.FOREST_LEAVES))
        args = argparse.Namespace(model="forest", maxnodes=build.FOREST_LEAVES)
        self.assertEqual(build.kind_leaves(args, "gbm"), build.BOOSTED_LEAVES)
        self.assertEqual(build.model_placeholders(["pe", "se"], {"pe": [1, 8]}, ["forest"]),
                         [("model_pe_c1.xml", "pe"), ("model_pe_c8.xml", "pe"), ("model_pe_forest.xml", "pe"),
                          ("model_pe_c1_forest.xml", "pe"), ("model_pe_c8_forest.xml", "pe"), ("model_se_forest.xml", "se")])
        self.assertEqual(build.model_placeholders(["pe"], {}, []), [])
        with tempfile.TemporaryDirectory() as db:
            for name in ("model_pe.xml", "model_pe_c4.xml", "model_se_gbm.xml", "model_pe_c1_forest.xml"):
                with open(os.path.join(db, name), "w") as fh:
                    fh.write("an earlier build's")
            build.refresh_model_placeholders(db, [("model_pe_c1_forest.xml", "pe"), ("model_se_forest.xml", "se")])
            self.assertEqual(sorted(os.listdir(db)), ["build_metadata.tsv", "model_pe.xml", "model_pe_c1_forest.xml",
                                                      "model_se_forest.xml"])
            self.assertEqual(model_pmml.model_kind(os.path.join(db, "model_se_forest.xml")), "placeholder")
            with open(os.path.join(db, "build_metadata.tsv")) as fh:
                self.assertEqual(fh.read(), build.METADATA_PLACEHOLDER)
            with open(os.path.join(db, "build_metadata.tsv"), "w") as fh:
                fh.write("gtdb_release\tr226\n")
            build.refresh_model_placeholders(db, [])
            with open(os.path.join(db, "build_metadata.tsv")) as fh:
                self.assertEqual(fh.read(), "gtdb_release\tr226\n", "an earlier build's is kept until the end")
        # The kind of a model as protal tells it (profiler::ModelKindOf): by its first segmentation.
        with tempfile.TemporaryDirectory() as tmp:
            for method, kind in (("modelChain", "gbm"), ("average", "forest"), ("majorityVote", "forest"), ("sum", "")):
                path = os.path.join(tmp, f"{method}.xml")
                with open(path, "w") as fh:
                    fh.write(f'<PMML><Header/><MiningModel><Segmentation multipleModelMethod="{method}">'
                             '<Segment/></Segmentation></MiningModel></PMML>')
                self.assertEqual(model_pmml.model_kind(path), kind, method)
        self.assertEqual(model_pmml.kind_model_file(model_pmml.gate_model_file("pb", 32), "gbm"), "model_PB_c32_gbm.xml")

    def test_the_wfa_traceback_reaches_both_builds(self):
        # --wfa-traceback: main() sets TRACEBACK_ARGS, which every protal --build carries (the training database's and the
        # finished one's), so that both record it in build_settings.tsv; unset, no option (protal's default, singletrack since 0.8.2), so that
        # earlier builds' commands and keys stay as they were.
        saved = build.TRACEBACK_ARGS[:]
        try:
            build.TRACEBACK_ARGS[:] = []
            self.assertNotIn("--wfa_traceback", build.build_command("protal", "/nonexistent/db", 4))
            build.TRACEBACK_ARGS[:] = ["--wfa_traceback", "singletrack"]
            for alleles in (None, build.FINAL_ALLELE_ARGS):
                command = build.build_command("protal", "/nonexistent/db", 4, "--compress_level", "3", alleles=alleles)
                self.assertEqual(command[command.index("--wfa_traceback") + 1], "singletrack")
                self.assertEqual(command[-2:], ["--compress_level", "3"])
        finally:
            build.TRACEBACK_ARGS[:] = saved

    def test_sparse_first_gates(self):
        # --sparse-first: a gate per read type (protal --sparse_first's form), one for all, or none; --test-c-gates: the
        # collector's --sparse_first_levels form, its default every read type's.
        self.assertEqual(build.sparse_first_spec(build.DEFAULT_SPARSE_FIRST), {"pe": 4, "se": 4, "pb": 16, "ont": 8})
        self.assertEqual(build.sparse_first_spec("2"), {"pe": 2, "se": 2, "pb": 2, "ont": 2})
        self.assertEqual(build.sparse_first_spec("none"), {})
        self.assertEqual(build.gates_text({"pe": 1, "pb": 8}), "pe:1,pb:8")
        self.assertEqual(build.gates_text({}), "none")
        for bad in ("3", "pe:64", "xx:1", "pe:1,", "pe"):
            with self.assertRaises(argparse.ArgumentTypeError, msg=bad):
                build.sparse_first_spec(bad)
        self.assertEqual(build.sparse_first_levels(build.DEFAULT_TEST_C_GATES),
                         {"pe": [0, 2, 16], "se": [0, 2, 16], "pb": [0, 4], "ont": [0, 4, 32]})
        self.assertEqual(build.sparse_first_levels("pb:16/4/4"), {"pb": [4, 16]})
        # --test-flex-indels: ONT's gates and the Ultima units' (key ultima) with --flex_indels too; none: no such rows.
        self.assertEqual(build.test_flex_spec(build.DEFAULT_TEST_FLEX_INDELS), {"ont": [0, 1, 4, 8], "ultima": [0, 1, 2, 4]})
        self.assertEqual(build.test_flex_spec("none"), {})
        with self.assertRaises(argparse.ArgumentTypeError):
            build.test_flex_spec("illumina:1")
        # --test-other-pick: auto (its default with --test-c-gates), each read type's database gate if 2 or more, with
        # the pick the database does not use (gates below 2 have none); given gates; none: no rows.
        default = build.sparse_first_levels(build.DEFAULT_TEST_C_GATES)
        gates = build.sparse_first_spec(build.DEFAULT_SPARSE_FIRST)
        self.assertEqual(build.test_stride_levels(build.test_stride_spec("auto"), default, gates),
                         {"pe": [4], "se": [4], "pb": [16], "ont": [8]})
        self.assertEqual(build.test_stride_levels("auto", {"pe": [0]}, {"pe": 1, "pb": 8}), {"pb": [8]})
        self.assertEqual(build.test_stride_levels("auto", None, gates), {})
        self.assertEqual(build.test_stride_levels(build.test_stride_spec("pb:8"), default, gates), {"pb": [8]})
        self.assertEqual((build.other_pick("stride"), build.other_pick("hash")), ("hash", "stride"))
        # --test-grid: auto (its default with --test-c-gates), each read type's gates with a model of 2 or more (the
        # ladder: --model-gates and the database's); given gates; none: no rows.
        ladder = {"pe": [1, 4, 8], "se": [1], "pb": [8, 16, 32], "ont": [0]}
        self.assertEqual(build.test_grid_levels(build.test_stride_spec("auto"), default, ladder),
                         {"pe": [4, 8], "pb": [8, 16, 32]})
        self.assertEqual(build.test_grid_levels("auto", None, ladder), {})
        self.assertEqual(build.test_grid_levels(build.test_stride_spec("pe:2"), default, ladder), {"pe": [2]})
        self.assertEqual(build.test_grid_levels(build.test_stride_spec("none"), default, ladder), {})
        # --sparse-first-pick: stride (the default) or hash for all, or read types and theirs (the others stride).
        self.assertEqual(build.sparse_first_pick_spec(build.DEFAULT_SPARSE_FIRST_PICK),
                         {"pe": "stride", "se": "stride", "pb": "stride", "ont": "stride"})
        picks = build.sparse_first_pick_spec("pb:hash")
        self.assertEqual(picks, {"pe": "stride", "se": "stride", "pb": "hash", "ont": "stride"})
        self.assertEqual(build.picks_text(picks), "pe:stride,se:stride,pb:hash,ont:stride")
        self.assertEqual(build.sparse_first_pick_spec("grid")["ont"], "grid")
        self.assertEqual(build.sparse_first_pick_spec("pb:grid")["pb"], "grid")
        for bad in ("strided", "pb:8", "xx:hash", "pe:"):
            with self.assertRaises(argparse.ArgumentTypeError, msg=bad):
                build.sparse_first_pick_spec(bad)
        # --model-gates: a ladder per read type, the database's gate added; none: the database's alone. A run gated at C
        # is profiled with the model of the gate nearest it (log2 C, off below 1, the lower of two).
        ladder = build.model_gate_ladder(build.model_gates_spec(build.DEFAULT_MODEL_GATES), gates, ["pe", "se", "pb", "ont"])
        self.assertEqual(ladder, {"pe": [1, 4, 8], "se": [1, 4, 8], "pb": [8, 16, 32], "ont": [1, 8, 16]})
        self.assertEqual(build.model_gate_ladder(build.model_gates_spec("pe:1/8"), {"pe": 2}, ["pe", "se"]),
                         {"pe": [1, 2, 8], "se": [0]})
        self.assertEqual(build.model_gates_spec("none"), {})
        self.assertEqual([build.nearest_gate(c, [1, 4, 8]) for c in (0, 1, 2, 4, 16, 32)], [1, 1, 1, 4, 8, 8])
        self.assertEqual([build.nearest_gate(c, [0, 4]) for c in (0, 1, 2, 8)], [0, 0, 4, 4])
        self.assertEqual(build.test_stride_spec("none"), {})
        for bad in ("pe:1", "ultima:4"):
            with self.assertRaises(argparse.ArgumentTypeError, msg=bad):
                build.test_stride_spec(bad)
        self.assertEqual(build.levels_text({"pe": [0, 2], "pb": [4]}), "pe:0/2,pb:4")
        for bad in ("pb:3", "pb", "xx:1", "pe:64"):
            with self.assertRaises(argparse.ArgumentTypeError, msg=bad):
                build.sparse_first_levels(bad)

    def test_flex_indels_support_from_the_full_help(self):
        # --test-flex-indels needs a protal that runs --flex_indels per sample: its --full_help lists the option and the
        # map's FLEX_INDELS column. An older one, which would ignore the column, does not; nor one that cannot run.
        with tempfile.TemporaryDirectory() as tmp:
            for n, (text, supported) in enumerate((("--flex_indels TYPES; map column FLEX_INDELS", True),
                                                   ("--flex_indels TYPES", False), ("--sparse_first C", False))):
                path = os.path.join(tmp, f"protal{n}")
                with open(path, "w") as fh:
                    fh.write(f"#!/bin/sh\necho '{text}'\n")
                os.chmod(path, 0o755)
                self.assertEqual(build.flex_indels_supported(path), supported, text)
            self.assertFalse(build.flex_indels_supported(os.path.join(tmp, "none")))

    def test_stride_support_from_the_full_help(self):
        # --test-stride needs a protal that runs the stride pick per sample: its --full_help lists the option and the
        # map's column (STRIDE_HELP). An older one, which would ignore the column, does not; nor one that cannot run.
        with tempfile.TemporaryDirectory() as tmp:
            for n, (text, supported) in enumerate(((" ".join(build.STRIDE_HELP), True), (build.STRIDE_HELP[0], False),
                                                   ("--sparse_first C; --stage1_probe stride", False))):
                path = os.path.join(tmp, f"protal{n}")
                with open(path, "w") as fh:
                    fh.write(f"#!/bin/sh\necho '{text}'\n")
                os.chmod(path, 0o755)
                self.assertEqual(build.stride_supported(path), supported, text)
            self.assertFalse(build.stride_supported(os.path.join(tmp, "none")))

    def test_the_genomes_of_the_strain_alleles_are_not_simulated(self):
        # The genomes that give strain alleles: a hash of the accession, the same in protal (tests/test_StrainAlleles.cpp
        # pins the same values); the genome table to simulate from loses them, but the representatives.
        from gtdb_to_protal_db import allele_genome
        self.assertTrue(allele_genome("GCF_000005845.2", 0.5))   # 0.278
        self.assertTrue(allele_genome("GCA_000001405.29", 0.5))  # 0.090
        self.assertFalse(allele_genome("GCA_900000000.1", 0.5))  # 0.696
        self.assertFalse(allele_genome("GCF_000195955.2", 0.5))  # 0.508
        self.assertTrue(allele_genome("GCF_000195955.2", 0.51))
        self.assertTrue(allele_genome("GCA_900000000.1", 1))
        self.assertFalse(allele_genome("GCF_000005845.2", 0))
        with tempfile.TemporaryDirectory() as root:
            table = os.path.join(root, "genomes.tsv")
            with open(table, "w") as fh:
                fh.write("GCF_000005845.2\td__Bacteria;s__A\t/a.fna\t10\n"      # an allele genome, but the representative
                         "RS_GCA_000001405.29\td__Bacteria;s__A\t/b.fna\t10\n"  # an allele genome
                         "GCA_900000000.1\td__Bacteria;s__A\t/c.fna\t10\n"      # a strain
                         "GCF_000195955.2\td__Bacteria;s__B\t/d.fna\t10\n")
            kept = os.path.join(root, "kept.tsv")
            out, removed = build.split_allele_genomes(table, {"GCF_000005845.2", "GCF_000195955.2"}, 0.5, kept)
            self.assertEqual((out, removed), (kept, 1))
            with open(kept) as fh:
                self.assertEqual([line.split("\t")[0] for line in fh], ["GCF_000005845.2", "GCA_900000000.1", "GCF_000195955.2"])
            # Without the representatives known every genome is judged; none to remove: the table itself.
            self.assertEqual(build.split_allele_genomes(table, None, 0.5, kept)[1], 2)
            self.assertEqual(build.split_allele_genomes(table, None, 0.05, kept), (table, 0))

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

    def test_options_checked_before_the_run(self):
        # Since 2026-10-15: a read type given twice (two trainers would write the same files hours in), a --holdout share
        # outside [0, 1) and a --depth-knob-read-types that names no read type stop the run before it starts.
        with tempfile.TemporaryDirectory() as root:
            for bad, said in ((["--read-types", "pe,se,pe"], "each once"), (["--holdout", "1"], "--holdout"),
                              (["--holdout", "-0.1"], "--holdout"), (["--depth-knob-read-types", "pe,xx"],
                                                                     "--depth-knob-read-types")):
                argv = ["build_gtdb_database.py", "--outdir", os.path.join(root, "out"), "--protal", "/nowhere/protal",
                        *bad]
                stderr = io.StringIO()
                with self.assertRaises(SystemExit) as stopped, contextlib.redirect_stderr(stderr), \
                        contextlib.redirect_stdout(io.StringIO()), unittest.mock.patch.object(sys, "argv", argv):
                    build.main()
                self.assertEqual(stopped.exception.code, 2, bad)
                self.assertIn(said, stderr.getvalue(), bad)
                self.assertFalse(os.path.exists(os.path.join(root, "out")), bad)

    def test_an_ignored_signal_stays_ignored(self):
        # Since 2026-10-15: SIGHUP ignored when the run starts (nohup) stays ignored; the others stop the run.
        saved = {sig: build.signal.getsignal(sig) for sig in (build.signal.SIGTERM, build.signal.SIGINT,
                                                               build.signal.SIGHUP)}
        try:
            build.signal.signal(build.signal.SIGHUP, build.signal.SIG_IGN)
            build.signal.signal(build.signal.SIGTERM, build.signal.SIG_DFL)
            build.install_signal_handlers()
            self.assertEqual(build.signal.getsignal(build.signal.SIGHUP), build.signal.SIG_IGN)
            self.assertEqual(build.signal.getsignal(build.signal.SIGTERM), build.on_signal)
        finally:
            for sig, handler in saved.items():
                build.signal.signal(sig, handler)

    def test_the_reports_commands_end_with_the_run(self):
        # Since 2026-10-15: the reports' commands (side_popen, in their thread) end with the jobs when the run stops,
        # with what they started (their process group).
        with tempfile.TemporaryDirectory() as root:
            pid_file = os.path.join(root, "child.pid")
            report = build.side_popen([sys.executable, "-c",
                                       "import subprocess, sys, time; "
                                       f"c = subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(60)']); "
                                       f"open({pid_file!r}, 'w').write(str(c.pid)); time.sleep(60)"])
            deadline = time.time() + 30
            while not os.path.isfile(pid_file) and time.time() < deadline:
                time.sleep(0.05)
            with open(pid_file) as fh:
                child = int(fh.read())
            began = time.time()
            build.stop_jobs()
            self.assertLess(time.time() - began, 30)
            self.assertIsNotNone(report.poll())
            deadline = time.time() + 10
            while time.time() < deadline:
                try:
                    os.kill(child, 0)
                except ProcessLookupError:
                    break
                with contextlib.suppress(ChildProcessError):
                    os.waitpid(child, os.WNOHANG)  # not ours: reaped by init
                time.sleep(0.05)
            else:
                self.fail("the report's own child still runs")
            self.assertEqual(build.SIDE_COMMANDS, [])

    def test_rows_at_other_gates_on_the_disk(self):
        # --test-c-gates: the samples of a simulation profiled again at other levels write a SAM, a profile and dumps each
        # time: KEPT_SHARE of its reads per level of its unit with the most (a paired-end point's pe and se units
        # together), in simulation_sizes' order and selection; none for a collection without levels.
        def command(*extra):
            return [sys.executable, "collect_training_data.py", "--db", "db", "--genome_table", "genomes.tsv", "-o", "out",
                    "--samples", "2", "--long_read_samples", "2", "--read_pairs", "1000", "--read_types", "pe,se,ont",
                    "--read_setups", "100:HS20:300:40", "--long_read_bases", "1000000", *extra]

        units = build.units_of(build.collector_args(command()[2:]))[1]
        ont = next(build.simulation_of(u) for u in units if u["type"] == "ont")
        plain = [("training data", command(), "training.log")]
        gated = [("test set", command("--sparse_first_levels", "pe:0/2,se:0,ont:4", "--stride_levels", "pe:2",
                                      "--flex_indels_levels", "ont:1/4", "--write_simulations", ont), "test.log")]
        sizes = build.simulation_sizes(gated)
        self.assertEqual(len(sizes), 2)
        # pe: 0, 2 and 2s (se only 0); ont: 4, 1f and 4f. Streamed, none (no --model_levels).
        self.assertEqual([c / (build.KEPT_SHARE * r) for _, r, c, _ in build.simulation_costs(gated)], [3, 3])
        self.assertEqual([s for _, _, _, s in build.simulation_costs(gated)], [0, 0])
        self.assertEqual([(largest, r) for largest, r, _, _ in build.simulation_costs(gated)], sizes)
        self.assertEqual([c for _, _, c, _ in build.simulation_costs(plain)], [0, 0])
        # --model_levels (--model-gates' other gates): every sample at them, a streamed simulation's too (pe: 1 and 8,
        # with 0 and 2: four levels written; ont: 8 with 4: two).
        models = [("training data", command("--sparse_first_levels", "pe:0/2,ont:4", "--model_levels", "pe:1/8,ont:8"),
                   "training.log")]
        self.assertEqual([(c / (build.KEPT_SHARE * r), s / (build.KEPT_SHARE * r))
                          for _, r, c, s in build.simulation_costs(models)], [(4, 2), (2, 1)])
        # Written whatever --stream-above says (--write_simulations) or streamable: the ONT simulation and the pe point's.
        written, streamable = build.simulation_costs(gated, "written"), build.simulation_costs(gated, "streamable")
        self.assertEqual(len(written), 1)
        self.assertEqual(len(streamable), 1)
        self.assertEqual(sorted(written + streamable), sorted(build.simulation_costs(gated)))
        self.assertEqual(build.simulation_bytes(units)[ont][1], written[0][1])

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

    def test_a_step_of_several_commands_in_one_log(self):
        # A step's commands one after the other in its log (the foreign scan and --add_tables, the read types' parity
        # checks): each appended one headed by its command line; a file given on the command line that is where the run
        # keeps it already is not copied onto itself.
        with tempfile.TemporaryDirectory() as root:
            log = os.path.join(root, "logs", "step.log")
            build.run([sys.executable, "-c", "print('first')"], log)
            build.run([sys.executable, "-c", "print('second')"], log, append=True)
            with open(log) as fh:
                lines = fh.read().splitlines()
            self.assertEqual(lines[0], "first")
            self.assertTrue(lines[1].startswith("--- ") and lines[1].endswith("print('second')"), lines)
            self.assertEqual(lines[2:], ["second"])
            kept = os.path.join(root, "heldout_species.txt")
            with open(kept, "w") as fh:
                fh.write("s__A a\n")
            build.copy_file(kept, kept)
            build.copy_file(kept, os.path.join(root, "copy.txt"))
            with open(os.path.join(root, "copy.txt")) as fh:
                self.assertEqual(fh.read(), "s__A a\n")

    def test_the_archive_to_share(self):
        # --share-logs: OUTDIR's logs, model_logs/ without the models, the build's metadata, the taxonomy and the tables
        # (shortened), under <name>/ as in OUTDIR; not the database, nor the rest of work/.
        with tempfile.TemporaryDirectory() as root:
            out = os.path.join(root, "r226_v18")
            files = {"console.log": "Ready protal database\n", "logs/convert.log": "converted\n",
                     "model_logs/summary.txt": "F1\n", "model_logs/trained_model.xml": "<PMML/>\n",
                     "model_logs/trained_model.joblib": "x", "model_logs/error_reads/pe/taxa.tsv.gz": "x",
                     "protal_db/build_metadata.tsv": "gtdb_release\tr226\n", "protal_db/database.protal": "x",
                     "work/internal_taxonomy.dmp": "1\t1\troot\n", "work/genomes.tsv": "x",
                     "work/training/training_data.tsv": "truth\tp\n1\t0.12345678901234567\n"}
            for path, text in files.items():
                os.makedirs(os.path.dirname(os.path.join(out, path)), exist_ok=True)
                with open(os.path.join(out, path), "w") as fh:
                    fh.write(text)
            # From the samples' disk (here the collection in work/): the design points' simulator logs and parameters,
            # in the points' order (p2 before p10), a long one cut to its ends; the protal runs' logs, which the build
            # gathers into logs/ only with --scratch; the in-silico strains' table.
            scratch = os.path.join(root, "scratch")
            points = os.path.join(out, "work", "training", "points")
            disk = {os.path.join(points, "p10", "simulate.log"): "p10 simulated\n",
                    os.path.join(points, "p2", "design.log"): "p2 design\n",
                    os.path.join(points, "p2", "sim", "run_params.tsv"): "param\tvalue\n",
                    os.path.join(points, "p3", "stream_pe.log"): "a" * (build.GATHERED_MAX + 10),
                    os.path.join(points, "p2", "sim", "reads", "x.log"): "not a point's log\n",
                    os.path.join(out, "work", "training", "profile_all", "run1", "protal.log"): "Sample s strain alleles\n",
                    os.path.join(scratch, "insilico_strains", "insilico_strains.tsv"): "strain\tcongener\n"}
            for path, text in disk.items():
                os.makedirs(os.path.dirname(path), exist_ok=True)
                with open(path, "w") as fh:
                    fh.write(text)
            with contextlib.redirect_stdout(io.StringIO()):
                build.share_archive(out, [("training", os.path.join(out, "work", "training")), ("test", None)],
                                    insilico_table=os.path.join(scratch, "insilico_strains", "insilico_strains.tsv"))
            with tarfile.open(os.path.join(out, "r226_v18_share.tar.gz")) as tar:
                members = {m.name for m in tar.getmembers() if m.isfile()}
                table = tar.extractfile("r226_v18/work/training/training_data.tsv").read().decode()
                simulations = tar.extractfile("r226_v18/logs/simulations_training.log").read().decode()
                runs = tar.extractfile("r226_v18/logs/protal_runs_training.log").read().decode()
            self.assertEqual(members, {"r226_v18/" + p for p in ("console.log", "logs/convert.log", "model_logs/summary.txt",
                                                                 "model_logs/error_reads/pe/taxa.tsv.gz",
                                                                 "protal_db/build_metadata.tsv",
                                                                 "work/internal_taxonomy.dmp",
                                                                 "work/training/training_data.tsv",
                                                                 "logs/simulations_training.log",
                                                                 "logs/protal_runs_training.log",
                                                                 "model_logs/insilico_strains.tsv")})
            self.assertEqual(table, "truth\tp\n1\t0.123456789\n")
            headers = re.findall(r"==> (.+) <==", simulations)
            self.assertEqual(headers, ["p2/design.log", "p2/sim/run_params.tsv", "p3/stream_pe.log", "p10/simulate.log"])
            self.assertIn("p10 simulated\n", simulations)
            self.assertIn(f"[... {build.GATHERED_MAX + 10 - 2 * build.GATHERED_KEEP} bytes left out ...]", simulations)
            self.assertLess(len(simulations), 2 * build.GATHERED_KEEP + 1000)
            self.assertEqual(runs, "==> run1 <==\nSample s strain alleles\n")
            self.assertFalse(os.path.exists(os.path.join(out, "work", "share_tables")))
            # The build gathered the protal runs' logs (--scratch): not again.
            os.makedirs(os.path.join(out, "logs"), exist_ok=True)
            with open(os.path.join(out, "logs", "protal_runs_training.log"), "w") as fh:
                fh.write("gathered\n")
            with contextlib.redirect_stdout(io.StringIO()):
                build.share_archive(out, [("training", os.path.join(out, "work", "training")), ("test", None)])
            with tarfile.open(os.path.join(out, "r226_v18_share.tar.gz")) as tar:
                names = [m.name for m in tar.getmembers() if m.isfile()]
                self.assertEqual(tar.extractfile("r226_v18/logs/protal_runs_training.log").read().decode(), "gathered\n")
            self.assertEqual(names.count("r226_v18/logs/protal_runs_training.log"), 1)
            self.assertNotIn("r226_v18/model_logs/insilico_strains.tsv", names)

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
        # --insilico-strains auto: the share of the one-genome species given an in-silico strain that keeps the real
        # strains at --min-real-strains of the simulated strains. Two species of two real strains (2/3 real each) and
        # eight one-genome species: real 4/3; without strains of the multi species, 8 * 0.5 * s <= 4/3 gives s = 1/3.
        reps = table({**{f"R{i}": (2, 0) for i in range(2)}, **{f"S{i}": (0, 0) for i in range(8)}})
        share, real, made, single = build.auto_insilico_share(path, reps, 0.0, 0.5)
        self.assertAlmostEqual(share, 1 / 3)
        self.assertAlmostEqual(real, (4 / 3) / 10)
        self.assertEqual((made, single), (0.0, 8))
        # A quarter of the multi species with an in-silico strain: each 1/4 drawn (k = 2), 2/4 real.
        share, real, made, single = build.auto_insilico_share(path, reps, 0.25, 0.5)
        real_draws = 2 * (0.75 * 2 / 3 + 0.25 * 2 / 4)
        made_draws = 2 * 0.25 / 4
        self.assertAlmostEqual(share, (real_draws - made_draws) / 4)
        self.assertAlmostEqual(made * 10, made_draws)
        # A stricter share leaves less room, a looser one never more than all; without real strains none.
        self.assertLess(build.auto_insilico_share(path, reps, 0.0, 0.8)[0], 1 / 3)
        self.assertEqual(build.auto_insilico_share(path, reps, 0.0, 0.01)[0], 1.0)
        reps = table({f"S{i}": (0, 0) for i in range(5)})
        self.assertEqual(build.auto_insilico_share(path, reps, 0.25, 0.5)[0], 0.0)
        # The option takes auto or a share.
        self.assertEqual(build.insilico_share_arg("auto"), "auto")
        self.assertEqual(build.insilico_share_arg("0.25"), 0.25)
        for bad in ("1.5", "some"):
            with self.assertRaises(Exception):
                build.insilico_share_arg(bad)

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

    def test_genome_facts_in_one_pass_and_the_cache(self):
        # One read of a genome gives its length and its contigs' names, as genome_length and trace_relatives.contig_names
        # count them; the cache keeps both by path, size and modification time, and the reports' contig cache gets the
        # names, so that neither a rebuild nor the reports read the genome again.
        import trace_relatives
        root = os.path.join(self.tmp.name, "facts")
        os.makedirs(root)
        plain, zipped = os.path.join(root, "a.fna"), os.path.join(root, "b.fna.gz")
        text = ">a one 123 ACGT\r\nACGTNNacgt\r\n\r\nRYK-*\n>second\nAC GT\n>\nACGT\n"  # 10 + 3 + 4 + 4 letters
        with open(plain, "w", newline="") as fh:
            fh.write(text)
        with gzip.open(zipped, "wt", newline="") as fh:
            fh.write("ACGT\n" + text)  # bases before the first header count too
        self.assertEqual(build.genome_facts(plain), (21, ["a", "second"]))
        self.assertEqual(build.genome_facts(zipped), (25, ["a", "second"]))
        self.assertEqual(build.genome_facts(plain)[1], [n for n in trace_relatives.contig_names(plain) if n])
        cache_path, contigs = os.path.join(root, "cache.tsv.gz"), os.path.join(root, "contigs.tsv.gz")
        rows = [["GA", "d__B;s__A", plain], ["GB", "d__B;s__B", zipped]]
        build.write_genome_table(rows, os.path.join(root, "genomes.tsv"), 2, build.GenomeCache(cache_path), contigs)
        with open(os.path.join(root, "genomes.tsv")) as fh:
            self.assertEqual(fh.read(), f"GA\td__B;s__A\t{plain}\t21\nGB\td__B;s__B\t{zipped}\t25\n")
        self.assertEqual(trace_relatives.genome_contigs([plain, zipped], 1, contigs),
                         {plain: ["a", "second"], zipped: ["a", "second"]})
        # The same file (size and mtime): its facts from the cache, whatever it holds now; a changed one is read again.
        stat = os.stat(plain)
        with open(plain, "w", newline="") as fh:
            fh.write(">x\n" + "G" * (stat.st_size - 4) + "\n")
        os.utime(plain, ns=(stat.st_atime_ns, stat.st_mtime_ns))
        cache = build.GenomeCache(cache_path)
        build.write_genome_table(rows[:1], os.path.join(root, "again.tsv"), 1, cache)
        with open(os.path.join(root, "again.tsv")) as fh:
            self.assertEqual(fh.read(), f"GA\td__B;s__A\t{plain}\t21\n")
        os.utime(plain, ns=(stat.st_atime_ns, stat.st_mtime_ns + 10**9))
        build.write_genome_table(rows[:1], os.path.join(root, "changed.tsv"), 1, cache)
        with open(os.path.join(root, "changed.tsv")) as fh:
            self.assertEqual(fh.read(), f"GA\td__B;s__A\t{plain}\t{stat.st_size - 4}\n")
        # Saved merged: both genomes, the changed one with its new facts.
        self.assertEqual(build.GenomeCache(cache_path).get(plain, build.file_stamps([plain])[plain]),
                         (stat.st_size - 4, ["x"]))
        self.assertEqual(build.GenomeCache(cache_path).get(zipped, build.file_stamps([zipped])[zipped]), (25, ["a", "second"]))

    def test_genome_rows_walk_each_folder_once(self):
        # The genome FASTAs of the release's species, by accession: each folder walked once; an accession with several
        # files takes the first folder's, there the first suffix's (.fna before .fna.gz ...); hidden files left out.
        gtdb = os.path.join(self.tmp.name, "walk")
        reps = os.path.join(gtdb, "genomic_files_reps", "gtdb_genomes_reps_r1", "GCF", "000")
        extra = os.path.join(self.tmp.name, "walk_extra")
        for folder in (reps, extra, os.path.join(extra, ".hidden")):
            os.makedirs(folder)
        with open(os.path.join(gtdb, "bac120_taxonomy_r1.tsv"), "w") as fh:
            for i in range(1, 5):
                fh.write(f"RS_GCF_00000000{i}.1\td__Bacteria;p__P;c__C;o__O;f__F;g__G;s__G s{i}\n")
        for name in ("GCF_000000001.1_genomic.fna.gz", "GCF_000000001.1_genomic.fna", "GCF_000000002.1.fa.gz",
                     "GCF_000000009.1.fna", ".GCF_000000003.1.fna"):
            open(os.path.join(reps, name), "w").close()
        for name in ("GCF_000000002.1.fna", "GCF_000000003.1.fasta", "GCF_000000004.1.txt"):
            open(os.path.join(extra, name), "w").close()
        open(os.path.join(extra, ".hidden", "GCF_000000004.1.fna"), "w").close()
        rows = build.genome_rows(gtdb, "1", [extra])
        self.assertEqual([(r[0], os.path.relpath(r[2], self.tmp.name)) for r in rows], [
            ("GCF_000000001.1", os.path.join("walk", "genomic_files_reps", "gtdb_genomes_reps_r1", "GCF", "000",
                                             "GCF_000000001.1_genomic.fna")),
            ("GCF_000000002.1", os.path.join("walk", "genomic_files_reps", "gtdb_genomes_reps_r1", "GCF", "000",
                                             "GCF_000000002.1.fa.gz")),
            ("GCF_000000003.1", os.path.join("walk_extra", "GCF_000000003.1.fasta"))])
        self.assertEqual(rows[0][1], "d__Bacteria;p__P;c__C;o__O;f__F;g__G;s__G s1")

    def test_the_bases_planned_before_the_in_silico_strains(self):
        # The genome store's room is estimated before the in-silico strains are made (in the background): the table's
        # genomes, and the share given of the one-genome species' genomes again.
        table = os.path.join(self.tmp.name, "planned.tsv")
        with open(table, "w") as fh:
            fh.write("A1\td__B;s__A\t/a1.fna\t100\nA2\td__B;s__A\t/a2.fna\t120\nB1\td__B;s__B\t/b1.fna\t1000\n"
                     "C1\td__B;s__C\t/c1.fna\t50\n")
        self.assertEqual(build.planned_bases(table, 1.0), 100 + 120 + 1000 + 50 + 1000 + 50)
        self.assertEqual(build.planned_bases(table, 0.5), 100 + 120 + 1000 + 50 + 525)
        self.assertEqual(build.planned_bases(table, 0), build.genome_bases(table))

    def test_c_gates_summary(self):
        # --test-c-gates: each read type's test samples at the database's gate (the test file's rows of the samples the
        # levels' tables have) and at its other levels (the trainer's extra_tests), with the alignment's seconds of those
        # samples from the collections' CPU tables (a later row of a sample wins; a sample missing leaves the seconds out),
        # each beside its baseline: the same samples at the database's gate. ont also with --flex_indels at the database's
        # gate 1 ("1f").
        root = self.tmp.name
        test = os.path.join(root, "test")
        scores = {"samples": 2, "taxa": 10, "present": 6, "TP": 5, "FP": 1, "FN": 1, "precision": 5 / 6,
                  "sensitivity": 5 / 6, "F1": 5 / 6}
        ont = dict(scores, F1=0.96)
        metrics = {"pe": {"extra_tests": {"c0": {"scores": dict(scores, FP=0, F1=10 / 11),
                                                 "test_file_same_samples": scores},
                                          "c4": {"scores": dict(scores, FN=3, F1=0.7), "test_file_same_samples": scores}}},
                   "pb": {"extra_tests": {"c4": {"scores": dict(scores, F1=0.9), "test_file_same_samples": scores},
                                          "c4s": {"scores": dict(scores, F1=0.92), "test_file_same_samples": scores},
                                          "c8s": {"scores": dict(scores, F1=0.85), "test_file_same_samples": scores}}},
                   "ont": {"extra_tests": {"c0": {"scores": dict(scores, F1=0.95), "test_file_same_samples": ont},
                                           "c1f": {"scores": dict(scores, F1=0.97), "test_file_same_samples": ont}}}}
        prefixes = {}
        for t, data in metrics.items():
            prefixes[t] = os.path.join(root, f"trained_model_{t}")
            with open(prefixes[t] + ".metrics.json", "w") as fh:
                json.dump(data, fh)
        for c, t, samples in ((0, "pe", ["s1", "s2"]), (4, "pe", ["s1", "s2"]), (0, "ont", ["o1"]), ("1f", "ont", ["o1"]),
                              (4, "pb", ["p1"]), ("4s", "pb", ["p1"]), ("8s", "pb", ["p1"])):
            os.makedirs(os.path.join(test, build.LEVELS_DIR, f"c{c}"), exist_ok=True)
            with open(os.path.join(test, build.LEVELS_DIR, f"c{c}", f"samples_{t}.txt"), "w") as fh:
                fh.writelines(s + "\n" for s in samples)
        cpu = [os.path.join(root, "a.tsv"), os.path.join(root, "b.tsv"), os.path.join(root, "missing.tsv")]
        rows = [("run 1", "aligning s1", 2, 10), ("run 1", "aligning s2", 2, 12), ("run 1", "aligning s1__c0", 3, 15),
                ("run 1", "aligning s2__c0", 3, 15), ("run 1", "aligning s1__c4", 1, 4), ("run 1", "profiling", 9, 9),
                ("run 1", "aligning o1", 5, 50), ("run 1", "aligning o1__c0", 9, 90), ("run 1", "aligning o1__c1f", 6, 60),
                ("run 1", "aligning p1", 2, 20), ("run 1", "aligning p1__c4", 3, 30), ("run 1", "aligning p1__c4s", 3, 31),
                ("run 1", "aligning p1__c8s", 2, 21)]
        for path, part in ((cpu[0], rows[:4]), (cpu[1], rows[4:] + [("run 2", "aligning s1", 1, 8)])):
            with open(path, "w") as fh:
                fh.write("run\tstage\twall_seconds\tcpu_seconds\tthreads\tcores_busy\n")
                fh.writelines(f"{r}\t{s}\t{w}\t{c}\t4\t\n" for r, s, w, c in part)
        out = os.path.join(root, "c_gates.tsv")
        lines = build.c_gates_summary(["pe", "pb", "ont"], prefixes, {"pe": [0, 4], "pb": [4, "4s", "8s"], "ont": [0, "1f"]},
                                      {"pe": 1, "pb": 8, "ont": 1}, test, cpu, out)
        self.assertEqual(lines, ["pe (2 samples): off F1 0.9091, 30 CPU s; 1 (database) F1 0.8333, 20 CPU s; 4 hash F1 "
                                 "0.7000",
                                 "pb (1 sample): 4 hash F1 0.9000, 30 CPU s; 4 stride F1 0.9200 (0.9000 by hash), 31 CPU s; "
                                 "8 hash (database) F1 0.8333, 20 CPU s; 8 stride F1 0.8500 (0.8333 by hash), 21 CPU s",
                                 "ont (1 sample): off F1 0.9500, 90 CPU s; 1 (database) F1 0.9600, 50 CPU s; 1 with "
                                 "--flex_indels F1 0.9700 (0.9600 at the database's gate), 60 CPU s"])
        table = list(csv.DictReader(open(out), delimiter="\t"))
        self.assertEqual(list(table[0]), build.C_GATES_COLUMNS)
        self.assertEqual([(r["read_type"], r["gate"], r["pick"], r["flex_indels"], r["database_gate"], r["model_gate"],
                           r["FP"], r["FN"], r["aligning_cpu_s"], r["aligning_wall_s"], r["baseline_F1"],
                           r["baseline_aligning_cpu_s"]) for r in table],
                         [("pe", "0", "-", "0", "0", "1", "0", "1", "30.0", "6.0", "0.833333", "20.0"),
                          ("pe", "1", "-", "0", "1", "1", "1", "1", "20.0", "3.0", "", ""),
                          ("pe", "4", "hash", "0", "0", "1", "1", "3", "", "", "0.833333", "20.0"),
                          ("pb", "4", "hash", "0", "0", "8", "1", "1", "30.0", "3.0", "0.833333", "20.0"),
                          ("pb", "4", "stride", "0", "0", "8", "1", "1", "31.0", "3.0", "0.833333", "20.0"),
                          ("pb", "8", "hash", "0", "1", "8", "1", "1", "20.0", "2.0", "", ""),
                          ("pb", "8", "stride", "0", "0", "8", "1", "1", "21.0", "2.0", "0.833333", "20.0"),
                          ("ont", "0", "-", "0", "0", "1", "1", "1", "90.0", "9.0", "0.96", "50.0"),
                          ("ont", "1", "-", "0", "1", "1", "1", "1", "50.0", "5.0", "", ""),
                          ("ont", "1", "-", "1", "0", "1", "1", "1", "60.0", "6.0", "0.96", "50.0")])
        self.assertEqual(table[0]["F1"], "0.909091")

    def test_c_gates_summary_with_models_at_other_gates(self):
        # --model-gates: a level is scored by the model of the gate nearest it (scored_by: its own, or the nearest
        # gate's trainer's extra_tests), its baseline by the database gate's model; the hash row beside the database's
        # stride pick at its gate. The present species without a row (unseen) are the trainer's, in the table.
        root = self.tmp.name
        test = os.path.join(root, "test")

        def scores(f1, unseen=0):
            return {"samples": 2, "taxa": 10, "present": 6, "unseen": unseen, "TP": 5, "FP": 1, "FN": 1, "precision": 5 / 6,
                    "sensitivity": 5 / 6, "F1": f1}
        base = scores(0.9)
        own = {"extra_tests": {f"c{c}": {"scores": scores(0.5), "test_file_same_samples": base}
                               for c in (0, 1, 2, 8, "4h")}}
        own["extra_tests"]["c4h"]["scores"] = scores(0.88, unseen=3)
        gate_models = {1: {"extra_tests": {"c0": {"scores": scores(0.91)}, "c1": {"scores": scores(0.92)},
                                           "c2": {"scores": scores(0.93)}}},
                       8: {"extra_tests": {"c8": {"scores": scores(0.85)}}}}
        prefixes = {"pe": os.path.join(root, "trained_model")}
        gate_prefixes = {"pe": {g: os.path.join(root, f"trained_model_c{g}") for g in gate_models}}
        for prefix, data in [(prefixes["pe"], own)] + [(gate_prefixes["pe"][g], m) for g, m in gate_models.items()]:
            with open(prefix + ".metrics.json", "w") as fh:
                json.dump(data, fh)
        levels = [0, 1, 2, 8, "4h"]
        for c in levels:
            os.makedirs(os.path.join(test, build.LEVELS_DIR, f"c{c}"), exist_ok=True)
            with open(os.path.join(test, build.LEVELS_DIR, f"c{c}", "samples_pe.txt"), "w") as fh:
                fh.write("s1\ns2\n")
        scored_by = {"pe": {c: build.nearest_gate(build.level_parts(c)[0], [1, 4, 8]) for c in levels}}
        self.assertEqual(scored_by["pe"], {0: 1, 1: 1, 2: 1, 8: 8, "4h": 4})
        out = os.path.join(root, "c_gates_models.tsv")
        lines = build.c_gates_summary(["pe"], prefixes, {"pe": levels}, {"pe": 4}, test, [], out, gate_prefixes, scored_by,
                                      {"pe": "stride"})
        self.assertEqual(lines, ["pe (2 samples): off (model 1) F1 0.9100; 1 (its model) F1 0.9200; 2 stride (model 1) F1 "
                                 "0.9300; 4 stride (database) F1 0.9000; 4 hash F1 0.8800 (0.9000 by stride); "
                                 "8 stride (its model) F1 0.8500"])
        table = list(csv.DictReader(open(out), delimiter="\t"))
        self.assertEqual([(r["gate"], r["pick"], r["database_gate"], r["model_gate"], r["unseen"], r["F1"], r["baseline_F1"])
                          for r in table],
                         [("0", "-", "0", "1", "0", "0.91", "0.9"), ("1", "-", "0", "1", "0", "0.92", "0.9"),
                          ("2", "stride", "0", "1", "0", "0.93", "0.9"), ("4", "stride", "1", "4", "0", "0.9", ""),
                          ("4", "hash", "0", "4", "3", "0.88", "0.9"), ("8", "stride", "0", "8", "0", "0.85", "0.9")])

    def test_c_gates_rows_over_the_samples_every_gate_profiled(self):
        # Since 2026-10-15: the samples streamed into protal are profiled again only at the models' gates (--model-gates),
        # so those gates' tables have more samples than the others'. Every row is over the samples all the levels without
        # --flex_indels have (compare_samples, which the trainers get as --extra-test-samples), its CPU too; a row with
        # --flex_indels over those of them it has; the database's row over them, from the first level without it.
        root = self.tmp.name
        test = os.path.join(root, "test_compare")
        levels = ["1f", 0, 1]  # 1: a model's gate, with the streamed s3
        for c, samples in (("1f", ["s1"]), (0, ["s1", "s2"]), (1, ["s1", "s2", "s3"])):
            os.makedirs(os.path.join(test, build.LEVELS_DIR, f"c{c}"), exist_ok=True)
            with open(os.path.join(test, build.LEVELS_DIR, f"c{c}", "samples_pe.txt"), "w") as fh:
                fh.writelines(s + "\n" for s in samples)
        self.assertEqual(build.compare_samples(test, "pe", levels), ["s1", "s2"])
        self.assertIsNone(build.compare_samples(test, "se", levels))

        def scores(f1):
            return {"samples": 2, "taxa": 10, "present": 6, "unseen": 0, "TP": 5, "FP": 1, "FN": 1, "F1": f1}
        prefix = os.path.join(root, "trained_model_compare")
        with open(prefix + ".metrics.json", "w") as fh:
            json.dump({"extra_tests": {"c0": {"scores": scores(0.8), "test_file_same_samples": scores(0.9)},
                                       "c1": {"scores": scores(0.7), "test_file_same_samples": scores(0.9)},
                                       "c1f": {"scores": scores(0.6), "test_file_same_samples": scores(0.5)}}}, fh)
        cpu = os.path.join(root, "compare_cpu.tsv")
        with open(cpu, "w") as fh:
            fh.write("run\tstage\twall_seconds\tcpu_seconds\tthreads\tcores_busy\n")
            for sample, c in (("s1", 1), ("s2", 2), ("s3", 100), ("s1__c0", 3), ("s2__c0", 4), ("s1__c1", 5),
                              ("s2__c1", 6), ("s3__c1", 500), ("s1__c1f", 7)):
                fh.write(f"run 1\taligning {sample}\t1\t{c}\t4\t\n")
        out = os.path.join(root, "c_gates_compare.tsv")
        lines = build.c_gates_summary(["pe"], {"pe": prefix}, {"pe": levels}, {"pe": 4}, test, [cpu], out,
                                      compare={"pe": build.compare_samples(test, "pe", levels)})
        self.assertEqual(lines, ["pe (2 samples): off F1 0.8000, 7 CPU s; 1 F1 0.7000, 11 CPU s; 4 hash (database) F1 "
                                 "0.9000, 3 CPU s; 1 with --flex_indels F1 0.6000 (0.5000 at the database's gate), 7 CPU s"])

    def test_cpu_logs(self):
        # The run's CPU use goes to logs only: a row per command that ended (Job.cpu_row: its wall and CPU seconds and the
        # cores it kept busy), and the cores the run kept busy over time (CpuTimeline), with the commands running.
        logs = os.path.join(self.tmp.name, "cpu")
        os.makedirs(logs)
        saved = build.Job.cpu_log
        build.Job.cpu_log = os.path.join(logs, "cpu_jobs.tsv")
        self.addCleanup(setattr, build.Job, "cpu_log", saved)
        timeline = build.CpuTimeline(os.path.join(logs, "cpu_timeline.tsv"), every=0.2)
        timeline.start()
        # Busy until it has used 0.6 s of CPU, however loaded the machine is.
        busy = [sys.executable, "-c", "import time\nwhile time.process_time() < 0.6: pass"]
        job = build.run(busy, os.path.join(logs, "busy.log"), label="keeping a core busy")
        timeline.stop()
        self.assertGreaterEqual(job.cpu, 0.6)
        rows = list(csv.DictReader(open(build.Job.cpu_log), delimiter="\t"))
        self.assertEqual(len(rows), 1)
        self.assertEqual((rows[0]["label"], rows[0]["log"], rows[0]["exit"], rows[0]["priority"]),
                         ("keeping a core busy", "busy.log", "0", "nice 0"))
        self.assertGreater(float(rows[0]["cores_busy"]), 0)
        self.assertLessEqual(float(rows[0]["cores_busy"]), 1.2)  # one thread
        timeline_rows = list(csv.DictReader(open(os.path.join(logs, "cpu_timeline.tsv")), delimiter="\t"))
        self.assertGreaterEqual(len(timeline_rows), 3)
        self.assertTrue(all(r["source"] in ("cgroup", "node") for r in timeline_rows))
        self.assertTrue(any("keeping a core busy" in r["running"] for r in timeline_rows))
        self.assertTrue(any(float(r["cores_busy"]) > 0 for r in timeline_rows))

    def test_a_job_that_need_not_stop_the_run(self):
        # Job(fatal=False): a check once the database is finished fails without stopping the build; its exit code kept.
        job = build.Job([sys.executable, "-c", "raise SystemExit(3)"], os.path.join(self.tmp.name, "check.log"),
                        fatal=False).finish()
        self.assertEqual(job.rc, 3)
        self.assertNotIn(job, build.Job.running)


class F1AuditTest(unittest.TestCase):
    """What build_gtdb_database.py does since the f1 audit (docs/claude/2026-10-15-f1-audit-since-0.6.0a): a test seed of
    each build's own (E2), the feature code's versions for the trainers (C1), the marker truth's F1 first in the summary
    (E3), and the scoring of the final check on the finished database (E1) and of the run without the allele scores (B1),
    database_check.py, on tables of a few rows."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def test_a_test_seed_per_build(self):
        # Drawn from the release, OUTDIR's name and the date; kept by a rerun into the OUTDIR; --test-seed replaces it.
        a = build.derive_test_seed("226", "/x/r226_v29", "2026-10-15", 1)
        self.assertEqual(a, build.derive_test_seed("226", "/y/r226_v29/", "2026-10-15", 1))
        self.assertNotEqual(a, build.derive_test_seed("226", "/x/r226_v30", "2026-10-15", 1))
        self.assertNotEqual(a, build.derive_test_seed("226", "/x/r226_v29", "2026-10-16", 1))
        for seed in (1, a, a + 5, 999_999):  # in range, and never within TEST_SEED_GAP of --seed
            value = build.derive_test_seed("226", "/x/r226_v29", "2026-10-15", seed)
            self.assertTrue(build.TEST_SEED_RANGE[0] <= value < build.TEST_SEED_RANGE[1], value)
            self.assertGreaterEqual(abs(value - seed), build.TEST_SEED_GAP)
        stages = build.Stages(os.path.join(self.tmp.name, "stages"))
        first, how = build.test_seed(stages, None, "226", "/x/r226_v29", 1)
        self.assertTrue(how.startswith("drawn for this build from r226, r226_v29 and "), how)
        again, how = build.test_seed(stages, None, "226", "/x/r226_v29", 1)
        self.assertEqual(again, first)
        self.assertTrue(how.startswith("kept from an earlier run into this OUTDIR (drawn"), how)
        given, how = build.test_seed(stages, 1001, "226", "/x/r226_v29", 1)
        self.assertEqual(given, 1001)
        self.assertIn(f"this OUTDIR's runs before had {first}", how)
        self.assertEqual(build.test_seed(stages, None, "226", "/x/r226_v29", 1)[0], 1001)

    def test_feature_versions_for_the_trainers(self):
        # protal --feature_versions saved and passed on (--feature-versions) where the trainer takes it; an older protal
        # or trainer is noted and skipped.
        new_protal, old_protal = (os.path.join(self.tmp.name, n) for n in ("protal_new", "protal_old"))
        with open(new_protal, "w") as fh:
            fh.write('#!/bin/sh\nif [ "$1" = "--feature_versions" ]; then\n'
                     "printf 'feature\\tversion\\tchanged\\tnote\\nfragments\\t2\\t2026-10-15\\tcounted from the "
                     "alignment\\nidentity\\t1\\t2026-08-10\\t\\n'\nexit 0\nfi\nexit 1\n")
        with open(old_protal, "w") as fh:
            fh.write("#!/bin/sh\necho \"Option 'feature_versions' does not exist\" >&2\nexit 1\n")
        for path in (new_protal, old_protal):
            os.chmod(path, 0o755)
        trainers = {}
        for name, option in (("new", "--feature-versions FILE"), ("old", "")):
            trainers[name] = os.path.join(self.tmp.name, f"trainer_{name}.py")
            with open(trainers[name], "w") as fh:
                fh.write(f"print('usage: trainer --truth-file T {option}')\n")
        path = os.path.join(self.tmp.name, "feature_versions.tsv")
        options, said = build.feature_versions(new_protal, path, trainers["new"])
        self.assertEqual(options, ["--feature-versions", path])
        self.assertIn("of 2 features", said)
        with open(path) as fh:
            self.assertEqual(fh.readline(), "feature\tversion\tchanged\tnote\n")
        options, said = build.feature_versions(new_protal, path, trainers["old"])
        self.assertEqual(options, [])
        self.assertIn("this trainer takes no --feature-versions", said)
        options, said = build.feature_versions(old_protal, os.path.join(self.tmp.name, "none.tsv"), trainers["new"])
        self.assertEqual(options, [])
        self.assertIn("has no --feature_versions", said)
        self.assertFalse(os.path.exists(os.path.join(self.tmp.name, "none.tsv")))

    def test_the_headline_is_the_marker_truth(self):
        # summary.txt's first lines: the trainer's "headline" (its held_out and test estimates' marker_truth_f1, with its
        # measure and notes) if its metrics have it, else the marker truth's F1 markers of the test set and with species
        # held out, else the F1 over the rows.
        rows = {"marker_truth": {"rows": [
            {"set": "species held out", "scenario": "(all)", "F1 markers": 0.91},
            {"set": "independent test set", "scenario": "(design)", "F1 markers": 0.88},
            {"set": "hold-out", "scenario": "gut", "F1 markers": 0.5}]},
            "evaluation": {"species": {"F1": 0.95}}, "test": {"this one": {"F1": 0.94}}}
        self.assertEqual(build.marker_headline(rows),
                         ("marker truth", [("independent test set", 0.88), ("held out", 0.91)], ""))
        head = {"knob": "knob curve", "measure": "marker truth", "scheme": "species and sample", "notes": ["a", "b"],
                "held_out": {"marker_truth_f1": 0.89, "rows_f1": 0.93}, "held_out_species": {"marker_truth_f1": 0.9},
                "test": {"marker_truth_f1": 0.87, "rows_f1": 0.92}}
        self.assertEqual(build.marker_headline({**rows, "headline": head}),
                         ("marker truth", [("independent test set", 0.87), ("held out", 0.89)], "a; b"))
        rows_only = {**head, "measure": "rows", "notes": "no gene positions",
                     "test": {"marker_truth_f1": None, "rows_f1": 0.92}}
        self.assertEqual(build.marker_headline({"headline": rows_only}),
                         ("rows", [("independent test set", 0.92), ("held out", 0.93)], "no gene positions"))
        self.assertEqual(build.marker_headline({k: v for k, v in rows.items() if k != "marker_truth"}),
                         ("rows", [("independent test set", 0.94), ("held out", 0.95)],
                          "no marker truth: collected without gene positions?"))
        self.assertEqual(build.marker_headline({})[1], [])
        prefixes = {t: os.path.join(self.tmp.name, f"model_{t}") for t in ("pe", "se")}
        with open(prefixes["pe"] + ".metrics.json", "w") as fh:
            json.dump({"headline": head, **rows}, fh)
        lines = build.headline_lines(["pe", "se"], prefixes)
        self.assertTrue(lines[0].startswith("Headline F1, over every species present"))
        self.assertEqual(lines[1:3], ["  pe: independent test set 0.8700, held out 0.8900 (marker truth; a; b)",
                                      "  se: - (no metrics)"])

    def write(self, name, header, rows, gz=False):
        path = os.path.join(self.tmp.name, name)
        with (gzip.open(path, "wt") if gz else open(path, "w")) as fh:
            fh.write("\t".join(header) + "\n")
            fh.writelines("\t".join(str(v) for v in row) + "\n" for row in rows)
        return path

    def calls_file(self, name, rows):
        """A PREFIX.calls.tsv.gz of rows (sample, taxon, truth, call, novel congener, marker reads, set)."""
        return self.write(name, ["meta_sample", "meta_scenario", "meta_novel_congener", "meta_marker_reads", "taxon",
                                 "taxon_name", "truth", "set", "p", "knob", "call"],
                          [(s, "", nc, m, 0, t, truth, st, 0.5, 0.5, call) for s, t, truth, call, nc, m, st in rows], gz=True)

    def table_file(self, name, rows, extra=()):
        """A collector table of rows (sample, taxon, truth, prediction, novel congener, marker reads, *extra values)."""
        return self.write(name, ["meta_sample", "meta_scenario", "meta_novel_congener", "meta_marker_reads", "truth",
                                 "prediction", "probability", "taxon", "taxon_name", *extra],
                          [(s, "", nc, m, truth, pred, 0.5, 0, t, *more) for s, t, truth, pred, nc, m, *more in rows])

    def test_the_final_check(self):
        import database_check
        heldout = os.path.join(self.tmp.name, "heldout_species.txt")
        with open(heldout, "w") as fh:
            fh.write("s__H\tspecies\ts__H\n")
        # The trainer's calls on the training database: s1 A found, B a false positive, s2 A missed; another sample and
        # the training rows are not compared.
        calls = self.calls_file("model.calls.tsv.gz", [
            ("s1", "s__A", 1, 1, 0, 5, "test"), ("s1", "s__B", 0, 1, 1, "", "test"), ("s2", "s__A", 1, 0, 0, 5, "test"),
            ("s2", "s__B", 0, 0, 0, "", "test"), ("s3", "s__A", 1, 0, 0, 5, "test"), ("t1", "s__A", 1, 0, 0, 5, "training")])
        # On the finished database: both A found, B rejected, the held-out H found (left out), D called (a row of its own).
        table = self.table_file("final.tsv", [
            ("s1", "s__A", 1, 1, 0, 5), ("s1", "s__B", 0, 0, 1, ""), ("s2", "s__A", 1, 1, 0, 5), ("s2", "s__B", 0, 0, 0, ""),
            ("s1", "s__H", 1, 1, 0, 9), ("s2", "s__D", 0, 1, 0, "")])
        # The training database missed C in s1 without a row (unseen), the finished database none.
        unseen = self.write("test.unseen.tsv", ["meta_design", "meta_sample", "meta_scenario", "taxon_name", "marker_reads"],
                            [("d", "s1", "", "s__C", 4), ("d", "s3", "", "s__C", 4), ("d", "s1", "", "s__H", 4)])
        out = os.path.join(self.tmp.name, "final_db_check.tsv")
        rows, missing = database_check.final_check(
            ["pe", "se"], {"pe": table}, {"pe": (unseen, None)}, {("pe", "forest"): calls}, {}, "forest",
            database_check.read_species(heldout), out)
        self.assertEqual(missing, ["se"])
        main, congener = rows
        self.assertEqual((main["rows_of"], main["calls"], main["samples"]), ("species not held out", "protal", 2))
        self.assertEqual([main[f"training_db_{c}"] for c in ("rows", "TP", "FP", "FN")], [4, 1, 1, 1])
        self.assertEqual([main[f"finished_db_{c}"] for c in ("rows", "TP", "FP", "FN")], [5, 2, 1, 0])
        self.assertAlmostEqual(main["training_db_F1"], 0.5)
        self.assertAlmostEqual(main["finished_db_F1"], 0.8)
        self.assertAlmostEqual(main["F1_change"], 0.3)
        self.assertEqual((main["training_db_marker_FN"], main["finished_db_marker_FN"]), (2, 0))  # C unseen in s1
        self.assertAlmostEqual(main["training_db_marker_F1"], 0.4)
        self.assertEqual([main[f] for f in database_check.FLIPS], [1, 1, 0, 1])
        # Without a held-out congener in the sample: s1's B left out on both.
        self.assertEqual((congener["training_db_rows"], congener["training_db_FP"], congener["finished_db_rows"]), (3, 0, 4))
        with open(out) as fh:
            self.assertEqual(fh.readline().rstrip("\n").split("\t"), database_check.FINAL_COLUMNS)
        kept = database_check.read_rows(out)
        self.assertEqual(kept[0]["training_db_TP"], 1)
        lines, warnings = database_check.final_check_lines(kept, missing, "forest")
        self.assertTrue(lines[0].startswith("pe forest (default), species not held out: F1 0.5000 -> 0.8000 (+0.3000;"),
                        lines[0])
        self.assertIn("no samples written to the disk of se", lines[-1])
        self.assertEqual(len(warnings), 1)  # the subset without congeners is no warning of its own
        self.assertIn("beyond a refit's noise of 0.002", warnings[0])

    def test_model_calls_and_the_allele_scores(self):
        import numpy as np
        import database_check
        from sklearn.ensemble import RandomForestClassifier
        # A forest of one split: present above identity 0.95, written as protal reads it, with a knob of 0.5.
        X = np.array([[10, 0.9], [10, 0.92], [10, 0.97], [10, 0.99]])
        forest = RandomForestClassifier(n_estimators=1, max_depth=1, bootstrap=False, random_state=0).fit(X, [0, 0, 1, 1])
        model = os.path.join(self.tmp.name, "model.xml")
        model_pmml.write_forest(forest, ["fragments", "identity"], model)
        # Normal run (the trainer's calls): A called in both samples, B called in s1 (a false positive), C missed.
        calls = self.calls_file("model.calls.tsv.gz", [
            ("s1", "s__A", 1, 1, 0, 5, "test"), ("s1", "s__B", 0, 1, 0, "", "test"), ("s2", "s__A", 1, 1, 0, 5, "test"),
            ("s2", "s__C", 1, 0, 0, 5, "test")])
        # Without the allele scores: A drops below 0.95 in s2 (missed), B too (rejected), C as before; A has alleles.
        table = self.table_file("no_alleles.tsv", [
            ("s1", "s__A", 1, 0, 0, 5, 10, 0.99, 0.6), ("s1", "s__B", 0, 0, 0, "", 10, 0.90, 0),
            ("s2", "s__A", 1, 0, 0, 5, 10, 0.93, 0.6), ("s2", "s__C", 1, 0, 0, 5, 10, 0.91, 0)],
            extra=("fragments", "identity", "allele_copy_share"))
        frame = database_check._frame(table)
        self.assertEqual(database_check.model_calls(frame, model).tolist(), [1, 0, 0, 0])
        out = os.path.join(self.tmp.name, "no_allele_scores.tsv")
        rows, missing = database_check.allele_check(["pe"], {"pe": table}, {"pe": None}, {("pe", "forest"): calls},
                                                    {("pe", "forest"): model}, out)
        self.assertEqual(missing, [])
        every, with_alleles, without = rows
        self.assertEqual((every["species"], with_alleles["species"], without["species"]),
                         ("all", "with strain alleles", "without strain alleles"))
        self.assertEqual([every[f"with_allele_scores_{c}"] for c in ("TP", "FP", "FN")], [2, 1, 1])
        self.assertEqual([every[f"without_allele_scores_{c}"] for c in ("TP", "FP", "FN")], [1, 0, 2])
        self.assertAlmostEqual(every["allele_scores_F1_gain"], 2 * 2 / 6 - 2 / 4)
        self.assertEqual((with_alleles["taxa"], with_alleles["with_allele_scores_TP"], with_alleles["without_allele_scores_FN"]),
                         (1, 2, 1))
        self.assertEqual((without["taxa"], without["with_allele_scores_FP"], without["without_allele_scores_FP"]), (2, 1, 0))
        lines = database_check.allele_check_lines(rows, missing)
        self.assertTrue(lines[0].startswith("pe forest, all (3 taxa): F1 0.5000 without the allele scores -> 0.6667 with them"),
                        lines[0])


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

    def test_holdout_keeps_species_complexes_whole(self):
        # With the species' clouds (protal --write_species_neighbours: each species' nearest congeners), congeners within
        # --holdout-complex-distance form a complex that is held out or kept whole; the draw is the old one without them.
        import collect_training_data as collect
        gtdb, db = shared_release()
        taxonomy = os.path.join(db, "internal_taxonomy.dmp")
        table = os.path.join(gtdb, "simulation", "genomes.tsv")
        _, ids = lineages.from_taxonomy(taxonomy) if False else (None, None)
        with tempfile.TemporaryDirectory() as tmp:
            # species_neighbours.tsv's format, by taxid: Mockella alpha and beta 0.08 apart, Fakibacter gamma alone.
            by_id, _ = __import__("lineages").from_taxonomy(taxonomy)
            taxid = {lin["species"]: tid for tid, lin in by_id.items() if "species" in lin}
            path = os.path.join(tmp, "species_clouds.tsv")
            with open(path, "w") as fh:
                fh.write("taxid\tneighbours\n")
                fh.write(f"{taxid['s__Mockella alpha']}\t{taxid['s__Mockella beta']}:0.08\n")
                fh.write(f"{taxid['s__Mockella beta']}\t{taxid['s__Mockella alpha']}:0.08\n")
                fh.write(f"{taxid['s__Fakibacter gamma']}\t\n")
            clouds = collect.read_clouds(path, taxonomy)
            self.assertEqual(clouds, {"s__Mockella alpha": [("s__Mockella beta", 0.08)],
                                      "s__Mockella beta": [("s__Mockella alpha", 0.08)], "s__Fakibacter gamma": []})
            self.assertEqual(build.species_complexes(clouds, 0.01), {})
            self.assertEqual(build.species_complexes(clouds, 0.1), {"s__Mockella alpha": 1, "s__Mockella beta": 1})
            # No complex within 0.01: the draw of before (one species of three at 0.34).
            plain = build.choose_holdout(table, taxonomy, 0.34, {}, 1.0, 3)
            self.assertEqual(plain, build.choose_holdout(table, taxonomy, 0.34, {}, 1.0, 3, clouds, 0.01))
            self.assertEqual(len(plain), 1)
            # Within 0.1, alpha and beta go together, whichever seed: both or neither.
            for seed in range(6):
                chosen = build.choose_holdout(table, taxonomy, 0.34, {}, 1.0, seed, clouds, 0.1)
                self.assertEqual(chosen, build.choose_holdout(table, taxonomy, 0.34, {}, 1.0, seed, clouds, 0.1))
                mockella = {"s__Mockella alpha", "s__Mockella beta"} & set(chosen)
                self.assertIn(len(mockella), (0, 2), chosen)
                self.assertTrue(chosen, chosen)
                # heldout_species.txt with the distance to the nearest kept congener and the complex; read back
                path = os.path.join(tmp, f"heldout{seed}.txt")
                build.write_holdout(path, chosen, clouds, 0.1)
                self.assertEqual(build.read_holdout(path), chosen)
                self.assertEqual(collect.read_novel(path), chosen)
                details = build.read_holdout_details(path)
                self.assertEqual(set(details), set(chosen))
                for s in chosen:
                    self.assertEqual(details[s], (None, "c1" if s in mockella else ""), s)  # the kept congener: none
                lines = build.describe_holdout(chosen, build.pool_species(table), details, 0.1)
                self.assertTrue(any("nearest kept congener" in line and "farther or none" in line for line in lines), lines)
                self.assertTrue(any(f"held out whole: {1 if mockella else 0}," in line for line in lines), lines)
            # One of the complex held out by hand: its twin kept, the distance to it recorded.
            torn = {"s__Mockella alpha": ("species", "s__Mockella alpha")}
            path = os.path.join(tmp, "torn.txt")
            build.write_holdout(path, torn, clouds, 0.1)
            self.assertEqual(build.read_holdout_details(path), {"s__Mockella alpha": (0.08, "c1")})
            self.assertEqual(build.cloud_band(0.08), "0.05-0.15")
            self.assertEqual(build.cloud_band(None), "farther or none")
            self.assertEqual(build.cloud_band(0.2), "farther or none")


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
        """[(gene id, header, sequence)] of a FASTA's bytes, one sequence line per record (the header's name taxid_geneid,
        a full reference's followed by the genome's accession)."""
        lines = data.split(b"\n")
        return [(int(lines[i][1:].split()[0].split(b"_")[1]), lines[i], lines[i + 1]) for i in range(0, len(lines) - 1, 2)]

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
        same_files(self, self.db, self.copy, "internal_taxonomy.dmp", "species_priors.tsv", "ncbi_taxids.tsv", "genome2tiid.tsv")
        # The genomes with a gene kept: the archaea have none of the bacterial genes chosen here.
        self.assertRegex(self.derived, r"gene_neighbours\.tsv derived from the 1[2-6] genomes of [6-8] species kept: \d+ lines "
                                       r"\(over the 6 genes kept\)")
        self.assertRegex(self.derived, r"without 0 species, 6 genes kept: \d+ representative sequences kept")

    def test_a_direct_conversion_gives_the_same_files(self):
        direct = os.path.join(self.tmp.name, "direct")
        run(CONVERT, "--gtdb", self.gtdb, "--outdir", direct, "--genes", self.spec)
        same_files(self, direct, self.copy, "reference.fna", "reference.map", "gene2geneid.tsv", "internal_taxonomy.dmp",
                   "species_priors.tsv", "ncbi_taxids.tsv")
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

        def build_db(folder, *extra):  # zstd level 1 without long-distance matching, 4 threads: the same content
            return subprocess.run([protal, "--build", "--no_profile", "-t", "4", "--no_bundle", "--compress_level", "1",
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
os.makedirs(os.path.join(out, "work"), exist_ok=True)
with open(os.path.join(out, "work", "internal_taxonomy.dmp"), "w") as fh:
    fh.write("taxid\tparent\tlevel\tname\trank\n" + "".join(f"{t}\t{p}\t-\t{name}\t{rank}\n" for t, p, name, rank in nodes))
with open(os.path.join(out, "model_logs", "summary.txt"), "w") as fh:
    fh.write("read type  evaluated on           knob  taxa  TP  FP  TN  FN  sensitivity  specificity  precision  F1      FP rate  FN rate  FP/sample\n"
             "pe         species held out       0.5   100   9   1   88  2   0.8182       0.9888       0.9000     0.8571  0.0112   0.1818   0.50\n"
             "pe         independent test set   0.5   50    5   0   44  1   0.8333       1.0000       1.0000     0.9091  0.0000   0.1667   0.00\n"
             "pe         soil: hold-out             0.5   40    3   2   33  2   0.6000       0.9429       0.6000     0.6000  0.0571   0.4000   1.00\n"
             "pe         marker truth: species held out  0.5   100   9   1   88  3   0.7500       0.9888       0.9000     0.8182  0.0112   0.2500   0.50\n"
             "pe         marker truth: independent test set  0.5   50    5   0   44  2   0.7143       1.0000       1.0000     0.8333  0.0000   0.2857   0.00\n"
             "pe         marker truth: soil hold-out  0.5   40    3   2   33  3   0.5000       0.9429       0.6000     0.5455  0.0571   0.5000   1.00\n")
if "--rank-genes" in args:
    with open(os.path.join(out, "model_logs", "gene_ranking.tsv"), "w") as fh:
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
                         ("3", os.path.join(self.out, "r226_full", "model_logs", "gene_ranking.tsv")))
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
            # (the model's own rows, not the scenario's or the marker truth's that follow them)
            self.assertEqual((r["pe_test_F1"], r["pe_test_FP_per_sample"], r["pe_heldout_F1"], r["se_test_F1"]),
                             ("0.9091", "0.00", "0.8571", ""))
            self.assertEqual((r["pe_test_marker_F1"], r["pe_heldout_marker_F1"]), ("0.8333", "0.8182"))
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
