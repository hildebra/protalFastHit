#!/usr/bin/env python3
"""test_collector.py - checks for collect_training_data.py (the training data's designs, the long reads' replay of the
paired-end communities, the simulations' scheduler, --follow, the worker processes), scenarios.py and hifi_reads.py.

A stand-in replaces protal; the long reads are made by simulate_metagenomes ($SIMULATE: those tests are skipped without
it, or fail under PROTAL_TESTS_REQUIRED=1), art_illumina is needed only for the host genome's paired-end reads (that
part is skipped without it). Needs numpy.

  python3 -m unittest scripts/mini_db/test_collector.py
"""

import argparse
import collections
import contextlib
import glob
import gzip
import io
import json
import math
import os
import random
import shutil
import socket
import subprocess
import sys
import tempfile
import threading
import time
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mini_db_fixtures  # noqa: E402,F401  (the scripts' folders on sys.path)
import prerequisites  # noqa: E402

import collect_training_data as collect  # noqa: E402
import compressed  # noqa: E402
import lineages  # noqa: E402
import scenarios  # noqa: E402


SIMULATOR = os.environ.get("SIMULATE", "")  # simulate_metagenomes, which makes the long reads


def need_simulator():
    """$SIMULATE (simulate_metagenomes): the test is skipped without it, or fails under PROTAL_TESTS_REQUIRED."""
    if not prerequisites.executable(SIMULATOR):
        prerequisites.missing("needs $SIMULATE (simulate_metagenomes) for the long reads")
    return SIMULATOR


def perfect_qshmm(path):
    """A pbsim3 qshmm model (QSHMM-*.model's format) whose reads have Q93 throughout, so no errors: one state per
    accuracy level 71-99, emitting Q93."""
    with open(path, "w") as fh:
        for level in range(71, 100):
            fh.write(f"{level} IP 1 1.0\n{level} EP 1 " + " ".join("1.0" if q == 93 else "0" for q in range(94)) +
                     f"\n{level} TP 1 1.0\n")


class Gate:
    """Jobs that record when they start and end, and wait at their start until the test opens the gate: what runs at
    once is then the scheduler's doing, not the threads' timing."""

    def __init__(self):
        self.cond = threading.Condition()
        self.open = False
        self.events, self.busy, self.most = [], collections.Counter(), collections.Counter()

    def job(self, name, group=None, wait=True, error=None, new=()):
        def run():
            with self.cond:
                self.events.append(("start", name))
                self.busy[group] += 1
                self.most[group] = max(self.most[group], self.busy[group])
                self.cond.notify_all()
                if wait:
                    self.cond.wait_for(lambda: self.open, timeout=60)
                self.busy[group] -= 1
                self.events.append(("end", name))
            return error, list(new)
        return run

    def started(self):
        return [name for kind, name in self.events if kind == "start"]

    def wait_until(self, condition):
        with self.cond:
            return self.cond.wait_for(condition, timeout=60)

    def release(self):
        with self.cond:
            self.open = True
            self.cond.notify_all()


class CollectorTest(unittest.TestCase):
    """collect_training_data.py's parts on their own."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def test_collector_designs(self):
        # Read setups (built-in and custom ART profiles), abundance models, long-read setups and units.
        opts = argparse.Namespace(read_setups="150:HSXt:350:50,150:file=/p1+/p2:350:50,250:MSv3:550:50",
                                  read_pairs="1000,5000:2", read_types=["pe", "se", "ont"], samples=4, long_read_samples=0,
                                  pb_setup="errhmm:ERRHMM-SEQUEL:15000:3000:0.999",
                                  ont_setup="qshmm:QSHMM-ONT-HQ:8000:6000:0.97", long_read_bases="1e6,2e6,3e6:5")
        points = collect.design_points(opts)
        self.assertEqual([p["name"] for p in points], ["rl150_HSXt_p1000", "rl150_HSXt_p5000", "rl150_custom1_p1000",
                                                       "rl150_custom1_p5000", "rl250_p1000", "rl250_p5000"])
        # DEPTH:SAMPLES gives a depth's points other samples than --samples.
        self.assertEqual([p["samples"] for p in points], [4, 2, 4, 2, 4, 2])
        # Setups of one length and profile are told apart by their fragments; a setup or depth given twice
        # would share a folder.
        same = argparse.Namespace(read_setups="150:HS25:350:50,150:HS25:500:80", read_pairs="1000", samples=1)
        self.assertEqual([p["name"] for p in collect.design_points(same)], ["rl150_HS25_f350-50_p1000", "rl150_HS25_f500-80_p1000"])
        for setups, pairs in (("150:HS25:350:50,150:HS25:350:50", "1000"), ("150:HS25:350:50", "1000,1000"),
                              ("150:HS25:350:50", "1000,5000:0"), ("150:HS25:350:50", "1000:x")):
            with self.assertRaises(SystemExit):
                collect.design_points(argparse.Namespace(read_setups=setups, read_pairs=pairs, samples=1))
        with self.assertRaises(SystemExit):
            collect.units_of(argparse.Namespace(**{**vars(opts), "long_read_bases": "1e6,1e6"}))
        with self.assertRaises(SystemExit):
            collect.art_profile_args("file=/nonexistent_r1.txt")
        self.assertEqual(collect.art_profile_args("HSXt"), ["--sequencer", "HSXt"])
        self.assertEqual(collect.abundance_args("lognormal:2.0"), ["--distribution", "poisson_lognormal", "--pln_sigma", "2.0"])
        self.assertEqual(collect.abundance_args("powerlaw:1.5"), ["--distribution", "power_law", "--alpha", "1.5"])
        self.assertEqual(collect.abundance_args(""), [])
        with self.assertRaises(SystemExit):
            collect.abundance_args("gamma:1")
        self.assertEqual(collect.parse_long_setup("qshmm:QSHMM-ONT-HQ:8000:6000:0.97")["length_mean"], 8000)
        hifi = collect.parse_long_setup("hifi:15000:3000:3")
        self.assertEqual((hifi["method"], hifi["length_mean"], hifi["length_sd"], hifi["q_sd"]), ("hifi", 15000, 3000, 3.0))
        with self.assertRaises(SystemExit):
            collect.parse_long_setup("art:x:1:1:1")
        _, units = collect.units_of(opts)
        self.assertEqual([u["type"] for u in units], ["pe"] * 6 + ["se"] * 6 + ["ont"] * 3)
        # A long-read point replays the communities of the paired-end points of the depth at its place (modulo the
        # depths), of every read setup: more samples than one point has (3:5 against 4 per point).
        ont = [u for u in units if u["type"] == "ont"]
        self.assertEqual([[p["name"] for p in u["communities"]] for u in ont],
                         [["rl150_HSXt_p1000", "rl150_custom1_p1000", "rl250_p1000"],
                          ["rl150_HSXt_p5000", "rl150_custom1_p5000", "rl250_p5000"],
                          ["rl150_HSXt_p1000", "rl150_custom1_p1000", "rl250_p1000"]])
        self.assertEqual([u["samples"] for u in ont], [4, 4, 5])
        self.assertEqual([u["samples"] for u in units if u["type"] == "se"], [4, 2, 4, 2, 4, 2])
        self.assertEqual(units[6]["name"], "rl150_HSXt_p1000_se")
        # No more long-read samples than the communities they replay (here 3 x 2 at 5,000 read pairs).
        with self.assertRaises(SystemExit):
            collect.units_of(argparse.Namespace(**{**vars(opts), "long_read_bases": "1e6,2e6:7"}))
        _, units = collect.units_of(argparse.Namespace(**{**vars(opts), "long_read_samples": 12}))
        self.assertEqual([u["samples"] for u in units if u["type"] == "ont"], [12, 6, 5])

    def test_long_read_replay(self):
        # pb/ont samples replay a paired-end point's communities: one simulate_metagenomes run per unit draws the
        # reads (a genome by relative abundance times length, a start uniform, either strand) until each sample's bases,
        # and makes one read of each (here by a qshmm model of Q93 throughout, so without errors).
        simulator = need_simulator()
        root = os.path.join(self.tmp.name, "longreads")
        point = os.path.join(root, "points", "rl150_p1000")
        os.makedirs(os.path.join(point, "sim"))
        genomes, sequences = {}, {}
        rng = random.Random(5)
        for name, length in (("GA", 20000), ("GB", 40000)):
            genomes[name] = os.path.join(root, name + ".fna.gz")
            sequences[name] = "".join(rng.choice("ACGT") for _ in range(length))
            with gzip.open(genomes[name], "wt") as fh:
                fh.write(f">{name}_contig\n" + "\n".join(sequences[name][i:i + 70] for i in range(0, length, 70)) + "\n")
        with open(os.path.join(point, "sim", "manifest.tsv"), "w") as fh:
            fh.write("sample\tgenome\tspecies\ttaxonomy\tgenome_length\tread_pairs\tvertical_coverage\trelative_abundance"
                     "\tfastq_r1\tfastq_r2\tfasta_path\tart_seed\n")
            for sample, abundances in (("rl150_p1000_s_1", (0.5, 0.5)), ("rl150_p1000_s_2", (0.8, 0.2))):
                for (name, path), a in zip(genomes.items(), abundances):
                    length = 20000 if name == "GA" else 40000
                    fh.write(f"{sample}\t{name}\tS {name}\td__B;s__S {name}\t{length}\t10\t1\t{a}\tr1\tr2\t{path}\t1\n")
        with open(os.path.join(point, "sim", "protal.meta"), "w") as fh:
            fh.write(f"#OUTPUT_DIR\t{point}/protal\n#INPUT_DIR\t{point}/sim/reads\n"
                     "#SAMPLEID\tFIRST\tSECOND\tSAM\tPREFIX\tPROFILE\tPROFILE_TRUTH\n"
                     "rl150_p1000_s_1\ta\tb\tc\td\te\t/truth/1\nrl150_p1000_s_2\ta\tb\tc\td\te\t/truth/2\n")
        models = os.path.join(root, "models")
        os.makedirs(models)
        perfect_qshmm(os.path.join(models, "PERFECT.model"))
        opts = argparse.Namespace(out=root, seed=1, simulator=simulator, pbsim="none", pbsim_models=models, samples=2)
        unit = {"type": "ont", "name": "ont_b300000", "bases": 300000,
                "setup": collect.parse_long_setup("qshmm:PERFECT:1000:0:0.97"), "samples": 2,
                "communities": [{"name": "rl150_p1000"}],
                "point": {"name": "ont_b300000", "read_length": "1000", "read_pairs": "300000"}}
        keys = {"ont_b300000": {"a": 1}}
        self.assertEqual(collect.simulate_long([(0, unit)], opts, 2, keys), {})
        self.assertTrue(collect.same_key(os.path.join(root, "points", "ont_b300000", "simulated.json"), keys["ont_b300000"]))
        self.assertFalse(os.path.exists(os.path.join(root, "points", "ont_b300000", "sim", "tmp")))
        rows, _ = collect.unit_map_rows(unit, opts)
        self.assertEqual([r["SAMPLEID"] for r in rows], ["ont_b300000_s_1", "ont_b300000_s_2"])
        self.assertEqual([r["PROFILE_TRUTH"] for r in rows], ["/truth/1", "/truth/2"])
        self.assertEqual({r["READ_TYPE"] for r in rows} | {r["SECOND"] for r in rows}, {"ont", "-"})
        complement = str.maketrans("ACGT", "TGCA")
        first = {}
        for row, (a, b) in zip(rows, ((0.5, 0.5), (0.8, 0.2))):
            # zstd by default (--read_compression), as the collector names it
            self.assertTrue(row["FIRST"].endswith(".fq.zst"), row["FIRST"])
            with open(row["FIRST"], "rb") as fh:
                self.assertEqual(fh.read(4), compressed.ZSTD_MAGIC)
            lines = compressed.read_text(row["FIRST"]).splitlines()
            first[row["SAMPLEID"]] = lines
            names, reads = lines[0::4], lines[1::4]
            self.assertEqual(names, [f"@g{n[2]}x_{i}" for i, n in enumerate(names, 1)], "g<genome>x_<n>, n in order")
            # The reads' bases reach the sample's; a read is 1 kb (the setup's mean, SD 0) or shorter at a
            # contig's end, from either strand of its genome.
            self.assertGreaterEqual(sum(map(len, reads)), 300000)
            self.assertLess(sum(map(len, reads)), 300000 + 1000)
            self.assertTrue(all(100 <= len(r) <= 1000 for r in reads))
            strands = collections.Counter()
            for name, read in zip(names, reads):
                genome = sequences["GA" if name.startswith("@g0x_") else "GB"]
                strands["+" if read in genome else "-" if read.translate(complement)[::-1] in genome else "?"] += 1
            self.assertEqual(strands["?"], 0)
            self.assertGreater(min(strands["+"], strands["-"]), len(reads) / 3)
            counts = {g: sum(n.startswith(f"@g{g}x_") for n in names) for g in (0, 1)}
            share_a = a * 20000 / (a * 20000 + b * 40000)
            self.assertEqual(counts[0] + counts[1], len(names))
            self.assertAlmostEqual(counts[0] / len(names), share_a, delta=0.1)
        self.assertTrue(collect.simulated(unit, opts))
        # The same seed, the same reads, on another number of slots; gzip (BGZF) instead of zstd, the same reads.
        self.assertEqual(collect.simulate_long([(0, unit)], opts, 1), {})
        for row in rows:
            self.assertEqual(compressed.read_text(row["FIRST"]).splitlines(), first[row["SAMPLEID"]])
        gzipped = argparse.Namespace(**vars(opts), read_compression="gzip")
        self.assertEqual(collect.simulate_long([(0, unit)], gzipped, 3), {})
        for row in collect.unit_map_rows(unit, gzipped)[0]:
            self.assertTrue(row["FIRST"].endswith(".fq.gz"))
            with gzip.open(row["FIRST"], "rt") as fh:
                self.assertEqual(fh.read().splitlines(), first[row["SAMPLEID"]])
        # A host share: the host genome (scenarios.prepare_host) last among the genomes, that share of the weight.
        host_fa = os.path.join(root, "host.fa")
        host_seq = "".join(rng.choice("ACGT") for _ in range(30000))
        with open(host_fa, "w") as fh:
            fh.write(">chr1\n" + host_seq + "\n")
        hosted = argparse.Namespace(**vars(opts), host_folder=scenarios.prepare_host(host_fa, os.path.join(root, "host")))
        unit_h = dict(unit, name="ont_host", host_share=0.5, point=dict(unit["point"], name="ont_host"))
        self.assertEqual(collect.simulate_long([(0, unit_h)], hosted, 2), {})
        for row in collect.unit_map_rows(unit_h, hosted)[0]:
            lines = compressed.read_text(row["FIRST"]).splitlines()
            host_reads = [r for n, r in zip(lines[0::4], lines[1::4]) if n.startswith("@g2x_")]
            self.assertAlmostEqual(len(host_reads) / len(lines[0::4]), 0.5, delta=0.12)
            self.assertTrue(all(r in host_seq or r.translate(complement)[::-1] in host_seq for r in host_reads))

    def test_scheduler(self):
        # The simulations' queue: the ready job of highest priority first, on its slots; a job waits for those it
        # comes after; a job's new jobs join; a failure starts nothing more and says why. The jobs record their starts;
        # what may run at once follows from the slots, not from the threads' timing.
        gate = Gate()
        gate.release()  # these jobs do not wait
        s = collect.Scheduler(1)  # one slot: one job at a time, by priority
        for name, priority in (("low", 1), ("high", 5), ("mid", 3)):
            s.add(name, gate.job(name), priority=priority)
        self.assertEqual(s.run(), {})
        self.assertEqual(gate.started(), ["high", "mid", "low"])
        gate = Gate()
        gate.release()
        s = collect.Scheduler(2)
        s.add("low", gate.job("low"), priority=1)
        s.add("high", gate.job("high"), priority=5)
        s.add("big", gate.job("big"), need=2, priority=3)
        s.add("after", gate.job("after", new=[{"name": "child", "run": gate.job("child")}]), after=["big"], priority=9)
        self.assertEqual(s.run(), {})
        started = gate.started()
        self.assertEqual(set(started[:2]), {"high", "low"})  # low fills the slot big cannot use yet
        self.assertEqual(started[2:], ["big", "after", "child"])
        # big, on both slots, ran alone: high and low had ended
        self.assertLess(max(gate.events.index(("end", n)) for n in ("high", "low")), gate.events.index(("start", "big")))
        gate = Gate()
        gate.release()
        s = collect.Scheduler(1)
        s.add("bad", gate.job("bad", error="it broke"), priority=2)
        s.add("next", gate.job("next"), priority=1)
        s.add("waits", gate.job("waits"), after=["bad"])
        failures = s.run()
        self.assertEqual(failures, {"bad": "it broke"})
        self.assertEqual(gate.started(), ["bad"])
        # The paired-end points' threads: all slots between them by their work, at least one each, at most 16 per
        # sample.
        points = [(0, {"name": "deep", "read_pairs": "10000000", "read_length": "150", "samples": 2}),
                  (1, {"name": "shallow", "read_pairs": "1000", "read_length": "150", "samples": 12}),
                  (2, {"name": "mid", "read_pairs": "500000", "read_length": "150", "samples": 12})]
        threads = collect.pe_threads(points, 30)
        self.assertEqual(sum(threads.values()), 30)
        self.assertGreater(threads["deep"], threads["mid"])
        self.assertGreaterEqual(threads["shallow"], 1)
        self.assertEqual(collect.pe_threads(points[:1], 64), {"deep": 32})
        # A group runs at most its limit of jobs at once (the long reads' drawings), others fill the slots beside it:
        # the first four jobs wait at the gate, two drawings (the limit) and two others; then the rest runs.
        gate = Gate()
        s = collect.Scheduler(4, {"draw": 2})
        for i in range(5):
            s.add(f"d{i}", gate.job(f"d{i}", "draw"), priority=5, group="draw")
            s.add(f"o{i}", gate.job(f"o{i}", "other"), priority=1)
        result = {}
        runner = threading.Thread(target=lambda: result.update(failures=s.run()))
        runner.start()
        self.assertTrue(gate.wait_until(lambda: len(gate.started()) == 4), gate.events)
        self.assertEqual(sorted(gate.started()), ["d0", "d1", "o0", "o1"])
        self.assertEqual((gate.busy["draw"], gate.busy["other"]), (2, 2))
        gate.release()
        runner.join(60)
        self.assertEqual(result.get("failures"), {})
        self.assertEqual(sorted(gate.started()), sorted([f"d{i}" for i in range(5)] + [f"o{i}" for i in range(5)]))
        self.assertEqual(gate.most["draw"], 2)

    def test_workers(self):
        # The collector's Python work (the host's paired-end fragments) runs in worker processes once they are
        # started, else on the calling thread: the same results.
        self.assertEqual(collect.Workers.call(len, "abc"), 3)
        with collect.Workers.started(2):
            self.assertIsNotNone(collect.Workers.pool)
            self.assertEqual(collect.Workers.call(len, "abcd"), 4)
        self.assertIsNone(collect.Workers.pool)

    def test_room_on_the_disk(self):
        # With space (folder, bytes to keep free), a job starts only while the disk has room for it and the jobs
        # running, and one that opens new work keeps the bytes free besides; finishing work may use them. a runs
        # until c has started beside it (or a minute has passed: then the order below fails).
        free, events, cond = [100], [], threading.Condition()
        usage = collections.namedtuple("usage", "total used free")

        def job(name, until=None):
            def run():
                with cond:
                    events.append(("start", name))
                    cond.notify_all()
                    if until:
                        cond.wait_for(lambda: ("start", until) in events, timeout=60)
                    events.append(("end", name))
                return None, []
            return run
        original = shutil.disk_usage
        shutil.disk_usage = lambda path: usage(1000, 1000 - free[0], free[0])
        try:
            s = collect.Scheduler(4, space=(self.tmp.name, 10))
            s.add("a", job("a", until="c"), disk=50, priority=3)  # 100 free: room for 50 + 10
            s.add("b", job("b"), disk=80, priority=2)  # 100 - 50 < 80 + 10: waits for a
            s.add("c", job("c"), disk=40, priority=1, opens=False)  # 100 - 50 >= 40: beside a
            self.assertEqual(s.run(), {})
        finally:
            shutil.disk_usage = original
        self.assertLess(events.index(("start", "c")), events.index(("end", "a")))
        self.assertLess(events.index(("end", "a")), events.index(("start", "b")))

    def test_follow(self):
        # --follow: units profiled as their simulations end, while the simulations go on (a protal run as soon as a
        # point is simulated, here), each point's reads removed once its pe and se units are; once the simulations
        # have ended, the rest. Profiled again (another protal) with the reads removed: simulated again first.
        root = os.path.join(self.tmp.name, "follow")
        os.makedirs(root)
        runs = os.path.join(root, "runs.txt")
        protal = os.path.join(root, "protal")
        with open(protal, "w") as fh:  # a stand-in: a dump of each sample of the map, which needs its reads
            fh.write("#!" + sys.executable + "\nimport sys\na = sys.argv[1:]\nm = a[a.index('--map') + 1]\n"
                     "rows = [l.rstrip('\\n').split('\\t') for l in open(m) if not l.startswith('#')]\n"
                     f"open({runs!r}, 'a').write(' '.join(r[0] for r in rows) + '\\n')\n"
                     "for r in rows:\n"
                     "    open(r[1]).close()\n"
                     "    open(r[5] + '.truth_annotated', 'w').write('taxon_name\\ttruth\\nt\\t1\\n')\n"
                     "    print('Write truth to: ' + r[5])\n")
        os.chmod(protal, 0o755)
        db = os.path.join(root, "db.protal")
        open(db, "w").close()
        opts = argparse.Namespace(out=root, db=db, protal=protal, threads=1, profile_block=1e-9, poll=0.05,
                                  protal_lock=os.path.join(root, "protal.lock"))
        points = [{"name": f"rl100_p{n}", "read_length": "100", "read_pairs": str(n), "samples": 2} for n in (10, 20)]
        units = [u for p in points for u in ({"type": "pe", "point": p, "name": p["name"], "samples": 2},
                                              {"type": "se", "point": p, "name": p["name"] + "_se", "samples": 2})]
        keys = {p["name"]: {"point": p["name"]} for p in points}

        def simulate(point):
            base, sim, _ = collect.point_dirs(point, opts)
            os.makedirs(os.path.join(sim, "reads"), exist_ok=True)
            with open(os.path.join(sim, "protal.meta"), "w") as fh:
                fh.write(f"#OUTPUT_DIR\t{base}/protal\n#INPUT_DIR\t{sim}/reads\n"
                         "#SAMPLEID\tFIRST\tSECOND\tSAM\tPREFIX\tPROFILE\tPROFILE_TRUTH\n")
                for s in (1, 2):
                    name = f"{point['name']}_s_{s}"
                    fh.write(f"{name}\t{name}_R1.fq.gz\t{name}_R2.fq.gz\t{name}.sam.gz\t{name}\t{name}.profile\t/t\n")
                    for r in (1, 2):
                        with gzip.open(os.path.join(sim, "reads", f"{name}_R{r}.fq.gz"), "wt") as out:
                            out.write("@r\nACGT\n+\nIIII\n")
            collect.write_key(os.path.join(base, "simulated.json"), keys[point["name"]])

        simulator = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(120)"])
        collect.write_key(collect.simulating_file(opts), {"pid": simulator.pid, "host": socket.gethostname()})
        try:
            simulate(points[0])
            follower = threading.Thread(target=collect.follow, args=(units, opts, keys, None))
            follower.start()
            key_of = collect.profile_keys(units, opts, keys)
            deadline = time.time() + 60
            while not all(collect.profiled(u, opts, key_of[u["name"]]) for u in units[:2]) and time.time() < deadline:
                time.sleep(0.05)
            self.assertTrue(collect.simulations_running(opts))
            simulate(points[1])  # the second point is simulated after the first is profiled
        finally:
            simulator.kill()
            simulator.wait()
        os.remove(collect.simulating_file(opts))
        follower.join(60)
        self.assertFalse(follower.is_alive())
        self.assertFalse(collect.simulations_running(opts))
        with open(runs) as fh:
            self.assertEqual(fh.read().splitlines(), ["rl100_p10_s_1 rl100_p10_s_2 rl100_p10_s_1_se rl100_p10_s_2_se",
                                                      "rl100_p20_s_1 rl100_p20_s_2 rl100_p20_s_1_se rl100_p20_s_2_se"])
        self.assertTrue(all(collect.profiled(u, opts, key_of[u["name"]]) for u in units))
        self.assertEqual(glob.glob(os.path.join(root, "points", "*", "sim", "reads", "*")), [])
        with open(os.path.join(root, "points", "rl100_p20", "sim", "reads_removed.txt")) as fh:
            self.assertEqual(len(fh.read().splitlines()), 4)
        # A rerun finds every unit profiled: no protal run, nothing simulated (simulate_again is not called).
        console = io.StringIO()
        with contextlib.redirect_stdout(console):
            collect.follow(units, opts, keys, None)
        self.assertIn("every design point is profiled, in 0 protal runs", console.getvalue())
        with open(runs) as fh:
            self.assertEqual(len(fh.read().splitlines()), 2)
        # Another protal: everything is to be profiled again, but the reads are gone and nothing simulates.
        with open(protal, "a") as fh:
            fh.write("# another version\n")
        again = []

        def simulate_again(names):
            again.append(sorted(names))
            for point in points:
                if point["name"] in names:
                    simulate(point)
        collect.follow(units, opts, keys, simulate_again)
        self.assertEqual(again, [["rl100_p10", "rl100_p20"]])
        key_of = collect.profile_keys(units, opts, keys)
        self.assertTrue(all(collect.profiled(u, opts, key_of[u["name"]]) for u in units))
        self.assertEqual(collect.reads_removed(units, opts, keys, key_of), set())
        with open(runs) as fh:
            self.assertEqual(len(fh.read().splitlines()), 3)  # the two points in one run: nothing simulates

    def test_long_read_templates(self):
        # simulate_metagenomes draws a long-read sample's templates from the contigs of 100 bases or more (plain or
        # gzipped FASTA, any case, CRLF), of the setup's lengths, and a failure says why; the setups are checked.
        simulator = need_simulator()
        root = os.path.join(self.tmp.name, "templates")
        os.makedirs(root)
        rng = random.Random(3)
        records = [(">a one", "".join(rng.choice("acgt") for _ in range(3000))), (">short", "C" * 99),
                   (">exact", "G" * 100), (">tiny", "ACGT")]
        text = "".join(f"{header}\r\n" + "\r\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)) + "\r\n\r\n"
                       for header, seq in records)
        plain, zipped = os.path.join(root, "plain.fna"), os.path.join(root, "zipped.fna.gz")
        with open(plain, "w", newline="") as fh:
            fh.write(text)
        with gzip.open(zipped, "wt", newline="") as fh:
            fh.write(text)
        models = os.path.join(root, "models")
        os.makedirs(models)
        perfect_qshmm(os.path.join(models, "PERFECT.model"))
        opts = argparse.Namespace(out=root, seed=1, simulator=simulator, pbsim="pbsim", pbsim_models=models,
                                  read_compression="gzip")

        def unit(name, fasta, setup="qshmm:PERFECT:1000:0:0.97", bases=50000):
            point = os.path.join(root, "points", "p_" + name, "sim")
            os.makedirs(point, exist_ok=True)
            with open(os.path.join(point, "manifest.tsv"), "w") as fh:
                fh.write("sample\tgenome\tspecies\ttaxonomy\tgenome_length\tread_pairs\tvertical_coverage\t"
                         f"relative_abundance\tfastq_r1\tfastq_r2\tfasta_path\tart_seed\np_{name}_s_1\tG\tS G\t"
                         f"d__B;s__S G\t3200\t10\t1\t1.0\tr1\tr2\t{fasta}\t1\n")
            with open(os.path.join(point, "protal.meta"), "w") as fh:
                fh.write(f"#OUTPUT_DIR\t{point}\n#INPUT_DIR\t{point}\n#SAMPLEID\tFIRST\tSECOND\tSAM\tPREFIX\tPROFILE\t"
                         f"PROFILE_TRUTH\np_{name}_s_1\ta\tb\tc\td\te\t/truth\n")
            return {"type": "ont", "name": name, "bases": bases, "setup": collect.parse_long_setup(setup), "samples": 1,
                    "communities": [{"name": "p_" + name}],
                    "point": {"name": name, "read_length": "1000", "read_pairs": str(bases)}}

        def reads_of(name):
            with gzip.open(os.path.join(root, "points", name, "sim", "reads", f"{name}_s_1.fq.gz"), "rt") as fh:
                return fh.read().splitlines()

        long_one, exact = records[0][1].upper(), records[2][1]
        for name, fasta in (("plain", plain), ("zipped", zipped)):
            self.assertEqual(collect.simulate_long([(0, unit(name, fasta))], opts, 2), {})
            reads = reads_of(name)[1::4]
            self.assertTrue(reads)
            for read in reads:  # from "a one" (up to 1 kb, cut at its end) or the whole of "exact", either strand
                back = read.translate(str.maketrans("ACGT", "TGCA"))[::-1]
                self.assertTrue(read in long_one or back in long_one or exact in (read, back), read[:20])
                self.assertTrue(100 <= len(read) <= 1000)
            self.assertFalse(os.path.exists(os.path.join(root, "points", name, "sim", "tmp")))
        # Only short sequences: nothing to simulate from, and the run says so; so does a simulator that fails.
        short = os.path.join(root, "short.fna")
        with open(short, "w") as fh:
            fh.write(">tiny\nACGT\n>short\n" + "A" * 99 + "\n")
        failures = collect.simulate_long([(0, unit("short", short))], opts, 1)
        self.assertIn("no sequence of 100 bases or more", "".join(failures.values()))
        failing = os.path.join(root, "failing")
        with open(failing, "w") as fh:
            fh.write("#!/bin/sh\necho 'ERROR: out of ideas'\nexit 3\n")
        os.chmod(failing, 0o755)
        failures = collect.simulate_long([(0, unit("failing", plain))], argparse.Namespace(**{**vars(opts), "simulator": failing}), 1)
        self.assertIn("failed (3): ERROR: out of ideas", "".join(failures.values()))
        # A hifi setup: no pbsim3 model needed; reads g<genome>x_<n> in the file's order, with qualities of HiFi reads.
        self.assertEqual(collect.simulate_long([(0, unit("hifi", plain, setup="hifi:1000:0:3"))], opts, 1), {})
        lines = reads_of("hifi")
        self.assertEqual([line[1:] for line in lines[0::4]], [f"g0x_{i}" for i in range(1, len(lines) // 4 + 1)])
        quality = [ord(c) - 33 for line in lines[3::4] for c in line]
        self.assertGreater(sorted(quality)[len(quality) // 2], 25)
        for bad in ("hifi:1000:0", "hifi:1000:0:30:3", "hifi:a:0:3", "errhmm:ERRHMM-SEQUEL:15000:3000:0.99",
                    "qshmm:M:a:0:0.9"):
            with self.assertRaises(SystemExit):
                collect.parse_long_setup(bad)

    def test_relation_to_novel_species(self):
        # An absent taxon is put down to a species the database lacks when that species is at least as close to
        # it as every present one.
        lin = lambda text: lineages.from_string(text)
        taxon = lin("d__B;p__P;c__C;o__O;f__F;g__G;s__G a")
        in_sample = {"s__G b": lin("d__B;p__P;c__C;o__O;f__F;g__G;s__G b"),
                     "s__H c": lin("d__B;p__P;c__C;o__O2;f__F2;g__H;s__H c")}
        self.assertEqual(collect.relation(taxon, in_sample, {}), ("genus", ""))
        self.assertEqual(collect.relation(taxon, in_sample, {"s__G b": ("species", "s__G b")}), ("genus", "species"))
        self.assertEqual(collect.relation(taxon, in_sample, {"s__H c": ("class", "c__C")}), ("genus", ""))
        self.assertEqual(collect.relation(taxon, {"s__H c": in_sample["s__H c"]}, {"s__H c": ("order", "o__O2")}),
                         ("class", "order"))
        self.assertEqual(collect.relation(taxon, {"s__X": lin("d__A;p__Q;s__X")}, {"s__X": ("phylum", "p__Q")}), ("none", ""))
        # a present taxon's closest other species (meta_neighbour_rank)
        self.assertEqual(collect.relation(taxon, {"s__H c": in_sample["s__H c"]}, {})[0], "class")
        self.assertEqual(collect.relation(taxon, {}, {})[0], "none")

    def test_sample_relations_as_relation(self):
        # SampleRelations gives what relation() gives, for every taxon of a sample at once (a soil scenario's 10,000
        # species would otherwise cost each taxon 10,000 comparisons).
        rng = random.Random(4)

        def lineage():
            names = [f"d__D{rng.randrange(2)}"]
            for prefix in "pcofg":
                names.append(f"{prefix}__{prefix.upper()}{rng.randrange(3)}")
            names.append(f"s__S{rng.randrange(4)}")
            if rng.random() < 0.1:  # a lineage without a rank
                del names[rng.randrange(1, 6)]
            return lineages.from_string(";".join(names))
        for _ in range(30):
            in_sample = {}
            for _ in range(rng.randrange(0, 25)):
                lin = lineage()
                in_sample.setdefault(lin.get("species", f"s__x{len(in_sample)}") + f" {len(in_sample)}", lin)
            names = list(in_sample)
            novel = {s: (rng.choice(["species", "genus", "family"]), "c") for s in names if rng.random() < 0.4}
            relations = collect.SampleRelations(in_sample, novel)
            for _ in range(20):
                taxon = lineage()
                self.assertEqual(relations.relation(taxon), collect.relation(taxon, in_sample, novel))
            for name in names:
                others = {s: lin for s, lin in in_sample.items() if s != name}
                self.assertEqual(relations.neighbour(in_sample[name], name), collect.relation(in_sample[name], others, {})[0])

    def test_collector_scenario_units(self):
        # The scenarios' points and units after the design's: a community point per scenario (its paired-end reads
        # at the community's part of the depth, the host's added later), and the units that draw their reads from
        # its communities, each of its own seed.
        opts = argparse.Namespace(read_setups="150:HSXt:350:50", read_pairs="1000", read_types=["pe", "se", "pb", "ont"],
                                  samples=2, long_read_samples=0, pb_setup="hifi:15000:3000:3",
                                  ont_setup="qshmm:QSHMM-ONT-HQ:8000:6000:0.97", long_read_bases="1e6",
                                  scenarios="host,gut:1", scenario_samples=3, scenario_file=None, seed=1)
        points, units = collect.units_of(opts)
        self.assertEqual([p["name"] for p in points], ["rl150_p1000", "sc_host_pe_p10000000", "sc_gut_pe_p20000000"])
        host = points[1]
        self.assertEqual((host["community_pairs"], host["host_pairs"], host["samples"]), ("1000000", 9000000, 3))
        # Each sample at its own depth, from half to twice the scenario's, the community's part and the host's of it.
        factors = scenarios.depth_factors(1, "host", 3, 2.0)
        pairs = scenarios.sample_depths(10_000_000, factors)
        self.assertEqual([c + h for c, h in zip(host["community_pairs_of"], host["host_pairs_of"])], pairs)
        self.assertEqual(host["community_pairs_of"], [round(p * 0.1) for p in pairs])
        self.assertTrue(all(5_000_000 <= p <= 20_000_000 for p in pairs))
        self.assertEqual(collect.community_pairs_of(points[0]), [1000, 1000])  # the design's: one depth
        names = [u["name"] for u in units if u.get("scenario")]
        self.assertEqual(names, ["sc_host_pe_p10000000", "sc_host_se_ultima_r10000000", "sc_host_pb_b3000000000",
                                 "sc_host_ont_b3000000000", "sc_gut_pe_p20000000", "sc_gut_pb_b6000000000",
                                 "sc_gut_ont_b6000000000"])
        ultima = units[names.index("sc_host_se_ultima_r10000000") + len(units) - len(names)]
        self.assertTrue(collect.drawn(ultima))
        self.assertEqual((ultima["type"], ultima["bases"], ultima["setup"]["method"], ultima["host_share"]),
                         ("se", 3_000_000_000, "ultima", 0.9))
        # sample s of every technology at the same factor: the same bases as the paired-end sample s
        self.assertEqual(ultima["bases_of"], scenarios.sample_depths(3_000_000_000, factors))
        self.assertEqual(collect.bases_of(ultima), ultima["bases_of"])
        self.assertEqual(collect.profile_dir(ultima, argparse.Namespace(out="/o")), "/o/points/sc_host_se_ultima_r10000000/protal")
        self.assertFalse(collect.drawn(units[1]))  # the design's se: the first reads of its paired-end point
        self.assertEqual(collect.parse_long_setup("ultima:300:40:25:2"),
                         {"method": "ultima", "model": None, "length_mean": 300, "length_sd": 40, "q_mean": 25.0,
                          "q_sd": 2.0, "ratio": "", "text": "ultima:300:40:25:2"})
        # Without paired-end reads collected, a scenario's communities come from a point without reads.
        _, units = collect.units_of(argparse.Namespace(**{**vars(opts), "read_types": ["pb"], "scenarios": "gut"}))
        community = [u for u in units if u.get("scenario")][0]["communities"][0]
        self.assertEqual((community["name"], community["reads"]), ("sc_gut_community", False))
        # --samples 0: the scenarios alone.
        points, units = collect.units_of(argparse.Namespace(**{**vars(opts), "samples": 0}))
        self.assertTrue(all(u.get("scenario") for u in units) and all(p.get("scenario") for p in points))
        command, error = collect.simulation_command(host, 0, argparse.Namespace(**{**vars(opts), "out": "/o", "seed": 1,
                                                                                 "simulator": "sim", "genome_table": "g"}),
                                                    4, {})
        self.assertIsNone(error)
        self.assertEqual(command[command.index("--genome_table") + 1], "/o/scenarios/host/genomes.tsv")
        self.assertEqual(command[command.index("--total_read_pairs") + 1], ",".join(map(str, host["community_pairs_of"])))
        self.assertEqual(command[command.index("--species_per_sample") + 1], "2-50")
        self.assertIn("power_law", command)
        self.assertNotIn("--test", command)
        # depth_spread 1: every sample at the scenario's depth, one value for the simulator
        path = os.path.join(self.tmp.name, "flat.json")
        with open(path, "w") as fh:
            json.dump({"host": {"depth_spread": 1}}, fh)
        flat_opts = argparse.Namespace(**{**vars(opts), "scenario_file": path, "out": "/o", "simulator": "sim",
                                          "genome_table": "g"})
        flat = collect.units_of(flat_opts)[0][1]
        self.assertNotIn("community_pairs_of", flat)
        command, _ = collect.simulation_command(flat, 0, flat_opts, 4, {})
        self.assertEqual(command[command.index("--total_read_pairs") + 1], "1000000")
        self.assertEqual(collect.host_pairs_of(flat), [9000000] * 3)
        # --unmapped_reads (all; READ_TYPE; READ_TYPE:design; READ_TYPE:SCENARIO): those units' map rows write unmapped
        # records (UNMAPPED_READS), their profile keys say so (profiled again once asked for); the others' rows count
        # them, their keys unchanged.
        self.assertEqual(collect.unmapped_spec("pe:gut, se:host,none"), {("pe", "gut"), ("se", "host")})
        self.assertEqual(collect.unmapped_spec("all"), {(t, None) for t in collect.READ_TYPES})
        self.assertEqual(collect.unmapped_spec("pb,pe:design"), {("pb", None), ("pe", "design")})
        for bad in ("gut", "xx:gut", "pe:"):
            with self.assertRaises(ValueError):
                collect.unmapped_spec(bad)

        def marked_units(spec):
            marked = argparse.Namespace(**{**vars(opts), "unmapped_units": collect.unmapped_spec(spec)})
            return marked, [u["name"] for u in collect.units_of(marked)[1] if collect.writes_unmapped(u, marked)]
        self.assertEqual(marked_units("pe:gut")[1], ["sc_gut_pe_p20000000"])
        self.assertEqual(marked_units("pe:design")[1], ["rl150_p1000"])
        self.assertEqual(marked_units("pb")[1], ["pb_b1e6", "sc_host_pb_b3000000000", "sc_gut_pb_b6000000000"])
        self.assertEqual(len(marked_units("all")[1]), len(collect.units_of(opts)[1]))
        marked, _ = marked_units("pe:gut")
        _, units = collect.units_of(marked)
        self.assertFalse(collect.writes_unmapped(units[0], opts))  # without the option
        sim = os.path.join(self.tmp.name, "unmapped", "points", "sc_gut_pe_p20000000", "sim")
        os.makedirs(sim)
        with open(os.path.join(sim, "protal.meta"), "w") as fh:
            fh.write(f"#OUTPUT_DIR\t{sim}\n#INPUT_DIR\t{sim}\n#SAMPLEID\tFIRST\tSECOND\tSAM\tPREFIX\tPROFILE\tPROFILE_TRUTH\n"
                     "s1\ta1.fq.zst\ta2.fq.zst\ts1.sam.zst\ts1\ts1.profile\t/t\n")
        marked.out = os.path.join(self.tmp.name, "unmapped")
        gut = next(u for u in units if u["name"] == "sc_gut_pe_p20000000")
        rows, _ = collect.unit_map_rows(gut, marked)
        self.assertEqual([r["UNMAPPED_READS"] for r in rows], ["write"])
        self.assertEqual([r["UNMAPPED_READS"] for r in collect.unit_map_rows(gut, argparse.Namespace(**{
            **vars(marked), "unmapped_units": set()}))[0]], ["count"])
        path = os.path.join(self.tmp.name, "unmapped", "samples.map")
        collect.write_map(path, rows + [{c: v for c, v in rows[0].items() if c != "UNMAPPED_READS"}])
        with open(path) as fh:  # a row without the column (an older collector's map) counts them
            self.assertEqual([line.rstrip("\n").split("\t")[-1] for line in fh][1:], ["UNMAPPED_READS", "write", "count"])
        keys = {u["name"]: {"k": 1} for u in units}
        key_of = collect.profile_keys(units, argparse.Namespace(**{**vars(marked), "db": sim, "protal": path}), keys)
        self.assertEqual(key_of["sc_gut_pe_p20000000"].get("unmapped_reads"), "write")
        self.assertNotIn("unmapped_reads", key_of["sc_host_pe_p10000000"])


class ScenariosTest(unittest.TestCase):
    """scenarios.py: the scenarios' presets and files, their genome tables and the host genome."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def test_scenario_definitions_and_tables(self):
        # The scenarios' presets, the selection of them, a file that changes or adds some, and a scenario's genome
        # table, which gives its share of species the database lacks.
        defs = scenarios.definitions()
        # The presets the build names (--scenarios gut, soil, soil_shallow, host) are there and well-formed; each reads
        # its communities with paired-end reads, and every other technology at the paired-end reads' bases (sample s of
        # every technology at the same bases); one has host reads.
        self.assertTrue({"gut", "soil", "soil_shallow", "host"} <= set(defs), list(defs))
        for name, d in defs.items():
            scenarios.check_definition(name, d)
            reads = {r["type"]: r for r in d["reads"]}
            pe = reads["pe"]
            bases = pe["depth"] * 2 * pe["length"]
            for kind in ("pb", "ont"):
                if kind in reads:
                    self.assertEqual(reads[kind]["depth"], bases, (name, kind))
            if "se" in reads:
                length = collect.parse_long_setup(reads["se"]["setup"])["length_mean"]
                self.assertEqual(reads["se"]["depth"] * length, bases, name)
        self.assertTrue(any(d["host_share"] > 0 for d in defs.values()))
        # Each sample's depth factor: log-uniform from 1/spread to spread, one in each equal part of the range, the
        # same for the same seed and name, other ones for another seed.
        factors = scenarios.depth_factors(1, "soil", 6, 2.0)
        self.assertEqual(factors, scenarios.depth_factors(1, "soil", 6, 2.0))
        self.assertNotEqual(factors, scenarios.depth_factors(1001, "soil", 6, 2.0))
        parts = sorted(int((math.log2(f) + 1) / 2 * 6) for f in factors)
        self.assertEqual(parts, list(range(6)))
        self.assertEqual(scenarios.depth_factors(1, "soil", 3, 1.0), [1.0, 1.0, 1.0])
        self.assertEqual(scenarios.depth_factors(1, "soil", 0, 2.0), [])
        self.assertEqual(scenarios.sample_depths(1000, [0.5, 1.999, 0.0001]), [500, 1999, 1])
        self.assertEqual(scenarios.selection("gut,soil:5", 3, defs), [("gut", 3), ("soil", 5)])
        self.assertEqual([n for n, _ in scenarios.selection("all", 2, defs)], list(defs))
        self.assertEqual(scenarios.selection("gut:0,host", 2, defs), [("host", 2)])
        self.assertEqual(scenarios.selection("gut:0,host", 2, defs, keep_zero=True), [("gut", 0), ("host", 2)])
        self.assertEqual(scenarios.selection("none", 2, defs), [])
        for bad in ("gut,gut", "mars", "gut:x"):
            with self.assertRaises(scenarios.ScenarioError):
                scenarios.selection(bad, 1, defs)
        path = os.path.join(self.tmp.name, "scenarios.json")
        with open(path, "w") as fh:
            json.dump({"soil": {"species": "40-50"},
                       "tiny": {"species": "3", "novel_share": 0.5, "abundance": "", "strains": "", "congeners": "0",
                                "host_share": 0, "reads": [{"type": "pe", "length": 100, "profile": "HS20",
                                                            "depth": 1000}]}}, fh)
        # A scenario larger than its table can hold is scaled down to it: its largest sample a TABLE_MARGIN-th of the
        # table, its smallest in proportion; one that fits is kept. (A soil of 9,000-11,000 species at 60% novel.)
        soil = {**defs["soil"], "species": "9000-11000", "novel_share": 0.6}
        self.assertEqual(scenarios.table_split(5400, 2600, 0.6), (1733, 2600))
        scaled, why = scenarios.fit_species(soil, 5400, 2600)
        self.assertEqual(scaled["species"], f"{round(9000 * 2888 / 11000)}-2888")  # 4333 / 1.5
        self.assertIn("scaled from 9000-11000 to 2363-2888 species per sample", why)
        self.assertEqual({k: v for k, v in scaled.items() if k != "species"}, {k: v for k, v in soil.items() if k != "species"})
        self.assertEqual(scenarios.fit_species(soil, 20000, 9000), (soil, None))
        self.assertEqual(scenarios.fit_species({**defs["gut"], "species": "350-450"}, 2, 1)[0]["species"], "1")
        with self.assertRaises(scenarios.ScenarioError):
            scenarios.fit_species(soil, 100, 0)
        changed = scenarios.definitions(path)
        self.assertEqual((changed["soil"]["species"], changed["soil"]["novel_share"]), ("40-50", defs["soil"]["novel_share"]))
        self.assertEqual(changed["tiny"]["reads"][0], {**scenarios.ILLUMINA, "length": 100, "profile": "HS20",
                                                       "depth": 1000})
        for wrong in ({"species": "3"}, {**changed["tiny"], "colour": 1}, {**changed["tiny"], "novel_share": 2},
                      {**changed["tiny"], "host_share": 1}, {**changed["tiny"], "reads": [{"type": "se", "depth": 9,
                                                                                         "setup": "hifi:1:1:1"}]},
                      {**changed["tiny"], "reads": [{"type": "pe", "depth": 9}] * 2},
                      {**changed["tiny"], "depth_spread": 0.5}, {**changed["tiny"], "depth_spread": "x"}):
            with self.assertRaises(scenarios.ScenarioError):
                scenarios.check_definition("x", wrong)
        # The table's species: all of the side short of the share, enough of the other.
        self.assertEqual(scenarios.table_plan(1000, 100, 0.05, 400), (1000, 53))
        self.assertEqual(scenarios.table_plan(1000, 100, 0.6, 150), (67, 100))
        self.assertEqual(scenarios.table_plan(1000, 100, 0.0, 400), (1000, 0))
        with self.assertRaises(scenarios.ScenarioError) as raised:
            scenarios.table_plan(1000, 100, 0.6, 400)
        self.assertIn("--rep_only_species", str(raised.exception))
        # A genome table of 40 species, 2 genomes each, 10 of them held out: a table of 20% held out takes the 30 the
        # database has and 30 x 0.2 / 0.8 = 7.5, 8, of the 10.
        table = os.path.join(self.tmp.name, "scenario_genomes.tsv")
        with open(table, "w") as fh:
            for s in range(40):
                for g in range(2):
                    fh.write(f"G{s}_{g}\td__Bacteria;p__P;c__C;o__O;f__F;g__G{s % 7};s__G{s % 7} sp{s}\t/x/{s}_{g}.fna\t1000\n")
        novel = {f"s__G{s % 7} sp{s}": ("species", "") for s in range(10)}
        out = os.path.join(self.tmp.name, "scenario_tables", "t.tsv")
        d = {**defs["gut"], "species": "10-20", "novel_share": 0.2}
        note = scenarios.scenario_table(table, novel, "t", d, 1, out)
        self.assertIn("38 species (8 the database lacks", note)
        with open(out) as fh:
            lines = fh.read().splitlines()
        species = collections.Counter(line.split("\t")[1].split(";")[-1] for line in lines)
        self.assertEqual(len(species), 38)
        self.assertEqual(set(species.values()), {2})  # every genome of a species taken
        self.assertEqual(sum(s in novel for s in species), 8)
        before = os.stat(out).st_mtime_ns
        time.sleep(0.01)
        scenarios.scenario_table(table, novel, "t", d, 1, out)  # the same: not written again
        self.assertEqual(os.stat(out).st_mtime_ns, before)
        with self.assertRaises(scenarios.ScenarioError):
            scenarios.scenario_table(table, novel, "t", {**d, "species": "60"}, 1, out)

    def test_host_genome(self):
        # A host genome as plain sequence, read by memory map: templates from its contigs of 1 kb or more; and its
        # paired-end reads by ART in amplicon mode, from fragments of it.
        rng = random.Random(7)
        chromosomes = {"chr1": "".join(rng.choice("ACGT") for _ in range(30000)),
                       "short": "ACGT" * 100, "chr2": "".join(rng.choice("acgt") for _ in range(20000))}
        fasta = os.path.join(self.tmp.name, "host.fna.gz")
        with gzip.open(fasta, "wt") as fh:
            for name, seq in chromosomes.items():
                fh.write(f">{name} a chromosome\n" + "\n".join(seq[i:i + 80] for i in range(0, len(seq), 80)) + "\n")
        folder = scenarios.prepare_host(fasta, os.path.join(self.tmp.name, "host"))
        self.assertEqual(scenarios.host_identity(folder)["bases"], 50000)  # "short" left out
        host = scenarios.Host.of(folder)
        self.assertIs(host, scenarios.Host.of(folder))
        draws = random.Random(1)
        upper = {k: v.upper() for k, v in chromosomes.items()}
        for _ in range(200):
            seq = host.draw(draws, 500).decode()
            self.assertTrue(seq in upper["chr1"] or seq in upper["chr2"])
            self.assertTrue(100 <= len(seq) <= 500)  # never a few bases at a contig's end
        self.assertTrue(all(len(host.draw(draws, 15000)) >= 100 for _ in range(500)))
        before = os.stat(os.path.join(folder, "host.seq")).st_mtime_ns
        scenarios.prepare_host(fasta, folder)  # the same FASTA: kept
        self.assertEqual(os.stat(os.path.join(folder, "host.seq")).st_mtime_ns, before)
        if not prerequisites.on_path("art_illumina"):
            prerequisites.missing("no art_illumina for the host's paired-end reads")
        tmp = os.path.join(self.tmp.name, "host_chunk")
        task = {"sample": "s", "host": folder, "chunk": 1, "art": "art_illumina", "art_args": ["-ss", "HS20"],
                "pairs": 500, "length": 100, "fragment_mean": 300.0, "fragment_sd": 40.0, "seed": 3,
                "tmp": os.path.join(tmp, "c1"), "r1": os.path.join(tmp, "r1.fq.gz"), "r2": os.path.join(tmp, "r2.fq.gz")}
        self.assertIsNone(scenarios.host_pe_chunk(task))
        with gzip.open(task["r1"], "rt") as a, gzip.open(task["r2"], "rt") as b:
            first, second = a.read().splitlines(), b.read().splitlines()
        self.assertEqual((len(first) // 4, len(second) // 4), (500, 500))
        self.assertEqual([n.split("/")[0] for n in first[0::4]], [n.split("/")[0] for n in second[0::4]])
        self.assertTrue(all(len(r) == 100 for r in first[1::4]))
        self.assertFalse(os.path.exists(task["tmp"]))
        both = "".join(upper[c] for c in ("chr1", "chr2"))
        back = str.maketrans("ACGT", "TGCA")
        near = sum(r[:30] in both or r[:30].translate(back)[::-1] in both for r in first[1::4])
        self.assertGreater(near, 400)  # reads of the host, but for their errors
        zstd_task = {**task, "r1": os.path.join(tmp, "r1.fq.zst"), "r2": os.path.join(tmp, "r2.fq.zst")}
        self.assertIsNone(scenarios.host_pe_chunk(zstd_task))
        with open(zstd_task["r1"], "rb") as fh:
            self.assertEqual(fh.read(4), compressed.ZSTD_MAGIC)
        self.assertEqual(compressed.read_text(zstd_task["r1"]).splitlines(), first)
        self.assertEqual(compressed.read_text(zstd_task["r2"]).splitlines(), second)


class HifiReadsTest(unittest.TestCase):
    """hifi_reads.py: the PacBio HiFi and Ultima reads the collector makes without pbsim3."""

    def test_hifi_reads(self):
        # hifi_reads.py: a read of each template whose qualities say how many errors it has (their expected errors
        # are its errors, overall), reads of Q50 at 5 kb, Q30 at 25 kb and Q20 at 50 kb (SD 3 around), only indels in
        # homopolymers, the same reads for the same seed.
        import numpy as np
        import hifi_reads
        self.assertEqual(list(hifi_reads.length_quality([1000, 5000, 15000, 25000, 37500, 50000, 75000])),
                         [50, 50, 40, 30, 25, 20, 10])
        rng = np.random.default_rng(3)
        seqs, groups = [], []
        for length, count in ((5000, 100), (15000, 100), (25000, 100), (50000, 40)):
            for _ in range(count):  # runs of random bases, a tenth of them homopolymers of 3 to 8
                runs = np.where(rng.random(length) < 0.9, rng.geometric(0.6, length), rng.integers(3, 9, length))
                seqs.append(np.repeat(np.frombuffer(b"ACGT", np.uint8)[rng.integers(0, 4, length)], runs)[:length].tobytes())
                groups.append(length)
        groups = np.array(groups)
        reads, stats = hifi_reads.mutate(seqs, np.random.default_rng(1), 3)
        self.assertEqual(len(reads), len(seqs))
        self.assertAlmostEqual(stats["events"].sum() / stats["expected"].sum(), 1.0, delta=0.08)
        self.assertGreaterEqual(stats["q"].min(), hifi_reads.Q_MIN)
        for length, q in ((5000, 50), (15000, 40), (25000, 30), (50000, 20)):
            self.assertAlmostEqual(float(stats["q"][groups == length].mean()), q, delta=1.5, msg=length)
            self.assertAlmostEqual(float(stats["q"][groups == length].std()), 3, delta=1.0, msg=length)
        self.assertEqual(stats["homopolymer_substitutions"], 0)
        for (read, quality), seq in zip(reads, seqs):
            self.assertEqual(len(read), len(quality))
            self.assertLess(abs(len(read) - len(seq)), 0.05 * len(seq))
            self.assertTrue(34 <= min(quality) and max(quality) <= 126)  # Q1 to Q93
        again, _ = hifi_reads.mutate(seqs, np.random.default_rng(1), 3)
        self.assertEqual(again, reads)
        # A read's errors follow its quality: about 10^(-Q/10) per base.
        errors = stats["events"] / np.array([len(s) for s in seqs])
        for length, q in ((25000, 30), (50000, 20)):
            self.assertAlmostEqual(float(np.median(errors[groups == length])), 10 ** (-q / 10), delta=0.4 * 10 ** (-q / 10))

    def test_flow_reads(self):
        # hifi_reads.py's flow model (Ultima reads of the scenarios): a read's bases average its quality (q_mean, SD
        # q_sd), its errors are what its qualities expect, homopolymers count from two bases and take no
        # substitution.
        import numpy as np
        import hifi_reads
        rng = np.random.default_rng(2)
        seqs = []
        for _ in range(3000):
            runs = np.where(rng.random(300) < 0.85, 1, rng.integers(2, 7, 300))
            seqs.append(np.repeat(np.frombuffer(b"ACGT", np.uint8)[rng.integers(0, 4, 300)], runs)[:300].tobytes())
        reads, stats = hifi_reads.mutate(seqs, np.random.default_rng(1), 2.0, q_mean=25.0)
        self.assertEqual(len(reads), len(seqs))
        self.assertAlmostEqual(float(stats["q"].mean()), 25, delta=0.2)
        self.assertAlmostEqual(float(stats["q"].std()), 2, delta=0.2)
        mean_phred = np.array([np.mean(np.frombuffer(q, np.uint8).astype(float) - 33) for _, q in reads])
        self.assertLess(abs(float(np.median(mean_phred - stats["q"]))), 0.5)  # the bases average the read's quality
        self.assertAlmostEqual(stats["events"].sum() / stats["expected"].sum(), 1.0, delta=0.08)
        self.assertEqual(stats["homopolymer_substitutions"], 0)
        self.assertNotEqual(hifi_reads.FLOW_MODEL, hifi_reads.MODEL)


if __name__ == "__main__":
    unittest.main()
