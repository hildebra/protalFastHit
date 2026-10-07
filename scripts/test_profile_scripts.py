#!/usr/bin/env python3
"""Tests of the scripts that read a protal run's profiles: recurrent_calls.py (species called thinly in several samples
beside the same abundant relative), prevalence_calls.py (calls adjusted by a species' prevalence across the run) and
protal_profile_utils (merging profiles into one table; shipped by the conda recipe and `just install`), on tiny
hand-made outputs whose expected numbers are worked out in the comments. Standard library only.

Run: python3 -m unittest scripts/test_profile_scripts.py
"""

import csv
import gzip
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
LINEAGE = "d__Bacteria;p__{p};c__C{p};o__O{p};f__{f};g__{g};s__{g} {s}"


def lineage(p, f, g, s):
    return LINEAGE.format(p=p, f=f, g=g, s=s)


def write_profile_log(path, rows):
    """A <sample>.profile.log as protal -o writes it: rows (taxid, name, lineage, hits, called, probability, abundance)."""
    with open(path, "w") as fh:
        fh.write("TaxID\tName\tLineage\tRepGenome\tAbundance\tPredicted\tProbability\tSummary\n")
        for taxid, name, lin, hits, called, p, abundance in rows:
            fh.write(f"{taxid}\t{name}\t{lin}\tGCF_{taxid}.1\t{abundance}\t{int(called)}\t{p}\t"
                     f"{{ Genes: 3, Total Hits: {hits}, VCOV: 0.1 }}\n")


def run(script, *args, cwd=None):
    return subprocess.run([sys.executable, os.path.join(HERE, script), *args], capture_output=True, text=True, cwd=cwd)


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


class RecurrentCallsTest(unittest.TestCase):
    """Species A (genus G1) is called on 2 and 1 hits in samples s1 and s2, each time beside its congener B (50 and 25
    hits): suspect, B its companion at genus rank. E (another genus of A's family) is thin in both too, but B has 30
    times E's hits in s1 only (E has 3 in s2, B fewer than 10 times that): no companion in every sample, not suspect.
    D, thin in s3 alone beside A (30 hits), is listed but not suspect; C, a species of another phylum, is never thin; X
    is not called."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        A, B = lineage("P1", "F1", "G1", "a"), lineage("P1", "F1", "G1", "b")
        C, D, E = lineage("P2", "F2", "G2", "c"), lineage("P1", "F1", "G1", "d"), lineage("P1", "F1", "G3", "e")
        run_dir = os.path.join(self.tmp.name, "run", "profiles")
        os.makedirs(run_dir)
        write_profile_log(os.path.join(run_dir, "s1.profile.log"), [
            ("11", "G1 a", A, 2, True, 0.61, 0.01), ("12", "G1 b", B, 50, True, 0.99, 0.3),
            ("13", "G2 c", C, 100, True, 0.99, 0.6), ("15", "G3 e", E, 1, True, 0.55, 0.01),
            ("19", "G9 x", lineage("P1", "F1", "G9", "x"), 1, False, 0.2, 0.0)])
        write_profile_log(os.path.join(run_dir, "s2.profile.log"), [
            ("11", "G1 a", A, 1, True, 0.52, 0.01), ("12", "G1 b", B, 25, True, 0.98, 0.9),
            ("15", "G3 e", E, 3, True, 0.7, 0.02)])
        write_profile_log(os.path.join(run_dir, "s3.profile.log"), [
            ("11", "G1 a", A, 30, True, 0.97, 0.5), ("13", "G2 c", C, 5, True, 0.9, 0.4),
            ("14", "G1 d", D, 1, True, 0.51, 0.01)])
        self.run_dir = os.path.join(self.tmp.name, "run")

    def test_thin_calls_beside_the_same_relative(self):
        out = os.path.join(self.tmp.name, "recurrent.tsv")
        result = run("recurrent_calls.py", self.run_dir, "-o", out)
        self.assertEqual(result.returncode, 0, result.stderr)
        rows = read_tsv(out)
        self.assertEqual([r["taxid"] for r in rows], ["11", "15", "14"])  # by thin calls, then taxid
        a, e, d = rows
        self.assertEqual((a["species"], a["samples_called"], a["thin_calls"], a["thin_samples"], a["hits"], a["probabilities"]),
                         ("G1 a", "3", "2", "s1,s2", "2,1", "0.61,0.52"))
        self.assertEqual((a["companion_taxid"], a["companion"], a["shared_rank"], a["companion_hits"], a["suspect"]),
                         ("12", "G1 b", "genus", "75", "1"))
        self.assertEqual((e["thin_calls"], e["companion_taxid"], e["suspect"]), ("2", "", "0"))
        self.assertEqual((d["thin_calls"], d["companion"], d["shared_rank"], d["companion_hits"], d["suspect"]),
                         ("1", "G1 a", "genus", "30", "0"))
        self.assertIn("3 samples; 3 species called on at most 3 hits in some sample, 1 of them", result.stdout)
        self.assertIn("G1 a (11): thin in 2 of 3 samples, beside G1 b (genus, 75 hits)", result.stdout)

    def test_options(self):
        # --min-samples 1: D is suspect too; --max-hits 1: each species thin in one sample (A in s2, E in s1, D in s3),
        # none suspect; --ratio 2: E's companion B in both samples (50 >= 2 x 1, 25 >= 2 x 3), so E is suspect too.
        out = os.path.join(self.tmp.name, "recurrent.tsv")
        self.assertEqual(run("recurrent_calls.py", self.run_dir, "-o", out, "--min-samples", "1").returncode, 0)
        self.assertEqual({r["taxid"]: r["suspect"] for r in read_tsv(out)}, {"11": "1", "15": "0", "14": "1"})
        self.assertEqual(run("recurrent_calls.py", self.run_dir, "-o", out, "--max-hits", "1").returncode, 0)
        self.assertEqual({r["taxid"]: (r["thin_samples"], r["suspect"]) for r in read_tsv(out)},
                         {"11": ("s2", "0"), "15": ("s1", "0"), "14": ("s3", "0")})
        self.assertEqual(run("recurrent_calls.py", self.run_dir, "-o", out, "--ratio", "2").returncode, 0)
        self.assertEqual({r["taxid"]: r["suspect"] for r in read_tsv(out)}, {"11": "1", "15": "1", "14": "0"})
        # One .profile.log given as a file; something else stops the script.
        log = os.path.join(self.run_dir, "profiles", "s1.profile.log")
        self.assertEqual(run("recurrent_calls.py", log, "-o", out).returncode, 0)
        self.assertEqual({r["thin_samples"] for r in read_tsv(out)}, {"s1"})
        self.assertNotEqual(run("recurrent_calls.py", os.path.join(self.tmp.name, "x.txt"), "-o", out).returncode, 0)
        empty = os.path.join(self.tmp.name, "empty")
        os.makedirs(empty)
        self.assertIn("no .profile.log files found", run("recurrent_calls.py", empty, "-o", out).stderr)


class PrevalenceCallsTest(unittest.TestCase):
    """Three samples. T is likely in s2 and s3 (0.9) and borderline in s1 (0.45); Q in s1 (0.55) and barely in s2
    (0.05); R in s1 alone. With the base rate 0.2 (--prior) and 2 pseudo-samples, T's prevalence for s1 is
    (0.9 + 0.9 + 2 x 0.2) / (2 + 2) = 0.55, odds 1.222 against the base's 0.25: a ratio of 4.9, capped at 4, so s1's T goes
    from odds 0.818 to 3.273, p 0.765957: called. Q's prevalence for s1 is (0.05 + 0 for s3 + 0.4) / 4 = 0.1125, a ratio
    of 0.507: by default (--direction up) not applied; with --direction both p 0.55 -> 0.382609, a call lost. R in one
    sample keeps its probability."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.run_dir = os.path.join(self.tmp.name, "run")
        os.makedirs(self.run_dir)
        T, Q, R = lineage("P1", "F1", "G1", "t"), lineage("P1", "F1", "G2", "q"), lineage("P2", "F2", "G3", "r")
        write_profile_log(os.path.join(self.run_dir, "s1.profile.log"), [
            ("1", "G1 t", T, 5, False, 0.45, 0.2), ("2", "G2 q", Q, 7, True, 0.55, 0.3), ("3", "G3 r", R, 9, True, 0.6, 0.5)])
        write_profile_log(os.path.join(self.run_dir, "s2.profile.log"), [
            ("1", "G1 t", T, 50, True, 0.9, 0.75), ("2", "G2 q", Q, 1, False, 0.05, 0.25)])
        write_profile_log(os.path.join(self.run_dir, "s3.profile.log"), [("1", "G1 t", T, 40, True, 0.9, 1.0)])

    def adjusted(self, *options):
        out = os.path.join(self.tmp.name, "prevalence.tsv")
        result = run("prevalence_calls.py", self.run_dir, "-o", out, "--prior", "0.2", *options)
        self.assertEqual(result.returncode, 0, result.stderr)
        return {(r["sample"], r["species"]): r for r in read_tsv(out)}, result.stdout

    def test_prevalent_species_boosted_by_default(self):
        rows, stdout = self.adjusted("--profiles", os.path.join(self.tmp.name, "profiles"))
        t = rows[("s1", "G1 t")]
        self.assertEqual((t["prevalence"], t["p_adjusted"], t["called"], t["called_adjusted"]), ("0.55", "0.765957", "0", "1"))
        q = rows[("s1", "G2 q")]
        self.assertEqual((q["prevalence"], q["p_adjusted"], q["called_adjusted"]), ("0.1125", "0.55", "1"))
        r = rows[("s1", "G3 r")]
        self.assertEqual((r["prevalence"], r["p_adjusted"], r["called_adjusted"]), ("", "0.6", "1"))
        self.assertEqual(len(rows), 6)
        self.assertIn("3 samples, 6 taxa with records, base rate 0.2000 (given); calls 4 -> 5: 1 gained, 0 lost", stdout)
        # The profiles with the adjusted calls, abundances renormalised over them: s1 now T, Q and R (0.2 + 0.3 + 0.5).
        with open(os.path.join(self.tmp.name, "profiles", "s1.profile")) as fh:
            profile = [line.rstrip("\n").split("\t") for line in fh]
        self.assertEqual([(rep, a) for rep, _, a in profile], [("GCF_1.1", "0.2"), ("GCF_2.1", "0.3"), ("GCF_3.1", "0.5")])
        with open(os.path.join(self.tmp.name, "profiles", "s2.profile")) as fh:
            self.assertEqual(fh.read(), f"GCF_1.1\t{lineage('P1', 'F1', 'G1', 't')}\t1\n")  # Q's call stays off

    def test_the_full_update_cuts_rare_species(self):
        rows, stdout = self.adjusted("--direction", "both", "--check")
        q = rows[("s1", "G2 q")]
        self.assertEqual((q["p_adjusted"], q["called"], q["called_adjusted"]), ("0.382609", "1", "0"))
        self.assertEqual(rows[("s1", "G1 t")]["p_adjusted"], "0.765957")
        self.assertIn("calls 4 -> 4: 1 gained, 1 lost", stdout)
        self.assertIn("taxa in one sample only: 1 of 3", stdout)
        # Only the cuts (down): T stays uncalled, Q is cut.
        rows, _ = self.adjusted("--direction", "down")
        self.assertEqual((rows[("s1", "G1 t")]["called_adjusted"], rows[("s1", "G2 q")]["called_adjusted"]), ("0", "0"))
        # No cap: T's ratio 4.889 in full.
        rows, _ = self.adjusted("--max-odds-ratio", "0")
        self.assertEqual(rows[("s1", "G1 t")]["p_adjusted"], "0.8")  # odds 0.818 x 4.889 = 4.0

    def test_base_rate_and_one_sample(self):
        out = os.path.join(self.tmp.name, "prevalence.tsv")
        result = run("prevalence_calls.py", self.run_dir, "-o", out)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("base rate 0.5750 (the run mean)", result.stdout)  # (0.45 + 0.55 + 0.6 + 0.9 + 0.05 + 0.9) / 6
        one = run("prevalence_calls.py", os.path.join(self.run_dir, "s1.profile.log"), "-o", out)
        self.assertNotEqual(one.returncode, 0)
        self.assertIn("prevalence needs two samples or more", one.stderr)


class ProfileUtilsTest(unittest.TestCase):
    """protal_profile_utils merge: profiles (representative, lineage, abundance) into one table of lineages by sample,
    sample names from the file names, gzipped or not; duplicate names stop it unless --resolve-samples tells them apart
    by their folders."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def merge(self, *args):
        return run("protal_profile_utils", "merge", *args)

    def test_merge(self):
        a, b = os.path.join(self.tmp.name, "s1.profile"), os.path.join(self.tmp.name, "s2.profile.gz")
        with open(a, "w") as fh:
            fh.write("GCF_1.1\td__B;s__A\t0.5\nGCF_2.1\td__B;s__B\t0.25\nGCF_9.1\td__B;s__B\t0.25\n\n")
        with gzip.open(b, "wt") as fh:
            fh.write("GCF_1.1\td__B;s__A\t0.25\nGCF_3.1\td__B;s__C\t0.75\n")
        result = self.merge("--input", os.path.join(self.tmp.name, "*.profile*"))
        self.assertEqual(result.returncode, 0, result.stderr)
        # The same lineage twice in a profile is summed; one a sample lacks is 0.
        self.assertEqual(result.stdout, "taxon\ts1\ts2\nd__B;s__A\t0.5\t0.25\nd__B;s__B\t0.50\t0\nd__B;s__C\t0\t0.75\n")

    def test_sample_names(self):
        for folder in ("run1", "run2"):
            os.makedirs(os.path.join(self.tmp.name, folder))
            with open(os.path.join(self.tmp.name, folder, "x.profile"), "w") as fh:
                fh.write(f"GCF_1.1\td__B;s__A\t{0.5 if folder == 'run1' else 1}\n")
        paths = [os.path.join(self.tmp.name, folder, "x.profile") for folder in ("run1", "run2")]
        result = self.merge("--input", *paths)
        self.assertEqual(result.returncode, 2)
        self.assertIn("Duplicate sample names detected", result.stderr)
        self.assertIn("Use --resolve-samples", result.stderr)
        result = self.merge("--input", *paths, "--resolve-samples")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout.splitlines()[0], "taxon\trun1_x\trun2_x")
        # Files of one run's naming: the common part of their names dropped, the rest kept as short as tells them apart.
        names = [os.path.join(self.tmp.name, f"study_{s}_R1.profile") for s in ("a", "b")]
        for name in names:
            with open(name, "w") as fh:
                fh.write("GCF_1.1\td__B;s__A\t1\n")
        result = self.merge("--input", *names, "--resolve-samples")
        self.assertEqual(result.stdout.splitlines()[0], "taxon\ta\tb")

    def test_malformed_profiles(self):
        bad = os.path.join(self.tmp.name, "bad.profile")
        for text, why in (("GCF_1.1\td__B;s__A\n", "expected 3 tab-delimited columns"),
                          ("GCF_1.1\t\t0.5\n", "empty taxonomy lineage"), ("GCF_1.1\td__B;s__A\tmany\n", "invalid abundance")):
            with open(bad, "w") as fh:
                fh.write(text)
            result = self.merge("--input", bad)
            self.assertEqual(result.returncode, 2, text)
            self.assertIn(why, result.stderr)
        self.assertNotEqual(self.merge("--input", os.path.join(self.tmp.name, "missing.profile")).returncode, 0)


if __name__ == "__main__":
    unittest.main()
