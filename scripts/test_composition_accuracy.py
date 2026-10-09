#!/usr/bin/env python3
"""Tests of composition_accuracy.py on a hand-made collection (collect_training_data.py's layout): a paired-end point of
one sample and a PacBio point that replays its community, the numbers worked out in the comments. Standard library only.

Run: python3 -m unittest scripts/test_composition_accuracy.py
"""
import csv
import gzip
import math
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import composition_accuracy as acc  # noqa: E402

LOG_HEADER = "Predicted\tProbability\tRepGenome\tLineage\tAbundance\tName\tTaxID\tSummary\tVCov\tGenomeSize\tGenomeFragments\tGeneCov0\n"
COMPOSITION_HEADER = "Sample\tScannedFragments\tScannedReads\tScannedBases\tFragmentBaseShare\tUnknownShare\n"


def write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)


def log_line(predicted, taxid, name, depth, size):
    return f"{predicted}\t0.9\tGCF\td__B;s__{name}\t0\ts__{name}\t{taxid}\t{{}}\t{depth}\t{size}\t0\t0\n"


class CompositionAccuracyTest(unittest.TestCase):
    """The community of p_s_1: A (2 Mb, 10,000 pairs, coverage 1.5), B (4 Mb, 10,000 pairs, 0.75) and C (1 Mb, 5,000
    pairs, 1.5), held out of the database: cells A 0.4, B 0.2, C 0.4, all cells' average genome 2.0 Mb. 25,000 pairs read
    (no host), 7.5 Mb, no overlapping mates. protal calls A at depth 1.44 (4% low, 2 Mb) and B at 0.75 (its size 4.4 Mb,
    10% high), not X, C's relative (0.2, 1.5 Mb). Explained: 2.88 + 3.3 = 6.18 Mb of 7.5 = 0.824 (truth 0.8); average
    genome 6.18 / 2.19 = 2.8219 Mb (truth of A and B 2.6667); unknown 1.32 Mb, 0.46777 genomes beside 2.19: 0.17600
    (truth 0.4: C's genome is smaller than the called species'); missing species 0.46777 / 1.095 and / 0.75."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        root = self.root = os.path.join(self.tmp.name, "test")
        sim, out = os.path.join(root, "points", "p", "sim"), os.path.join(root, "points", "p", "protal")
        lineage = "d__B;p__P;c__C;o__O;f__F;g__G;s__{}"
        write(os.path.join(sim, "manifest.tsv"),
              "sample\tgenome\tspecies\ttaxonomy\tgenome_length\tread_pairs\tvertical_coverage\trelative_abundance\n"
              f"p_s_1\tgA\tA\t{lineage.format('A')}\t2000000\t10000\t1.5\t0.4\n"
              f"p_s_1\tgB\tB\t{lineage.format('B')}\t4000000\t10000\t0.75\t0.2\n"
              f"p_s_1\tgC\tC\t{lineage.format('C')}\t1000000\t5000\t1.5\t0.4\n")
        write(os.path.join(sim, "protal.meta"), f"#OUTPUT_DIR\t{out}\n#INPUT_DIR\t{sim}\n"
              "#SAMPLEID\tFIRST\tSECOND\tSAM\tPREFIX\tPROFILE\tPROFILE_TRUTH\n"
              "p_s_1\tr1.fq\tr2.fq\tp_s_1.sam.zst\tp_s_1\tp_s_1.profile\t/truth\n")
        write(os.path.join(out, "profiles", "p_s_1.profile.log"), LOG_HEADER + log_line(1, 1, "A", 1.44, 2000000) +
              log_line(1, 2, "B", 0.75, 4400000) + log_line(0, 9, "X", 0.2, 1500000))
        write(os.path.join(out, "profiles", "p_s_1.profile.composition"),
              COMPOSITION_HEADER + "p_s_1\t25000\t50000\t7500000\t1\t0.176\n")
        # u: PacBio reads of p_s_1's community, 10 Mb drawn by bases (A 0.4, B 0.4, C 0.2 of them): A's depth 2.0, B's 1.0.
        usim, uout = os.path.join(root, "points", "u", "sim"), os.path.join(root, "points", "u", "protal")
        write(os.path.join(usim, "samples.tsv"), "sample\treads\ttruth\tcommunity\nu_s_1\tu.fq\t/truth\tp_s_1\n")
        write(os.path.join(uout, "profiles", "u_s_1.profile.log"), LOG_HEADER + log_line(1, 1, "A", 2.0, 2000000) +
              log_line(1, 2, "B", 1.0, 4000000))
        write(os.path.join(uout, "profiles", "u_s_1.profile.composition"), COMPOSITION_HEADER + "u_s_1\t100\t100\t10000000\t1\tNA\n")
        self.heldout = os.path.join(self.tmp.name, "heldout_species.txt")
        write(self.heldout, "s__C\tspecies\ts__C\n")

    def run_script(self, *extra):
        out = os.path.join(self.tmp.name, "out.tsv")
        lines = acc.main(["--test", self.root, "--heldout", self.heldout, "--out", out, *extra])
        with open(out) as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        self.assertEqual(len(rows), len(lines))
        return lines

    def test_the_profiles_own_calls(self):
        (line,) = self.run_script()
        self.assertEqual((line["species"], line["absent_species"], line["called"], line["true_called"], line["missing_true"]),
                         (3, 1, 2, 2, 1))
        self.assertAlmostEqual(line["host_share"], 0)
        self.assertAlmostEqual(line["explained"], 0.824)
        self.assertAlmostEqual(line["true_explained"], 0.8)
        self.assertAlmostEqual(line["true_explained_by_calls"], 0.8)
        unknown = (7.5e6 - 6.18e6) / (6.18e6 / 2.19)
        self.assertAlmostEqual(line["unknown"], unknown / (2.19 + unknown))
        self.assertAlmostEqual(line["unknown"], 0.176, places=3)  # as the profile's composition says
        self.assertAlmostEqual(line["true_unknown"], 0.4)
        self.assertAlmostEqual(line["ags"], 6.18e6 / 2.19, places=3)
        self.assertAlmostEqual(line["true_ags"], 2.0e6)
        self.assertAlmostEqual(line["true_ags_called"], 1.6e6 / 0.6, places=3)
        self.assertAlmostEqual(line["missing_at_median"], unknown / 1.095)
        self.assertAlmostEqual(line["missing_at_lowest"], unknown / 0.75)
        self.assertAlmostEqual(line["depth_ratio"], (0.96 + 1.0) / 2)
        self.assertAlmostEqual(line["size_ratio"], (1.0 + 1.1) / 2)

    def test_a_models_calls_and_drawn_reads(self):
        calls = os.path.join(self.tmp.name, "trained_model_pb.calls.tsv.gz")
        with gzip.open(calls, "wt") as fh:
            fh.write("meta_design\tmeta_sample\tmeta_read_type\tmeta_scenario\ttaxon\ttruth\tset\tp\tknob\tcall\n"
                     "u\tu_s_1\tpb\t\t1\t1\ttest\t0.9\t0.5\t1\n"
                     "u\tu_s_1\tpb\t\t2\t1\ttest\t0.8\t0.5\t1\n"
                     "p\tp_s_1\tpe\t\t1\t1\ttest\t0.9\t0.5\t1\n")  # another read type's: left out
        (line,) = self.run_script("--calls", calls, "--read-type", "pb")
        self.assertEqual((line["point"], line["sample"], line["read_type"]), ("u", "u_s_1", "pb"))
        self.assertAlmostEqual(line["host_share"], 0)  # p_s_1 read its community's pairs alone
        self.assertAlmostEqual(line["explained"], 0.8)  # 4 + 4 Mb of 10
        self.assertAlmostEqual(line["true_explained"], 0.8)
        self.assertAlmostEqual(line["depth_ratio"], 1.0)
        self.assertAlmostEqual(line["size_ratio"], 1.0)
        self.assertAlmostEqual(line["ags"], 8e6 / 3, places=3)
        self.assertAlmostEqual(line["unknown"], 0.75 / 3.75)  # 2 Mb left, 0.75 genomes beside 3
        # The paired-end sample with the model's calls: A and X (a false call), not B.
        calls_pe = os.path.join(self.tmp.name, "trained_model.calls.tsv.gz")
        with gzip.open(calls_pe, "wt") as fh:
            fh.write("meta_design\tmeta_sample\tmeta_read_type\ttaxon\tset\tcall\n"
                     "p\tp_s_1\tpe\t1\ttest\t1\np\tp_s_1\tpe\t2\ttest\t0\np\tp_s_1\tpe\t9\ttest\t1\n")
        (line,) = self.run_script("--calls", calls_pe)
        self.assertEqual((line["called"], line["true_called"], line["missing_true"]), (2, 1, 2))
        self.assertAlmostEqual(line["explained"], (2.88e6 + 0.3e6) / 7.5e6)
        self.assertAlmostEqual(line["true_explained_by_calls"], 0.4)
        self.assertAlmostEqual(line["true_unknown_given_calls"], 0.6)

    def test_the_species_table(self):
        # --species-out: A, B and C present (C held out and not in the profile: no taxid, no depth). X is in the profile
        # but neither present nor called (the profile's own calls are A and B): no row.
        species = os.path.join(self.tmp.name, "species.tsv.gz")
        self.run_script("--species-out", species)
        with gzip.open(species, "rt") as fh:
            rows = {r["species"]: r for r in csv.DictReader(fh, delimiter="\t")}
        self.assertEqual(sorted(rows), ["s__A", "s__B", "s__C"])
        a, b, c = rows["s__A"], rows["s__B"], rows["s__C"]
        self.assertEqual((a["taxid"], a["genomes"], a["present"], a["in_database"], a["called"]), ("1", "gA", "1", "1", "1"))
        self.assertAlmostEqual(float(a["true_depth"]), 1.5)
        self.assertAlmostEqual(float(a["depth"]), 1.44)
        self.assertAlmostEqual(float(a["true_genome_length"]), 2e6)
        self.assertAlmostEqual(float(b["genome_size"]), 4.4e6)
        self.assertAlmostEqual(float(b["true_cell_share"]), 0.2)
        self.assertEqual(b["true_read_pairs"], "10000")
        self.assertEqual((c["present"], c["in_database"], c["called"], c["taxid"], c["depth"]), ("1", "0", "0", "", "NA"))
        # With the model's calls X is called: a row of its own, absent (no truth), with its depth.
        calls_pe = os.path.join(self.tmp.name, "trained_model.calls.tsv.gz")
        with gzip.open(calls_pe, "wt") as fh:
            fh.write("meta_design\tmeta_sample\tmeta_read_type\ttaxon\tset\tcall\n"
                     "p\tp_s_1\tpe\t1\ttest\t1\np\tp_s_1\tpe\t2\ttest\t0\np\tp_s_1\tpe\t9\ttest\t1\n")
        self.run_script("--calls", calls_pe, "--species-out", species)
        with gzip.open(species, "rt") as fh:
            rows = {r["species"]: r for r in csv.DictReader(fh, delimiter="\t")}
        x = rows["s__X"]
        self.assertEqual((x["present"], x["called"], x["true_depth"], x["depth"]), ("0", "1", "NA", "0.2"))
        self.assertEqual((rows["s__B"]["called"], rows["s__B"]["depth"]), ("0", "0.75"))  # in the profile, not called

    def test_compose_matches_protal(self):
        # Composition.h's unit test (test_Composition.cpp, TheSharesAsComputedByHand): A depth 10 at 2 Mb, B 5 at 4 Mb, C 1
        # without a size; 100 Mb read, 10% the mates' overlap.
        r = acc.compose([(10, 2e6), (5, 4e6), (1, -1)], 400000, 100000000, 0.9, True)
        ags = (10 * 2e6 + 5 * 4e6) / 15
        attributed = 40e6 + ags
        unknown = (90e6 - attributed) / ags
        self.assertAlmostEqual(r["ags"], ags)
        self.assertAlmostEqual(r["explained"], attributed / 90e6)
        self.assertAlmostEqual(r["unknown"], unknown / (16 + unknown))
        self.assertEqual(r["quantiles"], [2e6, 2e6, 2e6, 4e6, 4e6])
        self.assertEqual(acc.compose([], 1000, 300000, 1, True)["unknown"], 1)
        self.assertTrue(math.isnan(acc.compose([], 1000, 300000, 1, False)["unknown"]))

    def test_summary_and_skipped_samples(self):
        # A sample without its composition file (a profile of protal before 2026-10-08) is skipped and said so.
        os.remove(os.path.join(self.root, "points", "p", "protal", "profiles", "p_s_1.profile.composition"))
        summary = os.path.join(self.tmp.name, "summary.txt")
        self.assertEqual(self.run_script("--summary", summary), [])
        with open(summary) as fh:
            self.assertIn("skipped: 1 (no composition", fh.read())


if __name__ == "__main__":
    unittest.main()
