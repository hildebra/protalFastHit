#!/usr/bin/env python3
"""test_ancestry_world.py - the true-positive test's pieces (docs/claude/2026-10-09-ancestry-true-positive-test/README.md):

  WorldTest    simulate_ancestry_world.py --quick: the release holds the database's genomes only, the roles sit where
               the trees put them, the planted genomes carry what they should, the samples are paired and balanced;
  OracleTest   ancestry_oracle.py on hand-made copies and reads: the sites, the consensus, a record's counts, the
               alleles a read is explained by and the polymorphic sites (the rules of AncestrySites.h and
               StrainAlleles.h, whose C++ tests these mirror);
  PipelineTest ancestry_truth_test.py run --quick end to end (needs $PROTAL and $SIMULATE; the models need scikit-learn
               in $PROTAL_TRAIN_PYTHON): protal's features equal the oracle's, the leaky control and the planted lines
               come out as built, the ancestry sites are the stem and the representative's lineage.

  python3 -m unittest scripts/mini_db/test_ancestry_world.py
  PROTAL=build/protal SIMULATE=build/simulate_metagenomes python3 -m unittest scripts/mini_db/test_ancestry_world.py
"""

import csv
import gzip
import os
import random
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.dirname(HERE)
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)
import prerequisites  # noqa: E402
import ancestry_oracle as oracle  # noqa: E402
import simulate_ancestry_world as world  # noqa: E402

GENERATOR = os.path.join(HERE, "simulate_ancestry_world.py")
RUNNER = os.path.join(SCRIPTS, "ancestry_truth_test.py")


def rows(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


class WorldTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.out = os.path.join(cls.tmp.name, "world")
        subprocess.run([sys.executable, GENERATOR, "--outdir", cls.out, "--quick", "--seed", "5"], check=True,
                       capture_output=True)
        sim = os.path.join(cls.out, "simulation")
        cls.genera = {int(r["genus"]): r for r in rows(os.path.join(sim, "genera.tsv"))}
        cls.roles = rows(os.path.join(sim, "roles.tsv"))
        cls.genomes = {r["accession"]: r for r in rows(os.path.join(sim, "all_genomes.tsv"))}
        cls.design = rows(os.path.join(cls.out, "samples", "design.tsv"))

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_release_holds_only_the_database_genomes(self):
        with open(os.path.join(self.out, "bac120_taxonomy_r226.tsv")) as fh:
            listed = {line.split("\t")[0][3:] for line in fh}
        query = {a for a, r in self.genomes.items() if r["in_release"] == "0"}
        self.assertTrue(query)
        self.assertFalse(listed & query)
        allele = {r["accession"] for r in self.roles if r["role"] == "allele"}
        self.assertTrue(allele <= listed)
        marker_dir = os.path.join(self.out, "genomic_files_all", "bac120_marker_genes_all_r226", "fna")
        first = sorted(os.listdir(marker_dir))[0]
        with open(os.path.join(marker_dir, first)) as fh:
            in_markers = {line[4:].strip() for line in fh if line.startswith(">")}
        self.assertFalse(in_markers & query)
        self.assertEqual(in_markers, listed)

    def test_roles_sit_where_the_trees_put_them(self):
        for r in self.roles:
            g = self.genera[int(r["genus"])]
            c_s, h_t, length = float(g["c_S"]), float(g["h_T"]), float(g["L"])
            if r["role"] == "Q_deep":
                self.assertAlmostEqual(float(r["mrca_rep"]), c_s, delta=2e-5)  # the crown: a long-diverged lineage
            if r["role"].startswith("N_") and not r["role"].startswith("N_imp"):
                b = float(r["b"])
                self.assertAlmostEqual(float(r["mrca_rep"]), h_t - b * length if b < 1 else c_s, delta=2e-5)
            if r["role"] == "Q_near":
                self.assertLess(float(r["mrca_nearest_allele"]), float(r["mrca_rep"]) + 1e-9)
            if r["class"] in ("Q", "N") and g["regime"] == "none":
                # without recombination a strain has the stem's derived states (a few repeat hits aside), and a
                # novel species keeps the ancestral base on the part of the stem it did not share
                stem = int(r["at_stem"])
                if r["class"] == "Q" and r["role"] != "Q_ils" and not r["role"].startswith("Q_imp"):
                    self.assertLess(stem, 0.02 * int(r["cons_sites"]))
                if r["role"] in ("N_0", "N_0.5"):
                    self.assertGreater(stem, 0.2 * int(r["cons_sites"]))

    def test_fixed_sites_tell_late_twins_from_deep_strains(self):
        for genus, g in self.genera.items():
            if g["regime"] != "none":
                continue
            by_role = {r["role"]: r for r in self.roles if int(r["genus"]) == genus}
            deep, twin = by_role["Q_deep"], by_role["N_0.8"]

            def fixed(r):
                return int(r["fixed_agree"]) / int(r["fixed_sites"])

            def every(r):
                return int(r["cons_agree"]) / int(r["cons_sites"])
            # A deep strain disagrees at the representative's lineage below the crown; the sites there count as
            # fixed unless an allele genome of the crown's other side makes them polymorphic.
            covered = int(deep["at_rep_lineage_poly"]) / int(deep["at_rep_lineage"])
            uncovered = int(deep["at_rep_lineage"]) - int(deep["at_rep_lineage_poly"])
            disagree = int(deep["fixed_sites"]) - int(deep["fixed_agree"])
            self.assertLessEqual(disagree, uncovered + 0.01 * int(deep["fixed_sites"]))
            self.assertGreater(disagree, 0.5 * uncovered)
            if covered > 0.99:
                self.assertGreater(fixed(deep), 0.97)
            self.assertAlmostEqual(fixed(twin), 0.8, delta=0.06 + (1 - covered) * 0.5)
            self.assertGreater(every(deep), every(twin))

    def test_planted_genomes(self):
        planted = [r for r in self.roles if r["role"].startswith(("A_a", "P_f"))]
        self.assertTrue(planted)
        for r in planted:
            if r["role"].startswith("A_a"):
                a = float(r["a"])
                self.assertAlmostEqual(int(r["cons_agree"]) / int(r["cons_sites"]), a, delta=0.01)
            self.assertGreaterEqual(int(r["diffs"]) / int(r["length"]), 0.0199)

    def test_imports(self):
        for r in self.roles:
            if r["role"] == "Q_imp0.2" and self.genera[int(r["genus"])]["regime"] == "none":
                self.assertEqual(int(r["genes_imported"]), round(0.2 * 120))
                self.assertGreater(int(r["imported_diffs"]), 0)

    def test_samples_are_paired_and_balanced(self):
        by_sample = {}
        for r in self.design:
            if r["genus"]:
                by_sample.setdefault((r["set"], r["sample"]), []).append(r)
        for (kind, _), focus in by_sample.items():
            self.assertEqual(sorted(int(r["genus"]) for r in focus), sorted(self.genera))
        for genus in self.genera:
            pe = [r for r in self.design if r["set"] == "pe" and r["genus"] == str(genus)]
            q = sum(r["class"] == "Q" for r in pe)
            self.assertEqual(q, len(pe) - len(pe) // 2)
        truth_dir = os.path.join(self.out, "samples", "pe", "truth")
        for (kind, sample), focus in by_sample.items():
            if kind != "pe":
                continue
            with open(os.path.join(truth_dir, sample + ".tsv")) as fh:
                present = {line.strip() for line in fh}
            for r in focus:
                target = self.genera[int(r["genus"])]["target"]
                self.assertEqual(target in present, r["class"] == "Q")
        manifest = rows(os.path.join(self.out, "samples", "pe", "manifest.tsv"))
        for column in ("sample", "genome", "species", "taxonomy", "genome_length", "read_pairs", "fasta_path"):
            self.assertIn(column, manifest[0])
            self.assertTrue(all(r[column] for r in manifest))

    def test_the_same_seed_makes_the_same_world(self):
        other = os.path.join(self.tmp.name, "again")
        subprocess.run([sys.executable, GENERATOR, "--outdir", other, "--quick", "--seed", "5"], check=True,
                       capture_output=True)
        for name in ("roles.tsv", "genera.tsv", "imports.tsv"):
            with open(os.path.join(self.out, "simulation", name)) as a, \
                 open(os.path.join(other, "simulation", name)) as b:
                self.assertEqual(a.read().replace(self.out, ""), b.read().replace(other, ""))


def random_gene(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def substitute_at(seq, positions, shift=1):
    s = list(seq)
    for p in positions:
        s[p] = "ACGT"[("ACGT".index(s[p]) + shift) % 4]
    return "".join(s)


class OracleTest(unittest.TestCase):
    def setUp(self):
        self.rng = random.Random(7)
        self.own = random_gene(self.rng, 900)

    def test_compare_finds_the_substitutions(self):
        changed = [30, 59, 300, 450, 700]
        other = substitute_at(self.own, changed)
        s = oracle.compare_copies(self.own, oracle.unique_kmers(self.own), other)
        self.assertEqual(s.positions, changed)
        self.assertEqual(s.bases, [oracle.CODE[other[p]] for p in changed])
        self.assertEqual(s.count(0, 900), 5)
        self.assertEqual(s.count(31, 59), 0)
        same = oracle.compare_copies(self.own, oracle.unique_kmers(self.own), self.own)
        self.assertEqual(same.compared, 900)
        self.assertEqual(same.positions, [])

    def test_consensus_keeps_the_sites_the_congeners_share(self):
        own = self.own
        congeners = [substitute_at(own, [100, 200, 800], 1), substitute_at(own, [100, 600, 800], 1),
                     substitute_at(own, [100, 800], 1), substitute_at(substitute_at(own, [100, 800], 1), [880], 2)]
        comparisons = []
        for k, c in enumerate(congeners):
            s = oracle.compare_copies(own, oracle.unique_kmers(own), c)
            s.congener = 10 + k
            comparisons.append(s)
        sites = oracle.consensus(len(own), comparisons)
        # 4 compared: a site needs all four on one other base (0.9 x 4 = 3.6); 880 is one congener's
        self.assertEqual(sites.positions, [100, 800])
        alone = oracle.consensus(len(own), comparisons[:1])
        self.assertEqual(alone.positions, [100, 200, 800])  # one congener: its differences
        two = oracle.consensus(len(own), comparisons[:2])
        self.assertEqual(two.positions, [100, 200, 800])  # two: the nearest's (the first of equal identity)

    def test_consensus_tolerates_one_congener_with_a_base_of_its_own(self):
        own = self.own
        # Six congeners: all differ at 100 (a site); at 300 five share a base and the first has a third (a site, from
        # six compared); at 400 five share a base and the first has the species' (no site); at 500 four share a base
        # and two have a third (no site).
        congeners = []
        for k in range(6):
            c = substitute_at(own, [100], 1)
            c = substitute_at(c, [300], 2 if k == 0 else 1)
            if k > 0:
                c = substitute_at(c, [400], 1)
            congeners.append(substitute_at(c, [500], 2 if k < 2 else 1))
        comparisons = []
        for k, c in enumerate(congeners):
            s = oracle.compare_copies(own, oracle.unique_kmers(own), c)
            s.congener = 20 + k
            comparisons.append(s)
        sites = oracle.consensus(len(own), comparisons)
        self.assertEqual(sites.positions, [100, 300])
        self.assertEqual(sites.bases[1], oracle.CODE[congeners[1][300]])  # the five's base
        # Without the first: five compared, all alike at 300 and 400; at 500 one of five has a base of its own, which
        # the nine-in-ten rule does not tolerate below six.
        five = oracle.consensus(len(own), comparisons[1:])
        self.assertEqual(five.positions, [100, 300, 400])
        self.assertTrue(oracle.agree(5, 1, 6))
        self.assertFalse(oracle.agree(5, 0, 6))
        self.assertFalse(oracle.agree(4, 2, 6))
        self.assertFalse(oracle.agree(4, 1, 5))
        self.assertTrue(oracle.agree(9, 0, 10))

    def test_the_family_consensus_polarises_a_small_genus(self):
        # One congener: its differences at 100 (the family carries the congener's base: a site), 200 (the family carries
        # the species' base: the congener's own change, no site) and 300 (a third base: no site); without the outgroup
        # all three, as before 2026-10-10.
        own = self.own
        congener = substitute_at(own, [100, 200, 300], 1)
        s = oracle.compare_copies(own, oracle.unique_kmers(own), congener)
        s.congener = 2
        outgroup = [oracle.CODE[b] for b in own]
        outgroup[100] = oracle.CODE[congener[100]]
        outgroup[300] = oracle.CODE[substitute_at(own, [300], 2)[300]]
        self.assertEqual(oracle.consensus(len(own), [s]).positions, [100, 200, 300])
        polarised = oracle.consensus(len(own), [s], outgroup)
        self.assertEqual(polarised.positions, [100])
        outgroup[200] = oracle.NO_BASE  # no consensus there: the fallback stands
        self.assertEqual(oracle.consensus(len(own), [s], outgroup).positions, [100, 200])

    def test_count_columns_weighs_mismatches_and_classifies_codons(self):
        # The case of tests/test_ColumnWeights.cpp: M K P G F M; a read with M->L at a conserved codon, a synonymous
        # AAA->AAG, a mismatch at a column without an estimate (not counted), and G->R at a conserved column.
        copy = "ATGAAACCCGGGTTTATG"
        columns = oracle.Columns(18)
        columns.within = [9, 9, 9, 2, 2, 2, 0, 0, 0, 12, 12, 12, 1, 1, 1, 9, 9, 9]
        columns.among = [5, 5, 5, 1, 1, 1, 0, 0, 0, 3, 3, 3, 1, 1, 1, 5, 5, 5]
        columns.aa = list(columns.within)
        read = list(copy)
        read[0], read[5], read[6], read[9] = "C", "G", "T", "A"
        counts = oracle.count_columns(columns, "1X4=1X1X2=1X8=", 1, "".join(read), copy)
        self.assertEqual(counts, [18, 15, 3 * 9 + 3 * 2 + 3 * 12 + 3 * 1 + 3 * 9, 3, 9 + 2 + 12, 5 + 1 + 3, 2, 1, 2, 2])
        self.assertEqual(oracle.count_columns(columns, "18M", 1, "*", copy)[1], 0)
        self.assertEqual(oracle.amino_acid(0, 3, 2), "M")
        self.assertEqual(oracle.amino_acid(3, 0, 0), "*")

    def test_count_walks_the_cigar(self):
        sites = oracle.Sites()
        sites.positions, sites.bases = [10, 20, 30, 40], [0, 1, 2, 3]
        seq = "A" * 60
        c = oracle.count(sites, "60M", 1, seq)
        self.assertEqual((c.sites, c.agree, c.congener), (4, 4, 0))
        # a mismatch at 20 with the congeners' base (C), one at 30 with another base
        read = list(seq)
        read[20], read[30] = "C", "T"
        c = oracle.count(sites, "20M1X9M1X29M", 1, "".join(read))
        self.assertEqual((c.sites, c.agree, c.congener), (4, 2, 1))
        # soft clip of 5: the read starts at reference 15
        c = oracle.count(sites, "5S40M", 16, "G" * 45)
        self.assertEqual((c.sites, c.agree), (3, 3))
        self.assertEqual(oracle.count(sites, "60M", 1, "*").sites, 0)

    def test_explain_and_best_allele(self):
        allele = oracle.parse_allele("0-200:10G,50T,120A")
        read = oracle.from_sam_record("10M1X39M1X49M", 1, "A" * 10 + "G" + "A" * 39 + "T" + "A" * 49)
        self.assertEqual([(d.pos, d.base) for d in read.diffs], [(10, 2), (50, 3)])
        self.assertEqual(oracle.explain(allele, read), (2, 0))  # 120 lies past the read
        idx, shift, explained = oracle.best_allele([allele], read)
        self.assertEqual((idx, shift, explained), (0, -2, 2))
        lacking = oracle.from_sam_record("130M", 1, "A" * 130)
        self.assertEqual(oracle.explain(allele, lacking), (0, 3))
        self.assertEqual(oracle.best_allele([allele], lacking)[0], -1)
        insertion = oracle.parse_allele("0-200:60i3")
        read = oracle.from_sam_record("62M3I40M", 1, "A" * 105)
        self.assertEqual(oracle.explain(insertion, read), (1, 0))  # within 4 bases, the same length

    def test_polymorphic_sites(self):
        alleles = [oracle.parse_allele("0-200:10G,50T"), oracle.parse_allele("20-200:50C,70d2")]
        read = oracle.from_sam_record("10M1X39M1X49M", 1, "A" * 10 + "G" + "A" * 39 + "G" + "A" * 49)
        poly = oracle.Polymorphism(alleles, read.begin, read.end)
        self.assertEqual([p for p, _, _ in poly.sites], [10, 50, 70, 71])
        sites, known, novel = oracle.count_sites(poly, read)
        self.assertEqual((sites, known, novel), (4, 1, 1))  # 10G known, 50G novel (T and C known), 70/71 the reference's
        self.assertEqual(poly.cover(10), 1)
        self.assertEqual(poly.cover(30), 2)
        self.assertEqual([oracle.fixed_weight(n) for n in (0, 1, 2, 4)], [0, 30, 40, 48])

    def test_note_record_counts_fixed_sites_apart(self):
        class Db:
            has_alleles = True
            alleles = {(1, 1): [oracle.parse_allele("0-100:20C")]}

            species_sites = oracle.Database.species_sites

            def __init__(self):
                self._sites = {}

            def sites(self, taxid, gene):
                s = oracle.Sites()
                s.positions, s.bases = [20, 40], [1, 1]
                return s
        e = oracle.Evidence()
        read = "A" * 20 + "C" + "A" * 19 + "C" + "A" * 39
        oracle.note_record(e, Db(), 1, 1, "20M1X19M1X39M", 1, read, 60, True)
        f = e.features()
        # site 20: polymorphic (the allele's C), the read the congener's base; site 40 fixed, the congener's base. The
        # ancestry counts take the sites the species' alleles share (since 2026-10-10): 40 alone; the fixed-site counts
        # all of them.
        self.assertEqual((e.ancestry_sites, e.ancestry_agree, e.ancestry_congener), (1, 0, 1))
        self.assertEqual(f["allele_explained_share"], 0.5)
        self.assertEqual(f["ancestry_fixed_gain"], 0.0)
        self.assertEqual(f["oracle_fixed_agreement"], 0.0)


@unittest.skipUnless(os.environ.get("PROTAL") and os.environ.get("SIMULATE") or prerequisites.REQUIRED,
                     "needs $PROTAL and $SIMULATE")
class PipelineTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        for name in ("PROTAL", "SIMULATE"):
            if not os.environ.get(name):
                prerequisites.missing(f"needs ${name}")
        cls.tmp = tempfile.TemporaryDirectory()
        cls.out = os.path.join(cls.tmp.name, "atp")
        command = [sys.executable, RUNNER, "run", "--outdir", cls.out, "--quick", "-t", "2", "--protal",
                   os.environ["PROTAL"], "--simulate", os.environ["SIMULATE"]]
        result = subprocess.run(command, capture_output=True, text=True)
        if result.returncode != 0:
            raise AssertionError(result.stdout[-3000:] + result.stderr[-3000:])
        cls.verdicts = {r["hypothesis"]: r for r in rows(os.path.join(cls.out, "report", "verdicts.tsv"))}

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def verdict(self, key):
        self.assertIn(key, self.verdicts)
        return self.verdicts[key]

    def test_features_equal_the_oracle(self):
        self.assertEqual(self.verdict("H1b/H2b oracle")["verdict"], "PASS", self.verdicts)

    def test_the_tables_are_the_trees(self):
        self.assertEqual(self.verdict("H1a table")["verdict"], "PASS", self.verdicts)
        # H3 compares protal's consensus rule with the trees' history: a finding (it recalls ~0.8 of the derived sites
        # with six diverse congeners), not a check of the code, which the oracle's equality is
        self.assertIn(self.verdict("H3 sites")["verdict"], ("PASS", "FAIL"))
        sites = rows(os.path.join(self.out, "report", "l0_sites.tsv"))
        self.assertTrue(sites)
        self.assertTrue(all(float(r["share_in_history"]) > 0.8 for r in sites if r["congeners"] == "6"))

    def test_controls_and_calibration(self):
        self.assertEqual(self.verdict("H1 leak")["verdict"], "PASS", self.verdicts)
        self.assertEqual(self.verdict("H2 calibration")["verdict"], "PASS", self.verdicts)


if __name__ == "__main__":
    unittest.main()
