#!/usr/bin/env python3
"""ancestry_sites.py on a synthetic genus: a strain's read sides with its species at the sites where the species
differs from its nearest congener and its private differences fall at the species' polymorphic sites; a novel
congener's read sides with the congener at half of the sites, at fixed ones.

    python3 -m unittest scripts/test_ancestry_sites.py
"""
import csv
import gzip
import os
import sys
import tempfile
import unittest

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ancestry_sites as an  # noqa: E402

LETTERS = "ACGT"


def mutate(seq, positions, rng):
    s = list(seq)
    for p in positions:
        s[p] = LETTERS[(LETTERS.index(s[p]) + 1 + rng.integers(0, 3)) % 4]
    return "".join(s)


def rows_of(path):
    with gzip.open(path, "rt") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


class AncestrySites(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        rng = np.random.default_rng(7)
        cls.tmp = tempfile.TemporaryDirectory()
        d = cls.tmp.name
        t = "".join(rng.choice(list(LETTERS), size=600))
        cls.sites2 = sorted(rng.choice(600, 24, replace=False).tolist())  # T vs its nearest congener T2
        t2 = mutate(t, cls.sites2, rng)
        sites3 = sorted(rng.choice(600, 48, replace=False).tolist())
        t3 = mutate(t, sites3, rng)
        with open(os.path.join(d, "reference.fna"), "w") as fh:
            for name, seq in (("11_5", t), ("12_5", t2), ("13_5", t3)):
                fh.write(f">{name}\n{seq}\n")
        # The species' own alleles (the full reference: every genome's copy under the species' taxid): the
        # representative's copy again, and one strain differing at the reads' private sites and eight more sites
        # outside the sites against T2.
        cls.private = [p for p in range(103, 247, 37) if p not in cls.sites2 and p not in sites3][:2]
        extra = [p for p in range(300, 600, 35) if p not in cls.sites2 and p not in cls.private][:8]
        cls.poly = sorted(cls.private + extra)
        allele = mutate(t, cls.poly, rng)
        # The converter names each record's genome after the name (since 2026-10-08): the report leaves a read's own
        # source genome out of the alleles.
        with open(os.path.join(d, "full_reference.fna"), "w") as fh:
            for name, seq in (("11_5 GCA_000000011.1", t), ("11_5 GCA_000000012.1", allele), ("12_5 GCA_000000021.1", t2),
                              ("13_5 GCA_000000031.1", t3)):
                fh.write(f">{name}\n{seq}\n")
        with open(os.path.join(d, "internal_taxonomy.dmp"), "w") as fh:
            fh.write("id\tparent_id\texternal_id\tname\trank\tlevel\trep_genome\n")
            fh.write("1\t1\t\td__Bacteria\tdomain\t0\t\n10\t1\t\tg__G\tgenus\t5\t\n")
            for tid, name in ((11, "s__G T"), (12, "s__G T2"), (13, "s__G T3"), (14, "s__G X")):
                fh.write(f"{tid}\t10\t\t{name}\tspecies\t6\tGCA_{tid}\n")
        with open(os.path.join(d, "heldout_species.txt"), "w") as fh:
            fh.write("s__G X\tspecies\n")
        # Reads of 150 bases from gene position 101 (1-based): the strain's with two private differences at polymorphic
        # sites; the novel congener's with T2's base at every other site in the window and the same two private ones.
        lo, hi = 100, 250
        window = [p for p in cls.sites2 if lo <= p < hi]
        cls.in_window = len(window)
        strain = mutate(t, cls.private, rng)[lo:hi]
        novel = list(t)
        for p in window[::2]:
            novel[p] = t2[p]
        novel = mutate("".join(novel), cls.private, rng)[lo:hi]
        cls.expected_alt = len(window[::2])
        # error_reads.py's layout: a sample's FP and FN files; the FN file holds the strain's read twice over (the
        # same record in both files is counted once) and a read of the same strain sequence from a genome the full
        # reference lacks (r4), the FP file the novel read and a secondary record. r1 comes from the allele's genome
        # (xg), so that allele is left out for it: its private differences are then not polymorphic; r4 keeps it.
        point = os.path.join(d, "training", "point")
        os.makedirs(point)
        strain_line = f"r1\t0\t11_5\t{lo + 1}\t60\t150M\t*\t0\t0\t{strain}\t*\txg:Z:GCA_000000012.1\txs:Z:s__G T\txe:Z:FN:11,source:11\n"
        other_line = f"r4\t0\t11_5\t{lo + 1}\t60\t150M\t*\t0\t0\t{strain}\t*\txg:Z:GCA_000000013.1\txs:Z:s__G T\txe:Z:FN:11,source:11\n"
        with open(os.path.join(point, "s1.FN.sam"), "w") as fh:
            fh.write("@HD\tVN:1.6\n@SQ\tSN:11_5\tLN:600\n" + strain_line + other_line)
        with open(os.path.join(point, "s1.FP.sam"), "w") as fh:
            fh.write("@HD\tVN:1.6\n@SQ\tSN:11_5\tLN:600\n" + strain_line)
            fh.write(f"r2\t0\t11_5\t{lo + 1}\t60\t150M\t*\t0\t0\t{novel}\t*\txg:Z:GCA_14\txs:Z:s__G X (not in the database)\txe:Z:FP:11\n")
            fh.write(f"r3\t256\t11_5\t{lo + 1}\t60\t150M\t*\t0\t0\t{novel}\t*\txs:Z:s__G X (not in the database)\txe:Z:FP:11\n")
        cls.out = os.path.join(d, "logs", "ancestry_sites", "pe")
        an.main(["--sams", d, "--reference", os.path.join(d, "reference.fna"),
                 "--full-reference", os.path.join(d, "full_reference.fna"),
                 "--taxonomy", os.path.join(d, "internal_taxonomy.dmp"), "--heldout", os.path.join(d, "heldout_species.txt"),
                 "--out", cls.out])

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_fragments(self):
        rows = rows_of(self.out + ".fragments.tsv.gz")
        self.assertEqual([r["qname"] for r in rows], ["r1", "r4", "r2"])  # the secondary record and the repeat left out
        strain, other_strain, novel = rows
        self.assertEqual(strain["sample"], "training:s1")
        self.assertEqual(strain["role"], "FN own")
        self.assertEqual(strain["relation"], "own")
        self.assertEqual(novel["role"], "FP")
        self.assertEqual(novel["relation"], "genus")
        self.assertEqual(int(strain["sites1"]), self.in_window)
        self.assertEqual(int(strain["agree1"]), self.in_window)
        self.assertEqual(int(strain["alt1"]), 0)
        self.assertEqual(int(novel["sites1"]), self.in_window)
        self.assertEqual(int(novel["alt1"]), self.expected_alt)
        self.assertEqual(int(novel["agree1"]), self.in_window - self.expected_alt)
        self.assertEqual(int(strain["mismatches"]), 2)
        self.assertAlmostEqual(float(strain["nearest_identity"]), 1 - 24 / 600, places=2)
        # The alleles: the novel read's and r4's private differences lie at polymorphic sites, the novel read's congener
        # bases at fixed ones; r1 comes from the allele's own genome, which is left out for it, so its differences are
        # at fixed sites. The sites against T2 in the window are all fixed (the allele differs elsewhere).
        self.assertAlmostEqual(float(strain["allele_divergence"]), 10 / 600, places=2)
        self.assertEqual(int(strain["poly_mismatches"]), 0)
        self.assertEqual(int(strain["nonpoly_mismatches"]), 2)
        self.assertEqual(int(strain["poly_covered"]), 0)
        self.assertEqual(int(other_strain["poly_mismatches"]), 2)
        self.assertEqual(int(other_strain["nonpoly_mismatches"]), 0)
        self.assertEqual(int(other_strain["poly_covered"]), 2)
        self.assertEqual(int(other_strain["fixed1"]), self.in_window)
        self.assertEqual(int(novel["poly_mismatches"]), 2)
        self.assertEqual(int(novel["nonpoly_mismatches"]), self.expected_alt)
        self.assertEqual(int(strain["fixed1"]), self.in_window)
        self.assertEqual(int(strain["fixed_agree1"]), self.in_window)
        self.assertEqual(int(novel["fixed_alt1"]), self.expected_alt)
        self.assertEqual(int(novel["poly_covered"]), 2)

    def test_taxa_and_summary(self):
        groups = {r["group"]: r for r in rows_of(self.out + ".taxa.tsv.gz")}
        self.assertEqual(set(groups), {"FN own", "FP genus"})
        self.assertEqual(groups["FN own"]["agree1"], groups["FN own"]["sites1"])
        self.assertLess(int(groups["FP genus"]["agree1"]), int(groups["FP genus"]["sites1"]))
        with open(self.out + ".auc.tsv") as fh:
            aucs = {(r["min_sites"], r["identity_band"], r["signal"]): r for r in csv.DictReader(fh, delimiter="\t")}
        self.assertEqual(aucs[("3", "all", "species_base_at_congener_sites")]["auc"], "1.0000")
        # Two congeners only: protal's consensus sites fall back to the nearest congener's.
        self.assertEqual(aucs[("3", "all", "species_base_at_consensus_sites")]["auc"], "1.0000")
        self.assertEqual(groups["FN own"]["sites4"], groups["FN own"]["sites1"])
        self.assertEqual(groups["FP genus"]["alt4"], groups["FP genus"]["alt1"])
        self.assertEqual(aucs[("3", "all", "fixed_site_identity")]["auc"], "1.0000")
        self.assertEqual(aucs[("10", "all", "identity")]["taxa"], "2")
        with open(self.out + ".summary.txt") as fh:
            summary = fh.read()
        self.assertIn("2 SAMs of 1 samples, 3 counted records", summary)
        self.assertIn("alleles: 1 copies of 1 taxon-genes", summary)
        self.assertIn("1 with their genome named", summary)
        self.assertIn("own source genome's allele was left out: FN own 1", summary)
        self.assertIn("N >= 10: 2 taxa (1 FN)", summary)

    def test_auc(self):
        self.assertAlmostEqual(an.auc([1, 1, 0, 0], [0.9, 0.8, 0.2, 0.1]), 1.0)
        self.assertAlmostEqual(an.auc([1, 0, 1, 0], [0.5, 0.5, 0.5, 0.5]), 0.5)
        self.assertAlmostEqual(an.auc([1, 1, 0, 0], [0.1, 0.2, 0.8, 0.9]), 0.0)


if __name__ == "__main__":
    unittest.main()
