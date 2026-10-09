#!/usr/bin/env python3
"""ancestry_sites.py on a synthetic genus: a strain's read sides with its species at the sites where the species
differs from its nearest congener and its private differences fall at the species' polymorphic sites; a novel
congener's read sides with the congener at half of the sites, at fixed ones. A run on several SAM folders writes at
each prefix what a run on that folder alone writes, byte for byte.

    python3 -m unittest scripts/test_ancestry_sites.py
"""
import contextlib
import csv
import gzip
import io
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

    def test_alleles_only_of_the_genomes_that_give_strain_alleles(self):
        # --allele-genome-share: the report takes alleles only from the genomes protal --build takes its strain alleles
        # from (gtdb_to_protal_db.allele_genome), the others being the ones the build simulates strains from.
        from gtdb_to_protal_db import allele_genome
        d = self.tmp.name
        reps = an.load_copies(os.path.join(d, "reference.fna"), {"11_5"})
        full = os.path.join(d, "full_reference.fna")
        everyone = an.load_alleles(full, {"11_5"}, reps, 6)
        self.assertEqual([genome for genome, _ in everyone["11_5"]], ["GCA_000000012.1"])
        self.assertEqual(an.load_alleles(full, {"11_5"}, reps, 6, 1e-12).get("11_5", []), [])
        half = an.load_alleles(full, {"11_5"}, reps, 6, 0.5).get("11_5", [])
        self.assertEqual(bool(half), allele_genome("GCA_000000012.1", 0.5))

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


OUTPUTS = (".fragments.tsv.gz", ".taxa.tsv.gz", ".auc.tsv", ".summary.txt")
X = "s__G X (not in the database)"  # the held-out congener the false positives' reads come from


def genome(t, k):
    """The accession of species t's k-th genome (k = 0: its representative)."""
    return f"GCA_{t:06d}{k:03d}.1"


def sam_line(qname, header, lo, seq, xs, xe, xg="", flag=0, mapq=60):
    """A record of error_reads.py's SAMs: 150 bases at gene position lo + 1, its source species, reasons and genome."""
    return (f"{qname}\t{flag}\t{header}\t{lo + 1}\t{mapq}\t150M\t*\t0\t0\t{seq}\t*\t" + (f"xg:Z:{xg}\t" if xg else "") +
            f"xs:Z:{xs}\txe:Z:{xe}\n")


def run(args):
    """ancestry_sites.main on args: its stdout."""
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        an.main(args)
    return out.getvalue()


class SeveralFolders(unittest.TestCase):
    """One run on two SAM folders (two read types' error reads) against a run on each alone: genus 10 has species 11-15
    in the training database and 16 held out (the false positives' reads come from it), more than --max-congeners 2,
    so that each folder's congeners are drawn at random; genus 20 has species 21 alone, without a congener. The
    folders' taxa overlap (11 in both, as its own taxon in one and as a congener in the other), their samples share a
    name and a record, and the full reference holds up to three other genomes' copies per species-gene (one twice),
    more than --max-alleles 2 for some."""

    @classmethod
    def setUpClass(cls):
        rng = np.random.default_rng(23)
        cls.tmp = tempfile.TemporaryDirectory()
        d = cls.tmp.name
        base = {g: "".join(rng.choice(list(LETTERS), size=600)) for g in (5, 6)}
        differences = {11: 0, 12: 15, 13: 30, 14: 45, 15: 60, 16: 20}
        copy = {(t, g): mutate(base[g], rng.choice(600, n, replace=False).tolist(), rng)
                for t, n in differences.items() for g in (5, 6)}
        copy[(21, 5)] = "".join(rng.choice(list(LETTERS), size=600))
        kept = [key for key in copy if key[0] != 16]
        with open(os.path.join(d, "reference.fna"), "w") as fh:
            for t, g in kept:
                fh.write(f">{t}_{g}\n{copy[(t, g)]}\n")
        # The species' other genomes: k = 1, 2, 3 per species-gene, 12_5's third the same as its first.
        allele = {}
        for (t, g), n in (((11, 5), 3), ((11, 6), 2), ((12, 5), 3), ((13, 6), 1), ((21, 5), 2)):
            for k in range(1, n + 1):
                allele[(t, g, k)] = mutate(copy[(t, g)], rng.choice(600, 4 * k, replace=False).tolist(), rng)
        allele[(12, 5, 3)] = allele[(12, 5, 1)]
        with open(os.path.join(d, "full_reference.fna"), "w") as fh:
            for t, g in kept:
                fh.write(f">{t}_{g} {genome(t, 0)}\n{copy[(t, g)]}\n")
            for k in (1, 2, 3):
                for (t, g, j), seq in allele.items():
                    if j == k:
                        fh.write(f">{t}_{g} {genome(t, k)}\n{seq}\n")
        with open(os.path.join(d, "internal_taxonomy.dmp"), "w") as fh:
            fh.write("id\tparent_id\texternal_id\tname\trank\tlevel\trep_genome\n")
            fh.write("1\t1\t\td__Bacteria\tdomain\t0\t\n10\t1\t\tg__G\tgenus\t5\t\n20\t1\t\tg__H\tgenus\t5\t\n")
            for t in (11, 12, 13, 14, 15):
                fh.write(f"{t}\t10\t\ts__G T{t}\tspecies\t6\t{genome(t, 0)}\n")
            fh.write(f"16\t10\t\ts__G X\tspecies\t6\t{genome(16, 0)}\n21\t20\t\ts__H T21\tspecies\t6\t{genome(21, 0)}\n")
        with open(os.path.join(d, "heldout_species.txt"), "w") as fh:
            fh.write("s__G X\tspecies\n")

        def own(t, g, lo, k=0):
            """A strain's read of species t: genome k's copy (the representative's for 0) with two private differences."""
            seq = allele[(t, g, k)] if k else copy[(t, g)]
            return mutate(seq, [lo + 17, lo + 101], rng)[lo:lo + 150]

        def novel(g, lo):
            return copy[(16, g)][lo:lo + 150]

        shared = sam_line("r1", "11_5", 100, own(11, 5, 100, 2), "s__G T11", "FN:11", genome(11, 2))
        cls.sams = {"pe": os.path.join(d, "error_reads", "pe"), "se": os.path.join(d, "error_reads", "se")}
        files = {
            ("pe", "s1.FN.sam"): [shared, sam_line("r2", "11_6", 200, own(11, 6, 200), "s__G T11", "FN:11", genome(11, 9))],
            ("pe", "s1.FP.sam"): [shared, sam_line("r3", "11_5", 100, novel(5, 100), X, "FP:11"),
                                  sam_line("r4", "13_6", 300, novel(6, 300), X, "FP:13"),
                                  sam_line("r5", "11_5", 100, novel(5, 100), X, "FP:11", flag=256)],
            ("pe", "s2.FN.sam"): [sam_line("r6", "14_5", 50, own(14, 5, 50), "s__G T14", "FN:14", genome(14, 1))],
            ("pe", "s2.FP.sam"): [sam_line("r7", "14_5", 50, novel(5, 50), X, "FP:14")],
            ("se", "s1.FN.sam"): [shared, sam_line("r8", "12_5", 120, own(12, 5, 120, 1), "s__G T12", "FN:12", genome(12, 1)),
                                  sam_line("r9", "15_6", 400, own(15, 6, 400), "s__G T15", "FN:15", genome(15, 1))],
            ("se", "s1.FP.sam"): [sam_line("r10", "12_5", 120, novel(5, 120), X, "FP:12"),
                                  sam_line("r11", "11_6", 200, novel(6, 200), X, "FP:11"),
                                  sam_line("r12", "21_5", 10, mutate(copy[(21, 5)], [30, 60, 90], rng)[10:160], X, "FP:21")],
            ("se", "s3.FP.sam"): [sam_line("r13", "12_5", 300, novel(5, 300), X, "FP:12", mapq=2),
                                  sam_line("r14", "12_5", 300, novel(5, 300), X, "FP:12")],
        }
        for (kind, name), lines in files.items():
            point = os.path.join(cls.sams[kind], "training", "point")
            os.makedirs(point, exist_ok=True)
            with open(os.path.join(point, name), "w") as fh:
                fh.write("@HD\tVN:1.6\n" + "".join(lines))
        cls.common = ["--reference", os.path.join(d, "reference.fna"), "--taxonomy", os.path.join(d, "internal_taxonomy.dmp"),
                      "--heldout", os.path.join(d, "heldout_species.txt"), "--max-congeners", "2"]
        cls.full = os.path.join(d, "full_reference.fna")

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def compare(self, name, extra):
        """Runs on each folder alone and on both (one after the other, and side by side): the same bytes at each
        prefix, and each folder's lines on stdout after a line naming its prefix."""
        d = os.path.join(self.tmp.name, name)
        kinds = ("pe", "se")
        alone = {}
        for kind in kinds:
            alone[kind] = run(["--sams", self.sams[kind], "--out", os.path.join(d, "alone", kind)] + self.common + extra)
            self.assertNotIn("==", alone[kind])
        for how, threads in (("serial", "1"), ("forked", "2")):
            prefix = {kind: os.path.join(d, how, kind) for kind in kinds}
            both = run(["--sams"] + [self.sams[k] for k in kinds] + ["--out"] + [prefix[k] for k in kinds] +
                       ["--threads", threads] + self.common + extra)
            starts = []
            for kind in kinds:
                for suffix in OUTPUTS:
                    with open(os.path.join(d, "alone", kind) + suffix, "rb") as fh:
                        expected = fh.read()
                    with open(prefix[kind] + suffix, "rb") as fh:
                        self.assertEqual(fh.read(), expected, f"{name}, {how}: {kind}{suffix}")
                block = f"== {prefix[kind]} (the SAMs of {self.sams[kind]}) ==\n" + alone[kind]
                self.assertIn(block, both, f"{name}, {how}: {kind}'s lines")
                starts.append(both.index(block))
            self.assertEqual(starts, sorted(starts))
        # Not vacuous: every counted record of each folder (the record both folders hold counted in each, the secondary
        # and the low-MAPQ ones in neither), both kinds of taxa, and the taxon without a congener in the second.
        self.assertEqual([r["qname"] for r in rows_of(os.path.join(d, "alone", "pe.fragments.tsv.gz"))],
                         ["r1", "r2", "r3", "r4", "r6", "r7"])
        self.assertEqual([r["qname"] for r in rows_of(os.path.join(d, "alone", "se.fragments.tsv.gz"))],
                         ["r1", "r8", "r9", "r10", "r11", "r12", "r14"])
        self.assertIn("taxa without a congener in the training database: 0", alone["pe"])
        self.assertIn("taxa without a congener in the training database: 1", alone["se"])
        self.assertIn("N >=  3:", alone["pe"])
        return alone

    def test_congener_sites_only(self):
        alone = self.compare("congeners", [])
        self.assertNotIn("alleles:", alone["pe"])

    def test_with_the_full_reference(self):
        alone = self.compare("alleles", ["--full-reference", self.full, "--max-alleles", "2"])
        # The folder's own taxon-genes only: 11_5 (its first two of three), 11_6 (2), 13_6 (1) and 14_5 (none), not
        # 12_5 (the other folder's, a congener's copy here).
        self.assertIn("alleles: 5 copies of 4 taxon-genes", alone["pe"])
        self.assertIn("own source genome's allele was left out: FN own 1", alone["pe"])
        # 11_5 (2 of 3), 11_6 (2), 12_5 (2 of 3), 21_5 (2); 15_6 none.
        self.assertIn("alleles: 8 copies of 5 taxon-genes", alone["se"])
        self.assertIn("own source genome's allele was left out: FN own 2", alone["se"])

    def test_with_a_share_of_the_genomes(self):
        self.compare("share", ["--full-reference", self.full, "--allele-genome-share", "0.5"])

    def test_sams_and_out_differ_in_count(self):
        err = io.StringIO()
        with contextlib.redirect_stderr(err), self.assertRaises(SystemExit) as raised:
            an.main(["--sams", self.sams["pe"], self.sams["se"], "--out", os.path.join(self.tmp.name, "one")] + self.common)
        self.assertEqual(raised.exception.code, 2)
        self.assertIn("--sams and --out take a value each per SAM folder, paired in order: 2 SAM folders, 1 output prefix",
                      err.getvalue())
        self.assertFalse(os.path.exists(os.path.join(self.tmp.name, "one.summary.txt")))
        with contextlib.redirect_stderr(err), self.assertRaises(SystemExit):
            an.main(["--sams", self.sams["pe"], self.sams["se"], "--out", "x", "./x"] + self.common)
        self.assertIn("--out names a prefix twice: ./x, x", err.getvalue())


if __name__ == "__main__":
    unittest.main()
