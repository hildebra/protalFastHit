#!/usr/bin/env python3
"""Tests of insilico_strains.py on a small synthetic genome: divergence reached, codon awareness, the table.

Run: python3 scripts/test_insilico_strains.py
"""

import gzip
import os
import re
import sys
import tempfile
import unittest

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import insilico_strains as ins  # noqa: E402

ACCESSION = re.compile(r"(?:RS_|GB_)?(GC[AF]_\d{9}\.\d+)")  # build_gtdb_database.py's


def coding_sequence(rng, codons):
    """A reading frame without stop codons (ATG first)."""
    sense = [i for i in range(64) if not ins.STOP[i]]
    picks = rng.choice(sense, size=codons)
    letters = "ACGT"
    return "ATG" + "".join(letters[c // 16] + letters[(c // 4) % 4] + letters[c % 4] for c in picks)


def reverse_complement(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


class InsilicoStrains(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        rng = np.random.default_rng(3)
        cls.tmp = tempfile.TemporaryDirectory()
        d = cls.tmp.name
        # One species of one genome: genes on both strands between random spacers; another species of two genomes.
        parts, cls.genes = [], []
        pos = 0
        for k in range(40):
            spacer = "".join(rng.choice(list("ACGT"), size=150))
            gene = coding_sequence(rng, 300)
            strand = "+" if k % 2 == 0 else "-"
            parts += [spacer, gene if strand == "+" else reverse_complement(gene)]
            start = pos + len(spacer) + 1
            cls.genes.append((start, start + len(gene) - 1, strand, k + 1))
            pos += len(spacer) + len(gene)
        cls.genome = "".join(parts)
        cls.single = os.path.join(d, "GCF_000000001.1_genomic.fna.gz")
        with gzip.open(cls.single, "wt") as fh:
            fh.write(">contig1 test\n" + "\n".join(cls.genome[i:i + 70] for i in range(0, len(cls.genome), 70)) + "\n")
        other = os.path.join(d, "other.fna")
        with open(other, "w") as fh:
            fh.write(">c\nACGTACGTACGT\n")
        cls.table = os.path.join(d, "genomes.tsv")
        with open(cls.table, "w") as fh:
            fh.write(f"GCF_000000001.1\td__Bacteria;s__Solo one\t{cls.single}\t{len(cls.genome)}\n")
            fh.write(f"GCF_000000002.1\td__Bacteria;s__Duo two\t{other}\t12\n")
            fh.write(f"GCF_000000003.1\td__Bacteria;s__Duo two\t{other}\t12\n")
        cls.positions = os.path.join(d, "gene_positions.tsv")
        with open(cls.positions, "w") as fh:
            fh.write("accession\ttaxid\tcontig\tcontig_length\tcircular\tgene\tstart\tend\tstrand\tplaced\tkmer_share\n")
            for s, e, strand, gene in cls.genes:
                fh.write(f"GCF_000000001.1\t1\tcontig1\t{len(cls.genome)}\t1\t{gene}\t{s}\t{e}\t{strand}\texact\t1\n")
        cls.out = os.path.join(d, "strains")
        cls.output = os.path.join(d, "genomes_simulated.tsv")
        ins.main(["--genome-table", cls.table, "--output", cls.output, "--out-dir", cls.out, "--positions",
                  cls.positions, "--ani", "96-96", "--marker-scale", "1", "-t", "1"])
        name = ins.strain_name("GCF_000000001.1")
        cls.strain = ins.read_fasta(os.path.join(cls.out, name + ".fna.gz"))[0][1].decode()

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_table(self):
        rows = ins.read_table(self.output)
        self.assertEqual(len(rows), 4)  # the three genomes and one strain: only the one-genome species gets one
        name = rows[-1][0]
        self.assertTrue(ins.is_insilico(name))
        self.assertIsNone(ACCESSION.search(name), "an accession pattern would take the strain for the representative")
        self.assertEqual(rows[-1][1], "d__Bacteria;s__Solo one")
        self.assertEqual(int(rows[-1][3]), len(self.genome))

    def test_contig_names(self):
        # Named after the strain, so that its reads (named by their contig) are told from the
        # representative's (trace_relatives.py).
        name = ins.strain_name("GCF_000000001.1")
        header = ins.read_fasta(os.path.join(self.out, name + ".fna.gz"))[0][0]
        self.assertTrue(header.startswith(name + "_contig1 "), header)

    def test_divergence(self):
        self.assertEqual(len(self.strain), len(self.genome))
        diff = np.frombuffer(self.strain.encode(), np.uint8) != np.frombuffer(self.genome.encode(), np.uint8)
        # Genes at the marker rate (scale 1: as the genome), the spacers at 1.3x a coding rate around it.
        self.assertAlmostEqual(diff.mean(), 0.04, delta=0.01)
        in_genes = np.zeros(len(self.genome), bool)
        for s, e, _, _ in self.genes:
            in_genes[s - 1:e] = True
        self.assertAlmostEqual(diff[in_genes].mean(), 0.04, delta=0.008)

    def test_codons(self):
        third = changed = new_stops = 0
        for s, e, strand, _ in self.genes:
            ref = self.genome[s - 1:e]
            alt = self.strain[s - 1:e]
            if strand == "-":
                ref, alt = reverse_complement(ref), reverse_complement(alt)
            for i in range(0, len(ref) - 2, 3):
                a, b = ref[i:i + 3], alt[i:i + 3]
                if a != b:
                    changed += 1
                    third += a[:2] == b[:2]
                    index = 16 * "ACGT".index(b[0]) + 4 * "ACGT".index(b[1]) + "ACGT".index(b[2])
                    new_stops += bool(ins.STOP[index])
        self.assertGreater(changed, 300)
        self.assertEqual(new_stops, 0)
        self.assertGreater(third / changed, 0.55)  # most differences on third positions, as a strain's

    def test_factors(self):
        strains = {f"s{k}": {g: 0.01 * (2.0 if g == 1 else 0.5 if g == 2 else 1.0) * (1 + k) for g in range(1, 21)}
                   for k in range(5)}
        factors = ins.conservation_factors(strains)
        self.assertAlmostEqual(factors[1], 2.0)
        self.assertAlmostEqual(factors[2], 0.5)
        self.assertAlmostEqual(factors[3], 1.0)
        self.assertEqual(len(ins.strain_divergences(strains)), 5)

    def test_codons_of_both_strands_at_one_base(self):
        # Where a minus-strand frame (sites 0-299, its first codon read from base 300) meets a plus-strand frame
        # (from base 300), both have a codon starting at base 300: two substitutions there were summed into one
        # codon index past 63 (the r226 v11 build's IndexError) and must be judged per codon.
        rng = np.random.default_rng(5)
        n = 900
        owner = np.where(np.arange(n) < 300, 1, 0)
        minus = owner == 1
        cpos = np.where(minus, (300 - np.arange(n)) % 3, (np.arange(n) - 300) % 3).astype(np.int8)
        for seed in range(40):
            codes = rng.integers(0, 4, n).astype(np.uint8)
            new, *_ = ins.mutate(codes, owner, cpos, minus, np.full(n, 0.6), 2, np.random.default_rng(seed), 3.0, 0.15)
            plus = new[300:900].astype(np.int64).reshape(-1, 3)
            index = 16 * plus[:, 0] + 4 * plus[:, 1] + plus[:, 2]
            before = codes[300:900].astype(np.int64).reshape(-1, 3)
            was_stop = ins.STOP[16 * before[:, 0] + 4 * before[:, 1] + before[:, 2]]
            self.assertFalse((ins.STOP[index] & ~was_stop).any(), seed)

    def test_substitutions_without_alignment(self):
        rng = np.random.default_rng(11)
        rep = rng.integers(0, 4, 900).astype(np.uint8)
        other = rep.copy()
        changed = np.sort(rng.choice(900, 30, replace=False))
        other[changed] = (other[changed] + 1 + rng.integers(0, 3, 30)) % 4
        positions, compared = ins.substitutions(rep, other)
        found, made = set(positions.tolist()), set(changed.tolist())
        self.assertTrue(found <= made)
        # Every substitution between the first and the last shared 12-mer; the ends past them are not compared.
        self.assertTrue({p for p in made if 24 <= p < 876} <= found)
        self.assertGreater(compared, 800)
        # An insertion in the other copy: the stretch after it leaves the main diagonal and is not compared, the
        # substitutions before it are still found.
        with_indel = np.concatenate([other[:450], rng.integers(0, 4, 7).astype(np.uint8), other[450:]])
        positions, compared = ins.substitutions(rep, with_indel)
        before = {p for p in changed.tolist() if p < 440}
        self.assertTrue(before <= set(positions.tolist()))
        self.assertLess(compared, 900)

    def test_spectrum_and_omega(self):
        # The in-silico strain of setUpClass was made at OMEGA_DEFAULT (no real strains to calibrate on): measured
        # against its representative on the placed genes, most of its substitutions are on third positions, and the
        # calibration recovers an omega near the one used; a more synonymous spectrum needs a smaller omega.
        name = ins.strain_name("GCF_000000001.1")
        strain_path = os.path.join(self.out, name + ".fna.gz")
        genes = [("contig1", gene, s, e, strand) for s, e, strand, gene in self.genes]
        third, subs, compared, copies = ins.spectrum_pair((strain_path, self.single, [(name + "_contig1",) + g[1:]
                                                                                       for g in genes], genes))
        self.assertGreater(subs, 300)
        self.assertGreater(compared, 30000)
        share = third / subs
        self.assertGreater(share, 0.55)
        rng = np.random.default_rng(2)
        omega = ins.calibrate_omega(share, copies, 3.0, rng)
        self.assertTrue(0.06 <= omega <= 0.35, omega)
        lower = ins.calibrate_omega(min(0.95, share + 0.1), copies, 3.0, rng)
        self.assertLess(lower, omega)
        self.assertEqual(ins.calibrate_omega(None, copies, 3.0, rng), ins.OMEGA_DEFAULT)

    def test_orfs_both_strands(self):
        codes = ins.CODE[np.frombuffer(self.genome.encode(), np.uint8)]
        found = ins.orfs(codes)
        for s, e, strand, _ in self.genes[:6]:
            self.assertTrue(any(f[2] == strand and f[0] <= s - 1 and f[1] >= e - 3 for f in found), (s, e, strand))


if __name__ == "__main__":
    unittest.main()
