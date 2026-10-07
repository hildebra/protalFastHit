#!/usr/bin/env python3
"""Tests of trace_relatives.py: reads traced to their genomes by the contigs in their names, and the script end to end
(where a held-out species' reads land, by the genes' conservation factors).

At GTDB r226 (v10) the trace found no gene: it took a read's genome from the part of its name before "_contig", as the
synthetic genomes of simulate_gtdb_release.py are named, while GTDB's genomes have NCBI's contig names
(docs/claude/2026-10-05-trace-relatives.md). Standard library only.

Run: python3 -m unittest scripts/test_trace_relatives.py
"""

import csv
import gzip
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import trace_relatives  # noqa: E402

LINEAGE = "d__Bacteria;p__P;c__C;o__O;f__F;g__G;s__G {}"
GENES = (1, 2, 3)


class TraceRelatives(unittest.TestCase):
    def world(self, ncbi_names=True, fasta_paths=True):
        """A sample of species A (in the database, taxid 1) and its congener B (held out at species rank), whose
        reads land on A; and C, held out at genus rank, whose reads are not compared. -> trace_relatives options."""
        d = tempfile.mkdtemp(dir=self.tmp.name)
        genomes = {"GCF_000000001.1": ("A", "NZ_CP000001.1"), "GCA_000000002.1": ("B", "JAAAAA010000001.1"),
                   "GCA_000000003.1": ("C", "NZ_CP000003.1")}
        point = os.path.join(d, "points", "rl150_p1000")
        os.makedirs(os.path.join(point, "sim"))
        os.makedirs(os.path.join(point, "protal", "alignments"))
        rows, contig = [], {}
        for acc, (sp, ncbi) in genomes.items():
            name = ncbi if ncbi_names else f"{acc}_contig1"
            contig[sp] = name
            fasta = os.path.join(d, f"{acc}.fna.gz")
            with gzip.open(fasta, "wt") as fh:
                fh.write(f">{name} {sp} chromosome\nACGT\n>{name}_plasmid\nACGT\n")
            rows.append({"sample": "rl150_p1000_s_1", "genome": acc, "taxonomy": LINEAGE.format(sp),
                         "vertical_coverage": "2.0", "fasta_path": fasta})
        columns = ["sample", "genome", "taxonomy", "vertical_coverage"] + (["fasta_path"] if fasta_paths else [])
        with open(os.path.join(point, "sim", "manifest.tsv"), "w") as fh:
            fh.write("\t".join(columns) + "\n")
            fh.writelines("\t".join(r[c] for c in columns) + "\n" for r in rows)
        with open(os.path.join(point, "protal", "alignments", "rl150_p1000_s_1.sam"), "w") as fh:
            fh.write("@HD\tVN:1.6\n")
            n = 0
            for gene in GENES:
                for sp, mapq in (("A", 60), ("A", 60), ("B", 60), ("B", 0), ("C", 60)):
                    n += 1
                    fh.write(f"{contig[sp]}-{n}\t0\t1_{gene}\t100\t{mapq}\t150M\t*\t0\t0\tA\tI\n")
        db = os.path.join(d, "training_db")
        os.makedirs(db)
        with open(os.path.join(db, "genome2tiid.tsv"), "w") as fh:
            fh.write(f"GCF_000000001.1\t1\tGCF_000000001.1\t{LINEAGE.format('A')}\n")
        with open(os.path.join(db, "gene_congeners.tsv"), "w") as fh:
            fh.write("geneid\twithin_factor\n" + "".join(f"{g}\t{0.5 * g}\n" for g in GENES))
        heldout = os.path.join(d, "heldout_species.txt")
        with open(heldout, "w") as fh:
            fh.write("s__G B\tspecies\tG\ns__G C\tgenus\tg__G\n")
        return trace_relatives.parse_args(["--points", os.path.join(d, "points"), "--db", db, "--heldout", heldout,
                                           "--out", os.path.join(d, "out"), "--threads", "2"])

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmp.cleanup()

    def check(self, result, problem):
        self.assertIsNone(problem)
        self.assertEqual(len(result["genes"]), len(GENES))
        for g in result["genes"]:
            self.assertEqual(g["own_records"], 2)
            self.assertEqual(g["relative_records"], 2)  # C's reads, held out at genus rank, are not counted
            self.assertAlmostEqual(g["relative_mapq_below_4"], 0.5)
            self.assertAlmostEqual(g["R"], 1.0)  # the same coverage, the same records per gene

    def test_ncbi_contig_names(self):
        result, problem = trace_relatives.trace(self.world())
        self.check(result, problem)
        self.assertEqual(result["unnamed"], 0)

    def test_synthetic_contig_names_without_fasta_paths(self):
        # An older manifest has no FASTA paths: the genome is the read name's part before "_contig".
        result, problem = trace_relatives.trace(self.world(ncbi_names=False, fasta_paths=False))
        self.check(result, problem)

    def test_nothing_traced_is_said(self):
        opts = self.world()
        os.remove(os.path.join(os.path.dirname(opts.db), "GCF_000000001.1.fna.gz"))
        for acc in ("GCA_000000002.1", "GCA_000000003.1"):
            os.remove(os.path.join(os.path.dirname(opts.db), f"{acc}.fna.gz"))
        result, problem = trace_relatives.trace(opts)
        self.assertIsNone(result)
        self.assertIn("none of their aligned records was traced", problem)

    def test_read_contig(self):
        self.assertEqual(trace_relatives.read_contig("NZ_CP012345.1-4321"), "NZ_CP012345.1")
        self.assertEqual(trace_relatives.read_contig("contig-with-dashes-7"), "contig-with-dashes")

    def test_reads_by_gene_factor(self):
        # The script end to end: where the reads of a species the training database lacks land, by the genes' factors
        # (gene_congeners.tsv), per unit coverage, with an older manifest's synthetic contig names; its table and report.
        tmp = self.tmp.name
        point = os.path.join(tmp, "points", "rl100_p1000")
        os.makedirs(os.path.join(point, "sim"))
        os.makedirs(os.path.join(point, "protal", "alignments"))
        lineage = "d__Bacteria;p__P;c__C;o__O;f__F;g__G;s__G {}"
        with open(os.path.join(point, "sim", "manifest.tsv"), "w") as fh:
            fh.write("sample\tgenome\tspecies\ttaxonomy\tvertical_coverage\n")
            fh.write(f"rl100_p1000_s_1\tGCF_1.1\tG a\t{lineage.format('a')}\t2.0\n")  # in the database (taxon 1)
            fh.write(f"rl100_p1000_s_1\tGCA_2.1\tG b\t{lineage.format('b')}\t1.0\n")  # held out: lands on taxon 3
        db = os.path.join(tmp, "training_db")
        os.makedirs(db)
        with open(os.path.join(db, "genome2tiid.tsv"), "w") as fh:
            for taxid, name in ((1, "a"), (3, "c")):
                fh.write(f"GCF_{taxid}.1\t{taxid}\tGCF_{taxid}.1\t{lineage.format(name)}\n")
        with open(os.path.join(db, "gene_congeners.tsv"), "w") as fh:
            fh.write("geneid\twithin_factor\tbetween_factor\tpairs\tspecies\tidentical_share\tnear_identical_share\n"
                     "1\t0.5\t0.4\t3\t2\t0\t0\n2\t1.5\t1.6\t3\t2\t0\t0\n")
        heldout = os.path.join(tmp, "heldout_species.txt")
        with open(heldout, "w") as fh:
            fh.write("s__G b\tspecies\ts__G b\n")
        records = []
        for i in range(8):  # the species' own: 4 records on each gene
            records.append((f"GCF_1.1_contig1-{i}", f"1_{1 + i % 2}", 60))
        for i, mapq in enumerate((60, 60, 0, 1)):  # the relative: 4 on the conserved gene, 2 of them ambiguous
            records.append((f"GCA_2.1_contig1-{i}", "3_1", mapq))
        records.append(("GCA_2.1_contig1-9", "3_2", 60))  # and 1 on the fast one
        records.append(("GCF_1.1_contig1-20", "3_1", 60))  # the species' own read on a congener: not counted
        with open(os.path.join(point, "protal", "alignments", "rl100_p1000_s_1.sam"), "w") as fh:
            fh.write("@HD\tVN:1.6\n")
            for name, ref, mapq in records:
                fh.write(f"{name}\t0\t{ref}\t1\t{mapq}\t100M\t*\t0\t0\t{'A' * 100}\t{'I' * 100}\n")
            fh.write("GCA_2.1_contig1-5\t256\t1_1\t1\t0\t100M\t*\t0\t0\t*\t*\n")  # secondary: not counted
        out = os.path.join(tmp, "logs", "relatives")
        command = [sys.executable, os.path.join(HERE, "trace_relatives.py"), "--points", os.path.join(tmp, "points"),
                   "--db", db, "--heldout", heldout, "--out", out]
        result = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        with open(out + ".tsv") as fh:
            genes = {r["geneid"]: r for r in csv.DictReader(fh, delimiter="\t")}
        # Per unit coverage: own 4 / 2 on each gene; the relative 4 / 1 on gene 1 (2 kept), 1 / 1 on gene 2.
        self.assertAlmostEqual(float(genes["1"]["R"]), 2.0)
        self.assertAlmostEqual(float(genes["1"]["R_kept"]), 1.0)
        self.assertAlmostEqual(float(genes["2"]["R"]), 0.5)
        self.assertAlmostEqual(float(genes["1"]["relative_mapq_below_4"]), 0.5)
        self.assertAlmostEqual(float(genes["1"]["relative_on_congener"]), 1.0)
        with open(out + ".txt") as fh:
            text = fh.read()
        self.assertIn("1 paired-end training samples", text)
        self.assertIn("| factor < 0.7 | 1 | 2.000 | 1.000 |", text)

        # Without factors there is nothing to trace by.
        os.remove(os.path.join(db, "gene_congeners.tsv"))
        result = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("Not traced: the training database has no gene conservation factors", result.stdout)


if __name__ == "__main__":
    unittest.main()
