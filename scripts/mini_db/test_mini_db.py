#!/usr/bin/env python3
"""test_mini_db.py - checks for simulate_gtdb_release.py + gtdb_to_protal_db.py.

Runs both scripts into a temporary directory and checks the invariants protal
relies on: reference.map byte offsets, a well-formed internal taxonomy whose
leaves are the reference taxids, determinism, and that sequence similarity
follows the taxonomy. Does not need a protal binary.

  python3 -m unittest scripts/mini_db/test_mini_db.py
"""

import gzip
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SIMULATE = os.path.join(HERE, "simulate_gtdb_release.py")
CONVERT = os.path.join(HERE, "gtdb_to_protal_db.py")


def run(*args):
    subprocess.run([sys.executable, *args], check=True, capture_output=True, text=True)


def read_fasta(path):
    records, header = {}, None
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                header = line[1:].split()[0]
                records[header] = ""
            else:
                records[header] += line
    return records


def kmers(seq, k=21):
    return {seq[i:i + k] for i in range(len(seq) - k + 1)}


class MiniDbTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.gtdb = os.path.join(cls.tmp.name, "gtdb")
        cls.db = os.path.join(cls.tmp.name, "db")
        run(SIMULATE, "--outdir", cls.gtdb, "--genome_length", "20000")
        run(CONVERT, "--gtdb", cls.gtdb, "--outdir", cls.db)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def taxonomy(self):
        with open(os.path.join(self.db, "internal_taxonomy.dmp")) as fh:
            next(fh)
            return {int(r[0]): r for r in (l.rstrip("\n").split("\t") for l in fh)}

    def test_reference_map_offsets(self):
        with open(os.path.join(self.db, "reference.fna"), "rb") as fh:
            data = fh.read()
        with open(os.path.join(self.db, "reference.map")) as fh:
            rows = [list(map(int, l.split("\t"))) for l in fh]
        self.assertEqual(len(rows), data.count(b">"))
        for taxid, geneid, start, end in rows:
            header_start = data.rindex(b">", 0, start)
            self.assertEqual(data[header_start:start], f">{taxid}_{geneid}\n".encode())
            self.assertEqual(data[end:end + 1], b"\n")
            self.assertRegex(data[start:end].decode(), r"^[ACGT]+$")

    def test_taxonomy_tree(self):
        tax = self.taxonomy()
        roots = [i for i, r in tax.items() if int(r[1]) == i]
        self.assertEqual(len(roots), 1)
        names = [r[3] for r in tax.values()]
        self.assertEqual(len(names), len(set(names)), "names must be unique (string_to_id)")
        for i, r in tax.items():
            if i != roots[0]:
                self.assertIn(int(r[1]), tax)
                self.assertEqual(int(r[5]), int(tax[int(r[1])][5]) + 1, f"level of {r[3]}")
        with open(os.path.join(self.db, "reference.map")) as fh:
            ref_taxids = {int(l.split("\t")[0]) for l in fh}
        parents = {int(r[1]) for i, r in tax.items() if i != roots[0]}
        for t in ref_taxids:
            self.assertEqual(tax[t][4], "species")
            self.assertTrue(tax[t][3].startswith("s__"))
            self.assertRegex(tax[t][6], r"^GCF_999\d{6}\.1$")
            self.assertNotIn(t, parents, "reference taxids must be leaves")

    def test_full_reference_covers_all_genomes(self):
        with open(os.path.join(self.gtdb, "simulation", "genomes.tsv")) as fh:
            n_genomes = sum(1 for _ in fh) - 1
        with open(os.path.join(self.db, "genome2tiid.tsv")) as fh:
            self.assertEqual(sum(1 for _ in fh), n_genomes)
        with open(os.path.join(self.db, "full_reference.fna")) as f, \
             open(os.path.join(self.db, "reference.fna")) as r:
            self.assertGreater(f.read().count(">"), 2 * r.read().count(">"))

    def test_deterministic(self):
        other = os.path.join(self.tmp.name, "gtdb_again")
        run(SIMULATE, "--outdir", other, "--genome_length", "20000")
        marker = os.path.join("genomic_files_all", "bac120_marker_genes_all_r226", "fna", "bac120_TIGR02013.fna")
        with open(os.path.join(self.gtdb, marker)) as a, open(os.path.join(other, marker)) as b:
            self.assertEqual(a.read(), b.read())

    def test_similarity_follows_taxonomy(self):
        seqs = read_fasta(os.path.join(self.gtdb, "genomic_files_all", "bac120_marker_genes_all_r226",
                                       "fna", "bac120_TIGR02013.fna"))

        def shared(a, b):
            ka, kb = kmers(seqs[a]), kmers(seqs[b])
            return len(ka & kb) / len(ka | kb)

        strain = shared("RS_GCF_999001001.1", "GB_GCA_999001002.1")      # alpha vs alpha
        congeneric = shared("RS_GCF_999001001.1", "RS_GCF_999002001.1")  # alpha vs beta
        phylum = shared("RS_GCF_999001001.1", "RS_GCF_999003001.1")      # alpha vs gamma
        self.assertGreater(strain, congeneric)
        self.assertGreater(congeneric, phylum)

    def test_no_internal_stops(self):
        for rec, seq in read_fasta(os.path.join(self.gtdb, "genomic_files_reps", "bac120_marker_genes_reps_r226",
                                                "faa", "bac120_TIGR02013.faa")).items():
            self.assertNotIn("*", seq, rec)


if __name__ == "__main__":
    unittest.main()
