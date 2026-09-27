#!/usr/bin/env python3
"""test_mini_db.py - checks for simulate_gtdb_release.py, gtdb_to_protal_db.py
and simulate_reads.py.

Runs the scripts into a temporary directory and checks the invariants protal
relies on: reference.map byte offsets, a well-formed internal taxonomy whose
leaves are the reference taxids, byte-identical output for a fixed seed, that
sequence similarity follows the taxonomy, and that simulated reads come from
the right genomes in the right proportions. Does not need a protal binary.

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
SIMULATE_READS = os.path.join(HERE, "simulate_reads.py")


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
        genome = os.path.join("genomic_files_reps", "gtdb_genomes_reps_r226", "database", "GCF", "999",
                              "001", "001", "GCF_999001001.1_genomic.fna.gz")
        for f in (os.path.join("genomic_files_all", "bac120_marker_genes_all_r226", "fna", "bac120_TIGR02013.fna"),
                  "bac120_metadata_r226.tsv.gz", genome):
            with open(os.path.join(self.gtdb, f), "rb") as a, open(os.path.join(other, f), "rb") as b:
                self.assertEqual(a.read(), b.read(), f"{f} differs between runs with the same seed")

    def simulate_reads(self, prefix, community, pairs=400, error_rate="0"):
        path = os.path.join(self.tmp.name, "community.tsv")
        with open(path, "w") as fh:
            fh.write("# comment\naccession\trelative_abundance\n")
            fh.writelines(f"{acc}\t{ab}\n" for acc, ab in community.items())
        out = os.path.join(self.tmp.name, prefix)
        run(SIMULATE_READS, "--genomes", os.path.join(self.gtdb, "simulation", "genomes.tsv"),
            "--community", path, "--out_prefix", out, "--pairs", str(pairs), "--error_rate", error_rate)
        return out

    def test_reads(self):
        community = {"GCA_999001002.1": 0.75, "GCA_999003003.1": 0.25}
        out = self.simulate_reads("reads", community)
        with open(out + ".truth.tsv") as fh:
            next(fh)
            truth = {f[0]: f for f in (l.rstrip("\n").split("\t") for l in fh)}
        self.assertEqual(set(truth), set(community))
        self.assertEqual(truth["GCA_999003003.1"][1], "s__Fakibacter gamma")
        pairs = {acc: int(f[5]) for acc, f in truth.items()}
        self.assertEqual(sum(pairs.values()), 400)
        # Read pairs follow abundance x genome length.
        weight = {acc: community[acc] * int(f[4]) for acc, f in truth.items()}
        for acc in community:
            self.assertAlmostEqual(pairs[acc], 400 * weight[acc] / sum(weight.values()), delta=1)

        with open(out + "_R1.fq") as f1, open(out + "_R2.fq") as f2:
            r1, r2 = f1.read().split("\n"), f2.read().split("\n")
        self.assertEqual(r1[0::4], r2[0::4], "mates share a read name")
        genomes = {}
        for acc in community:
            fasta = os.path.join(self.gtdb, "simulation", "genomes_nonreps", f"{acc}_genomic.fna.gz")
            genomes[acc] = "|".join(read_fasta(fasta).values())
        comp = str.maketrans("ACGT", "TGCA")
        for name, s1, s2 in zip(r1[0::4], r1[1::4], r2[1::4]):
            acc = name[1:].rsplit("-", 1)[0]
            self.assertIn(acc, community)
            # Error-free: both mates are genome substrings in FR orientation.
            g = genomes[acc]
            self.assertTrue(s1 in g or s1.translate(comp)[::-1] in g, name)
            self.assertTrue(s2 in g or s2.translate(comp)[::-1] in g, name)

        again = self.simulate_reads("reads_again", community)
        for suffix in ("_R1.fq", "_R2.fq", ".truth.tsv"):
            with open(out + suffix, "rb") as a, open(again + suffix, "rb") as b:
                self.assertEqual(a.read(), b.read(), f"{suffix} differs between runs with the same seed")

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
