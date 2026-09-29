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
LINEAGES = os.path.join(HERE, "gtdb_like_lineages.py")
STRAINS = os.path.join(HERE, "..", "gtdb_strain_genomes.py")


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

    def test_reference_order(self):
        def map_rows(db):
            with open(os.path.join(db, "reference.map")) as fh:
                return [tuple(map(int, l.split("\t")[:2])) for l in fh]
        rows = map_rows(self.db)  # default: by gene, then taxid
        self.assertEqual(rows, sorted(rows, key=lambda r: (r[1], r[0])))
        genome_db = os.path.join(self.tmp.name, "db_genome_order")
        run(CONVERT, "--gtdb", self.gtdb, "--outdir", genome_db, "--order", "genome")
        genome_rows = map_rows(genome_db)
        self.assertEqual(genome_rows, sorted(genome_rows))
        self.assertEqual(sorted(rows), sorted(genome_rows), "the same genes, only the order differs")

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

    def test_exclude_species(self):
        excluded = os.path.join(self.tmp.name, "excluded.txt")
        with open(excluded, "w") as fh:
            fh.write("# held out\nMockella beta\n")
        direct, copied = os.path.join(self.tmp.name, "db_direct"), os.path.join(self.tmp.name, "db_copied")
        run(CONVERT, "--gtdb", self.gtdb, "--outdir", direct, "--exclude_species", excluded)
        run(CONVERT, "--from_db", self.db, "--exclude_species", excluded, "--outdir", copied)
        taxids = {r[3]: r[0] for r in self.taxonomy().values() if r[4] == "species"}
        for f in ("reference.fna", "reference.map", "full_reference.fna", "internal_taxonomy.dmp"):
            with open(os.path.join(direct, f), "rb") as a, open(os.path.join(copied, f), "rb") as b:
                self.assertEqual(a.read(), b.read(), f"{f}: --from_db and --gtdb differ")
        with open(os.path.join(self.db, "internal_taxonomy.dmp"), "rb") as a, \
                open(os.path.join(direct, "internal_taxonomy.dmp"), "rb") as b:
            self.assertEqual(a.read(), b.read(), "the taxonomy must keep the species left out")
        with open(os.path.join(direct, "reference.fna")) as fh:
            kept = {line[1:].split("_")[0] for line in fh if line.startswith(">")}
        self.assertNotIn(taxids["s__Mockella beta"], kept)
        self.assertEqual(len(kept), len(taxids) - 1)
        # reference.map still points at each sequence line
        with open(os.path.join(direct, "reference.fna"), "rb") as fna, open(os.path.join(direct, "reference.map")) as fmap:
            data = fna.read()
            for line in fmap:
                tid, gid, start, end = line.split()
                header_end = int(start)
                self.assertEqual(data[data.rfind(b">", 0, header_end):header_end], f">{tid}_{gid}\n".encode())
                self.assertEqual(data[int(end):int(end) + 1], b"\n")
        with self.assertRaises(subprocess.CalledProcessError):
            with open(excluded, "w") as fh:
                fh.write("s__Nonexistent species\n")
            run(CONVERT, "--from_db", self.db, "--exclude_species", excluded, "--outdir", os.path.join(self.tmp.name, "x"))

    def test_strain_genomes(self):
        out = os.path.join(self.tmp.name, "strains")
        run(STRAINS, "--gtdb", self.gtdb, "-o", out, "--species", "2", "--per_species", "1", "--rep_only_species", "1")
        with open(os.path.join(out, "strain_accessions.txt")) as fh:
            accessions = fh.read().split()
        with open(os.path.join(out, "simulation_species.txt")) as fh:
            pool = fh.read().splitlines()
        self.assertEqual(len(accessions), 2)
        self.assertTrue(all(a.startswith("GCA_999") for a in accessions), "non-representatives, without RS_/GB_")
        self.assertEqual(len(pool), 3)
        self.assertTrue(all(s.startswith("s__") for s in pool))

    def test_gtdb_like_lineages(self):
        text = subprocess.run([sys.executable, LINEAGES, "--species", "300", "--archaea", "0.1", "--seed", "3"],
                              check=True, capture_output=True, text=True).stdout
        lineages = [l.split(";") for l in text.splitlines()]
        self.assertEqual(len(lineages), 300)
        self.assertTrue(all([r[:3] for r in l] == ["d__", "p__", "c__", "o__", "f__", "g__", "s__"] for l in lineages))
        self.assertEqual(sum(l[0] == "d__Archaea" for l in lineages), 30)
        self.assertEqual(len({l[6] for l in lineages}), 300)
        parent = {}
        for l in lineages:  # every name has one parent: names are unique across the tree
            for rank in range(1, 7):
                self.assertEqual(parent.setdefault(l[rank], l[rank - 1]), l[rank - 1], l[rank])
        genera = {}
        for l in lineages:
            genera[l[5]] = genera.get(l[5], 0) + 1
        self.assertGreater(sum(1 for n in genera.values() if n == 1), len(genera) / 3, "most genera have one species")
        self.assertGreater(max(genera.values()), 5)
        path = os.path.join(self.tmp.name, "lineages.txt")
        with open(path, "w") as fh:
            fh.write("\n".join(";".join(l) for l in lineages[:20]) + "\n")
        run(SIMULATE, "--outdir", os.path.join(self.tmp.name, "gtdb_like"), "--lineages", path,
            "--genome_length", "5000", "--genomes_per_species", "1")

    def test_divergence_ranges(self):
        def divergence(root):
            with open(os.path.join(root, "simulation", "divergence.tsv")) as fh:
                next(fh)
                return [(float(r[2]), float(r[3])) for r in (l.rstrip("\n").split("\t") for l in fh)]

        self.assertEqual(set(divergence(self.gtdb)), {(0.035, 0.005)})
        other = os.path.join(self.tmp.name, "gtdb_ranges")
        run(SIMULATE, "--outdir", other, "--genome_length", "20000",
            "--strain_divergence", "0.002-0.02", "--species_divergence", "0.01-0.04")
        rates = divergence(other)
        self.assertTrue(all(0.01 <= s <= 0.04 and 0.002 <= g <= 0.02 for s, g in rates))
        self.assertGreater(len({s for s, _ in rates}), 1)
        self.assertEqual(len({g for _, g in rates}), len(rates))

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
