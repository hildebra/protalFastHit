#!/usr/bin/env python3
"""test_mini_db.py - checks for simulate_gtdb_release.py, gtdb_to_protal_db.py, simulate_reads.py and
gtdb_like_lineages.py, the scripts that make the synthetic GTDB release and the mini database.

Runs the scripts into a temporary directory and checks the invariants protal relies on: reference.map byte offsets, a
well-formed internal taxonomy whose leaves are the reference taxids, byte-identical output for a fixed seed and any
number of threads, that sequence similarity follows the taxonomy, and that simulated reads come from the right genomes
in the right proportions. Needs numpy; no protal binary.

The scripts' other tests, beside this file (mini_db_fixtures.py holds what they share):
  test_downloads.py        download_gtdb.py, download_gtdb_releases.py (stand-ins for GTDB's and NCBI's servers)
  test_collector.py        collect_training_data.py, scenarios.py, hifi_reads.py
  test_gene_neighbours.py  gene_neighbours.py
  test_gtdb_build.py       build_gtdb_database.py, build_gtdb_releases.py, rank_genes.py, the reduced database; its
                           end-to-end build needs $PROTAL, $SIMULATE and scikit-learn

  python3 -m unittest scripts/mini_db/test_mini_db.py
"""

import collections
import csv
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mini_db_fixtures import (CONVERT, HERE, LINEAGES, SIMULATE, SIMULATE_READS, check_reference_map,  # noqa: E402
                              full_reference, kmers, read_fasta, run, same_files, shared_release)


class MiniDbTest(unittest.TestCase):
    """The default synthetic release and its conversion (shared_release), and what the scripts make of it with other
    options."""

    @classmethod
    def setUpClass(cls):
        cls.gtdb, cls.db = shared_release()
        cls.tmp = tempfile.TemporaryDirectory()
        cls.conversions = {("gene", 1): cls.db}

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def converted(self, order, threads):
        """The release converted with --order and -t, once for the class: the default (gene, 1) is the shared one."""
        if (order, threads) not in self.conversions:
            folder = os.path.join(self.tmp.name, f"db_{order}_t{threads}")
            run(CONVERT, "--gtdb", self.gtdb, "--outdir", folder, "--order", order, "-t", str(threads))
            self.conversions[(order, threads)] = folder
        return self.conversions[(order, threads)]

    def taxonomy(self):
        with open(os.path.join(self.db, "internal_taxonomy.dmp")) as fh:
            next(fh)
            return {int(r[0]): r for r in (l.rstrip("\n").split("\t") for l in fh)}

    def test_reference_map_offsets(self):
        check_reference_map(self, self.db)

    def test_reference_order(self):
        def map_rows(db):
            with open(os.path.join(db, "reference.map")) as fh:
                return [tuple(map(int, l.split("\t")[:2])) for l in fh]
        rows = map_rows(self.db)  # default: by gene, then taxid
        self.assertEqual(rows, sorted(rows, key=lambda r: (r[1], r[0])))
        genome_rows = map_rows(self.converted("genome", 1))
        self.assertEqual(genome_rows, sorted(genome_rows))
        self.assertEqual(sorted(rows), sorted(genome_rows), "the same genes, only the order differs")
        check_reference_map(self, self.converted("genome", 1))

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
        with open(os.path.join(self.db, "reference.fna")) as r:
            self.assertGreater(full_reference(self.db).count(b">"), 2 * r.read().count(">"))
        # zstd-compressed when the zstd command is there, as protal --build reads it.
        name = "full_reference.fna.zst" if shutil.which("zstd") else "full_reference.fna"
        self.assertEqual(sorted(f for f in os.listdir(self.db) if f.startswith("full_reference")), [name])

    def test_deterministic(self):
        other = os.path.join(self.tmp.name, "gtdb_again")
        run(SIMULATE, "--outdir", other, "--genome_length", "20000")
        genome = os.path.join("genomic_files_reps", "gtdb_genomes_reps_r226", "database", "GCF", "999",
                              "001", "001", "GCF_999001001.1_genomic.fna.gz")
        for f in (os.path.join("genomic_files_all", "bac120_marker_genes_all_r226", "fna", "bac120_TIGR02013.fna"),
                  "bac120_metadata_r226.tsv.gz", genome):
            with open(os.path.join(self.gtdb, f), "rb") as a, open(os.path.join(other, f), "rb") as b:
                self.assertEqual(a.read(), b.read(), f"{f} differs between runs with the same seed")

    def test_parallel_conversion(self):
        # The marker files are read by -t workers; the database must not depend on how many.
        for order in ("gene", "genome"):
            one, four = self.converted(order, 1), self.converted(order, 4)
            same_files(self, one, four, "reference.fna", "reference.map", "internal_taxonomy.dmp")
            self.assertEqual(full_reference(one), full_reference(four), f"{order} order, full reference")
            self.assertFalse(os.path.exists(os.path.join(four, ".convert_tmp")))

    def test_exclude_species(self):
        excluded = os.path.join(self.tmp.name, "excluded.txt")
        with open(excluded, "w") as fh:
            fh.write("# held out\nMockella beta\n")
        direct, copied = os.path.join(self.tmp.name, "db_direct"), os.path.join(self.tmp.name, "db_copied")
        run(CONVERT, "--gtdb", self.gtdb, "--outdir", direct, "--exclude_species", excluded)
        run(CONVERT, "--from_db", self.db, "--exclude_species", excluded, "--outdir", copied)
        taxids = {r[3]: r[0] for r in self.taxonomy().values() if r[4] == "species"}
        for f in ("reference.fna", "reference.map", "internal_taxonomy.dmp"):
            with open(os.path.join(direct, f), "rb") as a, open(os.path.join(copied, f), "rb") as b:
                self.assertEqual(a.read(), b.read(), f"{f}: --from_db and --gtdb differ")
        self.assertEqual(full_reference(direct), full_reference(copied), "full reference: --from_db and --gtdb differ")
        self.assertLess(len(full_reference(direct)), len(full_reference(self.db)))
        with open(os.path.join(self.db, "internal_taxonomy.dmp"), "rb") as a, \
                open(os.path.join(direct, "internal_taxonomy.dmp"), "rb") as b:
            self.assertEqual(a.read(), b.read(), "the taxonomy must keep the species left out")
        with open(os.path.join(direct, "reference.fna")) as fh:
            kept = {line[1:].split("_")[0] for line in fh if line.startswith(">")}
        self.assertNotIn(taxids["s__Mockella beta"], kept)
        self.assertEqual(len(kept), len(taxids) - 1)
        check_reference_map(self, direct)  # reference.map still points at each sequence line
        with self.assertRaises(subprocess.CalledProcessError):
            with open(excluded, "w") as fh:
                fh.write("s__Nonexistent species\n")
            run(CONVERT, "--from_db", self.db, "--exclude_species", excluded, "--outdir", os.path.join(self.tmp.name, "x"))

    def test_a_copy_takes_only_the_converted_files(self):
        # What a build of the folder wrote, or left when stopped, is not copied (its unique_kmers.tsv, of other
        # genes, stopped the copy's build), and what an earlier build of the copy left is removed.
        src, dst = os.path.join(self.tmp.name, "db_built"), os.path.join(self.tmp.name, "db_training")
        shutil.copytree(self.db, src)
        for name in ("unique_kmers.tsv", "gene_conservation.tsv", "index.prx.zst.partial", "database.protal.partial",
                     "database.protal", "build_metadata.tsv"):
            with open(os.path.join(src, name), "w") as fh:
                fh.write("left by a build\n")
        os.makedirs(dst)
        for name in ("unique_kmers.tsv", "index.prx.zst", "database.protal.partial", "model_ONT.xml"):
            with open(os.path.join(dst, name), "w") as fh:
                fh.write("left by an earlier build\n")
        excluded = os.path.join(self.tmp.name, "excluded_copy.txt")
        with open(excluded, "w") as fh:
            fh.write("Mockella beta\n")
        run(CONVERT, "--from_db", src, "--exclude_species", excluded, "--outdir", dst)
        from gtdb_to_protal_db import CONVERTED_FILES, full_reference_path
        expected = {"reference.fna", "reference.map", os.path.basename(full_reference_path(self.db))} | \
            {f for f in CONVERTED_FILES if os.path.isfile(os.path.join(self.db, f))}
        self.assertEqual(set(os.listdir(dst)), expected)
        # Converting a release anew removes the build outputs of the folder's earlier reference, too.
        run(CONVERT, "--gtdb", self.gtdb, "--outdir", src)
        for name in ("unique_kmers.tsv", "gene_conservation.tsv", "index.prx.zst.partial", "database.protal.partial",
                     "database.protal"):
            self.assertFalse(os.path.exists(os.path.join(src, name)), name)

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

    def test_gene_rates(self):
        # With --gene_rates categories, markers evolve at their category's speed (mean 1): congeneric species
        # stay closer at the ribosomal proteins than at the rest. Without it, nothing is written.
        self.assertFalse(os.path.exists(os.path.join(self.gtdb, "simulation", "gene_rates.tsv")))
        other = os.path.join(self.tmp.name, "gtdb_gene_rates")
        run(SIMULATE, "--outdir", other, "--genome_length", "20000", "--gene_rates", "categories")
        with open(os.path.join(other, "simulation", "gene_rates.tsv")) as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        rate = {r["marker"]: float(r["rate"]) for r in rows}
        self.assertAlmostEqual(sum(rate.values()) / len(rate), 1.0, places=3)
        by_category = collections.defaultdict(list)
        for r in rows:
            by_category[r["category"]].append(float(r["rate"]))
        means = {c: sum(v) / len(v) for c, v in by_category.items()}
        self.assertLess(means["ribosomal"], means["translation"])
        self.assertLess(means["translation"], means["other"])

        def similarity(category):  # Mockella alpha vs beta, over the category's bac120 markers
            shared = []
            for r in rows:
                if r["category"] != category:
                    continue
                path = os.path.join(other, "genomic_files_reps", "bac120_marker_genes_reps_r226", "fna",
                                    f"bac120_{r['marker']}.fna")
                if os.path.exists(path):
                    seqs = read_fasta(path)
                    a, b = seqs.get("RS_GCF_999001001.1"), seqs.get("RS_GCF_999002001.1")
                    if a and b:
                        ka, kb = kmers(a), kmers(b)
                        shared.append(len(ka & kb) / len(ka | kb))
            return sum(shared) / len(shared)

        self.assertGreater(similarity("ribosomal"), similarity("other"))

    def test_gene_rates_r226(self):
        # With --gene_rates r226, every gene evolves at the real r226 gene's speed (gene_rates_r226.tsv): strains at
        # its within-species factor, the lineage at its between-congener factor, each set's mean 1, archaea included.
        other = os.path.join(self.tmp.name, "gtdb_gene_rates_r226")
        run(SIMULATE, "--outdir", other, "--genome_length", "20000", "--gene_rates", "r226")
        with open(os.path.join(other, "simulation", "gene_rates.tsv")) as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        # A rate for every marker of both sets (some genes are in both), each once.
        released = set()
        for mset in ("bac120", "ar53"):
            folder = os.path.join(other, "genomic_files_all", f"{mset}_marker_genes_all_r226", "fna")
            released |= {(mset, name[len(mset) + 1:-len(".fna")]) for name in os.listdir(folder)}
        self.assertEqual(sorted((r["set"], r["marker"]) for r in rows), sorted(released))
        for mset in ("bac120", "ar53"):
            for column in ("strain_rate", "branch_rate"):
                values = [float(r[column]) for r in rows if r["set"] == mset]
                self.assertAlmostEqual(sum(values) / len(values), 1.0, places=3)
        strain = {(r["set"], r["name"]): float(r["strain_rate"]) for r in rows}
        # A ribosomal protein diverges slower than a tRNA synthetase, in bacteria and in archaea (uS9, which
        # --gene_rates categories took for a fast gene by its name).
        self.assertLess(strain[("bac120", "Ribosomal_S8")], strain[("bac120", "leuS_bact")])
        self.assertLess(strain[("ar53", "uS9_arch")], 1.0)
        # The table's gene ids follow the converter's order: the same as this world's database.
        import make_gene_rates
        with open(os.path.join(self.db, "gene2geneid.tsv")) as fh:
            ids = {f[0]: int(f[1]) for f in (line.split() for line in fh)}
        self.assertEqual(make_gene_rates.converter_ids(make_gene_rates.read_markers(
            os.path.join(HERE, "markers_r226.tsv"))), ids)
        # A table without a marker of the set stops the simulation.
        short = os.path.join(self.tmp.name, "short_rates.tsv")
        with open(os.path.join(HERE, "gene_rates_r226.tsv")) as fin, open(short, "w") as fout:
            fout.writelines(line for line in fin if "PF00380.20" not in line)
        result = subprocess.run([sys.executable, SIMULATE, "--outdir", os.path.join(self.tmp.name, "gtdb_short"),
                                 "--gene_rates", "r226", "--gene_rate_table", short], capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("lacks 1 markers", result.stderr)

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


class GtdbLikeLineagesTest(unittest.TestCase):
    """gtdb_like_lineages.py: lineages shaped like GTDB's, which simulate_gtdb_release.py --lineages takes."""

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
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "lineages.txt")
            with open(path, "w") as fh:
                fh.write("\n".join(";".join(l) for l in lineages[:20]) + "\n")
            out = os.path.join(tmp, "gtdb_like")
            run(SIMULATE, "--outdir", out, "--lineages", path, "--genome_length", "5000", "--genomes_per_species", "1")
            simulated = set()
            for mset in ("bac120", "ar53"):
                with open(os.path.join(out, f"{mset}_taxonomy_r226.tsv")) as fh:
                    simulated |= {line.rstrip("\n").split("\t")[1] for line in fh}
            self.assertEqual(simulated, {";".join(l) for l in lineages[:20]})


if __name__ == "__main__":
    unittest.main()
