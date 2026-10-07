#!/usr/bin/env python3
"""test_gene_neighbours.py - checks for gene_neighbours.py: where each marker gene lies in every genome to simulate
from, and the neighbours of its ends counted per clade (gene_neighbours.tsv, gene_positions.tsv), on synthetic releases
whose markers lie in operon-like clusters. Standard library only; no protal binary (tests/e2e has protal's side, its
GeneNeighboursTest).

  python3 -m unittest scripts/mini_db/test_gene_neighbours.py
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
from mini_db_fixtures import (CONVERT, GENE_NEIGHBOURS, SIMULATE, operon_release, read_fasta,  # noqa: E402
                              run)

import gene_neighbours as gn  # noqa: E402


class GeneNeighboursScriptTest(unittest.TestCase):
    """gene_neighbours.py on a release whose markers lie in operon-like clusters (operon_release: simulate_gtdb_release.py
    --operons), 4 species in one family, 2 in another and 2 archaea, two genomes each in 4 contigs: every gene of the
    database is found in every genome of its species where the simulator put it (a representative's exactly, the
    other strain's by its k-mer trace, a few bases off at most), and the lines of each family count exactly the
    neighbours that its species' marker positions give (each species the partner most of its genomes show there,
    worked out here independently); so do the lines derived from the positions alone, and those a training
    database's copy derives without a species."""

    MAX_GAP = 3000
    SLACK = 15  # bases a strain's gene placed by its trace may be off: its codon indels

    @classmethod
    def setUpClass(cls):
        cls.gtdb, cls.db, cls.output = operon_release()
        cls.positions = os.path.join(cls.db, "gene_positions.tsv")
        cls.tmp = tempfile.TemporaryDirectory()

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_genome_index_finds_what_find_finds(self):
        # place() looks genes up in an index of the genome's k-mers at a stride instead of scanning the genome with
        # find: it must find every occurrence find does, overlapping ones and those at contig ends too.
        import random
        rng = random.Random(9)
        contigs = [bytes(rng.choice(b"ACGT") for _ in range(n)) for n in (5000, 2000, 37, 3000)]
        genes = [contigs[0][100:1100], contigs[0][4900 - 63:], contigs[1][:63], contigs[1][5:67], contigs[3][17:3000],
                 contigs[2], b"A" * 40, contigs[0][2000:2062]]
        repeat = contigs[0][300:420]
        contigs[1] = contigs[1][:500] + repeat + repeat[:60] + repeat + contigs[1][500 + 300:]  # repeated, overlapping
        contigs[3] = contigs[3][:1000] + b"A" * 100 + contigs[3][1100:]
        genes += [repeat, repeat[:60] + repeat[:60], b"A" * 70, b"A" * 100]
        genome = gn.SEPARATOR.join(contigs)
        index = gn.GenomeIndex(genome)
        for seq in genes + [s.translate(gn.COMPLEMENT)[::-1] for s in genes] + [b"C" * 80, genes[0][:-1] + b"T"]:
            expected, pos = [], genome.find(seq)
            while pos >= 0:
                expected.append(pos)
                pos = genome.find(seq, pos + 1)
            self.assertEqual(index.occurrences(seq), expected, seq[:20])

    def table(self, name, key_columns):
        with open(os.path.join(self.db, name) if not os.path.isabs(name) else name) as fh:
            rows = [line.rstrip("\n").split("\t") for line in fh if not line.startswith("#")]
        header = rows[0]
        return [dict(zip(header, r)) for r in rows[1:]] if key_columns else rows

    def representatives(self):
        """{accession: (species taxid, family taxid)} of the representatives."""
        with open(os.path.join(self.db, "internal_taxonomy.dmp")) as fh:
            next(fh)
            nodes = {int(r[0]): r for r in (l.rstrip("\n").split("\t") for l in fh)}
        reps = {}
        for taxid, r in nodes.items():
            if r[4] == "species":
                parent = int(r[1])
                while nodes[parent][4] != "family":
                    parent = int(nodes[parent][1])
                reps[r[6]] = (taxid, parent)
        return reps

    def truth(self):
        """The simulator's marker positions in every genome, of the genes its species has in the database (its
        representative's): {accession: [(gene id, contig, start, end, strand)]} (1-based, inclusive); {accession:
        (species taxid, family taxid, representative accession)}; the contigs' lengths."""
        with open(os.path.join(self.db, "gene2geneid.tsv")) as fh:
            gene_id = dict(line.rstrip("\n").split("\t")[:2] for line in fh if line.strip())
        reps = self.representatives()
        lineage, lengths = {}, {}
        with open(os.path.join(self.gtdb, "simulation", "genomes.tsv")) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                lineage[r["accession"]] = r["gtdb_taxonomy"]
                lengths.update({name: len(seq) for name, seq in read_fasta(r["fasta_path"]).items()})
        rep_of = {lineage[a]: a for a in reps}
        genomes = {a: (*reps[rep_of[l]], rep_of[l]) for a, l in lineage.items()}
        raw = collections.defaultdict(list)
        with open(os.path.join(self.gtdb, "simulation", "marker_positions.tsv")) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                raw[r["accession"]].append((int(gene_id[r["marker"]]), r["contig"], int(r["start"]), int(r["end"]),
                                            r["strand"]))
        positions = {}
        for acc, (_, _, rep) in genomes.items():
            database = {g for g, *_ in raw[rep]}
            positions[acc] = [p for p in raw[acc] if p[0] in database]
        return positions, genomes, lengths

    def ends(self, genes, lengths):
        """One genome's informative gene ends: {(gene, end): (partner, partner end, gap)}, partner 0 for none."""
        by_contig = collections.defaultdict(list)
        for g, c, s, e, st in genes:
            by_contig[c].append((s, e, g, st))
        out = {}
        for contig, row in by_contig.items():
            row.sort()
            for i, (s, e, g, st) in enumerate(row):
                right, left = (3, 5) if st == "+" else (5, 3)
                if i + 1 < len(row) and row[i + 1][0] - 1 - e <= self.MAX_GAP:
                    n = row[i + 1]
                    out[(g, right)] = (n[2], 5 if n[3] == "+" else 3, n[0] - 1 - e)
                elif i + 1 < len(row) or lengths[contig] - e >= self.MAX_GAP:
                    out[(g, right)] = (0, 0, 0)
                if i > 0 and s - 1 - row[i - 1][1] <= self.MAX_GAP:
                    p = row[i - 1]
                    out[(g, left)] = (p[2], 3 if p[3] == "+" else 5, s - 1 - p[1])
                elif i > 0 or s - 1 >= self.MAX_GAP:
                    out[(g, left)] = (0, 0, 0)
        return out

    def expected(self, exclude=()):
        """The family lines the marker positions give, without the species of `exclude`: {(family, gene, end,
        partner, partner end): [the gap of each species]}, {(family, gene, end): informative species}. A species'
        partner at an end is the one most of its genomes informative there show, its representative's on a tie."""
        positions, genomes, lengths = self.truth()
        species = collections.defaultdict(dict)
        for acc, (taxid, family, rep) in genomes.items():
            if taxid not in exclude:
                species[(family, rep)][acc] = self.ends(positions[acc], lengths)
        lines, informative = collections.defaultdict(list), collections.Counter()
        for (family, rep), ends in species.items():
            for key in set().union(*ends.values()):
                seen = {acc: e[key] for acc, e in ends.items() if key in e}
                votes = collections.Counter(v[:2] for v in seen.values())
                tied = sorted(p for p, n in votes.items() if n == max(votes.values()))
                partner = seen[rep][:2] if rep in seen and seen[rep][:2] in tied else tied[0]
                gaps = sorted(v[2] for v in seen.values() if v[:2] == partner)
                lines[(family, *key, *partner)].append(gaps[len(gaps) // 2] if len(gaps) % 2 else
                                                       (gaps[len(gaps) // 2 - 1] + gaps[len(gaps) // 2]) / 2)
                informative[(family, *key)] += 1
        return lines, informative

    def check_family_lines(self, table, exclude=()):
        lines, informative = self.expected(exclude)
        families = {f for _, f in self.representatives().values()}
        rows = {}
        for r in self.table(table, True):
            key = tuple(int(r[c]) for c in ("clade", "gene", "end", "partner", "partner_end"))
            if key[0] in families:
                rows[key] = r
        self.assertEqual(set(rows), set(lines))
        for key, gaps in lines.items():
            r = rows[key]
            self.assertEqual(int(r["species"]), len(gaps), key)
            self.assertEqual(int(r["informative"]), informative[key[:3]], key)
            if key[3]:
                self.assertLessEqual(abs(int(r["gap_min"]) - min(gaps)), self.SLACK, key)
                self.assertLessEqual(abs(int(r["gap_max"]) - max(gaps)), self.SLACK, key)
        # Within the simulator's clusters, but where a genome lost the gene between two.
        gaps = sorted(int(r["gap_median"]) for r in rows.values() if r["partner"] != "0")
        self.assertLessEqual(gaps[len(gaps) // 2], 150)
        return rows

    def test_genes_are_found_where_the_simulator_put_them(self):
        positions, genomes, _ = self.truth()
        found = collections.defaultdict(dict)
        for r in self.table(self.positions, True):
            found[r["accession"]][int(r["gene"])] = r
        self.assertEqual(set(found), set(genomes))  # every genome used
        traced = 0
        for acc, genes in positions.items():
            self.assertEqual(set(found[acc]), {g for g, *_ in genes}, acc)  # each of its database genes, no other
            for g, contig, start, end, strand in genes:
                r = found[acc][g]
                self.assertEqual((r["contig"], r["strand"], r["circular"]), (contig, strand, "0"), (acc, g))
                if genomes[acc][2] == acc:
                    self.assertEqual((int(r["start"]), int(r["end"]), r["placed"]), (start, end, "exact"), (acc, g))
                else:
                    self.assertLessEqual(abs(int(r["start"]) - start), self.SLACK, (acc, g))
                    self.assertLessEqual(abs(int(r["end"]) - end), self.SLACK, (acc, g))
                    traced += r["placed"] == "trace"
        self.assertGreater(traced, 0.9 * sum(len(g) for a, g in positions.items() if genomes[a][2] != a))
        self.assertIn("16 genomes of 8 species", self.output)
        self.assertIn("(8 with their representative genome)", self.output)
        self.assertIn("circular: 0 with a whole replicon by its header, 0 representatives in one sequence", self.output)

    def test_family_lines_count_the_neighbours_of_its_species(self):
        rows = self.check_family_lines("gene_neighbours.tsv")
        # The families differ: a cluster that one of them broke up.
        by_family = collections.defaultdict(set)
        for family, gene, end, partner, partner_end in rows:
            if partner:
                by_family[family].add((gene, end, partner, partner_end))
        self.assertEqual(len(by_family), 3)
        self.assertGreater(len({frozenset(v) for v in by_family.values()}), 1)
        # The strains' contigs end elsewhere than their representatives': ends informative in the species that are
        # not in its representative.
        positions, genomes, lengths = self.truth()
        filled = 0
        for acc, (_, _, rep) in genomes.items():
            if acc != rep:
                filled += len(set(self.ends(positions[acc], lengths)) - set(self.ends(positions[rep], lengths)))
        self.assertGreater(filled, 0)

    def test_the_order_holds_both_families(self):
        rows = self.table("gene_neighbours.tsv", True)
        clades = collections.Counter(int(r["clade"]) for r in rows)
        families = {f for _, f in self.representatives().values()}
        self.assertTrue(families < set(clades))  # family, order, class, phylum, domain lines
        informative = max(int(r["informative"]) for r in rows)
        self.assertEqual(informative, 6)  # an end informative in all six species of the bacterial order
        self.assertIn("family: 3 clades", self.output)
        self.assertRegex(self.output, r"family: .* gene ends informative in fewer than 3 species \(they lean mostly on "
                                      r"the rank above\)")

    def test_the_positions_alone_give_the_same_table(self):
        derived = os.path.join(self.tmp.name, "derived.tsv")
        run(GENE_NEIGHBOURS, "--db", self.db, "--from_positions", self.positions, "--output", derived)
        with open(derived) as a, open(os.path.join(self.db, "gene_neighbours.tsv")) as b:
            self.assertEqual(a.read(), b.read())

    def test_a_training_copy_derives_the_table_without_its_species(self):
        src, dst = os.path.join(self.tmp.name, "db_positions"), os.path.join(self.tmp.name, "db_training")
        shutil.copytree(self.db, src)
        excluded = os.path.join(self.tmp.name, "excluded.txt")
        with open(excluded, "w") as fh:
            fh.write("Mockella s0\n")
        run(CONVERT, "--from_db", src, "--exclude_species", excluded, "--outdir", dst)
        taxid = next(t for t, _ in self.representatives().values() if self.species_name(t) == "s__Mockella s0")
        self.check_family_lines(os.path.join(dst, "gene_neighbours.tsv"), exclude={taxid})
        with open(self.positions) as fh:
            lines = fh.readlines()
        kept = [line for line in lines if line.startswith("#") or line.split("\t")[1] != str(taxid)]
        with open(os.path.join(dst, "gene_positions.tsv")) as fh:
            self.assertEqual(fh.readlines(), kept)
        self.assertLess(len(kept), len(lines))

    def species_name(self, taxid):
        with open(os.path.join(self.db, "internal_taxonomy.dmp")) as fh:
            next(fh)
            return next(r[3] for r in (l.rstrip("\n").split("\t") for l in fh) if int(r[0]) == taxid)


class SpeciesLinesTest(unittest.TestCase):
    """count_clades' lines of a species' own: six species of one family, five with genes 1, 2 and 3 in the order 1 2
    ... 3, one with 1 3 ... 2. Only the odd one gets lines of its own, at the three ends where its partner is rare in
    the family (1/6, below protal's expected share of 0.2); none for a family too small to judge, or without
    species lines."""

    @staticmethod
    def world(odd_order, species=6):
        nodes = {1: (1, "no rank", "root"), 2: (1, "domain", "d"), 10: (2, "phylum", "p"), 20: (10, "class", "c"),
                 30: (20, "order", "o"), 40: (30, "family", "f"), 50: (40, "genus", "g")}
        genomes = {}
        for i in range(species):
            taxid = 101 + i
            nodes[taxid] = (50, "species", f"s{i}")
            order = odd_order if i == species - 1 else (1, 2, 3)
            starts = (10000, 11100, 50000)
            placements = [(gene, "c1", 100000, False, start, start + 900, "+") for gene, start in zip(order, starts)]
            genomes[taxid] = {f"G{i}": placements}
        return genomes, nodes

    def lines(self, genomes, nodes, **kwargs):
        ranks = ["family", "order", "class", "phylum", "domain"]
        counts, informative, _ = gn.count_clades(genomes, nodes, ranks, {t: next(iter(g)) for t, g in genomes.items()},
                                                 3000, **kwargs)
        return {(clade, gene, end, *partner): (len(gaps), informative[(clade, gene, end)])
                for (clade, gene, end), partners in counts.items() for partner, gaps in partners.items()
                if nodes[clade][1] == "species"}

    def test_only_the_odd_species_gets_lines_where_its_partner_is_rare(self):
        genomes, nodes = self.world((1, 3, 2))
        self.assertEqual(self.lines(genomes, nodes), {(106, 1, 3, 3, 5): (1, 1), (106, 3, 5, 1, 3): (1, 1),
                                                      (106, 2, 5, 0, 0): (1, 1)})

    def test_none_without_them_or_in_a_family_too_small_to_judge(self):
        genomes, nodes = self.world((1, 3, 2))
        self.assertEqual(self.lines(genomes, nodes, species_lines=False), {})
        genomes, nodes = self.world((1, 3, 2), species=4)
        self.assertEqual(self.lines(genomes, nodes), {})

    def test_a_strain_s_own_partners_get_lines_too(self):
        # The odd species' representative has 1 3 ... 2, two other strains of it the family's 1 2 ... 3: the species
        # counts as the family does (its genomes' most common partner) in the clades, but the representative's
        # partners, which the reads of that genome show, get lines of the species', as do its 2' 5' end's
        # partners (gene 1 in two genomes, none in one).
        genomes, nodes = self.world((1, 3, 2))
        strain = [(gene, "c1", 100000, False, start, start + 900, "+") for gene, start in zip((1, 2, 3), (10000, 11100, 50000))]
        genomes[106]["G5b"] = strain
        genomes[106]["G5c"] = strain
        self.assertEqual(self.lines(genomes, nodes), {(106, 1, 3, 3, 5): (1, 1), (106, 3, 5, 1, 3): (1, 1),
                                                      (106, 2, 5, 0, 0): (1, 1)})


class CircularGeneNeighboursTest(unittest.TestCase):
    """gene_neighbours.py on a release whose genomes are each one sequence, short enough that a species' last
    cluster lies within --max_gap of its first across the origin: a representative in one sequence is circular,
    its last gene facing its first and every end informative; another strain's one sequence is not."""

    MAX_GAP = 3000

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        lineages = os.path.join(cls.tmp.name, "lineages.txt")
        with open(lineages, "w") as fh:
            fh.writelines(f"d__Bacteria;p__Simulatota;c__Simulatia;o__Simulales;f__Simulaceae;g__Mockella;"
                          f"s__Mockella s{i}\n" for i in range(3))
        cls.gtdb, cls.db = os.path.join(cls.tmp.name, "gtdb"), os.path.join(cls.tmp.name, "db")
        run(SIMULATE, "--outdir", cls.gtdb, "--lineages", lineages, "--operons", "--genome_length", "40000",
            "--genomes_per_species", "2", "--contigs", "1")
        run(CONVERT, "--gtdb", cls.gtdb, "--outdir", cls.db)
        cls.output = subprocess.run([sys.executable, GENE_NEIGHBOURS, "--db", cls.db, "--genome_table",
                                     os.path.join(cls.gtdb, "simulation", "genomes.tsv")],
                                    check=True, capture_output=True, text=True).stdout
        with open(os.path.join(cls.db, "gene_positions.tsv")) as fh:
            cls.positions = [r for r in csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t")]
        with open(os.path.join(cls.db, "gene_neighbours.tsv")) as fh:
            cls.lines = [r for r in csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t")]

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_a_representative_in_one_sequence_is_circular(self):
        self.assertIn("circular: 0 with a whole replicon by its header, 3 representatives in one sequence", self.output)
        circular = {(r["accession"], r["circular"]) for r in self.positions}
        self.assertEqual({c for a, c in circular if a.startswith("GCF_")}, {"1"})  # the representatives
        self.assertEqual({c for a, c in circular if a.startswith("GCA_")}, {"0"})  # the other strains

    def test_the_last_gene_faces_the_first_across_the_origin(self):
        by_genome = collections.defaultdict(list)
        for r in self.positions:
            if r["circular"] == "1":
                by_genome[r["accession"]].append(r)
        self.assertEqual(len(by_genome), 3)
        lines = {(r["gene"], r["end"], r["partner"], r["partner_end"]): r for r in self.lines}
        for genes in by_genome.values():
            genes.sort(key=lambda r: int(r["start"]))
            first, last = genes[0], genes[-1]
            wrap = int(last["contig_length"]) - int(last["end"]) + int(first["start"]) - 1
            self.assertLess(wrap, self.MAX_GAP)  # the release is short enough to test this
            last_end = "3" if last["strand"] == "+" else "5"
            first_end = "5" if first["strand"] == "+" else "3"
            for key in ((last["gene"], last_end, first["gene"], first_end), (first["gene"], first_end, last["gene"], last_end)):
                self.assertIn(key, lines)
                self.assertLessEqual(int(lines[key]["gap_min"]), wrap)
                self.assertGreaterEqual(int(lines[key]["gap_max"]), wrap)

    def test_neighbour_ends_of_a_circular_and_a_linear_contig(self):
        genes = [(1, 100, 1000, "+"), (2, 5000, 6000, "-"), (3, 9500, 9900, "+")]
        circular = gn.neighbour_ends([(g, "c", 10000, True, s, e, st) for g, s, e, st in genes], self.MAX_GAP)
        linear = gn.neighbour_ends([(g, "c", 10000, False, s, e, st) for g, s, e, st in genes], self.MAX_GAP)
        self.assertIn((1, 5, 3, 3, 200), circular)  # gene 1's 5' end faces gene 3's 3' end across the origin
        self.assertIn((3, 3, 1, 5, 200), circular)
        self.assertIn((1, 3, 0, 0, 0), circular)    # 4000 bases to gene 2: no neighbour
        self.assertEqual(len(circular), 6)          # every end informative
        self.assertEqual({(g, e) for g, e, *_ in linear}, {(1, 3), (2, 5), (2, 3), (3, 5)})  # not those near the ends
        with tempfile.NamedTemporaryFile("wb", suffix=".fna", delete=False) as fh:
            fh.write(b">c1 Mockella s0 chromosome, complete genome\nACGT\n>c2 a contig\nACGT\n>c3 [topology=circular]\nA\n")
        try:
            self.assertEqual([c for _, _, c in gn.read_contigs(fh.name)], [True, False, True])
        finally:
            os.remove(fh.name)


if __name__ == "__main__":
    unittest.main()
