#!/usr/bin/env python3
"""Tests of foreign_rates.py: the genome table's species and genera, the database's taxonomy, the counting of a scan's
SAM (each read's best record, MAPQ 4 or more, as a read of its copy's own species, another species of the genus or
another genus) and the table it writes, which src/SequenceUtils/ForeignRatesTable.h reads (tests/test_CongenerGaps.cpp
reads the same format). Standard library only.

Run: python3 -m unittest scripts/test_foreign_rates.py
"""

import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import foreign_rates  # noqa: E402

GENOMES = ("GCF_1.1\td__Bacteria;p__P;c__C;o__O;f__F;g__G;s__G alpha\t/x/a.fna\n"
           "GCA_2.1\td__Bacteria;p__P;c__C;o__O;f__F;g__G;s__G beta\t/x/b.fna\n"
           "GCA_3.1\td__Bacteria;p__P;c__C;o__O;f__F;g__H;s__H gamma\t/x/c.fna\n")
# id, parent, external id, name, rank, level, representative
TAXONOMY = ("1\t1\t-\troot\tno rank\t0\t\n"
            "10\t1\t-\tg__G\tgenus\t6\t\n"
            "11\t1\t-\tg__H\tgenus\t6\t\n"
            "100\t10\t-\ts__G alpha\tspecies\t7\tGCF_1.1\n"
            "101\t10\t-\ts__G beta\tspecies\t7\tGCA_2.1\n"
            "102\t11\t-\ts__H gamma\tspecies\t7\tGCA_3.1\n")


def record(qname, flag, rname, mapq):
    return f"{qname}\t{flag}\t{rname}\t1\t{mapq}\t150=\t*\t0\t0\t*\t*\tZU:i:0\tZT:i:0\n"


class ForeignRatesTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.genomes = os.path.join(self.tmp.name, "genomes.tsv")
        self.taxonomy = os.path.join(self.tmp.name, "internal_taxonomy.dmp")
        with open(self.genomes, "w") as fh:
            fh.write(GENOMES)
        with open(self.taxonomy, "w") as fh:
            fh.write(TAXONOMY)

    def tearDown(self):
        self.tmp.cleanup()

    def test_genomes_and_taxonomy(self):
        self.assertEqual(foreign_rates.read_genomes(self.genomes)["GCA_3.1"], ("s__H gamma", "g__H"))
        taxonomy = foreign_rates.read_taxonomy(self.taxonomy)
        self.assertEqual(taxonomy, {100: ("s__G alpha", "g__G"), 101: ("s__G beta", "g__G"), 102: ("s__H gamma", "g__H")})
        self.assertEqual(foreign_rates.last_rank("d__B;g__X;s__X y", "g__"), "g__X")
        self.assertEqual(foreign_rates.last_rank("d__B", "s__"), "")

    def test_each_reads_best_record_counts_for_its_copy(self):
        sam = os.path.join(self.tmp.name, "tiles.sam")
        with open(sam, "w") as fh:
            fh.write("@HD\tVN:1.6\n@SQ\tSN:100_5\tLN:900\n")
            fh.write(record("GCF_1.1:0:0", 0, "100_5", 60))     # its own species
            fh.write(record("GCF_1.1:0:500", 0, "100_5", 60))
            fh.write(record("GCA_2.1:3:1000", 0, "100_5", 30))  # another species of the genus
            fh.write(record("GCA_2.1:3:1000", 256, "101_5", 0))  # its secondary record: not counted again
            fh.write(record("GCA_3.1:0:0", 0, "100_5", 12))     # another genus
            fh.write(record("GCA_3.1:0:500", 0, "100_5", 2))    # MAPQ below 4: not counted
            fh.write(record("GCA_3.1:0:1000", 4, "*", 0))       # unmapped
            fh.write(record("GCX_9.9:0:0", 0, "100_5", 60))     # a genome the table lacks
            fh.write(record("GCA_2.1:0:0", 2048, "100_5", 60))  # a supplementary record first: the next one is the best
            fh.write(record("GCA_2.1:0:0", 0, "102_7", 60))
        counts, seen, counted, unknown = foreign_rates.count(sam, foreign_rates.read_genomes(self.genomes),
                                                             foreign_rates.read_taxonomy(self.taxonomy))
        self.assertEqual(dict(counts), {(100, 5): [4, 2, 1], (102, 7): [1, 1, 1]})
        self.assertEqual((seen, counted, unknown), (8, 5, 1))
        out = os.path.join(self.tmp.name, "foreign_rates.tsv")
        foreign_rates.write_table(out, {(102, 7): [1, 1, 1], (100, 5): [4, 2, 1], (100, 2): [70000, 0, 0]})
        with open(out) as fh:
            lines = fh.read().splitlines()
        self.assertTrue(lines[0].startswith("#"))
        self.assertEqual(lines[1:], ["taxid\tgene:reads:foreign:foreign_genus", "100\t2:65535:0:0,5:4:2:1", "102\t7:1:1:1"])


if __name__ == "__main__":
    unittest.main()
