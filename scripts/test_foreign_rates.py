#!/usr/bin/env python3
"""Tests of foreign_rates.py: the database's taxonomy, the tiler's sources file, the counting of a scan's SAM (each
read's best record, MAPQ 4 or more, as a read of its copy's own species, another species of the genus or another genus;
the read's source is the taxid its name starts with) and the table it writes, every copy of the sources listed, which
src/SequenceUtils/ForeignRatesTable.h reads (tests/test_CongenerGaps.cpp reads the same format). Standard library only.

Run: python3 -m unittest scripts/test_foreign_rates.py
"""

import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import foreign_rates  # noqa: E402

# id, parent, external id, name, rank, level, representative
TAXONOMY = ("1\t1\t-\troot\tno rank\t0\t\n"
            "10\t1\t-\tg__G\tgenus\t6\t\n"
            "11\t1\t-\tg__H\tgenus\t6\t\n"
            "100\t10\t-\ts__G alpha\tspecies\t7\tGCF_1.1\n"
            "101\t10\t-\ts__G beta\tspecies\t7\tGCA_2.1\n"
            "102\t11\t-\ts__H gamma\tspecies\t7\tGCA_3.1\n")
# simulate_metagenomes --tile_fasta's <tile_out>.sources.tsv: a header again when its gene's records are not contiguous
SOURCES = ("header\trecords\tkept\ttiles\n"
           "100_5\t30\t10\t40\n"
           "101_5\t1\t1\t4\n"
           "102_5\t2\t2\t8\n"
           "100_2\t1\t1\t3\n"
           "100_5\t1\t1\t4\n"
           "contig_x\t1\t1\t2\n")


def record(qname, flag, rname, mapq):
    return f"{qname}\t{flag}\t{rname}\t1\t{mapq}\t150=\t*\t0\t0\t*\t*\tZU:i:0\tZT:i:0\n"


class ForeignRatesTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.taxonomy = os.path.join(self.tmp.name, "internal_taxonomy.dmp")
        self.sources = os.path.join(self.tmp.name, "tiles.fq.zst.sources.tsv")
        with open(self.taxonomy, "w") as fh:
            fh.write(TAXONOMY)
        with open(self.sources, "w") as fh:
            fh.write(SOURCES)

    def tearDown(self):
        self.tmp.cleanup()

    def test_taxonomy_and_sources(self):
        taxonomy = foreign_rates.read_taxonomy(self.taxonomy)
        self.assertEqual(taxonomy, {100: ("s__G alpha", "g__G"), 101: ("s__G beta", "g__G"), 102: ("s__H gamma", "g__H")})
        sources, skipped = foreign_rates.read_sources(self.sources)
        self.assertEqual(dict(sources), {(100, 5): [31, 11, 44], (101, 5): [1, 1, 4], (102, 5): [2, 2, 8], (100, 2): [1, 1, 3]})
        self.assertEqual(skipped, 1)
        self.assertEqual(foreign_rates.copy_of("100_5"), (100, 5))
        self.assertIsNone(foreign_rates.copy_of("contig_x"))
        self.assertIsNone(foreign_rates.copy_of("100"))

    def test_each_reads_best_record_counts_for_its_copy(self):
        sam = os.path.join(self.tmp.name, "tiles.sam")
        with open(sam, "w") as fh:
            fh.write("@HD\tVN:1.6\n@SQ\tSN:100_5\tLN:900\n")
            fh.write(record("100_5:0:0", 0, "100_5", 60))      # its own species (another genome's copy, record 0)
            fh.write(record("100_5:7:250", 0, "100_5", 60))    # its own species, record 7
            fh.write(record("101_5:3:0", 0, "100_5", 30))      # another species of the genus
            fh.write(record("101_5:3:0", 256, "101_5", 0))     # its secondary record: not counted again
            fh.write(record("102_5:4:0", 0, "100_5", 12))      # another genus
            fh.write(record("102_5:4:250", 0, "100_5", 2))     # MAPQ below 4: not counted
            fh.write(record("102_5:4:500", 4, "*", 0))         # unmapped
            fh.write(record("999_5:9:0", 0, "100_5", 60))      # a species the taxonomy lacks
            fh.write(record("101_7:5:0", 2048, "100_5", 60))   # a supplementary record first: the next one is the best
            fh.write(record("101_7:5:0", 0, "102_7", 60))
        counts, seen, counted, unknown = foreign_rates.count(sam, foreign_rates.read_taxonomy(self.taxonomy))
        self.assertEqual(dict(counts), {(100, 5): [2, 2, 1], (102, 7): [0, 1, 1]})  # own, foreign, of them other genera
        self.assertEqual((seen, counted, unknown), (8, 5, 1))
        # The table's reads: the foreign reads plus one genome's worth of the copy's own (own over the header's records
        # tiled: 100_5 has 11, so its 2 own reads are 0 per genome; 100_2's 70000 over 1 record saturate); every copy of
        # the sources listed, with no reads where none landed.
        self.assertEqual(foreign_rates.table_reads(40, 8, 10), 12)
        self.assertEqual(foreign_rates.table_reads(0, 3, 0), 3)
        out = os.path.join(self.tmp.name, "foreign_rates.tsv")
        sources, _ = foreign_rates.read_sources(self.sources)
        listed = foreign_rates.write_table(out, {**counts, (100, 2): [70000, 0, 0], (101, 5): [40, 8, 2]}, sources)
        self.assertEqual(listed, 5)
        with open(out) as fh:
            lines = fh.read().splitlines()
        self.assertTrue(lines[0].startswith("#"))
        self.assertEqual(lines[1:], ["taxid\tgene:reads:foreign:foreign_genus", "100\t2:65535:0:0,5:2:2:1", "101\t5:48:8:2",
                                     "102\t5:0:0:0,7:1:1:1"])


if __name__ == "__main__":
    unittest.main()
