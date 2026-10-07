#!/usr/bin/env python3
"""Tests of error_read_features.py on the synthetic world of scripts/test_error_reads.py: a paired-end sample of
species A (true positive), B (present, not called: FN), C (present, no row: unseen) and E (held out of the training
database), and the absent taxon D, called on E's reads (FP).

Run: python3 docs/claude/2026-10-07-error-read-signatures/test_error_read_features.py
"""
import csv
import gzip
import os
import sys
import tempfile
import unittest
import zlib

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..", "..", "scripts")))
import compressed  # noqa: E402
import error_read_features as erf  # noqa: E402
import test_error_reads as world  # noqa: E402

TAXONOMY = ("id\tparent_id\texternal_id\tname\trank\tlevel\trep_genome\n"
            "11\t11\t0\tf__F\tfamily\t5\t\n10\t11\t0\tg__G\tgenus\t6\t\n"
            + "".join(f"{i}\t10\t0\ts__G {sp}\tspecies\t7\t{genome or ''}\n"
                      for i, (sp, (_, genome, _)) in enumerate(world.SPECIES.items(), 1)))


@unittest.skipUnless(world.HAVE_ZSTD, "needs the zstd command (or Python 3.14)")
class ErrorReadFeatures(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        maker = world.ErrorReads()
        maker.tmp = self.tmp
        self.args = maker.world()
        d = self.tmp.name
        with open(os.path.join(d, "internal_taxonomy.dmp"), "w") as fh:
            fh.write(TAXONOMY)
        with open(os.path.join(d, "heldout_species.txt"), "w") as fh:
            fh.write("s__G E\tspecies\ts__G E\n")
        for registry in (erf.TAXONOMY, erf.NAME_ID, erf.DB_SPECIES, erf.HELD_OUT):
            registry.clear()

    def tearDown(self):
        self.tmp.cleanup()

    def run_it(self, *extra):
        d = self.tmp.name
        out = os.path.join(d, "features")
        rc = erf.main(["--calls", os.path.join(d, "trained_model.calls.tsv.gz"), "--training",
                       os.path.join(d, "training"), "--db", os.path.join(d, "training_db"), "--taxonomy",
                       os.path.join(d, "internal_taxonomy.dmp"), "--heldout", os.path.join(d, "heldout_species.txt"),
                       "--samples", "design", "--threads", "1", "--out", out, *extra])
        self.assertEqual(rc, 0)
        base = os.path.join(out, "training", world.POINT, world.SAMPLE)
        with gzip.open(base + ".taxa.tsv.gz", "rt") as fh:
            taxa = {r["taxon"]: r for r in csv.DictReader(fh, delimiter="\t")}
        with gzip.open(base + ".records.tsv.gz", "rt") as fh:
            records = list(csv.DictReader(fh, delimiter="\t"))
        with compressed.open_read(base + ".sam.zst") as fh:
            sam = fh.read().decode()
        return out, taxa, records, sam

    def test_fp_and_fn_reads_and_no_unseen(self):
        out, taxa, records, sam = self.run_it()
        contig = {sp: c for sp, (_, _, c) in world.SPECIES.items() if c}
        frags = {r["qname"]: r["reasons"] for r in records}
        self.assertEqual(frags, {f"{contig['E']}-2": "FP:4", f"{contig['B']}-3": "FN:2,source:2",
                                 f"{contig['B']}-4": "source:2", f"{contig['B']}-5": "seeded:2,source:2"})
        self.assertEqual(len(records), 8)
        self.assertNotIn(contig["C"], sam)  # the unseen species' read is not taken
        self.assertNotIn(f"{contig['A']}-1\t", sam)
        fp = next(r for r in records if r["reasons"] == "FP:4" and r["taxon"] == "4" and r["mate"] == "1")
        self.assertEqual((fp["relation"], fp["source_status"], fp["source_taxid"], fp["taxon_call"]),
                         ("genus", "held_out", "5", "FP"))
        self.assertEqual((fp["identity"], fp["gene_len"], fp["zu"], fp["read_len"]), ("1.0", "900", "3", "150"))
        self.assertIn("xe:Z:FP:4", sam)
        self.assertIn("@SQ\tSN:4_1\tLN:900", sam)
        self.assertNotIn("protal failed candidates", sam)
        a, b, d = taxa["1"], taxa["2"], taxa["4"]
        self.assertEqual((a["class"], a["frags_any"], a["frags_best"], a["frags_counted"], a["frags_lowmapq"]),
                         ("TP", "4", "3", "2", "1"))
        self.assertEqual((a["src_own"], a["own_frags"], a["own_counted_here"], a["own_unaligned"], a["seeded_failed"]),
                         ("2", "3", "2", "1", "2"))
        self.assertEqual(a["zf_share"], "0.5")  # A-8 seeded on C too
        self.assertEqual((b["class"], b["frags_counted"], b["own_frags"], b["own_counted_here"],
                          b["own_lowmapq_elsewhere"], b["own_unaligned"], b["seeded_failed"]),
                         ("FN", "1", "3", "1", "1", "1", "1"))
        self.assertEqual((d["class"], d["frags_counted"], d["src_genus"], d["src_held_out"], d["other_taxon_share"],
                          d["both_mates_share"]), ("FP", "1", "1", "1", "1.0", "1.0"))
        with open(os.path.join(out, "summary.tsv")) as fh:
            summary = list(csv.DictReader(fh, delimiter="\t"))
        self.assertEqual((summary[0]["FP"], summary[0]["FN"], summary[0]["kept_fragments"],
                          summary[0]["not_consecutive"]), ("1", "1", "4", "0"))

    def test_cap_keeps_lowest_crc(self):
        contig = {sp: c for sp, (_, _, c) in world.SPECIES.items() if c}
        names = [f"{contig['B']}-{n}" for n in (3, 4, 5)]
        lowest = min(names, key=lambda q: (zlib.crc32(q.encode()), q))
        _, _, records, _ = self.run_it("--max-fragments", "1")
        kept = {r["qname"] for r in records}
        self.assertEqual(kept, {f"{contig['E']}-2", f"{contig['B']}-3", f"{contig['B']}-5", lowest})

    def test_drawn_long_reads(self):
        d = self.tmp.name
        out = os.path.join(d, "features_pb")
        self.assertEqual(erf.main(["--calls", os.path.join(d, "trained_model_pb.calls.tsv.gz"), "--training",
                                   os.path.join(d, "training"), "--db", os.path.join(d, "training_db"), "--taxonomy",
                                   os.path.join(d, "internal_taxonomy.dmp"), "--read-type", "pb", "--threads", "1",
                                   "--out", out]), 0)
        with gzip.open(os.path.join(out, "training", world.LONG_POINT, world.LONG_SAMPLE + ".records.tsv.gz"),
                       "rt") as fh:
            records = {r["qname"]: r for r in csv.DictReader(fh, delimiter="\t")}
        self.assertEqual({q: r["reasons"] for q, r in records.items()},
                         {"g1x_1": "FN:2,source:2", "g3x_2": "FP:4"})
        self.assertEqual((records["g3x_2"]["source_species"], records["g3x_2"]["relation"],
                          records["g1x_1"]["relation"], records["g1x_1"]["source_genome_kind"]),
                         ("s__G E", "genus", "own", "rep"))

    def test_cigar_counts(self):
        self.assertEqual(erf.cigar_counts(b"5S10M1X2I3M1D4M2H"), (17, 1, 2, 1, 2, 5, 2, 27))


if __name__ == "__main__":
    unittest.main()
