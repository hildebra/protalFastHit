#!/usr/bin/env python3
"""Tests of error_reads.py: the records of the reads behind a model's errors, taken from a sample's SAM.

A synthetic paired-end sample of species A (called: true positive), B (present, not called: FN), C (present, in the
training database, no row: unseen) and E (one the database lacks), and the absent taxon D, called on E's reads (FP);
and a PacBio sample of the same community, its reads named after their genome's place in the manifest. Standard library
only (and the zstd command, as compressed.py needs without Python 3.14).

Run: python3 -m unittest scripts/test_error_reads.py
"""

import csv
import gzip
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import compressed  # noqa: E402
import error_reads  # noqa: E402

LINEAGE = "d__Bacteria;p__P;c__C;o__O;f__F;g__G;s__G {}"
POINT, SAMPLE = "rl150_p1000", "rl150_p1000_s_1"
LONG_POINT, LONG_SAMPLE = "pb_b300000", "pb_b300000_s_1"
# species: (taxid in the training database or None, genome, contig); the manifest lists the genomes in this order
SPECIES = {"A": ("1", "GCF_000000001.1", "NZ_CP000001.1"), "B": ("2", "GCA_000000002.1", "JAAAAA010000001.1"),
           "C": ("3", "GCA_000000003.1", "NZ_CP000003.1"), "D": ("4", None, None),
           "E": (None, "GCA_000000005.1", "JAEEEE010000001.1")}
SAM_TEXT = """@HD\tVN:1.6
@SQ\tSN:1_1\tLN:900
@SQ\tSN:2_1\tLN:900
@SQ\tSN:4_1\tLN:900
@SQ\tSN:9_9\tLN:900
@CO\tprotal read type: pe
{A}-1\t65\t1_1\t10\t60\t150M\t=\t200\t0\tA\tI\tZU:i:3\tZT:i:0
{A}-1\t129\t1_1\t200\t60\t150M\t=\t10\t0\tA\tI\tZU:i:3\tZT:i:0
{E}-2\t65\t4_1\t10\t60\t150M\t=\t200\t0\tA\tI\tZU:i:3\tZT:i:0
{E}-2\t129\t4_1\t200\t60\t150M\t=\t10\t0\tA\tI\tZU:i:3\tZT:i:0
{E}-2\t321\t1_1\t10\t0\t150M\t9_9\t200\t0\tA\tI\tZU:i:0\tZT:i:0
{B}-3\t65\t2_1\t10\t60\t150M\t=\t200\t0\tA\tI\tZU:i:3\tZT:i:0
{B}-3\t129\t2_1\t200\t60\t150M\t=\t10\t0\tA\tI\tZU:i:3\tZT:i:0
{B}-4\t65\t1_1\t10\t2\t150M\t=\t200\t0\tA\tI\tZU:i:0\tZT:i:0
{B}-4\t129\t1_1\t200\t2\t150M\t=\t10\t0\tA\tI\tZU:i:0\tZT:i:0
{B}-5\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\tZU:i:0\tZT:i:0\tZF:Z:2,1
{C}-6\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\tZU:i:0\tZT:i:0\tZF:Z:3
{A}-7\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\tZU:i:0\tZT:i:0\tZF:Z:1
{A}-8\t65\t1_1\t10\t60\t150M\t=\t200\t0\tA\tI\tZU:i:3\tZT:i:0\tZF:Z:3
{A}-8\t129\t1_1\t200\t60\t150M\t=\t10\t0\tA\tI\tZU:i:3\tZT:i:0
"""
# The PacBio sample: g<i>x_<n>, i the genome's place in the manifest (A 0, B 1, C 2, E 3).
LONG_SAM = """@HD\tVN:1.6
@SQ\tSN:2_1\tLN:900
@CO\tprotal read type: pb
g1x_1\t0\t2_1\t10\t60\t900M\t*\t0\t0\tA\tI\tZU:i:3\tZT:i:0
g3x_2\t0\t4_1\t10\t60\t900M\t*\t0\t0\tA\tI\tZU:i:3\tZT:i:0
g0x_3\t0\t1_1\t10\t60\t900M\t*\t0\t0\tA\tI\tZU:i:3\tZT:i:0
g2x_4\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\tZU:i:0\tZT:i:0\tZF:Z:3
"""
COLUMNS = "meta_design\tmeta_sample\tmeta_read_type\tmeta_scenario\ttaxon\ttaxon_name\tdomain\ttruth\tset\tp\tknob\tcall\n"


class ErrorReads(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        error_reads.CONTIGS.clear()
        error_reads.DB_TAXID.clear()
        error_reads.DB_NAME.clear()

    def tearDown(self):
        self.tmp.cleanup()

    def world(self):
        """A training collection with the paired-end sample and a PacBio sample of its community, the calls of both and
        the training database. -> error_reads.py's arguments but --read-type and --calls."""
        d = self.tmp.name
        point = os.path.join(d, "training", "points", POINT)
        sim, protal = os.path.join(point, "sim"), os.path.join(point, "protal")
        os.makedirs(sim)
        os.makedirs(os.path.join(protal, "alignments"))
        with open(os.path.join(sim, "protal.meta"), "w") as fh:
            fh.write(f"#OUTPUT_DIR\t{protal}\n#INPUT_DIR\t{sim}\n#SAMPLEID\tFIRST\tSECOND\tSAM\tPREFIX\tPROFILE\tPROFILE_TRUTH\n"
                     f"{SAMPLE}\tr1.fq.zst\tr2.fq.zst\t{SAMPLE}.sam.zst\t{SAMPLE}\t{SAMPLE}.profile\t/truth\n")
        columns = ["sample", "genome", "taxonomy", "read_pairs", "fasta_path"]
        with open(os.path.join(sim, "manifest.tsv"), "w") as fh:
            fh.write("\t".join(columns) + "\n")
            for sp, (_, genome, contig) in SPECIES.items():
                if genome is None:
                    continue
                fasta = os.path.join(d, genome + ".fna.gz")
                with gzip.open(fasta, "wt") as out:
                    out.write(f">{contig} {sp} chromosome\nACGT\n")
                fh.write("\t".join([SAMPLE, genome, LINEAGE.format(sp), "100", fasta]) + "\n")
        with compressed.open_write(os.path.join(protal, "alignments", SAMPLE + ".sam.zst")) as fh:
            fh.write(SAM_TEXT.format(**{sp: c for sp, (_, _, c) in SPECIES.items() if c}).encode())
        long_point = os.path.join(d, "training", "points", LONG_POINT)
        os.makedirs(os.path.join(long_point, "sim"))
        os.makedirs(os.path.join(long_point, "protal", "alignments"))
        with open(os.path.join(long_point, "sim", "samples.tsv"), "w") as fh:
            fh.write("sample\treads\ttruth\tcommunity\n" f"{LONG_SAMPLE}\tr.fq.zst\t/truth\t{SAMPLE}\n")
        with gzip.open(os.path.join(long_point, "protal", "alignments", LONG_SAMPLE + ".sam.gz"), "wt") as fh:
            fh.write(LONG_SAM)
        db = os.path.join(d, "training_db")
        os.makedirs(db)
        with open(os.path.join(db, "genome2tiid.tsv"), "w") as fh:
            for sp, (taxid, genome, _) in SPECIES.items():
                if taxid:
                    fh.write(f"{genome or 'GCF_000000004.1'}\t{taxid}\t{genome}\t{LINEAGE.format(sp)}\n")
        with gzip.open(os.path.join(d, "trained_model.calls.tsv.gz"), "wt") as fh:
            fh.write(COLUMNS)
            for sp, truth, p, call in (("A", 1, 0.9, 1), ("B", 1, 0.3, 0), ("D", 0, 0.8, 1)):
                fh.write(f"{POINT}\t{SAMPLE}\tpe\t\t{SPECIES[sp][0]}\ts__G {sp}\tBacteria\t{truth}\ttraining\t{p}\t0.7\t{call}\n")
            # A test sample, whose collection is not given: said, and left out; a scenario's, left out by --samples.
            fh.write(f"{POINT}\t{SAMPLE}\tpe\t\t1\ts__G A\tBacteria\t1\ttest\t0.9\t0.7\t1\n")
            fh.write(f"sc_gut_pe_p1\tsc_gut_pe_p1_s_1\tpe\tgut\t1\ts__G A\tBacteria\t1\ttraining\t0.9\t0.7\t1\n")
        with gzip.open(os.path.join(d, "trained_model_pb.calls.tsv.gz"), "wt") as fh:
            fh.write(COLUMNS)
            for sp, truth, p, call in (("A", 1, 0.9, 1), ("B", 1, 0.3, 0), ("D", 0, 0.8, 1)):
                fh.write(f"{LONG_POINT}\t{LONG_SAMPLE}\tpb\t\t{SPECIES[sp][0]}\ts__G {sp}\tBacteria\t{truth}\ttraining\t"
                         f"{p}\t0.6\t{call}\n")
        return ["--training", os.path.join(d, "training"), "--db", db, "--threads", "1"]

    def run_pe(self, args):
        d = self.tmp.name
        self.assertEqual(error_reads.main(args + ["--calls", os.path.join(d, "trained_model.calls.tsv.gz"),
                                                  "--samples", "design", "--out", os.path.join(d, "out")]), 0)
        return os.path.join(d, "out")

    def test_the_reads_of_the_errors(self):
        out = self.run_pe(self.world())
        with open(os.path.join(out, "summary.tsv")) as fh:
            summary = list(csv.DictReader(fh, delimiter="\t"))
        self.assertEqual(len(summary), 1)
        row = summary[0]
        self.assertEqual((row["set"], row["point"], row["scenario"], row["sample"]), ("training", POINT, "design", SAMPLE))
        self.assertEqual((row["FP"], row["FN"], row["unseen"]), ("1", "1", "1"))
        # E's read on D (FP), B's three (FN: on itself, on A, seeded on B and A but aligned nowhere), C's (unseen,
        # aligned nowhere) and A's read 8, which seeded on C; not A's reads 1 and 7.
        self.assertEqual((row["fragments"], row["records"], row["unknown_source"]), ("6", "11", "0"))
        self.assertEqual(row["sam"], os.path.join("training", POINT, SAMPLE + ".sam.zst"))
        lines = compressed.read_text(os.path.join(out, row["sam"])).splitlines()
        header = [line for line in lines if line.startswith("@")]
        records = [line.split("\t") for line in lines if not line.startswith("@")]
        # The genes the records name, 9_9 too (a mate's RNEXT); none other.
        self.assertEqual([h.split("\t")[1] for h in header if h.startswith("@SQ")], ["SN:1_1", "SN:2_1", "SN:4_1", "SN:9_9"])
        self.assertIn("@CO\tprotal read type: pe", header)
        self.assertTrue(any(h.startswith("@CO\terror_reads.py:") and "1 FP, 1 FN and 1 unseen taxa" in h for h in header))
        self.assertEqual(len(records), 11)
        tags = {}
        for r in records:
            tags.setdefault(r[0].rsplit("-", 1)[1], {t[:2]: t[5:] for t in r[11:] if t[:2] in ("xg", "xs", "xe")})
        self.assertEqual(tags["2"], {"xg": "GCA_000000005.1", "xs": "s__G E (not in the database)", "xe": "FP:4"})
        self.assertEqual(tags["3"], {"xg": "GCA_000000002.1", "xs": "s__G B", "xe": "FN:2,source:2"})
        self.assertEqual(tags["4"]["xe"], "source:2")
        self.assertEqual(tags["5"]["xe"], "seeded:2,source:2")
        self.assertEqual(tags["6"], {"xg": "GCA_000000003.1", "xs": "s__G C", "xe": "seeded:3,source:3"})
        self.assertEqual(tags["8"], {"xg": "GCF_000000001.1", "xs": "s__G A", "xe": "seeded:3"})
        self.assertNotIn("1", tags)
        self.assertNotIn("7", tags)
        # The unmapped records are kept whole: FLAG 4, ZF.
        unmapped = {r[0].rsplit("-", 1)[1]: r for r in records if r[1] == "4"}
        self.assertEqual(sorted(unmapped), ["5", "6"])
        self.assertIn("ZF:Z:2,1", unmapped["5"])

        with open(os.path.join(out, "training", POINT, SAMPLE + ".taxa.tsv")) as fh:
            taxa = {r["taxon_name"]: r for r in csv.DictReader(fh, delimiter="\t")}
        self.assertEqual(list(taxa), ["s__G D", "s__G B", "s__G C"])  # FP, FN, unseen
        fp, fn, unseen = taxa["s__G D"], taxa["s__G B"], taxa["s__G C"]
        self.assertEqual((fp["error"], fp["taxid"], fp["p"], fp["knob"]), ("FP", "4", "0.8", "0.7"))
        self.assertEqual((fp["fragments_on_taxon"], fp["best_on_taxon"], fp["best_on_taxon_mapq4"]), ("1", "1", "1"))
        self.assertEqual(fp["fragments_on_taxon_from"], "s__G E (not in the database):1")
        self.assertEqual(fp["own_fragments"], "")
        self.assertEqual((fn["error"], fn["genomes"], fn["read_pairs"]), ("FN", "GCA_000000002.1", "100"))
        self.assertEqual((fn["own_fragments"], fn["own_best_on_taxon"], fn["own_best_on_taxon_mapq4"],
                          fn["own_best_elsewhere"], fn["own_unaligned"]), ("3", "1", "1", "1", "1"))
        self.assertEqual(fn["own_best_elsewhere_on"], "s__G A:1")
        self.assertEqual((fn["fragments_on_taxon"], fn["seeded_not_aligned"]), ("1", "1"))
        self.assertEqual((unseen["error"], unseen["taxid"], unseen["p"]), ("unseen", "3", ""))
        self.assertEqual((unseen["own_fragments"], unseen["own_unaligned"], unseen["fragments_on_taxon"]), ("1", "1", "0"))
        self.assertEqual(unseen["seeded_not_aligned"], "2")  # C's own read and A's read 8

    def test_drawn_reads_by_their_genome_s_place(self):
        # A PacBio sample: its reads' sources by g<i>x_, i the genome's place among its community's in the manifest.
        d = self.tmp.name
        args = self.world() + ["--calls", os.path.join(d, "trained_model_pb.calls.tsv.gz"), "--read-type", "pb",
                               "--out", os.path.join(d, "pb")]
        self.assertEqual(error_reads.main(args), 0)
        with open(os.path.join(d, "pb", "summary.tsv")) as fh:
            row = next(csv.DictReader(fh, delimiter="\t"))
        self.assertEqual((row["point"], row["sample"], row["FP"], row["FN"], row["unseen"]),
                         (LONG_POINT, LONG_SAMPLE, "1", "1", "1"))
        lines = compressed.read_text(os.path.join(d, "pb", row["sam"])).splitlines()
        tags = {line.split("\t")[0]: line.split("\txg:Z:")[1] for line in lines if not line.startswith("@")}
        self.assertEqual(tags, {"g1x_1": "GCA_000000002.1\txs:Z:s__G B\txe:Z:FN:2,source:2",
                                "g3x_2": "GCA_000000005.1\txs:Z:s__G E (not in the database)\txe:Z:FP:4",
                                "g2x_4": "GCA_000000003.1\txs:Z:s__G C\txe:Z:seeded:3,source:3"})
        with open(os.path.join(d, "pb", "training", LONG_POINT, LONG_SAMPLE + ".taxa.tsv")) as fh:
            fn = next(r for r in csv.DictReader(fh, delimiter="\t") if r["error"] == "FN")
        self.assertEqual((fn["genomes"], fn["read_pairs"], fn["own_fragments"], fn["own_best_on_taxon"]),
                         ("GCA_000000002.1", "", "1", "1"))  # no read pairs of a drawn sample's genomes

    def test_a_contig_two_genomes_share_names_neither(self):
        # An in-silico strain of an earlier insilico_strains.py kept its representative's contig names: a read of such
        # a contig has no known source.
        args = self.world()
        fasta = os.path.join(self.tmp.name, "GCA_000000005.1.fna.gz")
        with gzip.open(fasta, "wt") as fh:
            fh.write(f">{SPECIES['B'][2]} E's copy\nACGT\n>{SPECIES['E'][2]}\nACGT\n")
        out = self.run_pe(args)
        lines = compressed.read_text(os.path.join(out, "training", POINT, SAMPLE + ".sam.zst")).splitlines()
        sources = {line.split("\t")[0].rsplit("-", 1)[1]: line.split("\txs:Z:")[1].split("\t")[0]
                   for line in lines if not line.startswith("@")}
        self.assertEqual(sources["3"], "?")  # B's read 3: on B (FN), its source unknown
        self.assertNotIn("4", sources)  # B's read on A: taken for its source only, which is unknown now
        self.assertEqual(sources["2"], "s__G E (not in the database)")

    def test_the_contig_cache(self):
        # The genomes' contig names cached by path, size and time: read once for trace_relatives.py and this script.
        args = self.world()
        cache = os.path.join(self.tmp.name, "contigs.tsv.gz")
        self.run_pe(args + ["--contig-cache", cache])
        with gzip.open(cache, "rt") as fh:
            cached = {line.split("\t")[0]: line.rstrip("\n").split("\t")[2] for line in fh}
        self.assertEqual(len(cached), 4)
        self.assertEqual(cached[os.path.join(self.tmp.name, "GCA_000000002.1.fna.gz")], SPECIES["B"][2])
        # From the cache: a FASTA unchanged is not read again (its cached names stand), one changed is.
        import trace_relatives
        path = os.path.join(self.tmp.name, "GCA_000000002.1.fna.gz")
        stat = os.stat(path)
        with gzip.open(cache, "at") as fh:
            fh.write(f"{path}\t{stat.st_size}:{stat.st_mtime_ns}\tcached_name other_name\n")
        self.assertEqual(trace_relatives.genome_contigs([path], 1, cache)[path], ["cached_name", "other_name"])
        with gzip.open(path, "wt") as fh:
            fh.write(">changed\nACGTACGT\n")
        os.utime(path, ns=(stat.st_atime_ns, stat.st_mtime_ns + 10 ** 9))
        self.assertEqual(trace_relatives.genome_contigs([path], 1, cache)[path], ["changed"])


if __name__ == "__main__":
    unittest.main()
