#!/usr/bin/env python3
"""Tests of the strain test's report scripts (scripts/strain_test, `just strain-dbcounts` and `just strain-report`):
db_gene_counts.py (each species' marker genes in the database) and strain_report.py (the QC report of a strain run and
qcmsa), on a tiny hand-made strain output whose expected numbers and checks are worked out in the comments. Standard
library only.

Run: python3 -m unittest scripts/test_strain_scripts.py
"""

import csv
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
DB_GENE_COUNTS = os.path.join(HERE, "strain_test", "db_gene_counts.py")
STRAIN_REPORT = os.path.join(HERE, "strain_test", "strain_report.py")
ALPHA, GAMMA, NOVEL = "s__Mockella_alpha", "s__Fakibacter_gamma", "s__Novel_x"


def write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)


def run(script, *args):
    return subprocess.run([sys.executable, script, *args], capture_output=True, text=True)


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


class StrainScriptsTest(unittest.TestCase):
    """A strain run of three species. Mockella alpha: 3 samples x 2 genes in its meta table (gene 1 covered in all,
    gene 2 in s1 only), protal's MSA of the reference and the 3 samples over 2 genes, every variant position of its 2
    samples' SNP stats filtered, and qcmsa's output (2 samples, 1 gene, 5 of 1000 sites). Fakibacter gamma and Novel x:
    a meta table, no MSA. The database: alpha with two genomes of 5 and 3 marker rows, gamma with one of 2 rows, another
    species; Novel x is not in it."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        d = self.tmp.name
        self.strains, self.db = os.path.join(d, "run", "strains"), os.path.join(d, "db")
        lineage = "d__Bacteria;p__P;c__C;o__O;f__F;g__{};s__{}"
        write(os.path.join(self.db, "genome2tiid.tsv"),
              f"GCF_1.1\t1\tGCF_1.1\t{lineage.format('Mockella', 'Mockella alpha')}\n"
              f"GCA_2.1\t2\tGCF_1.1\t{lineage.format('Mockella', 'Mockella alpha')}\n"
              f"GCF_3.1\t3\tGCF_3.1\t{lineage.format('Fakibacter', 'Fakibacter gamma')}\n"
              f"GCF_4.1\t4\tGCF_4.1\t{lineage.format('Otherella', 'Otherella one')}\n")
        write(os.path.join(self.db, "unique_kmers.tsv"),
              "".join(f"{tiid}\t{gene}\t100\t50\n" for tiid, genes in ((1, 5), (2, 3), (3, 2), (4, 7)) for gene in range(genes)))
        meta = "sample\tgene_id\thcov\tmean_vcov_nonzero\tvertical_coverage\tmulti_rate_vcov2\n"
        write(os.path.join(self.strains, ALPHA + ".meta.tsv"), meta + "".join(
            f"{s}\t{g}\t{0.9 if g == 1 or s == 's1' else 0.2}\t5\t4\t0.01\n" for s in ("s1", "s2", "s3") for g in (1, 2)))
        for species in (GAMMA, NOVEL):
            write(os.path.join(self.strains, species + ".meta.tsv"), meta + "s1\t3\t0.9\t5\t4\t0\n")
        write(os.path.join(self.strains, ALPHA + ".raw.msa.fna"),
              "".join(f">{name}\nACGTACGT\n" for name in (ALPHA + "_reference", "s1", "s2", "s3")))
        write(os.path.join(self.strains, ALPHA + ".raw.partition.txt"), "DNA, 1 = 1-4\nDNA, 2 = 5-8\n")
        stats = "\t".join(("sample", "snps_retained", "variants_filtered_qual_sum", "variants_filtered_obs_cov",
                           "positions_below_min_cov", "positions_no_coverage", "total_variant_positions")) + "\n"
        write(os.path.join(self.strains, ALPHA + ".snp_stats.tsv"),
              stats + f"{ALPHA}_reference\t9\t0\t0\t0\t0\t9\ns1\t0\t3\t2\t1\t0\t5\ns2\t0\t2\t1\t0\t4\t3\n")
        write(os.path.join(self.strains, ALPHA + ".qcmsa_summary.tsv"),
              "section\tkey\tvalue\treason\nparam\tmax_mrate\t0.1\t\ncount\tsites_in\t1000\t\ncount\tsites_kept\t5\t\n"
              "count\tgenes_filtered\t1\t\ncount\tsamples_filtered\t1\t\ncount\toutlier_cells\t0\t\n"
              "gene_filtered\t2\t0.4\tmrate\nsample_filtered\ts3\t0.3\tmrate\n")
        write(os.path.join(self.strains, ALPHA + ".msa.fna"),
              "".join(f">{name}\nACGT\n" for name in (ALPHA + "_reference", "s1", "s2")))
        write(os.path.join(self.strains, ALPHA + ".partition.txt"), "DNA, 1 = 1-4\n")
        # protal's run log with the M3 thresholds: gene 1 passes in 3 samples, more than 1; gene 2 in 1, not more.
        write(os.path.join(self.tmp.name, "run", "protal_run.log"),
              "gene min hcov frac: 0.5\ngene min mean depth: 1\nmsa min samples: 1\n")
        self.counts = os.path.join(self.tmp.name, "run", "db_gene_counts.tsv")

    def test_db_gene_counts(self):
        # Each species of the strain run (spaces and underscores alike): its genome with the most rows of
        # unique_kmers.tsv, one row per marker gene; a species the database lacks gets empty fields.
        result = run(DB_GENE_COUNTS, "--db", self.db, "--strains", self.strains, "--out", self.counts)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual([(r["species"], r["genome_id"], r["db_markers"]) for r in read_tsv(self.counts)],
                         [(GAMMA, "3", "2"), (ALPHA, "1", "5"), (NOVEL, "", "")])
        self.assertIn("(3 species)", result.stderr)
        # A database packed into database.protal has to be unpacked first.
        os.remove(os.path.join(self.db, "unique_kmers.tsv"))
        write(os.path.join(self.db, "database.protal"), "")
        result = run(DB_GENE_COUNTS, "--db", self.db, "--strains", self.strains, "--out", self.counts)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("protal --unpack_db", result.stderr)

    def test_report(self):
        self.assertEqual(run(DB_GENE_COUNTS, "--db", self.db, "--strains", self.strains, "--out", self.counts).returncode, 0)
        out = os.path.join(self.tmp.name, "report")
        result = run(STRAIN_REPORT, "--strains", self.strains, "--out", out)
        self.assertEqual(result.returncode, 0, result.stderr)
        rows = {r["species"]: r for r in read_tsv(os.path.join(out, "summary.tsv"))}
        self.assertEqual(sorted(rows), sorted((ALPHA, GAMMA, NOVEL)))
        alpha = rows[ALPHA]
        # The reference row is no sample; the counts before (protal) and after qcmsa, and qcmsa's own.
        self.assertEqual({k: alpha[k] for k in ("meta_samples", "meta_genes", "before_seqs", "before_sample_seqs",
                                                "before_genes", "after_seqs", "after_sample_seqs", "after_genes",
                                                "sites_in", "sites_kept", "genes_filtered", "samples_filtered",
                                                "outlier_cells", "qcmsa_ran")},
                         {"meta_samples": "3", "meta_genes": "2", "before_seqs": "4", "before_sample_seqs": "3",
                          "before_genes": "2", "after_seqs": "3", "after_sample_seqs": "2", "after_genes": "1",
                          "sites_in": "1000", "sites_kept": "5", "genes_filtered": "1", "samples_filtered": "1",
                          "outlier_cells": "0", "qcmsa_ran": "True"})
        self.assertEqual((rows[GAMMA]["before_seqs"], rows[GAMMA]["qcmsa_ran"]), ("None", "False"))
        with open(os.path.join(out, "report.md")) as fh:
            report = fh.read()
        self.assertIn("- Species with a protal MSA: **1**", report)
        # Alpha: all 8 variant positions of its 2 samples filtered (the reference's row is no sample); 2 of the 5 marker
        # genes of its genome hit; the M3 thresholds of the run log drop gene 2; qcmsa keeps 5 of 1000 sites and 3
        # sequences. Gamma and Novel x: no MSA.
        expected = [f"[{ALPHA}]: ALL 8 variant positions filtered by SNP filters (M1) across 2 samples",
                    f"[{ALPHA}]: only 2/5 marker genes in the DB genome got any read hits (3 unhit, 60%)",
                    f"[{ALPHA}]: protal M3 gene-coverage filter dropped 1/2 observed genes (50%)",
                    f"[{ALPHA}]: qcmsa site cleanup kept only 5/1000 sites (0.5%)",
                    f"[{ALPHA}]: after qcmsa only 3 sequences remain (<4)",
                    f"[{GAMMA}]: no MSA produced by protal", f"[{NOVEL}]: no MSA produced by protal"]
        self.assertIn(f"## Automated checks: 0 FAIL, {len(expected)} WARN", report)
        for line in expected:
            self.assertIn(f"- **WARN** {line}", report)
        self.assertIn("(0 FAIL, 7 WARN)", result.stderr)
        with open(os.path.join(out, "report.html")) as fh:
            page = fh.read()
        self.assertIn("<svg", page)  # the plots and alpha's heatmaps, inline
        self.assertIn(ALPHA, page)
        # qcmsa keeping more genes than protal's MSA had is a failure, and the exit code says so; thresholds given on
        # the command line replace the log's (no gene dropped at 0).
        write(os.path.join(self.strains, ALPHA + ".partition.txt"), "DNA, 1 = 1-4\nDNA, 2 = 5-8\nDNA, 3 = 9-12\n")
        result = run(STRAIN_REPORT, "--strains", self.strains, "--out", out, "--gene-min-hcov-frac", "0")
        self.assertEqual(result.returncode, 1, result.stderr)
        with open(os.path.join(out, "report.md")) as fh:
            report = fh.read()
        self.assertIn(f"- **FAIL** [{ALPHA}]: qcmsa kept MORE genes (3) than input (2)", report)
        self.assertNotIn("M3 gene-coverage filter dropped", report)


if __name__ == "__main__":
    unittest.main()
