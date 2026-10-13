#!/usr/bin/env python3
"""Tests of the benchmark on public mock communities: download_mocks.py (choosing a run's files, cutting a FASTQ while
it streams, MD5 checks, resuming), score_mocks.py (scores of hand-made profiles against each kind of truth, replicate
pairs, the protal runs of a stand-in binary) and prepare_tables.py (reading .xlsx, writing truth tables). A local HTTP
server stands in for ENA. Standard library only.

Run: python3 -m unittest scripts/mock_benchmark/test_mock_benchmark.py
"""
import csv
import functools
import gzip
import hashlib
import http.server
import io
import json
import math
import os
import stat
import subprocess
import sys
import tempfile
import threading
import unittest
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import download_mocks as dl  # noqa: E402
import prepare_tables as prep  # noqa: E402
import score_mocks as sc  # noqa: E402


def fastq(lengths, prefix="r"):
    return "".join(f"@{prefix}{i}\n{'A' * n}\n+\n{'I' * n}\n" for i, n in enumerate(lengths)).encode()


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


class Server:
    """A local HTTP server over a folder."""

    def __init__(self, folder):
        handler = functools.partial(http.server.SimpleHTTPRequestHandler, directory=folder)
        handler.log_message = lambda *a: None
        self.httpd = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
        self.url = f"http://127.0.0.1:{self.httpd.server_address[1]}/"
        self.thread = threading.Thread(target=self.httpd.serve_forever, daemon=True)
        self.thread.start()

    def stop(self):
        self.httpd.shutdown()
        self.httpd.server_close()


class DownloadTest(unittest.TestCase):
    def setUp(self):
        dl.RETRY_WAIT = 0
        self.tmp = tempfile.TemporaryDirectory()
        self.served = os.path.join(self.tmp.name, "served")
        os.makedirs(self.served)
        self.files = {
            "SRR1_1.fastq.gz": gzip.compress(fastq([150] * 10, "a")),
            "SRR1_2.fastq.gz": gzip.compress(fastq([150] * 10, "b")),
            "ERR2.fastq.gz": gzip.compress(fastq([100, 200, 300, 400, 500], "l")),
            "noeol.fastq.gz": gzip.compress(fastq([50, 60]).rstrip(b"\n")),
            "broken.fastq.gz": gzip.compress(fastq([50, 60])[:-20]),
        }
        for name, data in self.files.items():
            with open(os.path.join(self.served, name), "wb") as fh:
                fh.write(data)
        self.server = Server(self.served)
        self.out = os.path.join(self.tmp.name, "out")
        os.makedirs(self.out)

    def tearDown(self):
        self.server.stop()
        self.tmp.cleanup()

    def meta(self, run, files, reads, bases, field="fastq"):
        md5s = [hashlib.md5(self.files[f]).hexdigest() for f in files]
        sizes = [str(len(self.files[f])) for f in files]
        other = "submitted" if field == "fastq" else "fastq"
        return {"run_accession": run, "library_layout": "PAIRED", "read_count": str(reads), "base_count": str(bases),
                f"{field}_ftp": ";".join(self.server.url + f for f in files), f"{field}_md5": ";".join(md5s),
                f"{field}_bytes": ";".join(sizes), f"{other}_ftp": "", f"{other}_md5": "", f"{other}_bytes": ""}

    def kept(self, path):
        with gzip.open(path, "rb") as fh:
            return fh.read()

    def test_run_files(self):
        meta = {"run_accession": "SRR9", "fastq_ftp": "ftp.x/SRR9.fastq.gz;ftp.x/SRR9_1.fastq.gz;ftp.x/SRR9_2.fastq.gz",
                "fastq_md5": "a;b;c", "fastq_bytes": "1;2;3", "submitted_ftp": "", "submitted_md5": "",
                "submitted_bytes": ""}
        first, second = dl.run_files(meta, "pe")
        self.assertEqual(first, [("https://ftp.x/SRR9_1.fastq.gz", "b", 2)])
        self.assertEqual(second, [("https://ftp.x/SRR9_2.fastq.gz", "c", 3)])
        submitted = {"run_accession": "ERR9", "fastq_ftp": "", "fastq_md5": "", "fastq_bytes": "",
                     "submitted_ftp": "ftp.x/M_S1_R1_001.fastq.gz;ftp.x/M_S1_R2_001.fastq.gz",
                     "submitted_md5": "m1;m2", "submitted_bytes": "5;6"}
        first, second = dl.run_files(submitted, "pe")
        self.assertTrue(first[0][0].endswith("R1_001.fastq.gz") and second[0][0].endswith("R2_001.fastq.gz"))
        single = dict(submitted, submitted_ftp="ftp.x/x.pod5.tar.gz;ftp.x/ont.fastq.gz", submitted_md5="m;n",
                      submitted_bytes="1;2")
        self.assertEqual(dl.run_files(single, "ont"), ([("https://ftp.x/ont.fastq.gz", "n", 2)], []))
        with self.assertRaises(ValueError):
            dl.run_files(dict(single, submitted_ftp="ftp.x/x.pod5.tar.gz"), "ont")

    def test_copy_records_cuts_by_reads_and_bases(self):
        path = os.path.join(self.out, "a.fq.gz")
        reads, bases, more = dl.copy_records(self.server.url + "SRR1_1.fastq.gz", path, max_reads=3)
        self.assertEqual((reads, bases, more), (3, 450, True))
        self.assertEqual(self.kept(path), fastq([150] * 3, "a"))
        # the record that passes the limit is kept: 100 + 200 >= 250
        reads, bases, more = dl.copy_records(self.server.url + "ERR2.fastq.gz", path, max_bases=250)
        self.assertEqual((reads, bases, more), (2, 300, True))
        reads, bases, more = dl.copy_records(self.server.url + "ERR2.fastq.gz", path)
        self.assertEqual((reads, bases, more), (5, 1500, False))
        reads, _, _ = dl.copy_records(self.server.url + "noeol.fastq.gz", path)
        self.assertEqual(reads, 2)
        self.assertEqual(self.kept(path), fastq([50, 60]))
        with self.assertRaises(IOError):
            dl.copy_records(self.server.url + "broken.fastq.gz", path)

    def test_copy_whole_checks_md5(self):
        path = os.path.join(self.out, "w.fq.gz")
        good = hashlib.md5(self.files["ERR2.fastq.gz"]).hexdigest()
        self.assertEqual(dl.copy_whole(self.server.url + "ERR2.fastq.gz", good, path), len(self.files["ERR2.fastq.gz"]))
        os.remove(path)
        with self.assertRaises(IOError):
            dl.copy_whole(self.server.url + "ERR2.fastq.gz", "0" * 32, path)
        self.assertFalse(os.path.exists(path) or os.path.exists(path + ".part"))

    def test_fetch_sample_cuts_resumes_and_copies_whole(self):
        metas = {"SRR1": self.meta("SRR1", ["SRR1_1.fastq.gz", "SRR1_2.fastq.gz"], 10, 3000),
                 "ERR2": self.meta("ERR2", ["ERR2.fastq.gz"], 5, 1500)}
        sample = {"set": "s", "sample": "s_SRR1", "runs": "SRR1", "read_type": "pe"}
        row = dl.fetch_sample(sample, metas, self.out, max_pairs=4, max_bases=None)
        self.assertEqual((row["reads"], row["bases"], row["cut"]), (4, 1200, 1))
        first, second = row["first"].split(","), row["second"].split(",")
        self.assertEqual(self.kept(first[0]), fastq([150] * 4, "a"))
        self.assertEqual(self.kept(second[0]), fastq([150] * 4, "b"))
        with open(os.path.join(self.out, "s", "s_SRR1", "SRR1.json")) as fh:
            self.assertEqual(json.load(fh)["reads"], 4)
        self.server.stop()  # a finished run is not fetched again
        self.server = Server(self.served)
        metas["SRR1"] = dict(metas["SRR1"], fastq_ftp="http://127.0.0.1:1/gone_1.fastq.gz;http://127.0.0.1:1/gone_2.fastq.gz")
        self.assertEqual(dl.fetch_sample(sample, metas, self.out, max_pairs=4, max_bases=None), row)
        # a long-read run within the limit is copied whole, as ENA has it
        metas["ERR2"] = self.meta("ERR2", ["ERR2.fastq.gz"], 5, 1500)  # the new server's address
        long_sample = {"set": "s", "sample": "s_ERR2", "runs": "ERR2", "read_type": "ont"}
        row = dl.fetch_sample(long_sample, metas, self.out, max_pairs=4, max_bases=10_000)
        self.assertEqual((row["reads"], row["bases"], row["cut"], row["second"]), (5, 1500, 0, "-"))
        with open(row["first"], "rb") as fh:
            self.assertEqual(fh.read(), self.files["ERR2.fastq.gz"])

    def test_estimate(self):
        metas = {"SRR1": self.meta("SRR1", ["SRR1_1.fastq.gz", "SRR1_2.fastq.gz"], 10, 3000)}
        table = dl.estimate([{"set": "s", "sample": "x", "runs": "SRR1", "read_type": "pe"}], metas, 5, None)
        size = len(self.files["SRR1_1.fastq.gz"]) + len(self.files["SRR1_2.fastq.gz"])
        n, runs, bases, ena, kept = table["s"]
        self.assertEqual((n, runs, bases, ena), (1, 1, 3000, size))
        self.assertAlmostEqual(kept, size * 0.5 * 1.15)


def write_profile(path, rows, unknown):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        for species, share in rows:
            fh.write(f"GCF_{species}.1\td__Bacteria;p__P;c__C;o__O;f__F;g__G;s__{species}\t{share}\n")
        fh.write(f"?\t?\t{unknown}\n")


def write_truth(folder, name, kind, rows):
    with open(os.path.join(folder, name + ".tsv"), "w") as fh:
        fh.write(f"#name\t{name}\n#kind\t{kind}\n#source\ttest\n#basis\ttest\nspecies\tshare\tmembers\n")
        for species, share in rows:
            fh.write(f"{species}\t{share}\tm\n")


class ScoreTest(unittest.TestCase):
    """Truth T: A 0.6, B 0.3, ? 0.1 (placed: A 2/3, B 1/3). Sample s1 calls A 0.5 and C 0.1 (placed: A 5/6, C 1/6):
    TP A, FP C, FN B; L1 = 1/6 + 1/3 + 1/6 = 2/3. Sample s2 matches the truth; s1 and s2 differ by Bray-Curtis 1/3
    and share 1 of their 3 species (Jaccard 1/3)."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        root = self.tmp.name
        self.truth = os.path.join(root, "truth")
        os.makedirs(self.truth)
        write_truth(self.truth, "t", "abundance", [("s__A", 0.6), ("s__B", 0.3), ("?", 0.1)])
        write_truth(self.truth, "closed", "closed", [("s__A", "NA"), ("s__B", "NA")])
        write_truth(self.truth, "drop", "dropout", [("s__X", "NA"), ("s__Y", "0")])
        write_truth(self.truth, "negative", "negative", [])
        profiles = os.path.join(root, "out", "protal", "set1", "profiles")
        write_profile(os.path.join(profiles, "s1.profile"), [("A", 0.5), ("C", 0.1)], 0.4)
        write_profile(os.path.join(profiles, "s2.profile"), [("A", 0.6), ("B", 0.3)], 0.1)
        write_profile(os.path.join(profiles, "s3.profile"), [("D", 0.01)], 0.99)
        write_profile(os.path.join(profiles, "s4.profile"), [("X", 0.5), ("Y", 0.2)], 0.3)
        write_profile(os.path.join(profiles, "s5.profile"), [("A", 0.5), ("E", 0.2)], 0.3)
        write_profile(os.path.join(profiles, "s6.profile"), [("F", 0.5)], 0.5)
        self.samples = os.path.join(root, "samples.tsv")
        with open(self.samples, "w") as fh:
            fh.write("set\tsample\truns\tread_type\ttruth\tgroup\tnote\n")
            for name, truth, group in (("s1", "t", "g1"), ("s2", "t", "g1"), ("s3", "negative", "-"),
                                       ("s4", "drop", "-"), ("s5", "closed", "-"), ("s6", "-", "-"),
                                       ("s7", "t", "g1")):  # s7 has no profile
                fh.write(f"set1\t{name}\tR\tpe\t{truth}\t{group}\t-\n")
        self.out = os.path.join(root, "out")

    def tearDown(self):
        self.tmp.cleanup()

    def test_scores(self):
        result = subprocess.run([sys.executable, os.path.join(HERE, "score_mocks.py"), "score", "--reads", self.out,
                                 "--samples", self.samples, "--truth", self.truth], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        rows = {r["sample"]: r for r in read_tsv(os.path.join(self.out, "scores", "samples.tsv"))}
        self.assertEqual(sorted(rows), ["s1", "s2", "s3", "s4", "s5", "s6"])
        s1 = rows["s1"]
        self.assertEqual((s1["called"], s1["tp"], s1["fp"], s1["fn"]), ("2", "1", "1", "1"))
        for column, expected in (("precision", 0.5), ("recall", 0.5), ("f1", 0.5), ("recall_1pct", 0.5),
                                 ("fp_share", 1 / 6), ("l1", 2 / 3), ("bray_curtis", 1 / 3), ("unknown", 0.4),
                                 ("unknown_expected", 0.1)):
            self.assertAlmostEqual(float(s1[column]), expected, places=5, msg=column)
        self.assertEqual(s1["false_species"], "s__C")
        self.assertAlmostEqual(float(rows["s2"]["l1"]), 0.0)
        self.assertEqual((rows["s3"]["kind"], rows["s3"]["fp"], rows["s3"]["false_species"]), ("negative", "1", "s__D"))
        self.assertEqual((rows["s4"]["absent_called"], rows["s4"]["present_missed"]), ("s__Y", "-"))
        self.assertEqual((rows["s5"]["tp"], rows["s5"]["fp"], rows["s5"]["fn"], rows["s5"]["l1"]), ("1", "1", "1", "NA"))
        self.assertEqual((rows["s6"]["kind"], rows["s6"]["called"], rows["s6"]["tp"]), ("-", "1", "NA"))
        pairs = read_tsv(os.path.join(self.out, "scores", "replicates.tsv"))
        self.assertEqual(len(pairs), 1)
        self.assertAlmostEqual(float(pairs[0]["bray_curtis"]), 1 / 3, places=5)
        self.assertAlmostEqual(float(pairs[0]["jaccard"]), 1 / 3, places=5)
        species = [r for r in read_tsv(os.path.join(self.out, "scores", "species.tsv")) if r["sample"] == "s1"]
        self.assertEqual({r["species"]: r["status"] for r in species}, {"s__A": "TP", "s__B": "FN", "s__C": "FP"})
        summary = read_tsv(os.path.join(self.out, "scores", "summary.tsv"))
        abundance = next(r for r in summary if r["kind"] == "abundance")
        self.assertEqual((abundance["samples"], abundance["f1"]), ("2", "0.75"))
        dropout = next(r for r in summary if r["kind"] == "dropout")
        self.assertEqual(dropout["dropouts_called"], "1")

    def test_run_with_a_stand_in_protal(self):
        reads = os.path.join(self.tmp.name, "mocks")
        os.makedirs(os.path.join(reads, "reads"))
        with open(os.path.join(reads, "reads", "downloaded.tsv"), "w") as fh:
            fh.write("set\tsample\tread_type\tfirst\tsecond\treads\tbases\tcut\n")
            fh.write("a\ta_1\tpe\t/r/a_1_1.fq.gz\t/r/a_1_2.fq.gz\t10\t3000\t0\n")
            fh.write("a\ta_2\tont\t/r/a_2.fq.gz\t-\t5\t9000\t1\n")
            fh.write("b\tb_1\tpe\t/r/b_1_1.fq.gz\t/r/b_1_2.fq.gz\t10\t3000\t0\n")
        fake = os.path.join(self.tmp.name, "fake_protal")
        calls = os.path.join(self.tmp.name, "calls.txt")
        with open(fake, "w") as fh:
            fh.write(f"""#!{sys.executable}
import os, sys
args = sys.argv[1:]
with open({calls!r}, "a") as log:
    log.write(" ".join(args) + "\\n")
lines = open(args[args.index("--map") + 1]).read().splitlines()
out = lines[0].split("\\t")[1]
os.makedirs(os.path.join(out, "profiles"), exist_ok=True)
for line in lines[2:]:
    sample = line.split("\\t")[0]
    open(os.path.join(out, "profiles", sample + ".profile"), "w").write("?\\t?\\t1\\n")
""")
        os.chmod(fake, os.stat(fake).st_mode | stat.S_IEXEC)
        command = [sys.executable, os.path.join(HERE, "score_mocks.py"), "run", "--reads", reads, "--db", "/db",
                   "--protal", fake, "-t", "3", "--sets", "a"]
        result = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        with open(os.path.join(reads, "protal", "a", "map.tsv")) as fh:
            lines = fh.read().splitlines()
        self.assertEqual(lines[1], "#SAMPLEID\tPREFIX\tFIRST\tSECOND\tREAD_TYPE")
        self.assertEqual(lines[2:], ["a_1\ta_1\t/r/a_1_1.fq.gz\t/r/a_1_2.fq.gz\tpe", "a_2\ta_2\t/r/a_2.fq.gz\t-\tont"])
        with open(calls) as fh:
            called = fh.read().splitlines()
        self.assertEqual(len(called), 1)
        self.assertIn("--db /db", called[0])
        self.assertIn("-t 3 --no_strains --force", called[0])
        self.assertTrue(os.path.exists(os.path.join(reads, "protal", "a", "profiles", "a_2.profile")))
        result = subprocess.run(command, capture_output=True, text=True)  # profiles there: protal is not run again
        self.assertEqual(result.returncode, 0, result.stderr)
        with open(calls) as fh:
            self.assertEqual(len(fh.read().splitlines()), 1)
        # two chunks of set b: chunk 1 has no sample (b has one), chunk 0 runs it into its own folder
        for chunk in ("1", "0"):
            result = subprocess.run(command[:-1] + ["b", "--chunks", "2", "--chunk", chunk], capture_output=True,
                                    text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(os.path.exists(os.path.join(reads, "protal", "b", "chunk_0", "profiles", "b_1.profile")))


class PrepareTest(unittest.TestCase):
    def test_read_xlsx(self):
        buffer = io.BytesIO()
        main = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
        rels = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
        with zipfile.ZipFile(buffer, "w") as z:
            z.writestr("xl/workbook.xml", f'<workbook xmlns="{main}" xmlns:r="{rels}"><sheets>'
                                          f'<sheet name="S1" sheetId="1" r:id="rId1"/></sheets></workbook>')
            z.writestr("xl/_rels/workbook.xml.rels",
                       '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
                       '<Relationship Id="rId1" Target="worksheets/sheet1.xml"/></Relationships>')
            z.writestr("xl/sharedStrings.xml", f'<sst xmlns="{main}"><si><t>name</t></si><si><t>E. coli</t></si></sst>')
            z.writestr("xl/worksheets/sheet1.xml",
                       f'<worksheet xmlns="{main}"><sheetData>'
                       '<row r="1"><c r="A1" t="s"><v>0</v></c><c r="C1" t="inlineStr"><is><t>share</t></is></c></row>'
                       '<row r="2"><c r="A2" t="s"><v>1</v></c><c r="C2"><v>0.25</v></c></row>'
                       '</sheetData></worksheet>')
        self.assertEqual(prep.read_xlsx(buffer.getvalue()), {"S1": [["name", "", "share"], ["E. coli", "", "0.25"]]})

    def test_write_truth(self):
        with tempfile.TemporaryDirectory() as folder:
            prep.write_truth(folder, "x", "abundance", "src", "genome copies",
                             [("a1", "s__A", 3.0), ("a2", "s__A", 1.0), ("b", "s__B", 0.5), ("y", None, 0.5)])
            with open(os.path.join(folder, "x.tsv")) as fh:
                lines = fh.read().splitlines()
            self.assertEqual(lines[:2], ["#name\tx", "#kind\tabundance"])
            self.assertEqual(lines[4:], ["species\tshare\tmembers", "s__A\t0.8\ta1; a2", "s__B\t0.1\tb", "?\t0.1\ty"])
            kind, rows, unknown = sc.read_truth(os.path.join(folder, "x.tsv"))
            self.assertEqual((kind, rows, unknown), ("abundance", {"s__A": 0.8, "s__B": 0.1}, 0.1))
            prep.write_truth(folder, "d", "dropout", "src", "dropout", [("x", "s__X", None), ("y", "s__Y", 0)])
            kind, rows, _ = sc.read_truth(os.path.join(folder, "d.tsv"))
            self.assertTrue(math.isnan(rows["s__X"]))
            self.assertEqual((kind, rows["s__Y"]), ("dropout", 0.0))

    def test_read_types_and_microbench_names(self):
        def run(model, layout="PAIRED"):
            return {"instrument_model": model, "library_layout": layout}
        self.assertEqual([prep.read_type(run(m)) for m in ("Illumina NovaSeq 6000", "PromethION", "Sequel II",
                                                           "Ion Torrent S5", "DNBSEQ-T7")],
                         ["pe", "ont", "pb", "se", "pe"])
        self.assertEqual(prep.read_type(run("Illumina HiSeq 2000", "SINGLE")), "se")
        self.assertEqual(prep.microbench_sample({"sample_alias": "ZymoGut metagenome", "sample_title": "x"}), "ZymoGut")
        self.assertEqual(prep.microbench_sample({"sample_alias": "", "sample_title": "ZymoHMW"}), "ZymoHMW")


if __name__ == "__main__":
    unittest.main()
