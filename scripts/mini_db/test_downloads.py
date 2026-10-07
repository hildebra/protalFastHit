#!/usr/bin/env python3
"""test_downloads.py - checks for download_gtdb.py and download_gtdb_releases.py, offline.

The synthetic GTDB release (mini_db_fixtures.shared_release) is served by a local HTTP server laid out as GTDB's mirror
and NCBI's FTP server are, and NCBI's `datasets` is a stand-in script: nothing goes to the internet (every download is
given --mirror, --ftp_url and --datasets). Standard library only.

  python3 -m unittest scripts/mini_db/test_downloads.py
"""

import argparse
import collections
import contextlib
import csv
import gzip
import hashlib
import http.server
import io
import json
import os
import random
import shutil
import subprocess
import sys
import tempfile
import threading
import time
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mini_db_fixtures import (CONVERT, DOWNLOAD, RELEASES_DOWNLOAD, FakeNcbi, assembly_folder, fake_datasets,  # noqa: E402
                              run, shared_release)

import download_gtdb as dl  # noqa: E402


class DownloadTest(unittest.TestCase):
    """download_gtdb.py and download_gtdb_releases.py on the synthetic release, served by the stand-ins."""

    @classmethod
    def setUpClass(cls):
        cls.gtdb, cls.db = shared_release()
        cls.tmp = tempfile.TemporaryDirectory()
        cls.ncbi = FakeNcbi(cls.gtdb, os.path.join(cls.tmp.name, "served"))

    @classmethod
    def tearDownClass(cls):
        cls.ncbi.stop()
        cls.tmp.cleanup()

    def test_download_gtdb(self):
        ncbi = self.ncbi
        # GCA_999002003.1 is the PacBio assembly among the strains of Mockella beta, so it is the one picked,
        # and NCBI no longer has it.
        env = dict(ncbi.env, FAKE_SUPPRESSED="GCA_999002003.1", FAKE_LONG="GCA_999002003.1")
        out = os.path.join(self.tmp.name, "inputs")
        command = [sys.executable, DOWNLOAD, "-o", out, *ncbi.download_options, "--species", "3", "--per_species", "1",
                   "--rep_only_species", "0", "--batch", "4", "-t", "2", "--host_genome", "none"]
        first = subprocess.run(command, env=env, capture_output=True, text=True)
        self.assertEqual(first.returncode, 0, first.stdout + first.stderr)

        release = os.path.join(out, "release")
        for path in ("bac120_taxonomy_r226.tsv.gz", "bac120_metadata_r226.tsv.gz", "VERSION.txt",
                     "genomic_files_reps/bac120_marker_genes_reps_r226", "genomic_files_all/bac120_marker_genes_all_r226"):
            self.assertTrue(os.path.exists(os.path.join(release, path)), path)
        self.assertFalse(os.path.exists(os.path.join(release, "genomic_files_reps", "bac120_marker_genes_reps_r226.tar.gz")),
                         "archives are removed once extracted")
        genomes = sorted(os.listdir(os.path.join(out, "genomes")))
        # 3 species with a strain each and their 3 representatives, less the suppressed strain
        self.assertEqual(len(genomes), 5, genomes)
        self.assertEqual(sum(g.startswith("GCF_") for g in genomes), 3)
        with gzip.open(os.path.join(out, "genomes", genomes[0]), "rt") as fh:
            self.assertTrue(fh.readline().startswith(">"))
        with open(os.path.join(out, "missing.txt")) as fh:
            self.assertEqual(fh.read().split(), ["GCA_999002003.1"])
        with open(os.path.join(out, "simulation_species.txt")) as fh:
            self.assertEqual(len(fh.read().split("\n")) - 1, 3)
        with open(os.path.join(out, "download.json")) as fh:
            state = json.load(fh)
        self.assertEqual(state["release"]["version"], "226.0")
        self.assertEqual((state["genomes"]["delivered"], state["genomes"]["missing"]), (5, 1))
        self.assertEqual(state["genomes"]["strain_category"], {"isolate": 2})
        with open(os.path.join(out, "genomes.tsv")) as fh:
            header = next(fh).rstrip("\n").split("\t")
            table = [dict(zip(header, line.rstrip("\n").split("\t"))) for line in fh]
        self.assertEqual({r["role"] for r in table}, {"strain", "representative"})
        self.assertEqual({r["genome_category"] for r in table}, {"isolate"})
        self.assertTrue(all(r["sequencing_tech"] == "Illumina HiSeq" for r in table if r["role"] == "strain"), table)
        with open(os.path.join(out, "ncbi_info.tsv")) as fh:  # every candidate strain was asked for, once
            asked = [line.split("\t")[0] for line in fh][1:]
        self.assertIn("GCA_999002003.1", asked)
        self.assertEqual(len(asked), len(set(asked)))

        # A rerun downloads nothing it has and asks NCBI nothing it knows, and the release converts as downloaded.
        # A genome of an earlier choice is moved out of genomes/, which a build reads whole.
        with gzip.open(os.path.join(out, "genomes", "GCA_999999999.1.fna.gz"), "wt") as fh:
            fh.write(">stray\nACGT\n")
        calls = os.path.join(self.tmp.name, "datasets_rerun_calls.txt")
        again = subprocess.run(command, env=dict(env, FAKE_LOG=calls), capture_output=True, text=True)
        self.assertEqual(again.returncode, 0, again.stdout + again.stderr)
        self.assertIn("1 genomes that this choice does not include moved to", again.stdout)
        self.assertEqual(os.listdir(os.path.join(out, "genomes_unused")), ["GCA_999999999.1.fna.gz"])
        self.assertNotIn("GCA_999999999.1.fna.gz", os.listdir(os.path.join(out, "genomes")))
        self.assertIn("marker_genes_reps_r226.tar.gz: already there", again.stdout)
        self.assertIn("5 already there, 1 to download", again.stdout)
        self.assertIn(f"sequencing technology of {len(asked)} candidate strains: {len(asked)} known, 0 to ask", again.stdout)
        with open(calls) as fh:
            self.assertNotIn("summary", fh.read())
        # Genomes NCBI no longer has (withdrawn): datasets packs nothing to fetch for them. A rerun asking only for such
        # genomes took that for NCBI being down (r226, 2026-10-05: 52 such genomes, 13 failed requests in a row, no
        # host genome fetched); they are missing, and the rest of the run goes on.
        gone = os.path.join(self.tmp.name, "inputs_withdrawn")
        shutil.copytree(out, gone)
        lost = sorted(os.listdir(os.path.join(gone, "genomes")))[:4]
        for name in lost:
            os.remove(os.path.join(gone, "genomes", name))
        withdrawn = [name[:-len(".fna.gz")] for name in lost] + ["GCA_999002003.1"]
        rerun = [gone if c == out else c for c in command]
        rerun[rerun.index("--batch") + 1] = "1"  # 5 requests of one, as many as stop a run that fails them all
        result = subprocess.run(rerun, env=dict(env, FAKE_SUPPRESSED="", FAKE_WITHDRAWN=",".join(withdrawn)),
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertNotIn("failed in a row", result.stderr)
        self.assertIn("Inputs for GTDB r226", result.stdout)
        with open(os.path.join(gone, "missing.txt")) as fh:
            self.assertEqual(sorted(fh.read().split()), sorted(withdrawn))
        run(CONVERT, "--gtdb", release, "--outdir", os.path.join(self.tmp.name, "db_downloaded"))
        with open(os.path.join(self.db, "reference.fna"), "rb") as a, \
                open(os.path.join(self.tmp.name, "db_downloaded", "reference.fna"), "rb") as b:
            self.assertEqual(a.read(), b.read())

        missing = subprocess.run([sys.executable, DOWNLOAD, "-o", out + "_x", "--release", "999", "--mirror", ncbi.url,
                                  "--ftp_url", ncbi.ftp_url, "--no_genomes"], capture_output=True, text=True)
        self.assertNotEqual(missing.returncode, 0)

        # With NCBI down, a few failed requests in a row stop the run (not ~2n of them halving every batch),
        # in the technology lookup or, without it, in the download, and it does not say its inputs are ready.
        for extra, failed in (([], "datasets summary failed (1)"), (["--no_tech_lookup"], "datasets download failed (1)")):
            calls = os.path.join(self.tmp.name, "datasets_calls.txt")
            if os.path.exists(calls):
                os.remove(calls)
            down = subprocess.run([sys.executable, DOWNLOAD, "-o", out + "_down", *ncbi.download_options, "--species", "3",
                                   "--per_species", "1", "--rep_only_species", "0", "--batch", "1", "-t", "2", *extra],
                                  env=dict(env, FAKE_DOWN="1", FAKE_LOG=calls), capture_output=True, text=True)
            self.assertNotEqual(down.returncode, 0, down.stdout)
            self.assertIn(f"NCBI requests failed in a row, the last: {failed}: Error: Gateway Timeout", down.stderr)
            self.assertNotIn("Inputs for GTDB", down.stdout)
            with open(calls) as fh:
                self.assertEqual(len(fh.readlines()), 5)  # --batch 1: 1 .bit_length() + 4

    def test_download_releases(self):
        # download_gtdb_releases.py: download_gtdb.py per release into INPUTS/gtdb_r<release>, the options it does
        # not know passed on, and a summary of what each release's inputs hold.
        ncbi = self.ncbi
        out = os.path.join(self.tmp.name, "inputs_releases")
        result = subprocess.run([sys.executable, RELEASES_DOWNLOAD, "-o", out, "--releases", "226", "-t", "2",
                                 *ncbi.download_options, "--species", "2", "--per_species", "1", "--rep_only_species", "0",
                                 "--batch", "4", "--host_genome", "none"],
                                env=ncbi.env, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertTrue(os.path.isfile(os.path.join(out, "gtdb_r226", "download.json")))
        self.assertTrue(os.path.isfile(os.path.join(out, "download_r226.log")))
        with open(os.path.join(out, "download_summary.tsv")) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            rows = [dict(zip(header, line.rstrip("\n").split("\t"))) for line in fh]
        self.assertEqual(len(rows), 1)
        row = rows[0]
        self.assertEqual((row["release"], row["version"], row["status"]), ("226", "226.0", "ok"))
        self.assertEqual(row["folder"], os.path.join(out, "gtdb_r226"))
        with open(os.path.join(self.gtdb, "bac120_taxonomy_r226.tsv")) as fh:
            lineages = [line.rstrip("\n").split("\t")[1] for line in fh]
        self.assertEqual(int(row["genomes_in_taxonomy"]), len(lineages))
        self.assertEqual(int(row["species"]), len(set(lineages)))
        self.assertEqual(int(row["genera"]), len({l.rsplit(";", 1)[0] for l in lineages}))
        self.assertEqual(int(row["phyla"]), len({";".join(l.split(";")[:2]) for l in lineages}))
        self.assertEqual(int(row["bacteria_species"]) + int(row["archaea_species"]), int(row["species"]))
        self.assertEqual(int(row["genomes_delivered"]), 4)  # 2 species, a strain and the representative each
        self.assertEqual(row["genomes_missing"], "0")
        self.assertEqual(row["simulation_species"], "2")
        self.assertGreater(int(row["marker_files"]), 0)
        self.assertRegex(row["size_gb"], r"^\d+\.\d$")
        self.assertRegex(row["download_time"], r"^\d+:\d\d:\d\d$")
        self.assertIn("Summary: " + os.path.join(out, "download_summary.tsv"), result.stdout)
        self.assertRegex(result.stdout, r"\nrelease  version  status  genomes_in_taxonomy")
        failing = subprocess.run([sys.executable, RELEASES_DOWNLOAD, "-o", out + "_x", "--releases", "999",
                                  *ncbi.download_options, "--no_genomes"], env=ncbi.env, capture_output=True, text=True)
        self.assertEqual(failing.returncode, 1)
        self.assertIn("failed: r999", failing.stdout)

    def test_download_direct(self):
        # Genomes come straight from the (stand-in) FTP server, in parallel, without datasets; the ones it does
        # not have go through datasets.
        ncbi = FakeNcbi(self.gtdb, os.path.join(self.tmp.name, "served_direct"), self, assembly_names=True)
        mirror = ncbi.mirror
        with open(os.path.join(self.gtdb, "simulation", "genomes.tsv")) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                a = row["accession"]
                folder = os.path.join(mirror, "ftp", a[:3], a[4:7], a[7:10], a[10:13], assembly_folder(a))
                os.makedirs(folder)
                shutil.copy(row["fasta_path"], os.path.join(folder, assembly_folder(a) + "_genomic.fna.gz"))
        # The host genome (the human one by default) where NCBI keeps it, a small stand-in.
        host = os.path.join(mirror, "ftp", "GCF", "009", "914", "755", "GCF_009914755.1_T2T-CHM13v2.0")
        os.makedirs(host)
        with gzip.open(os.path.join(host, "GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz"), "wt") as fh:
            fh.write(">NC_060925.1 chromosome 1\n" + "ACGTTGCA" * 400 + "\n")
        datasets = ncbi.datasets
        env = ncbi.env

        def download(out, *extra, **more):
            return subprocess.run([sys.executable, DOWNLOAD, "-o", out, "--mirror", ncbi.url, "--ftp_url", ncbi.ftp_url,
                                   "--species", "3", "--per_species", "1", "--rep_only_species", "0", "--connections", "4",
                                   *extra], env=dict(env, **more), capture_output=True, text=True)

        # Without datasets at all (and no lookup of the sequencing technology, which needs it).
        out = os.path.join(self.tmp.name, "inputs_direct")
        first = download(out, "--datasets", "/nonexistent/datasets", "--no_tech_lookup")
        self.assertEqual(first.returncode, 0, first.stdout + first.stderr)
        self.assertIn("fetching 6 genomes from", first.stdout)
        self.assertIn("6 fetched directly", first.stdout)
        genomes = sorted(os.listdir(os.path.join(out, "genomes")))
        self.assertEqual(len(genomes), 6, genomes)
        self.assertFalse([g for g in genomes if g.endswith(".part")])
        with open(os.path.join(out, "genomes.tsv")) as fh:
            rows = [r for r in csv.DictReader(fh, delimiter="\t")]
        with open(os.path.join(self.gtdb, "simulation", "genomes.tsv")) as fh:
            source = {r["accession"]: r["fasta_path"] for r in csv.DictReader(fh, delimiter="\t")}
        for r in rows:  # the file is NCBI's, byte for byte
            with open(os.path.join(out, "genomes", r["accession"] + ".fna.gz"), "rb") as a, open(source[r["accession"]], "rb") as b:
                self.assertEqual(a.read(), b.read())
        with open(os.path.join(out, "download.json")) as fh:
            state = json.load(fh)
        self.assertEqual(state["genomes"]["fetched_directly"], 6)
        # The host genome, gzipped as served, for the scenarios with host reads (build_gtdb_database.py finds it).
        self.assertEqual(state["host"]["path"], os.path.join("host", "GCF_009914755.1.fna.gz"))
        with gzip.open(os.path.join(out, state["host"]["path"]), "rt") as fh:
            self.assertTrue(fh.readline().startswith(">NC_060925.1"))
        again = download(out, "--datasets", "/nonexistent/datasets", "--no_tech_lookup")
        self.assertEqual(again.returncode, 0, again.stdout + again.stderr)
        self.assertIn("6 already there, 0 to download", again.stdout)
        self.assertIn("host genome: " + os.path.join(out, "host", "GCF_009914755.1.fna.gz") + " is there", again.stdout)
        # A host the server lacks, without datasets, is only noted: the GTDB inputs serve every build without it.
        missing_host = download(out + "_nohost", "--datasets", "/nonexistent/datasets", "--no_tech_lookup",
                                "--host_genome", "GCF_000000001.1_none")
        self.assertEqual(missing_host.returncode, 0, missing_host.stdout + missing_host.stderr)
        self.assertIn("host genome GCF_000000001.1: not downloaded (not found)", missing_host.stdout)
        with open(os.path.join(out + "_nohost", "download.json")) as fh:
            self.assertNotIn("host", json.load(fh))

        # A genome the server lacks (404) goes through datasets, and the others do not.
        lacking = next(r["accession"] for r in rows if r["role"] == "strain")
        a = lacking
        os.remove(os.path.join(mirror, "ftp", a[:3], a[4:7], a[7:10], a[10:13], assembly_folder(a),
                               assembly_folder(a) + "_genomic.fna.gz"))
        calls = os.path.join(self.tmp.name, "datasets_direct_calls.txt")
        out2 = os.path.join(self.tmp.name, "inputs_direct2")
        second = download(out2, "--datasets", datasets, "--no_tech_lookup", FAKE_LOG=calls)
        self.assertEqual(second.returncode, 0, second.stdout + second.stderr)
        self.assertIn("5 fetched directly; 1 left for datasets (1 not found)", second.stdout)
        self.assertEqual(len(os.listdir(os.path.join(out2, "genomes"))), 6)
        with open(calls) as fh:
            asked = fh.read().splitlines()
        self.assertEqual(len([c for c in asked if c.startswith("download genome accession")]), 1, asked)
        self.assertEqual(len([c for c in asked if c.startswith("rehydrate")]), 1, asked)
        self.assertEqual(len(asked), 2, asked)


class DownloadFunctionsTest(unittest.TestCase):
    """download_gtdb.py's parts on their own: fetching a file, choosing the strains, reading GTDB's metadata, asking
    NCBI."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def serve_handler(self, handler):
        server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
        threading.Thread(target=server.serve_forever, daemon=True).start()
        self.addCleanup(server.server_close)
        self.addCleanup(server.shutdown)
        return f"http://127.0.0.1:{server.server_port}"

    def test_fetch_file(self):
        self.assertEqual(dl.ftp_url("http://x", "GCF_000005845.2", "ASM584v2"),
                         "http://x/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz")
        self.assertEqual(dl.ftp_url("http://x", "GCA_947500805.1", "ATCC 21022 (v1)#2"),
                         "http://x/GCA/947/500/805/GCA_947500805.1_ATCC_21022__v1__2/GCA_947500805.1_ATCC_21022__v1__2_genomic.fna.gz")
        good = gzip.compress(b">g\nACGT\n" * 1000)
        requests = collections.Counter()

        class Handler(http.server.BaseHTTPRequestHandler):
            protocol_version = "HTTP/1.1"

            def log_message(self, *args):
                pass

            def do_GET(self):
                name = self.path.lstrip("/")
                requests[name] += 1
                if name.startswith("gone"):
                    self.send_error(404)
                    return
                if name.startswith("forbidden"):
                    self.send_error(403)
                    return
                if name.startswith("busy") or (name.startswith("flaky") and requests[name] <= 2):
                    self.send_error(503)
                    return
                if name.startswith("later") and requests[name] == 1:
                    self.send_response(503)
                    self.send_header("Retry-After", "1")
                    self.send_header("Content-Length", "0")
                    self.end_headers()
                    return
                body = good
                if name.startswith("bad"):
                    body = b"this is not gzip" * 50
                self.send_response(200)
                self.send_header("Content-Length", str(len(body) + (500 if name.startswith("short") else 0)))
                self.end_headers()
                self.wfile.write(body)
                self.close_connection = True

        base = self.serve_handler(Handler)
        folder = os.path.join(self.tmp.name, "fetch_file")
        os.makedirs(folder)
        wait, dl.RETRY_WAIT = dl.RETRY_WAIT, 0.0
        self.addCleanup(setattr, dl, "RETRY_WAIT", wait)
        dest = os.path.join(folder, "x.fna.gz")
        self.assertEqual(dl.fetch_file(f"{base}/ok.gz", dest), "")
        with open(dest, "rb") as fh:
            self.assertEqual(fh.read(), good)
        self.assertEqual(dl.fetch_file(f"{base}/flaky.gz", dest), "")  # 503 twice, then the file
        self.assertEqual(requests["flaky.gz"], 3)
        self.assertEqual(dl.fetch_file(f"{base}/gone.gz", dest), "not found")
        self.assertEqual(requests["gone.gz"], 1, "a missing file is not asked for again")
        self.assertEqual(dl.fetch_file(f"{base}/forbidden.gz", dest), "HTTP 403")
        self.assertEqual(requests["forbidden.gz"], 1)
        self.assertEqual(dl.fetch_file(f"{base}/busy.gz", dest), "HTTP 503")
        self.assertEqual(requests["busy.gz"], dl.ATTEMPTS)
        started = time.time()  # the server says when to come back
        self.assertEqual(dl.fetch_file(f"{base}/later.gz", dest), "")
        self.assertGreaterEqual(time.time() - started, 1.0)
        self.assertEqual(requests["later.gz"], 2)
        self.assertIn("BadGzipFile", dl.fetch_file(f"{base}/bad.gz", dest))  # the right length, not gzip
        short = dl.fetch_file(f"{base}/short.gz", dest)
        self.assertTrue(short, "a file shorter than announced is no delivery")
        self.assertEqual(sorted(os.listdir(folder)), ["x.fna.gz"], "no .part file is left behind")

        # NCBI not answering: after FAILS_IN_A_ROW failures the rest is left for datasets, not asked for.
        opts = argparse.Namespace(connections=1)
        before = sum(requests.values())
        urls = [(f"G{i}", f"{base}/busy{i}.gz") for i in range(dl.FAILS_IN_A_ROW + 5)]
        delivered, failed = dl.fetch_direct(opts, urls, folder)
        self.assertEqual((len(delivered), len(failed)), (0, len(urls)))
        self.assertEqual(sum(1 for why in failed.values() if why.startswith("not asked for")), 5)
        self.assertEqual(sum(requests.values()) - before, dl.ATTEMPTS * dl.FAILS_IN_A_ROW)
        # Missing files are no sign of trouble, however many.
        urls = [(f"M{i}", f"{base}/gone{i}.gz") for i in range(dl.FAILS_IN_A_ROW + 5)] + [("OK", f"{base}/ok.gz")]
        delivered, failed = dl.fetch_direct(opts, urls, folder)
        self.assertEqual(delivered, {"OK"})
        self.assertEqual(set(failed.values()), {"not found"})

    def test_strain_choice(self):
        # Strains are taken by quality: isolates before MAGs, complete assemblies before drafts, long-read
        # assemblies before the others; ties are drawn at random; only the best candidates by what the metadata
        # says are asked for their sequencing technology.
        rows = []

        def genome(acc, species, category="isolate", level="contig", contigs=50, rep=False, completeness=99.0):
            lineage = f"d__Bacteria;p__P;c__C;o__O;f__F;g__G;s__{species}"
            rows.append(dl.Genome(acc, acc.replace("GCF_", "GCA_"), "s__" + species, lineage, rep, completeness, 1.0,
                                  category, level, contigs))

        genome("GCF_1.1", "A", rep=True)
        genome("GCA_2.1", "A", "mag", "complete genome", 1)
        genome("GCA_3.1", "A", "isolate", "contig", 40)
        genome("GCA_4.1", "A", "isolate", "complete genome", 2)
        genome("GCA_5.1", "A", "isolate", "complete genome", 2)
        genome("GCA_6.1", "A", "isolate", "scaffold", 8)
        genome("GCA_7.1", "A", "isolate", "complete genome", 1, completeness=70.0)  # fails the filter
        genome("GCF_8.1", "B", rep=True)
        genome("GCA_9.1", "B", "mag", "contig", 300)
        genome("GCA_10.1", "B", "mag", "contig", 30)
        genome("GCA_11.1", "B", "sag", "contig", 400)
        opts = argparse.Namespace(seed=1, species=2, per_species=2, rep_only_species=0, min_completeness=90.0,
                                  max_contamination=5.0, tech_candidates=3)
        asked = []

        def tech(accessions):
            asked.append(accessions)
            return {a: "Oxford Nanopore MinION" if a == "GCA_5.1" else "Illumina" for a in accessions}

        picked, chosen, _, _, representative, candidates = dl.pick(rows, opts, tech)
        by_species = {s: [g.accession for g, _ in picked if g.species == s] for s in chosen}
        # A: the shortlist is the 3 best isolates (two complete, one scaffold), the Nanopore one first
        self.assertEqual(by_species["s__A"], ["GCA_5.1", "GCA_4.1"])
        # B: a SAG before MAGs, then the MAG with fewer contigs
        self.assertEqual(by_species["s__B"], ["GCA_11.1", "GCA_10.1"])
        self.assertEqual(sorted(asked[0]), ["GCA_10.1", "GCA_11.1", "GCA_4.1", "GCA_5.1", "GCA_6.1", "GCA_9.1"])
        self.assertEqual(len(asked), 1, "one lookup for all species")
        self.assertEqual(dict((g.accession, t) for g, t in picked)["GCA_5.1"], "Oxford Nanopore MinION")
        self.assertEqual(sorted(representative), ["s__A", "s__B"])
        self.assertNotIn("GCA_7.1", [g.accession for g in candidates["s__A"]])
        # Without a lookup, the same quality order (level, then contigs) decides
        picked, *_ = dl.pick(rows, opts, None)
        self.assertEqual(sorted(g.accession for g, _ in picked if g.species == "s__A"), ["GCA_4.1", "GCA_5.1"])
        # Genomes that tie are drawn at random, not by name: another seed picks others
        for i in range(12):
            genome(f"GCA_{100 + i}.1", "C", "isolate", "complete genome", 1)
        genome("GCF_99.1", "C", rep=True)
        opts.species = 3
        drawn = set()
        for seed in range(1, 9):
            opts.seed = seed
            picked, *_ = dl.pick(rows, opts, None)
            drawn.add(tuple(g.accession for g, _ in picked if g.species == "s__C"))
        self.assertGreater(len(drawn), 3)
        # Only the genomes of a proGenomes table are taken as strains
        opts.seed = 1
        picked, *_ = dl.pick(rows, opts, None, allowed={"GCA_3.1", "GCA_9.1", "GCA_5.1"})
        self.assertEqual(sorted(g.accession for g, _ in picked), ["GCA_3.1", "GCA_5.1", "GCA_9.1"])

    def test_species_choice_kept_on_rerun(self):
        # A rerun with more species with strains keeps the earlier choice (its strain species and their strains) and
        # takes the new strain species first among the earlier representative-only ones, whose genomes are there.
        rows = []
        for i in range(40):
            name = f"s__S{i:02d}"
            lineage = f"d__Bacteria;p__P;c__C;o__O;f__F;g__G;{name}"
            rows.append(dl.Genome(f"GCF_{i}.1", f"GCA_{i}.1", name, lineage, True, 99.0, 1.0, "isolate", "complete genome", 1))
            if i < 30:  # 30 species with two candidate strains each, the second ranked lower (more contigs)
                for j in range(2):
                    rows.append(dl.Genome(f"GCA_{i}0{j}.1", f"GCA_{i}0{j}.1", name, lineage, False, 99.0, 1.0, "isolate",
                                          "contig", 10 + j))
        opts = argparse.Namespace(seed=1, species=10, per_species=1, rep_only_species=10, min_completeness=90.0,
                                  max_contamination=5.0, tech_candidates=3)
        picked, chosen, rep_only, *_ = dl.pick(rows, opts)
        self.assertEqual((len(chosen), len(rep_only)), (10, 10))
        self.assertEqual(dl.pick(rows, opts)[1], chosen)  # the same seed, the same choice
        # The earlier run took each species' second-ranked strain (a tie broken otherwise, say): it stays.
        previous = {"strains": {s: [f"GCA_{int(s[4:])}01.1"] for s in chosen}, "rep_only": set(rep_only)}
        opts.species = 20
        picked2, chosen2, rep_only2, *_ = dl.pick(rows, opts, None, None, previous)
        self.assertEqual(len(chosen2), 20)
        self.assertTrue(set(chosen) <= set(chosen2))
        self.assertEqual({g.species: g.accession for g, _ in picked2 if g.species in chosen},
                         {s: previous["strains"][s][0] for s in chosen})
        eligible_before = {s for s in rep_only if int(s[4:]) < 30}  # earlier representative-only species with strains
        new = set(chosen2) - set(chosen)
        self.assertTrue(eligible_before <= new if len(eligible_before) <= 10 else new <= eligible_before)
        # Its other representative-only species stay so; the pool holds them all still.
        self.assertTrue((set(rep_only) - set(chosen2)) <= set(rep_only2))
        # A new folder takes the seed's order: the first of it with strains.
        self.assertEqual(dl.previous_choice(os.path.join(self.tmp.name, "no_such_folder")), {"strains": {}, "rep_only": set()})

    def test_read_metadata_quality_columns(self):
        # GTDB's columns for the quality of a genome (names of metadata_field_desc.tsv), with none, missing
        # and odd values.
        folder = os.path.join(self.tmp.name, "meta_quality")
        os.makedirs(folder)
        header = ["accession", "checkm2_completeness", "checkm2_contamination", "gtdb_representative", "gtdb_taxonomy",
                  "ncbi_assembly_level", "ncbi_genome_category", "ncbi_contig_count", "contig_count",
                  "ncbi_genbank_assembly_accession", "ncbi_organism_name"]
        lineage = "d__Bacteria;p__P;c__C;o__O;f__F;g__G;s__G a"
        rows = [["RS_GCF_000000001.1", "99.5", "0.5", "t", lineage, "Complete Genome", "none", "2", "2",
                 "GCA_000000001.1", "Gus a"],
                ["GB_GCA_000000002.1", "95", "1", "f", lineage, "Contig", "derived from metagenome", "none", "310",
                 "GCA_000000002.1", "uncultured Gus"],
                ["GB_GCA_000000003.1", "91", "none", "f", lineage, "", "derived from single cell", "", "",
                 "none", ""]]
        with open(os.path.join(folder, "bac120_metadata_r226.tsv"), "w") as fh:
            fh.write("\t".join(header) + "\n" + "".join("\t".join(r) + "\n" for r in rows))
        a, b, c = dl.read_metadata(folder, "226")
        self.assertEqual((a.accession, a.genbank, a.is_rep, a.category, a.level, a.contigs),
                         ("GCF_000000001.1", "GCA_000000001.1", True, "isolate", "complete genome", 2))
        self.assertEqual((b.accession, b.category, b.level, b.contigs), ("GCA_000000002.1", "mag", "contig", 310))
        self.assertEqual((c.category, c.level, c.contigs, c.genbank, c.contamination),
                         ("sag", "", None, "GCA_000000003.1", 0.0))
        self.assertLess(dl.quality_rank(a), dl.quality_rank(b))
        self.assertLess(dl.quality_rank(b, True), dl.quality_rank(b))

    def test_ncbi_summary_and_categories(self):
        self.assertEqual([dl.genome_category(c, o) for c, o in (
            ("none", "Escherichia coli"), ("", ""), ("derived from metagenome", ""), ("derived from single cell", "x"),
            ("none", "uncultured Pseudoalteromonas sp."), ("none", "bacterium metagenome bin 3"),
            ("derived from environmental sample", ""))], ["isolate", "isolate", "mag", "sag", "mag", "mag", "mag"])
        for text in ("PacBio Sequel", "Oxford Nanopore MinION; Illumina", "Illumina HiSeq; PacBio RS II", "PacBio HiFi"):
            self.assertTrue(dl.LONG_READS.search(text), text)
        for text in ("Illumina HiSeq 2500", "454 GS FLX Titanium", "Sanger dideoxy sequencing", "Ion Torrent"):
            self.assertFalse(dl.LONG_READS.search(text), text)
        datasets = fake_datasets(os.path.join(self.tmp.name, "datasets_summary"))
        opts = argparse.Namespace(datasets=datasets)
        work = os.path.join(self.tmp.name, "summary_work")
        for snake in ("", "1"):  # datasets writes assemblyInfo or assembly_info, by version
            os.environ["FAKE_LONG"], os.environ["FAKE_SNAKE"] = "GCA_1.1", snake
            try:
                found, why = dl.ncbi_summary(opts, ["GCA_1.1", "GCA_2.1"], work)
            finally:
                del os.environ["FAKE_LONG"], os.environ["FAKE_SNAKE"]
            self.assertEqual((found, why), ({"GCA_1.1": "PacBio Sequel; Illumina HiSeq", "GCA_2.1": "Illumina HiSeq"}, ""))
        os.environ["FAKE_DOWN"] = "1"
        try:
            found, why = dl.ncbi_summary(opts, ["GCA_1.1"], work)
        finally:
            del os.environ["FAKE_DOWN"]
        self.assertEqual((found, why), ({}, "datasets summary failed (1): Error: Gateway Timeout"))

    def test_download_goes_on_where_a_connection_dropped(self):
        data = random.Random(5).randbytes(50000)
        md5 = hashlib.md5(data).hexdigest()
        requests = []

        class Server(http.server.BaseHTTPRequestHandler):
            drop = None  # bytes of each response sent before the connection closes; None: all

            def log_message(self, *args):
                pass

            def do_GET(self):  # Range requests as GTDB's server answers them (206, 416 from the end on)
                requests.append(self.headers.get("Range"))
                start = int(self.headers["Range"].split("=")[1].split("-")[0]) if self.headers.get("Range") else 0
                if start >= len(data):
                    self.send_response(416)
                    self.send_header("Content-Range", f"bytes */{len(data)}")
                    self.send_header("Content-Length", "0")
                    self.end_headers()
                    return
                self.send_response(206 if start else 200)
                if start:
                    self.send_header("Content-Range", f"bytes {start}-{len(data) - 1}/{len(data)}")
                self.send_header("Content-Length", str(len(data) - start))
                self.end_headers()
                self.wfile.write(data[start:] if Server.drop is None else data[start:start + Server.drop])
                self.close_connection = True

        url = self.serve_handler(Server) + "/big.tar.gz"
        dest = os.path.join(self.tmp.name, "big.tar.gz")
        quiet = contextlib.redirect_stdout(io.StringIO())
        # The connection drops after 12,000 bytes each time, without an error: each request goes on from there.
        Server.drop = 12000
        with quiet:
            dl.download(url, dest, md5)
        with open(dest, "rb") as fh:
            self.assertEqual(fh.read(), data)
        self.assertEqual(requests, [None, "bytes=12000-", "bytes=24000-", "bytes=36000-", "bytes=48000-"])
        # A whole .part, left by a run stopped before it was checked: the server has nothing after its end (416).
        os.replace(dest, dest + ".part")
        requests.clear()
        Server.drop = None
        with quiet:
            dl.download(url, dest, md5)
        self.assertEqual(requests, ["bytes=50000-"])
        self.assertTrue(os.path.isfile(dest) and not os.path.exists(dest + ".part"))
        # Requests that bring nothing new, three in a row: give up, keeping what came for a rerun.
        with open(dest + ".part", "wb") as fh:
            fh.write(data[:1000])
        Server.drop = 0
        with quiet, self.assertRaises(SystemExit) as stopped:
            dl.download(url, os.path.join(self.tmp.name, "big.tar.gz"), md5)
        self.assertIn("a rerun goes on from", str(stopped.exception))
        self.assertEqual(os.path.getsize(dest + ".part"), 1000)


if __name__ == "__main__":
    unittest.main()
