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

import collections
import contextlib
import csv
import functools
import gzip
import hashlib
import http.server
import io
import json
import os
import random
import shutil
import signal
import subprocess
import sys
import tarfile
import tempfile
import threading
import time
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SIMULATE = os.path.join(HERE, "simulate_gtdb_release.py")
CONVERT = os.path.join(HERE, "gtdb_to_protal_db.py")
SIMULATE_READS = os.path.join(HERE, "simulate_reads.py")
LINEAGES = os.path.join(HERE, "gtdb_like_lineages.py")
DOWNLOAD = os.path.join(HERE, "..", "download_gtdb.py")

# A stand-in for NCBI's `datasets`: `download genome accession` writes the list into the zip,
# `rehydrate` copies the synthetic release's genomes; accessions in $FAKE_SUPPRESSED fail a request, and
# every request fails with $FAKE_DOWN set (NCBI down). `summary genome accession` reports the sequencing
# technology: PacBio for the accessions in $FAKE_LONG, Illumina for the others (keys in snake_case with
# $FAKE_SNAKE). $FAKE_LOG, if set, gets a line per call.
FAKE_DATASETS = r'''#!/usr/bin/env python3
import csv, gzip, json, os, shutil, sys, zipfile
args = sys.argv[1:]
if os.environ.get("FAKE_LOG"):
    open(os.environ["FAKE_LOG"], "a").write(" ".join(args[:3]) + "\n")
if os.environ.get("FAKE_DOWN"):
    sys.exit("Error: Gateway Timeout")
if args[:3] == ["summary", "genome", "accession"]:
    snake = bool(os.environ.get("FAKE_SNAKE"))
    for a in open(args[args.index("--inputfile") + 1]).read().split():
        tech = "PacBio Sequel; Illumina HiSeq" if a in os.environ.get("FAKE_LONG", "").split(",") else "Illumina HiSeq"
        info = {"assembly_level": "Contig", "sequencing_tech": tech} if snake else {"assemblyLevel": "Contig", "sequencingTech": tech}
        print(json.dumps({"accession": a, "assembly_info" if snake else "assemblyInfo": info}))
elif args[:3] == ["download", "genome", "accession"]:
    accessions = open(args[args.index("--inputfile") + 1]).read().split()
    if set(accessions) & set(os.environ.get("FAKE_SUPPRESSED", "").split(",")):
        sys.exit("Error: some accessions are not valid")
    with zipfile.ZipFile(args[args.index("--filename") + 1], "w") as z:
        z.writestr("ncbi_dataset/fetch.txt", "\n".join(accessions))
elif args[0] == "rehydrate":
    folder = args[args.index("--directory") + 1]
    paths = {r["accession"]: r["fasta_path"] for r in csv.DictReader(open(os.environ["FAKE_TABLE"]), delimiter="\t")}
    for a in open(os.path.join(folder, "ncbi_dataset", "fetch.txt")).read().split():
        os.makedirs(os.path.join(folder, "ncbi_dataset", "data", a), exist_ok=True)
        with gzip.open(paths[a]) as fin, open(os.path.join(folder, "ncbi_dataset", "data", a, a + "_ASM1v1_genomic.fna"), "wb") as fout:
            shutil.copyfileobj(fin, fout)
'''


def assembly_folder(accession):
    """The name of a synthetic genome's folder on the stand-in FTP server (ASM<digits>v1 is its assembly name)."""
    return f"{accession}_ASM{accession[4:13]}v1"


def gtdb_mirror(release, root, assembly_names=False):
    """release, a synthetic GTDB release, laid out as GTDB's server has it under root/release226/226.0. With
    assembly_names, its metadata has the ncbi_assembly_name column (the synthetic one has none)."""
    base = os.path.join(root, "release226", "226.0")
    os.makedirs(os.path.join(base, "genomic_files_reps"))
    os.makedirs(os.path.join(base, "genomic_files_all"))
    for mset in ("bac120", "ar53"):
        with open(os.path.join(release, f"{mset}_taxonomy_r226.tsv"), "rb") as fin, \
                gzip.open(os.path.join(base, f"{mset}_taxonomy_r226.tsv.gz"), "wb") as fout:
            fout.write(fin.read())
        with open(os.path.join(release, f"{mset}_metadata_r226.tsv.gz"), "rb") as fin, \
                open(os.path.join(base, f"{mset}_metadata_r226.tsv.gz"), "wb") as fout:
            fout.write(fin.read())
        if assembly_names:
            path = os.path.join(base, f"{mset}_metadata_r226.tsv.gz")
            with gzip.open(path, "rt") as fh:
                lines = fh.read().splitlines()
            column = lines[0].split("\t").index("accession")
            out = [lines[0] + "\tncbi_assembly_name"]
            for line in lines[1:]:
                accession = line.split("\t")[column]
                accession = accession[3:] if accession[:3] in ("RS_", "GB_") else accession
                out.append(f"{line}\tASM{accession[4:13]}v1")
            with gzip.open(path, "wt", newline="\n") as fh:
                fh.write("\n".join(out) + "\n")
        for kind in ("reps", "all"):
            name = f"{mset}_marker_genes_{kind}_r226"
            with tarfile.open(os.path.join(base, f"genomic_files_{kind}", name + ".tar.gz"), "w:gz") as tar:
                tar.add(os.path.join(release, f"genomic_files_{kind}", name), arcname=name)
    with open(os.path.join(base, "VERSION.txt"), "w") as fh:
        fh.write("v226.0\n")
    with open(os.path.join(base, "MD5SUM.txt"), "w") as fh:
        for folder, _dirs, names in os.walk(base):
            for name in sorted(names):
                if name != "MD5SUM.txt":
                    path = os.path.join(folder, name)
                    with open(path, "rb") as f:
                        fh.write(f"{hashlib.md5(f.read()).hexdigest()}  ./{os.path.relpath(path, base)}\n")


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

    def test_parallel_conversion(self):
        # The marker files are read by -t workers; the database must not depend on how many.
        for order in ("gene", "genome"):
            one, four = (os.path.join(self.tmp.name, f"db_{order}_t{t}") for t in (1, 4))
            run(CONVERT, "--gtdb", self.gtdb, "--outdir", one, "--order", order, "-t", "1")
            run(CONVERT, "--gtdb", self.gtdb, "--outdir", four, "--order", order, "-t", "4")
            for f in ("reference.fna", "reference.map", "full_reference.fna", "internal_taxonomy.dmp"):
                with open(os.path.join(one, f), "rb") as a, open(os.path.join(four, f), "rb") as b:
                    self.assertEqual(a.read(), b.read(), f"{order} order, {f}")
            self.assertFalse(os.path.exists(os.path.join(four, ".convert_tmp")))

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

    def test_a_copy_takes_only_the_converted_files(self):
        # What a build of the folder wrote, or left when stopped, is not copied (its unique_kmers.tsv, of other
        # genes, stopped the copy's build), and what an earlier build of the copy left is removed.
        src, dst = os.path.join(self.tmp.name, "db_built"), os.path.join(self.tmp.name, "db_training")
        shutil.copytree(self.db, src)
        for name in ("unique_kmers.tsv", "index.prx.zst.partial", "database.protal.partial", "database.protal",
                     "build_metadata.tsv"):
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
        sys.path.insert(0, HERE)
        from gtdb_to_protal_db import CONVERTED_FILES
        expected = {"reference.fna", "reference.map", "full_reference.fna"} | \
            {f for f in CONVERTED_FILES if os.path.isfile(os.path.join(self.db, f))}
        self.assertEqual(set(os.listdir(dst)), expected)
        # Converting a release anew removes the build outputs of the folder's earlier reference, too.
        run(CONVERT, "--gtdb", self.gtdb, "--outdir", src)
        for name in ("unique_kmers.tsv", "index.prx.zst.partial", "database.protal.partial", "database.protal"):
            self.assertFalse(os.path.exists(os.path.join(src, name)), name)

    def test_download_gtdb(self):
        mirror = os.path.join(self.tmp.name, "mirror")
        gtdb_mirror(self.gtdb, mirror)
        class Quiet(http.server.SimpleHTTPRequestHandler):
            def log_message(self, *args):
                pass

        server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), functools.partial(Quiet, directory=mirror))
        threading.Thread(target=server.serve_forever, daemon=True).start()
        self.addCleanup(server.server_close)
        self.addCleanup(server.shutdown)
        datasets = os.path.join(self.tmp.name, "datasets")
        with open(datasets, "w") as fh:
            fh.write(FAKE_DATASETS)
        os.chmod(datasets, 0o755)
        # GCA_999002003.1 is the PacBio assembly among the strains of Mockella beta, so it is the one picked,
        # and NCBI no longer has it.
        env = dict(os.environ, FAKE_TABLE=os.path.join(self.gtdb, "simulation", "genomes.tsv"),
                   FAKE_SUPPRESSED="GCA_999002003.1", FAKE_LONG="GCA_999002003.1")
        out = os.path.join(self.tmp.name, "inputs")
        command = [sys.executable, DOWNLOAD, "-o", out, "--mirror", f"http://127.0.0.1:{server.server_port}",
                   "--datasets", datasets, "--species", "3", "--per_species", "1", "--rep_only_species", "0",
                   "--batch", "4", "-t", "2"]
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
        run(CONVERT, "--gtdb", release, "--outdir", os.path.join(self.tmp.name, "db_downloaded"))
        with open(os.path.join(self.db, "reference.fna"), "rb") as a, \
                open(os.path.join(self.tmp.name, "db_downloaded", "reference.fna"), "rb") as b:
            self.assertEqual(a.read(), b.read())

        missing = subprocess.run([sys.executable, DOWNLOAD, "-o", out + "_x", "--release", "999", "--mirror",
                                  f"http://127.0.0.1:{server.server_port}", "--no_genomes"], capture_output=True, text=True)
        self.assertNotEqual(missing.returncode, 0)

        # With NCBI down, a few failed requests in a row stop the run (not ~2n of them halving every batch),
        # in the technology lookup or, without it, in the download, and it does not say its inputs are ready.
        for extra, failed in (([], "datasets summary failed (1)"), (["--no_tech_lookup"], "datasets download failed (1)")):
            calls = os.path.join(self.tmp.name, "datasets_calls.txt")
            if os.path.exists(calls):
                os.remove(calls)
            down = subprocess.run([sys.executable, DOWNLOAD, "-o", out + "_down", "--mirror",
                                   f"http://127.0.0.1:{server.server_port}", "--datasets", datasets, "--species", "3",
                                   "--per_species", "1", "--rep_only_species", "0", "--batch", "1", "-t", "2", *extra],
                                  env=dict(env, FAKE_DOWN="1", FAKE_LOG=calls), capture_output=True, text=True)
            self.assertNotEqual(down.returncode, 0, down.stdout)
            self.assertIn(f"NCBI requests failed in a row, the last: {failed}: Error: Gateway Timeout", down.stderr)
            self.assertNotIn("Inputs for GTDB", down.stdout)
            with open(calls) as fh:
                self.assertEqual(len(fh.readlines()), 5)  # --batch 1: 1 .bit_length() + 4

    def test_download_direct(self):
        # Genomes come straight from the (stand-in) FTP server, in parallel, without datasets; the ones it does
        # not have go through datasets.
        mirror = os.path.join(self.tmp.name, "mirror_direct")
        gtdb_mirror(self.gtdb, mirror, assembly_names=True)
        with open(os.path.join(self.gtdb, "simulation", "genomes.tsv")) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                a = row["accession"]
                folder = os.path.join(mirror, "ftp", a[:3], a[4:7], a[7:10], a[10:13], assembly_folder(a))
                os.makedirs(folder)
                shutil.copy(row["fasta_path"], os.path.join(folder, assembly_folder(a) + "_genomic.fna.gz"))

        class Quiet(http.server.SimpleHTTPRequestHandler):
            def log_message(self, *args):
                pass

        server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), functools.partial(Quiet, directory=mirror))
        threading.Thread(target=server.serve_forever, daemon=True).start()
        self.addCleanup(server.server_close)
        self.addCleanup(server.shutdown)
        datasets = os.path.join(self.tmp.name, "datasets_direct")
        with open(datasets, "w") as fh:
            fh.write(FAKE_DATASETS)
        os.chmod(datasets, 0o755)
        port = server.server_port
        env = dict(os.environ, FAKE_TABLE=os.path.join(self.gtdb, "simulation", "genomes.tsv"))

        def download(out, *extra, **more):
            return subprocess.run([sys.executable, DOWNLOAD, "-o", out, "--mirror", f"http://127.0.0.1:{port}",
                                   "--ftp_url", f"http://127.0.0.1:{port}/ftp", "--species", "3", "--per_species", "1",
                                   "--rep_only_species", "0", "--connections", "4", *extra],
                                  env=dict(env, **more), capture_output=True, text=True)

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
            self.assertEqual(json.load(fh)["genomes"]["fetched_directly"], 6)
        again = download(out, "--datasets", "/nonexistent/datasets", "--no_tech_lookup")
        self.assertEqual(again.returncode, 0, again.stdout + again.stderr)
        self.assertIn("6 already there, 0 to download", again.stdout)

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

    def test_fetch_file(self):
        sys.path.insert(0, os.path.join(HERE, ".."))
        import argparse
        import download_gtdb as dl
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

        server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), Handler)
        threading.Thread(target=server.serve_forever, daemon=True).start()
        self.addCleanup(server.server_close)
        self.addCleanup(server.shutdown)
        base = f"http://127.0.0.1:{server.server_port}"
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
        sys.path.insert(0, os.path.join(HERE, ".."))
        import argparse
        import download_gtdb as dl
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

    def test_read_metadata_quality_columns(self):
        # GTDB's columns for the quality of a genome (names of metadata_field_desc.tsv), with none, missing
        # and odd values.
        sys.path.insert(0, os.path.join(HERE, ".."))
        import download_gtdb as dl
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
        sys.path.insert(0, os.path.join(HERE, ".."))
        import argparse
        import download_gtdb as dl
        self.assertEqual([dl.genome_category(c, o) for c, o in (
            ("none", "Escherichia coli"), ("", ""), ("derived from metagenome", ""), ("derived from single cell", "x"),
            ("none", "uncultured Pseudoalteromonas sp."), ("none", "bacterium metagenome bin 3"),
            ("derived from environmental sample", ""))], ["isolate", "isolate", "mag", "sag", "mag", "mag", "mag"])
        for text in ("PacBio Sequel", "Oxford Nanopore MinION; Illumina", "Illumina HiSeq; PacBio RS II", "PacBio HiFi"):
            self.assertTrue(dl.LONG_READS.search(text), text)
        for text in ("Illumina HiSeq 2500", "454 GS FLX Titanium", "Sanger dideoxy sequencing", "Ion Torrent"):
            self.assertFalse(dl.LONG_READS.search(text), text)
        datasets = os.path.join(self.tmp.name, "datasets_summary")
        with open(datasets, "w") as fh:
            fh.write(FAKE_DATASETS)
        os.chmod(datasets, 0o755)
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
        sys.path.insert(0, os.path.join(HERE, ".."))
        import download_gtdb
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

        server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), Server)
        threading.Thread(target=server.serve_forever, daemon=True).start()
        self.addCleanup(server.server_close)
        self.addCleanup(server.shutdown)
        url, dest = f"http://127.0.0.1:{server.server_port}/big.tar.gz", os.path.join(self.tmp.name, "big.tar.gz")
        quiet = contextlib.redirect_stdout(io.StringIO())
        # The connection drops after 12,000 bytes each time, without an error: each request goes on from there.
        Server.drop = 12000
        with quiet:
            download_gtdb.download(url, dest, md5)
        with open(dest, "rb") as fh:
            self.assertEqual(fh.read(), data)
        self.assertEqual(requests, [None, "bytes=12000-", "bytes=24000-", "bytes=36000-", "bytes=48000-"])
        # A whole .part, left by a run stopped before it was checked: the server has nothing after its end (416).
        os.replace(dest, dest + ".part")
        requests.clear()
        Server.drop = None
        with quiet:
            download_gtdb.download(url, dest, md5)
        self.assertEqual(requests, ["bytes=50000-"])
        self.assertTrue(os.path.isfile(dest) and not os.path.exists(dest + ".part"))
        # Requests that bring nothing new, three in a row: give up, keeping what came for a rerun.
        with open(dest + ".part", "wb") as fh:
            fh.write(data[:1000])
        Server.drop = 0
        with quiet, self.assertRaises(SystemExit) as stopped:
            download_gtdb.download(url, os.path.join(self.tmp.name, "big.tar.gz"), md5)
        self.assertIn("a rerun goes on from", str(stopped.exception))
        self.assertEqual(os.path.getsize(dest + ".part"), 1000)

    def test_lineages_and_clade_holdout(self):
        # The training scripts read lineages from internal_taxonomy.dmp; build_gtdb_database.py holds out whole
        # clades and then single species.
        sys.path.insert(0, os.path.join(HERE, ".."))
        import build_gtdb_database as build
        import collect_training_data as collect
        import lineages
        taxonomy = os.path.join(self.db, "internal_taxonomy.dmp")
        by_id, ids = lineages.from_taxonomy(taxonomy)
        with open(os.path.join(self.db, "genome2tiid.tsv")) as fh:
            for accession, taxid, _rep, lineage in (line.rstrip("\n").split("\t") for line in fh):
                self.assertEqual(by_id[taxid], lineages.from_string(lineage), accession)
        table = os.path.join(self.gtdb, "simulation", "genomes.tsv")
        pool = build.pool_species(table)
        species = {lin["species"]: lin for lin in by_id.values() if "species" in lin}
        families = {}
        for name, lin in species.items():
            families.setdefault(lin["family"], set()).add(name)
        eligible = {f for f, members in families.items() if len(members & pool) >= 2}
        self.assertTrue(eligible, "the synthetic release needs a family with two species")
        chosen = build.choose_holdout(table, taxonomy, 0.5, {"family": 1}, 1.0, 3)
        self.assertEqual(chosen, build.choose_holdout(table, taxonomy, 0.5, {"family": 1}, 1.0, 3))
        held_family = {clade for rank, clade in chosen.values() if rank == "family"}
        self.assertEqual(len(held_family), 1)
        family = held_family.pop()
        self.assertIn(family, eligible)
        self.assertEqual({s for s, (rank, _) in chosen.items() if rank == "family"}, families[family])
        rest = sorted(pool - families[family])
        self.assertEqual(sum(rank == "species" for rank, _ in chosen.values()),
                         sum(round(0.5 * len([s for s in rest if species[s]["domain"] == d]))
                             for d in {species[s]["domain"] for s in rest}))
        # a clade holding more than the share, or inside one taken before, is not drawn
        self.assertEqual(build.choose_holdout(table, taxonomy, 0, {"family": 1}, 0.0, 3), {})
        phylum_first = build.choose_holdout(table, taxonomy, 0, {"phylum": 1, "family": 5}, 1.0, 3)
        phylum = next(c for r, c in phylum_first.values() if r == "phylum")
        for s, (rank, clade) in phylum_first.items():
            if rank == "family":
                self.assertNotEqual(species[s]["phylum"], phylum)
        # each rank alone takes a whole clade of that rank; the synthetic release has one clade with two
        # species, from phylum to genus, so with every rank asked for, the first (phylum) takes it
        for rank in build.CLADE_RANKS:
            alone = build.choose_holdout(table, taxonomy, 0, {rank: 1}, 1.0, 3)
            taken = {c for r, c in alone.values()}
            self.assertEqual({r for r, _ in alone.values()}, {rank})
            self.assertEqual(len(taken), 1)
            clade = taken.pop()
            self.assertEqual(set(alone), {s for s, lin in species.items() if lin.get(rank) == clade})
            self.assertEqual(collect.novel_clades(alone, table, 1), {rank: [clade]})
        every = build.choose_holdout(table, taxonomy, 0, build.parse_clades("phylum:1,class:1,order:1,family:1,genus:2"), 1.0, 3)
        self.assertEqual({r for r, _ in every.values()}, {"phylum"})
        # heldout_species.txt round trip, and the collector's reading of it
        path = os.path.join(self.tmp.name, "heldout.txt")
        with open(path, "w") as fh:
            fh.write("".join(f"{s}\t{r}\t{c}\n" for s, (r, c) in sorted(chosen.items())) + "s__Alone one\n")
        self.assertEqual(build.read_holdout(path), {**chosen, "s__Alone one": ("species", "s__Alone one")})
        self.assertEqual(collect.read_novel(path), build.read_holdout(path))
        self.assertEqual(build.parse_clades("phylum:2,class:4"), {"phylum": 2, "class": 4})
        self.assertEqual(build.parse_clades("none"), {})

    def test_collector_designs(self):
        # Read setups (built-in and custom ART profiles), abundance models, long-read setups and units.
        sys.path.insert(0, os.path.join(HERE, ".."))
        import argparse
        import collect_training_data as collect
        opts = argparse.Namespace(read_setups="150:HSXt:350:50,150:file=/p1+/p2:350:50,250:MSv3:550:50",
                                  read_pairs="1000,5000", read_types=["pe", "se", "ont"],
                                  pb_setup="errhmm:ERRHMM-SEQUEL:15000:3000:0.999",
                                  ont_setup="qshmm:QSHMM-ONT-HQ:8000:6000:0.97", long_read_bases="1e6,2e6,3e6")
        points = collect.design_points(opts)
        self.assertEqual([p["name"] for p in points], ["rl150_HSXt_p1000", "rl150_HSXt_p5000", "rl150_custom1_p1000",
                                                       "rl150_custom1_p5000", "rl250_p1000", "rl250_p5000"])
        # Setups of one length and profile are told apart by their fragments; a setup or depth given twice
        # would share a folder.
        same = argparse.Namespace(read_setups="150:HS25:350:50,150:HS25:500:80", read_pairs="1000")
        self.assertEqual([p["name"] for p in collect.design_points(same)], ["rl150_HS25_f350-50_p1000", "rl150_HS25_f500-80_p1000"])
        for setups, pairs in (("150:HS25:350:50,150:HS25:350:50", "1000"), ("150:HS25:350:50", "1000,1000")):
            with self.assertRaises(SystemExit):
                collect.design_points(argparse.Namespace(read_setups=setups, read_pairs=pairs))
        with self.assertRaises(SystemExit):
            collect.units_of(argparse.Namespace(**{**vars(opts), "long_read_bases": "1e6,1e6"}))
        with self.assertRaises(SystemExit):
            collect.art_profile_args("file=/nonexistent_r1.txt")
        self.assertEqual(collect.art_profile_args("HSXt"), ["--sequencer", "HSXt"])
        self.assertEqual(collect.abundance_args("lognormal:2.0"), ["--distribution", "poisson_lognormal", "--pln_sigma", "2.0"])
        self.assertEqual(collect.abundance_args("powerlaw:1.5"), ["--distribution", "power_law", "--alpha", "1.5"])
        self.assertEqual(collect.abundance_args(""), [])
        with self.assertRaises(SystemExit):
            collect.abundance_args("gamma:1")
        self.assertEqual(collect.parse_long_setup("qshmm:QSHMM-ONT-HQ:8000:6000:0.97")["length_mean"], 8000)
        with self.assertRaises(SystemExit):
            collect.parse_long_setup("art:x:1:1:1")
        _, units = collect.units_of(opts)
        self.assertEqual([u["type"] for u in units], ["pe"] * 6 + ["se"] * 6 + ["ont"] * 3)
        self.assertEqual([u["community"]["name"] for u in units if u["type"] == "ont"],
                         ["rl150_HSXt_p1000", "rl150_HSXt_p5000", "rl150_custom1_p1000"])
        self.assertEqual(units[6]["name"], "rl150_HSXt_p1000_se")

    def test_long_read_replay(self):
        # pb/ont samples replay a paired-end point's communities: each genome gets its share of the bases by
        # relative abundance times length (a stand-in pbsim writes depth x length / 1000 reads of 1 kb).
        sys.path.insert(0, os.path.join(HERE, ".."))
        import argparse
        import collect_training_data as collect
        root = os.path.join(self.tmp.name, "longreads")
        point = os.path.join(root, "points", "rl150_p1000")
        os.makedirs(os.path.join(point, "sim"))
        genomes = {}
        for name, length in (("GA", 20000), ("GB", 40000)):
            genomes[name] = os.path.join(root, name + ".fna.gz")
            with gzip.open(genomes[name], "wt") as fh:
                fh.write(f">{name}_contig\n" + "ACGT" * (length // 4) + "\n")
        with open(os.path.join(point, "sim", "manifest.tsv"), "w") as fh:
            fh.write("sample\tgenome\tspecies\ttaxonomy\tgenome_length\tread_pairs\tvertical_coverage\trelative_abundance"
                     "\tfastq_r1\tfastq_r2\tfasta_path\tart_seed\n")
            for sample, abundances in (("rl150_p1000_s_1", (0.5, 0.5)), ("rl150_p1000_s_2", (0.8, 0.2))):
                for (name, path), a in zip(genomes.items(), abundances):
                    length = 20000 if name == "GA" else 40000
                    fh.write(f"{sample}\t{name}\tS {name}\td__B;s__S {name}\t{length}\t10\t1\t{a}\tr1\tr2\t{path}\t1\n")
        with open(os.path.join(point, "sim", "protal.meta"), "w") as fh:
            fh.write(f"#OUTPUT_DIR\t{point}/protal\n#INPUT_DIR\t{point}/sim/reads\n"
                     "#SAMPLEID\tFIRST\tSECOND\tSAM\tPREFIX\tPROFILE\tPROFILE_TRUTH\n"
                     "rl150_p1000_s_1\ta\tb\tc\td\te\t/truth/1\nrl150_p1000_s_2\ta\tb\tc\td\te\t/truth/2\n")
        fake = os.path.join(root, "pbsim")
        with open(fake, "w") as fh:
            fh.write("#!" + sys.executable + "\nimport sys\na = sys.argv[1:]\nget = lambda k: a[a.index(k) + 1]\n"
                     "seq = ''.join(l.strip() for l in open(get('--genome')) if not l.startswith('>'))\n"
                     "n = round(float(get('--depth')) * len(seq) / 1000)\n"
                     "with open(get('--prefix') + '_0001.fastq', 'w') as out:\n"
                     "    for i in range(n):\n"
                     "        out.write('@%s1_%d\\n%s\\n+\\n%s\\n' % (get('--id-prefix'), i, seq[:1000], 'I' * 1000))\n")
        os.chmod(fake, 0o755)
        models = os.path.join(root, "models")
        os.makedirs(models)
        open(os.path.join(models, "FAKE.model"), "w").close()
        opts = argparse.Namespace(out=root, seed=1, pbsim=fake, pbsim_models=models, samples=2)
        unit = {"type": "ont", "name": "ont_b60000", "bases": 60000,
                "setup": collect.parse_long_setup("qshmm:FAKE:1000:0:0.97"), "community": {"name": "rl150_p1000"},
                "point": {"name": "ont_b60000", "read_length": "1000", "read_pairs": "60000"}}
        self.assertIsNone(collect.simulate_long(unit, 0, opts, 2))
        rows, _ = collect.unit_map_rows(unit, opts)
        self.assertEqual([r["SAMPLEID"] for r in rows], ["ont_b60000_s_1", "ont_b60000_s_2"])
        self.assertEqual([r["PROFILE_TRUTH"] for r in rows], ["/truth/1", "/truth/2"])
        self.assertEqual({r["READ_TYPE"] for r in rows} | {r["SECOND"] for r in rows}, {"ont", "-"})
        for row, (a, b) in zip(rows, ((0.5, 0.5), (0.8, 0.2))):
            names = [line.strip() for i, line in enumerate(gzip.open(row["FIRST"], "rt")) if i % 4 == 0]
            self.assertEqual(len(names), len(set(names)), "read names must be unique within a sample")
            counts = {g: sum(n.startswith(f"@g{g}x") for n in names) for g in (0, 1)}
            share_a = a * 20000 / (a * 20000 + b * 40000)
            self.assertAlmostEqual(counts[0], round(60 * share_a), delta=1)
            self.assertAlmostEqual(counts[1], round(60 * (1 - share_a)), delta=1)
        self.assertTrue(collect.simulated(unit, opts))

    def test_long_read_short_contigs(self):
        # pbsim3 stops at a reference sequence under 100 bases ("Reference is too short"), as assemblies of
        # metagenomes have: it gets the longer sequences only, and a failure says why.
        sys.path.insert(0, os.path.join(HERE, ".."))
        import collect_training_data as collect
        root = os.path.join(self.tmp.name, "shortcontigs")
        os.makedirs(root)
        fake = os.path.join(root, "pbsim")
        with open(fake, "w") as fh:
            fh.write("#!" + sys.executable + "\nimport sys\na = sys.argv[1:]\nlengths = []\n"
                     "for l in open(a[a.index('--genome') + 1]):\n"
                     "    if l.startswith('>'): lengths.append(0)\n"
                     "    else: lengths[-1] += len(l.strip())\n"
                     "if min(lengths) < 100:\n"
                     "    print('ERROR: Reference is too short. Acceptable length >= 100.')\n"
                     "    sys.exit(255)\n"
                     "open(a[a.index('--prefix') + 1] + '_0001.fastq', 'w').write('@r\\nACGT\\n+\\nIIII\\n')\n")
        os.chmod(fake, 0o755)
        records = [(">a one", "ACGT" * 75), (">short", "ACGT" * 24 + "ACG"), (">exact", "A" * 100), (">tiny", "ACGT")]
        text = "".join(f"{header}\r\n" + "\r\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)) + "\r\n\r\n"
                       for header, seq in records)
        plain, zipped = os.path.join(root, "plain.fna"), os.path.join(root, "zipped.fna.gz")
        with open(plain, "w", newline="") as fh:
            fh.write(text)
        with gzip.open(zipped, "wt", newline="") as fh:
            fh.write(text)

        def task(name, fasta, pbsim=fake):
            return {"genome": name, "fasta": fasta, "tmp": os.path.join(root, "tmp", name), "pbsim": pbsim,
                    "setup": collect.parse_long_setup("qshmm:FAKE:1000:0:0.97"), "model": "FAKE", "depth": 1.0,
                    "seed": 1, "id_prefix": "g0x"}

        for name, fasta in (("plain", plain), ("zipped", zipped)):
            fastqs, error = collect.long_read_genome(task(name, fasta))
            self.assertIsNone(error)
            self.assertEqual(len(fastqs), 1)
            with open(os.path.join(root, "tmp", name, "genome.fna")) as fh:
                kept = [line.strip() for line in fh if line.startswith(">")]
            self.assertEqual(kept, [">a one", ">exact"])
        # What the fix is for: pbsim3 fails on the file as it is.
        self.assertEqual(subprocess.run([fake, "--genome", plain, "--prefix", os.path.join(root, "x")],
                                        capture_output=True).returncode, 255)
        # Only short sequences: nothing to simulate from. A pbsim failure carries its last line of output.
        short = os.path.join(root, "short.fna")
        with open(short, "w") as fh:
            fh.write(">tiny\nACGT\n")
        _, error = collect.long_read_genome(task("short", short))
        self.assertIn("no sequence of 100 bases or more", error)
        failing = os.path.join(root, "failing")
        with open(failing, "w") as fh:
            fh.write("#!/bin/sh\necho 'ERROR: out of ideas'\nexit 3\n")
        os.chmod(failing, 0o755)
        _, error = collect.long_read_genome(task("failing", plain, failing))
        self.assertIn("pbsim failed (3)", error)
        self.assertIn("ERROR: out of ideas", error)

    def test_relation_to_novel_species(self):
        # An absent taxon is put down to a species the database lacks when that species is at least as close to
        # it as every present one.
        sys.path.insert(0, os.path.join(HERE, ".."))
        import collect_training_data as collect
        import lineages
        lin = lambda text: lineages.from_string(text)
        taxon = lin("d__B;p__P;c__C;o__O;f__F;g__G;s__G a")
        in_sample = {"s__G b": lin("d__B;p__P;c__C;o__O;f__F;g__G;s__G b"),
                     "s__H c": lin("d__B;p__P;c__C;o__O2;f__F2;g__H;s__H c")}
        self.assertEqual(collect.relation(taxon, in_sample, {}), ("genus", ""))
        self.assertEqual(collect.relation(taxon, in_sample, {"s__G b": ("species", "s__G b")}), ("genus", "species"))
        self.assertEqual(collect.relation(taxon, in_sample, {"s__H c": ("class", "c__C")}), ("genus", ""))
        self.assertEqual(collect.relation(taxon, {"s__H c": in_sample["s__H c"]}, {"s__H c": ("order", "o__O2")}),
                         ("class", "order"))
        self.assertEqual(collect.relation(taxon, {"s__X": lin("d__A;p__Q;s__X")}, {"s__X": ("phylum", "p__Q")}), ("none", ""))
        # a present taxon's closest other species (meta_neighbour_rank)
        self.assertEqual(collect.relation(taxon, {"s__H c": in_sample["s__H c"]}, {})[0], "class")
        self.assertEqual(collect.relation(taxon, {}, {})[0], "none")

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


BUILD = os.path.join(HERE, "..", "build_gtdb_database.py")


@unittest.skipUnless(os.environ.get("PROTAL") and os.environ.get("SIMULATE") and shutil.which("art_illumina"),
                     "needs $PROTAL, $SIMULATE (simulate_metagenomes) and art_illumina")
class GtdbBuildTest(unittest.TestCase):
    """build_gtdb_database.py end to end, on a synthetic GTDB-like release of 60 species downloaded from a
    fake GTDB mirror and a fake NCBI: the database is built and its pe and se models trained; a rerun skips
    the conversion and both builds; another seed (other species held out) rebuilds the training database
    only, and simulates and profiles again; a build that fails in the background stops the run at once; a
    run stopped with SIGTERM leaves no command running. The trainer needs scikit-learn: $PROTAL_TRAIN_PYTHON
    (default: this Python). About five minutes."""

    @classmethod
    def setUpClass(cls):
        cls.python = os.environ.get("PROTAL_TRAIN_PYTHON", sys.executable)
        if subprocess.run([cls.python, "-c", "import joblib, numpy, pandas, sklearn"], capture_output=True).returncode:
            raise unittest.SkipTest(f"{cls.python} cannot import scikit-learn, joblib, numpy and pandas ($PROTAL_TRAIN_PYTHON)")
        cls.tmp = tempfile.TemporaryDirectory()
        work = cls.tmp.name
        lineages = os.path.join(work, "lineages.txt")
        with open(lineages, "w") as fh:
            subprocess.run([sys.executable, LINEAGES, "--species", "60", "--archaea", "0.1", "--seed", "1"], stdout=fh, check=True)
        gtdb = os.path.join(work, "gtdb")
        run(SIMULATE, "--outdir", gtdb, "--lineages", lineages, "--genomes_per_species", "3", "--genome_length", "40000",
            "--strain_divergence", "0.002-0.015", "--species_divergence", "0.015-0.04", "--seed", "3")
        gtdb_mirror(gtdb, os.path.join(work, "mirror"))

        class Quiet(http.server.SimpleHTTPRequestHandler):
            def log_message(self, *args):
                pass

        server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), functools.partial(Quiet, directory=os.path.join(work, "mirror")))
        threading.Thread(target=server.serve_forever, daemon=True).start()
        datasets = os.path.join(work, "datasets")
        with open(datasets, "w") as fh:
            fh.write(FAKE_DATASETS)
        os.chmod(datasets, 0o755)
        cls.inputs = os.path.join(work, "inputs")
        subprocess.run([sys.executable, DOWNLOAD, "-o", cls.inputs, "--mirror", f"http://127.0.0.1:{server.server_port}",
                        "--datasets", datasets, "--species", "15", "--per_species", "2", "--rep_only_species", "10",
                        "--batch", "8", "-t", "2"], env=dict(os.environ, FAKE_TABLE=os.path.join(gtdb, "simulation", "genomes.tsv")),
                       check=True, capture_output=True)
        server.shutdown()
        server.server_close()

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def build(self, out, *extra, protal=None, wait=True):
        command = [self.python, BUILD, "--inputs", self.inputs, "--outdir", os.path.join(self.tmp.name, out),
                   "--protal", protal or os.environ["PROTAL"], "--simulator", os.environ["SIMULATE"], "-t", "2",
                   "--samples", "2", "--read-pairs", "1000,4000", "--read-setups", "100:HS20:300:40",
                   "--species-per-sample", "6-8", "--archaea", "1", "--holdout-max-share", "0.2",
                   "--holdout-clades", "family:1,genus:1", "--read-types", "pe,se", "--test-samples", "1",
                   "--test-read-pairs", "2000", "--ntree", "16", "--evaluation", "basic", *extra]
        if not wait:
            return subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        return subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=3000)

    def text(self, *path):
        with open(os.path.join(self.tmp.name, *path)) as fh:
            return fh.read()

    def test_a_build_and_rerun(self):
        first = self.build("out")
        self.assertEqual(first.returncode, 0, first.stdout[-3000:])
        self.assertIn("Ready protal database", first.stdout)
        for path in ("protal_db/database.protal", "training_db/database.protal", "model_logs/summary.txt",
                     ".stages/convert.json", ".stages/protal_db.json", ".stages/training_db.json"):
            self.assertTrue(os.path.isfile(os.path.join(self.tmp.name, "out", path)), path)

        # A rerun converts and builds nothing, and the collector reuses its samples and dumps.
        again = self.build("out")
        self.assertEqual(again.returncode, 0, again.stdout[-3000:])
        self.assertIn("protal_db was built by an earlier run from the same release and protal; kept", again.stdout)
        self.assertIn("training_db was built by an earlier run with the same species left out; kept", again.stdout)
        self.assertNotIn("Built ", again.stdout)
        self.assertNotRegex(self.text("out", "training_data.log"), "simulating|profiling")

        # Other species held out (another seed): only the training database is built again, from the release
        # converted anew (the finished database's build consumed the converted files), and every point is
        # simulated and profiled again rather than mixed into the table.
        other = self.build("out", "--seed", "2")
        self.assertEqual(other.returncode, 0, other.stdout[-3000:])
        self.assertIn("protal_db was built by an earlier run from the same release and protal; kept", other.stdout)
        self.assertRegex(other.stdout, r"Built \S+training_db in \d+ s")
        self.assertNotRegex(other.stdout, r"Built \S+protal_db")
        self.assertIn("was simulated from other inputs (or by an older collector): simulating it again",
                      self.text("out", "training_data.log"))
        self.assertFalse(os.path.exists(os.path.join(self.tmp.name, "out", ".converted")))

    def test_b_a_failed_background_build_stops_the_run(self):
        # The finished database's build fails in the background, 2 s in: the run stops then, not after the
        # collection, the training and the parity checks.
        failing = os.path.join(self.tmp.name, "failing_protal")
        with open(failing, "w") as fh:
            fh.write(f'#!/bin/sh\ncase "$*" in *--build*protal_db*) sleep 2; exit 3;; esac\nexec {os.environ["PROTAL"]} "$@"\n')
        os.chmod(failing, 0o755)
        started = time.time()
        result = self.build("out_fail", protal=failing)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Command failed (3)", result.stdout)
        self.assertIn("index_and_package.log", result.stdout)
        self.assertNotIn("Collected the training data", result.stdout)
        self.assertLess(time.time() - started, 600)

    def test_c_sigterm_stops_every_command(self):
        run_ = self.build("out_term", wait=False)
        log = os.path.join(self.tmp.name, "out_term", "training_data.log")
        deadline = time.time() + 900
        while not os.path.exists(log) and run_.poll() is None and time.time() < deadline:
            time.sleep(1)
        self.assertIsNone(run_.poll(), "the run ended before the collection")
        time.sleep(3)  # the collector has started the simulator
        run_.send_signal(signal.SIGTERM)
        output, _ = run_.communicate(timeout=120)
        self.assertNotEqual(run_.returncode, 0)
        self.assertIn("Stopped by SIGTERM", output)
        out = os.path.join(self.tmp.name, "out_term")
        left = []
        for pid in os.listdir("/proc"):
            try:
                with open(f"/proc/{pid}/cmdline", "rb") as fh:
                    if out.encode() in fh.read():
                        left.append(pid)
            except OSError:
                pass
        self.assertEqual(left, [], "commands of the stopped run are still running")


if __name__ == "__main__":
    unittest.main()
