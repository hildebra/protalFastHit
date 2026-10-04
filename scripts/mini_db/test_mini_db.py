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
import glob
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
GENE_NEIGHBOURS = os.path.join(HERE, "gene_neighbours.py")
SIMULATE_READS = os.path.join(HERE, "simulate_reads.py")
LINEAGES = os.path.join(HERE, "gtdb_like_lineages.py")
DOWNLOAD = os.path.join(HERE, "..", "download_gtdb.py")
RANK_GENES = os.path.join(HERE, "..", "rank_genes.py")
RELEASES_DOWNLOAD = os.path.join(HERE, "..", "download_gtdb_releases.py")
RELEASES_BUILD = os.path.join(HERE, "..", "build_gtdb_releases.py")

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


def full_reference(folder):
    """The content of a converted folder's full reference (full_reference.fna.zst, or .fna without zstd)."""
    sys.path.insert(0, HERE)
    from gtdb_to_protal_db import full_reference_path, read_full_reference
    with read_full_reference(full_reference_path(folder)) as fh:
        return fh.read()


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
            one, four = (os.path.join(self.tmp.name, f"db_{order}_t{t}") for t in (1, 4))
            run(CONVERT, "--gtdb", self.gtdb, "--outdir", one, "--order", order, "-t", "1")
            run(CONVERT, "--gtdb", self.gtdb, "--outdir", four, "--order", order, "-t", "4")
            for f in ("reference.fna", "reference.map", "internal_taxonomy.dmp"):
                with open(os.path.join(one, f), "rb") as a, open(os.path.join(four, f), "rb") as b:
                    self.assertEqual(a.read(), b.read(), f"{order} order, {f}")
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
        sys.path.insert(0, HERE)
        from gtdb_to_protal_db import CONVERTED_FILES, full_reference_path
        expected = {"reference.fna", "reference.map", os.path.basename(full_reference_path(self.db))} | \
            {f for f in CONVERTED_FILES if os.path.isfile(os.path.join(self.db, f))}
        self.assertEqual(set(os.listdir(dst)), expected)
        # Converting a release anew removes the build outputs of the folder's earlier reference, too.
        run(CONVERT, "--gtdb", self.gtdb, "--outdir", src)
        for name in ("unique_kmers.tsv", "gene_conservation.tsv", "index.prx.zst.partial", "database.protal.partial",
                     "database.protal"):
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

    def test_download_releases(self):
        # download_gtdb_releases.py: download_gtdb.py per release into INPUTS/gtdb_r<release>, the options it does
        # not know passed on, and a summary of what each release's inputs hold.
        mirror = os.path.join(self.tmp.name, "mirror_releases")
        gtdb_mirror(self.gtdb, mirror)

        class Quiet(http.server.SimpleHTTPRequestHandler):
            def log_message(self, *args):
                pass

        server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), functools.partial(Quiet, directory=mirror))
        threading.Thread(target=server.serve_forever, daemon=True).start()
        self.addCleanup(server.server_close)
        self.addCleanup(server.shutdown)
        datasets = os.path.join(self.tmp.name, "datasets_releases")
        with open(datasets, "w") as fh:
            fh.write(FAKE_DATASETS)
        os.chmod(datasets, 0o755)
        env = dict(os.environ, FAKE_TABLE=os.path.join(self.gtdb, "simulation", "genomes.tsv"))
        out = os.path.join(self.tmp.name, "inputs_releases")
        result = subprocess.run([sys.executable, RELEASES_DOWNLOAD, "-o", out, "--releases", "226", "-t", "2", "--mirror",
                                 f"http://127.0.0.1:{server.server_port}", "--datasets", datasets, "--species", "2",
                                 "--per_species", "1", "--rep_only_species", "0", "--batch", "4"],
                                env=env, capture_output=True, text=True)
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
        failing = subprocess.run([sys.executable, RELEASES_DOWNLOAD, "-o", out + "_x", "--releases", "999", "--mirror",
                                  f"http://127.0.0.1:{server.server_port}", "--datasets", datasets, "--no_genomes"],
                                 env=env, capture_output=True, text=True)
        self.assertEqual(failing.returncode, 1)
        self.assertIn("failed: r999", failing.stdout)

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

    def test_gene_conservation_in_the_build_metadata(self):
        # build_metadata.tsv records what protal --build said of the genes' conservation factors.
        sys.path.insert(0, os.path.join(HERE, ".."))
        import build_gtdb_database as build
        log = os.path.join(self.tmp.name, "index_and_package.log")
        with open(log, "w") as fh:
            fh.write("Uniqueness check took 1s\nGene conservation: factors 0.26-3.3 for 168 genes, from 120 species "
                     "(65559 copies compared): /data/db/gene_conservation.tsv\nGene conservation took 1s\n")
        self.assertEqual(build.gene_conservation_summary(log),
                         "factors 0.26-3.3 for 168 genes, from 120 species (65559 copies compared)")
        with open(log, "w") as fh:
            fh.write("Gene conservation: no factors (0 species with other genomes' copies of their genes in x.fna, 0 of "
                     "them with enough genes that differ from the representative's): every gene keeps the whole margin\n")
        self.assertTrue(build.gene_conservation_summary(log).startswith("no factors (0 species"))
        with open(log, "w") as fh:
            fh.write("Run build took 5s\n")
        self.assertEqual(build.gene_conservation_summary(log), "none (this protal does not estimate them)")
        self.assertEqual(build.gene_conservation_summary(log + ".missing"), "unknown (no build log)")

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

    def test_genome_table_lengths(self):
        # build_gtdb_database.py gives the simulator each genome's length (it would read every genome for it at
        # each design point): letters outside header lines, gzipped or not; a given table without lengths gets a
        # copy with them, one with them is taken as it is.
        sys.path.insert(0, os.path.join(HERE, ".."))
        import build_gtdb_database as build
        root = os.path.join(self.tmp.name, "lengths")
        os.makedirs(root)
        plain, zipped = os.path.join(root, "a.fna"), os.path.join(root, "b.fna.gz")
        text = ">a one 123 ACGT\r\nACGTNNacgt\r\n\r\nRYK-*\n>second\nAC GT\n"  # 10 + 3 + 4 letters
        with open(plain, "w", newline="") as fh:
            fh.write(text)
        with gzip.open(zipped, "wt", newline="") as fh:
            fh.write(text * 2)
        self.assertEqual(build.genome_length(plain), 17)
        self.assertEqual(build.genome_length(zipped), 34)
        given = os.path.join(root, "given.tsv")
        with open(given, "w") as fh:
            fh.write(f"name\ttaxonomy\tfasta_path\nGA\td__B;s__A\t{plain}\nGB\td__B;s__B\t{zipped}\n")
        copy = build.with_lengths(given, os.path.join(root, "genomes.tsv"), 2)
        with open(copy) as fh:
            self.assertEqual(fh.read(), f"GA\td__B;s__A\t{plain}\t17\nGB\td__B;s__B\t{zipped}\t34\n")
        self.assertEqual(build.with_lengths(copy, os.path.join(root, "other.tsv"), 2), copy)

    def test_collector_designs(self):
        # Read setups (built-in and custom ART profiles), abundance models, long-read setups and units.
        sys.path.insert(0, os.path.join(HERE, ".."))
        import argparse
        import collect_training_data as collect
        opts = argparse.Namespace(read_setups="150:HSXt:350:50,150:file=/p1+/p2:350:50,250:MSv3:550:50",
                                  read_pairs="1000,5000:2", read_types=["pe", "se", "ont"], samples=4, long_read_samples=0,
                                  pb_setup="errhmm:ERRHMM-SEQUEL:15000:3000:0.999",
                                  ont_setup="qshmm:QSHMM-ONT-HQ:8000:6000:0.97", long_read_bases="1e6,2e6,3e6:5")
        points = collect.design_points(opts)
        self.assertEqual([p["name"] for p in points], ["rl150_HSXt_p1000", "rl150_HSXt_p5000", "rl150_custom1_p1000",
                                                       "rl150_custom1_p5000", "rl250_p1000", "rl250_p5000"])
        # DEPTH:SAMPLES gives a depth's points other samples than --samples.
        self.assertEqual([p["samples"] for p in points], [4, 2, 4, 2, 4, 2])
        # Setups of one length and profile are told apart by their fragments; a setup or depth given twice
        # would share a folder.
        same = argparse.Namespace(read_setups="150:HS25:350:50,150:HS25:500:80", read_pairs="1000", samples=1)
        self.assertEqual([p["name"] for p in collect.design_points(same)], ["rl150_HS25_f350-50_p1000", "rl150_HS25_f500-80_p1000"])
        for setups, pairs in (("150:HS25:350:50,150:HS25:350:50", "1000"), ("150:HS25:350:50", "1000,1000"),
                              ("150:HS25:350:50", "1000,5000:0"), ("150:HS25:350:50", "1000:x")):
            with self.assertRaises(SystemExit):
                collect.design_points(argparse.Namespace(read_setups=setups, read_pairs=pairs, samples=1))
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
        hifi = collect.parse_long_setup("hifi:15000:3000:3")
        self.assertEqual((hifi["method"], hifi["length_mean"], hifi["length_sd"], hifi["q_sd"]), ("hifi", 15000, 3000, 3.0))
        with self.assertRaises(SystemExit):
            collect.parse_long_setup("art:x:1:1:1")
        _, units = collect.units_of(opts)
        self.assertEqual([u["type"] for u in units], ["pe"] * 6 + ["se"] * 6 + ["ont"] * 3)
        # A long-read point replays the communities of the paired-end points of the depth at its place (modulo the
        # depths), of every read setup: more samples than one point has (3:5 against 4 per point).
        ont = [u for u in units if u["type"] == "ont"]
        self.assertEqual([[p["name"] for p in u["communities"]] for u in ont],
                         [["rl150_HSXt_p1000", "rl150_custom1_p1000", "rl250_p1000"],
                          ["rl150_HSXt_p5000", "rl150_custom1_p5000", "rl250_p5000"],
                          ["rl150_HSXt_p1000", "rl150_custom1_p1000", "rl250_p1000"]])
        self.assertEqual([u["samples"] for u in ont], [4, 4, 5])
        self.assertEqual([u["samples"] for u in units if u["type"] == "se"], [4, 2, 4, 2, 4, 2])
        self.assertEqual(units[6]["name"], "rl150_HSXt_p1000_se")
        # No more long-read samples than the communities they replay (here 3 x 2 at 5,000 read pairs).
        with self.assertRaises(SystemExit):
            collect.units_of(argparse.Namespace(**{**vars(opts), "long_read_bases": "1e6,2e6:7"}))
        _, units = collect.units_of(argparse.Namespace(**{**vars(opts), "long_read_samples": 12}))
        self.assertEqual([u["samples"] for u in units if u["type"] == "ont"], [12, 6, 5])

    @staticmethod
    def fake_templ_pbsim(path):
        """A stand-in for pbsim3 --strategy templ: one read per template, its sequence as it is, named r_<n>
        after its place in the file, into <prefix>.fq.gz; it fails without --strategy templ."""
        with open(path, "w") as fh:
            fh.write("#!" + sys.executable + "\nimport gzip, sys\na = sys.argv[1:]\nget = lambda k: a[a.index(k) + 1]\n"
                     "if get('--strategy') != 'templ': sys.exit(2)\n"
                     "seqs = [l.strip() for l in open(get('--template')) if not l.startswith('>')]\n"
                     "with gzip.open(get('--prefix') + '.fq.gz', 'wt') as out:\n"
                     "    for i, s in enumerate(seqs, 1):\n"
                     "        out.write('@%s_%d\\n%s\\n+\\n%s\\n' % (get('--id-prefix'), i, s, 'I' * len(s)))\n")
        os.chmod(path, 0o755)

    def test_long_read_replay(self):
        # pb/ont samples replay a paired-end point's communities: the collector draws the reads (a genome by
        # relative abundance times length, a start uniform, either strand) until the sample's bases, and pbsim3
        # makes one read of each (a stand-in here, which keeps the sequence).
        sys.path.insert(0, os.path.join(HERE, ".."))
        import argparse
        import random
        import collect_training_data as collect
        root = os.path.join(self.tmp.name, "longreads")
        point = os.path.join(root, "points", "rl150_p1000")
        os.makedirs(os.path.join(point, "sim"))
        genomes, sequences = {}, {}
        rng = random.Random(5)
        for name, length in (("GA", 20000), ("GB", 40000)):
            genomes[name] = os.path.join(root, name + ".fna.gz")
            sequences[name] = "".join(rng.choice("ACGT") for _ in range(length))
            with gzip.open(genomes[name], "wt") as fh:
                fh.write(f">{name}_contig\n" + "\n".join(sequences[name][i:i + 70] for i in range(0, length, 70)) + "\n")
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
        self.fake_templ_pbsim(fake)
        models = os.path.join(root, "models")
        os.makedirs(models)
        open(os.path.join(models, "FAKE.model"), "w").close()
        opts = argparse.Namespace(out=root, seed=1, pbsim=fake, pbsim_models=models, samples=2)
        unit = {"type": "ont", "name": "ont_b300000", "bases": 300000,
                "setup": collect.parse_long_setup("qshmm:FAKE:1000:0:0.97"), "samples": 2,
                "communities": [{"name": "rl150_p1000"}],
                "point": {"name": "ont_b300000", "read_length": "1000", "read_pairs": "300000"}}
        keys = {"ont_b300000": {"a": 1}}
        self.assertEqual(collect.simulate_long([(0, unit)], opts, 2, keys), {})
        self.assertTrue(collect.same_key(os.path.join(root, "points", "ont_b300000", "simulated.json"), keys["ont_b300000"]))
        self.assertFalse(os.path.exists(os.path.join(root, "points", "ont_b300000", "sim", "tmp")))
        rows, _ = collect.unit_map_rows(unit, opts)
        self.assertEqual([r["SAMPLEID"] for r in rows], ["ont_b300000_s_1", "ont_b300000_s_2"])
        self.assertEqual([r["PROFILE_TRUTH"] for r in rows], ["/truth/1", "/truth/2"])
        self.assertEqual({r["READ_TYPE"] for r in rows} | {r["SECOND"] for r in rows}, {"ont", "-"})
        complement = str.maketrans("ACGT", "TGCA")
        first = {}
        for row, (a, b) in zip(rows, ((0.5, 0.5), (0.8, 0.2))):
            with gzip.open(row["FIRST"], "rt") as fh:
                lines = fh.read().splitlines()
            first[row["SAMPLEID"]] = lines
            names, reads = lines[0::4], lines[1::4]
            self.assertEqual(len(names), len(set(names)), "read names must be unique within a sample")
            # The reads' bases reach the sample's; a read is 1 kb (the setup's mean, SD 0) or shorter at a
            # contig's end, from either strand of its genome.
            self.assertGreaterEqual(sum(map(len, reads)), 300000)
            self.assertLess(sum(map(len, reads)), 300000 + 1000)
            self.assertTrue(all(100 <= len(r) <= 1000 for r in reads))
            strands = collections.Counter()
            for name, read in zip(names, reads):
                genome = sequences["GA" if name.startswith("@g0x_") else "GB"]
                strands["+" if read in genome else "-" if read.translate(complement)[::-1] in genome else "?"] += 1
            self.assertEqual(strands["?"], 0)
            self.assertGreater(min(strands["+"], strands["-"]), len(reads) / 3)
            counts = {g: sum(n.startswith(f"@g{g}x_") for n in names) for g in (0, 1)}
            share_a = a * 20000 / (a * 20000 + b * 40000)
            self.assertEqual(counts[0] + counts[1], len(names))
            self.assertAlmostEqual(counts[0] / len(names), share_a, delta=0.1)
        self.assertTrue(collect.simulated(unit, opts))
        # The same seed, the same reads.
        self.assertEqual(collect.simulate_long([(0, unit)], opts, 1), {})
        for row in rows:
            with gzip.open(row["FIRST"], "rt") as fh:
                self.assertEqual(fh.read().splitlines(), first[row["SAMPLEID"]])
        # A sample of more bases than --long_read_chunk is simulated in chunks side by side, joined into one gzip
        # file: every read named once (chunk c of k names every k-th from c + 1), the sample's bases, reads of its
        # genomes; the same reads on any number of slots, and its chunks' files gone.
        chunked = argparse.Namespace(**vars(opts), long_read_chunk=70000)
        self.assertEqual(collect.simulate_long([(0, unit)], chunked, 3), {})
        for row in rows:
            with gzip.open(row["FIRST"], "rt") as fh:
                lines = fh.read().splitlines()
            names, reads = lines[0::4], lines[1::4]
            self.assertNotEqual(lines, first[row["SAMPLEID"]])
            ids = [int(n.split("x_")[1]) for n in names]
            self.assertEqual(len(set(ids)), len(ids))
            self.assertEqual({i % 5 for i in ids}, set(range(5)))  # 5 chunks, each its own every 5th name
            self.assertGreaterEqual(sum(map(len, reads)), 300000)
            self.assertLess(sum(map(len, reads)), 300000 + 5 * 1000)
            for name, read in zip(names, reads):
                genome = sequences["GA" if name.startswith("@g0x_") else "GB"]
                self.assertTrue(read in genome or read.translate(complement)[::-1] in genome)
            with gzip.open(row["FIRST"], "rt") as fh:
                chunked_first = fh.read()
            self.assertEqual(collect.simulate_long([(0, unit)], chunked, 1), {})
            with gzip.open(row["FIRST"], "rt") as fh:
                self.assertEqual(fh.read(), chunked_first)
        self.assertFalse(os.path.exists(os.path.join(root, "points", "ont_b300000", "sim", "tmp")))
        tasks = collect.long_read_chunks({"bases": 300000, "seed": 7, "tmp": "/t", "out": "/o"}, 70000)
        self.assertEqual([t["bases"] for t in tasks], [60000] * 4 + [60000])
        self.assertEqual(len({t["seed"] for t in tasks}), 5)
        self.assertEqual(collect.long_read_chunks({"bases": 300000, "seed": 7}, 0)[0]["seed"], 7)
        self.assertEqual(len(collect.long_read_chunks({"bases": 300000, "seed": 7}, 300000)), 1)

    def test_scheduler(self):
        # The simulations' queue: the ready job of highest priority first, on its slots; a job waits for those it
        # comes after; a job's new jobs join; a failure starts nothing more and says why.
        sys.path.insert(0, os.path.join(HERE, ".."))
        import collect_training_data as collect
        order, lock = [], threading.Lock()

        def job(name, seconds=0.05, new=(), error=None):
            def run():
                with lock:
                    order.append(name)
                time.sleep(seconds)
                return error, list(new)
            return run
        s = collect.Scheduler(2)
        s.add("low", job("low"), priority=1)
        s.add("high", job("high"), priority=5)
        s.add("big", job("big"), need=2, priority=3)
        s.add("after", job("after", new=[{"name": "child", "run": job("child")}]), after=["big"], priority=9)
        self.assertEqual(s.run(), {})
        self.assertEqual(order[:3], ["high", "low", "big"])  # low fills the slot big cannot use yet
        self.assertLess(order.index("big"), order.index("after"))
        self.assertLess(order.index("after"), order.index("child"))
        self.assertEqual(sorted(order), ["after", "big", "child", "high", "low"])
        order.clear()
        s = collect.Scheduler(1)
        s.add("bad", job("bad", error="it broke"), priority=2)
        s.add("next", job("next"), priority=1)
        s.add("waits", job("waits"), after=["bad"])
        failures = s.run()
        self.assertEqual(failures, {"bad": "it broke"})
        self.assertEqual(order, ["bad"])
        # The paired-end points' threads: all slots between them by their work, at least one each, at most 16 per
        # sample.
        points = [(0, {"name": "deep", "read_pairs": "10000000", "read_length": "150", "samples": 2}),
                  (1, {"name": "shallow", "read_pairs": "1000", "read_length": "150", "samples": 12}),
                  (2, {"name": "mid", "read_pairs": "500000", "read_length": "150", "samples": 12})]
        threads = collect.pe_threads(points, 30)
        self.assertEqual(sum(threads.values()), 30)
        self.assertGreater(threads["deep"], threads["mid"])
        self.assertGreaterEqual(threads["shallow"], 1)
        self.assertEqual(collect.pe_threads(points[:1], 64), {"deep": 32})

    def test_long_read_templates(self):
        # A long-read sample's reads come from the contigs of 100 bases or more (pbsim3's shortest read; it stops
        # at a shorter reference sequence), have the setup's lengths, and a failure says why.
        sys.path.insert(0, os.path.join(HERE, ".."))
        import random
        import collect_training_data as collect
        root = os.path.join(self.tmp.name, "templates")
        os.makedirs(root)
        rng = random.Random(3)
        records = [(">a one", "".join(rng.choice("acgt") for _ in range(3000))), (">short", "C" * 99),
                   (">exact", "G" * 100), (">tiny", "ACGT")]
        text = "".join(f"{header}\r\n" + "\r\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)) + "\r\n\r\n"
                       for header, seq in records)
        plain, zipped = os.path.join(root, "plain.fna"), os.path.join(root, "zipped.fna.gz")
        with open(plain, "w", newline="") as fh:
            fh.write(text)
        with gzip.open(zipped, "wt", newline="") as fh:
            fh.write(text)
        expected = [seq.upper().encode() for _, seq in records]
        self.assertEqual(collect.read_contigs(plain), expected)
        self.assertEqual(collect.read_contigs(zipped), expected)
        fake = os.path.join(root, "pbsim")
        self.fake_templ_pbsim(fake)

        def task(name, fasta, pbsim=fake, bases=50000, setup="qshmm:FAKE:1000:0:0.97"):
            return {"sample": name, "out": os.path.join(root, name + ".fq.gz"), "bases": bases, "pbsim": pbsim,
                    "setup": collect.parse_long_setup(setup), "model": "FAKE", "seed": 1, "tmp": os.path.join(root, "tmp", name),
                    "genomes": [{"genome": "G", "fasta": fasta, "weight": 1.0}]}

        long_one, exact = records[0][1].upper(), records[2][1]
        for name, fasta in (("plain", plain), ("zipped", zipped)):
            self.assertIsNone(collect.long_read_sample(task(name, fasta)))
            with gzip.open(os.path.join(root, name + ".fq.gz"), "rt") as fh:
                reads = fh.read().splitlines()[1::4]
            self.assertTrue(reads)
            for read in reads:  # from "a one" (up to 1 kb, cut at its end) or the whole of "exact", either strand
                back = read.translate(str.maketrans("ACGT", "TGCA"))[::-1]
                self.assertTrue(read in long_one or back in long_one or exact in (read, back), read[:20])
                self.assertTrue(100 <= len(read) <= 1000)
        # Only short sequences: nothing to simulate from. A pbsim failure carries its last line of output, and
        # a pbsim that makes fewer reads than templates is caught.
        short = os.path.join(root, "short.fna")
        with open(short, "w") as fh:
            fh.write(">tiny\nACGT\n>short\n" + "A" * 99 + "\n")
        self.assertIn("no sequence of 100 bases or more", collect.long_read_sample(task("short", short)))
        failing = os.path.join(root, "failing")
        with open(failing, "w") as fh:
            fh.write("#!/bin/sh\necho 'ERROR: out of ideas'\nexit 3\n")
        os.chmod(failing, 0o755)
        error = collect.long_read_sample(task("failing", plain, failing))
        self.assertIn("pbsim failed (3)", error)
        self.assertIn("ERROR: out of ideas", error)
        lossy = os.path.join(root, "lossy")
        with open(lossy, "w") as fh:
            fh.write("#!" + sys.executable + "\nimport gzip, sys\na = sys.argv[1:]\n"
                     "with gzip.open(a[a.index('--prefix') + 1] + '.fq.gz', 'wt') as out:\n"
                     "    out.write('@r_1\\nACGT\\n+\\nIIII\\n')\n")
        os.chmod(lossy, 0o755)
        self.assertIn("pbsim made 1 reads of", collect.long_read_sample(task("lossy", plain, lossy)))
        # Read lengths: the setup's gamma distribution between 100 and 1,000,000; the mean when the SD is 0.
        rng = random.Random(1)
        lengths = [collect.read_length(rng, 8000, 6000) for _ in range(20000)]
        self.assertTrue(all(100 <= n <= 1000000 for n in lengths))
        self.assertAlmostEqual(sum(lengths) / len(lengths), 8000, delta=200)
        self.assertEqual(collect.read_length(rng, 1000, 0), 1000)
        # A hifi setup makes the reads with hifi_reads.py, no pbsim3 needed: one of each template, named after it,
        # with qualities of HiFi reads.
        missing = os.path.join(root, "no_pbsim")
        self.assertIsNone(collect.long_read_sample(task("hifi", plain, missing, setup="hifi:1000:0:3")))
        with gzip.open(os.path.join(root, "hifi.fq.gz"), "rt") as fh:
            lines = fh.read().splitlines()
        # Named after their templates, g<genome>x_<n> in the order drawn; the templates go once the reads are there.
        self.assertEqual([line[1:] for line in lines[0::4]], [f"g0x_{i}" for i in range(1, len(lines) // 4 + 1)])
        self.assertFalse(os.path.exists(os.path.join(root, "tmp", "hifi")))
        self.assertFalse(os.path.exists(os.path.join(root, "tmp", "plain")), "a pbsim3 sample's temporary files too")
        quality = [ord(c) - 33 for line in lines[3::4] for c in line]
        self.assertGreater(sorted(quality)[len(quality) // 2], 25)
        for bad in ("hifi:1000:0", "hifi:1000:0:30:3", "hifi:a:0:3"):
            with self.assertRaises(SystemExit):
                collect.parse_long_setup(bad)

    def test_hifi_reads(self):
        # hifi_reads.py: a read of each template whose qualities say how many errors it has (their expected errors
        # are its errors, overall), reads of Q50 at 5 kb, Q30 at 25 kb and Q20 at 50 kb (SD 3 around), only indels in
        # homopolymers, the same reads for the same seed.
        sys.path.insert(0, os.path.join(HERE, ".."))
        import numpy as np
        import hifi_reads
        self.assertEqual(list(hifi_reads.length_quality([1000, 5000, 15000, 25000, 37500, 50000, 75000])),
                         [50, 50, 40, 30, 25, 20, 10])
        rng = np.random.default_rng(3)
        seqs, groups = [], []
        for length, count in ((5000, 100), (15000, 100), (25000, 100), (50000, 40)):
            for _ in range(count):  # runs of random bases, a tenth of them homopolymers of 3 to 8
                runs = np.where(rng.random(length) < 0.9, rng.geometric(0.6, length), rng.integers(3, 9, length))
                seqs.append(np.repeat(np.frombuffer(b"ACGT", np.uint8)[rng.integers(0, 4, length)], runs)[:length].tobytes())
                groups.append(length)
        groups = np.array(groups)
        reads, stats = hifi_reads.mutate(seqs, np.random.default_rng(1), 3)
        self.assertEqual(len(reads), len(seqs))
        self.assertAlmostEqual(stats["events"].sum() / stats["expected"].sum(), 1.0, delta=0.08)
        self.assertGreaterEqual(stats["q"].min(), hifi_reads.Q_MIN)
        for length, q in ((5000, 50), (15000, 40), (25000, 30), (50000, 20)):
            self.assertAlmostEqual(float(stats["q"][groups == length].mean()), q, delta=1.5, msg=length)
            self.assertAlmostEqual(float(stats["q"][groups == length].std()), 3, delta=1.0, msg=length)
        self.assertEqual(stats["homopolymer_substitutions"], 0)
        for (read, quality), seq in zip(reads, seqs):
            self.assertEqual(len(read), len(quality))
            self.assertLess(abs(len(read) - len(seq)), 0.05 * len(seq))
            self.assertTrue(34 <= min(quality) and max(quality) <= 126)  # Q1 to Q93
        again, _ = hifi_reads.mutate(seqs, np.random.default_rng(1), 3)
        self.assertEqual(again, reads)
        # A read's errors follow its quality: about 10^(-Q/10) per base.
        errors = stats["events"] / np.array([len(s) for s in seqs])
        for length, q in ((25000, 30), (50000, 20)):
            self.assertAlmostEqual(float(np.median(errors[groups == length])), 10 ** (-q / 10), delta=0.4 * 10 ** (-q / 10))

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


class GeneNeighboursTest(unittest.TestCase):
    """gene_neighbours.py on a release whose markers lie in operon-like clusters (simulate_gtdb_release.py
    --operons), 4 species in one family and 2 in another, two genomes each in 4 contigs: every gene of the
    database is found in every genome of its species where the simulator put it (a representative's exactly, the
    other strain's by its k-mer trace, a few bases off at most), and the lines of each family count exactly the
    neighbours that its species' marker positions give (each species the partner most of its genomes show there,
    worked out here independently); so do the lines derived from the positions alone, and those a training
    database's copy derives without a species."""

    MAX_GAP = 3000
    SLACK = 15  # bases a strain's gene placed by its trace may be off: its codon indels

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        lineages = os.path.join(cls.tmp.name, "lineages.txt")
        order = "d__Bacteria;p__Simulatota;c__Simulatia;o__Simulales"
        with open(lineages, "w") as fh:
            fh.writelines(f"{order};f__Simulaceae;g__Mockella;s__Mockella s{i}\n" for i in range(4))
            fh.writelines(f"{order};f__Otheraceae;g__Otherella;s__Otherella s{i}\n" for i in range(2))
        cls.gtdb = os.path.join(cls.tmp.name, "gtdb")
        cls.db = os.path.join(cls.tmp.name, "db")
        run(SIMULATE, "--outdir", cls.gtdb, "--lineages", lineages, "--operons", "--operon_breaks", "0.5",
            "--genome_length", "400000", "--genomes_per_species", "2", "--contigs", "4")
        run(CONVERT, "--gtdb", cls.gtdb, "--outdir", cls.db)
        cls.positions = os.path.join(cls.tmp.name, "positions.tsv")
        cls.output = subprocess.run([sys.executable, GENE_NEIGHBOURS, "--db", cls.db, "--genome_table",
                                     os.path.join(cls.gtdb, "simulation", "genomes.tsv"), "--positions", cls.positions,
                                     "-t", "2"], check=True, capture_output=True, text=True).stdout

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_genome_index_finds_what_find_finds(self):
        # place() looks genes up in an index of the genome's k-mers at a stride instead of scanning the genome with
        # find: it must find every occurrence find does, overlapping ones and those at contig ends too.
        sys.path.insert(0, os.path.join(HERE))
        import gene_neighbours as gn
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
        self.assertIn("12 genomes of 6 species", self.output)
        self.assertIn("(6 with their representative genome)", self.output)
        self.assertIn("circular: 0 with a whole replicon by its header, 0 representatives in one sequence", self.output)

    def test_family_lines_count_the_neighbours_of_its_species(self):
        rows = self.check_family_lines("gene_neighbours.tsv")
        # The two families differ: a cluster that one of them broke up.
        by_family = collections.defaultdict(set)
        for family, gene, end, partner, partner_end in rows:
            if partner:
                by_family[family].add((gene, end, partner, partner_end))
        self.assertEqual(len(by_family), 2)
        a, b = by_family.values()
        self.assertTrue(a - b or b - a)
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
        self.assertEqual(informative, 6)  # an end informative in all six species of the order
        self.assertIn("family: 2 clades", self.output)
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
        shutil.copy(self.positions, os.path.join(src, "gene_positions.tsv"))
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
        sys.path.insert(0, HERE)
        import gene_neighbours as gn
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
        sys.path.insert(0, HERE)
        import gene_neighbours as gn
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


@unittest.skipUnless(os.environ.get("PROTAL") and os.environ.get("SIMULATE") and shutil.which("art_illumina"),
                     "needs $PROTAL, $SIMULATE (simulate_metagenomes) and art_illumina")
class TraceRelativesTest(unittest.TestCase):
    """trace_relatives.py: where the reads of a species the training database lacks land, by the genes' factors."""

    def test_reads_by_gene_factor(self):
        with tempfile.TemporaryDirectory() as tmp:
            point = os.path.join(tmp, "points", "rl100_p1000")
            os.makedirs(os.path.join(point, "sim"))
            os.makedirs(os.path.join(point, "protal", "alignments"))
            lineage = "d__Bacteria;p__P;c__C;o__O;f__F;g__G;s__G {}"
            with open(os.path.join(point, "sim", "manifest.tsv"), "w") as fh:
                fh.write("sample\tgenome\tspecies\ttaxonomy\tvertical_coverage\n")
                fh.write(f"rl100_p1000_s_1\tGCF_1.1\tG a\t{lineage.format('a')}\t2.0\n")  # in the database (taxon 1)
                fh.write(f"rl100_p1000_s_1\tGCA_2.1\tG b\t{lineage.format('b')}\t1.0\n")  # held out: lands on taxon 3
            db = os.path.join(tmp, "training_db")
            os.makedirs(db)
            with open(os.path.join(db, "genome2tiid.tsv"), "w") as fh:
                for taxid, name in ((1, "a"), (3, "c")):
                    fh.write(f"GCF_{taxid}.1\t{taxid}\tGCF_{taxid}.1\t{lineage.format(name)}\n")
            with open(os.path.join(db, "gene_congeners.tsv"), "w") as fh:
                fh.write("geneid\twithin_factor\tbetween_factor\tpairs\tspecies\tidentical_share\tnear_identical_share\n"
                         "1\t0.5\t0.4\t3\t2\t0\t0\n2\t1.5\t1.6\t3\t2\t0\t0\n")
            heldout = os.path.join(tmp, "heldout_species.txt")
            with open(heldout, "w") as fh:
                fh.write("s__G b\tspecies\ts__G b\n")
            records = []
            for i in range(8):  # the species' own: 4 records on each gene
                records.append((f"GCF_1.1_contig1-{i}", f"1_{1 + i % 2}", 60))
            for i, mapq in enumerate((60, 60, 0, 1)):  # the relative: 4 on the conserved gene, 2 of them ambiguous
                records.append((f"GCA_2.1_contig1-{i}", "3_1", mapq))
            records.append(("GCA_2.1_contig1-9", "3_2", 60))  # and 1 on the fast one
            records.append(("GCF_1.1_contig1-20", "3_1", 60))  # the species' own read on a congener: not counted
            with open(os.path.join(point, "protal", "alignments", "rl100_p1000_s_1.sam"), "w") as fh:
                fh.write("@HD\tVN:1.6\n")
                for name, ref, mapq in records:
                    fh.write(f"{name}\t0\t{ref}\t1\t{mapq}\t100M\t*\t0\t0\t{'A' * 100}\t{'I' * 100}\n")
                fh.write("GCA_2.1_contig1-5\t256\t1_1\t1\t0\t100M\t*\t0\t0\t*\t*\n")  # secondary: not counted
            out = os.path.join(tmp, "logs", "relatives")
            result = subprocess.run([sys.executable, os.path.join(HERE, "..", "trace_relatives.py"), "--points",
                                     os.path.join(tmp, "points"), "--db", db, "--heldout", heldout, "--out", out],
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            with open(out + ".tsv") as fh:
                genes = {r["geneid"]: r for r in csv.DictReader(fh, delimiter="\t")}
            # Per unit coverage: own 4 / 2 on each gene; the relative 4 / 1 on gene 1 (2 kept), 1 / 1 on gene 2.
            self.assertAlmostEqual(float(genes["1"]["R"]), 2.0)
            self.assertAlmostEqual(float(genes["1"]["R_kept"]), 1.0)
            self.assertAlmostEqual(float(genes["2"]["R"]), 0.5)
            self.assertAlmostEqual(float(genes["1"]["relative_mapq_below_4"]), 0.5)
            self.assertAlmostEqual(float(genes["1"]["relative_on_congener"]), 1.0)
            with open(out + ".txt") as fh:
                text = fh.read()
            self.assertIn("1 paired-end training samples", text)
            self.assertIn("| factor < 0.7 | 1 | 2.000 | 1.000 |", text)

            # Without factors there is nothing to trace by.
            os.remove(os.path.join(db, "gene_congeners.tsv"))
            result = subprocess.run([sys.executable, os.path.join(HERE, "..", "trace_relatives.py"), "--points",
                                     os.path.join(tmp, "points"), "--db", db, "--heldout", heldout, "--out", out],
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertIn("Not traced: the training database has no gene conservation factors", result.stdout)


class BuildOptionsTest(unittest.TestCase):
    def test_leaves_per_read_type(self):
        # build_gtdb_database.py --maxnodes: N for the read types not named, TYPE:N for one; 256 without either.
        sys.path.insert(0, os.path.join(HERE, ".."))
        import build_gtdb_database as build
        self.assertEqual([build.max_leaves("512,pb:128,ont:128", t) for t in ("pe", "se", "pb", "ont")],
                         [512, 512, 128, 128])
        self.assertEqual(build.max_leaves("64", "ont"), 64)
        self.assertEqual(build.max_leaves("se:32", "se"), 32)
        self.assertEqual(build.max_leaves("se:32", "pe"), 256)
        for bad in ("512,xx:4", "pb:", "many"):
            with self.assertRaises(ValueError):
                build.max_leaves(bad, "pe")


class BinaryCheckTest(unittest.TestCase):
    """The commit a build records for --version (protal_commit.cmake), and build_gtdb_database.py's check at its start
    that protal and the simulator were built from the source its scripts are at, on a throwaway checkout."""

    def setUp(self):
        if not shutil.which("git"):
            self.skipTest("no git")
        self.tmp = tempfile.TemporaryDirectory()
        self.source = os.path.join(self.tmp.name, "checkout")
        os.makedirs(os.path.join(self.source, "src"))
        self.write("CMakeLists.txt", "project(protal VERSION 0.7.3)\n")
        self.write("src/a.h", "int a;\n")
        self.git("init", "-q")
        self.commit("one")
        self.built = self.git("rev-parse", "HEAD")

    def tearDown(self):
        self.tmp.cleanup()

    def write(self, path, text, mode="w"):
        with open(os.path.join(self.source, path), mode) as fh:
            fh.write(text)

    def git(self, *arguments):
        return subprocess.run(["git", "-C", self.source, "-c", "user.name=test", "-c", "user.email=test@example.org",
                               *arguments], check=True, capture_output=True, text=True).stdout.strip()

    def commit(self, message):
        self.git("add", "-A")
        self.git("commit", "-q", "-m", message)

    def check(self, said, code=0):
        # A stand-in binary that says `said` to --version.
        sys.path.insert(0, os.path.join(HERE, ".."))
        import build_gtdb_database as build
        binary = os.path.join(self.tmp.name, "protal")
        with open(binary, "w") as fh:
            fh.write(f"#!/bin/sh\necho '{said}'\nexit {code}\n")
        os.chmod(binary, 0o755)
        return build.build_check(binary, "protal", self.source)

    def test_the_commit_in_the_build(self):
        if not shutil.which("cmake"):
            self.skipTest("no cmake")
        header = os.path.join(self.tmp.name, "protal_commit.h")

        def written():
            subprocess.run(["cmake", f"-DSOURCE_DIR={self.source}", f"-DOUTPUT={header}", "-P",
                            os.path.join(HERE, "..", "..", "protal_commit.cmake")], check=True, capture_output=True)
            with open(header) as fh:
                return fh.read()

        def expected(commit, modified):
            return f'#pragma once\n#define PROTAL_GIT_COMMIT "{commit}"\n#define PROTAL_GIT_MODIFIED {modified}\n'

        self.assertEqual(written(), expected(self.built, 0))
        self.write("README.md", "docs\n")  # not what the binaries are built from
        self.assertEqual(written(), expected(self.built, 0))
        self.write("src/b.h", "int b;\n")  # a new file in src/ is
        self.assertEqual(written(), expected(self.built, 1))
        shutil.rmtree(os.path.join(self.source, ".git"))
        self.assertEqual(written(), expected("", 0))

    def test_binaries_of_another_source_stop_the_run(self):
        self.assertEqual(self.check(f"protal v0.7.3 (commit {self.built})"), (None, None))
        # Commits that leave src/, lib/ and the build files alone keep the binary current.
        self.write("README.md", "docs\n")
        self.commit("docs")
        self.assertEqual(self.check(f"protal v0.7.3 (commit {self.built})"), (None, None))
        problem, _ = self.check("protal v0.7.1")
        self.assertIn("is v0.7.1, these scripts v0.7.3", problem)
        problem, _ = self.check("unknown option(s): --version", 1)  # simulate_metagenomes before 0.7.3
        self.assertIn("--version gives no version (1: unknown option", problem)
        problem, _ = self.check(f"protal v0.7.3 (commit {'0' * 40})")
        self.assertIn(f"was built from commit 0000000000, which {self.source} does not have", problem)
        self.write("src/a.h", "int a2;\n")
        problem, _ = self.check(f"protal v0.7.3 (commit {self.built})")
        self.assertIn(f"was built from commit {self.built[:10]}, and src/, lib/ or the build files have changed since "
                      "(uncommitted changes): rebuild it", problem)
        self.commit("two")
        problem, _ = self.check(f"protal v0.7.3 (commit {self.built})")
        self.assertIn("have changed since (1 commit): rebuild it", problem)
        self.assertIsNone(self.check(f"protal v0.7.3 (commit {self.git('rev-parse', 'HEAD')})")[0])

    def test_binaries_that_cannot_say_are_noted(self):
        problem, note = self.check("protal v0.7.3")
        self.assertIsNone(problem)
        self.assertIn("does not say which commit it was built from", note)
        problem, note = self.check(f"protal v0.7.3 (commit {self.built}, with uncommitted changes)")
        self.assertIsNone(problem)
        self.assertIn("was built with uncommitted changes to its source", note)
        shutil.rmtree(os.path.join(self.source, ".git"))
        problem, note = self.check(f"protal v0.7.3 (commit {self.built})")
        self.assertIsNone(problem)
        self.assertIn("is no git checkout", note)


def read_table(path):
    """A tab-separated table with a header line (# lines skipped) -> [row dict]."""
    with open(path) as fh:
        rows = [line.rstrip("\n").split("\t") for line in fh if not line.startswith("#")]
    return [dict(zip(rows[0], r)) for r in rows[1:]]


class GeneSubsetTest(unittest.TestCase):
    """A database of a subset of the marker genes. The converter's --genes (with --from_db, and with --gtdb the
    same files) keeps the genes named, by GTDB marker id or protal gene id, with their ids in reference.fna,
    reference.map, the full reference and gene2geneid.tsv, and counts the gene neighbours anew over them: a gene
    whose neighbour in the genome is left out faces the nearest gene kept, as a read of the reduced database
    meets it. With $PROTAL: protal --build of such a folder lists only its genes in unique_kmers.tsv and
    gene_conservation.tsv; rank_genes.py ranks a full build's genes and writes a gene list that the converter and
    --build_gene_subset take; and --build_gene_subset refuses a folder whose neighbours were counted over every
    gene."""

    MAX_GAP = 3000

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        lineages = os.path.join(cls.tmp.name, "lineages.txt")
        order = "d__Bacteria;p__Simulatota;c__Simulatia;o__Simulales"
        with open(lineages, "w") as fh:
            fh.writelines(f"{order};f__Simulaceae;g__Mockella;s__Mockella s{i}\n" for i in range(4))
            fh.writelines(f"{order};f__Otheraceae;g__Otherella;s__Otherella s{i}\n" for i in range(2))
            fh.writelines(f"d__Archaea;p__Archota;c__Archia;o__Archales;f__Archaceae;g__Archella;s__Archella s{i}\n"
                          for i in range(2))
        cls.gtdb = os.path.join(cls.tmp.name, "gtdb")
        cls.db = os.path.join(cls.tmp.name, "db")
        run(SIMULATE, "--outdir", cls.gtdb, "--lineages", lineages, "--operons", "--operon_breaks", "0.5",
            "--genome_length", "200000", "--genomes_per_species", "2", "--contigs", "2")
        run(CONVERT, "--gtdb", cls.gtdb, "--outdir", cls.db)
        run(GENE_NEIGHBOURS, "--db", cls.db, "--genome_table", os.path.join(cls.gtdb, "simulation", "genomes.tsv"), "-t", "2")
        sys.path.insert(0, HERE)
        from gtdb_to_protal_db import read_gene_ids
        cls.gene_ids = read_gene_ids(os.path.join(cls.db, "gene2geneid.tsv"))
        cls.markers = {gid: marker for marker, gid in cls.gene_ids.items()}
        # The subset: every other gene along the longest contig of one representative (the genes in between are
        # left out, so their neighbours in the subset's tables are the next genes kept): the first two by their
        # marker id, the third by the id without its version, the rest by protal gene id.
        with open(os.path.join(cls.db, "internal_taxonomy.dmp")) as fh:
            next(fh)
            cls.nodes = {int(r[0]): r for r in (l.rstrip("\n").split("\t") for l in fh)}
        cls.taxid, cls.rep = next((t, r[6]) for t, r in cls.nodes.items() if r[3] == "s__Mockella s0")
        by_contig = collections.defaultdict(list)
        for r in read_table(os.path.join(cls.db, "gene_positions.tsv")):
            if r["accession"] == cls.rep:
                by_contig[r["contig"]].append((int(r["start"]), int(r["end"]), int(r["gene"])))
        along = sorted(max(by_contig.values(), key=len))
        # Three genes in a row within the neighbour gap: the first and the third are kept, the middle one left
        # out, so that in the subset the first faces the third; four more genes from elsewhere on the contig.
        first = next(i for i in range(len(along) - 2)
                     if along[i + 1][0] - along[i][1] <= cls.MAX_GAP and along[i + 2][0] - along[i + 1][1] <= cls.MAX_GAP
                     and along[i + 2][0] - along[i][1] <= cls.MAX_GAP)
        others = [k for k in range(len(along)) if k not in (first, first + 1, first + 2)]
        indices = [first, first + 2] + others[::max(1, len(others) // 4)][:4]
        cls.subset = sorted(along[k][2] for k in indices)
        tokens = [cls.markers[g] for g in cls.subset[:2]] + [cls.markers[cls.subset[2]].split(".")[0]] + \
            [str(g) for g in cls.subset[3:]]
        cls.spec = ",".join(tokens)
        cls.copy = os.path.join(cls.tmp.name, "subset")
        cls.derived = subprocess.run([sys.executable, CONVERT, "--from_db", cls.db, "--genes", cls.spec, "--outdir", cls.copy],
                                     check=True, capture_output=True, text=True).stderr

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    @staticmethod
    def map_rows(folder):
        with open(os.path.join(folder, "reference.map")) as fh:
            return [tuple(map(int, line.split("\t"))) for line in fh]

    @staticmethod
    def records(data):
        """[(gene id, header, sequence)] of a FASTA's bytes, one sequence line per record."""
        lines = data.split(b"\n")
        return [(int(lines[i][1:].split(b"_")[1]), lines[i], lines[i + 1]) for i in range(0, len(lines) - 1, 2)]

    def same_file(self, a, b, name):
        with open(os.path.join(a, name), "rb") as x, open(os.path.join(b, name), "rb") as y:
            self.assertEqual(x.read(), y.read(), name)

    def test_the_copy_holds_the_genes_named_with_their_ids(self):
        self.assertEqual(len(self.subset), 6)
        rows = self.map_rows(self.copy)
        self.assertEqual(sorted({r[1] for r in rows}), self.subset)
        self.assertEqual([r[:2] for r in rows], [r[:2] for r in self.map_rows(self.db) if r[1] in self.subset])
        with open(os.path.join(self.copy, "reference.fna"), "rb") as fh:
            data = fh.read()
        for tid, gid, start, end in rows:
            self.assertEqual(data[data.rfind(b">", 0, start):start], f">{tid}_{gid}\n".encode())
            self.assertEqual(data[end:end + 1], b"\n")
            self.assertRegex(data[start:end].decode(), r"^[ACGT]+$")
        with open(os.path.join(self.db, "reference.fna"), "rb") as fh:
            self.assertEqual(self.records(data), [r for r in self.records(fh.read()) if r[0] in self.subset])
        self.assertEqual(self.records(full_reference(self.copy)),
                         [r for r in self.records(full_reference(self.db)) if r[0] in self.subset])
        with open(os.path.join(self.copy, "gene2geneid.tsv")) as fh:
            kept = {m: int(g) for m, g in (line.rstrip("\n").split("\t") for line in fh)}
        self.assertEqual(sorted(kept.values()), self.subset)
        self.assertEqual(kept, {m: g for m, g in self.gene_ids.items() if g in self.subset})
        for name in ("internal_taxonomy.dmp", "species_priors.tsv", "genome2tiid.tsv"):
            self.same_file(self.db, self.copy, name)
        # The genomes with a gene kept: the archaea have none of the bacterial genes chosen here.
        self.assertRegex(self.derived, r"gene_neighbours\.tsv derived from the 1[2-6] genomes of [6-8] species kept: \d+ lines "
                                       r"\(over the 6 genes kept\)")
        self.assertRegex(self.derived, r"without 0 species, 6 genes kept: \d+ representative sequences kept")

    def test_a_direct_conversion_gives_the_same_files(self):
        direct = os.path.join(self.tmp.name, "direct")
        run(CONVERT, "--gtdb", self.gtdb, "--outdir", direct, "--genes", self.spec)
        for name in ("reference.fna", "reference.map", "gene2geneid.tsv", "internal_taxonomy.dmp", "species_priors.tsv"):
            self.same_file(direct, self.copy, name)
        self.assertEqual(full_reference(direct), full_reference(self.copy))
        # A file lists the genes too (its first column), with or without the species left out.
        listed = os.path.join(self.tmp.name, "genes.txt")
        with open(listed, "w") as fh:
            fh.write("# the subset\n" + "".join(f"{t}\tcomment\n" for t in self.spec.split(",")))
        excluded = os.path.join(self.tmp.name, "excluded.txt")
        with open(excluded, "w") as fh:
            fh.write("s__Otherella s1\n")
        both = os.path.join(self.tmp.name, "both")
        run(CONVERT, "--from_db", self.db, "--genes", listed, "--exclude_species", excluded, "--outdir", both)
        rows = self.map_rows(both)
        self.assertEqual(sorted({r[1] for r in rows}), self.subset)
        self.assertEqual(len({r[0] for r in rows}), len({r[0] for r in self.map_rows(self.copy)}) - 1)
        self.same_file(both, self.copy, "gene2geneid.tsv")
        for bad in ("PF99999.1", "0", str(max(self.gene_ids.values()) + 1), "#"):
            with self.assertRaises(subprocess.CalledProcessError, msg=bad):
                run(CONVERT, "--from_db", self.db, "--genes", bad, "--outdir", os.path.join(self.tmp.name, "bad"))

    def test_neighbours_are_counted_over_the_genes_kept(self):
        subset = set(self.subset)
        kept = read_table(os.path.join(self.copy, "gene_neighbours.tsv"))
        self.assertTrue(kept)
        self.assertTrue(all(int(r["gene"]) in subset and int(r["partner"]) in subset | {0} for r in kept))
        positions = read_table(os.path.join(self.copy, "gene_positions.tsv"))
        self.assertEqual({int(r["gene"]) for r in positions}, subset)
        self.assertEqual({r["accession"] for r in positions},
                         {r["accession"] for r in read_table(os.path.join(self.db, "gene_positions.tsv")) if int(r["gene"]) in subset})
        # What the representative's genome gives among the genes kept (neighbour_ends on its placements) is in
        # the copy's table, under the species or one of its clades; a gene whose neighbour was left out faces the
        # next gene kept, farther away, where the full table had the gene left out.
        import gene_neighbours as gn
        ancestors, node = set(), self.taxid
        while node not in ancestors:
            ancestors.add(node)
            node = int(self.nodes[node][1])
        species, settings = gn.read_positions(os.path.join(self.db, "gene_positions.tsv"), genes=subset)
        max_gap = int(settings.get("max_gap", self.MAX_GAP))
        expected = gn.neighbour_ends(species[self.taxid][self.rep], max_gap)
        lines = {(int(r["clade"]), int(r["gene"]), int(r["end"]), int(r["partner"]), int(r["partner_end"])) for r in kept}
        for gene, end, partner, partner_end, gap in expected:
            self.assertTrue(any((clade, gene, end, partner, partner_end) in lines for clade in ancestors),
                            f"gene {gene} end {end} facing {partner} (end {partner_end}, {gap} bases)")
        whole, _ = gn.read_positions(os.path.join(self.db, "gene_positions.tsv"))
        before = {(g, e): (p, gap) for g, e, p, _, gap in gn.neighbour_ends(whole[self.taxid][self.rep], max_gap)}
        bridged = [(g, e, p, gap) for g, e, p, _, gap in expected if p and before[(g, e)][0] not in subset]
        self.assertTrue(bridged, "no gene end faces a gene kept where a gene left out was")
        self.assertTrue(all(gap > before[(g, e)][1] for g, e, p, gap in bridged))

    def test_protal_builds_the_subset_and_ranks_a_full_build(self):
        protal = os.environ.get("PROTAL", "")
        if not os.access(protal, os.X_OK):
            self.skipTest("$PROTAL names no protal binary")
        from gtdb_to_protal_db import full_reference_path

        def build(folder, *extra):
            return subprocess.run([protal, "--build", "--no_profile", "-t", "2", "--no_bundle", "--db", folder, "--reference",
                                   os.path.join(folder, "reference.fna"), "--full_reference", full_reference_path(folder), *extra],
                                  capture_output=True, text=True)

        def genes_of(folder, name):
            with open(os.path.join(folder, name)) as fh:
                return sorted({int(line.split("\t")[1]) for line in fh if line[0].isdigit() and line.split("\t")[1].isdigit()})

        subset = os.path.join(self.tmp.name, "subset_built")
        shutil.copytree(self.copy, subset)
        result = build(subset)
        self.assertEqual(result.returncode, 0, result.stdout[-2000:] + result.stderr[-2000:])
        self.assertEqual(genes_of(subset, "unique_kmers.tsv"), self.subset)
        if os.path.isfile(os.path.join(subset, "gene_conservation.tsv")):  # too few genes per species for factors here
            self.assertTrue(set(genes_of(subset, "gene_conservation.tsv")) <= set(self.subset))
        else:
            self.assertIn("Gene conservation: no factors", result.stdout)
        self.assertRegex(result.stdout, r"Gene neighbours: \d+ rules of \d+ clades from 1[2-6] genomes")
        # Every gene, built and ranked: a table of all genes with the domains' columns, the 3 chosen as a gene
        # list, one of them reserved for archaea (a gene of their marker set alone is rare over all species).
        every = os.path.join(self.tmp.name, "every_gene")
        run(CONVERT, "--from_db", self.db, "--outdir", every)
        self.same_file(every, self.db, "reference.map")
        result = build(every)
        self.assertEqual(result.returncode, 0, result.stdout[-2000:] + result.stderr[-2000:])
        table, listed = os.path.join(self.tmp.name, "ranking.tsv"), os.path.join(self.tmp.name, "best.txt")
        ranking = subprocess.run([sys.executable, RANK_GENES, "--db", every, "-o", table, "--top", "3", "--subset", listed],
                                 check=True, capture_output=True, text=True)
        self.assertRegex(ranking.stderr, r"best\.txt: 3 genes: \S+, \S+, \S+; scores [0-9.]+ down to [0-9.]+; in half the "
                                         r"species or more of: bacteria \d, archaea \d")
        rows = read_table(table)
        self.assertEqual(len(rows), len({r[1] for r in self.map_rows(every)}))  # the reference's genes
        self.assertEqual([int(r["rank"]) for r in rows], list(range(1, len(rows) + 1)))
        scores = [float(r["score"]) for r in rows]
        self.assertEqual(scores, sorted(scores, reverse=True))
        self.assertTrue(all(0 <= float(r[c]) <= 1 for r in rows
                            for c in ("score", "prevalence", "unique_share", "bacteria_prevalence", "archaea_prevalence")))
        self.assertTrue(all(self.gene_ids[r["marker"]] == int(r["gene_id"]) for r in rows))
        self.assertTrue(all(int(r["species"]) <= 8 and float(r["mean_length"]) > 0 for r in rows))
        self.assertGreater(scores[0], 0)
        archaeal = {int(r["gene_id"]) for r in rows if float(r["archaea_prevalence"]) >= 0.5}
        bacterial = {int(r["gene_id"]) for r in rows if float(r["bacteria_prevalence"]) >= 0.5}
        self.assertTrue(archaeal - bacterial, "genes of the archaeal marker set alone")
        self.assertTrue(all(float(r["score"]) < 0.5 for r in rows if int(r["gene_id"]) in archaeal - bacterial))
        with open(listed) as fh:
            lines = [line.rstrip("\n") for line in fh]
        best = [int(line) for line in lines if not line.startswith("#")]
        self.assertEqual(len(best), 3)
        self.assertTrue(set(best) & archaeal, "a gene archaea have")
        self.assertTrue(set(best) & bacterial, "a gene bacteria have")
        self.assertEqual(best, sorted(best, key=lambda g: next(int(r["rank"]) for r in rows if int(r["gene_id"]) == g)))
        self.assertEqual(sum(line.startswith("# rank ") for line in lines), 3)
        self.assertTrue(any("chosen for archaea" in line for line in lines))
        # Without the domains' share, the 3 best by the overall score are bacterial genes only.
        plain = os.path.join(self.tmp.name, "plain.txt")
        subprocess.run([sys.executable, RANK_GENES, "--db", every, "--top", "3", "--per-domain", "0", "--subset", plain],
                       check=True, capture_output=True, text=True)
        with open(plain) as fh:
            self.assertEqual([int(line) for line in fh if not line.startswith("#")], [int(r["gene_id"]) for r in rows[:3]])
        self.assertFalse({int(r["gene_id"]) for r in rows[:3]} & (archaeal - bacterial))
        from_list = os.path.join(self.tmp.name, "from_list")
        run(CONVERT, "--from_db", self.db, "--genes", listed, "--outdir", from_list)
        self.assertEqual(sorted({r[1] for r in self.map_rows(from_list)}), sorted(best))
        # protal --build_gene_subset: refused with a neighbours table counted over every gene; without one, the
        # database gets the listed genes' k-mers, rows and factors only, the other genes stay in reference.fna.
        with_neighbours = os.path.join(self.tmp.name, "with_neighbours")
        run(CONVERT, "--from_db", self.db, "--outdir", with_neighbours)
        result = build(with_neighbours, "--build_gene_subset", listed)
        self.assertEqual(result.returncode, 8)
        self.assertIn("Cannot build with --build_gene_subset: the folder has gene_neighbours.tsv", result.stderr)
        for name in ("gene_neighbours.tsv", "gene_positions.tsv"):
            os.remove(os.path.join(with_neighbours, name))
        result = build(with_neighbours, "--build_gene_subset", listed)
        self.assertEqual(result.returncode, 0, result.stdout[-2000:] + result.stderr[-2000:])
        self.assertEqual(genes_of(with_neighbours, "unique_kmers.tsv"), sorted(best))
        if os.path.isfile(os.path.join(with_neighbours, "gene_conservation.tsv")):
            self.assertTrue(set(genes_of(with_neighbours, "gene_conservation.tsv")) <= set(best))
        self.assertEqual(self.map_rows(with_neighbours), self.map_rows(self.db))  # every gene is still in the reference


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
                   "--test-read-pairs", "2000", "--ntree", "16", "--evaluation", "basic", "--progress-every", "5",
                   *extra]
        if not wait:
            return subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        return subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=3000)

    def text(self, *path):
        with open(os.path.join(self.tmp.name, *path)) as fh:
            return fh.read()

    def test_a_build_and_rerun(self):
        # The samples are simulated on a scratch disk of their own (--scratch), the tables copied to OUTDIR.
        scratch = ("--scratch", os.path.join(self.tmp.name, "scratch"))
        first = self.build("out", *scratch)
        self.assertEqual(first.returncode, 0, first.stdout[-3000:])
        self.assertIn("Ready protal database", first.stdout)
        for path in ("protal_db/database.protal", "model_logs/summary.txt",
                     ".stages/convert.json", ".stages/protal_db.json", ".stages/training_db.json"):
            self.assertTrue(os.path.isfile(os.path.join(self.tmp.name, "out", path)), path)
        # The training database, read only by the collections and the parity check, is built on the scratch disk.
        self.assertTrue(os.path.isfile(os.path.join(self.tmp.name, "scratch", "training_db", "database.protal")))
        self.assertFalse(os.path.exists(os.path.join(self.tmp.name, "out", "training_db")))
        # full_reference.fna, which only the builds read, is gone once they are done.
        for path in ("out/protal_db/full_reference.fna", "scratch/training_db/full_reference.fna"):
            for name in (path, path + ".zst"):
                self.assertFalse(os.path.exists(os.path.join(self.tmp.name, name)), name)
        self.assertRegex(first.stdout, r"\n\[[^]]+\]     built protal_db in the background in \d+:\d\d:\d\d.*; full_reference\.fna"
                         + (r"\.zst" if shutil.which("zstd") else "") + r" removed \([\d.]+ [MG]B\)")
        # The genes' conservation factors are in the database, and their summary in build_metadata.tsv.
        metadata = dict(line.rstrip("\n").split("\t", 1) for line in open(os.path.join(self.tmp.name, "out", "protal_db",
                                                                                       "build_metadata.tsv")))
        self.assertRegex(metadata["gene_conservation"], r"^factors [0-9.]+-[0-9.]+ for \d+ genes, from \d+ species")
        self.assertEqual(metadata["classifier_previous_procedure"], "not compared")  # without --previous-procedure
        # The gene neighbours' features and the relatives' four by the references' distance (not all of them, by
        # default), trained on samples with congener groups, the models calling at their knob curves (the other seed's build below trains the relatives features and the
        # calls at a target share of false calls).
        self.assertEqual(metadata["classifier_features"], "normalized+adjacency+distance+depth+divergence+unfiltered")
        self.assertIn("gene copies", metadata["suspect_copies"])  # the build looked for suspect copies
        self.assertIn("; congeners 0.25:2-5", metadata["classifier_training_design"])
        commands = [open(p).read() for p in glob.glob(os.path.join(self.tmp.name, "**", "run_params.tsv"), recursive=True)]
        self.assertTrue(commands, "the simulators' run_params.tsv")  # beside the samples, in --scratch
        self.assertTrue(all("--congener_groups 0.25:2-5" in c for c in commands))
        self.assertEqual(metadata["classifier_call_mode"], "curve")
        self.assertNotIn("model_pe_false_calls", metadata)
        self.assertNotIn("--fdr-calls", self.text("out", "classifier_training.log"))
        # What the conservation features rest on, on the release's genomes: how the genes differ between congeners
        # (protal --build), and where the held-out species' reads land (trace_relatives.py).
        self.assertRegex(metadata["gene_congeners"], r"^\d+ pairs of species of \d+ genera")
        # The gene neighbours of every genome to simulate from: the frequencies and the positions in the database,
        # the training database's frequencies derived anew without the species it leaves out.
        self.assertRegex(first.stdout, r"; the genes placed in \d+ genomes \(\d+ exactly, \d+ by their k-mer trace\) "
                                       r"and their neighbours counted in \d+:\d\d:\d\d[^(]* \(gene_neighbours\.log\)")
        self.assertRegex(metadata["gene_neighbours"], r"^\d+ rules of \d+ clades from \d+ genomes")
        self.assertRegex(metadata["gene_positions"], r"^\d+ genes \(\d+ placed by their k-mer trace\) in \d+ genomes of "
                                                     r"\d+ species, \d+ read as circular")
        self.assertRegex(self.text("out", "training_db.log"), r"gene_neighbours\.tsv derived from the \d+ genomes of "
                                                             r"\d+ species kept")
        for name in ("gene_congeners.tsv", "relatives_by_gene_conservation.txt"):
            self.assertTrue(os.path.isfile(os.path.join(self.tmp.name, "out", "model_logs", name)), name)
        # Each step says when it starts and, indented, when it ended, how long it took and what it made; a stage
        # running for a while (5 s here, --progress-every) how it is doing; the collector how far the simulations
        # and protal are.
        self.assertRegex(first.stdout, r"\[\d\d:\d\d:\d\d \+\d+:\d\d:\d\d\] 4/8 training data \(training_data\.log\): "
                                       r"4 pe, 4 se samples\n")
        # The models go into the database in one rewrite, each with its knob curve over depth (the trainer's
        # --depth-knobs for every read type), with up to 512 leaves per tree for short reads (--maxnodes); the
        # converter logs its steps' times.
        self.assertTrue(os.path.isfile(os.path.join(self.tmp.name, "out", "final_package.log")))
        self.assertEqual(glob.glob(os.path.join(self.tmp.name, "out", "final_package_*.log")), [])
        self.assertEqual(self.text("out", "final_package.log").count("Models for read types:"), 1)
        self.assertEqual(metadata["classifier_depth_knobs"], "pe,se")
        self.assertEqual(metadata["classifier_max_leaves"], "pe:512,se:512")
        self.assertIn("--maxnodes 512", self.text("out", "classifier_training_se.log"))
        self.assertRegex(self.text("out", "convert.log"), r"spooled the representatives' marker genes \(\d+ species\): [\d.]+ s")
        self.assertRegex(self.text("out", "convert.log"), r"joined them into full_reference\.fna(\.zst)?: [\d.]+ s")
        self.assertRegex(first.stdout, r"\n\[[^]]+\]     collected in \d+:\d\d:\d\d.*; taxa present/absent: pe \d+/\d+, "
                                       r"se \d+/\d+;")
        self.assertRegex(first.stdout, r"\n\[[^]]+\] 8/8 adding the pe, se models to protal_db")
        self.assertRegex(first.stdout, r": \d+:\d\d:\d\d so far")
        # Both collections simulate in the background from the holdout on, during the builds, and profile after.
        self.assertRegex(first.stdout, r"\n\[[^]]+\]     simulating the training data and the independent test set "
                                       r"meanwhile, in the background \(training_data_simulation\.log, "
                                       r"test_data_simulation\.log\)\n")
        self.assertLess(first.stdout.index("simulating the training data and"), first.stdout.index("built training_db in"))
        self.assertRegex(first.stdout, r"\n\[[^]]+\]     simulated the training data in the background in \d+:\d\d:\d\d")
        simulation = self.text("out", "training_data_simulation.log")
        self.assertRegex(simulation, r"2 of 2 paired-end design points, \d+:\d\d:\d\d in all")
        self.assertIn("profiling left to a run without --simulate_only", simulation)
        self.assertNotIn("protal profiled", simulation)
        collection = self.text("out", "training_data.log")
        self.assertNotRegex(collection, "simulating")
        # The test set's samples (1 pe, 1 se) are profiled in the training data's protal run: the database loads once.
        self.assertIn("and 2 samples of another collection in one protal run", collection)
        self.assertRegex(collection, r"protal profiled 10 samples in \d+:\d\d:\d\d")
        self.assertNotIn("protal profiled", self.text("out", "test_data.log"))
        self.assertIn("profiling the independent test set's samples in the same protal run", first.stdout)
        self.assertRegex(self.text("out", "test_data.log"), r"\d+ taxa in \S+training_data\.tsv: \d+ present")
        # The genome table has each genome's length, as the simulator counts it when the table has none.
        with open(os.path.join(self.tmp.name, "out", "genomes.tsv")) as fh:
            table = [line.rstrip("\n").split("\t") for line in fh]
        self.assertTrue(table and all(len(f) == 4 and f[3].isdigit() for f in table))
        three = os.path.join(self.tmp.name, "three_columns.tsv")
        with open(three, "w") as fh:
            fh.write("".join("\t".join(f[:3]) + "\n" for f in table))
        subprocess.run([os.environ["SIMULATE"], "--genome_table", three, "-o", os.path.join(self.tmp.name, "lengths"),
                        "-n", "3", "--total_read_pairs", "100", "--species_per_sample", "6", "--seed", "4", "--test"],
                       check=True, capture_output=True)
        with open(os.path.join(self.tmp.name, "lengths", "manifest.tsv")) as fh:
            counted = {row["genome"]: row["genome_length"] for row in csv.DictReader(fh, delimiter="\t")}
        self.assertTrue(counted)
        self.assertEqual(counted, {f[0]: f[3] for f in table if f[0] in counted})
        self.assertTrue(os.path.isfile(os.path.join(self.tmp.name, "out", "training", "training_data_se.tsv")))
        self.assertFalse(os.path.exists(os.path.join(self.tmp.name, "out", "training", "points")))
        self.assertTrue(os.path.isdir(os.path.join(self.tmp.name, "scratch", "training", "points")))
        self.assertRegex(first.stdout, r"The run took at most [\d.]+ [MG]B on \S+scratch")

        # A rerun converts and builds nothing, and the collector reuses its samples and dumps.
        again = self.build("out", *scratch)
        self.assertEqual(again.returncode, 0, again.stdout[-3000:])
        self.assertIn("protal_db was built by an earlier run from the same release and protal; kept", again.stdout)
        self.assertIn("training_db was built by an earlier run with the same species left out; kept", again.stdout)
        self.assertNotRegex(again.stdout, r"built (protal|training)_db")
        self.assertNotRegex(self.text("out", "training_data_simulation.log"), "simulating")
        self.assertNotRegex(self.text("out", "training_data.log"), "simulating|profiling")

        # Other species held out (another seed): only the training database is built again, from the release
        # converted anew (the finished database's build consumed the converted files), and every point is
        # simulated and profiled again rather than mixed into the table.
        # Its models with the relatives features and calls at a target share of false calls: trained, checked for
        # parity with protal, and in the database.
        other = self.build("out", "--seed", "2", "--features", "normalized+adjacency+relatives", "--call-mode", "fdr",
                           *scratch)
        self.assertEqual(other.returncode, 0, other.stdout[-3000:])
        self.assertIn("protal_db was built by an earlier run from the same release and protal; kept", other.stdout)
        self.assertRegex(other.stdout, r"built training_db in \d+:\d\d:\d\d")
        self.assertNotRegex(other.stdout, r"built protal_db")
        self.assertIn("was simulated from other inputs (or by an older collector): simulating it again",
                      self.text("out", "training_data_simulation.log"))
        self.assertFalse(os.path.exists(os.path.join(self.tmp.name, "out", ".converted")))
        metadata = dict(line.rstrip("\n").split("\t", 1) for line in open(os.path.join(self.tmp.name, "out", "protal_db",
                                                                                       "build_metadata.tsv")))
        self.assertEqual(metadata["classifier_features"], "normalized+adjacency+relatives")
        self.assertEqual(metadata["classifier_call_mode"], "fdr")
        for t in ("pe", "se"):
            self.assertRegex(metadata[f"model_{t}_false_calls"], r"^target [0-9.]+ per sample; species held out F1 [0-9.]+")
            log = self.text("out", "classifier_training" + ("" if t == "pe" else "_se") + ".log")
            self.assertIn("--fdr-calls", log)
            self.assertIn("em_own_share", log)  # among the model's features

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
        self.assertNotIn("collected in", result.stdout)
        self.assertLess(time.time() - started, 600)

    def test_c_sigterm_stops_every_command(self):
        run_ = self.build("out_term", wait=False)
        log = os.path.join(self.tmp.name, "out_term", "training_data_simulation.log")
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

    def test_d_a_reduced_database_of_the_most_distinctive_genes(self):
        # --n-genes 3: the release converted whole into .converted, the genes ranked from a full build of the
        # training database, both database folders derived for the 3 best, the models trained on them, the
        # whole conversion gone at the end; a rerun ranks, derives and builds nothing; --genes names the genes
        # instead, and both databases are built again from the samples already simulated.
        scratch = ("--scratch", os.path.join(self.tmp.name, "scratch_genes"))
        first = self.build("out_genes", "--n-genes", "3", *scratch)
        self.assertEqual(first.returncode, 0, first.stdout[-3000:])
        self.assertRegex(first.stdout, r"\n\[[^]]+\] 2/9 the release \(convert\.log\): converted whole into \.converted")
        self.assertRegex(first.stdout, r"\n\[[^]]+\] 3/9 marker genes \(gene_subset\.txt\): the 3 most distinctive by "
                                       r"prevalence x unique k-mer share\n")
        self.assertRegex(first.stdout, r"ranked the \d+ genes from a full build of the training database \(every gene, "
                                       r"\d+ species left out\), built in \d+:\d\d:\d\d[^(]*\(gene_ranking_build\.log\): "
                                       r"gene_ranking\.tsv")
        self.assertRegex(first.stdout, r"the 3 best of \d+: \S+, \S+, \S+; scores [0-9.]+ down to [0-9.]+; in half the "
                                       r"species or more of: bacteria [1-3], archaea [1-3]\n")
        self.assertRegex(first.stdout, r"protal_db's files derived for these genes in \d+:\d\d:\d\d[^(]*\(protal_db_files\.log\)")
        self.assertRegex(first.stdout, r"\n\[[^]]+\] 4/9 training database ")
        self.assertIn("Ready protal database", first.stdout)
        self.assertIn("(marker genes: 3 of ", first.stdout)
        out = os.path.join(self.tmp.name, "out_genes")
        for path in ("gene_ranking.tsv", "gene_subset.txt", "gene_ranking_build.log", "gene_ranking_files.log",
                     "protal_db_files.log", "model_logs/gene_ranking.tsv", "model_logs/gene_subset.txt",
                     "protal_db/database.protal", ".stages/gene_ranking.json", ".stages/convert.json"):
            self.assertTrue(os.path.isfile(os.path.join(out, path)), path)
        self.assertFalse(os.path.exists(os.path.join(out, ".converted")))
        self.assertFalse(os.path.exists(os.path.join(self.tmp.name, "scratch_genes", "ranking_db")))
        with open(os.path.join(out, "gene_subset.txt")) as fh:
            chosen = [int(line) for line in fh if not line.startswith("#")]
        self.assertEqual(len(chosen), 3)
        with open(os.path.join(out, "gene_ranking.tsv")) as fh:
            ranking = [line.rstrip("\n").split("\t") for line in fh]
        self.assertEqual(ranking[0][:4], ["rank", "gene_id", "marker", "score"])
        self.assertIn(int(ranking[1][1]), chosen)  # the best overall, and the best of each domain
        self.assertTrue(set(chosen) <= {int(r[1]) for r in ranking[1:]})
        with open(os.path.join(out, "protal_db", "build_metadata.tsv")) as fh:
            metadata = dict(line.rstrip("\n").split("\t", 1) for line in fh)
        self.assertRegex(metadata["marker_genes"], r"^3 of \d+, the most distinctive by prevalence x unique k-mer share "
                                                   r"\(scripts/rank_genes\.py, 1 per domain, from a full build of the training "
                                                   r"database\): \S+, \S+, \S+; in half the species or more of bacteria \d, archaea \d$")
        # One of the three is a gene archaea have (the archaeal marker set's genes score low over all species).
        header = ranking[0]
        archaeal = {int(r[1]) for r in ranking[1:] if float(r[header.index("archaea_prevalence")]) >= 0.5}
        self.assertTrue(set(chosen) & archaeal, "no gene archaea have among the three chosen")
        with open(os.path.join(out, "gene_subset.txt")) as fh:
            self.assertIn("chosen for archaea", fh.read())
        self.assertRegex(self.text("out_genes", "training_db.log"), r"derived from the \d+ genomes of \d+ species kept: "
                                                                   r"\d+ lines \(over the 3 genes kept\)")
        self.assertIn(", 3 genes kept: ", self.text("out_genes", "protal_db_files.log"))
        self.assertNotIn("genes kept", self.text("out_genes", "gene_ranking_files.log"))  # every gene, species left out
        again = self.build("out_genes", "--n-genes", "3", *scratch)
        self.assertEqual(again.returncode, 0, again.stdout[-3000:])
        self.assertIn("ranked by an earlier run from a full build of the training database; kept", again.stdout)
        self.assertIn("protal_db was built by an earlier run from these genes; kept", again.stdout)
        self.assertIn("was built by an earlier run with the same species left out and the same genes; kept", again.stdout)
        self.assertNotRegex(again.stdout, r"built (protal|training|ranking)_db|derived for these genes")
        self.assertNotRegex(self.text("out_genes", "training_data_simulation.log"), "simulating")
        # Other genes, named: the 4th and 5th of the ranking by marker id, the best by gene id.
        spec = ",".join([ranking[4][2], ranking[5][2].split(".")[0], ranking[1][1]])
        named = self.build("out_genes", "--genes", spec, *scratch)
        self.assertEqual(named.returncode, 0, named.stdout[-3000:])
        self.assertRegex(named.stdout, r"\n\[[^]]+\] 3/9 marker genes \(gene_subset\.txt\): the ones listed \(--genes\)\n")
        self.assertRegex(named.stdout, r"3 of the \d+ marker genes: \S+, \S+, \S+\n")
        self.assertRegex(named.stdout, r"built training_db in \d+:\d\d:\d\d")
        self.assertRegex(named.stdout, r"built protal_db in the background in \d+:\d\d:\d\d")
        self.assertNotIn("ranking", named.stdout)
        with open(os.path.join(out, "gene_subset.txt")) as fh:
            self.assertEqual(sorted(int(line) for line in fh if not line.startswith("#")),
                             sorted(int(ranking[i][1]) for i in (1, 4, 5)))
        with open(os.path.join(out, "protal_db", "build_metadata.tsv")) as fh:
            metadata = dict(line.rstrip("\n").split("\t", 1) for line in fh)
        self.assertRegex(metadata["marker_genes"], r"^3 of \d+ \(--genes\): \S+, \S+, \S+$")
        self.assertNotRegex(self.text("out_genes", "training_data_simulation.log"), "simulating")
        self.assertFalse(os.path.exists(os.path.join(out, ".converted")))

    def test_e_releases_full_and_reduced(self):
        # build_gtdb_releases.py: the full database with --rank-genes (the genes ranked from its training database,
        # unpacked), then the reduced one from that ranking (no ranking database built), each in its folder, with
        # the summary of both; a rerun keeps them.
        root = os.path.join(self.tmp.name, "releases")
        os.makedirs(root, exist_ok=True)
        if not os.path.exists(os.path.join(root, "gtdb_r226")):
            os.symlink(self.inputs, os.path.join(root, "gtdb_r226"))
        out = os.path.join(self.tmp.name, "dbs")
        command = [self.python, RELEASES_BUILD, "--inputs", root, "--outdir", out, "--variants", "n3,full", "--n-genes", "3",
                   "--protal", os.environ["PROTAL"], "--simulator", os.environ["SIMULATE"], "-t", "2",
                   "--scratch", os.path.join(self.tmp.name, "scratch_releases"),
                   "--samples", "2", "--read-pairs", "1000,4000", "--read-setups", "100:HS20:300:40",
                   "--species-per-sample", "6-8", "--archaea", "1", "--holdout-max-share", "0.2",
                   "--holdout-clades", "family:1,genus:1", "--read-types", "pe,se", "--test-samples", "1",
                   "--test-read-pairs", "2000", "--ntree", "16", "--evaluation", "basic"]
        result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=3000)
        self.assertEqual(result.returncode, 0, result.stdout[-3000:])
        self.assertLess(result.stdout.index("r226, full: building"), result.stdout.index("r226, n3: building"))
        full_log, reduced_log = self.text("dbs", "r226_full.log"), self.text("dbs", "r226_n3.log")
        self.assertRegex(full_log, r"ranked the \d+ genes of training_db in \d+:\d\d:\d\d[^(]*\(--rank-genes, gene_ranking\.log\): "
                                   r"gene_ranking\.tsv; the 12 a reduced database would take: ")
        self.assertTrue(os.path.isfile(os.path.join(out, "r226_full", "gene_ranking.tsv")))
        self.assertTrue(os.path.isfile(os.path.join(out, "r226_full", "model_logs", "gene_ranking.tsv")))
        self.assertFalse(os.path.exists(os.path.join(self.tmp.name, "scratch_releases", "r226_full", "ranking_files")))
        self.assertIn("ranked by " + os.path.join(out, "r226_full", "gene_ranking.tsv"), reduced_log)
        self.assertNotIn("gene_ranking_build.log", reduced_log)
        self.assertFalse(os.path.exists(os.path.join(out, "r226_n3.log.partial")))
        with open(os.path.join(out, "build_summary.tsv")) as fh:
            rows = [dict(zip(*pair)) for pair in zip([fh.readline().rstrip("\n").split("\t")] * 2,
                                                     [line.rstrip("\n").split("\t") for line in fh])]
        self.assertEqual([(r["release"], r["variant"], r["status"]) for r in rows], [("226", "full", "ok"), ("226", "n3", "ok")])
        full, reduced = rows
        self.assertEqual(full["genes"], "all")
        self.assertEqual(reduced["genes"], "3")
        self.assertIn("1 per domain", reduced["marker_genes"])
        for row in rows:
            self.assertRegex(row["protal_version"], r"^protal v[0-9.]+")
            self.assertEqual(row["species"], "60")  # the release's species, not only those simulated from
            self.assertTrue(int(row["genera"]) <= int(row["species"]) and int(row["phyla"]) >= 1, row)
            self.assertEqual(int(row["bacteria_species"]) + int(row["archaea_species"]), 60)
            self.assertGreater(int(row["archaea_species"]), 0)
            self.assertRegex(row["genomes_simulated"], r"^\d+ genomes of \d+ species$")
            self.assertRegex(row["species_held_out"], r"^\d+$")
            self.assertRegex(row["database_gb"], r"^\d+\.\d\d$")
            self.assertRegex(row["build_time"], r"^\d+:\d\d:\d\d$")
            self.assertRegex(row["build_peak_gb"], r"^\d+\.\d+$")
            self.assertRegex(row["profiling_peak_gb"], r"^\d+\.\d+$")
            for t in ("pe", "se"):
                self.assertRegex(row[f"{t}_test_F1"], r"^[0-9.]+$|^-$")
                self.assertRegex(row[f"{t}_heldout_F1"], r"^[0-9.]+$")
                self.assertRegex(row[f"{t}_heldout_FP_per_sample"], r"^[0-9.]+$")
            for t in ("pb", "ont"):
                self.assertEqual(row[f"{t}_test_F1"], "")
        self.assertLess(float(reduced["database_gb"]), float(full["database_gb"]) + 0.01)
        summary = self.text("dbs", "build_summary.txt")
        self.assertIn("GTDB r226, full: ok; protal", summary)
        self.assertRegex(summary, r"taxa: 60 species \(\d+ bacteria, \d+ archaea\), \d+ genera, \d+ families, \d+ orders, "
                                  r"\d+ classes, \d+ phyla")
        self.assertRegex(summary, r"\n  pe: independent test F1 [0-9.-]+, FP/sample [0-9.-]+, sensitivity [0-9.-]+, "
                                  r"precision [0-9.-]+; species held out F1 [0-9.]+, FP/sample [0-9.]+")
        again = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=600)
        self.assertEqual(again.returncode, 0, again.stdout[-3000:])
        self.assertEqual(again.stdout.count("is built and trained; kept"), 2)
        self.assertNotIn("building", again.stdout)


if __name__ == "__main__":
    unittest.main()
