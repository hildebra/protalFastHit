"""mini_db_fixtures.py - what the tests of the mini database's scripts and of the GTDB build scripts share.

The paths of the scripts, the synthetic GTDB release (simulate_gtdb_release.py) and its conversion, made once per process
for every test class that reads them (shared_release); a release whose markers lie in operon-like clusters, with its
gene positions (operon_release); the stand-ins for GTDB's and NCBI's servers and NCBI's `datasets`; and checks and
readers of the converter's files. Imported by test_mini_db.py, test_downloads.py, test_collector.py,
test_gene_neighbours.py and test_gtdb_build.py, which run from the repository root as

  python3 -m unittest scripts/mini_db/test_mini_db.py
"""

import atexit
import functools
import gzip
import hashlib
import http.server
import os
import subprocess
import sys
import tarfile
import tempfile
import threading

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.dirname(HERE)
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)
import prerequisites  # noqa: E402,F401  (scripts/prerequisites.py, for the test files)

SIMULATE = os.path.join(HERE, "simulate_gtdb_release.py")
CONVERT = os.path.join(HERE, "gtdb_to_protal_db.py")
GENE_NEIGHBOURS = os.path.join(HERE, "gene_neighbours.py")
SIMULATE_READS = os.path.join(HERE, "simulate_reads.py")
LINEAGES = os.path.join(HERE, "gtdb_like_lineages.py")
DOWNLOAD = os.path.join(SCRIPTS, "download_gtdb.py")
RANK_GENES = os.path.join(SCRIPTS, "rank_genes.py")
RELEASES_DOWNLOAD = os.path.join(SCRIPTS, "download_gtdb_releases.py")
RELEASES_BUILD = os.path.join(SCRIPTS, "build_gtdb_releases.py")
BUILD = os.path.join(SCRIPTS, "build_gtdb_database.py")

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
    # Withdrawn genomes: in the package's report, not in its fetch.txt, which a package of only such has none of.
    have = [a for a in accessions if a not in os.environ.get("FAKE_WITHDRAWN", "").split(",")]
    with zipfile.ZipFile(args[args.index("--filename") + 1], "w") as z:
        z.writestr("ncbi_dataset/data/assembly_data_report.jsonl", "")
        if have:
            z.writestr("ncbi_dataset/fetch.txt", "\n".join(have))
elif args[0] == "rehydrate":
    folder = args[args.index("--directory") + 1]
    paths = {r["accession"]: r["fasta_path"] for r in csv.DictReader(open(os.environ["FAKE_TABLE"]), delimiter="\t")}
    for a in open(os.path.join(folder, "ncbi_dataset", "fetch.txt")).read().split():
        os.makedirs(os.path.join(folder, "ncbi_dataset", "data", a), exist_ok=True)
        with gzip.open(paths[a]) as fin, open(os.path.join(folder, "ncbi_dataset", "data", a, a + "_ASM1v1_genomic.fna"), "wb") as fout:
            shutil.copyfileobj(fin, fout)
'''


def run(*args):
    """A script of this repository with this Python; CalledProcessError if it fails."""
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


def read_table(path):
    """A tab-separated table with a header line (# lines skipped) -> [row dict]."""
    with open(path) as fh:
        rows = [line.rstrip("\n").split("\t") for line in fh if not line.startswith("#")]
    return [dict(zip(rows[0], r)) for r in rows[1:]]


def full_reference(folder):
    """The content of a converted folder's full reference (full_reference.fna.zst, or .fna without zstd)."""
    from gtdb_to_protal_db import full_reference_path, read_full_reference
    with read_full_reference(full_reference_path(folder)) as fh:
        return fh.read()


def kmers(seq, k=21):
    return {seq[i:i + k] for i in range(len(seq) - k + 1)}


def map_rows(folder):
    """reference.map's rows: (taxid, gene id, start, end)."""
    with open(os.path.join(folder, "reference.map")) as fh:
        return [tuple(map(int, line.split("\t"))) for line in fh]


def check_reference_map(test, folder):
    """reference.map points at each sequence line of reference.fna: the header ">taxid_geneid" just before its start,
    a newline at its end, bases only between, one row per record. -> the rows."""
    with open(os.path.join(folder, "reference.fna"), "rb") as fh:
        data = fh.read()
    rows = map_rows(folder)
    test.assertEqual(len(rows), data.count(b">"))
    for taxid, geneid, start, end in rows:
        test.assertEqual(data[data.rindex(b">", 0, start):start], f">{taxid}_{geneid}\n".encode())
        test.assertEqual(data[end:end + 1], b"\n")
        test.assertRegex(data[start:end].decode(), r"^[ACGT]+$")
    return rows


def same_files(test, a, b, *names):
    for name in names:
        with open(os.path.join(a, name), "rb") as x, open(os.path.join(b, name), "rb") as y:
            test.assertEqual(x.read(), y.read(), name)


def scratch_dir(prefix):
    """A temporary folder removed when the process ends."""
    tmp = tempfile.TemporaryDirectory(prefix=prefix)
    atexit.register(tmp.cleanup)
    return tmp.name


@functools.lru_cache(maxsize=None)
def shared_release():
    """The default synthetic release (simulate_gtdb_release.py --genome_length 20000: Mockella alpha and beta, Fakibacter
    gamma, three genomes each) and its conversion, made once per process: -> (gtdb, db). Read only: a test that changes
    files copies them first."""
    root = scratch_dir("mini_db_release_")
    gtdb, db = os.path.join(root, "gtdb"), os.path.join(root, "db")
    run(SIMULATE, "--outdir", gtdb, "--genome_length", "20000")
    run(CONVERT, "--gtdb", gtdb, "--outdir", db)
    return gtdb, db


@functools.lru_cache(maxsize=None)
def operon_release():
    """A release whose markers lie in operon-like clusters (--operons, half of the clusters broken in a family), 4
    species of the family Simulaceae, 2 of Otheraceae and 2 archaea, two genomes each in 4 contigs of 300 kb in all;
    converted, and gene_neighbours.py run on it (gene_neighbours.tsv and gene_positions.tsv in the folder), made once per
    process: -> (gtdb, db, gene_neighbours.py's console output). Read only."""
    root = scratch_dir("mini_db_operons_")
    lineages = os.path.join(root, "lineages.txt")
    order = "d__Bacteria;p__Simulatota;c__Simulatia;o__Simulales"
    with open(lineages, "w") as fh:
        fh.writelines(f"{order};f__Simulaceae;g__Mockella;s__Mockella s{i}\n" for i in range(4))
        fh.writelines(f"{order};f__Otheraceae;g__Otherella;s__Otherella s{i}\n" for i in range(2))
        fh.writelines(f"d__Archaea;p__Archota;c__Archia;o__Archales;f__Archaceae;g__Archella;s__Archella s{i}\n"
                      for i in range(2))
    gtdb, db = os.path.join(root, "gtdb"), os.path.join(root, "db")
    run(SIMULATE, "--outdir", gtdb, "--lineages", lineages, "--operons", "--operon_breaks", "0.5",
        "--genome_length", "300000", "--genomes_per_species", "2", "--contigs", "4")
    run(CONVERT, "--gtdb", gtdb, "--outdir", db)
    output = subprocess.run([sys.executable, GENE_NEIGHBOURS, "--db", db, "--genome_table",
                             os.path.join(gtdb, "simulation", "genomes.tsv"), "-t", "2"],
                            check=True, capture_output=True, text=True).stdout
    return gtdb, db, output


# ---- stand-ins for GTDB's and NCBI's servers -------------------------------------------------------------------

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


class QuietHandler(http.server.SimpleHTTPRequestHandler):
    def log_message(self, *args):
        pass


def serve(directory):
    """A local HTTP server of a folder's files (GTDB's mirror, NCBI's FTP server) -> (its URL, a function that stops
    it)."""
    server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), functools.partial(QuietHandler, directory=directory))
    threading.Thread(target=server.serve_forever, daemon=True).start()

    def stop():
        server.shutdown()
        server.server_close()
    return f"http://127.0.0.1:{server.server_port}", stop


def fake_datasets(path):
    """The `datasets` stand-in (FAKE_DATASETS) as an executable at path -> path."""
    with open(path, "w") as fh:
        fh.write(FAKE_DATASETS)
    os.chmod(path, 0o755)
    return path


class FakeNcbi:
    """A synthetic release served as GTDB's mirror is (gtdb_mirror) and as NCBI's FTP server (URL + "/ftp", which
    holds nothing unless a test puts genomes there), and the `datasets` stand-in, in folder: url, ftp_url, datasets,
    env (FAKE_TABLE for the stand-in) and download_options (--mirror, --ftp_url and --datasets for download_gtdb.py:
    nothing goes to NCBI). stop() ends the server; a TestCase passes itself to have it stopped at its cleanup."""

    def __init__(self, release, folder, test=None, assembly_names=False):
        self.mirror = os.path.join(folder, "mirror")
        gtdb_mirror(release, self.mirror, assembly_names)
        self.url, self.stop = serve(self.mirror)
        if test is not None:
            test.addCleanup(self.stop)
        self.ftp_url = self.url + "/ftp"
        self.datasets = fake_datasets(os.path.join(folder, "datasets"))
        self.env = dict(os.environ, FAKE_TABLE=os.path.join(release, "simulation", "genomes.tsv"))
        self.download_options = ["--mirror", self.url, "--ftp_url", self.ftp_url, "--datasets", self.datasets]
