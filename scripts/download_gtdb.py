#!/usr/bin/env python3
"""Download what building a protal database and training its model need, on a node with internet.

Run it once per GTDB release on a node that can reach the internet; the folder it writes serves
every later build of that release (build_gtdb_database.py --inputs), and a rerun downloads only
what is missing or does not match its checksum.

    python3 scripts/download_gtdb.py -o /shared/protal_inputs/gtdb_r226             # GTDB r226
    python3 scripts/download_gtdb.py -o /shared/protal_inputs/gtdb_r220 --release 220
    python3 scripts/build_gtdb_database.py --inputs /shared/protal_inputs/gtdb_r226 --outdir OUT ...

Two parts:
- GTDB's files (data.gtdb.ecogenomic.org, --mirror): taxonomy, metadata, and the marker genes of
  the representatives and of all genomes, checked against the release's MD5SUM.txt and extracted.
  Releases 207 and later (bac120 and ar53 marker sets); the newest point release (e.g. 214.1) unless
  --release names one. With --rep_genomes gtdb also GTDB's archive of all representative genomes
  (127 GB for r226). Unless --no_msa also GTDB-Tk's alignments of the representatives' marker
  proteins (kept packed) and the trees, for the column weights (2.2 GB for r226).
- Genomes to simulate training samples from, from NCBI (the `datasets` CLI): GTDB distributes the
  representatives' genomes only, which are the database's own references; real samples hold other
  strains. From the metadata, per domain in GTDB's proportions: --species species with
  non-representative genomes that pass the quality filters (CheckM2 completeness and
  contamination), up to --per_species of those each, drawn at random; and --rep_only_species more
  species simulated from their representative only. With --rep_genomes ncbi (the default) the
  representatives of these species come from NCBI too, so the 127 GB archive is not needed. The
  simulator draws a species uniformly, then one of its genomes: a species with k strains is
  simulated from a strain k/(k+1) of the time. The species are taken in a random order fixed by the
  seed, and a rerun keeps what the folder has: a larger --species adds strain species (the earlier
  representative-only ones first) and downloads only their strains.

Which strains of a species (those passing the CheckM2 filters) are taken is not left to chance: the
best by, in this order, (1) isolate genomes before single-cell and metagenome-assembled ones (GTDB's
ncbi_genome_category, and an "uncultured" or "metagenome" organism name), (2) assembly level (complete
genome, chromosome, scaffold, contig), (3) a PacBio or Nanopore assembly (the sequencing technology at
NCBI, asked for the --tech_candidates best candidates of each species), (4) fewer contigs. Genomes
that tie are drawn at random, so a species with thousands of complete genomes is not always
represented by the same two. A species whose genomes are all MAGs keeps its best MAGs. --progenomes
limits the strains to the genomes of a proGenomes table (isolate genomes that passed CheckM2 and GUNC).
The representatives are GTDB's and stay as they are, as the database's references.

The genomes are fetched from NCBI's FTP server (--ftp_url) directly, --connections at a time: the
URL of a genome is built from its accession and the assembly name in GTDB's metadata, the file is
kept as NCBI compresses it, and is checked against its announced length and read through gzip. What
that does not deliver (an assembly NCBI renamed or withdrew, a server that fails) goes through the
`datasets` CLI, batch by batch, which is slower: it asks NCBI per file, one batch after the other.

Written to OUT/:
  release/                GTDB's files as build_gtdb_database.py --gtdb reads them
  genomes/                <accession>.fna.gz of the NCBI genomes
  genomes.tsv             accession, species, role (representative or strain), lineage, CheckM2
                          completeness and contamination, genome category (isolate, sag or mag),
                          assembly level, contig count and sequencing technology at NCBI (strains
                          only), of the genomes in genomes/
  ncbi_info.tsv           the sequencing technology NCBI gives for the candidate strains (kept for reruns)
  simulation_species.txt  the species to simulate from (build_gtdb_database.py --simulate-species)
  missing.txt             accessions NCBI did not deliver
  host/                   <accession>.fna.gz, the host genome of the scenarios with host reads
                          (build_gtdb_database.py --scenarios host): by default the human genome
                          (T2T-CHM13v2.0, 0.9 GB), kept gzipped as NCBI serves it (3.1 GB unpacked)
  download.json           the release, the options, the counts and the checksums of what is there
"""

import argparse
import collections
import concurrent.futures
import datetime
import gzip
import hashlib
import http.client
import json
import os
import random
import re
import shutil
import subprocess
import sys
import threading
import time
import urllib.error
import urllib.request
import zlib

MIRROR = "https://data.gtdb.ecogenomic.org/releases"
FTP = "https://ftp.ncbi.nlm.nih.gov/genomes/all"
MARKER_SETS = ("bac120", "ar53")
USER_AGENT = "protal-download_gtdb (https://protal.earlham.ac.uk)"
ATTEMPTS = 5  # tries of a genome before it is given up on
RETRY_WAIT = 2.0  # seconds before the second try of a genome, doubled for each further try (or the server's Retry-After)
FAILS_IN_A_ROW = 25  # direct fetches that fail in a row (not counting missing files) stop the direct fetching
# Host genomes by name (--host_genome): NCBI accession and assembly name. The human one is the complete T2T assembly
# (every centromere and the rDNA arrays, no alternative haplotypes, chrY of HG002), whose reads are what a host-
# dominated sample holds, rather than GRCh38's, whose gaps and alternative loci are not.
HOST_GENOMES = {"human": ("GCF_009914755.1", "T2T-CHM13v2.0")}


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("-o", "--out", required=True, help="folder for this release's inputs, reused by later builds")
    p.add_argument("--release", default="226",
                   help="GTDB release: 226 (default), another from 207 on (e.g. 220), or a point release (214.1)")
    p.add_argument("--mirror", default=MIRROR, help=f"base URL of GTDB's releases (default {MIRROR})")
    p.add_argument("--rep_genomes", choices=["ncbi", "gtdb"], default="ncbi",
                   help="representatives' genomes of the simulated species from NCBI (default), or GTDB's archive "
                        "of all representative genomes (127 GB for r226)")
    p.add_argument("--no_genomes", action="store_true", help="GTDB's files only, no genomes from NCBI")
    p.add_argument("--no_msa", action="store_true",
                   help="leave out GTDB's alignments of the representatives' marker proteins and its trees (2.2 GB at "
                        "r226), from which protal --build takes the column weights' columns and ancestral sequences; "
                        "without them it aligns the proteins itself")
    p.add_argument("--dry_run", action="store_true", help="list GTDB's files to download, with their sizes, and stop")
    p.add_argument("--keep_archives", action="store_true",
                   help="keep the downloaded .tar.gz archives after extracting them (default: removed)")
    p.add_argument("--species", type=int, default=16000,
                   help="species to download strains of (default 16000; 6000 before 2026-10-05, when the r226 v12 build "
                        "had 6000 such species against 18978 simulated from their representative and an in-silico strain, "
                        "and its models missed real strains 2.4 times as often as in-silico ones: "
                        "docs/claude/2026-10-05-r226-v12-scenarios). A rerun with more keeps the species and strains it "
                        "has, and gives strains to its representative-only species first")
    p.add_argument("--per_species", type=int, default=2, help="strains per species at most (default 2)")
    p.add_argument("--rep_only_species", type=int, default=9000,
                   help="further species, simulated from their representative only (and an in-silico strain, "
                        "build_gtdb_database.py --insilico-strains) (default 9000: with --species, 25,000 species, what "
                        "the soil scenarios need at full size; 2000 before 2026-10-05)")
    p.add_argument("--min_completeness", type=float, default=90.0, help="CheckM2 completeness, %% (default 90)")
    p.add_argument("--max_contamination", type=float, default=5.0, help="CheckM2 contamination, %% (default 5)")
    p.add_argument("--tech_candidates", type=int, default=30,
                   help="candidate strains per species whose sequencing technology is asked of NCBI, the best by "
                        "category and assembly level (default 30); PacBio and Nanopore assemblies are preferred "
                        "among them")
    p.add_argument("--no_tech_lookup", action="store_true",
                   help="do not ask NCBI for sequencing technologies: no preference for long-read assemblies")
    p.add_argument("--progenomes",
                   help="a proGenomes ANI-clustering table (pg4_ANI_clustering.tsv.gz of progenomes.embl.de/download.cgi: "
                        "a cluster and the GenBank accessions of its genomes), or any file listing accessions: "
                        "only strains it lists are taken")
    p.add_argument("--host_genome", default="human",
                   help="the host genome of the scenarios with host reads (build_gtdb_database.py --scenarios host), "
                        "into OUT/host, gzipped: human (default: T2T-CHM13v2.0, GCF_009914755.1, 0.9 GB), "
                        "ACCESSION_ASSEMBLYNAME of another NCBI assembly (e.g. GCF_000001405.40_GRCh38.p14), or none")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--connections", type=int, default=8,
                   help="genomes fetched at a time straight from NCBI's FTP server (default 8: NCBI answered HTTP 503 "
                        "to 32 at a time); 0: only through the datasets CLI, a batch at a time")
    p.add_argument("--ftp_url", default=FTP, help=f"NCBI's genome archive (default {FTP})")
    p.add_argument("--datasets", default="datasets", help="NCBI datasets binary (default: on PATH)")
    p.add_argument("--batch", type=int, default=500, help="genomes per NCBI request (default 500)")
    p.add_argument("-t", "--threads", type=int, default=8, help="parallel downloads and compression (default 8)")
    return p.parse_args(argv)


# ---- GTDB's files -----------------------------------------------------------------------------------------

def fetch_text(url):
    try:
        with urllib.request.urlopen(url, timeout=120) as response:
            return response.read().decode()
    except OSError as e:
        sys.exit(f"cannot read {url}: {e}")


def resolve_release(opts):
    """(release number as in file names, e.g. '226'; version folder, e.g. '226.0')."""
    text = opts.release.lower().removeprefix("r").removeprefix("release")
    number = text.split(".")[0]
    if not number.isdigit():
        sys.exit(f"--release {opts.release}: expected a number such as 226 or 214.1")
    if int(number) < 207:
        sys.exit(f"GTDB r{number} predates the ar53 marker set (it has ar122); protal databases need release 207 or later")
    if "." in text:
        return number, text
    listing = fetch_text(f"{opts.mirror}/release{number}/")
    versions = sorted({v for v in re.findall(rf'href="({number}\.\d+)/?"', listing)}, key=lambda v: int(v.split(".")[1]))
    if not versions:
        sys.exit(f"no GTDB release {number} at {opts.mirror}/release{number}/")
    return number, versions[-1]


def release_files(md5sums, number, rep_genomes, msa=True):
    """The files of the release to download: [(path in the release, md5, extract?)]."""
    wanted = []

    def one_of(*names):
        for name in names:
            if name in md5sums:
                return name
        return None

    for mset in MARKER_SETS:
        taxonomy = one_of(f"{mset}_taxonomy_r{number}.tsv.gz", f"{mset}_taxonomy_r{number}.tsv")
        metadata = one_of(f"{mset}_metadata_r{number}.tsv.gz", f"{mset}_metadata_r{number}.tar.gz")
        reps = one_of(f"genomic_files_reps/{mset}_marker_genes_reps_r{number}.tar.gz")
        every = one_of(f"genomic_files_all/{mset}_marker_genes_all_r{number}.tar.gz")
        missing = [what for what, name in (("taxonomy", taxonomy), ("metadata", metadata),
                                           ("representatives' marker genes", reps), ("all marker genes", every)) if not name]
        if missing:
            sys.exit(f"GTDB r{number} has no {mset} " + ", ".join(missing) + (
                " (releases before 207 use the ar122 marker set; protal reads bac120 and ar53)" if mset == "ar53" else ""))
        wanted += [(taxonomy, False), (metadata, metadata.endswith(".tar.gz")), (reps, True), (every, True)]
    if rep_genomes == "gtdb":
        genomes = one_of(f"genomic_files_reps/gtdb_genomes_reps_r{number}.tar.gz")
        if not genomes:
            sys.exit(f"GTDB r{number} has no genomic_files_reps/gtdb_genomes_reps_r{number}.tar.gz")
        wanted.append((genomes, True))
    # The species clusters (each representative's ANI circumscription radius, its cluster's intra-species ANI and
    # size): optional, for the converter's species_priors.tsv.
    clusters = one_of(f"auxillary_files/sp_clusters_r{number}.tsv", f"auxillary_files/sp_clusters_r{number}.tsv.gz")
    if clusters:
        wanted.append((clusters, False))
    # GTDB-Tk's alignments of the representatives' marker proteins and the trees inferred from them (since 2026-10-10,
    # unless --no_msa): the converter writes them into the database folder for protal --build's column weights,
    # which then take the family columns from GTDB's alignment and each genus's and family's ancestral sequence from
    # the tree (2.2 GB at r226, read packed).
    if msa:
        for mset in MARKER_SETS:
            for name in (f"genomic_files_reps/{mset}_msa_marker_genes_reps_r{number}.tar.gz", f"{mset}_r{number}.tree"):
                if name in md5sums:
                    wanted.append((name, False))
    if "VERSION.txt" in md5sums:
        wanted.append(("VERSION.txt", False))
    return [(name, md5sums[name], extract) for name, extract in wanted]


def md5_of(path):
    digest = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 22), b""):
            digest.update(block)
    return digest.hexdigest()


def expected_size(response, start):
    """The whole file's size by a response to a request from byte start: from Content-Range (bytes a-b/total)
    or Content-Length; None if it does not say."""
    total = re.match(r"bytes \d+-\d+/(\d+)", response.headers.get("Content-Range", ""))
    if total:
        return int(total.group(1))
    length = response.headers.get("Content-Length")
    return start + int(length) if length and length.isdigit() else None


def download(url, dest, md5, attempts=3):
    """dest from url, resuming a partial download (dest.part, of this run or an earlier one), until its MD5 is
    md5. A connection that drops returns short data without an error: the next request goes on from there.
    `attempts` requests in a row that bring nothing new give up; a .part that turns out complete (the server
    has nothing after its end, 416) is checked and kept."""
    part, name = dest + ".part", os.path.basename(dest)
    failed = 0
    while failed < attempts:
        start = os.path.getsize(part) if os.path.exists(part) else 0
        request = urllib.request.Request(url, headers={"Range": f"bytes={start}-"} if start else {})
        size = None
        try:
            with urllib.request.urlopen(request, timeout=300) as response:
                resumed = start and response.status == 206
                size = expected_size(response, start if resumed else 0)
                with open(part, "ab" if resumed else "wb") as fh:
                    shutil.copyfileobj(response, fh, 1 << 22)
        except urllib.error.HTTPError as e:
            if e.code != 416 or not start:
                failed += 1
                print(f"  {name}: {e} (attempt {failed} of {attempts})", flush=True)
                continue
            size = start  # nothing after the .part's end: it is whole, or longer than the file
        except (OSError, http.client.HTTPException) as e:
            failed += 1 if (os.path.getsize(part) if os.path.exists(part) else 0) == start else 0
            print(f"  {name}: {e}; going on from byte {os.path.getsize(part) if os.path.exists(part) else 0}", flush=True)
            continue
        have = os.path.getsize(part)
        if size is not None and have < size:
            failed += 1 if have == start else 0
            print(f"  {name}: {have} of {size} bytes (the connection dropped); going on", flush=True)
            continue
        if md5_of(part) == md5:
            os.replace(part, dest)
            return
        failed += 1
        print(f"  {name}: checksum mismatch, downloading again (attempt {failed} of {attempts})", flush=True)
        os.remove(part)
    sys.exit(f"could not download {url} with MD5 {md5}" +
             (f"; {part} holds what came, which a rerun goes on from" if os.path.exists(part) else ""))


def extracted_path(archive):
    """What the converter reads from an archive: bac120_marker_genes_reps_r226.tar.gz -> the folder
    bac120_marker_genes_reps_r226, bac120_metadata_r207.tar.gz -> bac120_metadata_r207.tsv."""
    base = archive[:-len(".tar.gz")]
    return base + ".tsv" if "_metadata_" in os.path.basename(base) else base


def extract(archive, into):
    command = ["tar", "-xf", archive, "-C", into]
    command[1:1] = ["-I", "pigz"] if shutil.which("pigz") else ["-z"]
    if subprocess.run(command).returncode != 0:
        sys.exit(f"extracting {archive} failed")
    if not os.path.exists(extracted_path(archive)):
        sys.exit(f"{archive} did not extract to {extracted_path(archive)}, which the converter reads")


def get_release(opts, state):
    number, version = resolve_release(opts)
    base = f"{opts.mirror}/release{number}/{version}"
    release_dir = os.path.join(opts.out, "release")
    os.makedirs(release_dir, exist_ok=True)
    md5sums = {}
    for line in fetch_text(f"{base}/MD5SUM.txt").splitlines():
        parts = line.split()
        if len(parts) == 2:
            md5sums[parts[1].removeprefix("./")] = parts[0]
    files = state.setdefault("files", {})
    print(f"GTDB r{number} ({version}) from {base}", flush=True)
    if opts.dry_run:
        total = 0
        for name, _md5, _unpack in release_files(md5sums, number, opts.rep_genomes, not opts.no_msa):
            with urllib.request.urlopen(urllib.request.Request(f"{base}/{name}", method="HEAD"), timeout=120) as r:
                size = int(r.headers.get("Content-Length", 0))
            total += size
            print(f"  {name}: {size / 1e9:.2f} GB", flush=True)
        print(f"  total {total / 1e9:.1f} GB from GTDB; genomes from NCBI come on top (about 4 MB each)")
        sys.exit(0)
    for name, md5, unpack in release_files(md5sums, number, opts.rep_genomes, not opts.no_msa):
        dest = os.path.join(release_dir, name)
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        done = files.get(name, {})
        there = os.path.exists(extracted_path(dest)) if unpack else os.path.exists(dest)
        if done.get("md5") == md5 and there:
            print(f"  {name}: already there", flush=True)
            continue
        if not (os.path.exists(dest) and md5_of(dest) == md5):
            print(f"  {name}: downloading", flush=True)
            download(f"{base}/{name}", dest, md5)
        entry = {"md5": md5, "bytes": os.path.getsize(dest)}
        if unpack:
            print(f"  {name}: extracting", flush=True)
            extract(dest, os.path.dirname(dest))
            entry["extracted"] = True
            if not opts.keep_archives:
                os.remove(dest)
        files[name] = entry
        save_state(opts, state)
    state["release"] = {"number": number, "version": version, "url": base}
    return number


# ---- genomes from NCBI ------------------------------------------------------------------------------------

Genome = collections.namedtuple(
    "Genome", "accession genbank species lineage is_rep completeness contamination category level contigs name",
    defaults=("",))
# category: isolate, sag (single-cell amplified) or mag (metagenome-assembled); level: NCBI's assembly level
# in lower case ("" if unknown); contigs: the number of contigs (None if unknown); genbank: the GenBank
# accession (GCA_) of the assembly, which proGenomes lists; name: the assembly name ("" if unknown), which
# is part of the genome's address on NCBI's FTP server.

LEVELS = ("complete genome", "chromosome", "scaffold", "contig")
CATEGORIES = ("isolate", "sag", "mag")
LONG_READS = re.compile(r"pacbio|pacific bio|nanopore|oxford|minion|gridion|promethion|flongle|sequel|revio|smrt|hifi",
                        re.I)


def genome_category(category, organism):
    """isolate, sag or mag, from GTDB's ncbi_genome_category ("none", or what the genome is derived from:
    "derived from metagenome", "derived from single cell", ...) and, where that says none, an organism
    name that marks an uncultured or metagenomic genome."""
    text = (category or "").strip().lower()
    if "single cell" in text:
        return "sag"
    if text not in ("", "none", "na", "n/a"):
        return "mag"
    return "mag" if re.search(r"uncultured|metagenom", organism or "", re.I) else "isolate"


def read_metadata(gtdb, release):
    """Every genome as a Genome."""
    genomes = []
    for mset in MARKER_SETS:
        path = next((os.path.join(gtdb, f"{mset}_metadata_r{release}{ext}") for ext in (".tsv", ".tsv.gz")
                     if os.path.exists(os.path.join(gtdb, f"{mset}_metadata_r{release}{ext}"))), None)
        if path is None:
            continue
        with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            col = {name: i for i, name in enumerate(header)}
            completeness = col.get("checkm2_completeness", col.get("checkm_completeness"))
            contamination = col.get("checkm2_contamination", col.get("checkm_contamination"))
            for name in ("accession", "gtdb_taxonomy", "gtdb_representative"):
                if name not in col:
                    sys.exit(f"{path} has no {name} column")

            def field(f, *names):  # the first of these columns that has a value
                for name in names:
                    if name in col and col[name] < len(f) and f[col[name]] not in ("", "none", "None", "na", "NA"):
                        return f[col[name]]
                return ""

            for line in fh:
                f = line.rstrip("\n").split("\t")
                lineage = f[col["gtdb_taxonomy"]]
                accession = f[col["accession"]]
                accession = accession[3:] if accession[:3] in ("RS_", "GB_") else accession
                genbank = field(f, "ncbi_genbank_assembly_accession") or accession.replace("GCF_", "GCA_", 1)
                contigs = field(f, "ncbi_contig_count", "contig_count")
                genomes.append(Genome(
                    accession, genbank, lineage.split(";")[-1], lineage, f[col["gtdb_representative"]] == "t",
                    float(f[completeness]) if completeness is not None and f[completeness] not in ("", "none") else 100.0,
                    float(f[contamination]) if contamination is not None and f[contamination] not in ("", "none") else 0.0,
                    genome_category(field(f, "ncbi_genome_category"), field(f, "ncbi_organism_name")),
                    field(f, "ncbi_assembly_level").lower(),
                    int(float(contigs)) if contigs else None,
                    field(f, "ncbi_assembly_name")))
    if not genomes:
        sys.exit(f"no bac120/ar53_metadata_r{release}.tsv[.gz] in {gtdb}")
    return genomes


def fragmentation(contigs):
    """0 (a few contigs) to 4 (many), 5 if unknown."""
    if contigs is None:
        return 5
    return next((i for i, limit in enumerate((5, 20, 100, 500)) if contigs <= limit), 4)


def quality_rank(genome, long_read=False):
    """What makes a strain better for simulating samples from, as a tuple that is smaller for a better genome:
    isolate before single-cell before metagenome-assembled, assembly level, a long-read assembly, fewer
    contigs."""
    level = LEVELS.index(genome.level) if genome.level in LEVELS else len(LEVELS)
    return (CATEGORIES.index(genome.category), level, 0 if long_read else 1, fragmentation(genome.contigs))


def read_accessions(path):
    """The GenBank and RefSeq accessions (GCA_/GCF_) in a file, plain or gzipped."""
    with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as fh:
        return set(re.findall(r"GC[AF]_\d{9}\.\d+", fh.read()))


def previous_choice(out):
    """What an earlier run into the folder chose: {"strains": {species: [strain accessions delivered]}, "rep_only":
    {species simulated from their representative only}}, from its genomes.tsv and simulation_species.txt (empty for a
    new folder)."""
    strains, pool = collections.defaultdict(list), set()
    path = os.path.join(out, "genomes.tsv")
    if os.path.isfile(path):
        with open(path) as fh:
            header = next(fh, "").rstrip("\n").split("\t")
            if "species" in header and "role" in header:
                species, role, accession = header.index("species"), header.index("role"), header.index("accession")
                for fields in (line.rstrip("\n").split("\t") for line in fh):
                    if len(fields) > role and fields[role] == "strain":
                        strains[fields[species]].append(fields[accession])
    path = os.path.join(out, "simulation_species.txt")
    if os.path.isfile(path):
        with open(path) as fh:
            pool = {line.strip() for line in fh if line.strip()}
    return {"strains": dict(strains), "rep_only": pool - set(strains)}


def pick(genomes, opts, technology=None, allowed=None, previous=None):
    """(strains [(Genome, sequencing technology)], species with strains, species simulated from their
    representative only, {species: domain}, {species: representative Genome}, {species: candidate strains}).
    technology(accessions) -> {accession: sequencing technology string} is asked for the best candidates of the
    species that get strains; allowed: accessions, if only these genomes may be strains.

    Species are taken per domain in one random order of the domain's species (the seed's): the first --species of
    those with candidate strains, then --rep_only_species of the others. previous (previous_choice): what an earlier
    run into the folder chose, kept as far as the counts allow: its strain species first, then, for more strain
    species, its representative-only ones (whose representatives are there), and its strains of a species before the
    species' other candidates. So a rerun with a larger --species downloads only the new strains, instead of a new
    draw of every species and genome."""
    previous = previous or {"strains": {}, "rep_only": set()}
    candidates = collections.defaultdict(list)
    domain, lineage, representative = {}, {}, {}
    for g in genomes:
        domain[g.species], lineage[g.species] = g.lineage.split(";")[0], g.lineage
        if g.is_rep:
            representative[g.species] = g
        elif g.completeness >= opts.min_completeness and g.contamination <= opts.max_contamination \
                and (allowed is None or g.genbank in allowed or g.accession in allowed):
            candidates[g.species].append(g)
    with_strains = collections.defaultdict(list)
    without = collections.defaultdict(list)
    for species in sorted(domain):
        (with_strains if candidates.get(species) else without)[domain[species]].append(species)
    earlier_strains, earlier_rep_only = set(previous["strains"]), set(previous["rep_only"])
    chosen, rep_only = [], []
    for d in sorted(set(domain.values())):
        share = sum(1 for s in domain if domain[s] == d) / len(domain)
        order = with_strains[d] + without[d]
        order.sort()
        random.Random(f"{opts.seed}:{d}").shuffle(order)  # a larger count takes the species after
        eligible = [s for s in order if candidates.get(s)]
        # The earlier strain species, then the earlier representative-only ones, then the rest (a stable sort keeps
        # the random order within each).
        eligible.sort(key=lambda s: 0 if s in earlier_strains else 1 if s in earlier_rep_only else 2)
        picked_here = eligible[:min(len(eligible), round(opts.species * share))]
        chosen += sorted(picked_here)
        taken = set(picked_here)
        rest = [s for s in order if s not in taken]
        rest.sort(key=lambda s: 0 if s in earlier_rep_only or s in earlier_strains else 1)
        rep_only += rest[:min(len(rest), round(opts.rep_only_species * share))]
    # Genomes that tie are drawn at random: shuffled once (by the species' own seed), then sorted by rank, which keeps
    # their order. A species' shortlist is its best candidates by what the metadata says (category, level); NCBI is
    # asked for their sequencing technology in one go, then the shortlist is ranked again with it. An earlier run's
    # strains of the species come first, while they are still candidates.
    shortlists, kept = {}, {}
    for species in sorted(chosen):
        options = sorted(candidates[species], key=lambda g: g.accession)
        random.Random(f"{opts.seed}:{species}").shuffle(options)
        options.sort(key=lambda g: quality_rank(g)[:2])
        shortlists[species] = options[:max(opts.per_species, opts.tech_candidates)] if technology else options
        before = set(previous["strains"].get(species, ()))
        kept[species] = [g for g in candidates[species] if g.accession in before]
    asked = {g.accession for s in shortlists.values() for g in s} | {g.accession for k in kept.values() for g in k}
    tech = technology(sorted(asked)) if technology else {}
    picked = []
    for species in sorted(chosen):
        ranked = sorted(shortlists[species], key=lambda g: quality_rank(g, bool(LONG_READS.search(tech.get(g.accession, "")))))
        ranked = kept[species] + [g for g in ranked if g not in kept[species]]
        picked += [(g, tech.get(g.accession, "")) for g in ranked[:opts.per_species]]
    return picked, sorted(chosen), sorted(rep_only), domain, representative, candidates


def gzip_into(source, dest):
    with open(source, "rb") as fin, gzip.open(dest + ".part", "wb", compresslevel=6) as fout:
        shutil.copyfileobj(fin, fout, 1 << 22)
    os.replace(dest + ".part", dest)


def ftp_url(base, accession, name):
    """Where NCBI's FTP server has an assembly: base/GCA/000/005/845/GCA_000005845.2_ASM584v2/
    GCA_000005845.2_ASM584v2_genomic.fna.gz (characters of the name other than letters, digits, . _ - become _)."""
    folder = f"{accession}_{re.sub('[^A-Za-z0-9._-]', '_', name)}"
    digits = accession[4:13]
    return f"{base}/{accession[:3]}/{digits[0:3]}/{digits[3:6]}/{digits[6:9]}/{folder}/{folder}_genomic.fna.gz"


def fetch_file(url, dest, attempts=None):
    """Downloads a gzip file to dest: "" if it arrived whole (the length it announced, and gzip reads it to its
    end), else why not. A missing file ("not found") or another refusal (HTTP 4xx but 429) is not asked for
    again; other failures are (ATTEMPTS times), after growing waits (RETRY_WAIT), or as long as the server
    asks for (Retry-After, at most a minute)."""
    why = ""
    part = dest + ".part"
    asked_to_wait = 0.0
    for attempt in range(attempts or ATTEMPTS):
        if attempt:
            time.sleep(max(RETRY_WAIT * 2 ** (attempt - 1), asked_to_wait))
        asked_to_wait = 0.0
        try:
            request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(request, timeout=60) as response, open(part, "wb") as fh:
                expected = response.headers.get("Content-Length")
                shutil.copyfileobj(response, fh, 1 << 20)
            if expected is not None and os.path.getsize(part) != int(expected):
                why = f"{os.path.getsize(part)} of {expected} bytes arrived"
                continue
            with gzip.open(part, "rb") as gz:  # CRC and length of the compressed stream
                while gz.read(1 << 22):
                    pass
            os.replace(part, dest)
            return ""
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return "not found"
            why = f"HTTP {e.code}"
            if 400 <= e.code < 500 and e.code != 429:
                break
            retry_after = e.headers.get("Retry-After", "") if e.headers else ""
            if retry_after.isdigit():
                asked_to_wait = min(60.0, float(retry_after))
        except (OSError, EOFError, http.client.HTTPException, zlib.error) as e:  # also timeouts, resets, bad gzip
            why = f"{type(e).__name__}: {e}"
    if os.path.exists(part):
        os.remove(part)
    return why


def fetch_direct(opts, wanted, folder):
    """Downloads wanted [(accession, url)] into folder as <accession>.fna.gz, --connections at a time. Returns
    ({accession delivered}, {accession: why not}). Missing files are no sign of trouble (an assembly NCBI renamed
    or withdrew goes through datasets), but FAILS_IN_A_ROW other failures in a row mean that NCBI does not
    answer: the rest is left undone, for datasets."""
    delivered, failed = set(), {}
    lock = threading.Lock()
    in_a_row = [0]
    started = time.time()
    size = [0]

    def one(item):
        accession, url = item
        if in_a_row[0] >= FAILS_IN_A_ROW:
            return accession, "not asked for: too many failures in a row"
        dest = os.path.join(folder, accession + ".fna.gz")
        why = fetch_file(url, dest)
        with lock:
            if why == "":
                in_a_row[0] = 0
                size[0] += os.path.getsize(dest)
            elif why != "not found":
                in_a_row[0] += 1
        return accession, why

    step = max(1, min(500, len(wanted) // 10))
    with concurrent.futures.ThreadPoolExecutor(max(1, opts.connections)) as executor:
        for n, (accession, why) in enumerate(executor.map(one, wanted), 1):
            (failed.__setitem__(accession, why) if why else delivered.add(accession))
            if n % step == 0 or n == len(wanted):
                seconds = max(time.time() - started, 1e-6)
                print(f"  {len(delivered)} of {n} fetched directly ({size[0] / 1e6:.0f} MB, "
                      f"{size[0] / 1e6 / seconds:.1f} MB/s, {n / seconds:.1f} genomes/s)", flush=True)
    return delivered, failed


def ncbi_batch(opts, accessions, work):
    """Downloads accessions with datasets into work: {accession: FASTA path} of those delivered, and why the
    request failed ("" if it did not). A package without fetch.txt is datasets' answer that it has none of the
    genomes (withdrawn or suppressed assemblies): none delivered, and no failure (rehydrate would fail on it, and
    a rerun asking only for such genomes took NCBI for down)."""
    shutil.rmtree(work, ignore_errors=True)
    os.makedirs(work)
    listing = os.path.join(work, "accessions.txt")
    with open(listing, "w") as fh:
        fh.writelines(a + "\n" for a in accessions)
    zip_path, unpacked = os.path.join(work, "batch.zip"), os.path.join(work, "batch")
    steps = ([opts.datasets, "download", "genome", "accession", "--inputfile", listing, "--include", "genome",
              "--dehydrated", "--filename", zip_path],
             ["unzip", "-o", "-q", zip_path, "-d", unpacked],
             [opts.datasets, "rehydrate", "--directory", unpacked, "--max-workers", str(max(1, min(opts.threads, 30)))])
    log_path = os.path.join(work, "log.txt")
    with open(log_path, "w") as log:
        for command in steps:
            if command[1] == "rehydrate" and not os.path.isfile(os.path.join(unpacked, "ncbi_dataset", "fetch.txt")):
                return {}, ""  # nothing to fetch: datasets has none of these genomes
            rc = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT).returncode
            if rc != 0:
                log.close()
                with open(log_path, errors="replace") as fh:
                    tail = " | ".join(line.strip() for line in fh.read().strip().splitlines()[-3:])
                return {}, f"{os.path.basename(command[0])} {command[1]} failed ({rc}): {tail}"
    found = {}
    for root, _dirs, names in os.walk(unpacked):
        for name in names:
            m = re.match(r"(GC[AF]_\d{9}\.\d+)_.*_genomic\.fna$", name)
            if m and m.group(1) in accessions:
                found[m.group(1)] = os.path.join(root, name)
    return found, ""


def require_datasets(opts):
    if not shutil.which(opts.datasets):
        sys.exit(f"the NCBI datasets CLI ({opts.datasets}) is not on PATH: install it (conda: ncbi-datasets-cli; see "
                 "envs/protal-db-build.yaml)")


def ncbi_requests(opts, items, request):
    """Yields (batch, found, requests left, the last failure) for the batches of items that request(batch) ->
    (found, why the request failed, "" if it did not) answers. A request that fails is halved, down to single
    items, to find the ones NCBI refuses; one such item fails at most ~log2(--batch) requests in a row. More
    failures in a row mean that NCBI (or datasets) is not working, which halving would only ask about ~2n
    times: the run stops. A request that succeeds is not halved, whatever it delivered: what it lacks NCBI does
    not have."""
    batches = [items[i:i + opts.batch] for i in range(0, len(items), opts.batch)]
    in_a_row, most_in_a_row, error = 0, max(1, opts.batch).bit_length() + 4, ""
    while batches:
        batch = batches.pop(0)
        found, why = request(batch)
        if why:
            in_a_row, error = in_a_row + 1, why
            if in_a_row >= most_in_a_row:
                sys.exit(f"{in_a_row} NCBI requests failed in a row, the last: {why}. Is NCBI reachable from "
                         f"here, and does {opts.datasets} work? A rerun downloads what is still missing")
        else:
            in_a_row = 0
        if why and not found and len(batch) > 1:  # a failed request: halve it, down to single items
            batches[:0] = [batch[:len(batch) // 2], batch[len(batch) // 2:]]
            continue
        yield batch, found, len(batches), error


def ncbi_summary(opts, accessions, work):
    """The sequencing technology (a string, "" if NCBI has none) of accessions, from `datasets summary genome
    accession`: ({accession: technology} of those NCBI answered, and why the request failed ("" if it did not))."""
    os.makedirs(work, exist_ok=True)
    listing = os.path.join(work, "accessions.txt")
    with open(listing, "w") as fh:
        fh.writelines(a + "\n" for a in accessions)
    result = subprocess.run([opts.datasets, "summary", "genome", "accession", "--inputfile", listing,
                             "--as-json-lines"], capture_output=True, text=True, errors="replace")
    if result.returncode != 0:
        tail = " | ".join(line.strip() for line in (result.stderr or result.stdout).strip().splitlines()[-3:])
        return {}, f"datasets summary failed ({result.returncode}): {tail}"
    found = {}
    for line in result.stdout.splitlines():
        try:
            report = json.loads(line)
        except ValueError:
            continue
        accession = flat_get(report, "accession")
        if isinstance(accession, str) and accession in accessions:
            tech = flat_get(report, "assemblyinfo", "sequencingtech")
            found[accession] = " ".join(tech.split()) if isinstance(tech, str) else ""
    return found, ""


def flat_get(report, *path):
    """report[path[0]][path[1]]..., the keys compared in lower case without underscores: datasets writes
    assemblyInfo/sequencingTech or assembly_info/sequencing_tech, by version."""
    for name in path:
        if not isinstance(report, dict):
            return None
        report = next((v for k, v in report.items() if k.replace("_", "").lower() == name), None)
    return report


def technologies(opts, accessions):
    """{accession: sequencing technology at NCBI} of those NCBI answered for, from OUT/ncbi_info.tsv and what
    it is asked for (and then added to it, so that a rerun asks for nothing it knows)."""
    path = os.path.join(opts.out, "ncbi_info.tsv")
    known = {}
    if os.path.isfile(path):
        with open(path) as fh:
            next(fh, None)
            known = {f[0]: f[1] if len(f) > 1 else "" for f in (line.rstrip("\n").split("\t") for line in fh) if f[0]}
    todo = [a for a in sorted(set(accessions)) if a not in known]
    print(f"sequencing technology of {len(set(accessions))} candidate strains: {len(known.keys() & set(accessions))} "
          f"known, {len(todo)} to ask NCBI for", flush=True)
    if todo:
        require_datasets(opts)
        work = os.path.join(opts.out, "ncbi_summary")
        new = not os.path.isfile(path)
        with open(path, "a") as fh:
            if new:
                fh.write("accession\tsequencing_tech\n")
            for batch, found, left, _ in ncbi_requests(opts, todo, lambda b: ncbi_summary(opts, b, work)):
                fh.writelines(f"{a}\t{t}\n" for a, t in found.items())
                fh.flush()
                known.update(found)
                print(f"  {len(found)} of {len(batch)} answered; {left} requests left", flush=True)
        shutil.rmtree(work, ignore_errors=True)
    return {a: known[a] for a in set(accessions) if a in known}


def count_text(counter, order=()):
    """'3 isolate, 1 mag' from a Counter, the keys in `order` first."""
    keys = [k for k in order if counter.get(k)] + sorted(k for k in counter if k not in order)
    return ", ".join(f"{counter[k]} {k or 'unknown'}" for k in keys) or "none"


def get_genomes(opts, state, release):
    gtdb = os.path.join(opts.out, "release")
    genomes = read_metadata(gtdb, release)
    allowed = None
    if opts.progenomes:
        allowed = read_accessions(opts.progenomes)
        if not allowed:
            sys.exit(f"{opts.progenomes} lists no GCA_/GCF_ accessions: is it proGenomes' pg4_ANI_clustering.tsv.gz?")
        print(f"proGenomes: {len(allowed)} genomes in {opts.progenomes}; strains are taken from these only", flush=True)
    previous = previous_choice(opts.out)
    strains, chosen, rep_only, domain, representative, candidates = pick(
        genomes, opts, None if opts.no_tech_lookup else lambda accessions: technologies(opts, accessions), allowed,
        previous)
    pool = sorted(chosen + rep_only)
    eligible = sum(1 for s in domain if candidates.get(s))
    kept = sum(1 for s in chosen if s in previous["strains"])
    converted = sum(1 for s in chosen if s in previous["rep_only"])
    print(f"species: {len(chosen)} with strains (--species {opts.species}; {eligible} of GTDB's {len(domain)} species have "
          f"a strain passing the filters" + (f"; {kept} kept from the earlier choice, {converted} of its representative-"
                                             "only species given strains" if previous["strains"] or previous["rep_only"] else "")
          + f"), {len(rep_only)} from their representative only (--rep_only_species {opts.rep_only_species})", flush=True)
    if len(chosen) < round(opts.species * 0.99):
        print(f"  only {len(chosen)} species have strains passing the filters (--min_completeness, --max_contamination"
              f"{', --progenomes' if allowed else ''}): fewer than --species asks", flush=True)
    # (accession, species, role, lineage, CheckM2 completeness and contamination, category, level, contigs, technology)
    wanted = [(g.accession, g.species, "strain", g.lineage, g.completeness, g.contamination, g.category, g.level,
               "" if g.contigs is None else g.contigs, tech) for g, tech in strains]
    if opts.rep_genomes == "ncbi":
        wanted += [(g.accession, g.species, "representative", g.lineage, g.completeness, g.contamination, g.category,
                    g.level, "" if g.contigs is None else g.contigs, "")
                   for g in (representative[s] for s in pool if s in representative)]
    folder = os.path.join(opts.out, "genomes")
    os.makedirs(folder, exist_ok=True)
    with open(os.path.join(opts.out, "simulation_species.txt"), "w") as fh:
        fh.writelines(s + "\n" for s in pool)
    have = {name[:-len(".fna.gz")] for name in os.listdir(folder) if name.endswith(".fna.gz")}
    todo = [w[0] for w in wanted if w[0] not in have]
    print(f"genomes: {len(wanted)} wanted ({len(strains)} strains of {len(chosen)} species"
          + (f", the representatives of {len(pool)} species" if opts.rep_genomes == "ncbi" else "")
          + f"), {len(wanted) - len(todo)} already there, {len(todo)} to download", flush=True)
    # Straight from NCBI's FTP server, in parallel; what that does not deliver goes through datasets.
    fetched = 0
    if todo and opts.connections > 0 and opts.ftp_url:
        names = {g.accession: g.name for g in genomes}
        direct = [(a, ftp_url(opts.ftp_url.rstrip("/"), a, names[a])) for a in todo if names.get(a)]
        if direct:
            print(f"fetching {len(direct)} genomes from {opts.ftp_url}, {opts.connections} at a time", flush=True)
            delivered_direct, failed_direct = fetch_direct(opts, direct, folder)
            fetched = len(delivered_direct)
            todo = [a for a in todo if a not in delivered_direct]
            reasons = collections.Counter(failed_direct.values())
            print(f"  {fetched} fetched directly" + (f"; {len(todo)} left for datasets ("
                  + "; ".join(f"{n} {why}" for why, n in reasons.most_common(3)) + ")" if todo else ""), flush=True)
            if reasons["HTTP 503"] or reasons["HTTP 429"]:
                print(f"  NCBI limited the requests: use fewer --connections (now {opts.connections})", flush=True)
    if todo:
        require_datasets(opts)
    missing = []
    work = os.path.join(opts.out, "ncbi_batch")
    error = ""
    with concurrent.futures.ProcessPoolExecutor(max(1, opts.threads)) as pool_executor:
        for batch, found, left, error in ncbi_requests(opts, todo, lambda b: ncbi_batch(opts, b, work)):
            list(pool_executor.map(gzip_into, found.values(), [os.path.join(folder, a + ".fna.gz") for a in found]))
            missing += [a for a in batch if a not in found]
            print(f"  {len(found)} of {len(batch)} delivered; {left} requests left", flush=True)
    shutil.rmtree(work, ignore_errors=True)
    have = {name[:-len(".fna.gz")] for name in os.listdir(folder) if name.endswith(".fna.gz")}
    if wanted and not any(w[0] in have for w in wanted):
        sys.exit(f"NCBI delivered none of the {len(wanted)} genomes" + (f"; the last failure: {error}" if error else ""))
    # build_gtdb_database.py reads every genome in genomes/, so those of an earlier choice (other options, or a
    # choice by chance before strains were ranked by quality) are moved out, not deleted.
    stray = sorted(have - {w[0] for w in wanted})
    if stray:
        unused = os.path.join(opts.out, "genomes_unused")
        os.makedirs(unused, exist_ok=True)
        for accession in stray:
            os.replace(os.path.join(folder, accession + ".fna.gz"), os.path.join(unused, accession + ".fna.gz"))
        print(f"genomes: {len(stray)} genomes that this choice does not include moved to {unused} (a build reads all "
              f"of {folder}); delete them to free the space", flush=True)
    with open(os.path.join(opts.out, "genomes.tsv"), "w") as fh:
        fh.write("accession\tspecies\trole\tlineage\tcheckm2_completeness\tcheckm2_contamination\tgenome_category\t"
                 "assembly_level\tcontig_count\tsequencing_tech\n")
        fh.writelines("\t".join(map(str, w)) + "\n" for w in wanted if w[0] in have)
    with open(os.path.join(opts.out, "missing.txt"), "w") as fh:
        fh.writelines(a + "\n" for a in sorted(w[0] for w in wanted if w[0] not in have))
    delivered = [w for w in wanted if w[0] in have]
    per_species = collections.Counter(w[1] for w in delivered if w[2] == "strain")
    other = sum(k / (k + 1) for k in per_species.values()) / max(1, len(pool))
    # How good the strains are, and what choice there was: species whose candidates were all MAGs, and
    # representatives that are MAGs where an isolate strain exists.
    picked_strains = [w for w in delivered if w[2] == "strain"]
    known_tech = [w for w in picked_strains if w[9]]
    long_read = [w for w in known_tech if LONG_READS.search(w[9])]
    mag_only = sum(1 for s in chosen if all(g.category == "mag" for g in candidates[s]))
    mag_reps = sum(1 for s in pool if s in representative and representative[s].category == "mag"
                   and any(g.category == "isolate" for g in candidates.get(s, ())))
    state["genomes"] = {"wanted": len(wanted), "delivered": len(delivered), "missing": len(wanted) - len(delivered),
                        "strains": sum(per_species.values()), "species_with_strains": len(per_species),
                        "simulation_species": len(pool), "other_strain_share": round(other, 3),
                        "by_domain": dict(collections.Counter(domain[s] for s in pool)),
                        "strain_category": dict(collections.Counter(w[6] for w in picked_strains)),
                        "strain_assembly_level": dict(collections.Counter(w[7] or "unknown" for w in picked_strains)),
                        "strains_with_known_technology": len(known_tech), "strains_long_read": len(long_read),
                        "species_with_mags_only": mag_only, "mag_representatives_with_isolate_strains": mag_reps,
                        "moved_to_genomes_unused": len(stray), "fetched_directly": fetched}
    print(f"genomes: {len(delivered)} of {len(wanted)} in {folder} ({len(wanted) - len(delivered)} missing, see "
          f"missing.txt); a simulated species is another strain than the representative {100 * other:.0f}% of the time",
          flush=True)
    print(f"strains: {count_text(collections.Counter(w[6] for w in picked_strains), CATEGORIES)}; assembly level "
          f"{count_text(collections.Counter(w[7] for w in picked_strains), LEVELS)}; {len(long_read)} of "
          f"{len(known_tech)} with a known sequencing technology are PacBio or Nanopore assemblies; {mag_only} of "
          f"{len(chosen)} species have MAGs only; {mag_reps} representatives are MAGs of species that have isolate "
          "strains (left as GTDB's references)", flush=True)


def get_host(opts, state):
    """The host genome (--host_genome) as OUT/host/<accession>.fna.gz: from NCBI's FTP server, kept gzipped as it is
    served (and checked through gzip), or through datasets, whose plain FASTA is gzipped here. Kept from an earlier run
    when there. state["host"]: accession, assembly, path (relative to OUT) and size."""
    spec = (opts.host_genome or "").strip()
    if spec.lower() in ("", "none"):
        return
    if spec in HOST_GENOMES:
        accession, name = HOST_GENOMES[spec]
    else:
        m = re.fullmatch(r"(GC[AF]_\d{9}\.\d+)_(.+)", spec)
        if not m:
            sys.exit(f"--host_genome {spec!r}: one of {', '.join(HOST_GENOMES)}, ACCESSION_ASSEMBLYNAME (e.g. "
                     "GCF_000001405.40_GRCh38.p14), or none")
        accession, name = m.groups()
    folder = os.path.join(opts.out, "host")
    os.makedirs(folder, exist_ok=True)
    dest = os.path.join(folder, accession + ".fna.gz")
    if os.path.isfile(dest):
        print(f"host genome: {dest} is there ({os.path.getsize(dest) / 1e9:.2f} GB)", flush=True)
    else:
        why = "no FTP server (--ftp_url)"
        if opts.connections > 0 and opts.ftp_url:
            url = ftp_url(opts.ftp_url.rstrip("/"), accession, name)
            print(f"host genome: fetching {accession} ({name}) from {url}", flush=True)
            why = fetch_file(url, dest)
        if why and shutil.which(opts.datasets):
            print(f"host genome: not fetched directly ({why}); asking datasets", flush=True)
            work = os.path.join(opts.out, "host_batch")
            found, error = ncbi_batch(opts, [accession], work)
            if accession in found:
                gzip_into(found[accession], dest)
                why = ""
            else:
                why = error or f"{why}; datasets did not deliver it either"
            shutil.rmtree(work, ignore_errors=True)
        if why:  # the GTDB inputs serve every build without it: only the scenarios with host reads need it
            print(f"host genome {accession}: not downloaded ({why}). Only the scenarios with host reads "
                  "(build_gtdb_database.py --scenarios host) need it: rerun, or give the build --host-genome", flush=True)
            state.pop("host", None)
            return
        print(f"host genome: {dest} ({os.path.getsize(dest) / 1e9:.2f} GB)", flush=True)
    state["host"] = {"name": spec, "accession": accession, "assembly": name,
                     "path": os.path.relpath(dest, opts.out), "bytes": os.path.getsize(dest)}


# ---- main -------------------------------------------------------------------------------------------------

def load_state(opts):
    path = os.path.join(opts.out, "download.json")
    if os.path.exists(path):
        with open(path) as fh:
            return json.load(fh)
    return {}


def save_state(opts, state):
    path = os.path.join(opts.out, "download.json")
    with open(path + ".part", "w") as fh:
        json.dump(state, fh, indent=1)
    os.replace(path + ".part", path)


def main(argv=None):
    opts = parse_args(argv)
    os.makedirs(opts.out, exist_ok=True)
    state = load_state(opts)
    release = get_release(opts, state)
    state["options"] = {k: v for k, v in vars(opts).items() if k not in ("out", "datasets", "threads")}
    save_state(opts, state)
    if not opts.no_genomes:
        get_genomes(opts, state, release)
        get_host(opts, state)
    state["updated"] = datetime.datetime.now().isoformat(timespec="seconds")
    save_state(opts, state)
    print(f"Inputs for GTDB r{release}: {opts.out} (build_gtdb_database.py --inputs {opts.out})", flush=True)


if __name__ == "__main__":
    main()
