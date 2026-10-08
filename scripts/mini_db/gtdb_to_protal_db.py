#!/usr/bin/env python3
"""gtdb_to_protal_db.py - turn an (extracted) GTDB release into the input files
of a protal database folder, ready for `protal --build`.

Reads from the GTDB release directory:
  {bac120,ar53}_taxonomy_r<R>.tsv[.gz]     accession -> GTDB lineage
  {bac120,ar53}_metadata_r<R>.tsv[.gz]     species representatives, CheckM quality and
                                           genome sizes (optional; without it, genomes
                                           present in the *_marker_genes_reps_* files
                                           are the reps, and sizes are unknown)
  genomic_files_reps/*_marker_genes_reps_r<R>/**.fna[.gz]   marker genes of reps
  genomic_files_all/*_marker_genes_all_r<R>/**.fna[.gz]     marker genes of all
                                           genomes (optional, for unique k-mers)
The marker gene tarballs must be extracted first. One FASTA per marker; the
marker id (PFxxxxx.x / TIGRxxxxx) is taken from the file name and the genome
accession from the record header.

Writes to <outdir>:
  reference.fna          representative marker genes, header >taxid_geneid, one
                         sequence line per record; ordered by gene, then taxid (taxids follow
                         the taxonomy), so that a gene's copies in related species are
                         neighbours: zstd compresses this ~2x better than genome by genome
                         when related genomes are far apart in the file (--order genome for
                         the old order). protal finds genes through reference.map, in any order.
  reference.map          taxid, geneid, start byte, end byte of each sequence line
  internal_taxonomy.dmp  id, parent_id, external_id, name, rank, level, rep_genome
                         (species are the leaves; their ids are the reference taxids)
  full_reference.fna.zst marker genes of all genomes, header >taxid_geneid of the
                         genome's species (only if genomic_files_all is present);
                         zstd-compressed by the zstd command (86 GB raw at r226), a frame
                         per marker file (or gene, in a --from_db copy) and a seek table,
                         so that protal --build reads it on several threads;
                         full_reference.fna without the zstd command
  gene2geneid.tsv        marker id -> geneid
  genome2tiid.tsv        accession, species taxid, species rep accession, lineage
  species_priors.tsv     per species: the representative's markers, CheckM quality, GTDB's species
                         cluster, and (since 2026-10-08) the species' genome size (the mean of its
                         genomes' sizes corrected by CheckM), the representative's assembly size and
                         the bases of its marker genes; --priors_only writes this file alone
  model_pe.xml           copy of --model (the profiler's random forest for paired-end reads;
                         models for other read types: protal --add_model, see README)

Then build the index with
  protal --build --no_profile --db <outdir> --reference <outdir>/reference.fna \\
         --full_reference <outdir>/full_reference.fna
which writes database.protal (index.prx.zst and reference.fna.zst with --no_bundle; raw files with
--no_compress); for --full_reference full_reference.fna, protal reads its .zst.

Each marker's genes are spooled to <outdir>/.convert_tmp (or <--tmp>/.convert_tmp) and sorted one gene
at a time, so memory stays at about one gene's sequences per worker; -t reads the marker files in
parallel (the output is the same for any -t). The workers compress their chunks of the full reference
themselves (zstd frames, joined as they are, then the seek table), and the log gives each step's time.

Usage:
  gtdb_to_protal_db.py --gtdb <release dir> --outdir <db dir> [--release 226] [--model FILE]
      [--exclude_species FILE] [--genes LIST]
  gtdb_to_protal_db.py --from_db <db dir> [--exclude_species FILE] [--genes LIST] --outdir <copy>
  gtdb_to_protal_db.py --gtdb <release dir> --priors_only --outdir <dir>   (then protal --add_tables)

--genes keeps a subset of the marker genes (a reduced database: less memory, fewer hits): GTDB
marker ids (PF00380.20, TIGR00001; PF00380 without its version matches any) or protal gene ids
(the numbers of gene2geneid.tsv), comma-separated or in a file (one per line, the first column, #
comments; e.g. the subset file of scripts/rank_genes.py). The genes keep their ids; the other
genes are left out of reference.fna, reference.map, the full reference and gene2geneid.tsv, and
the gene neighbours of a --from_db copy are counted anew over the genes kept (from
gene_positions.tsv: the gene a read meets next is the nearest one of the subset).

--exclude_species leaves the marker genes of those species out and keeps them in the taxonomy
with their taxids, for a training database: reads of the species left out land on relatives, as
those of species GTDB lacks do in real samples. --from_db makes such a copy of a folder this script
wrote, before protal --build packs it, without reading the release again.
"""

import argparse
import contextlib
import glob
import gzip
import heapq
import multiprocessing
import os
import re
import shutil
import struct
import subprocess
import sys
import time

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FULL_REFERENCE = "full_reference.fna"
# zstd settings of full_reference.fna.zst (86 GB raw at GTDB r226): level 6 with a 128 MB long-distance window
# packed a GTDB-like world's 12.8x at ~60 MB/s per thread, and decompresses at ~800 MB/s
# (docs/claude/2026-10-01-scratch-space). protal --build reads it.
ZSTD_OPTIONS = ("-6", "--long=27")
MARKER_SETS = ("bac120", "ar53")
RANKS = [("d", "domain"), ("p", "phylum"), ("c", "class"), ("o", "order"),
         ("f", "family"), ("g", "genus"), ("s", "species")]
ACCESSION_RE = re.compile(r"(GC[AF]_\d{9}\.\d+)")
MARKER_RE = re.compile(r"(PF\d{5}\.\d+|TIGR\d{5})")
FASTA_RE = re.compile(r"\.(fna|fa|fasta)(\.gz)?$")


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def find_one(directory, stem):
    """<stem>.tsv, else <stem>.tsv.gz, else None."""
    for ext in (".tsv", ".tsv.gz"):
        if os.path.exists(os.path.join(directory, stem + ext)):
            return os.path.join(directory, stem + ext)
    return None


def detect_release(gtdb):
    tags = sorted({m.group(1) for f in os.listdir(gtdb)
                   for m in [re.match(r"(?:bac120|ar53)_taxonomy_r(\d+)", f)] if m})
    if len(tags) != 1:
        sys.exit(f"Cannot detect the GTDB release in {gtdb} (found: {tags or 'none'}); pass --release")
    return tags[0]


def normalize_accession(text):
    """'RS_GCF_000005845.2 some description' -> 'GCF_000005845.2'."""
    m = ACCESSION_RE.search(text)
    if m:
        return m.group(1)
    token = text.split()[0] if text.split() else text
    return token[3:] if token[:3] in ("RS_", "GB_") else token


def allele_genome(accession, share):
    """Whether the genome `accession` (normalize_accession's form, which the full reference names after each record's
    name) may give strain alleles: the top 53 bits of its FNV-1a 64 hash as a fraction below `share`. protal --build
    (--allele_genome_share, StrainAlleles.h AlleleGenome) and ancestry_sites.py take alleles from these genomes,
    build_gtdb_database.py simulates strains from the others, so that no simulated strain is its species' own allele."""
    if share >= 1:
        return True
    if share <= 0:
        return False
    h = 1469598103934665603
    for byte in accession.encode():
        h = ((h ^ byte) * 1099511628211) & 0xFFFFFFFFFFFFFFFF
    return (h >> 11) / float(1 << 53) < share


def read_taxonomy(gtdb, rel):
    lineage = {}
    for mset in MARKER_SETS:
        path = find_one(gtdb, f"{mset}_taxonomy_r{rel}")
        if not path:
            continue
        with open_text(path) as fh:
            for line in fh:
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 2 and fields[1].startswith("d__"):
                    lineage[normalize_accession(fields[0])] = fields[1].strip()
    if not lineage:
        sys.exit(f"No {'/'.join(MARKER_SETS)}_taxonomy_r{rel}.tsv[.gz] found in {gtdb}")
    return lineage


def read_metadata(gtdb, rel):
    """-> ({accession: (completeness, contamination)}, {accession: (genome size, species representative)}) from the
    metadata: CheckM's columns (checkm2_* where the release has them, else checkm_*), genome_size and
    gtdb_genome_representative (the accession itself where the column is missing); a missing or empty value is None.
    Every genome of the release, not only the representatives. Empty dicts without metadata or columns."""
    quality, sizes = {}, {}

    def number(fields, col):
        try:
            return float(fields[col]) if col is not None and fields[col] not in ("", "none", "N/A", "NA") else None
        except (ValueError, IndexError):
            return None

    for mset in MARKER_SETS:
        path = find_one(gtdb, f"{mset}_metadata_r{rel}")
        if not path:
            continue
        with open_text(path) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            if "accession" not in header:
                continue
            acc_col = header.index("accession")
            cols = []
            for what in ("completeness", "contamination"):
                col = next((header.index(c) for c in (f"checkm2_{what}", f"checkm_{what}") if c in header), None)
                cols.append(col)
            size_col = header.index("genome_size") if "genome_size" in header else None
            rep_col = header.index("gtdb_genome_representative") if "gtdb_genome_representative" in header else None
            for line in fh:
                fields = line.rstrip("\n").split("\t")
                if len(fields) <= acc_col:
                    continue
                acc = normalize_accession(fields[acc_col])
                if any(c is not None for c in cols):
                    quality[acc] = tuple(number(fields, col) for col in cols)
                size = number(fields, size_col)
                if size is not None and size > 0:
                    rep = fields[rep_col] if rep_col is not None and rep_col < len(fields) else ""
                    sizes[acc] = (size, normalize_accession(rep) if rep not in ("", "none", "N/A", "NA") else acc)
    return quality, sizes


# A genome's size is taken as its assembly's corrected by CheckM: size * (100 - contamination) / completeness
# (Rodriguez-Gijon et al. 2022, Nat Ecol Evol). A species' size is the mean over its high-quality genomes (MIMAG:
# completeness >= 90, contamination <= 5), else over its medium-quality ones (>= 50, <= 10), else the mean of its
# genomes' assembly sizes as they are (a release without CheckM's columns, or no genome good enough).
SIZE_QUALITY_TIERS = ((90.0, 5.0), (50.0, 10.0))


def species_genome_sizes(quality, sizes):
    """-> {species representative: (genome size, genomes averaged)} from read_metadata's dicts."""
    by_rep = {}
    for acc, (size, rep) in sizes.items():
        by_rep.setdefault(rep, []).append((size, *quality.get(acc, (None, None))))
    result = {}
    for rep, genomes in by_rep.items():
        chosen = None
        for min_comp, max_cont in SIZE_QUALITY_TIERS:
            corrected = [size * (100.0 - cont) / comp for size, comp, cont in genomes
                         if comp is not None and cont is not None and comp >= min_comp and cont <= max_cont]
            if corrected:
                chosen = corrected
                break
        if chosen is None:
            chosen = [size for size, _, _ in genomes]
        result[rep] = (sum(chosen) / len(chosen), len(chosen))
    return result


SP_CLUSTERS_COLUMNS = {"radius": "ANI circumscription radius", "mean_ani": "Mean intra-species ANI",
                       "min_ani": "Min intra-species ANI", "genomes": "No. clustered genomes"}


def read_sp_clusters(gtdb, rel):
    """-> {representative accession: {radius, mean_ani, min_ani, genomes}} from GTDB's sp_clusters_r<rel>.tsv (in
    auxillary_files/, GTDB's spelling, or beside the metadata), None where a value is N/A; {} without the file."""
    path = None
    for folder in (os.path.join(gtdb, "auxillary_files"), os.path.join(gtdb, "auxiliary_files"), gtdb):
        path = find_one(folder, f"sp_clusters_r{rel}") if os.path.isdir(folder) else None
        if path:
            break
    if not path:
        return {}
    clusters = {}
    with open_text(path) as fh:
        header = [h.strip() for h in fh.readline().rstrip("\n").split("\t")]
        try:
            acc_col = header.index("Representative genome")
            cols = {key: header.index(name) for key, name in SP_CLUSTERS_COLUMNS.items()}
        except ValueError:
            sys.exit(f"{path} lacks the columns 'Representative genome', " + ", ".join(f"'{n}'" for n in SP_CLUSTERS_COLUMNS.values()))
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= max(cols.values()):
                continue
            row = {}
            for key, col in cols.items():
                try:
                    row[key] = float(fields[col]) if fields[col] not in ("", "N/A", "NA", "none") else None
                except ValueError:
                    row[key] = None
            clusters[normalize_accession(fields[acc_col])] = row
    return clusters


def read_representatives(gtdb, rel):
    """-> set of representative accessions from the metadata, or None if there is none."""
    reps, found = set(), False
    for mset in MARKER_SETS:
        path = find_one(gtdb, f"{mset}_metadata_r{rel}")
        if not path:
            continue
        found = True
        with open_text(path) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            try:
                acc_col, rep_col = header.index("accession"), header.index("gtdb_representative")
            except ValueError:
                sys.exit(f"{path} lacks the accession/gtdb_representative columns")
            for line in fh:
                fields = line.rstrip("\n").split("\t")
                if len(fields) > rep_col and fields[rep_col] == "t":
                    reps.add(normalize_accession(fields[acc_col]))
    return reps if found else None


def marker_files(gtdb, subdir, kind, rel):
    """-> [(marker id, path)] of the nucleotide marker files, sorted by set then marker."""
    files = []
    for mset in MARKER_SETS:
        root = os.path.join(gtdb, subdir, f"{mset}_marker_genes_{kind}_r{rel}")
        if not os.path.isdir(root):
            if os.path.exists(root + ".tar.gz"):
                sys.exit(f"Extract {root}.tar.gz first")
            continue
        found = []
        for dirpath, _dirs, names in os.walk(root):
            if os.path.basename(dirpath) == "faa":
                continue
            for name in names:
                if FASTA_RE.search(name):
                    m = MARKER_RE.search(name)
                    marker = m.group(1) if m else FASTA_RE.sub("", name).split("_", 1)[-1]
                    found.append((marker, os.path.join(dirpath, name)))
        files += sorted(found)
    return files


def read_fasta(path):
    header, chunks = None, []
    with open_text(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line[1:], []
            elif line:
                chunks.append(line)
    if header is not None:
        yield header, "".join(chunks)


_WORK = {}  # lineage, reps, tmp, gene_ids, taxid, drop: set before workers fork (_parallel)


def _parallel(threads, function, items):
    """function over items, in threads forked workers (which see _WORK as it is now), in order."""
    if threads <= 1:
        return list(map(function, items))
    with multiprocessing.get_context("fork").Pool(threads) as pool:
        return pool.map(function, items, chunksize=1)


def _spool_representatives(item):
    """One marker: the representatives' records of all its files (a genome's first copy), spooled
    to rep_<gene id>.tsv; -> (marker, [(accession, bases of its copy)], skipped counts, the genomes with more copies)."""
    marker, paths = item
    lineage, reps = _WORK["lineage"], _WORK["reps"]
    counts = {"not in taxonomy": 0, "not a representative": 0, "duplicate": 0}
    kept, accessions, duplicated = set(), [], set()
    with open(os.path.join(_WORK["tmp"], f"rep_{_WORK['gene_ids'][marker]}.tsv"), "w", newline="\n") as out:
        for path in paths:
            for header, seq in read_fasta(path):
                acc = normalize_accession(header)
                if acc not in lineage:
                    counts["not in taxonomy"] += 1
                elif reps is not None and acc not in reps:
                    counts["not a representative"] += 1
                elif acc in kept:
                    counts["duplicate"] += 1
                    duplicated.add(acc)  # a single-copy marker twice in one genome: CheckM's contamination signature
                else:
                    kept.add(acc)
                    accessions.append((acc, len(seq)))
                    out.write(f"{acc}\t{seq.upper()}\n")
    return marker, accessions, counts, sorted(duplicated)


def _write_gene(marker):
    """One gene: its records with their species' taxids (without the species left out), sorted by
    taxid, as gene_<gene id>.fna; -> [(taxid, sequence length)]."""
    gid = _WORK["gene_ids"][marker]
    lineage, taxid, drop = _WORK["lineage"], _WORK["taxid"], _WORK["drop"]
    records = []
    with open(os.path.join(_WORK["tmp"], f"rep_{gid}.tsv")) as fh:
        for line in fh:
            acc, seq = line.rstrip("\n").split("\t")
            tid = taxid[lineage[acc].split(";")[-1]]
            if tid not in drop:
                records.append((tid, seq))
    records.sort(key=lambda r: r[0])
    with open(os.path.join(_WORK["tmp"], f"gene_{gid}.fna"), "wb") as out:
        for tid, seq in records:
            out.write(f">{tid}_{gid}\n{seq}\n".encode())
    return [(tid, len(seq)) for tid, seq in records]


def _gene_records(path, gid):
    """(taxid, gene id, sequence) of a gene_<gene id>.fna, in its order (by taxid)."""
    with open(path, "rb") as fh:
        for header in fh:
            yield int(header[1:].split(b"_", 1)[0]), gid, fh.readline().rstrip(b"\n")


def _write_full_reference(item):
    """One marker's files of all genomes: records of the database's species (not those left out), a
    genome's first copy, one chunk per file, zstd-compressed (a frame of its own: the chunks joined are
    the compressed full reference) when there is the zstd command (_WORK["zstd"]); -> {file index:
    (chunk, records)}."""
    marker, files = item
    gid = _WORK["gene_ids"][marker]
    lineage, taxid, drop = _WORK["lineage"], _WORK["taxid"], _WORK["drop"]
    seen, result = set(), {}
    for index, path in files:
        chunk, n, size = os.path.join(_WORK["tmp"], f"full_{index}.fna" + (".zst" if _WORK["zstd"] else "")), 0, 0
        with (zstd_writer(chunk, 1) if _WORK["zstd"] else open(chunk, "wb")) as fh:
            for header, seq in read_fasta(path):
                acc = normalize_accession(header)
                sp = lineage.get(acc, "").split(";")[-1]
                if sp not in taxid or taxid[sp] in drop or acc in seen:
                    continue
                seen.add(acc)
                # The genome's accession after the name, so that a report can leave a read's own genome out of a
                # species' alleles (ancestry_sites.py); protal and the tiler read the first word.
                record = f">{taxid[sp]}_{gid} {acc}\n{seq.upper()}\n".encode()
                fh.write(record)
                n += 1
                size += len(record)
        result[index] = (chunk, n, size)
    return result


def seek_table(frames):
    """zstd's seekable format's seek table (a skippable frame) for frames given as (compressed size, content size):
    protal --build finds the full reference's frames through it and decompresses them on several threads; the zstd
    command skips it. None if a frame is too large for the table's 4-byte entries."""
    if any(c >= 1 << 32 or d >= 1 << 32 for c, d in frames):
        return None
    entries = b"".join(struct.pack("<II", c, d) for c, d in frames)
    return (struct.pack("<II", 0x184D2A5E, len(entries) + 9) + entries +
            struct.pack("<IBI", len(frames), 0, 0x8F92EAB1))


def read_species_list(path):
    """Species names, one per line (first tab-separated field; 's__' optional; '#' comments)."""
    names = set()
    with open(path) as fh:
        for line in fh:
            name = line.rstrip("\n").split("\t")[0].strip()
            if name and not name.startswith("#"):
                names.add(name if name.startswith("s__") else "s__" + name)
    return names


def species_taxids(taxonomy_rows, names, what):
    """Taxids of the species named in `names` among (taxid, name, rank) rows; exits naming any not found."""
    found = {name: tid for tid, name, rank in taxonomy_rows if rank == "species" and name in names}
    missing = sorted(names - set(found))
    if missing:
        sys.exit(f"{len(missing)} species of {what} are not in the taxonomy: " + ", ".join(missing[:10])
                 + (" ..." if len(missing) > 10 else ""))
    return set(found.values())


# The files of a database folder this script writes besides reference.fna, reference.map and
# full_reference.fna, and those build_gtdb_database.py puts beside them: the models (model_pmml.MODEL_FILES)
# and gene_neighbours.tsv (gene_neighbours.py; per clade). A copy without some species derives
# gene_neighbours.tsv anew from gene_positions.tsv, without their genomes, when the folder has one.
CONVERTED_FILES = ("internal_taxonomy.dmp", "gene2geneid.tsv", "genome2tiid.tsv", "gene_neighbours.tsv",
                   "species_priors.tsv", "model_pe.xml", "model_se.xml", "model_PB.xml", "model_ONT.xml")
GENE_POSITIONS = "gene_positions.tsv"
# What protal --build writes into the folder (and leaves there when stopped: .partial files), stale once
# the folder's reference is written anew; it would stop the next build (unique_kmers.tsv of other genes)
# or shadow the new files (database.protal).
BUILD_OUTPUTS = ("index.prx", "index.prx.zst", "reference.fna.zst", "unique_kmers.tsv", "gene_conservation.tsv", "gene_congeners.tsv",
                 "database.protal", "build_metadata.tsv")


def full_reference_path(folder):
    """A converted folder's full reference: full_reference.fna.zst, or full_reference.fna (written so without the
    zstd command); None if it has neither."""
    for name in (FULL_REFERENCE + ".zst", FULL_REFERENCE):
        if os.path.isfile(os.path.join(folder, name)):
            return os.path.join(folder, name)
    return None


def remove_full_reference(folder):
    """Removes both variants of a folder's full reference; the bytes removed."""
    removed = 0
    for name in (FULL_REFERENCE + ".zst", FULL_REFERENCE):
        path = os.path.join(folder, name)
        if os.path.isfile(path):
            removed += os.path.getsize(path)
            os.remove(path)
    return removed


@contextlib.contextmanager
def zstd_writer(path, threads):
    """A binary file to write `path` through the zstd command (ZSTD_OPTIONS, up to 8 of `threads`): written as
    path.partial, renamed once complete."""
    partial = path + ".partial"
    process = subprocess.Popen([shutil.which("zstd") or "zstd", "-q", "-f", f"-T{max(1, min(threads, 8))}", *ZSTD_OPTIONS,
                                "-o", partial], stdin=subprocess.PIPE)
    try:
        yield process.stdin
    finally:
        process.stdin.close()
        rc = process.wait()
    if rc:
        sys.exit(f"zstd failed with exit code {rc} writing {partial}")
    os.replace(partial, path)


class FramedZstdFile:
    """A file written as zstd frames, each through a zstd command of its own (ZSTD_OPTIONS, up to 8 of `threads`), with
    a seek table at the end (seek_table): a frame ends where the writer calls new_frame (at whole records), so that
    protal --build reads the frames on several threads. Written as path.partial, renamed by close."""

    def __init__(self, path, threads):
        self.path, self.threads = path, threads
        self.fh = open(path + ".partial", "wb")
        self.process, self.frames, self.size, self.start = None, [], 0, 0

    def write(self, data):
        if self.process is None:
            self.start = self.fh.seek(0, os.SEEK_END)  # the commands write through the same file
            self.process = subprocess.Popen([shutil.which("zstd") or "zstd", "-q", "-c", f"-T{max(1, min(self.threads, 8))}",
                                             *ZSTD_OPTIONS], stdin=subprocess.PIPE, stdout=self.fh)
        self.process.stdin.write(data)
        self.size += len(data)

    def new_frame(self):
        if self.process is None:
            return
        self.process.stdin.close()
        rc = self.process.wait()
        if rc:
            sys.exit(f"zstd failed with exit code {rc} writing {self.path}.partial")
        self.frames.append((self.fh.seek(0, os.SEEK_END) - self.start, self.size))
        self.process, self.size = None, 0

    def close(self):
        self.new_frame()
        table = seek_table(self.frames)
        if table:
            self.fh.write(table)
        self.fh.close()
        os.replace(self.path + ".partial", self.path)


class PlainFile:
    """full_reference.fna without the zstd command: the frames FramedZstdFile would make are of no use."""

    def __init__(self, path):
        self.path, self.fh = path, open(path + ".partial", "wb")

    def write(self, data):
        self.fh.write(data)

    def new_frame(self):
        pass

    def close(self):
        self.fh.close()
        os.replace(self.path + ".partial", self.path)


@contextlib.contextmanager
def full_reference_writer(folder, threads):
    """A file to write a folder's full reference to (write(bytes), new_frame() between genes): full_reference.fna.zst
    in zstd frames with a seek table (FramedZstdFile), or full_reference.fna without the zstd command, written as a
    .partial file and renamed once complete; the variant of an earlier conversion goes first."""
    zstd = shutil.which("zstd")
    remove_full_reference(folder)
    target = os.path.join(folder, FULL_REFERENCE + (".zst" if zstd else ""))
    if not zstd:
        sys.stderr.write(f"Note: no zstd command, so {target} is written uncompressed\n")
    out = FramedZstdFile(target, threads) if zstd else PlainFile(target)
    yield out
    out.close()


@contextlib.contextmanager
def read_full_reference(path):
    """A full reference, plain or .zst (through the zstd command), as a binary file of its lines."""
    if not path.endswith(".zst"):
        with open(path, "rb") as fh:
            yield fh
        return
    zstd = shutil.which("zstd")
    if not zstd:
        sys.exit(f"Reading {path} needs the zstd command")
    process = subprocess.Popen([zstd, "-q", "-dc", "--long=31", path], stdout=subprocess.PIPE)
    try:
        yield process.stdout
    finally:
        process.stdout.close()
        rc = process.wait()
    if rc:
        sys.exit(f"zstd failed with exit code {rc} reading {path}")


def clear_build_outputs(folder):
    for name in os.listdir(folder) if os.path.isdir(folder) else ():
        path = os.path.join(folder, name)
        if (name in BUILD_OUTPUTS or name.endswith(".partial")) and os.path.isfile(path):
            os.remove(path)


def read_gene_ids(path):
    """gene2geneid.tsv of a converted folder -> {marker id: gene id}."""
    gene_ids = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2 and f[1].isdigit():
                gene_ids[f[0]] = int(f[1])
    return gene_ids


def read_gene_list(spec, gene_ids, what="--genes"):
    """The gene ids of a --genes value: GTDB marker ids (PF00380.20, TIGR00001; PF00380 without its version
    matches any) or protal gene ids (the numbers of gene_ids, {marker: gene id}), comma-separated or in a file
    (one per line, its first tab- or space-separated column, # comments). -> sorted gene ids; exits on an unknown
    id, an empty list or a marker id that matches several."""
    if os.path.isfile(spec):
        items = []
        with open(spec) as fh:
            for line in fh:
                token = line.strip().split("\t")[0].split(" ")[0]
                if token and not token.startswith("#"):
                    items.append(token)
    else:
        items = [t.strip() for t in spec.split(",") if t.strip()]
    if not items:
        sys.exit(f"{what}: no gene ids in {spec}")
    by_stem = {}
    for marker, gid in gene_ids.items():
        by_stem.setdefault(marker.split(".")[0], []).append((marker, gid))
    genes = set()
    for token in items:
        if token.isdigit():
            if int(token) not in gene_ids.values():
                sys.exit(f"{what}: gene id {token} is not in gene2geneid.tsv (ids 1-{max(gene_ids.values(), default=0)})")
            genes.add(int(token))
        elif token in gene_ids:
            genes.add(gene_ids[token])
        elif token in by_stem:
            if len(by_stem[token]) > 1:
                sys.exit(f"{what}: {token} matches several markers ({', '.join(m for m, _ in by_stem[token])}): name one")
            genes.add(by_stem[token][0][1])
        else:
            sys.exit(f"{what}: {token} is not a marker of gene2geneid.tsv (a GTDB marker id such as PF00380.20 or "
                     "TIGR00001, or a protal gene id)")
    return sorted(genes)


def derive_db(src, dst, names=(), genes=None, threads=1):
    """Copy the converted database folder src (before protal --build) to dst, leaving out the marker genes
    of the species `names` and, with `genes` (gene ids), every gene but those, from reference.fna, reference.map
    and the full reference. The taxonomy keeps the species with the same taxids: truth files still name them,
    and a model trained on dst applies to src; the genes keep their ids (gene2geneid.tsv lists the ones kept).
    The gene neighbours are counted anew from gene_positions.tsv over the species and genes kept. Only the
    files this script (and build_gtdb_database.py) wrote are copied, not what a build of src wrote or left
    there; dst's own build outputs are removed."""
    fna = os.path.join(src, "reference.fna")
    if not os.path.isfile(fna):
        sys.exit(f"{src} has no reference.fna: derive the copy before protal --build packs the database")
    with open(os.path.join(src, "internal_taxonomy.dmp")) as fh:
        next(fh)
        rows = [(f[0], f[3], f[4]) for f in (line.rstrip("\n").split("\t") for line in fh)]
    drop = species_taxids(rows, names, "--exclude_species") if names else set()
    keep = set(genes) if genes is not None else None
    os.makedirs(dst, exist_ok=True)
    clear_build_outputs(dst)
    for name in CONVERTED_FILES:
        if os.path.isfile(os.path.join(src, name)):
            shutil.copyfile(os.path.join(src, name), os.path.join(dst, name))
        elif os.path.isfile(os.path.join(dst, name)):
            os.remove(os.path.join(dst, name))
    if keep is not None and os.path.isfile(os.path.join(src, "gene2geneid.tsv")):
        with open(os.path.join(src, "gene2geneid.tsv")) as fin, \
                open(os.path.join(dst, "gene2geneid.tsv"), "w", newline="\n") as fout:
            for line in fin:
                f = line.rstrip("\n").split("\t")
                if len(f) < 2 or not f[1].isdigit() or int(f[1]) in keep:
                    fout.write(line)
    if os.path.isfile(os.path.join(src, GENE_POSITIONS)):
        # The gene neighbours of the genomes of the species kept, among the genes kept: the left-out species'
        # gene order is not known to the copy, as an organism the database lacks is not; and a read meets the
        # nearest gene of the subset next, not a gene the copy lacks.
        import gene_neighbours  # here: it imports this module
        lines, genomes, kept_species = gene_neighbours.derive(
            os.path.join(src, GENE_POSITIONS), os.path.join(dst, "internal_taxonomy.dmp"),
            os.path.join(dst, gene_neighbours.FILE_NAME), {int(t) for t in drop}, os.path.join(dst, GENE_POSITIONS), keep)
        sys.stderr.write(f"{gene_neighbours.FILE_NAME} derived from the {genomes} genomes of {kept_species} species "
                         f"kept: {lines} lines" + (f" (over the {len(keep)} genes kept)" if keep is not None else "") + "\n")
    elif os.path.isfile(os.path.join(dst, GENE_POSITIONS)):
        os.remove(os.path.join(dst, GENE_POSITIONS))

    def records(path):
        with open(path) as fh:
            for header in fh:
                seq = fh.readline()
                if not header.startswith(">") or not seq:
                    sys.exit(f"{path}: expected a header and one sequence line per record")
                tid, gid = header[1:].split()[0].split("_", 1)
                yield header, seq, tid, gid

    def wanted(tid, gid):
        return tid not in drop and (keep is None or int(gid) in keep)

    kept = dropped = offset = 0
    genes_kept = set()
    with open(os.path.join(dst, "reference.fna"), "w", newline="\n") as out, \
            open(os.path.join(dst, "reference.map"), "w", newline="\n") as fmap:
        for header, seq, tid, gid in records(fna):
            if not wanted(tid, gid):
                dropped += 1
                continue
            out.write(header + seq)
            start = offset + len(header)
            fmap.write(f"{tid}\t{gid}\t{start}\t{start + len(seq) - 1}\n")
            offset = start + len(seq)
            kept += 1
            genes_kept.add(int(gid))
    if keep is not None and genes_kept != keep:
        missing = sorted(keep - genes_kept)
        sys.exit(f"--genes: {fna} has no sequence of gene{'s' if len(missing) > 1 else ''} {', '.join(map(str, missing))}")
    full = full_reference_path(src)
    full_kept = 0
    if full:
        with read_full_reference(full) as fin, full_reference_writer(dst, threads) as out:
            frame_gene = None  # a frame per gene (the full reference is gene by gene), whose copies stay together
            for header in fin:
                seq = fin.readline()
                if not header.startswith(b">") or not seq:
                    sys.exit(f"{full}: expected a header and one sequence line per record")
                tid, gid = header[1:].split()[0].split(b"_", 1)  # the name; the genome's accession follows it
                if wanted(tid.decode(), gid.decode()):
                    if gid != frame_gene:
                        out.new_frame()
                        frame_gene = gid
                    out.write(header)
                    out.write(seq)
                    full_kept += 1
    else:
        remove_full_reference(dst)  # of an earlier copy; protal --build would take it
    sys.stderr.write(f"{src} -> {dst} without {len(drop)} species" +
                     (f", {len(keep)} genes kept" if keep is not None else "") +
                     f": {kept} representative sequences kept, {dropped} left out; full reference {full_kept} sequences\n")


def build_taxonomy(species_lineages):
    """Species (sorted by lineage) get ids 1..S, then root, then higher ranks.
    -> (rows for internal_taxonomy.dmp, {species name: taxid})"""
    ordered = sorted(species_lineages.items(), key=lambda kv: kv[1])
    ids, rows = {}, []
    for i, (species, _lin) in enumerate(ordered, 1):
        ids[species] = i
    root_id = len(ordered) + 1
    next_id = root_id + 1
    prefix_id = {}
    for level, _ in enumerate(RANKS[:-1], 1):
        for prefix in sorted({";".join(lin.split(";")[:level]) for _, lin in ordered}):
            prefix_id[prefix] = next_id
            next_id += 1

    rows.append((root_id, root_id, 0, "root", "no rank", 0, ""))
    for prefix, tid in sorted(prefix_id.items(), key=lambda kv: kv[1]):
        parts = prefix.split(";")
        parent = prefix_id.get(";".join(parts[:-1]), root_id)
        rows.append((tid, parent, 0, parts[-1], RANKS[len(parts) - 1][1], len(parts), ""))
    return rows, ids, prefix_id, root_id


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--gtdb", help="extracted GTDB release directory")
    ap.add_argument("--from_db", help="instead of --gtdb: a folder this script wrote (before protal --build), "
                                      "copied without the species of --exclude_species")
    ap.add_argument("--outdir", required=True, help="protal database directory to write")
    ap.add_argument("--release", help="release number, e.g. 226 (default: detected)")
    ap.add_argument("--model", default=os.path.join(SCRIPT_DIR, "..", "random_forest.xml"),
                    help="random forest PMML for paired-end reads, copied to <outdir>/model_pe.xml")
    ap.add_argument("--order", choices=("gene", "genome"), default="gene",
                    help="reference.fna record order: by gene, then taxid (compresses better), or by taxid, then gene")
    ap.add_argument("-t", "--threads", type=int, default=1,
                    help="marker files read in parallel (default 1); the output does not depend on it")
    ap.add_argument("--tmp", help="where the marker genes are spooled, in a .convert_tmp folder in it, removed when "
                                  "done (default: OUTDIR), e.g. a node's own disk")
    ap.add_argument("--exclude_species", help="file of species (s__Genus species, one per line) whose marker genes "
                                              "are left out; the taxonomy keeps them, with the taxids they have "
                                              "with all species (a training database with species held out)")
    ap.add_argument("--genes", help="keep only these marker genes (a reduced database): GTDB marker ids (PF00380.20, "
                                    "TIGR00001; PF00380 matches any version) or protal gene ids, comma-separated or "
                                    "one per line in a file (first column, # comments). The genes keep their ids; "
                                    "with --from_db the gene neighbours are counted anew over them")
    ap.add_argument("--priors_only", action="store_true",
                    help="write species_priors.tsv only (into --outdir, nothing else there is touched): the priors with "
                         "the genome sizes for a database converted from this release before 2026-10-08, to store with "
                         "protal --add_tables (the taxids are the same for a release, with or without --exclude_species "
                         "and --genes)")
    args = ap.parse_args()
    if bool(args.gtdb) == bool(args.from_db):
        ap.error("give --gtdb or --from_db")
    if args.priors_only and not args.gtdb:
        ap.error("--priors_only reads the release: give --gtdb")
    exclude = read_species_list(args.exclude_species) if args.exclude_species else set()
    if args.from_db:
        # Without --exclude_species and --genes: a plain copy of the converted files (a folder to build apart).
        genes = None
        if args.genes:
            table = os.path.join(args.from_db, "gene2geneid.tsv")
            if not os.path.isfile(table):
                sys.exit(f"--genes: {args.from_db} has no gene2geneid.tsv")
            genes = read_gene_list(args.genes, read_gene_ids(table))
        derive_db(args.from_db, args.outdir, exclude, genes, args.threads)
        return

    clock = [time.time()]

    def phase(what):
        """Logs how long the step that ended took (convert.log of build_gtdb_database.py)."""
        now = time.time()
        sys.stderr.write(f"{what}: {now - clock[0]:.1f} s\n")
        sys.stderr.flush()
        clock[0] = now

    rel = args.release or detect_release(args.gtdb)
    lineage = read_taxonomy(args.gtdb, rel)
    reps = read_representatives(args.gtdb, rel)
    phase(f"read the taxonomy and metadata ({len(lineage)} genomes)")

    rep_files = marker_files(args.gtdb, "genomic_files_reps", "reps", rel)
    if not rep_files:
        sys.exit(f"No *_marker_genes_reps_r{rel} marker files under {args.gtdb}/genomic_files_reps")
    all_files = marker_files(args.gtdb, "genomic_files_all", "all", rel)
    phase(f"found the marker files ({len(rep_files)} of representatives, {len(all_files)} of all genomes)")

    # Gene ids follow the order of first appearance (bac120 markers, then ar53-only ones).
    gene_ids = {}
    for marker, _ in rep_files + all_files:
        gene_ids.setdefault(marker, len(gene_ids) + 1)
    # --genes: the database holds these genes only, with the ids they have among all markers (a copy of the
    # full conversion derived with --from_db --genes gives the same files). The representatives' files of
    # every marker are still read: species_priors.tsv counts the markers found over the whole set.
    keep = set(read_gene_list(args.genes, gene_ids)) if args.genes else None
    if keep is not None:
        all_files = [(marker, path) for marker, path in all_files if gene_ids[marker] in keep]

    # The marker files are read in parallel (--threads), one marker at a time per worker, and spooled
    # to a temporary folder (--tmp): the representatives' genes are never all in memory.
    tmp = os.path.join(args.tmp or args.outdir, ".convert_tmp")
    shutil.rmtree(tmp, ignore_errors=True)
    os.makedirs(tmp)
    _WORK.update(lineage=lineage, reps=reps, tmp=tmp, gene_ids=gene_ids, zstd=bool(shutil.which("zstd")))

    # --- representative marker genes -------------------------------------------------
    skipped = {"not in taxonomy": 0, "not a representative": 0, "duplicate": 0}
    by_marker = {}
    for marker, path in rep_files:
        by_marker.setdefault(marker, []).append(path)
    species_rep = {}       # species name -> rep accession
    markers_of = {}        # rep accession -> markers found, markers with more than one copy, bases of the markers
    for marker, accessions, counts, duplicated in _parallel(args.threads, _spool_representatives, list(by_marker.items())):
        for key, n in counts.items():
            skipped[key] += n
        for acc, bases in accessions:
            species = lineage[acc].split(";")[-1]
            other = species_rep.setdefault(species, acc)
            if other != acc:
                sys.exit(f"Species {species} has two representatives with marker genes: {other}, {acc}")
            found = markers_of.setdefault(acc, [0, 0, 0])
            found[0] += 1
            found[2] += bases
        for acc in duplicated:
            markers_of.setdefault(acc, [0, 0, 0])[1] += 1
    phase(f"spooled the representatives' marker genes ({len(species_rep)} species)")
    quality, sizes = read_metadata(args.gtdb, rel)
    genome_sizes = species_genome_sizes(quality, sizes)
    clusters = read_sp_clusters(args.gtdb, rel)

    species_lineage = {sp: lineage[acc] for sp, acc in species_rep.items()}
    for sp, lin in species_lineage.items():
        if len(lin.split(";")) != len(RANKS):
            sys.exit(f"Lineage of {sp} does not have the 7 GTDB ranks: {lin}")
    rows, taxid, prefix_id, root_id = build_taxonomy(species_lineage)
    for sp, lin in species_lineage.items():
        parent = prefix_id[";".join(lin.split(";")[:-1])]
        rows.append((taxid[sp], parent, 0, sp, "species", len(RANKS), species_rep[sp]))
    rows.sort()
    drop = species_taxids([(str(taxid[sp]), sp, "species") for sp in taxid], exclude, "--exclude_species")
    drop = {int(t) for t in drop}

    os.makedirs(args.outdir, exist_ok=True)
    out = lambda name: os.path.join(args.outdir, name)
    if not args.priors_only:
        clear_build_outputs(args.outdir)
        with open(out("internal_taxonomy.dmp"), "w", newline="\n") as fh:
            fh.write("id\tparent_id\texternal_id\tname\trank\tlevel\trep_genome\n")
            for row in rows:
                fh.write("\t".join(map(str, row)) + "\n")

    # species_priors.tsv: what GTDB knows of each species' representative and cluster before any read, for the
    # model's prior features (protal's SpeciesPriors.h; -1: unknown). Duplicated single-copy markers in the
    # representative are CheckM's contamination signature; a contaminating contig's genes put every present
    # organism's reads on the species (docs/claude/2026-10-03-false-positive-anatomy). Since 2026-10-08 also the
    # species' genome size (species_genome_sizes), the representative's assembly size and the bases of its marker genes
    # (of every marker, as `markers`: a reduced database holds fewer), and their share of the assembly, for the profile's
    # composition (protal's Composition.h: the reads the called species explain, the sample's unknown share).
    def number(v):
        return "-1" if v is None else (f"{v:g}" if isinstance(v, float) else str(v))
    with_dups = with_quality = with_cluster = with_size = 0
    with open(out("species_priors.tsv"), "w", newline="\n") as fh:
        fh.write("taxid\trep_genome\tmarkers\tduplicate_markers\tcheckm_completeness\tcheckm_contamination"
                 "\tani_radius\tmean_intra_ani\tmin_intra_ani\tclustered_genomes"
                 "\tgenome_size\tsized_genomes\trep_genome_size\tmarker_bases\tmarker_share\n")
        for sp in sorted(taxid, key=lambda s: taxid[s]):
            acc = species_rep[sp]
            found, dups, marker_bases = markers_of.get(acc, [0, 0, 0])
            comp, cont = quality.get(acc, (None, None))
            cl = clusters.get(acc, {})
            size, sized = genome_sizes.get(acc, (None, 0))
            rep_size = sizes[acc][0] if acc in sizes else None
            with_dups += dups > 0
            with_quality += comp is not None
            with_cluster += bool(cl)
            with_size += size is not None
            fh.write("\t".join([str(taxid[sp]), acc, str(found), str(dups), number(comp), number(cont), number(cl.get("radius")),
                                 number(cl.get("mean_ani")), number(cl.get("min_ani")), number(cl.get("genomes")),
                                 "-1" if size is None else str(round(size)), str(sized),
                                 "-1" if rep_size is None else str(round(rep_size)), str(marker_bases),
                                 f"{marker_bases / rep_size:.4g}" if rep_size else "-1"]) + "\n")
    phase(f"species priors ({with_dups} representatives with a duplicated marker, CheckM quality for {with_quality}, "
          f"species clusters for {with_cluster}, genome sizes for {with_size} of {len(taxid)} species"
          + ("" if clusters else "; no sp_clusters file") + ("" if sizes else "; no genome_size in the metadata") + ")")
    if args.priors_only:
        shutil.rmtree(tmp, ignore_errors=True)
        sys.stderr.write(f"GTDB r{rel} -> {out('species_priors.tsv')} ({len(taxid)} species; store it in a database built "
                         f"from this release with protal --add_tables {out('species_priors.tsv')} --db DB)\n")
        return

    # Each gene's records, sorted by taxid, then all genes in gene id order (--order gene), or merged
    # by taxid, then gene (--order genome).
    _WORK.update(taxid=taxid, drop=drop)
    markers = sorted((m for m in by_marker if keep is None or gene_ids[m] in keep), key=lambda m: gene_ids[m])
    if keep is not None and len(markers) < len(keep):
        sys.exit(f"--genes: no representative has gene{'s' if len(keep) - len(markers) > 1 else ''} "
                 f"{', '.join(str(g) for g in sorted(keep) if g not in {gene_ids[m] for m in markers})}")
    phase("wrote the taxonomy")
    lengths = dict(zip(markers, _parallel(args.threads, _write_gene, markers)))
    n_records = sum(len(v) for v in lengths.values())
    phase(f"sorted each gene's representative sequences ({n_records})")
    with open(out("reference.fna"), "wb") as fna, open(out("reference.map"), "w", newline="\n") as fmap:
        offset = 0

        def put(tid, gid, seq_length):
            nonlocal offset
            offset += len(f">{tid}_{gid}\n")
            fmap.write(f"{tid}\t{gid}\t{offset}\t{offset + seq_length}\n")
            offset += seq_length + 1

        if args.order == "gene":
            for marker in markers:
                with open(os.path.join(tmp, f"gene_{gene_ids[marker]}.fna"), "rb") as chunk:
                    shutil.copyfileobj(chunk, fna, 1 << 22)
                for tid, seq_length in lengths[marker]:
                    put(tid, gene_ids[marker], seq_length)
        else:
            streams = [_gene_records(os.path.join(tmp, f"gene_{gene_ids[m]}.fna"), gene_ids[m]) for m in markers]
            for tid, gid, seq in heapq.merge(*streams):
                fna.write(f">{tid}_{gid}\n".encode() + seq + b"\n")
                put(tid, gid, len(seq))

    with open(out("gene2geneid.tsv"), "w", newline="\n") as fh:
        for marker, gid in gene_ids.items():
            if keep is None or gid in keep:
                fh.write(f"{marker}\t{gid}\n")

    with open(out("genome2tiid.tsv"), "w", newline="\n") as fh:
        for acc in sorted(lineage):
            sp = lineage[acc].split(";")[-1]
            if sp in taxid:
                fh.write(f"{acc}\t{taxid[sp]}\t{species_rep[sp]}\t{lineage[acc]}\n")
    phase("wrote reference.fna, reference.map, gene2geneid.tsv and genome2tiid.tsv")

    # --- all genomes, for the unique k-mer check -----------------------------------------
    n_full = 0
    remove_full_reference(args.outdir)  # of an earlier conversion; protal --build would take it
    if all_files:
        # One worker per marker (duplicates are per genome and marker), one chunk per file, compressed by the
        # worker when there is zstd (each chunk a zstd frame: joined in file order they are the compressed file,
        # so the ~90 GB of GTDB r226 are never written or read uncompressed), joined in file order.
        all_by_marker = {}
        for index, (marker, path) in enumerate(all_files):
            all_by_marker.setdefault(marker, []).append((index, path))
        written = {}
        for chunks in _parallel(args.threads, _write_full_reference, list(all_by_marker.items())):
            written.update(chunks)
        phase(f"wrote the full reference's chunks ({len(all_files)}" + (", zstd-compressed)" if _WORK["zstd"] else ")"))
        target = out(FULL_REFERENCE + (".zst" if _WORK["zstd"] else ""))
        if not _WORK["zstd"]:
            sys.stderr.write(f"Note: no zstd command, so {target} is written uncompressed\n")
        frames = []  # (compressed size, content size) of each chunk, a zstd frame each
        with open(target + ".partial", "wb") as fh:
            for index in range(len(all_files)):
                chunk, count, size = written[index]
                frames.append((os.path.getsize(chunk), size))
                with open(chunk, "rb") as part:
                    shutil.copyfileobj(part, fh, 1 << 22)
                os.remove(chunk)  # not kept until all are joined
                n_full += count
            table = seek_table(frames) if _WORK["zstd"] else None
            if table:  # protal --build decompresses the frames on several threads
                fh.write(table)
        os.replace(target + ".partial", target)
        phase(f"joined them into {os.path.basename(target)}")
    shutil.rmtree(tmp, ignore_errors=True)
    full = full_reference_path(args.outdir)

    if os.path.exists(args.model):
        shutil.copyfile(args.model, out("model_pe.xml"))
    else:
        sys.stderr.write(f"Warning: model {args.model} not found; add model_pe.xml before profiling\n")

    sys.stderr.write(
        f"GTDB r{rel} -> {args.outdir}\n"
        f"  species:            {len(taxid)} (taxids 1..{len(taxid)}, root {root_id})"
        + (f", {len(drop)} without marker genes (--exclude_species)" if drop else "") + "\n"
        f"  marker genes:       {len(markers)} ids" + (f" of {len(gene_ids)} (--genes)" if keep is not None else "")
        + f", {n_records} representative sequences\n"
        f"  full reference:     {n_full} sequences" +
        (f" ({os.path.basename(full)}, {os.path.getsize(full) / 1e9:.2f} GB)" if full else " (no genomic_files_all)") + "\n"
        f"  skipped rep records: " + ", ".join(f"{k}: {v}" for k, v in skipped.items()) + "\n")


if __name__ == "__main__":
    main()
