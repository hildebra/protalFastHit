#!/usr/bin/env python3
"""gtdb_to_protal_db.py - turn an (extracted) GTDB release into the input files
of a protal database folder, ready for `protal --build`.

Reads from the GTDB release directory:
  {bac120,ar53}_taxonomy_r<R>.tsv[.gz]     accession -> GTDB lineage
  {bac120,ar53}_metadata_r<R>.tsv[.gz]     species representatives (optional;
                                           without it, genomes present in the
                                           *_marker_genes_reps_* files are the reps)
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
                         zstd-compressed by the zstd command (86 GB raw at r226),
                         full_reference.fna without one
  gene2geneid.tsv        marker id -> geneid
  genome2tiid.tsv        accession, species taxid, species rep accession, lineage
  model_pe.xml           copy of --model (the profiler's random forest for paired-end reads;
                         models for other read types: protal --add_model, see README)

Then build the index with
  protal --build --no_profile --db <outdir> --reference <outdir>/reference.fna \\
         --full_reference <outdir>/full_reference.fna
which writes index.prx.zst and replaces reference.fna by reference.fna.zst unless --no_compress
(for --full_reference full_reference.fna, protal reads its .zst).

Each marker's genes are spooled to <outdir>/.convert_tmp and sorted one gene at a time, so memory
stays at about one gene's sequences per worker; -t reads the marker files in parallel (the output is
the same for any -t).

Usage:
  gtdb_to_protal_db.py --gtdb <release dir> --outdir <db dir> [--release 226] [--model FILE]
      [--exclude_species FILE]
  gtdb_to_protal_db.py --from_db <db dir> --exclude_species FILE --outdir <training db dir>

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
import subprocess
import sys

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
    to rep_<gene id>.tsv; -> (marker, their accessions, skipped counts)."""
    marker, paths = item
    lineage, reps = _WORK["lineage"], _WORK["reps"]
    counts = {"not in taxonomy": 0, "not a representative": 0, "duplicate": 0}
    kept, accessions = set(), []
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
                else:
                    kept.add(acc)
                    accessions.append(acc)
                    out.write(f"{acc}\t{seq.upper()}\n")
    return marker, accessions, counts


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
    genome's first copy, one chunk per file; -> {file index: (chunk, records)}."""
    marker, files = item
    gid = _WORK["gene_ids"][marker]
    lineage, taxid, drop = _WORK["lineage"], _WORK["taxid"], _WORK["drop"]
    seen, result = set(), {}
    for index, path in files:
        chunk, n = os.path.join(_WORK["tmp"], f"full_{index}.fna"), 0
        with open(chunk, "w", newline="\n") as fh:
            for header, seq in read_fasta(path):
                acc = normalize_accession(header)
                sp = lineage.get(acc, "").split(";")[-1]
                if sp not in taxid or taxid[sp] in drop or acc in seen:
                    continue
                seen.add(acc)
                fh.write(f">{taxid[sp]}_{gid}\n{seq.upper()}\n")
                n += 1
        result[index] = (chunk, n)
    return result


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
# and gene_neighbours.tsv (gene_neighbours.py; per clade, so it holds for a copy without some species).
CONVERTED_FILES = ("internal_taxonomy.dmp", "gene2geneid.tsv", "genome2tiid.tsv", "gene_neighbours.tsv",
                   "model_pe.xml", "model_se.xml", "model_PB.xml", "model_ONT.xml")
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
def full_reference_writer(folder, threads):
    """A binary file to write a folder's full reference to: full_reference.fna.zst through the zstd command
    (ZSTD_OPTIONS, up to 8 of `threads`), or full_reference.fna without one. Written as a .partial file,
    renamed once complete; the variant of an earlier conversion goes first."""
    zstd = shutil.which("zstd")
    remove_full_reference(folder)
    target = os.path.join(folder, FULL_REFERENCE + (".zst" if zstd else ""))
    partial = target + ".partial"
    if not zstd:
        sys.stderr.write(f"Note: no zstd command, so {target} is written uncompressed\n")
        with open(partial, "wb") as fh:
            yield fh
    else:
        process = subprocess.Popen([zstd, "-q", "-f", f"-T{max(1, min(threads, 8))}", *ZSTD_OPTIONS, "-o", partial],
                                   stdin=subprocess.PIPE)
        try:
            yield process.stdin
        finally:
            process.stdin.close()
            rc = process.wait()
        if rc:
            sys.exit(f"zstd failed with exit code {rc} writing {partial}")
    os.replace(partial, target)


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


def exclude_from_db(src, dst, names, threads=1):
    """Copy the converted database folder src (before protal --build) to dst, leaving out the marker genes
    of the species `names` from reference.fna, reference.map and the full reference. The taxonomy keeps
    them with the same taxids: truth files still name them, and a model trained on dst applies to src.
    Only the files this script (and build_gtdb_database.py) wrote are copied, not what a build of src
    wrote or left there; dst's own build outputs are removed."""
    fna = os.path.join(src, "reference.fna")
    if not os.path.isfile(fna):
        sys.exit(f"{src} has no reference.fna: exclude species before protal --build packs the database")
    with open(os.path.join(src, "internal_taxonomy.dmp")) as fh:
        next(fh)
        rows = [(f[0], f[3], f[4]) for f in (line.rstrip("\n").split("\t") for line in fh)]
    drop = species_taxids(rows, names, "--exclude_species")
    os.makedirs(dst, exist_ok=True)
    clear_build_outputs(dst)
    for name in CONVERTED_FILES:
        if os.path.isfile(os.path.join(src, name)):
            shutil.copyfile(os.path.join(src, name), os.path.join(dst, name))
        elif os.path.isfile(os.path.join(dst, name)):
            os.remove(os.path.join(dst, name))

    def records(path):
        with open(path) as fh:
            for header in fh:
                seq = fh.readline()
                if not header.startswith(">") or not seq:
                    sys.exit(f"{path}: expected a header and one sequence line per record")
                yield header, seq, header[1:].split("_", 1)[0]

    kept = dropped = offset = 0
    with open(os.path.join(dst, "reference.fna"), "w", newline="\n") as out, \
            open(os.path.join(dst, "reference.map"), "w", newline="\n") as fmap:
        for header, seq, tid in records(fna):
            if tid in drop:
                dropped += 1
                continue
            gid = header[1:].rstrip("\n").split("_", 1)[1]
            out.write(header + seq)
            start = offset + len(header)
            fmap.write(f"{tid}\t{gid}\t{start}\t{start + len(seq) - 1}\n")
            offset = start + len(seq)
            kept += 1
    full = full_reference_path(src)
    full_kept = 0
    if full:
        drop_bytes = {tid.encode() for tid in drop}
        with read_full_reference(full) as fin, full_reference_writer(dst, threads) as out:
            for header in fin:
                seq = fin.readline()
                if not header.startswith(b">") or not seq:
                    sys.exit(f"{full}: expected a header and one sequence line per record")
                if header[1:].split(b"_", 1)[0] not in drop_bytes:
                    out.write(header)
                    out.write(seq)
                    full_kept += 1
    else:
        remove_full_reference(dst)  # of an earlier copy; protal --build would take it
    sys.stderr.write(f"{src} -> {dst} without {len(drop)} species: {kept} representative sequences kept, "
                     f"{dropped} left out; full reference {full_kept} sequences\n")


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
    ap.add_argument("--exclude_species", help="file of species (s__Genus species, one per line) whose marker genes "
                                              "are left out; the taxonomy keeps them, with the taxids they have "
                                              "with all species (a training database with species held out)")
    args = ap.parse_args()
    if bool(args.gtdb) == bool(args.from_db):
        ap.error("give --gtdb or --from_db")
    exclude = read_species_list(args.exclude_species) if args.exclude_species else set()
    if args.from_db:
        if not exclude:
            ap.error("--from_db needs --exclude_species")
        exclude_from_db(args.from_db, args.outdir, exclude, args.threads)
        return

    rel = args.release or detect_release(args.gtdb)
    lineage = read_taxonomy(args.gtdb, rel)
    reps = read_representatives(args.gtdb, rel)

    rep_files = marker_files(args.gtdb, "genomic_files_reps", "reps", rel)
    if not rep_files:
        sys.exit(f"No *_marker_genes_reps_r{rel} marker files under {args.gtdb}/genomic_files_reps")
    all_files = marker_files(args.gtdb, "genomic_files_all", "all", rel)

    # Gene ids follow the order of first appearance (bac120 markers, then ar53-only ones).
    gene_ids = {}
    for marker, _ in rep_files + all_files:
        gene_ids.setdefault(marker, len(gene_ids) + 1)

    # The marker files are read in parallel (--threads), one marker at a time per worker, and spooled
    # to a temporary folder: the representatives' genes are never all in memory.
    tmp = os.path.join(args.outdir, ".convert_tmp")
    shutil.rmtree(tmp, ignore_errors=True)
    os.makedirs(tmp)
    _WORK.update(lineage=lineage, reps=reps, tmp=tmp, gene_ids=gene_ids)

    # --- representative marker genes -------------------------------------------------
    skipped = {"not in taxonomy": 0, "not a representative": 0, "duplicate": 0}
    by_marker = {}
    for marker, path in rep_files:
        by_marker.setdefault(marker, []).append(path)
    species_rep = {}       # species name -> rep accession
    for marker, accessions, counts in _parallel(args.threads, _spool_representatives, list(by_marker.items())):
        for key, n in counts.items():
            skipped[key] += n
        for acc in accessions:
            species = lineage[acc].split(";")[-1]
            other = species_rep.setdefault(species, acc)
            if other != acc:
                sys.exit(f"Species {species} has two representatives with marker genes: {other}, {acc}")

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
    clear_build_outputs(args.outdir)
    out = lambda name: os.path.join(args.outdir, name)

    with open(out("internal_taxonomy.dmp"), "w", newline="\n") as fh:
        fh.write("id\tparent_id\texternal_id\tname\trank\tlevel\trep_genome\n")
        for row in rows:
            fh.write("\t".join(map(str, row)) + "\n")

    # Each gene's records, sorted by taxid, then all genes in gene id order (--order gene), or merged
    # by taxid, then gene (--order genome).
    _WORK.update(taxid=taxid, drop=drop)
    markers = sorted(by_marker, key=lambda m: gene_ids[m])
    lengths = dict(zip(markers, _parallel(args.threads, _write_gene, markers)))
    n_records = sum(len(v) for v in lengths.values())
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
            fh.write(f"{marker}\t{gid}\n")

    with open(out("genome2tiid.tsv"), "w", newline="\n") as fh:
        for acc in sorted(lineage):
            sp = lineage[acc].split(";")[-1]
            if sp in taxid:
                fh.write(f"{acc}\t{taxid[sp]}\t{species_rep[sp]}\t{lineage[acc]}\n")

    # --- all genomes, for the unique k-mer check -----------------------------------------
    n_full = 0
    if all_files:
        # One worker per marker (duplicates are per genome and marker), one chunk per file, joined in
        # file order.
        all_by_marker = {}
        for index, (marker, path) in enumerate(all_files):
            all_by_marker.setdefault(marker, []).append((index, path))
        written = {}
        for chunks in _parallel(args.threads, _write_full_reference, list(all_by_marker.items())):
            written.update(chunks)
        with full_reference_writer(args.outdir, args.threads) as fh:
            for index in range(len(all_files)):
                chunk, count = written[index]
                with open(chunk, "rb") as part:
                    shutil.copyfileobj(part, fh, 1 << 22)
                os.remove(chunk)  # at GTDB size the chunks hold ~90 GB: not kept until all are joined
                n_full += count
    else:
        remove_full_reference(args.outdir)  # of an earlier conversion; protal --build would take it
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
        f"  marker genes:       {len(gene_ids)} ids, {n_records} representative sequences\n"
        f"  full reference:     {n_full} sequences" +
        (f" ({os.path.basename(full)}, {os.path.getsize(full) / 1e9:.2f} GB)" if full else " (no genomic_files_all)") + "\n"
        f"  skipped rep records: " + ", ".join(f"{k}: {v}" for k, v in skipped.items()) + "\n")


if __name__ == "__main__":
    main()
