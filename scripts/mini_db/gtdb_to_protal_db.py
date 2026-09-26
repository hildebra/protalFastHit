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
                         sequence line per record
  reference.map          taxid, geneid, start byte, end byte of each sequence line
  internal_taxonomy.dmp  id, parent_id, external_id, name, rank, level, rep_genome
                         (species are the leaves; their ids are the reference taxids)
  full_reference.fna     marker genes of all genomes, header >taxid_geneid of the
                         genome's species (only if genomic_files_all is present)
  gene2geneid.tsv        marker id -> geneid
  genome2tiid.tsv        accession, species taxid, species rep accession, lineage
  model.xml              copy of --model (the profiler's random forest)

Then build the index with
  protal --build --no_profile --db <outdir> --reference <outdir>/reference.fna \\
         --full_reference <outdir>/full_reference.fna
which writes index.prx.zst and replaces reference.fna by reference.fna.zst unless --no_compress.

Intended for small/sparse releases: representative sequences are held in memory.

Usage:
  gtdb_to_protal_db.py --gtdb <release dir> --outdir <db dir> [--release 226] [--model FILE]
"""

import argparse
import glob
import gzip
import os
import re
import shutil
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
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
    ap.add_argument("--gtdb", required=True, help="extracted GTDB release directory")
    ap.add_argument("--outdir", required=True, help="protal database directory to write")
    ap.add_argument("--release", help="release number, e.g. 226 (default: detected)")
    ap.add_argument("--model", default=os.path.join(SCRIPT_DIR, "..", "random_forest.xml"),
                    help="random forest PMML copied to <outdir>/model.xml")
    args = ap.parse_args()

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

    # --- representative marker genes -------------------------------------------------
    rep_seqs = {}          # (accession, marker) -> sequence
    skipped = {"not in taxonomy": 0, "not a representative": 0, "duplicate": 0}
    for marker, path in rep_files:
        for header, seq in read_fasta(path):
            acc = normalize_accession(header)
            if acc not in lineage:
                skipped["not in taxonomy"] += 1
            elif reps is not None and acc not in reps:
                skipped["not a representative"] += 1
            elif (acc, marker) in rep_seqs:
                skipped["duplicate"] += 1
            else:
                rep_seqs[(acc, marker)] = seq.upper()

    species_rep = {}       # species name -> rep accession
    for acc, _marker in rep_seqs:
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

    os.makedirs(args.outdir, exist_ok=True)
    out = lambda name: os.path.join(args.outdir, name)

    with open(out("internal_taxonomy.dmp"), "w", newline="\n") as fh:
        fh.write("id\tparent_id\texternal_id\tname\trank\tlevel\trep_genome\n")
        for row in rows:
            fh.write("\t".join(map(str, row)) + "\n")

    records = sorted(((taxid[lineage[acc].split(";")[-1]], gene_ids[marker], seq)
                      for (acc, marker), seq in rep_seqs.items()))
    with open(out("reference.fna"), "wb") as fna, open(out("reference.map"), "w", newline="\n") as fmap:
        offset = 0
        for tid, gid, seq in records:
            header = f">{tid}_{gid}\n".encode()
            offset += len(header)
            fna.write(header + seq.encode() + b"\n")
            fmap.write(f"{tid}\t{gid}\t{offset}\t{offset + len(seq)}\n")
            offset += len(seq) + 1

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
        seen = set()
        with open(out("full_reference.fna"), "w", newline="\n") as fh:
            for marker, path in all_files:
                for header, seq in read_fasta(path):
                    acc = normalize_accession(header)
                    sp = lineage.get(acc, "").split(";")[-1]
                    if sp not in taxid or (acc, marker) in seen:
                        continue
                    seen.add((acc, marker))
                    fh.write(f">{taxid[sp]}_{gene_ids[marker]}\n{seq.upper()}\n")
                    n_full += 1

    if os.path.exists(args.model):
        shutil.copyfile(args.model, out("model.xml"))
    else:
        sys.stderr.write(f"Warning: model {args.model} not found; add model.xml before profiling\n")

    sys.stderr.write(
        f"GTDB r{rel} -> {args.outdir}\n"
        f"  species:            {len(taxid)} (taxids 1..{len(taxid)}, root {root_id})\n"
        f"  marker genes:       {len(gene_ids)} ids, {len(records)} representative sequences\n"
        f"  full reference:     {n_full} sequences" + ("" if all_files else " (no genomic_files_all)") + "\n"
        f"  skipped rep records: " + ", ".join(f"{k}: {v}" for k, v in skipped.items()) + "\n")


if __name__ == "__main__":
    main()
