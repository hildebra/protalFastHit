#!/usr/bin/env python3
"""Where the marker reads of long-read samples go in protal: for the rowless-ONT-species study
(docs/claude/2026-10-15-ont-rowless-species).

The reads must be named with their origin, g<i>x_<n>|<contig>|<start>|<end>|<errors>|<strand> (simulate_metagenomes
with the study's local patch read_origin_names.patch and PROTAL_SIM_READ_ORIGIN=1; the same reads as without it;
start/end 0-based half-open on the contig, errors the simulator's substitutions, insertions and deleted bases in the
read). For each sample:
  - each read's overlap with the marker genes of its genome (gene_positions.tsv; an in-silico strain has its
    representative's genes on contigs named insilico_<accession>_<contig>): the most bases it shares with one gene,
    and whether it holds that gene whole; a marker read shares --min-overlap bases or more (the simulator's rule);
  - its fate in protal's SAM (written with --unmapped_reads, so a read that seeded but aligned nowhere has an unmapped
    record): none (no record: it seeded on nothing, or the gate dropped it), unaligned (an unmapped record), own (its
    primary record on its own species, MAPQ >= 4, as the profiler counts it), own_lowmapq (own species, MAPQ < 4), or
    elsewhere (its primary record on another species);
  - the species with a row: those protal profiled (the sample's .profile lines, or --rows: a file of taxids).
Writes OUT.reads.tsv.gz (every marker read), OUT.species.tsv (every species of the sample with marker reads: its reads
by fate, and whether it has a row) and prints a summary.

Usage: trace_marker_reads.py --reads R.fq.zst --sam S.sam[.gz|.zst] --genomes genomes.tsv --sample NAME
                             --gene-positions gene_positions.tsv --rows ROWS
                             [--heldout heldout_species.txt --taxonomy internal_taxonomy.dmp] --out PREFIX
"""
import argparse
import bisect
import collections
import gzip
import io
import os
import subprocess
import sys

MIN_MAPQ = 4  # Profiler::m_min_mapq


def open_text(path):
    if path.endswith(".zst"):
        p = subprocess.Popen(["zstd", "-dcq", path], stdout=subprocess.PIPE)
        return io.TextIOWrapper(p.stdout, encoding="ascii", errors="replace")
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def representative(genome):
    """insilico_GCA_000018565_1 -> GCA_000018565.1; others unchanged."""
    if not genome.startswith("insilico_"):
        return genome
    acc = genome[len("insilico_"):]
    head, _, version = acc.rpartition("_")
    return f"{head}.{version}"


def load_genes(path):
    """{accession: {contig: (starts, ends, genes, max_ends)}} 0-based half-open, and {accession: taxid}."""
    rows = collections.defaultdict(lambda: collections.defaultdict(list))
    taxid_of = {}
    with open(path) as fh:
        header = None
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if header is None:
                header = {k: i for i, k in enumerate(f)}
                continue
            acc, contig = f[header["accession"]], f[header["contig"]]
            start, end = int(f[header["start"]]) - 1, int(f[header["end"]])
            rows[acc][contig].append((start, end, f[header["gene"]]))
            taxid_of[acc] = f[header["taxid"]]
    genes = {}
    for acc, contigs in rows.items():
        genes[acc] = {}
        for contig, items in contigs.items():
            items.sort()
            starts = [s for s, _, _ in items]
            ends = [e for _, e, _ in items]
            names = [g for _, _, g in items]
            genes[acc][contig] = (starts, ends, names)
    return genes, taxid_of


def overlap(genes, genome, contig, start, end):
    """(the most bases the read shares with one marker gene, that gene, whether the read holds it whole, the gene's
    length) on its genome's contig."""
    acc = representative(genome)
    contigs = genes.get(acc)
    if not contigs:
        return 0, "", False, 0
    if genome.startswith("insilico_"):
        prefix = genome + "_"
        if contig.startswith(prefix):
            contig = contig[len(prefix):]
    entry = contigs.get(contig)
    if not entry:
        return 0, "", False, 0
    starts, ends, names = entry
    best = (0, "", False, 0)
    # genes are short (<= a few kb): scan those starting before the read's end, back to a few kb before its start
    i = bisect.bisect_left(starts, end)
    j = i - 1
    while j >= 0 and starts[j] > start - 20000:
        s, e = starts[j], ends[j]
        o = min(end, e) - max(start, s)
        if o > best[0]:
            best = (o, names[j], s >= start and e <= end, e - s)
        j -= 1
    return best


def read_genomes(path, sample):
    """The sample's genomes in the simulator's order (g<i>): genomes.tsv rows of the sample, host last."""
    out = []
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        h = {k: i for i, k in enumerate(head)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[h["sample"]] == sample:
                out.append(f[h["genome"]])
    return out


def cigar_ref_span(cigar):
    span, num = 0, ""
    for c in cigar:
        if c.isdigit():
            num += c
        else:
            if c in "MDN=X":
                span += int(num)
            num = ""
    return span


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reads", required=True)
    ap.add_argument("--sam", required=True)
    ap.add_argument("--genomes", required=True, help="the simulator's --long_genomes table")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--gene-positions", required=True)
    ap.add_argument("--rows", required=True, help="the taxids with a row (one per line)")
    ap.add_argument("--heldout", help="the species held out of the profiling database (model_logs/heldout_species.txt: "
                                      "names in the first column)")
    ap.add_argument("--taxonomy", help="the build's internal_taxonomy.dmp (id and name), to name the held-out species")
    ap.add_argument("--min-overlap", type=int, default=50)
    ap.add_argument("--out", required=True)
    o = ap.parse_args()

    genes, taxid_of = load_genes(o.gene_positions)
    genomes = read_genomes(o.genomes, o.sample)
    with open(o.rows) as fh:
        rows = {line.strip() for line in fh if line.strip()}
    heldout = set()
    if o.heldout and o.taxonomy:
        with open(o.heldout) as fh:
            names = {line.split("\t", 1)[0].strip() for line in fh if line.strip()}
        with open(o.taxonomy) as fh:
            head = fh.readline().rstrip("\n").split("\t")
            h = {k: i for i, k in enumerate(head)}
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if f[h["name"]] in names:
                    heldout.add(f[h["id"]])

    # marker reads from the FASTQ names (origin) and their qualities
    marker = {}
    n_reads = 0
    with open_text(o.reads) as fh:
        while True:
            name = fh.readline()
            if not name:
                break
            seq = fh.readline().rstrip("\n")
            fh.readline()
            qual = fh.readline().rstrip("\n")
            n_reads += 1
            f = name[1:].rstrip("\n").split("|")
            if len(f) < 5:
                continue
            gi = int(f[0][1:f[0].index("x")])
            genome = genomes[gi] if gi < len(genomes) else "host"
            contig, start, end = f[1], int(f[2]), int(f[3])
            errors = int(f[4]) if len(f) >= 6 else -1  # the simulator's errors in the read (patch v2)
            ov, gene, whole, glen = overlap(genes, genome, contig, start, end)
            if ov < min(o.min_overlap, end - start):
                continue
            meanq = sum(ord(c) - 33 for c in qual) / max(len(qual), 1)
            marker[name[1:].rstrip("\n")] = {
                "genome": genome, "taxid": taxid_of.get(representative(genome), "?"), "length": len(seq),
                "overlap": ov, "gene": gene, "whole": whole, "gene_length": glen, "meanq": meanq,
                "error_rate": errors / max(end - start, 1) if errors >= 0 else -1,
                "fate": "none", "mapq": -1, "on": "", "span": 0, "nm": -1, "records": 0}

    # their records
    with open_text(o.sam) as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            f = line.rstrip("\n").split("\t")
            r = marker.get(f[0])
            if r is None:
                continue
            r["records"] += 1
            flag = int(f[1])
            if flag & 4:
                if r["fate"] == "none":
                    r["fate"] = "unaligned"
                continue
            if flag & 0x904:
                continue
            taxid = f[2].split("_", 1)[0]
            mapq = int(f[4])
            own = taxid == r["taxid"]
            r["fate"] = ("own" if mapq >= MIN_MAPQ else "own_lowmapq") if own else "elsewhere"
            r["mapq"], r["on"], r["span"] = mapq, f[2], cigar_ref_span(f[5])
            for t in f[11:]:
                if t.startswith("NM:i:"):
                    r["nm"] = int(t[5:])

    with gzip.open(o.out + ".reads.tsv.gz", "wt") as out:
        cols = ["sample", "read", "genome", "taxid", "heldout", "length", "meanq", "error_rate", "gene", "gene_length",
                "overlap", "whole", "fate", "mapq", "on", "span", "nm", "records", "row"]
        out.write("\t".join(cols) + "\n")
        for name, r in marker.items():
            out.write("\t".join(str(x) for x in [o.sample, name, r["genome"], r["taxid"], int(r["taxid"] in heldout),
                                                 r["length"], f"{r['meanq']:.1f}", f"{r['error_rate']:.4f}",
                                                 r["gene"], r["gene_length"],
                                                 r["overlap"], int(r["whole"]), r["fate"], r["mapq"], r["on"],
                                                 r["span"], r["nm"], r["records"], int(r["taxid"] in rows)]) + "\n")

    by_species = collections.defaultdict(collections.Counter)
    for r in marker.values():
        c = by_species[r["taxid"]]
        c["marker_reads"] += 1
        c[r["fate"]] += 1
        c["whole" if r["whole"] else "partial"] += 1
    fates = ["own", "own_lowmapq", "elsewhere", "unaligned", "none"]
    with open(o.out + ".species.tsv", "w") as out:
        out.write("\t".join(["sample", "taxid", "heldout", "row", "marker_reads", "whole", "partial"] + fates) + "\n")
        for taxid, c in sorted(by_species.items(), key=lambda kv: -kv[1]["marker_reads"]):
            out.write("\t".join([o.sample, taxid, str(int(taxid in heldout)), str(int(taxid in rows)),
                                 str(c["marker_reads"]), str(c["whole"]), str(c["partial"])] +
                                [str(c[f]) for f in fates]) + "\n")

    in_db = [t for t in by_species if t not in heldout]
    rowless = [t for t in in_db if t not in rows]
    print(f"{o.sample}: {n_reads} reads, {len(marker)} marker reads of {len(by_species)} species "
          f"({len(by_species) - len(in_db)} held out); {len(rowless)} database species with marker reads have no row")
    total = collections.Counter(r["fate"] for r in marker.values() if r["taxid"] not in heldout)
    print("  database species' marker reads by fate: " + ", ".join(f"{k} {total[k]}" for k in fates))


if __name__ == "__main__":
    main()
