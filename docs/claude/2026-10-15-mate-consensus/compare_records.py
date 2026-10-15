#!/usr/bin/env python3
"""The primary records of two runs of the same samples side by side (e.g. --mate_guidance guide against consensus):
per fragment each mate's taxon and MAPQ, and, for simulated reads, whether the taxon is the species the fragment was
simulated from (its name: <contig>-<number>; the contig's genome from the build's genome_contigs.tsv.gz, its species
from genomes_simulated.tsv; a record's species from genome2tiid.tsv).

A mate "counts" when its MAPQ is 4 or more (the profiler's cutoff: only those records add to a taxon's fragments).

  per sample (records.tsv): fragments, those with any record in either run, those whose records differ (a mate's
  taxon, or only a MAPQ), those that cross the cutoff; for simulated samples the counted mates on the right and on
  the wrong species in each run, and the fragments whose mates were written on two taxa in each run.

usage: compare_records.py OUT.tsv LABEL_A LABEL_B DATABASE_FOLDER [--work WORK] SAMPLE:SAM_A:SAM_B [...]
  WORK: the build's work folder (genome_contigs.tsv.gz, genomes_simulated.tsv); without it the samples are real.
Standard library only (zstdcat for .zst SAMs).
"""
import argparse
import collections
import gzip
import os
import subprocess
import sys

CUTOFF = 4


def open_sam(path):
    if path.endswith(".zst"):
        proc = subprocess.Popen(["zstdcat", path], stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
        return proc.stdout
    return open(path)


def species_of_taxonomy(taxonomy):
    return taxonomy.rstrip().split(";")[-1]


def load_db_species(folder):
    """tiid -> species name."""
    species = {}
    with open(os.path.join(folder, "genome2tiid.tsv")) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 4:
                species[f[1]] = species_of_taxonomy(f[3])
    return species


def load_truth(work):
    """contig -> species name of the genome it is from."""
    taxonomy = {}
    with open(os.path.join(work, "genomes_simulated.tsv")) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2:
                taxonomy[f[0]] = species_of_taxonomy(f[1])
    contigs = {}
    with gzip.open(os.path.join(work, "genome_contigs.tsv.gz"), "rt") as fh:
        for line in fh:
            path, _, names = line.rstrip("\n").split("\t")
            genome = os.path.basename(path)
            for suffix in (".fna.gz", ".fa.gz", ".fasta.gz", ".fna", ".fa", ".fasta"):
                if genome.endswith(suffix):
                    genome = genome[:-len(suffix)]
                    break
            if genome not in taxonomy:
                continue
            for name in names.split(" "):
                contigs[name] = taxonomy[genome]
    return contigs


def read_records(path):
    """qname -> [taxid1, mapq1, taxid2, mapq2] of the fragment's primary records (None where a mate has none)."""
    fragments = {}
    for line in open_sam(path):
        if line.startswith("@"):
            continue
        f = line.split("\t", 6)
        flag = int(f[1])
        if flag & 0x904:  # unmapped, secondary, supplementary
            continue
        entry = fragments.get(f[0])
        if entry is None:
            entry = fragments[f[0]] = [None, None, None, None]
        taxid = f[2].split("_")[0]
        mate = 2 if flag & 0x80 else 0
        entry[mate] = taxid
        entry[mate + 1] = int(f[4])
    return fragments


def main():
    p = argparse.ArgumentParser()
    p.add_argument("out")
    p.add_argument("label_a")
    p.add_argument("label_b")
    p.add_argument("database")
    p.add_argument("--work", default="")
    p.add_argument("samples", nargs="+")
    args = p.parse_args()
    species = load_db_species(args.database)
    truth = load_truth(args.work) if args.work else None
    a, b = args.label_a, args.label_b
    columns = ["sample", "fragments_with_records", "differ", "taxon_differs", "mapq_differs_only", "cross_cutoff",
               f"split_taxa_{a}", f"split_taxa_{b}"]
    if truth is not None:
        columns += [f"counted_right_{a}", f"counted_wrong_{a}", f"counted_right_{b}", f"counted_wrong_{b}",
                    "moved_to_right", "moved_to_wrong", f"fragments_right_{a}", f"fragments_wrong_{a}",
                    f"fragments_right_{b}", f"fragments_wrong_{b}"]
    rows = []
    for spec in args.samples:
        name, path_a, path_b = spec.split(":")
        ra, rb = read_records(path_a), read_records(path_b)
        c = collections.Counter()
        for qname in set(ra) | set(rb):
            ea = ra.get(qname, [None] * 4)
            eb = rb.get(qname, [None] * 4)
            c["fragments_with_records"] += 1
            if ea != eb:
                c["differ"] += 1
                if (ea[0], ea[2]) != (eb[0], eb[2]):
                    c["taxon_differs"] += 1
                else:
                    c["mapq_differs_only"] += 1
                counted = lambda e, m: e[m] is not None and e[m + 1] >= CUTOFF
                if any(counted(ea, m) != counted(eb, m) for m in (0, 2)):
                    c["cross_cutoff"] += 1
            for label, e in ((a, ea), (b, eb)):
                if e[0] is not None and e[2] is not None and e[0] != e[2]:
                    c[f"split_taxa_{label}"] += 1
            if truth is None:
                continue
            own = truth.get(qname.rsplit("-", 1)[0])
            if own is None:
                continue
            for label, e in ((a, ea), (b, eb)):
                right = wrong = 0
                for m in (0, 2):
                    if e[m] is not None and e[m + 1] >= CUTOFF:
                        if species.get(e[m]) == own:
                            right += 1
                        else:
                            wrong += 1
                c[f"counted_right_{label}"] += right
                c[f"counted_wrong_{label}"] += wrong
                # A fragment counts for the taxa of its counted mates (once per taxon).
                taxa = {species.get(e[m]) for m in (0, 2) if e[m] is not None and e[m + 1] >= CUTOFF}
                c[f"fragments_right_{label}"] += own in taxa
                c[f"fragments_wrong_{label}"] += len(taxa - {own})
            if ea != eb:
                before = {species.get(ea[m]) for m in (0, 2) if ea[m] is not None and ea[m + 1] >= CUTOFF}
                after = {species.get(eb[m]) for m in (0, 2) if eb[m] is not None and eb[m + 1] >= CUTOFF}
                c["moved_to_right"] += (own in after) and (own not in before)
                c["moved_to_wrong"] += len(after - {own} - before)
        row = {"sample": name}
        row.update({k: c.get(k, 0) for k in columns[1:]})
        rows.append(row)
        print(name, " ".join(f"{k}={row[k]}" for k in columns[1:]), file=sys.stderr, flush=True)
    total = {"sample": "all"}
    total.update({k: sum(r[k] for r in rows) for k in columns[1:]})
    rows.append(total)
    with open(args.out, "w") as fh:
        fh.write("\t".join(columns) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[k]) for k in columns) + "\n")


if __name__ == "__main__":
    main()
