#!/usr/bin/env python3
"""foreign_rates.tsv for a protal database: per gene copy, how many reads of a tiled scan of the genomes at hand landed on
it, and how many of them came from another species, and from another genus.

simulate_metagenomes --tiles cuts error-free reads every --stride bases from every genome of the genome table (the
species of --exclude left out: a training database's held-out species, which its models are trained to find missing),
protal aligns them against the database once (single-end, no profile), and each read whose best record (the first
without FLAG 0x904) has MAPQ 4 or more, as the profiler counts reads, counts for that record's copy (RNAME taxid_gene):
from its own species (the read's genome's species, by the GTDB lineage of the genome table, is the copy's species in the
database's internal_taxonomy.dmp), from another species of the genus, or from another genus. A copy that other species'
reads reach (a transferred gene, a contaminating contig of the reference, a gene conserved across its genus) is weaker
evidence of its species: the profiler's "foreign" features (src/SequenceUtils/ForeignRatesTable.h). protal --add_tables
stores the table in the database (docs/claude/2026-10-07-error-read-signatures, section 7 and the follow-up).

    python3 scripts/foreign_rates.py --db DB --genome-table genomes.tsv --exclude heldout_species.txt \\
        --out DB/foreign_rates.tsv -t 32
    protal --add_tables DB/foreign_rates.tsv --db DB
"""
import argparse
import collections
import csv
import glob
import os
import shutil
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compressed  # noqa: E402

MIN_MAPQ = 4  # the profiler's (Profiler::m_min_mapq)
SATURATED = 65535  # the table's counts are 16-bit in a run (ForeignRatesTable.h saturates)


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--db", required=True, help="the database (a folder with database.protal or its separate files)")
    p.add_argument("--genome-table", required=True,
                   help="genomes at hand: accession, GTDB lineage, FASTA path (as simulate_metagenomes takes)")
    p.add_argument("--taxonomy", help="the database's internal_taxonomy.dmp (default: DB/internal_taxonomy.dmp)")
    p.add_argument("--exclude", help="species to leave out of the scan, first column (heldout_species.txt)")
    p.add_argument("--out", required=True, help="the table to write (foreign_rates.tsv)")
    p.add_argument("--length", type=int, default=150, help="read length (default 150)")
    p.add_argument("--stride", type=int, default=500, help="a read every this many bases (default 500)")
    p.add_argument("-t", "--threads", type=int, default=1)
    p.add_argument("--protal", default="protal")
    p.add_argument("--simulate", default="simulate_metagenomes")
    p.add_argument("--genome-store", help="simulate_metagenomes --genome_store (genomes decoded once)")
    p.add_argument("--workdir", help="where the tiles and the SAM go (default: OUT's folder/foreign_rates_scan)")
    p.add_argument("--keep", action="store_true", help="keep the tiles and the SAM")
    return p.parse_args(argv)


def last_rank(lineage, prefix):
    """The lineage's field of a rank ("s__", "g__"), or ""."""
    for field in reversed(lineage.split(";")):
        field = field.strip()
        if field.startswith(prefix):
            return field
    return ""


def read_genomes(path):
    """{accession: (species, genus)} by the GTDB lineage of the genome table (s__ and g__ fields)."""
    out = {}
    with open(path, newline="") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) < 2 or not row[0] or row[0].startswith("#"):
                continue
            out[row[0]] = (last_rank(row[1], "s__"), last_rank(row[1], "g__"))
    return out


def read_taxonomy(path):
    """{taxid: (species name, genus name)} of the database's species (internal_taxonomy.dmp: id, parent, external id,
    name, rank, ...)."""
    parent, name, rank = {}, {}, {}
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if len(row) < 5 or not row[0].isdigit():
                continue
            parent[row[0]], name[row[0]], rank[row[0]] = row[1], row[3], row[4]

    def genus_of(taxid):
        for _ in range(64):
            taxid = parent.get(taxid)
            if taxid is None:
                return ""
            if rank.get(taxid) == "genus":
                return name[taxid]
            if parent.get(taxid) == taxid:
                return ""
        return ""
    return {int(t): (name[t], genus_of(t)) for t in name if rank[t] == "species"}


def run(command, log):
    print("+ " + " ".join(command), flush=True)
    began = time.time()
    with open(log, "w") as fh:
        done = subprocess.run(command, stdout=fh, stderr=subprocess.STDOUT)
    if done.returncode != 0:
        with open(log) as fh:
            sys.stderr.write(fh.read()[-4000:])
        sys.exit(f"{command[0]} failed (exit {done.returncode}); its log: {log}")
    print(f"  {time.time() - began:.0f} s", flush=True)


def count(sam, genomes, taxonomy):
    """{(taxid, gene): [reads, foreign, foreign_genus]} of the counted reads, and the reads seen, counted and of an
    unknown genome."""
    counts = collections.defaultdict(lambda: [0, 0, 0])
    seen, counted, unknown = 0, 0, 0
    last = None
    with compressed.open_read(sam) as fh:
        for raw in fh:
            if raw.startswith(b"@"):
                continue
            f = raw.split(b"\t", 6)
            if len(f) < 6:
                continue
            qname = f[0]
            if qname != last:
                last, done = qname, False
                seen += 1
            flag = int(f[1])
            if done or flag & 0x904:
                continue
            done = True  # the read's best record
            if f[2] == b"*" or int(f[4]) < MIN_MAPQ:
                continue
            rname = f[2].decode()
            taxid, _, gene = rname.partition("_")
            accession = qname.decode().split(":", 1)[0]
            source = genomes.get(accession)
            target = taxonomy.get(int(taxid))
            if source is None or target is None:
                unknown += 1
                continue
            counted += 1
            c = counts[(int(taxid), int(gene))]
            c[0] += 1
            if source[0] != target[0]:
                c[1] += 1
                if source[1] != target[1]:
                    c[2] += 1
    return counts, seen, counted, unknown


def write_table(path, counts):
    """foreign_rates.tsv as ForeignRatesTable.h reads it: per species taxid<TAB>gene:reads:foreign:foreign_genus,...
    (genes ascending); counts saturated at 16 bits."""
    by_taxid = collections.defaultdict(list)
    for (taxid, gene), (reads, foreign, foreign_genus) in counts.items():
        by_taxid[taxid].append((gene, min(reads, SATURATED), min(foreign, SATURATED), min(foreign_genus, SATURATED)))
    partial = path + ".partial"
    with open(partial, "w") as out:
        out.write("# protal foreign rates: reads of a tiled scan per gene copy, of them from other species and genera\n")
        out.write("taxid\tgene:reads:foreign:foreign_genus\n")
        for taxid in sorted(by_taxid):
            entries = sorted(by_taxid[taxid])
            out.write(f"{taxid}\t" + ",".join(f"{g}:{r}:{f}:{fg}" for g, r, f, fg in entries) + "\n")
    os.replace(partial, path)


def main(argv=None):
    opts = parse_args(argv)
    taxonomy_path = opts.taxonomy or os.path.join(opts.db, "internal_taxonomy.dmp")
    if not os.path.isfile(taxonomy_path):
        sys.exit(f"no taxonomy at {taxonomy_path}: give the database's internal_taxonomy.dmp with --taxonomy")
    workdir = opts.workdir or os.path.join(os.path.dirname(os.path.abspath(opts.out)), "foreign_rates_scan")
    os.makedirs(workdir, exist_ok=True)
    tiles = os.path.join(workdir, "tiles.fq.zst")
    command = [opts.simulate, "--tiles", f"{opts.length}:{opts.stride}", "--genome_table", opts.genome_table,
               "--tile_out", tiles, "-t", str(opts.threads)]
    if opts.exclude:
        command += ["--tile_exclude", opts.exclude]
    if opts.genome_store:
        command += ["--genome_store", opts.genome_store]
    run(command, os.path.join(workdir, "tiles.log"))
    for old in glob.glob(os.path.join(workdir, "tiles*.sam*")):
        os.remove(old)
    run([opts.protal, "--db", opts.db, "-1", tiles, "--read_type", "se", "--no_profile", "-o", workdir, "--prefix", "tiles",
         "-t", str(opts.threads), "--force"], os.path.join(workdir, "protal.log"))
    sams = sorted(glob.glob(os.path.join(workdir, "**", "tiles*.sam*"), recursive=True))
    sams = [s for s in sams if not s.endswith(".partial")]
    if len(sams) != 1:
        sys.exit(f"expected one SAM of the tiles in {workdir}, found {sams}")
    genomes = read_genomes(opts.genome_table)
    taxonomy = read_taxonomy(taxonomy_path)
    counts, seen, counted, unknown = count(sams[0], genomes, taxonomy)
    write_table(opts.out, counts)
    copies = len(counts)
    foreign = sum(1 for c in counts.values() if c[1] > 0)
    foreign_genus = sum(1 for c in counts.values() if c[2] > 0)
    print(f"Foreign rates: {seen} reads with a record, {counted} counted (MAPQ >= {MIN_MAPQ}), {unknown} of an unknown "
          f"genome or taxon; {copies} gene copies reached, {foreign} by another species' reads, {foreign_genus} by "
          f"another genus's: {opts.out}")
    if not opts.keep:
        shutil.rmtree(workdir, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
