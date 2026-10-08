#!/usr/bin/env python3
"""foreign_rates.tsv for a protal database: per gene copy, how many reads of a tiled scan of the full reference's marker
genes landed on it, and how many of them came from another species, and from another genus.

The full reference (full_reference.fna.zst, the converter's: every genome's marker genes under >taxid_geneid of its
species; a training database's lacks the species held out of it) is cut into error-free reads every --stride bases of at
most --per-header copies of each species' gene (simulate_metagenomes --tiles --tile_fasta --tile_per_header), protal
aligns them against the database once (single-end, no profile), and each read whose best record (the first without FLAG
0x904) has MAPQ 4 or more, as the profiler counts reads, counts for that record's copy (RNAME taxid_gene): from its own
species (the read's name starts with the taxid of the genome's species), from another species of the genus (by the
database's internal_taxonomy.dmp), or from another genus. Every copy of the full reference gets a row, reached by a read
or not: a copy that other species' reads reach (a gene conserved across its genus, a strain of a congener that crosses
the species boundary on it) is weaker evidence of its species, the profiler's "foreign" features
(src/SequenceUtils/ForeignRatesTable.h). protal --add_tables stores the table in the database.

The table's reads of a copy are its foreign reads plus one genome's worth of its own species' reads (own reads over the
own records tiled), not every own read: own reads scale with the species' genome count (up to --per-header), and a
share of foreign over all reads would then be ten times higher for a one-genome species than for one with ten, the
cluster-size rule the simulation cannot test (docs/claude/2026-10-03-r226-v9-evaluation, why the priors are opt-in).
The profiler's foreign / (reads + 1) is thus foreign / (foreign + own per genome + 1), whatever the genome count.

Until 2026-10-08 the scan tiled the genomes the training samples are drawn from, and only copies a read reached were
listed: the features then told the models which species the simulation could draw (docs/claude/2026-10-07-congener-gaps,
"The foreign features leak"). Scanning every species' own copies alike, and listing every copy, is what removes that.

    python3 scripts/foreign_rates.py --db DB --full-reference DB/full_reference.fna.zst --taxonomy internal_taxonomy.dmp \\
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
PER_HEADER = 10  # genomes per species and gene tiled (--per-header): E. coli's thousands would swamp the alignment
STRIDE = 250  # a marker gene of 1 kb gives four reads of 150 bases (--stride)


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--db", required=True, help="the database (a folder with database.protal or its separate files)")
    p.add_argument("--full-reference", required=True,
                   help="full_reference.fna[.zst]: every genome's marker genes, >taxid_geneid of its species (a training "
                        "database's lacks its held-out species)")
    p.add_argument("--taxonomy", help="the database's internal_taxonomy.dmp (default: DB/internal_taxonomy.dmp)")
    p.add_argument("--out", required=True, help="the table to write (foreign_rates.tsv)")
    p.add_argument("--length", type=int, default=150, help="read length (default 150)")
    p.add_argument("--stride", type=int, default=STRIDE, help=f"a read every this many bases of a copy (default {STRIDE})")
    p.add_argument("--per-header", type=int, default=PER_HEADER,
                   help=f"at most this many copies of each species' gene tiled, the full reference's first (default "
                        f"{PER_HEADER}; 0: all)")
    p.add_argument("-t", "--threads", type=int, default=1)
    p.add_argument("--protal", default="protal")
    p.add_argument("--simulate", default="simulate_metagenomes")
    p.add_argument("--workdir", help="where the tiles and the SAM go (default: OUT's folder/foreign_rates_scan)")
    p.add_argument("--keep", action="store_true", help="keep the tiles and the SAM")
    return p.parse_args(argv)


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


def read_sources(path):
    """{(taxid, gene): [records seen, records kept, tiles]} of the tiler's <tiles>.sources.tsv (a header's lines summed:
    the tiler lists a header again when its gene's records are not contiguous); headers that are no taxid_gene are
    skipped."""
    sources = collections.defaultdict(lambda: [0, 0, 0])
    skipped = 0
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if len(row) < 4 or row[0] == "header" or row[0].startswith("#"):
                continue
            key = copy_of(row[0])
            if key is None:
                skipped += 1
                continue
            s = sources[key]
            s[0] += int(row[1])
            s[1] += int(row[2])
            s[2] += int(row[3])
    return sources, skipped


def copy_of(header):
    """(taxid, gene) of a taxid_gene header, or None."""
    taxid, sep, gene = header.partition("_")
    if not sep or not taxid.isdigit() or not gene.isdigit():
        return None
    return int(taxid), int(gene)


def run(command, log, summary=None):
    """Runs the command with its output to log; prints its summary line (the last of the log starting with `summary`)
    and how long it took."""
    print("+ " + " ".join(command), flush=True)
    began = time.time()
    with open(log, "w") as fh:
        done = subprocess.run(command, stdout=fh, stderr=subprocess.STDOUT)
    if done.returncode != 0:
        with open(log) as fh:
            sys.stderr.write(fh.read()[-4000:])
        sys.exit(f"{command[0]} failed (exit {done.returncode}); its log: {log}")
    if summary:
        with open(log, errors="replace") as fh:
            lines = [line.rstrip("\n") for line in fh if line.startswith(summary)]
        if lines:
            print("  " + lines[-1], flush=True)
    print(f"  {time.time() - began:.0f} s", flush=True)


def count(sam, taxonomy):
    """{(taxid, gene): [own, foreign, foreign_genus]} of the counted reads (from the copy's own species, from other
    species, of them from other genera), and the reads seen, counted and of an unknown source or target taxon. A read's
    name is <taxid_gene>:<record>:<start> (the tiler's): its source is the species taxid."""
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
            source_key = copy_of(qname.decode().split(":", 1)[0])
            source = taxonomy.get(source_key[0]) if source_key else None
            target = taxonomy.get(int(taxid))
            if source is None or target is None:
                unknown += 1
                continue
            counted += 1
            c = counts[(int(taxid), int(gene))]
            if source_key[0] != int(taxid):
                c[1] += 1
                if source[1] != target[1]:
                    c[2] += 1
            else:
                c[0] += 1
    return counts, seen, counted, unknown


def table_reads(own, foreign, kept):
    """A copy's reads for the table: its foreign reads plus one genome's worth of its own (own reads over the `kept`
    records of its header tiled; 0 without such records), so that the species' genome count drops out of
    foreign / (reads + 1)."""
    return foreign + (round(own / kept) if kept > 0 else 0)


def write_table(path, counts, sources=None):
    """foreign_rates.tsv as ForeignRatesTable.h reads it: per species taxid<TAB>gene:reads:foreign:foreign_genus,...
    (genes ascending), reads as table_reads; counts saturated at 16 bits. `sources` (read_sources: {(taxid, gene):
    [seen, kept, tiles]}) gives each header's records tiled, and every copy of it is listed, counted or not (with no
    reads). Returns the copies listed."""
    sources = sources or {}
    by_taxid = collections.defaultdict(list)
    listed = set()
    for (taxid, gene), (own, foreign, foreign_genus) in counts.items():
        kept = sources.get((taxid, gene), (0, 0, 0))[1]
        by_taxid[taxid].append((gene, min(table_reads(own, foreign, kept), SATURATED), min(foreign, SATURATED),
                                min(foreign_genus, SATURATED)))
        listed.add((taxid, gene))
    for taxid, gene in sources:
        if (taxid, gene) not in listed:
            by_taxid[taxid].append((gene, 0, 0, 0))
            listed.add((taxid, gene))
    partial = path + ".partial"
    with open(partial, "w") as out:
        out.write("# protal foreign rates: reads of a tiled scan of the full reference's marker genes per gene copy (the "
                  "foreign ones plus one genome's worth of the copy's own), of them from other species and genera; every "
                  "copy listed, reached or not\n")
        out.write("taxid\tgene:reads:foreign:foreign_genus\n")
        for taxid in sorted(by_taxid):
            entries = sorted(by_taxid[taxid])
            out.write(f"{taxid}\t" + ",".join(f"{g}:{r}:{f}:{fg}" for g, r, f, fg in entries) + "\n")
    os.replace(partial, path)
    return len(listed)


def main(argv=None):
    opts = parse_args(argv)
    taxonomy_path = opts.taxonomy or os.path.join(opts.db, "internal_taxonomy.dmp")
    if not os.path.isfile(taxonomy_path):
        sys.exit(f"no taxonomy at {taxonomy_path}: give the database's internal_taxonomy.dmp with --taxonomy")
    if not os.path.isfile(opts.full_reference):
        sys.exit(f"no full reference at {opts.full_reference}")
    workdir = opts.workdir or os.path.join(os.path.dirname(os.path.abspath(opts.out)), "foreign_rates_scan")
    os.makedirs(workdir, exist_ok=True)
    tiles = os.path.join(workdir, "tiles.fq.zst")
    run([opts.simulate, "--tiles", f"{opts.length}:{opts.stride}", "--tile_fasta", opts.full_reference, "--tile_per_header",
         str(max(0, opts.per_header)), "--tile_out", tiles, "-t", str(opts.threads)], os.path.join(workdir, "tiles.log"),
        summary="Tiles:")
    for old in glob.glob(os.path.join(workdir, "tiles*.sam*")):
        os.remove(old)
    run([opts.protal, "--db", opts.db, "-1", tiles, "--read_type", "se", "--no_profile", "-o", workdir, "--prefix", "tiles",
         "-t", str(opts.threads), "--force"], os.path.join(workdir, "protal.log"))
    sams = sorted(glob.glob(os.path.join(workdir, "**", "tiles*.sam*"), recursive=True))
    sams = [s for s in sams if not s.endswith(".partial")]
    if len(sams) != 1:
        sys.exit(f"expected one SAM of the tiles in {workdir}, found {sams}")
    taxonomy = read_taxonomy(taxonomy_path)
    sources, skipped = read_sources(tiles + ".sources.tsv")
    counts, seen, counted, unknown = count(sams[0], taxonomy)
    listed = write_table(opts.out, counts, sources)
    reached = sum(1 for c in counts.values() if c[0] > 0 or c[1] > 0)
    foreign = sum(1 for c in counts.values() if c[1] > 0)
    foreign_genus = sum(1 for c in counts.values() if c[2] > 0)
    genomes = sum(s[1] for s in sources.values())
    print(f"Foreign rates: {seen} reads with a record, {counted} counted (MAPQ >= {MIN_MAPQ}), {unknown} of an unknown "
          f"taxon, from {genomes} gene copies of {len(sources)} species' genes (at most {opts.per_header or 'all'} each"
          f"{f', {skipped} headers that are no taxid_gene skipped' if skipped else ''}); {listed} gene copies listed, "
          f"{reached} reached by a counted read, {foreign} by another species' reads, {foreign_genus} by another genus's: "
          f"{opts.out}")
    if not opts.keep:
        shutil.rmtree(workdir, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
