#!/usr/bin/env python3
"""The alignments of the reads behind a presence model's false positives and false negatives, with their features,
from the SAMs a database build left on its scratch disk: what tells the reads of an FP taxon from those of a present
one, and where an FN taxon's reads went.

Written for the r226 v15 build, whose paired-end error reads were lost (error_reads.py failed for pe), to run on the
cluster node that holds the samples. It reads only the build's SAMs, the trainer's calls and the training database's
species, and takes no reads of "unseen" species (present species without a row). One pass over each sample's SAM
(protal writes a read's records one after the other) gives, per sample:

  OUT/<set>/<point>/<sample>.sam.zst         the records of the fragments behind the sample's FP and FN taxa, QUAL *:
                                             FP:<taxid>      aligned to an FP taxon
                                             FN:<taxid>      aligned to an FN taxon
                                             seeded:<taxid>  seeded on an FN taxon but did not align to it (ZF)
                                             source:<taxid>  simulated from a genome of an FN species, wherever it went
                                             at most --max-fragments fragments per taxon and reason (the lowest CRC-32
                                             of the read name, as error_reads.py takes them), each record tagged
                                             xg:Z:<source genome> xs:Z:<source species> xe:Z:<its reasons>
  OUT/<set>/<point>/<sample>.records.tsv.gz  one line per record of those fragments: its alignment (CIGAR counts,
                                             identity, clips, gene ends, mate), its tags (ZU, ZT, ZA, ZF, ZR), its
                                             taxon's call, and the read's source and how it relates to the taxon (own
                                             species, same genus, same family, other)
  OUT/<set>/<point>/<sample>.taxa.tsv.gz     one line per row of the sample (TP, FP, FN and TN alike, from all its
                                             records): the fragments the profiler counts for the taxon (best record on
                                             it at MAPQ 4 or more) and where they came from, means and shares of their
                                             alignment features, and where the reads of the taxon's own species went
and OUT/summary.tsv (one line per sample), OUT/failed.txt (the traceback of each sample that failed; the others go on).

    python3 docs/claude/2026-10-07-error-read-signatures/error_read_features.py --build OUTDIR --scratch SCRATCH \\
        --read-type pe --out OUTDIR/error_read_features/pe --threads 16

--build gives the trainer's calls (OUTDIR/trained_model[_<type>].calls.tsv.gz, or under model_logs/), the taxonomy
(internal_taxonomy.dmp) and heldout_species.txt; --scratch the collections (SCRATCH/training, SCRATCH/test), the
training database (SCRATCH/training_db) and the genomes' contig names read before (SCRATCH/genome_contigs.tsv.gz).
Each can be given on its own instead. Standard library only, and the scripts of this checkout (error_reads.py,
trace_relatives.py, compressed.py: the zstd command or Python 3.14).
"""
import argparse
import collections
import concurrent.futures
import csv
import gzip
import heapq
import multiprocessing
import os
import re
import sys
import time
import traceback
import zlib

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.normpath(os.path.join(HERE, "..", "..", "..", "scripts"))
sys.path.insert(0, SCRIPTS)
import compressed  # noqa: E402
import error_reads  # noqa: E402
import trace_relatives  # noqa: E402

MIN_MAPQ = 4  # Profiler::m_min_mapq: a fragment counts for the taxon of its best record at this MAPQ or more
CLIP = 10  # a record clipped by at least this many bases at one end is "clipped"
MAX_FRAGMENTS = 20  # as error_reads.py
CIGAR = re.compile(rb"(\d+)([MIDNSHP=X])")
RECENT = 4096  # read names remembered to tell a read whose records are not one after the other

RECORD_COLUMNS = [
    "set", "scenario", "point", "sample", "qname", "reasons",
    "source_genome", "source_species", "source_taxid", "source_genome_kind", "source_status", "source_call",
    "frag_records", "frag_aligned_taxa", "frag_best_taxon", "frag_best_mapq", "frag_zf_n",
    "flag", "mate", "primary", "proper", "mate_unmapped", "unmapped",
    "taxon", "gene", "pos", "gene_len", "mapq", "read_len", "aln_query", "aln_ref", "matches", "mismatches", "ins",
    "dels", "gap_opens", "clip_left", "clip_right", "identity", "at_gene_start", "at_gene_end", "rnext", "pnext",
    "tlen", "zu", "zt", "za", "za_n", "za_min", "za_ties", "za_source", "zf", "zf_n", "zf_source", "zr",
    "taxon_call", "taxon_p", "knob", "relation"]
# The taxa table: the row's call, then what its fragments are (counted: best record on the taxon at MAPQ >= 4).
TAXA_COLUMNS = [
    "set", "scenario", "point", "sample", "taxon", "taxon_name", "truth", "call", "class", "p", "knob",
    "frags_any", "frags_best", "frags_counted", "frags_lowmapq",
    "src_own", "src_genus", "src_family", "src_other", "src_host", "src_unknown", "src_held_out", "src_strain",
    "src_insilico", "src_species", "src_top_share",
    "identity_mean", "identity_min", "identity_sd", "mapq_mean", "mapq_max_share", "zu_mean", "zt_mean", "zu0_share",
    "zt0_share", "za_share", "za_tie_share", "za_close_share", "za_called_share", "zf_share", "zf_mean",
    "zf_called_share", "clip_share", "edge_share",
    "both_mates_share", "split_genes_share", "other_taxon_share", "proper_share", "mate_unmapped_share", "zr_share",
    "secondary_share", "aln_mean", "genes", "fragments_per_gene_max",
    "own_frags", "own_counted_here", "own_lowmapq_here", "own_counted_elsewhere", "own_lowmapq_elsewhere",
    "own_unaligned", "seeded_failed"]
SUMMARY_COLUMNS = ["set", "point", "scenario", "sample", "rows", "FP", "FN", "fragments", "records",
                   "aligned_records", "unknown_source", "not_consecutive", "kept_fragments", "kept_records",
                   "sam_bytes", "out_bytes", "seconds"]

# Set before the worker processes fork (main): the taxonomy {id: (parent, name, rank, representative genome)},
# {species name: id}, the held-out species and the species with genes in the training database (names).
TAXONOMY, NAME_ID, HELD_OUT, DB_SPECIES = {}, {}, set(), set()


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--build", help="the build's OUTDIR (calls, internal_taxonomy.dmp, heldout_species.txt)")
    p.add_argument("--scratch", help="the build's --scratch (training, test, training_db, genome_contigs.tsv.gz)")
    p.add_argument("--calls", help="the trainer's PREFIX.calls.tsv.gz (default: from --build and --read-type)")
    p.add_argument("--training", help="the training collection (default SCRATCH/training)")
    p.add_argument("--test", help="the test collection (default SCRATCH/test)")
    p.add_argument("--db", help="the training database's folder (default SCRATCH/training_db)")
    p.add_argument("--taxonomy", help="internal_taxonomy.dmp (default OUTDIR's, else the database's)")
    p.add_argument("--heldout", help="heldout_species.txt (default OUTDIR's)")
    p.add_argument("--contig-cache",
                   help="the genomes' contig names read before (default SCRATCH/genome_contigs.tsv.gz)")
    p.add_argument("--read-type", default="pe", choices=("pe", "se", "pb", "ont"), help="default pe")
    p.add_argument("--samples", default="all", help="all (default), or design and scenario names, comma-separated")
    p.add_argument("--limit", type=int, default=0, help="only this many samples, the smallest SAMs first (a check)")
    p.add_argument("--max-fragments", type=int, default=MAX_FRAGMENTS,
                   help=f"fragments per taxon and reason in the SAMs and records tables (default {MAX_FRAGMENTS}; 0: "
                        "all); the taxa tables count all")
    p.add_argument("--qualities", action="store_true", help="keep the records' QUAL (default *)")
    p.add_argument("--threads", type=int, default=4, help="samples at once (default 4)")
    p.add_argument("--out", required=True, help="output folder")
    opts = p.parse_args(argv)

    def pick(given, *candidates):
        if given:
            return given
        return next((c for c in candidates if c and os.path.exists(c)), None)

    b, s = opts.build, opts.scratch
    prefix = "trained_model" + ("" if opts.read_type == "pe" else "_" + opts.read_type)
    opts.calls = pick(opts.calls, b and os.path.join(b, prefix + ".calls.tsv.gz"),
                      b and os.path.join(b, "model_logs", prefix + ".calls.tsv.gz"))
    opts.training = pick(opts.training, s and os.path.join(s, "training"))
    opts.test = pick(opts.test, s and os.path.join(s, "test"))
    opts.db = pick(opts.db, s and os.path.join(s, "training_db"), b and os.path.join(b, "training_db"))
    opts.taxonomy = pick(opts.taxonomy, b and os.path.join(b, "internal_taxonomy.dmp"),
                         opts.db and os.path.join(opts.db, "internal_taxonomy.dmp"))
    opts.heldout = pick(opts.heldout, b and os.path.join(b, "heldout_species.txt"),
                        b and os.path.join(b, "model_logs", "heldout_species.txt"))
    opts.contig_cache = pick(opts.contig_cache, s and os.path.join(s, "genome_contigs.tsv.gz"))
    missing = [name for name, value in (("--calls", opts.calls), ("--db", opts.db), ("--taxonomy", opts.taxonomy))
               if not value]
    if missing:
        p.error(f"not found: {', '.join(missing)} (give --build and --scratch, or each)")
    if not opts.training and not opts.test:
        p.error("no collection: give --scratch, or --training/--test")
    opts.samples = {x.strip() for x in opts.samples.split(",") if x.strip()}
    return opts


def read_taxonomy(path):
    """{id: (parent, name, rank, representative genome)} of an internal_taxonomy.dmp, and {species name: id}."""
    nodes, names = {}, {}
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader, None)
        for row in reader:
            if len(row) < 5:
                continue
            nodes[row[0]] = (row[1], row[3], row[4], row[6] if len(row) > 6 else "")
            if row[4] == "species":
                names.setdefault(row[3], row[0])
    return nodes, names


def ancestor(taxid, rank):
    """The id of a taxon's ancestor of a rank (the taxon itself if it has that rank), or None."""
    seen = 0
    while taxid in TAXONOMY and seen < 64:
        parent, _, r, _ = TAXONOMY[taxid]
        if r == rank:
            return taxid
        if parent == taxid:
            return None
        taxid, seen = parent, seen + 1
    return None


def relation(source_taxid, taxid, source_species):
    """How a read's source relates to the taxon it aligned to: own, genus, family, other, host or unknown."""
    if source_species == "host":
        return "host"
    if not source_taxid or taxid not in TAXONOMY:
        return "unknown"
    if source_taxid == taxid:
        return "own"
    for rank in ("genus", "family"):
        a = ancestor(source_taxid, rank)
        if a is not None and a == ancestor(taxid, rank):
            return rank
    return "other"


def cigar_counts(cigar):
    """(matches, mismatches, insertions, deletions, gap opens, clip left, clip right, query bases incl. clips) of a
    CIGAR whose M are matches (protal writes X for mismatches; = counts as M)."""
    m = x = i = d = opens = 0
    ops = CIGAR.findall(cigar)
    left = right = 0
    k = 0
    while k < len(ops) and ops[k][1] in b"SH":
        left += int(ops[k][0])
        k += 1
    j = len(ops)
    while j > k and ops[j - 1][1] in b"SH":
        right += int(ops[j - 1][0])
        j -= 1
    for n, op in ops[k:j]:
        n = int(n)
        if op == b"M" or op == b"=":
            m += n
        elif op == b"X":
            x += n
        elif op == b"I":
            i += n
            opens += 1
        elif op == b"D":
            d += n
            opens += 1
    return m, x, i, d, opens, left, right, left + right + m + x + i


class Record:
    """A SAM record's fields that the features need."""
    __slots__ = ("line", "flag", "unmapped", "taxon", "gene", "pos", "mapq", "m", "x", "i", "d", "opens", "left",
                 "right", "read_len", "rnext", "pnext", "tlen", "zu", "zt", "za", "zf", "zr", "gene_len")

    def __init__(self, line, gene_len):
        f = line.rstrip(b"\r\n").split(b"\t")
        self.line = line
        self.flag = int(f[1])
        self.unmapped = bool(self.flag & 4) or f[2] == b"*"
        self.zu = self.zt = self.zr = 0
        self.za, self.zf = None, ()
        for tag in f[11:]:
            key = tag[:2]
            if key == b"ZU":
                self.zu = int(tag[5:])
            elif key == b"ZT":
                self.zt = int(tag[5:])
            elif key == b"ZF":  # taxid or taxid:gene, comma-separated
                self.zf = tuple(t.split(b":", 1)[0].decode() for t in tag[5:].split(b",") if t)
            elif key == b"ZA":
                self.za = tag[5:].decode()
            elif key == b"ZR":
                self.zr = int(tag[5:])
        if self.unmapped:
            self.taxon = self.gene = None
            return
        taxon, _, gene = f[2].partition(b"_")
        self.taxon, self.gene = taxon.decode(), gene.decode()
        self.pos, self.mapq = int(f[3]), int(f[4])
        (self.m, self.x, self.i, self.d, self.opens, self.left, self.right,
         self.read_len) = cigar_counts(f[5])
        self.rnext, self.pnext, self.tlen = f[6].decode(), int(f[7]), int(f[8])
        self.gene_len = gene_len.get(f[2], 0)

    @property
    def aln_query(self):
        return self.m + self.x + self.i

    @property
    def aln_ref(self):
        return self.m + self.x + self.d

    @property
    def identity(self):
        n = self.m + self.x + self.i + self.d
        return self.m / n if n else 0.0

    @property
    def at_start(self):
        return self.pos == 1

    @property
    def at_end(self):
        return bool(self.gene_len) and self.pos + self.aln_ref - 1 >= self.gene_len

    def alternatives(self):
        """[(taxid, edits more)] of the ZA tag."""
        if not self.za or self.za == "*":
            return []
        out = []
        for alt in self.za.split(","):
            taxid, _, more = alt.partition(":")
            try:
                out.append((taxid, int(more)))
            except ValueError:
                out.append((taxid, 99))
        return out


class Taxon:
    """The counts of a taxon's fragments, and sums of the features of those counted for it (c["identity"] and
    c["identity2"]: the sums of their identities and squared identities)."""
    __slots__ = ("c", "identity_min", "sources", "genes")

    def __init__(self):
        self.c = collections.Counter()
        self.identity_min = None
        self.sources = collections.Counter()
        self.genes = collections.Counter()


def mean(total, n):
    return round(total / n, 6) if n else ""


def taxa_rows(meta, calls, names, taxa, own, failed):
    """The taxa table's lines: one per row of the sample."""
    for taxid, (cls, row) in calls.items():
        t = taxa.get(taxid) or Taxon()
        c, n = t.c, t.c["counted"]
        sd = ""
        if n > 1:
            sd = round(max(0.0, (c["identity2"] - c["identity"] ** 2 / n) / (n - 1)) ** 0.5, 6)
        o = own.get(taxid, collections.Counter())
        top = max(t.sources.values()) if t.sources else 0
        yield dict(meta, taxon=taxid, taxon_name=names.get(taxid, row.get("taxon_name", "")), truth=row["truth"],
                   call=row["call"], **{"class": cls}, p=row.get("p", ""), knob=row.get("knob", ""),
                   frags_any=c["any"], frags_best=c["best"], frags_counted=n, frags_lowmapq=c["lowmapq"],
                   src_own=c["src_own"], src_genus=c["src_genus"], src_family=c["src_family"],
                   src_other=c["src_other"], src_host=c["src_host"], src_unknown=c["src_unknown"],
                   src_held_out=c["src_held_out"], src_strain=c["src_strain"], src_insilico=c["src_insilico"],
                   src_species=len(t.sources), src_top_share=mean(top, n),
                   identity_mean=mean(c["identity"], n),
                   identity_min=round(t.identity_min, 6) if t.identity_min is not None else "",
                   identity_sd=sd, mapq_mean=mean(c["mapq"], n), mapq_max_share=mean(c["mapq_max"], n),
                   zu_mean=mean(c["zu"], n), zt_mean=mean(c["zt"], n), zu0_share=mean(c["zu0"], n),
                   zt0_share=mean(c["zt0"], n), za_share=mean(c["za"], n), za_tie_share=mean(c["za_tie"], n),
                   za_close_share=mean(c["za_close"], n), za_called_share=mean(c["za_called"], n),
                   zf_share=mean(c["zf"], n), zf_mean=mean(c["zf_n"], n), zf_called_share=mean(c["zf_called"], n),
                   clip_share=mean(c["clip"], n), edge_share=mean(c["edge"], n),
                   both_mates_share=mean(c["both_mates"], n), split_genes_share=mean(c["split_genes"], n),
                   other_taxon_share=mean(c["other_taxon"], n), proper_share=mean(c["proper"], n),
                   mate_unmapped_share=mean(c["mate_unmapped"], n), zr_share=mean(c["zr"], n),
                   secondary_share=mean(c["secondary"], n), aln_mean=mean(c["aln"], n), genes=len(t.genes),
                   fragments_per_gene_max=max(t.genes.values()) if t.genes else 0,
                   own_frags=o["frags"], own_counted_here=o["counted_here"], own_lowmapq_here=o["lowmapq_here"],
                   own_counted_elsewhere=o["counted_elsewhere"], own_lowmapq_elsewhere=o["lowmapq_elsewhere"],
                   own_unaligned=o["unaligned"], seeded_failed=failed.get(taxid, 0))


def features(job):
    """One sample: its taxa table, and the records of its FP and FN fragments (SAM and table). -> its summary row."""
    began = time.time()
    which, point, scenario, sample, rows, sam, manifest, drawn, contigs, out_dir, cap, qualities = job
    error_reads.CONTIGS.clear()
    error_reads.CONTIGS.update(contigs)
    source = error_reads.source_finder(manifest, drawn)
    calls = {}
    for row in rows:
        truth, call = row["truth"].lower() in ("1", "true"), row["call"] == "1"
        calls[row["taxon"]] = ("TP" if truth and call else "FN" if truth else "FP" if call else "TN", row)
    fp = {t for t, (cls, _) in calls.items() if cls == "FP"}
    fn = {t for t, (cls, _) in calls.items() if cls == "FN"}
    # The taxa the sample calls: a read whose other candidates (ZA) or failed seeds (ZF) are on them may be theirs.
    called = {t for t, (_, row) in calls.items() if row["call"] == "1"}
    names = {t: TAXONOMY[t][1] for t in calls if t in TAXONOMY}

    genome_kind = {}

    def kind_of_genome(genome, sid, species):
        if genome is None:
            return "unknown"
        if species == "host" or genome == "host":
            return "host"
        k = genome_kind.get(genome)
        if k is None:
            rep = TAXONOMY.get(sid, ("", "", "", ""))[3] if sid else ""
            k = genome_kind[genome] = ("insilico" if genome.startswith("insilico_") else
                                       "rep" if rep and genome == rep else "strain")
        return k

    def status_of(species):
        if not species:
            return "unknown"
        if species == "host":
            return "host"
        if species in HELD_OUT:
            return "held_out"
        return "db" if species in DB_SPECIES else "absent"

    relations = {}

    def related(sid, taxid, species):
        key = (sid, taxid, species == "host")
        r = relations.get(key)
        if r is None:
            r = relations[key] = relation(sid, taxid, species)
        return r

    header, gene_len = [], {}
    taxa = collections.defaultdict(Taxon)
    own = collections.defaultdict(collections.Counter)
    failed = collections.Counter()
    heaps = collections.defaultdict(list)  # reason -> max-heap of (-CRC, qname): the fragments it keeps
    held = {}  # qname -> [fragment, heaps holding it]
    stats = collections.Counter()
    recent = collections.OrderedDict()

    def fragment(qname, lines):
        recs = [Record(line, gene_len) for line in lines]
        stats["fragments"] += 1
        stats["records"] += len(recs)
        aligned = [r for r in recs if not r.unmapped]
        stats["aligned_records"] += len(aligned)
        genome, species = source(qname)
        stats["unknown_source"] += genome is None
        sid = NAME_ID.get(species) if species and species != "host" else None
        best = next((r for r in aligned if not r.flag & 0x904), None) or (aligned[0] if aligned else None)
        zf = set()
        for r in recs:
            zf.update(r.zf)
        on = collections.defaultdict(list)
        for r in aligned:
            on[r.taxon].append(r)
        for taxid in on:
            taxa[taxid].c["any"] += 1
        for taxid in zf:
            failed[taxid] += 1
        if best is not None:
            t = taxa[best.taxon]
            t.c["best"] += 1
            if best.mapq >= MIN_MAPQ:
                mine = on[best.taxon]
                c = t.c
                c["counted"] += 1
                n = sum(r.m + r.x + r.i + r.d for r in mine)
                identity = sum(r.m for r in mine) / n if n else 0.0
                c["identity"] += identity
                c["identity2"] += identity * identity
                if t.identity_min is None or identity < t.identity_min:
                    t.identity_min = identity
                c["mapq"] += best.mapq
                c["mapq_max"] += best.mapq >= 60
                zu, zt = sum(r.zu for r in mine), sum(r.zt for r in mine)
                c["zu"] += zu
                c["zt"] += zt
                c["zu0"] += zu == 0
                c["zt0"] += zt == 0
                alts = best.alternatives()
                c["za"] += bool(alts)
                c["za_tie"] += any(more == 0 for _, more in alts)
                c["za_close"] += any(more <= 1 for _, more in alts)
                c["za_called"] += any(t in called and t != best.taxon for t, _ in alts)
                c["zf"] += bool(zf)
                c["zf_n"] += len(zf)
                c["zf_called"] += any(x in called and x != best.taxon for x in zf)
                c["clip"] += any(r.left >= CLIP or r.right >= CLIP for r in mine)
                c["edge"] += any(r.at_start or r.at_end for r in mine)
                c["both_mates"] += (any(r.flag & 0x40 for r in mine) and any(r.flag & 0x80 for r in mine))
                c["split_genes"] += len({r.gene for r in mine}) > 1
                c["other_taxon"] += len(on) > 1
                c["proper"] += bool(best.flag & 2)
                c["mate_unmapped"] += bool(best.flag & 1 and best.flag & 8)
                c["zr"] += any(r.zr for r in mine)
                c["secondary"] += any(r.flag & 0x900 for r in mine)
                c["aln"] += sum(r.aln_query for r in mine)
                t.genes[best.gene] += 1
                rel = related(sid, best.taxon, species)
                c["src_" + rel] += 1
                status = status_of(species)
                c["src_held_out"] += status == "held_out"
                kind = kind_of_genome(genome, sid, species)
                c["src_strain"] += kind == "strain"
                c["src_insilico"] += kind == "insilico"
                t.sources[species or "?"] += 1
            else:
                t.c["lowmapq"] += 1
        if sid:
            o = own[sid]
            o["frags"] += 1
            if best is None:
                o["unaligned"] += 1
            elif best.taxon == sid:
                o["counted_here" if best.mapq >= MIN_MAPQ else "lowmapq_here"] += 1
            else:
                o["counted_elsewhere" if best.mapq >= MIN_MAPQ else "lowmapq_elsewhere"] += 1

        reasons = set()
        for taxid in on:
            if taxid in fp:
                reasons.add("FP:" + taxid)
            elif taxid in fn:
                reasons.add("FN:" + taxid)
        reasons.update("seeded:" + t for t in zf & fn)
        if sid in fn:
            reasons.add("source:" + sid)
        if not reasons:
            return
        item = (-zlib.crc32(qname), qname)
        entry = [(recs, reasons, genome, species, sid, best, len(zf)), 0]
        for reason in reasons:
            h = heaps[reason]
            if cap and len(h) >= cap:
                if item <= h[0]:
                    continue
                _, dropped = heapq.heapreplace(h, item)
                held[dropped][1] -= 1
                if held[dropped][1] == 0:
                    del held[dropped]
            else:
                heapq.heappush(h, item)
            if qname not in held:
                held[qname] = entry
            held[qname][1] += 1

    with compressed.open_read(sam) as fh:
        current, lines = None, []
        for line in fh:
            if line.startswith(b"@"):
                if current is None:
                    header.append(line)
                    if line.startswith(b"@SQ\t"):
                        f = dict(x.split(b":", 1) for x in line.rstrip(b"\r\n").split(b"\t")[1:] if b":" in x)
                        if b"SN" in f and b"LN" in f:
                            gene_len[f[b"SN"]] = int(f[b"LN"])
                continue
            if not line.strip():
                continue
            qname = line.split(b"\t", 1)[0]
            if qname != current:
                if lines:
                    fragment(current, lines)
                    recent[current] = None
                    if len(recent) > RECENT:
                        recent.popitem(last=False)
                if qname in recent:
                    stats["not_consecutive"] += 1
                current, lines = qname, []
            lines.append(line)
        if lines:
            fragment(current, lines)

    meta = {"set": which, "scenario": scenario or "design", "point": point, "sample": sample}
    os.makedirs(out_dir, exist_ok=True)
    taxa_path = os.path.join(out_dir, sample + ".taxa.tsv.gz")
    with gzip.open(taxa_path + ".partial", "wt", newline="", compresslevel=6) as fh:
        writer = csv.DictWriter(fh, TAXA_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(taxa_rows(meta, calls, names, taxa, own, failed))
    os.replace(taxa_path + ".partial", taxa_path)

    kept = sorted(held.items(), key=lambda kv: kv[0])
    genes, n_records = set(), 0
    sam_path = os.path.join(out_dir, sample + ".sam.zst")
    records_path = os.path.join(out_dir, sample + ".records.tsv.gz")
    with compressed.open_write(sam_path + ".partial") as out, \
            gzip.open(records_path + ".partial", "wt", newline="", compresslevel=6) as table:
        for qname, ((recs, _, _, _, _, _, _), _) in kept:
            for r in recs:
                if not r.unmapped:
                    genes.add(f"{r.taxon}_{r.gene}".encode())
                    if r.rnext not in ("=", "*"):
                        genes.add(r.rnext.encode())
        lines = [h for h in header if not h.startswith(b"@SQ\t") or
                 h.split(b"\t", 2)[1][3:].rstrip(b"\r\n") in genes]
        lines = [h for h in lines if not h.startswith(b"@CO\tprotal failed candidates")]
        counts = collections.Counter(cls for cls, _ in calls.values())
        lines.append((f"@CO\terror_read_features.py: the records of the fragments behind the {counts['FP']} FP and "
                      f"{counts['FN']} FN taxa of {sample} ({which} set, {scenario or 'design'}), at most "
                      f"{cap or 'all'} per taxon and reason; xg and xs: the read's source genome and species, xe: why it was "
                      f"taken\n").encode())
        out.write(b"".join(lines))
        writer = csv.DictWriter(table, RECORD_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for qname, ((recs, reasons, genome, species, sid, best, zf_n), _) in kept:
            why = ",".join(sorted(reasons))
            status, kind = status_of(species), kind_of_genome(genome, sid, species)
            source_call = (calls[sid][0] if sid in calls else "held_out" if status == "held_out" else
                           "unseen" if status == "db" else status)
            frag = {"qname": qname.decode(), "reasons": why, "source_genome": genome or "",
                    "source_species": species or "", "source_taxid": sid or "", "source_genome_kind": kind,
                    "source_status": status, "source_call": source_call, "frag_records": len(recs),
                    "frag_aligned_taxa": len({r.taxon for r in recs if not r.unmapped}),
                    "frag_best_taxon": best.taxon if best else "", "frag_best_mapq": best.mapq if best else "",
                    "frag_zf_n": zf_n}
            tags = f"\txg:Z:{genome or '?'}\txs:Z:{species or '?'}\txe:Z:{why}".encode()
            for r in recs:
                f = r.line.rstrip(b"\r\n").split(b"\t")
                if not qualities and len(f) > 10:
                    f[10] = b"*"
                out.write(b"\t".join(f) + tags + b"\n")
                n_records += 1
                row = dict(meta, **frag, flag=r.flag, mate=1 if r.flag & 0x40 else 2 if r.flag & 0x80 else 0,
                           primary=int(not r.flag & 0x900), proper=int(bool(r.flag & 2)),
                           mate_unmapped=int(bool(r.flag & 1 and r.flag & 8)), unmapped=int(r.unmapped),
                           zu=r.zu, zt=r.zt, zf=",".join(r.zf), zf_n=len(r.zf),
                           zf_source=int(bool(sid) and sid in r.zf), zr=r.zr)
                if not r.unmapped:
                    alts = r.alternatives()
                    cls, call_row = calls.get(r.taxon, ("", {}))
                    row.update(taxon=r.taxon, gene=r.gene, pos=r.pos, gene_len=r.gene_len, mapq=r.mapq,
                               read_len=r.read_len, aln_query=r.aln_query, aln_ref=r.aln_ref, matches=r.m,
                               mismatches=r.x, ins=r.i, dels=r.d, gap_opens=r.opens, clip_left=r.left,
                               clip_right=r.right, identity=round(r.identity, 6), at_gene_start=int(r.at_start),
                               at_gene_end=int(r.at_end), rnext=r.rnext, pnext=r.pnext, tlen=r.tlen,
                               za=r.za or "", za_n=len(alts), za_min=min((m for _, m in alts), default=""),
                               za_ties=sum(m == 0 for _, m in alts), za_source=int(any(t == sid for t, _ in alts)),
                               taxon_call=cls, taxon_p=call_row.get("p", ""), knob=call_row.get("knob", ""),
                               relation=related(sid, r.taxon, species))
                writer.writerow(row)
    os.replace(sam_path + ".partial", sam_path)
    os.replace(records_path + ".partial", records_path)
    counts = collections.Counter(cls for cls, _ in calls.values())
    out_bytes = sum(os.path.getsize(p) for p in (taxa_path, sam_path, records_path))
    return dict(meta, rows=len(calls), FP=counts["FP"], FN=counts["FN"], fragments=stats["fragments"],
                records=stats["records"], aligned_records=stats["aligned_records"],
                unknown_source=stats["unknown_source"], not_consecutive=stats["not_consecutive"],
                kept_fragments=len(kept), kept_records=n_records, sam_bytes=os.path.getsize(sam),
                out_bytes=out_bytes, seconds=round(time.time() - began, 1))


def safe(job):
    """features(job), or (sample, traceback) if it failed."""
    try:
        return features(job)
    except Exception:  # noqa: BLE001 - one sample's failure is reported, the others go on
        return (job[3], traceback.format_exc())


def main(argv=None):
    opts = parse_args(argv)
    began = time.time()
    nodes, names = read_taxonomy(opts.taxonomy)
    TAXONOMY.update(nodes)
    NAME_ID.update(names)
    by_name, _ = error_reads.database_species(opts.db)
    if opts.heldout:
        with open(opts.heldout) as fh:
            HELD_OUT.update(line.split("\t")[0].strip() for line in fh if line.strip())
    DB_SPECIES.update(n for n in by_name if n not in HELD_OUT)
    print(f"taxonomy {opts.taxonomy}: {len(NAME_ID)} species; training database {opts.db}: {len(DB_SPECIES)} species "
          f"with genes, {len(HELD_OUT)} held out ({opts.heldout or 'no list'})", flush=True)
    calls = error_reads.read_calls(opts.calls, opts.samples, opts.read_type)
    print(f"calls {opts.calls}: {sum(len(r) for r in calls.values())} rows of {len(calls)} {opts.read_type} samples",
          flush=True)
    collections_ = {"training": opts.training, "test": opts.test}
    jobs, missing = [], []
    for (which, point, sample), rows in sorted(calls.items()):
        collection = collections_.get(which)
        files = error_reads.sample_files(collection, point, sample, opts.read_type) if collection else None
        if files is None:
            missing.append(f"{sample} ({which})")
            continue
        sam, manifest, drawn = files
        jobs.append([which, point, rows[0].get("meta_scenario") or "", sample, rows, sam, manifest, drawn, None,
                     os.path.join(opts.out, which, point), opts.max_fragments, opts.qualities])
    if missing:
        print(f"no SAM for {len(missing)} samples (not here, or removed): {', '.join(missing[:6])}"
              f"{' ...' if len(missing) > 6 else ''}", flush=True)
    if not jobs:
        print(f"no {opts.read_type} samples with their SAMs: nothing to do", flush=True)
        return 1
    if opts.limit:
        jobs = sorted(jobs, key=lambda job: os.path.getsize(job[5]))[:opts.limit]
    fastas = sorted({r.get("fasta_path") for job in jobs if not job[7] for r in job[6] if r.get("fasta_path")})
    contigs = trace_relatives.genome_contigs(fastas, opts.threads, opts.contig_cache) if fastas else {}
    for job in jobs:  # each worker gets the contig names of its sample's genomes only
        job[8] = {} if job[7] else {r.get("fasta_path"): contigs.get(r.get("fasta_path"), ()) for r in job[6]}
    del contigs
    print(f"{len(jobs)} samples, the contigs of {len(fastas)} genomes in {time.time() - began:.0f} s", flush=True)
    jobs.sort(key=lambda job: -os.path.getsize(job[5]))  # the largest SAMs first
    done, failures = [], []
    if opts.threads > 1 and len(jobs) > 1 and "fork" in multiprocessing.get_all_start_methods():
        with concurrent.futures.ProcessPoolExecutor(min(opts.threads, len(jobs)),
                                                    mp_context=multiprocessing.get_context("fork")) as pool:
            for result in pool.map(safe, jobs):
                (failures if isinstance(result, tuple) else done).append(result)
    else:
        for result in map(safe, jobs):
            (failures if isinstance(result, tuple) else done).append(result)
    os.makedirs(opts.out, exist_ok=True)
    done.sort(key=lambda r: (r["set"], r["point"], r["sample"]))
    with open(os.path.join(opts.out, "summary.tsv"), "w", newline="") as fh:
        writer = csv.DictWriter(fh, SUMMARY_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(done)
    with open(os.path.join(opts.out, "failed.txt"), "w") as fh:
        for sample, text in failures:
            fh.write(f"== {sample}\n{text}\n")
    if failures:
        print(f"{len(failures)} samples failed (tracebacks in {os.path.join(opts.out, 'failed.txt')}): "
              f"{', '.join(s for s, _ in failures[:6])}\n{failures[0][1]}", flush=True)
    print(f"{len(done)} {opts.read_type} samples: {sum(r['FP'] for r in done)} FP and {sum(r['FN'] for r in done)} FN "
          f"taxa; {sum(r['kept_fragments'] for r in done)} fragments kept of {sum(r['fragments'] for r in done)}; "
          f"{sum(r['out_bytes'] for r in done) / 1e6:.1f} MB in {opts.out}; "
          f"{sum(r['not_consecutive'] for r in done)} reads whose records were not one after the other; "
          f"{sum(r['unknown_source'] for r in done)} fragments of unknown source; {time.time() - began:.0f} s",
          flush=True)
    return 0 if done else 1


if __name__ == "__main__":
    sys.exit(main())
