#!/usr/bin/env python3
"""The reads behind a presence model's errors in the samples it was trained and tested on, from the samples' SAMs.

build_gtdb_database.py --error-reads (default all) has the collectors profile the samples with an unmapped record for
every read that seeded on taxa but aligned nowhere (collect_training_data.py --unmapped_reads, the map's
UNMAPPED_READS: the non-hits), keeps their SAMs, and runs this script once the models are trained. The trainer's
PREFIX.calls.tsv.gz holds the call of every row at the model's knob, as protal calls by default: set "training" (the
design's samples and the scenarios' hold-in samples) with species held out, so that no model saw the taxon's species,
and set "test" (the independent test set and the scenarios' hold-out samples) by the final model. A sample's errors are
  FP      absent taxa called
  FN      present taxa not called
  unseen  species of the sample that the training database has but that protal profiled no reads to (no row: the
          model never scored them; protal's own count of false negatives includes them)
and its SAM gives every record of the reads (of paired-end reads: of the fragments) that
  - align to an FP or FN taxon (RNAME <taxid>_<gene>),
  - seeded on an FN or unseen taxon but did not align to it (its ZF tag; an unmapped record if they aligned nowhere), or
  - come from a genome of an FN or unseen species. ART names a read after its contig (<contig>-<n>; single-end reads
    are the paired-end reads' first), whose genome the FASTAs of the sample's manifest.tsv tell (trace_relatives.py);
    the collector names a drawn read (PacBio, Nanopore, Ultima) g<i>x_<n>, i the genome's place among its community's
    in the manifest (the host's after them).
Each record gains its read's source and why it was taken: xg:Z:<genome>, xs:Z:<species> (" (not in the database)" for a
species the training database lacks), xe:Z:FP:<taxid>,FN:<taxid>,seeded:<taxid>,source:<taxid> (lower-case tags: the
SAM specification leaves them to users). A read that seeded on nothing has no record: most of a genome's reads are
outside its marker genes. The simulations are seeded: the collector replays a sample's reads byte for byte.

Written to OUT/<training|test>/<design point>/:
  <sample>.sam.zst   those records, the header's @SQ lines cut to the genes they name
  <sample>.taxa.tsv  each error taxon: its score and knob; for an FN or unseen species its genomes and the reads (pairs)
                     simulated from them (paired-end samples), its fragments with a record, those whose best record is
                     on itself (at MAPQ 4 or more, as the profiler counts them), elsewhere (where, by taxon) or none; for
                     any, the fragments with a record on it and their sources, and those that seeded on it but failed to
                     align
and OUT/summary.tsv (one line per sample).

    python3 scripts/error_reads.py --calls OUT/trained_model.calls.tsv.gz --training OUT/training --test OUT/test \\
        --db OUT/training_db --read-type pe --out OUT/model_logs/error_reads/pe
"""
import argparse
import collections
import concurrent.futures
import csv
import glob
import gzip
import multiprocessing
import os
import re
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "mini_db"))
import collect_training_data as collect  # noqa: E402
import compressed  # noqa: E402
import trace_relatives  # noqa: E402

MIN_MAPQ = trace_relatives.MIN_MAPQ  # the profiler's MAPQ floor (Profiler::m_min_mapq)
SETS = ("training", "test")
DESIGN = "design"  # --samples: the design's samples (no meta_scenario)
NOT_IN_DB = " (not in the database)"
TOP = 5  # taxa listed in a taxa.tsv cell
DRAWN_NAME = re.compile(rb"g(\d+)x_")  # collect_training_data.draw_templates
TAXA_COLUMNS = ["sample", "set", "error", "taxid", "taxon_name", "p", "knob", "genomes", "read_pairs",
                "own_fragments", "own_best_on_taxon", "own_best_on_taxon_mapq4", "own_best_elsewhere", "own_unaligned",
                "own_best_elsewhere_on", "fragments_on_taxon", "best_on_taxon", "best_on_taxon_mapq4",
                "fragments_on_taxon_from", "seeded_not_aligned"]
SUMMARY_COLUMNS = ["set", "point", "scenario", "sample", "FP", "FN", "unseen", "fragments", "records", "unknown_source",
                   "sam", "sam_bytes", "source_sam_bytes", "seconds"]

# Set before the worker processes fork (main): {FASTA path: [contig names]}, {species: taxid} and {taxid: species}
# of the training database.
CONTIGS, DB_TAXID, DB_NAME = {}, {}, {}


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--calls", required=True, help="the trainer's PREFIX.calls.tsv.gz")
    p.add_argument("--training", help="the collection of the training samples (collect_training_data.py -o; "
                                      "build_gtdb_database.py: OUT/training or SCRATCH/training)")
    p.add_argument("--test", help="the collection of the test samples (OUT/test)")
    p.add_argument("--db", required=True, help="the training database's folder (genome2tiid.tsv: its species)")
    p.add_argument("--read-type", default="pe", choices=list(collect.READ_TYPES), help="the samples' read type")
    p.add_argument("--samples", default="all",
                   help="whose samples: all (default), or design and scenario names, comma-separated")
    p.add_argument("--out", required=True, help="output folder")
    p.add_argument("--threads", type=int, default=4, help="samples at once, and genome FASTAs read at once (default 4)")
    p.add_argument("--contig-cache", help="a file of the genomes' contig names read before (trace_relatives.py "
                                          "--contig-cache), joined by those read here")
    opts = p.parse_args(argv)
    opts.samples = {s.strip() for s in opts.samples.split(",") if s.strip()}
    if not opts.training and not opts.test:
        p.error("give --training or --test (or both)")
    return opts


def species_of(lineage):
    """s__Genus species, as protal names a species, of a GTDB lineage."""
    return lineage.strip().split(";")[-1].strip()


def database_species(db):
    """({species: taxid}, {taxid: species}) of a database's genome2tiid.tsv (accession, species taxid, representative,
    lineage)."""
    by_name, by_id = {}, {}
    with open(os.path.join(db, "genome2tiid.tsv")) as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) >= 4 and row[1]:
                name = species_of(row[3])
                by_name.setdefault(name, row[1])
                by_id.setdefault(row[1], name)
    return by_name, by_id


def read_calls(path, samples, read_type):
    """{(set, design point, sample): [rows]} of the calls of a read type's samples (samples: {"all"}, or design and
    scenario names)."""
    out = collections.defaultdict(list)
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("meta_read_type", read_type) != read_type or row.get("set") not in SETS:
                continue
            if "all" in samples or (row.get("meta_scenario") or DESIGN) in samples:
                out[(row["set"], row["meta_design"], row["meta_sample"])].append(row)
    return out


def sam_of(folder, sample):
    """The SAM protal wrote for a sample in a folder (.sam.zst, .sam.gz or .sam), or None."""
    found = [p for p in glob.glob(os.path.join(folder, glob.escape(sample) + ".sam*"))
             if os.path.basename(p)[len(sample):] in (".sam", ".sam.gz", ".sam.zst")]
    return found[0] if found else None


def sample_files(collection, point, sample, read_type):
    """(SAM path, manifest rows of the sample's community in the manifest's order, whether its reads are drawn) of a
    sample of a collection, as the collector lays them out; None for a sample whose SAM is not there."""
    points = os.path.join(collection, "points")
    drawn = os.path.join(points, point, "sim", "samples.tsv")
    if os.path.isfile(drawn):  # long reads and Ultima reads: drawn from a community sample's genomes
        with open(drawn) as fh:
            next(fh)
            community = next((f[3] for f in (line.rstrip("\n").split("\t") for line in fh) if f[0] == sample), None)
        sam = sam_of(os.path.join(points, point, "protal", "alignments"), sample)
        if community is None or sam is None:
            return None
        sim = os.path.join(points, community.rsplit("_s_", 1)[0], "sim")
        return sam, [r for r in collect.manifest_rows(sim) if r["sample"] == community], True
    if read_type == "se" and point.endswith("_se") and sample.endswith("_se"):  # a paired-end point's first reads
        base = point[:-len("_se")]
        sam = sam_of(os.path.join(points, base, "protal_se", "alignments"), sample)
        sim = os.path.join(points, base, "sim")
        if sam is None or not os.path.isfile(os.path.join(sim, "manifest.tsv")):
            return None
        return sam, [r for r in collect.manifest_rows(sim) if r["sample"] == sample[:-len("_se")]], False
    sim = os.path.join(points, point, "sim")
    meta = os.path.join(sim, "protal.meta")
    if not os.path.isfile(meta):
        return None
    _, rows, _ = collect.map_rows(meta)
    sam = next((r["SAM"] for r in rows if r["SAMPLEID"] == sample), None)
    if sam is None or not os.path.isfile(sam):
        return None
    return sam, [r for r in collect.manifest_rows(sim) if r["sample"] == sample], False


def errors_of(rows, manifest):
    """The errors of a sample: {taxid: (FP or FN, its row)} of its rows, and {taxid: species} of the training
    database's species in its manifest without a row (unseen)."""
    errors = {}
    for row in rows:
        truth, call = row["truth"].lower() in ("1", "true"), row["call"] == "1"
        if truth != call:
            errors[row["taxon"]] = ("FN" if truth else "FP", row)
    profiled = {row["taxon"] for row in rows}
    unseen = {}
    for species in dict.fromkeys(species_of(r.get("taxonomy", "")) for r in manifest):
        taxid = DB_TAXID.get(species)
        if taxid and taxid not in profiled:
            unseen[taxid] = species
    return errors, unseen


class Fragment:
    """What a fragment's records say: the taxa it aligned to, its best record's (primary: no flag 0x904) taxon and
    MAPQ, the taxa it seeded on but did not align to (ZF), its records."""
    __slots__ = ("aligned", "best", "best_mapq", "failed", "records")

    def __init__(self):
        self.aligned, self.best, self.best_mapq, self.failed, self.records = set(), None, 0, set(), 0


def fields_of(line):
    """(QNAME, FLAG, RNAME, MAPQ, ZF taxids) of a SAM record (bytes)."""
    f = line.split(b"\t", 5)
    zf = line.rfind(b"\tZF:Z:")
    failed = line[zf + 6:].rstrip(b"\r\n").split(b",") if zf >= 0 else ()
    return f[0], int(f[1]), f[2], int(f[4]), failed


def records(path):
    """The records of a SAM (bytes lines, the header's left out)."""
    with compressed.open_read(path) as fh:
        for line in fh:
            if not line.startswith(b"@") and line.strip():
                yield line


def source_finder(manifest, drawn):
    """qname -> (genome, species) of the read's source, or (None, None): by the contig ART named the read after, or by
    a drawn read's genome index."""
    species_of_genome = {r["genome"]: species_of(r.get("taxonomy", "")) for r in manifest}
    if drawn:
        genomes = [r["genome"] for r in manifest] + ["host"]

        def source(qname):
            m = DRAWN_NAME.match(qname)
            if not m or int(m.group(1)) >= len(genomes):
                return None, None
            genome = genomes[int(m.group(1))]
            return genome, species_of_genome.get(genome, "host")
        return source
    # The genome of each contig of the sample (a contig name two genomes share: neither).
    genome_of, shared = {}, set()
    for r in manifest:
        for contig in CONTIGS.get(r.get("fasta_path"), ()):
            key = contig.encode()
            if genome_of.setdefault(key, r["genome"]) != r["genome"]:
                shared.add(key)
    for key in shared:
        del genome_of[key]

    def source(qname):
        genome = genome_of.get(qname.rsplit(b"-", 1)[0])
        return genome, (species_of_genome.get(genome) if genome else None)
    return source


def extract(job):
    """One sample: its errors, the records of their reads into its SAM, its taxa table. -> its summary row."""
    began = time.time()
    which, point, scenario, sample, rows, sam, manifest, drawn, out_dir = job
    errors, unseen = errors_of(rows, manifest)
    reason_of = {t.encode(): e for t, (e, _) in errors.items()}  # FP or FN: the records on the taxon
    seeded = {t.encode() for t, (e, _) in errors.items() if e == "FN"} | {t.encode() for t in unseen}
    # The species whose reads are taken wherever they went: the FN and unseen ones, by name (protal's and GTDB's).
    sources = {row.get("taxon_name", ""): t for t, (e, row) in errors.items() if e == "FN"}
    sources.update({species: t for t, species in unseen.items()})
    source = source_finder(manifest, drawn)
    source_reason = {r["genome"]: "source:" + sources[species_of(r.get("taxonomy", ""))] for r in manifest
                     if species_of(r.get("taxonomy", "")) in sources}

    # Pass 1: the fragments to take, and why.
    why = collections.defaultdict(set)
    for line in records(sam):
        qname, _, rname, _, failed = fields_of(line)
        if rname != b"*":
            taxid = rname.split(b"_", 1)[0]
            if taxid in reason_of:
                why[qname].add(reason_of[taxid] + ":" + taxid.decode())
        for taxid in failed:
            if taxid in seeded:
                why[qname].add("seeded:" + taxid.decode())
        genome, _ = source(qname)
        if genome in source_reason:
            why[qname].add(source_reason[genome])

    # Pass 2: their records, tagged, into a file of their own; the genes they name; what each fragment did.
    os.makedirs(out_dir, exist_ok=True)
    target = os.path.join(out_dir, sample + ".sam.zst")
    with compressed.open_read(sam) as fh:
        header = []
        for line in fh:
            if not line.startswith(b"@"):
                break
            header.append(line)
    genes, fragments, n_records, unknown = set(), {}, 0, 0
    with tempfile.TemporaryFile(dir=out_dir) as body:
        for line in records(sam):
            qname, flag, rname, mapq, failed = fields_of(line)
            reasons = why.get(qname)
            if reasons is None:
                continue
            n_records += 1
            fragment = fragments.get(qname)
            if fragment is None:
                fragment = fragments[qname] = Fragment()
            fragment.records += 1
            fragment.failed.update(failed)
            if rname != b"*":
                genes.add(rname)
                taxid = rname.split(b"_", 1)[0]
                fragment.aligned.add(taxid)
                if not flag & 0x904 and fragment.best is None:
                    fragment.best, fragment.best_mapq = taxid, mapq
                rnext = line.split(b"\t", 7)[6]
                if rnext not in (b"=", b"*"):
                    genes.add(rnext)
            genome, species = source(qname)
            unknown += genome is None
            tags = f"\txg:Z:{genome or '?'}\txs:Z:{source_name(species)}\txe:Z:{','.join(sorted(reasons))}"
            body.write(line.rstrip(b"\r\n") + tags.encode() + b"\n")
        body.seek(0)
        kept = [h for h in header if not h.startswith(b"@SQ") or h.split(b"\t", 2)[1][3:].rstrip(b"\r\n") in genes]
        counts = collections.Counter(e for e, _ in errors.values())
        note = (f"@CO\terror_reads.py: the records of the reads of the model's errors in {sample} ({which} set, "
                f"{scenario or DESIGN}): {counts['FP']} FP, {counts['FN']} FN and {len(unseen)} unseen taxa; xg and xs: "
                f"the read's source genome and species, xe: why it was taken\n").encode()
        with compressed.open_write(target + ".partial") as out:
            out.write(b"".join(kept) + note)
            while True:
                chunk = body.read(1 << 20)
                if not chunk:
                    break
                out.write(chunk)
    os.replace(target + ".partial", target)
    write_taxa(os.path.join(out_dir, sample + ".taxa.tsv"), sample, which, errors, unseen, manifest, fragments, source,
               drawn)
    return {"set": which, "point": point, "scenario": scenario or DESIGN, "sample": sample, "FP": counts["FP"],
            "FN": counts["FN"], "unseen": len(unseen), "fragments": len(fragments), "records": n_records,
            "unknown_source": unknown, "sam": os.path.relpath(target, os.path.dirname(os.path.dirname(out_dir))),
            "sam_bytes": os.path.getsize(target), "source_sam_bytes": os.path.getsize(sam),
            "seconds": round(time.time() - began, 1)}


def source_name(species):
    """A read's source species as the outputs name it."""
    return "?" if not species else species + ("" if species in DB_TAXID or species == "host" else NOT_IN_DB)


def top(counter, name=str):
    """The most frequent keys of a counter, "name:count" comma-separated."""
    return ",".join(f"{name(k)}:{n}" for k, n in counter.most_common(TOP))


def write_taxa(path, sample, which, errors, unseen, manifest, fragments, source, drawn):
    """<sample>.taxa.tsv: a line per error taxon (TAXA_COLUMNS), FP first, then FN and unseen, by name."""
    def taxon_name(taxid):
        text = taxid.decode()
        return DB_NAME.get(text, text)

    genomes_of = collections.defaultdict(list)
    for r in manifest:
        genomes_of[species_of(r.get("taxonomy", ""))].append(r)
    own = collections.defaultdict(list)  # source species -> its fragments
    on = collections.defaultdict(list)   # taxid -> (source species, fragment) of the fragments aligned to it
    failed = collections.Counter()       # taxid -> fragments that seeded on it but did not align to it
    for qname, fragment in fragments.items():
        _, species = source(qname)
        own[species].append(fragment)
        for taxid in fragment.aligned:
            on[taxid].append((species, fragment))
        for taxid in fragment.failed - fragment.aligned:
            failed[taxid] += 1
    entries = [(taxid, error, row.get("taxon_name", ""), row.get("p", ""), row.get("knob", ""))
               for taxid, (error, row) in errors.items()]
    entries += [(taxid, "unseen", species, "", "") for taxid, species in unseen.items()]
    order = {"FP": 0, "FN": 1, "unseen": 2}
    with open(path + ".partial", "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(TAXA_COLUMNS)
        for taxid, error, name, p, knob in sorted(entries, key=lambda e: (order[e[1]], e[2])):
            key = taxid.encode()
            on_taxon = on.get(key, [])
            best_on = [f for _, f in on_taxon if f.best == key]
            row = [sample, which, error, taxid, name, p, knob]
            if error == "FP":
                row += [""] * 8
            else:
                genomes, mine = genomes_of.get(name, []), own.get(name, [])
                best_here = [f for f in mine if f.best == key]
                elsewhere = collections.Counter(f.best for f in mine if f.best is not None and f.best != key)
                pairs = "" if drawn else sum(int(r.get("read_pairs") or 0) for r in genomes)
                row += [",".join(r["genome"] for r in genomes), pairs, len(mine), len(best_here),
                        sum(f.best_mapq >= MIN_MAPQ for f in best_here), sum(elsewhere.values()),
                        sum(f.best is None for f in mine), top(elsewhere, taxon_name)]
            row += [len(on_taxon), len(best_on), sum(f.best_mapq >= MIN_MAPQ for f in best_on),
                    top(collections.Counter(source_name(s) for s, _ in on_taxon)), failed.get(key, 0)]
            writer.writerow(row)
    os.replace(path + ".partial", path)


def main(argv=None):
    opts = parse_args(argv)
    began = time.time()
    by_name, by_id = database_species(opts.db)
    DB_TAXID.update(by_name)
    DB_NAME.update(by_id)
    calls = read_calls(opts.calls, opts.samples, opts.read_type)
    collections_ = {"training": opts.training, "test": opts.test}
    jobs, missing = [], []
    for (which, point, sample), rows in sorted(calls.items()):
        collection = collections_[which]
        files = sample_files(collection, point, sample, opts.read_type) if collection else None
        if files is None:
            missing.append(f"{sample} ({which})")
            continue
        sam, manifest, drawn = files
        jobs.append((which, point, rows[0].get("meta_scenario") or "", sample, rows, sam, manifest, drawn,
                     os.path.join(opts.out, which, point)))
    if missing:
        print(f"no SAM for {len(missing)} samples (not collected here, or removed): {', '.join(missing[:6])}"
              f"{' ...' if len(missing) > 6 else ''}", flush=True)
    if not jobs:
        print(f"no {opts.read_type} samples in {opts.calls} with their SAMs", flush=True)
        return 0
    fastas = [r.get("fasta_path") for job in jobs if not job[7] for r in job[6]]
    CONTIGS.update(trace_relatives.genome_contigs(fastas, opts.threads, opts.contig_cache))
    print(f"the contigs of {len(CONTIGS)} genomes in {time.time() - began:.0f} s", flush=True)
    jobs.sort(key=lambda job: -os.path.getsize(job[5]))  # the largest SAMs first, so that the workers end together
    methods = multiprocessing.get_all_start_methods()
    if opts.threads > 1 and len(jobs) > 1 and "fork" in methods:  # the workers share CONTIGS as forked
        with concurrent.futures.ProcessPoolExecutor(min(opts.threads, len(jobs)),
                                                    mp_context=multiprocessing.get_context("fork")) as pool:
            done = list(pool.map(extract, jobs))
    else:
        done = [extract(job) for job in jobs]
    done.sort(key=lambda r: (SETS.index(r["set"]), r["point"], r["sample"]))
    os.makedirs(opts.out, exist_ok=True)
    with open(os.path.join(opts.out, "summary.tsv"), "w", newline="") as fh:
        writer = csv.DictWriter(fh, SUMMARY_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(done)
    by_group = collections.defaultdict(collections.Counter)
    for row in done:
        by_group[(row["set"], row["scenario"])].update(samples=1, FP=row["FP"], FN=row["FN"], unseen=row["unseen"],
                                                        fragments=row["fragments"], bytes=row["sam_bytes"],
                                                        unknown=row["unknown_source"], records=row["records"])
    for (which, scenario), c in sorted(by_group.items()):
        print(f"{opts.read_type} {which} {scenario}: {c['samples']} samples, {c['FP']} FP, {c['FN']} FN, {c['unseen']} "
              f"unseen; {c['fragments']} fragments, {c['records']} records ({c['unknown']} of unknown source), "
              f"{c['bytes'] / 1e6:.1f} MB", flush=True)
    print(f"error reads of {len(done)} {opts.read_type} samples: {sum(r['FP'] for r in done)} FP, "
          f"{sum(r['FN'] for r in done)} FN and {sum(r['unseen'] for r in done)} unseen taxa, "
          f"{sum(r['fragments'] for r in done)} fragments, {sum(r['sam_bytes'] for r in done) / 1e6:.1f} MB in "
          f"{opts.out}, {time.time() - began:.0f} s", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
