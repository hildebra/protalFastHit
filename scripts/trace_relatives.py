#!/usr/bin/env python3
"""Where the reads of species the training database lacks land, by how fast the genes diverge.

build_gtdb_database.py simulates its training samples from genomes of every species, and aligns them to a training
database that leaves some species out. The reads of a species left out (whose genus stays) can only land on its
congeners. This script follows them read by read: a simulated read is named after the contig it came from
(<contig>-<n>), whose genome the genomes' FASTAs in the samples' manifest.tsv tell (their headers, read once per
genome; GTDB's contigs are named by NCBI, not by their genome), and its primary record names the taxon and gene it
aligned to (RNAME <taxid>_<geneid>). With the genomes being real (GTDB's), it tells how a relative's reads spread over the genes of
real congeners, by the genes' conservation factors (gene_conservation.tsv, here from the training database's
gene_congeners.tsv), which the model's conservation features rest on (docs/claude/2026-10-01-conservation-pattern).

    python3 scripts/trace_relatives.py --points OUT/training/points --db OUT/training_db \\
        --heldout OUT/heldout_species.txt --out OUT/model_logs/relatives_by_gene_conservation

For each gene, per unit of its source genomes' coverage (the samples' manifest.tsv), the primary records of
  own       reads of species in the database, on their own species (the baseline of how the gene takes reads)
  relative  reads of species held out at species rank, on any taxon
and R, the relative's records per coverage over the own ones' (so that gene length and how well a gene takes reads
cancel), before and after the profiler's MAPQ filter (MAPQ 4: Profiler::m_min_mapq); the share of the relative's
records on a congener of its source, below MAPQ 4, and on the taxon its genome's reads of that gene hit most. Paired-
end samples only: single-end samples are their first reads, and long reads (pbsim) are not named by their genome.

Writes OUT.txt (the medians by class of factor, and the rank correlation of the factor with R) and OUT.tsv (each
gene). A sample past --max-records records is read only that far.
"""
import argparse
import collections
import concurrent.futures
import contextlib
import csv
import fcntl
import glob
import gzip
import io
import os
import statistics
import subprocess
import sys

MIN_MAPQ = 4  # Profiler::m_min_mapq
CLASSES = [("factor < 0.7", 0, 0.7), ("0.7 - 1", 0.7, 1.0), ("1 - 1.4", 1.0, 1.4), ("factor >= 1.4", 1.4, float("inf"))]
SHARED = object()  # a contig name of genomes of two species in one sample


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--points", required=True, help="the collector's points folder (collect_training_data.py -o DIR: DIR/points)")
    p.add_argument("--db", required=True, help="the training database's folder: genome2tiid.tsv, and gene_congeners.tsv "
                                               "or gene_conservation.tsv")
    p.add_argument("--heldout", required=True, help="heldout_species.txt: species <tab> rank it was held out at")
    p.add_argument("--out", required=True, help="output prefix: PREFIX.txt and PREFIX.tsv")
    p.add_argument("--max-records", type=int, default=5_000_000, help="records read per sample at most (default 5M)")
    p.add_argument("--threads", type=int, default=8, help="genome FASTAs read at once for their contig names (default 8)")
    p.add_argument("--contig-cache", help="a file of the genomes' contig names read before, joined by those read here "
                                          "(error_reads.py reads it too)")
    return p.parse_args(argv)


def species_of(lineage):
    return lineage.strip().split(";")[-1].strip().removeprefix("s__")


def genus_of(lineage):
    parts = lineage.strip().split(";")
    return parts[5] if len(parts) > 5 else ""


def factors_of(db):
    """geneid -> conservation factor, from gene_congeners.tsv (within_factor) or gene_conservation.tsv."""
    for name, column in (("gene_congeners.tsv", "within_factor"), ("gene_conservation.tsv", "factor")):
        path = os.path.join(db, name)
        if os.path.isfile(path):
            with open(path) as fh:
                rows = csv.DictReader((line for line in fh if not line.startswith("#")), delimiter="\t")
                return {r["geneid"]: float(r[column]) for r in rows if r.get(column) not in (None, "", "NA")}
    return {}


def sam_lines(path):
    """The lines of a SAM file (.sam, .sam.gz, or .sam.zst through the zstd command)."""
    if path.endswith(".zst"):
        process = subprocess.Popen(["zstd", "-q", "-dc", path], stdout=subprocess.PIPE)
        try:
            yield from io.TextIOWrapper(process.stdout, encoding="ascii", errors="replace")
        finally:
            process.kill()
            process.wait()
        return
    with (gzip.open(path, "rt", encoding="ascii", errors="replace") if path.endswith(".gz")
          else open(path, encoding="ascii", errors="replace")) as fh:
        yield from fh


def paired_samples(points):
    """(manifest rows of the point, sample, SAM path) of the paired-end samples under `points`, but those of the
    collector's scenarios (sc_*: deep samples, host reads not named by a contig)."""
    for manifest in sorted(glob.glob(os.path.join(points, "*", "sim", "manifest.tsv"))):
        point = os.path.dirname(os.path.dirname(manifest))
        if os.path.basename(point).startswith("sc_"):
            continue
        with open(manifest) as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        for sam in sorted(glob.glob(os.path.join(point, "protal", "alignments", "*.sam*"))):
            if sam.endswith(".err"):
                continue
            sample = os.path.basename(sam).split(".sam")[0]
            yield [r for r in rows if r.get("sample") == sample], sample, sam


def contig_names(path):
    """The first words of a genome FASTA's headers (its contigs, which name the reads), .gz or plain; an empty
    list for a file that cannot be read."""
    try:
        opener = gzip.open if path.endswith(".gz") else open
        with opener(path, "rb") as fh:
            return [line[1:].split(None, 1)[0].decode(errors="replace") for line in fh if line.startswith(b">") and
                    len(line) > 2]
    except OSError:
        return []


def genome_contigs(paths, threads, cache=None):
    """{FASTA path: [contig names]} of the genomes, read on `threads` processes. cache: a file (TSV, gzipped) of the
    names read before, by path, size and modification time, which the names read here join, so that a build reads its
    ~40,000 genomes once for trace_relatives.py and error_reads.py."""
    paths = sorted(set(p for p in paths if p))
    if not paths:
        return {}
    known, stamps = {}, {}
    for path in paths:
        try:
            st = os.stat(path)
            stamps[path] = f"{st.st_size}:{st.st_mtime_ns}"
        except OSError:
            stamps[path] = ""
    if cache and os.path.isfile(cache):
        with cache_lock(cache, fcntl.LOCK_SH), gzip.open(cache, "rt") as fh:
            for line in fh:
                path, stamp, names = (line.rstrip("\n").split("\t") + ["", ""])[:3]
                if stamps.get(path) == stamp and stamp:
                    known[path] = names.split(" ") if names else []  # a contig's name has no blank
    todo = [p for p in paths if p not in known]
    if todo:
        with concurrent.futures.ProcessPoolExecutor(max(1, min(threads, len(todo)))) as pool:
            known.update(zip(todo, pool.map(contig_names, todo, chunksize=16)))
        if cache:
            # One gzip member written at once under the lock: runs side by side (build_gtdb_database.py's reports) add
            # to the cache without mixing their writes, nor reading one half written.
            os.makedirs(os.path.dirname(os.path.abspath(cache)), exist_ok=True)
            text = "".join(f"{path}\t{stamps[path]}\t{' '.join(known[path])}\n" for path in todo if stamps[path])
            with cache_lock(cache, fcntl.LOCK_EX), open(cache, "ab") as fh:
                fh.write(gzip.compress(text.encode(), compresslevel=1))
    return {p: known[p] for p in paths}


@contextlib.contextmanager
def cache_lock(cache, how):
    """The contig cache's lock (CACHE.lock, flock): shared to read the cache, exclusive to add to it."""
    with open(cache + ".lock", "a") as fh:
        fcntl.flock(fh, how)
        try:
            yield
        finally:
            fcntl.flock(fh, fcntl.LOCK_UN)


def read_contig(name):
    """The contig a simulated read came from, by its name (<contig>-<n>; the aligner cut a mate's /1 or /2)."""
    return name.rsplit("-", 1)[0]


def trace(opts):
    factor = factors_of(opts.db)
    if not factor:
        return None, "the training database has no gene conservation factors (gene_congeners.tsv, gene_conservation.tsv)"
    heldout = {}
    with open(opts.heldout) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2 and fields[0]:
                heldout[fields[0].removeprefix("s__")] = fields[1]
    taxon_lineage = {}
    with open(os.path.join(opts.db, "genome2tiid.tsv")) as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) >= 4:
                taxon_lineage[row[1]] = row[3]
    coverage = collections.Counter()
    counts = collections.defaultdict(collections.Counter)  # (class, gene) -> counts
    spread = collections.defaultdict(collections.Counter)  # (sample, relative genome, gene) -> {taxid: kept records}
    samples = 0
    paired = list(paired_samples(opts.points))
    contigs = genome_contigs([r.get("fasta_path") for rows, _, _ in paired for r in rows], opts.threads,
                             opts.contig_cache)
    unnamed = ambiguous = 0
    for rows, sample, sam in paired:
        source, contig_source = {}, {}
        for r in rows:
            sp = species_of(r.get("taxonomy") or r.get("species", ""))
            cls = "relative" if heldout.get(sp) == "species" else ("own" if sp not in heldout else None)
            if cls:
                source[r["genome"]] = (cls, r.get("taxonomy", ""))
                coverage[cls] += float(r.get("vertical_coverage") or 0)
            for contig in contigs.get(r.get("fasta_path"), []):
                # Every genome of the sample (those of species held out at other ranks are not compared: class
                # None). A contig name two genomes share (an in-silico strain of an earlier insilico_strains.py
                # kept its representative's) is kept if both are of one species, else no read of it is traced.
                other = contig_source.setdefault(contig, (cls, r.get("taxonomy", ""), r["genome"]))
                if other is not SHARED and species_of(other[1]) != sp:
                    contig_source[contig] = SHARED
        if not source:
            continue
        samples += 1
        lines = sam_lines(sam)
        for n, line in enumerate(lines):
            if n >= opts.max_records:
                break
            if line.startswith("@"):
                continue
            f = line.split("\t", 5)
            if len(f) < 5 or int(f[1]) & 0x904 or f[2] == "*":
                continue
            if contig_source:
                contig = read_contig(f[0])
                found = contig_source.get(contig)
                if found is None or found is SHARED:
                    unnamed += found is None
                    ambiguous += found is SHARED
                    continue
                cls, src, acc = found
                if cls is None:
                    continue
            else:  # a manifest without FASTA paths: simulate_gtdb_release.py names its contigs <accession>_contigN
                acc = f[0].split("_contig")[0]
                if acc not in source:
                    unnamed += 1
                    continue
                cls, src = source[acc]
            taxid, _, gene = f[2].partition("_")
            hit = taxon_lineage.get(taxid)
            if hit is None:
                continue
            own = species_of(hit) == species_of(src)
            if cls == "own" and not own:
                continue
            c = counts[(cls, gene)]
            kept = int(f[4]) >= MIN_MAPQ
            c["records"] += 1
            c["kept"] += kept
            c["congener"] += genus_of(hit) == genus_of(src) and not own
            if cls == "relative" and kept:
                spread[(sample, acc, gene)][taxid] += 1
        lines.close()
    if not coverage["own"] or not coverage["relative"]:
        return None, (f"{samples} paired-end samples, but no reads of species held out at species rank and of species "
                      "in the database to compare")
    if not counts:
        return None, (f"{samples} paired-end samples, but none of their aligned records was traced to a genome of "
                      f"them ({unnamed} records of no contig of their genomes, {ambiguous} of contigs two species share)")
    top = collections.defaultdict(list)
    for (_, _, gene), hits in spread.items():
        top[gene].append(max(hits.values()) / sum(hits.values()))
    genes = []
    for gene in sorted({g for _, g in counts}, key=lambda g: (factor.get(g, 0), g)):
        o, r = counts[("own", gene)], counts[("relative", gene)]
        if gene not in factor or not o["records"] or not o["kept"]:
            continue
        per = lambda c, key, cls: c[key] / coverage[cls]
        share = lambda key: r[key] / r["records"] if r["records"] else float("nan")
        genes.append(dict(geneid=gene, factor=factor[gene], own_records=o["records"], relative_records=r["records"],
                          R=per(r, "records", "relative") / per(o, "records", "own"),
                          R_kept=per(r, "kept", "relative") / per(o, "kept", "own"),
                          relative_on_congener=share("congener"), relative_mapq_below_4=1 - share("kept"),
                          own_mapq_below_4=1 - o["kept"] / o["records"],
                          relative_on_top_taxon=statistics.mean(top[gene]) if top[gene] else float("nan")))
    return dict(samples=samples, coverage=coverage, genes=genes, unnamed=unnamed, ambiguous=ambiguous), None


def median(values):
    values = [v for v in values if v == v]
    return statistics.median(values) if values else float("nan")


def spearman(x, y):
    def ranks(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(v):
            j = i
            while j + 1 < len(v) and v[order[j + 1]] == v[order[i]]:
                j += 1
            for k in range(i, j + 1):
                r[order[k]] = (i + j) / 2
            i = j + 1
        return r
    if len(x) < 3:
        return float("nan")
    rx, ry = ranks(x), ranks(y)
    mx, my = statistics.mean(rx), statistics.mean(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = (sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry)) ** 0.5
    return num / den if den else float("nan")


def report(result):
    genes = result["genes"]
    cov = result["coverage"]
    lines = ["Where the reads of species held out of the training database land, by the gene's conservation factor",
             "(trace_relatives.py; docs/claude/2026-10-01-conservation-pattern)", "",
             f"{result['samples']} paired-end training samples; source coverage summed: species in the database "
             f"{cov['own']:.1f}, species held out at species rank {cov['relative']:.1f}; records of no contig of the "
             f"sample's genomes {result.get('unnamed', 0)}, of contigs two species share {result.get('ambiguous', 0)}. "
             "R: the held-out species' "
             "records per unit coverage over a species' own on the same gene; medians over the genes of each class.", "",
             "| gene factor | genes | R, all records | R, MAPQ >= 4 | on a congener | MAPQ < 4 | own MAPQ < 4 | "
             "kept on the most-hit taxon |", "|---|---|---|---|---|---|---|---|"]
    for label, low, high in CLASSES:
        rows = [g for g in genes if low <= g["factor"] < high]
        if rows:
            lines.append(f"| {label} | {len(rows)} | " + " | ".join(
                f"{median(g[k] for g in rows):.3f}" for k in ("R", "R_kept", "relative_on_congener",
                                                              "relative_mapq_below_4", "own_mapq_below_4",
                                                              "relative_on_top_taxon")) + " |")
    lines += ["", f"Spearman correlation of the factor with R over {len(genes)} genes: "
                  f"{spearman([g['factor'] for g in genes], [g['R'] for g in genes]):+.3f} (all records), "
                  f"{spearman([g['factor'] for g in genes], [g['R_kept'] for g in genes]):+.3f} (MAPQ >= 4)"]
    return "\n".join(lines) + "\n"


def main(argv=None):
    opts = parse_args(argv)
    result, problem = trace(opts)
    os.makedirs(os.path.dirname(os.path.abspath(opts.out)), exist_ok=True)
    if problem:
        text = "Not traced: " + problem + "\n"
    else:
        text = report(result)
        with open(opts.out + ".tsv", "w") as fh:
            columns = list(result["genes"][0]) if result["genes"] else ["geneid"]
            fh.write("\t".join(columns) + "\n")
            for g in result["genes"]:
                fh.write("\t".join(f"{g[c]:.4g}" if isinstance(g[c], float) else str(g[c]) for c in columns) + "\n")
    with open(opts.out + ".txt", "w") as fh:
        fh.write(text)
    print(text, end="")


if __name__ == "__main__":
    sys.exit(main())
