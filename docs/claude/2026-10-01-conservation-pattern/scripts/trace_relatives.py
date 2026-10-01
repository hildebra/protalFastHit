#!/usr/bin/env python3
"""trace_relatives.py [BENCH_DIR] - where the reads of a species the database lacks land, by how fast the gene evolves.

On the v0.7.1 benchmark world (BENCH_DIR, default ~/bench071): the 0.7.1 pipeline's paired-end training samples
(V071/training/points/rl*), which protal aligned to the training database (V071/training_db: 290 of the 765 species
held out, 119 of them as single species whose genus stays). A read's name starts with the accession of the genome it
was simulated from; its primary records name the taxon and marker gene they align to (RNAME <taxid>_<geneid>). The
gene's true rate is the simulator's (world/full/simulation/gene_rates.tsv: the factor that scales the divergence of
every branch at that marker, within and between species alike).

For each gene, per unit of its source genomes' coverage (manifest.tsv), the primary records of
  own      reads of species in the database, on their own species (the baseline: alignability, gene length)
  relative reads of species held out at species rank (their genus is in the database), on any taxon
and of the relative's records: the share on a congener of the source, the share with MAPQ below 4 (dropped by the
profiler, Profiler::m_min_mapq) and below 10 (low_mapq_share), and the mean identity. The relative's records per unit
coverage over the own ones' (R) says on which genes the relative's reads land, relative to how the gene takes reads
anyway. Writes BENCH_DIR/results_conservation/trace_relatives.md and genes.tsv.
"""
import collections
import csv
import glob
import gzip
import io
import os
import re
import statistics
import subprocess
import sys

B = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/bench071")
W = os.path.join(B, "world", "full", "simulation")
DB = os.path.join(B, "V071", "training_db")
OUT = os.path.join(B, "results_conservation")
MIN_MAPQ, LOW_MAPQ = 4, 10
CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def species_of(lineage):
    return lineage.strip().split(";")[-1].removeprefix("s__")


def genus_of(lineage):
    return lineage.strip().split(";")[5]


def read_sam(path):
    proc = subprocess.Popen(["zstd", "-dc", path], stdout=subprocess.PIPE)
    for line in io.TextIOWrapper(proc.stdout, encoding="ascii"):
        if not line.startswith("@"):
            yield line.split("\t", 6)
    proc.wait()


def identity(cigar):
    m = d = 0
    for n, op in CIGAR.findall(cigar):
        if op in "M=":
            m += int(n)
        elif op in "XID":
            d += int(n)
    return m / (m + d) if m + d else 0.0


def main():
    heldout = {}
    with open(os.path.join(B, "V071", "heldout_species.txt")) as fh:
        for line in fh:
            name, rank = line.rstrip("\n").split("\t")[:2]
            heldout[name.removeprefix("s__")] = rank
    taxon_lineage = {}  # training-db taxid -> lineage
    with open(os.path.join(DB, "genome2tiid.tsv")) as fh:
        for row in csv.reader(fh, delimiter="\t"):
            taxon_lineage[row[1]] = row[3]
    marker = {}
    with open(os.path.join(DB, "gene2geneid.tsv")) as fh:
        for row in csv.reader(fh, delimiter="\t"):
            marker[row[1]] = row[0]
    rate = {}
    with open(os.path.join(W, "gene_rates.tsv")) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            rate[row["marker"]] = float(row["rate"])

    coverage = collections.Counter()  # class -> summed source coverage
    counts = collections.defaultdict(collections.Counter)  # (class, gene) -> counter
    idsum = collections.defaultdict(float)
    spread = collections.defaultdict(collections.Counter)  # (sample, relative genome, gene) -> {taxid: kept records}
    samples = 0
    for point in sorted(glob.glob(os.path.join(B, "V071", "training", "points", "rl*"))):
        manifest = os.path.join(point, "sim", "manifest.tsv")
        with open(manifest) as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        for sam in sorted(glob.glob(os.path.join(point, "protal", "alignments", "*.sam.zst"))):
            sample = os.path.basename(sam)[:-len(".sam.zst")]
            samples += 1
            source_class = {}
            for r in rows:
                if r["sample"] != sample:
                    continue
                sp = species_of(r["taxonomy"])
                cls = "relative" if heldout.get(sp) == "species" else ("own" if sp not in heldout else None)
                if cls:
                    source_class[r["genome"]] = (cls, r["taxonomy"])
                    coverage[cls] += float(r["vertical_coverage"])
            for f in read_sam(sam):
                flag = int(f[1])
                if flag & 0x904 or f[2] == "*":
                    continue
                acc = f[0].split("_contig")[0]
                if acc not in source_class:
                    continue
                cls, src = source_class[acc]
                taxid, gene = f[2].split("_", 1)
                hit = taxon_lineage.get(taxid)
                if hit is None:
                    continue
                own = species_of(hit) == species_of(src)
                if cls == "own" and not own:
                    continue
                key = (cls, gene)
                c = counts[key]
                mapq = int(f[4])
                c["records"] += 1
                c["congener"] += genus_of(hit) == genus_of(src) and not own
                c["mapq_lt4"] += mapq < MIN_MAPQ
                c["mapq_lt10"] += mapq < LOW_MAPQ
                if mapq >= MIN_MAPQ:
                    c["kept"] += 1
                    if cls == "relative":
                        spread[(sample, acc, gene)][taxid] += 1
                idsum[key] += identity(f[5])

    top_share = collections.defaultdict(list)  # gene -> share of a relative's kept records on its most-hit taxon
    for (_, _, gene), hits in spread.items():
        top_share[gene].append(max(hits.values()) / sum(hits.values()))
    os.makedirs(OUT, exist_ok=True)
    genes = sorted({g for _, g in counts}, key=lambda g: rate.get(marker.get(g, ""), 0))
    table = []
    for g in genes:
        o, r = counts[("own", g)], counts[("relative", g)]
        if not o["records"]:
            continue
        e_own = o["records"] / coverage["own"]
        e_rel = r["records"] / coverage["relative"]
        table.append(dict(gene=g, marker=marker.get(g, "?"), rate=rate.get(marker.get(g, ""), float("nan")),
                          own_per_cov=e_own, relative_per_cov=e_rel, R=e_rel / e_own,
                          R_kept=((r["kept"] / coverage["relative"]) / (o["kept"] / coverage["own"])
                                  if o["kept"] else float("nan")),
                          top_taxon_share=statistics.mean(top_share[g]) if top_share[g] else float("nan"),
                          congener_share=r["congener"] / r["records"] if r["records"] else float("nan"),
                          mapq_lt4=r["mapq_lt4"] / r["records"] if r["records"] else float("nan"),
                          mapq_lt10=r["mapq_lt10"] / r["records"] if r["records"] else float("nan"),
                          own_mapq_lt4=o["mapq_lt4"] / o["records"],
                          relative_identity=idsum[("relative", g)] / r["records"] if r["records"] else float("nan"),
                          own_identity=idsum[("own", g)] / o["records"]))
    with open(os.path.join(OUT, "genes.tsv"), "w") as fh:
        cols = list(table[0])
        fh.write("\t".join(cols) + "\n")
        for t in table:
            fh.write("\t".join(f"{t[c]:.4g}" if isinstance(t[c], float) else str(t[c]) for c in cols) + "\n")

    def med(rows, key):
        vals = [r[key] for r in rows if r[key] == r[key]]
        return statistics.median(vals) if vals else float("nan")

    classes = [("rate < 0.7", lambda x: x < 0.7), ("0.7-1", lambda x: 0.7 <= x < 1), ("1-1.4", lambda x: 1 <= x < 1.4),
               ("rate >= 1.4", lambda x: x >= 1.4)]
    lines = ["# Where a missing species' reads land, by the gene's rate (trace_relatives.py)", "",
             f"{samples} paired-end training samples of the 0.7.1 pipeline; source coverage summed: own "
             f"{coverage['own']:.1f}, relative {coverage['relative']:.1f}. Medians over the genes of each class.", "",
             "| gene rate | genes | own records / coverage | relative records / coverage | R (relative / own) | "
             "on a congener | MAPQ < 4 | MAPQ < 10 | own MAPQ < 4 | R of the records kept (MAPQ >= 4) | "
             "kept records on the most-hit taxon | relative identity | own identity |",
             "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for label, test in classes:
        rows = [t for t in table if test(t["rate"])]
        if not rows:
            continue
        lines.append(f"| {label} | {len(rows)} | " + " | ".join(
            f"{med(rows, k):.3f}" for k in ("own_per_cov", "relative_per_cov", "R", "congener_share", "mapq_lt4",
                                           "mapq_lt10", "own_mapq_lt4", "R_kept", "top_taxon_share",
                                           "relative_identity", "own_identity")) + " |")
    xs = [t["rate"] for t in table if t["R"] == t["R"]]
    ys = [t["R"] for t in table if t["R"] == t["R"]]
    if len(xs) > 2:
        lines += ["", f"Spearman correlation of the gene's rate with R over {len(xs)} genes: "
                      f"{spearman(xs, ys):+.3f}"]
    text = "\n".join(lines) + "\n"
    with open(os.path.join(OUT, "trace_relatives.md"), "w") as fh:
        fh.write(text)
    print(text)


def spearman(x, y):
    def ranks(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        for rank, i in enumerate(order):
            r[i] = rank
        return r
    rx, ry = ranks(x), ranks(y)
    mx, my = statistics.mean(rx), statistics.mean(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = (sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry)) ** 0.5
    return num / den if den else float("nan")


if __name__ == "__main__":
    main()
