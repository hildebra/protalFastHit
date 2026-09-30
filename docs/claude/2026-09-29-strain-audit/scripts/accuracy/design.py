#!/usr/bin/env python3
"""design.py - write simulate_metagenomes --from_manifest manifests for the accuracy experiment.

Run A (depth):   42 samples x 8 species. Species k in sample j gets genome index (j//7 + k) % 6
                 (0 = the representative itself, the control; 1..5 = strains GCA_999xxx002..006)
                 at depth D[(j + k) % 7], D = 1,2,3,5,10,20,50x. Every (genome, depth) pair of
                 every species occurs exactly once; one strain per species per sample.
Run B (mix):     12 samples x 8 species, two strains of each species per sample. Minor fraction
                 index (m + k) % 6 over 0.5,0.3,0.2,0.1,0.05,0.02; total depth 40x (m<6) or 20x.
                 Major strain = strain (m + k) % 5 + 2, minor = next strain.
Run C (frag):    7 samples x 8 species, one strain each at D[(j + k) % 7], strain (j + k) % 5 + 2;
                 the same manifest is simulated with short (200+-20) and long (600+-50) fragments.
Read pairs = depth * genome_length / (2 * 150).
"""
import os, sys

W = os.path.expanduser("~/audit5/world/gtdb_r226")
OUT = os.path.expanduser("~/audit5/accuracy/design")
DEPTHS = [1, 2, 3, 5, 10, 20, 50]
MINOR = [0.5, 0.3, 0.2, 0.1, 0.05, 0.02]
RL = 150


def load():
    g = {}
    with open(os.path.join(W, "simulation/genomes.tsv")) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = dict(zip(hdr, line.rstrip("\n").split("\t")))
            f["species"] = f["gtdb_taxonomy"].split(";s__")[-1]
            g[f["accession"]] = f
    species = sorted({f["species"] for f in g.values()}, key=lambda s: min(a for a in g if g[a]["species"] == s))
    by_sp = {s: sorted(a for a in g if g[a]["species"] == s) for s in species}
    # rep first (GCF_...001), then GCA ...002..006
    for s in by_sp:
        by_sp[s].sort(key=lambda a: (not a.startswith("GCF"), a))
    return g, species, by_sp


COLS = ["sample", "genome", "species", "taxonomy", "genome_length", "read_pairs", "vertical_coverage",
        "relative_abundance", "fasta_path", "art_seed"]


def write(path, rows):
    # relative abundance by read pairs within the sample
    tot = {}
    for r in rows:
        tot[r["sample"]] = tot.get(r["sample"], 0) + r["read_pairs"]
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        for r in rows:
            r["relative_abundance"] = "%.6g" % (r["read_pairs"] / tot[r["sample"]])
            fh.write("\t".join(str(r[c]) for c in COLS) + "\n")


def row(g, sample, acc, depth, seed):
    f = g[acc]
    L = int(f["genome_length"])
    rp = max(1, int(round(depth * L / (2 * RL))))
    return dict(sample=sample, genome=acc, species=f["species"], taxonomy=f["gtdb_taxonomy"], genome_length=L,
                read_pairs=rp, vertical_coverage="%.4f" % (rp * 2 * RL / L), fasta_path=f["fasta_path"],
                art_seed=seed)


def main():
    g, species, by_sp = load()
    seed = 1000
    # Run A
    rows = []
    for j in range(42):
        for k, s in enumerate(species):
            acc = by_sp[s][(j // 7 + k) % 6]
            seed += 1
            rows.append(row(g, "A%02d" % j, acc, DEPTHS[(j + k) % 7], seed))
    write(os.path.join(OUT, "A.manifest.tsv"), rows)
    # Run B
    rows = []
    for m in range(12):
        total = 40 if m < 6 else 20
        for k, s in enumerate(species):
            fmin = MINOR[(m + k) % 6]
            major = by_sp[s][1 + (m + k) % 5]
            minor = by_sp[s][1 + (m + k + 1) % 5]
            seed += 1
            rows.append(row(g, "B%02d" % m, major, total * (1 - fmin), seed))
            seed += 1
            rows.append(row(g, "B%02d" % m, minor, total * fmin, seed))
    write(os.path.join(OUT, "B.manifest.tsv"), rows)
    # Run C
    rows = []
    for j in range(7):
        for k, s in enumerate(species):
            acc = by_sp[s][1 + (j + k) % 5]
            seed += 1
            rows.append(row(g, "C%02d" % j, acc, DEPTHS[(j + k) % 7], seed))
    write(os.path.join(OUT, "C.manifest.tsv"), rows)
    print("species order:", species, file=sys.stderr)


if __name__ == "__main__":
    main()
