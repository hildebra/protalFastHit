#!/usr/bin/env python3
"""make_strains.py [BENCH_DIR] - the strain benchmark's samples: strains of known tree and known SNPs.

From the v0.7.1 benchmark world (BENCH_DIR, default ~/bench071: world/full/simulation), SPECIES species of the
databases (10 bacterial, 2 archaeal, none of the 135 unknown ones), each grown into STRAINS strains from its
representative genome (the databases' reference) along a random coalescent tree (Kingman): root-to-tip depth
0.4% to 1.2% substitutions per site at a gene of rate 1, every marker gene at its rate (world gene_rates.tsv:
ribosomal proteins slowest), substitutions only, so that a strain's genes keep the representative's coordinates and
its SNPs against the reference are known exactly. Sample j holds strain j of every species (8 samples: an MSA of 8
rows per species, whose true tree is the strain tree), at a coverage of the species' own (2 to 40x, log-uniform) times 0.4-1.6
per sample, and 15 background species of the world (any, unknown ones included) at 0.5-5x.

Reads: paired-end 2x150 (art_illumina HSXt, 350 +- 50 bp fragments); PacBio HiFi-like (pbsim3 errhmm
ERRHMM-SEQUEL, 15 +- 3 kb, 99.9%) and Nanopore (pbsim3 qshmm QSHMM-ONT-HQ, 8 +- 6 kb, 97%, 39:24:36) of the same
genomes at the same coverages. pbsim3's errhmm model gives every PacBio base quality 0 ('!'), which protal's SNP
filters reject; the PacBio reads get Q30 ('?', the accuracy of 99.9% they are simulated at) instead, as HiFi reads
have qualities of Q20 and more (make_strains.py --pb-qualities only rewrites existing ones so). Writes BENCH_DIR/strains/: genomes/, trees.tsv (species, Newick with branch lengths),
strains.tsv (species, taxonomy, representative, strain, sample, genome, coverage), reads/<sample>_R1.fq.gz, _R2,
_pb.fq.gz, _ont.fq.gz.
"""
import csv
import glob
import math
import gzip
import os
import random
import shutil
import subprocess
import sys

import numpy as np

B = os.path.expanduser(next((a for a in sys.argv[1:] if not a.startswith("--")), "~/bench071"))
W = os.path.join(B, "world", "full", "simulation")
OUT = os.path.join(B, "strains")
ENV = os.path.expanduser("~/micromamba/envs/protal-db-build")
ART = os.path.join(ENV, "bin", "art_illumina")
PBSIM = os.path.join(ENV, "bin", "pbsim")
SPECIES_BACTERIA, SPECIES_ARCHAEA, STRAINS, BACKGROUND = 10, 2, 8, 15
SEED = 2026
BASES = np.frombuffer(b"ACGT", dtype=np.uint8)


def read_fasta(path):
    opener = gzip.open if path.endswith(".gz") else open
    contigs, name, seq = [], None, []
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    contigs.append((name, "".join(seq)))
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip().upper())
    if name is not None:
        contigs.append((name, "".join(seq)))
    return contigs


def write_fasta(path, contigs):
    with open(path, "w") as fh:
        for name, seq in contigs:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")


def coalescent(rng, n, depth):
    """A random Kingman tree of n leaves, root height `depth`: (children of each node, height of each node, branch
    length to the parent); leaves 0..n-1, the root last."""
    nodes = list(range(n))
    height = {i: 0.0 for i in range(n)}
    children = {}
    t, nxt = 0.0, n
    while len(nodes) > 1:
        k = len(nodes)
        t += rng.expovariate(k * (k - 1) / 2)
        a, b = rng.sample(nodes, 2)
        nodes.remove(a)
        nodes.remove(b)
        children[nxt] = (a, b)
        height[nxt] = t
        nodes.append(nxt)
        nxt += 1
    root = nodes[0]
    scale = depth / height[root]
    for k in height:
        height[k] *= scale
    parent = {c: p for p, cs in children.items() for c in cs}
    branch = {k: height[parent[k]] - height[k] for k in parent}
    return children, height, branch, root


def newick(children, branch, root, names):
    def node(k):
        if k not in children:
            return f"{names[k]}:{branch[k]:.6f}"
        a, b = children[k]
        inner = f"({node(a)},{node(b)})"
        return inner + (f":{branch[k]:.6f}" if k in branch else "")
    return node(root) + ";"


def mutate(nrng, seq_arr, rate_arr, length):
    """seq_arr (uint8 ACGT) with each base substituted with probability rate * length."""
    hits = np.nonzero(nrng.random(len(seq_arr)) < rate_arr * length)[0]
    out = seq_arr.copy()
    for i in hits:
        choices = BASES[BASES != out[i]]
        out[i] = choices[nrng.integers(0, 3)]
    return out


def hifi_qualities(path):
    """Rewrites the FASTQ at path (gzip) with every base quality Q30."""
    tmp = path + ".partial"
    with gzip.open(path, "rt") as src, gzip.open(tmp, "wt", compresslevel=3) as dst:
        for i, line in enumerate(src):
            dst.write("?" * len(line.rstrip("\n")) + "\n" if i % 4 == 3 else line)
    os.replace(tmp, path)


def pbsim_model(name):
    for folder in (sorted(glob.glob(os.path.join(ENV, "share", "pbsim*", "data"))) +
                   sorted(glob.glob(os.path.join(ENV, "share", "pbsim*"))) + [os.path.join(ENV, "data")]):
        for candidate in (name, name + ".model"):
            path = os.path.join(folder, candidate)
            if os.path.isfile(path):
                return path
    sys.exit(f"pbsim3 model {name} not found under {ENV}/share or {ENV}/data")


def main():
    rng = random.Random(SEED)
    nrng = np.random.default_rng(SEED)
    with open(os.path.join(B, "world", "unknown.txt")) as fh:
        unknown = {line.strip().split(";")[-1].removeprefix("s__") for line in fh if line.strip()}
    reps = []
    with open(os.path.join(W, "genomes.tsv")) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["gtdb_representative"] == "t":
                reps.append(row)
    species_of = lambda r: r["gtdb_taxonomy"].split(";")[-1].removeprefix("s__")
    known = [r for r in reps if species_of(r) not in unknown]
    bacteria = [r for r in known if r["gtdb_taxonomy"].startswith("d__Bacteria")]
    archaea = [r for r in known if r["gtdb_taxonomy"].startswith("d__Archaea")]
    chosen = rng.sample(bacteria, SPECIES_BACTERIA) + rng.sample(archaea, SPECIES_ARCHAEA)
    chosen_species = {species_of(r) for r in chosen}
    background_pool = [r for r in reps if species_of(r) not in chosen_species]
    rate = {}
    with open(os.path.join(W, "gene_rates.tsv")) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            rate[row["marker"]] = float(row["rate"])
    markers = {}
    with open(os.path.join(W, "marker_positions.tsv")) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            markers.setdefault(row["accession"], []).append(row)

    shutil.rmtree(OUT, ignore_errors=True)
    for d in ("genomes", "reads", "tmp"):
        os.makedirs(os.path.join(OUT, d))
    samples = [f"str_s{j + 1}" for j in range(STRAINS)]
    strains_rows, trees = [], []
    for si, rep in enumerate(chosen):
        sp = species_of(rep)
        contigs = read_fasta(rep["fasta_path"])
        arrays, rates = [], []
        for name, seq in contigs:
            arr = np.frombuffer(seq.encode(), dtype=np.uint8).copy()
            r = np.ones(len(arr))
            for m in markers.get(rep["accession"], []):
                if m["contig"] == name:
                    r[int(m["start"]) - 1:int(m["end"])] = rate.get(m["marker"], 1.0)
            r[~np.isin(arr, BASES)] = 0  # Ns stay
            arrays.append(arr)
            rates.append(r)
        depth = rng.uniform(0.004, 0.012)
        children, height, branch, root = coalescent(rng, STRAINS, depth)
        seqs = {root: arrays}
        order = sorted(children, key=lambda k: -height[k])  # parents before children
        for k in order:
            for c in children[k]:
                seqs[c] = [mutate(nrng, a, r, branch[c]) for a, r in zip(seqs[k], rates)]
        names = {j: f"STR{si + 1:02d}_{j + 1}" for j in range(STRAINS)}
        trees.append((sp, newick(children, branch, root, names), f"{depth:.4f}"))
        cov = math.exp(rng.uniform(math.log(2), math.log(40)))
        for j in range(STRAINS):
            path = os.path.join(OUT, "genomes", names[j] + ".fna")
            write_fasta(path, [(f"{names[j]}_contig{i + 1}", arr.tobytes().decode())
                               for i, arr in enumerate(seqs[j])])
            strains_rows.append(dict(species=sp, taxonomy=rep["gtdb_taxonomy"], representative=rep["accession"],
                                     strain=names[j], sample=samples[j], genome=path,
                                     coverage=f"{cov * rng.uniform(0.4, 1.6):.2f}", role="strain"))
        print(f"{sp}: 8 strains, root depth {depth:.4f}, coverage about {cov:.0f}x", flush=True)
    for j, sample in enumerate(samples):
        for r in rng.sample(background_pool, BACKGROUND):
            path = os.path.join(OUT, "tmp", f"{sample}_{r['accession']}.fna")
            write_fasta(path, read_fasta(r["fasta_path"]))
            strains_rows.append(dict(species=species_of(r), taxonomy=r["gtdb_taxonomy"], representative=r["accession"],
                                     strain=r["accession"], sample=sample, genome=path,
                                     coverage=f"{rng.uniform(0.5, 5):.2f}", role="background"))
    with open(os.path.join(OUT, "strains.tsv"), "w") as fh:
        cols = list(strains_rows[0])
        fh.write("\t".join(cols) + "\n")
        fh.writelines("\t".join(str(r[c]) for c in cols) + "\n" for r in strains_rows)
    with open(os.path.join(OUT, "trees.tsv"), "w") as fh:
        fh.write("species\troot_depth\tnewick\n")
        fh.writelines(f"{sp}\t{depth}\t{nwk}\n" for sp, nwk, depth in trees)

    pb_model, ont_model = pbsim_model("ERRHMM-SEQUEL"), pbsim_model("QSHMM-ONT-HQ")
    for j, sample in enumerate(samples):
        rows = [r for r in strains_rows if r["sample"] == sample]
        tmp = os.path.join(OUT, "tmp", sample)
        os.makedirs(tmp, exist_ok=True)
        r1, r2 = [], []
        for k, r in enumerate(rows):
            prefix = os.path.join(tmp, f"pe{k}_")
            subprocess.run([ART, "-ss", "HSXt", "-i", r["genome"], "-p", "-l", "150", "-f", r["coverage"], "-m", "350",
                            "-s", "50", "-na", "-rs", str(SEED * 100 + j * 50 + k), "-o", prefix],
                           check=True, stdout=subprocess.DEVNULL)
            r1.append(prefix + "1.fq")
            r2.append(prefix + "2.fq")
        for parts, mate in ((r1, "R1"), (r2, "R2")):
            with gzip.open(os.path.join(OUT, "reads", f"{sample}_{mate}.fq.gz"), "wt", compresslevel=3) as out:
                for part in parts:
                    with open(part) as fh:
                        shutil.copyfileobj(fh, out)
        for kind, method, model, extra in (
                ("pb", "errhmm", pb_model, ["--length-mean", "15000", "--length-sd", "3000", "--accuracy-mean", "0.999"]),
                ("ont", "qshmm", ont_model, ["--length-mean", "8000", "--length-sd", "6000", "--accuracy-mean", "0.97",
                                             "--difference-ratio", "39:24:36"])):
            fastqs = []
            for k, r in enumerate(rows):
                prefix = os.path.join(tmp, f"{kind}{k}")
                subprocess.run([PBSIM, "--strategy", "wgs", "--method", method, f"--{method}", model, "--depth",
                                r["coverage"], "--genome", r["genome"], "--prefix", prefix,
                                "--seed", str(SEED * 100 + j * 50 + k), *extra],
                               check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                fastqs += sorted(glob.glob(prefix + "_*.fq*") + glob.glob(prefix + "_*.fastq*"))
            with gzip.open(os.path.join(OUT, "reads", f"{sample}_{kind}.fq.gz"), "wt", compresslevel=3) as out:
                for part in fastqs:
                    opener = gzip.open if part.endswith(".gz") else open
                    with opener(part, "rt") as fh:
                        shutil.copyfileobj(fh, out)
            if kind == "pb":
                hifi_qualities(os.path.join(OUT, "reads", f"{sample}_{kind}.fq.gz"))
        shutil.rmtree(tmp)
        print(f"{sample}: reads written", flush=True)
    print("done")


if __name__ == "__main__":
    if "--pb-qualities" in sys.argv:
        for path in sorted(glob.glob(os.path.join(OUT, "reads", "*_pb.fq.gz"))):
            hifi_qualities(path)
            print(f"{path}: Q30")
    else:
        main()
