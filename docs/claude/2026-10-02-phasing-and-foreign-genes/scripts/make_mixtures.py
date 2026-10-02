#!/usr/bin/env python3
"""make_mixtures.py SCRIPTS_DIR [BENCH_DIR] - long-read samples of strain mixtures, for the phasing benchmark.

From the strain benchmark of 2026-10-02 (BENCH_DIR, default ~/bench071; ../2026-10-02-v072-benchmark/scripts/
make_strains.py: 12 species of 8 strains each, of known SNPs and tree, strains.tsv): samples whose species each hold two
or three of their strains at known shares, and, as in that benchmark, samples of one strain each:

  pure1..pure8   strain j of every species
  mix1..mix4     strains (1, 2), (3, 4), (5, 6), (7, 8) at 70:30
  mix5, mix6     strains (1, 5), (2, 6) at 50:50
  mix7           strains 3, 7, 8 at 50:30:20
  mix8           strains 4, 5 at 85:15

A species' total coverage in every sample is the mean of its 8 strains' coverages in that benchmark (2 to 40x), each
sample also has the 15 background genomes of the benchmark's sample of the same number (pure j / mix j: str_s<j>).
Reads as protal's collector draws them (collect_training_data.long_read_sample, from SCRIPTS_DIR): PacBio HiFi by
hifi_reads.py (15 +- 3 kb, Q by length + SD 3) and Nanopore by pbsim3 (QSHMM-ONT-HQ, 8 +- 6 kb, 97%, 39:24:36), one
seed per sample and kind. Writes BENCH_DIR/phasing/: reads/<sample>_<pb|ont>.fq.gz, truth.tsv (sample, species,
strain, share, coverage, genome).
"""
import csv
import os
import sys
from concurrent.futures import ProcessPoolExecutor

SCRIPTS = os.path.expanduser(sys.argv[1])
B = os.path.expanduser(sys.argv[2] if len(sys.argv) > 2 else "~/bench071")
S = os.path.join(B, "strains")
OUT = os.path.join(B, "phasing")
ENV = os.path.expanduser("~/micromamba/envs/protal-db-build")
sys.path.insert(0, SCRIPTS)
import collect_training_data as collect  # noqa: E402

DESIGNS = [(f"pure{j}", [(j, 1.0)]) for j in range(1, 9)] + [
    ("mix1", [(1, .7), (2, .3)]), ("mix2", [(3, .7), (4, .3)]), ("mix3", [(5, .7), (6, .3)]),
    ("mix4", [(7, .7), (8, .3)]), ("mix5", [(1, .5), (5, .5)]), ("mix6", [(2, .5), (6, .5)]),
    ("mix7", [(3, .5), (7, .3), (8, .2)]), ("mix8", [(4, .85), (5, .15)])]
SETUPS = {"pb": "hifi:15000:3000:3", "ont": "qshmm:QSHMM-ONT-HQ:8000:6000:0.97:39/24/36"}


def genome_length(path):
    with open(path) as fh:
        return sum(len(line.strip()) for line in fh if not line.startswith(">"))


def run(task):
    return task["out"], collect.long_read_sample(task)


def main():
    rows = list(csv.DictReader(open(os.path.join(S, "strains.tsv")), delimiter="\t"))
    strains = {}  # species -> strain number -> row
    coverage = {}
    for r in rows:
        if r["role"] == "strain":
            strains.setdefault(r["species"], {})[int(r["strain"].rsplit("_", 1)[1])] = r
    for species, of in strains.items():
        coverage[species] = sum(float(r["coverage"]) for r in of.values()) / len(of)
    background = {}
    for r in rows:
        if r["role"] == "background":
            background.setdefault(r["sample"], []).append(r)
    os.makedirs(os.path.join(OUT, "reads"), exist_ok=True)
    truth, tasks = [], []
    for index, (sample, mix) in enumerate(DESIGNS):
        number = int("".join(c for c in sample if c.isdigit()))
        genomes = []
        for species in sorted(strains):
            for strain, share in mix:
                r = strains[species][strain]
                cov = coverage[species] * share
                genomes.append((r["genome"], cov))
                truth.append(dict(sample=sample, species=species, strain=r["strain"], share=share,
                                  coverage=f"{cov:.3f}", genome=r["genome"]))
        for r in background[f"str_s{number}"]:
            genomes.append((r["genome"], float(r["coverage"])))
        weighted = [{"genome": os.path.basename(g), "fasta": g, "weight": cov * genome_length(g)} for g, cov in genomes]
        bases = int(sum(w["weight"] for w in weighted))
        for kind, setup in SETUPS.items():
            out = os.path.join(OUT, "reads", f"{sample}_{kind}.fq.gz")
            if os.path.exists(out):
                continue
            tasks.append({"sample": sample, "tmp": os.path.join(OUT, "tmp", f"{sample}_{kind}"), "out": out,
                          "seed": 7000 + 10 * index + (kind == "ont"), "setup": collect.parse_long_setup(setup),
                          "bases": bases, "genomes": weighted, "pbsim": os.path.join(ENV, "bin", "pbsim"),
                          "model": os.path.join(ENV, "data", "QSHMM-ONT-HQ.model")})
    with open(os.path.join(OUT, "truth.tsv"), "w") as fh:
        cols = list(truth[0])
        fh.write("\t".join(cols) + "\n")
        fh.writelines("\t".join(str(t[c]) for c in cols) + "\n" for t in truth)
    with ProcessPoolExecutor(max_workers=4) as pool:
        for out, error in pool.map(run, tasks):
            print(out, error or "ok", flush=True)


if __name__ == "__main__":
    main()
