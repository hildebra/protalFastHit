#!/usr/bin/env python3
"""Pick other strains of GTDB species than their representatives, to simulate training data from, and
download them from NCBI.

GTDB distributes whole genomes of its species representatives only, and protal's database holds
their marker genes: reads simulated from a representative match the reference better than those of
the strains in real samples, which differ from it by up to a few percent. GTDB's metadata
(bac120/ar53_metadata_r<R>.tsv.gz) lists every genome placed in a species, with its NCBI accession
and quality. This picks, per domain in GTDB's proportions:
- --species species with non-representative genomes that pass the quality filters, and up to
  --per_species of those genomes each (drawn at random), to download;
- --rep_only_species more species, simulated from their representative only (species with one
  genome are a large part of GTDB, and strains close to the reference occur, too).

Written to OUT/:
  strain_accessions.txt   NCBI accessions of the genomes to download
  strain_genomes.tsv      accession, species, lineage, CheckM2 completeness and contamination
  simulation_species.txt  both sets of species: build_gtdb_database.py --simulate-species
  genomes/                with --download: the FASTAs (NCBI datasets CLI), for --extra-genomes

The simulator draws species uniformly from its genome table, then one of their genomes, so a species
with k downloaded strains is simulated from a strain k/(k+1) of the time.

  python3 scripts/gtdb_strain_genomes.py --gtdb GTDB_DIR -o strains --download
  python3 scripts/build_gtdb_database.py --gtdb GTDB_DIR --outdir OUT \\
      --extra-genomes strains/genomes --simulate-species strains/simulation_species.txt ...
"""

import argparse
import collections
import gzip
import os
import random
import shutil
import subprocess
import sys


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--gtdb", required=True, help="extracted GTDB release directory (with the *_metadata_r<R> files)")
    p.add_argument("-o", "--out", required=True)
    p.add_argument("--release", help="GTDB release number (default: detected)")
    p.add_argument("--species", type=int, default=4000, help="species to download strains of (default 4000)")
    p.add_argument("--per_species", type=int, default=2, help="strains per species at most (default 2)")
    p.add_argument("--rep_only_species", type=int, default=1000,
                   help="further species, simulated from their representative only (default 1000: with the "
                        "other defaults, a simulated species is another strain about 53%% of the time)")
    p.add_argument("--min_completeness", type=float, default=90.0, help="CheckM2 completeness, %% (default 90)")
    p.add_argument("--max_contamination", type=float, default=5.0, help="CheckM2 contamination, %% (default 5)")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--download", action="store_true", help="download the genomes with the NCBI datasets CLI")
    p.add_argument("--datasets", default="datasets", help="NCBI datasets binary (default: on PATH)")
    return p.parse_args(argv)


def detect_release(gtdb):
    found = {name.split("_metadata_r")[1].split(".")[0] for name in os.listdir(gtdb)
             if "_metadata_r" in name and name.split("_")[0] in ("bac120", "ar53")}
    if len(found) != 1:
        sys.exit(f"cannot detect one GTDB release in {gtdb} (found {sorted(found) or 'none'}); pass --release")
    return found.pop()


def read_metadata(gtdb, release):
    """Every genome: (NCBI accession, species, lineage, is representative, completeness, contamination)."""
    genomes = []
    for mset in ("bac120", "ar53"):
        path = next((os.path.join(gtdb, f"{mset}_metadata_r{release}{ext}") for ext in (".tsv", ".tsv.gz")
                     if os.path.exists(os.path.join(gtdb, f"{mset}_metadata_r{release}{ext}"))), None)
        if path is None:
            continue
        with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            col = {name: i for i, name in enumerate(header)}
            completeness = col.get("checkm2_completeness", col.get("checkm_completeness"))
            contamination = col.get("checkm2_contamination", col.get("checkm_contamination"))
            for name in ("accession", "gtdb_taxonomy", "gtdb_representative"):
                if name not in col:
                    sys.exit(f"{path} has no {name} column")
            for line in fh:
                f = line.rstrip("\n").split("\t")
                lineage = f[col["gtdb_taxonomy"]]
                accession = f[col["accession"]]
                genomes.append((accession[3:] if accession[:3] in ("RS_", "GB_") else accession, lineage.split(";")[-1],
                                lineage, f[col["gtdb_representative"]] == "t",
                                float(f[completeness]) if completeness is not None and f[completeness] not in ("", "none") else 100.0,
                                float(f[contamination]) if contamination is not None and f[contamination] not in ("", "none") else 0.0))
    if not genomes:
        sys.exit(f"no bac120/ar53_metadata_r{release}.tsv[.gz] in {gtdb}")
    return genomes


def pick(genomes, opts):
    rng = random.Random(opts.seed)
    strains = collections.defaultdict(list)
    domain, lineage = {}, {}
    for acc, species, lin, is_rep, comp, cont in genomes:
        domain[species], lineage[species] = lin.split(";")[0], lin
        if not is_rep and comp >= opts.min_completeness and cont <= opts.max_contamination:
            strains[species].append((acc, comp, cont))
    with_strains = collections.defaultdict(list)
    without = collections.defaultdict(list)
    for species in sorted(domain):
        (with_strains if strains.get(species) else without)[domain[species]].append(species)
    chosen, rep_only = [], []
    for d in sorted(set(domain.values())):
        share = sum(1 for s in domain if domain[s] == d) / len(domain)
        picked_here = set(rng.sample(with_strains[d], min(len(with_strains[d]), round(opts.species * share))))
        chosen += sorted(picked_here)
        rest = [s for s in with_strains[d] if s not in picked_here] + without[d]
        rep_only += rng.sample(rest, min(len(rest), round(opts.rep_only_species * share)))
    picked = []
    for species in sorted(chosen):
        options = sorted(strains[species])
        for acc, comp, cont in rng.sample(options, min(opts.per_species, len(options))):
            picked.append((acc, species, lineage[species], comp, cont))
    return picked, sorted(chosen), sorted(rep_only), domain


def download(opts, accessions_file):
    datasets = shutil.which(opts.datasets)
    if not datasets:
        sys.exit(f"the NCBI datasets CLI ({opts.datasets}) is not on PATH: install it "
                 "(https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/) "
                 "or download strain_accessions.txt another way")
    zip_path = os.path.join(opts.out, "strains.zip")
    genomes = os.path.join(opts.out, "genomes")
    for command in ([datasets, "download", "genome", "accession", "--inputfile", accessions_file,
                     "--include", "genome", "--dehydrated", "--filename", zip_path],
                    ["unzip", "-o", "-q", zip_path, "-d", genomes],
                    [datasets, "rehydrate", "--directory", genomes]):
        print("+ " + " ".join(command), flush=True)
        if subprocess.run(command).returncode != 0:
            sys.exit(f"failed: {' '.join(command)}")
    found = sum(1 for _root, _dirs, files in os.walk(genomes) for f in files if f.endswith((".fna", ".fna.gz")))
    print(f"{found} genome FASTAs in {genomes}")


def main(argv=None):
    opts = parse_args(argv)
    release = opts.release or detect_release(opts.gtdb)
    genomes = read_metadata(opts.gtdb, release)
    picked, chosen, rep_only, domain = pick(genomes, opts)
    os.makedirs(opts.out, exist_ok=True)
    accessions = os.path.join(opts.out, "strain_accessions.txt")
    with open(accessions, "w") as fh:
        fh.writelines(acc + "\n" for acc, *_ in picked)
    with open(os.path.join(opts.out, "strain_genomes.tsv"), "w") as fh:
        fh.write("accession\tspecies\tlineage\tcheckm2_completeness\tcheckm2_contamination\n")
        fh.writelines(f"{acc}\t{sp}\t{lin}\t{comp}\t{cont}\n" for acc, sp, lin, comp, cont in picked)
    with open(os.path.join(opts.out, "simulation_species.txt"), "w") as fh:
        fh.writelines(s + "\n" for s in sorted(chosen + rep_only))
    counts = collections.Counter(domain[s] for s in chosen)
    per = collections.Counter(sp for _, sp, *_ in picked)
    other = sum(k / (k + 1) for k in per.values()) / max(1, len(chosen) + len(rep_only))
    print(f"GTDB r{release}: {len(genomes)} genomes of {len(domain)} species")
    print(f"strains: {len(picked)} genomes of {len(chosen)} species ("
          + ", ".join(f"{d} {n}" for d, n in sorted(counts.items())) + f") -> {accessions}")
    print(f"simulation pool: {len(chosen) + len(rep_only)} species ({len(rep_only)} from the representative only); "
          f"a simulated species is another strain than the representative {100 * other:.0f}% of the time")
    if opts.download:
        download(opts, accessions)
    else:
        print("download them with --download (NCBI datasets CLI), or: datasets download genome accession "
              f"--inputfile {accessions} --include genome --dehydrated --filename strains.zip; unzip strains.zip "
              "-d genomes; datasets rehydrate --directory genomes")


if __name__ == "__main__":
    main()
