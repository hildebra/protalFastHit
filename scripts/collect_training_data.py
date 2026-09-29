#!/usr/bin/env python3
"""Simulate metagenomes and profile them with protal, to train its presence model.

For every design point, a read setup (length, ART profile, fragment size) and a sequencing depth,
simulate_metagenomes draws random communities from a genome table and protal profiles them
against a database, knowing the true species. The training dumps of all samples
(<profile>.truth_annotated: every taxon protal saw, its features and whether it was present) are
joined into one table, with meta_* columns saying where each row comes from (meta_domain from
the genome table's lineages; with --novel_species and --taxonomy, which rows share a genus with a
species the database lacks, and which present species were simulated from another genome than
the database's reference). Train on it with

    python3 scripts/random_forest_cmdline.py --truth-file OUT/training_data.tsv --output-prefix OUT/model

A genome table with species the database lacks gives the negatives that matter most: a relative
the database has picks up their reads. Archaea (--archaea) have fewer marker genes than bacteria
and need to be in the training data, too. Points already done are skipped, so a run can be resumed.

usage: collect_training_data.py --db DB --genome_table genomes.tsv -o OUT [options]
"""

import argparse
import collections
import csv
import glob
import os
import random
import subprocess
import sys


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--db", required=True, help="protal database")
    p.add_argument("--genome_table", required=True,
                   help="simulate_metagenomes genome table (accession, GTDB taxonomy, FASTA path)")
    p.add_argument("-o", "--out", required=True, help="output directory")
    p.add_argument("--protal", default="protal", help="protal binary (default: protal on PATH)")
    p.add_argument("--simulator", default="simulate_metagenomes", help="simulate_metagenomes binary")
    p.add_argument("--samples", type=int, default=4, help="samples per design point (default: 4)")
    p.add_argument("--read_pairs", default="5000,20000,100000,500000",
                   help="comma-separated read pairs per sample, one design point each")
    p.add_argument("--read_setups", default="100:HS20:300:40,150:HS25:350:50,250:MSv3:550:50",
                   help="comma-separated LENGTH:ART_PROFILE:FRAGMENT_MEAN:FRAGMENT_SD, one design point each")
    p.add_argument("--species_per_sample", default="5-30", help="species per sample, N or MIN-MAX")
    p.add_argument("--archaea", type=int, default=0, help="archaeal species per sample (default: 0)")
    p.add_argument("--congeners", type=int, default=0,
                   help="species of one genus in every sample of a design point, the genus drawn per point among "
                        "those with that many species (default: 0). Species are otherwise drawn uniformly, so among "
                        "many genera relatives hardly ever share a sample, while in real samples they often do")
    p.add_argument("-t", "--threads", type=int, default=4)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--novel_species",
                   help="species the database lacks, one per line (e.g. those a training database leaves out): "
                        "meta_novel_species counts them in each sample, meta_novel_congener marks the taxa of their "
                        "genera, which their reads land on")
    p.add_argument("--taxonomy", help="internal_taxonomy.dmp of the database: meta_rep_genome says whether a present "
                                      "species was simulated from its representative genome (the database's "
                                      "reference, 1) or from another strain (0)")
    return p.parse_args(argv)


META_COLUMNS = ["meta_design", "meta_sample", "meta_read_length", "meta_read_pairs", "meta_domain",
                "meta_novel_species", "meta_novel_congener", "meta_rep_genome"]


def large_genera(genome_table, size):
    """Genera with at least `size` species in the genome table, sorted."""
    species = collections.defaultdict(set)
    with open(genome_table) as fh:
        for line in fh:
            lineage = next((f for f in line.rstrip("\n").split("\t") if f.startswith("d__") and ";s__" in f), None)
            if lineage:
                ranks = lineage.split(";")
                species[ranks[5][3:]].add(ranks[-1])
    return sorted(g for g, s in species.items() if len(s) >= size)


def genus_of(species):
    """The genus of a GTDB species name: s__Bacillus_A cereus -> Bacillus_A."""
    return species[3:].split(" ")[0] if species.startswith("s__") else ""


def read_list(path):
    with open(path) as fh:
        names = {line.rstrip("\n").split("\t")[0].strip() for line in fh}
    return {n if n.startswith("s__") else "s__" + n for n in names if n and not n.startswith("#")}


def representatives(taxonomy):
    """species name -> accession of its representative genome, from internal_taxonomy.dmp."""
    with open(taxonomy) as fh:
        header = next(fh).rstrip("\n").split("\t")
        name, rank, rep = header.index("name"), header.index("rank"), header.index("rep_genome")
        return {f[name]: f[rep] for f in (line.rstrip("\n").split("\t") for line in fh) if f[rank] == "species"}


def simulated_genomes(point_dir):
    """sample -> {species name: [genome, ...]} from the simulator's manifest."""
    genomes = {}
    with open(os.path.join(point_dir, "sim", "manifest.tsv")) as fh:
        header = next(fh).rstrip("\n").split("\t")
        for line in fh:
            row = dict(zip(header, line.rstrip("\n").split("\t")))
            species = row["taxonomy"].split(";")[-1]
            genomes.setdefault(row["sample"], {}).setdefault(species, []).append(row["genome"])
    return genomes


def species_domains(genome_table):
    """GTDB species name (s__Genus species, as protal names taxa) -> domain, from the lineages in the table."""
    domains = {}
    with open(genome_table) as fh:
        for line in fh:
            for field in line.rstrip("\n").split("\t"):
                if field.startswith("d__") and ";s__" in field:
                    ranks = field.split(";")
                    domains[ranks[-1].strip()] = ranks[0][3:]
                    break
    return domains


def run(command, log):
    with open(log, "w") as fh:
        rc = subprocess.run(command, stdout=fh, stderr=subprocess.STDOUT).returncode
    if rc != 0:
        sys.exit(f"{command[0]} failed with exit code {rc}; see {log}")


def design_points(opts):
    points = []
    for setup in opts.read_setups.split(","):
        length, profile, fragment_mean, fragment_sd = setup.split(":")
        for pairs in opts.read_pairs.split(","):
            points.append({"name": f"rl{length}_p{pairs}", "read_length": length, "sequencer": profile,
                           "fragment_mean": fragment_mean, "fragment_sd": fragment_sd, "read_pairs": pairs})
    return points


def simulate_and_profile(point, index, opts):
    base = os.path.join(opts.out, "points", point["name"])
    sim, profiles = os.path.join(base, "sim"), os.path.join(base, "protal")
    meta = os.path.join(sim, "protal.meta")
    if not os.path.isfile(meta):
        command = [opts.simulator, "--genome_table", opts.genome_table, "-o", sim, "-n", str(opts.samples),
                   "--sample_prefix", point["name"] + "_s", "--total_read_pairs", point["read_pairs"],
                   "--species_per_sample", opts.species_per_sample, "--read_length", point["read_length"],
                   "--sequencer", point["sequencer"], "--fragment_mean", point["fragment_mean"],
                   "--fragment_stdev", point["fragment_sd"], "--seed", str(opts.seed + index),
                   "-t", str(opts.threads), "--protal_metafile", profiles]
        if opts.archaea > 0:
            command += ["--taxon", f"d__Archaea:{opts.archaea}"]
        if opts.congeners > 0:
            genera = large_genera(opts.genome_table, opts.congeners)
            if not genera:
                sys.exit(f"no genus in {opts.genome_table} has {opts.congeners} species (--congeners)")
            command += ["--genus", "g__" + random.Random(opts.seed * 1000 + index).choice(genera) + f":{opts.congeners}"]
        if opts.archaea > 0 or opts.congeners > 0:
            command += ["--pick_random_demand_if_fail"]
        os.makedirs(sim, exist_ok=True)
        run(command, os.path.join(base, "simulate.log"))
    dumps = sorted(glob.glob(os.path.join(profiles, "**", "*.truth_annotated"), recursive=True))
    if len(dumps) < opts.samples:
        run([opts.protal, "--db", opts.db, "--map", meta, "-t", str(opts.threads), "--no_strains", "--no_qcmsa"],
            os.path.join(base, "protal.log"))
        dumps = sorted(glob.glob(os.path.join(profiles, "**", "*.truth_annotated"), recursive=True))
    if len(dumps) != opts.samples:
        sys.exit(f"{point['name']}: expected {opts.samples} training dumps in {profiles}, found {len(dumps)}")
    return dumps


def main(argv=None):
    opts = parse_args(argv)
    os.makedirs(opts.out, exist_ok=True)
    domains = species_domains(opts.genome_table)
    novel = read_list(opts.novel_species) if opts.novel_species else set()
    reps = representatives(opts.taxonomy) if opts.taxonomy else {}
    header, rows = None, 0
    table = os.path.join(opts.out, "training_data.tsv")
    with open(table + ".partial", "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        for index, point in enumerate(design_points(opts)):
            present = absent = strains = congeners = 0
            dumps = simulate_and_profile(point, index, opts)
            genomes = simulated_genomes(os.path.join(opts.out, "points", point["name"]))
            for dump in dumps:
                sample = os.path.basename(dump).split(".profile")[0]
                in_sample = genomes.get(sample, {})
                novel_here = [s for s in in_sample if s in novel]
                novel_genera = {genus_of(s) for s in novel_here}
                with open(dump, newline="") as fh:
                    reader = csv.reader(fh, delimiter="\t")
                    dump_header = next(reader)
                    if header is None:
                        header = dump_header
                        writer.writerow(META_COLUMNS + header)
                    elif dump_header != header:
                        sys.exit(f"{dump} has other columns than the dumps before it (another protal version?)")
                    name, truth = header.index("taxon_name"), header.index("truth")
                    for row in reader:
                        taxon = row[name]
                        is_present = row[truth].lower() in ("1", "true")
                        congener = genus_of(taxon) in novel_genera
                        rep = ""
                        if is_present and reps.get(taxon) and taxon in in_sample:
                            rep = "1" if all(g == reps[taxon] for g in in_sample[taxon]) else "0"
                        writer.writerow([point["name"], sample, point["read_length"], point["read_pairs"],
                                         domains.get(taxon, "unknown"), len(novel_here), int(congener), rep] + row)
                        rows += 1
                        present += is_present
                        absent += not is_present
                        strains += rep == "0"
                        congeners += congener and not is_present
            print(f"{point['name']}: {present} present ({strains} from other genomes than the representative) and "
                  f"{absent} absent taxa ({congeners} congeners of species the database lacks) in {opts.samples} "
                  "samples", flush=True)
    os.replace(table + ".partial", table)
    print(f"{rows} taxa in {table}")


if __name__ == "__main__":
    main()
